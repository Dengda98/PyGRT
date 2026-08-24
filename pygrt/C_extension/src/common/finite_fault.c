/**
 * @file   finite_fault.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * Coulomb 格式有限断层：读入、衍生量与几何剖分
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include "grt/common/finite_fault.h"
#include "grt/common/util.h"

#define COULOMB_HEADER_MIN_TOKENS 11
#define COULOMB_HEADER_MAX_TOKENS 12
#define COULOMB_HEADER_TOKEN_SIZE 32

/**
 * 检查有限断层的走向长度、倾角和深度范围
 *
 * @param[in]  f       有限断层结构体
 * @param[in]  where   错误位置描述
 */
static void check_fault_geometry(const FINITE_FAULT *f, const char *where)
{
    if(hypot(f->east_end - f->east_begin, f->north_end - f->north_begin) <= 0.0){
        GRTRaiseError("%s: fault along-strike length must be positive.", where);
    }
    if(f->dip <= 0.0 || f->dip > 90.0){
        GRTRaiseError("%s: dip (%.6g deg) must be in (0, 90].", where, f->dip);
    }
    if(!(f->bot > f->top)){
        GRTRaiseError("%s: bot (%.6g km) must be greater than top (%.6g km).", where, f->bot, f->top);
    }
}


/**
 * 由 Coulomb 分量建立有限断层的走向、滑动角和滑动量
 *
 * @param[in,out]  f   有限断层结构体
 */
static void set_fault_derived(FINITE_FAULT *f)
{
    f->strike = 1.0 / DEG1 * atan2(f->east_end - f->east_begin, f->north_end - f->north_begin);
    if((f->right_lateral == 0.0) && (f->reverse == 0.0)){
        f->rake = GRT_FINITE_FAULT_UNDEFINED_RAKE;
    } else {
        f->rake = 1.0 / DEG1 * atan2(f->reverse, -f->right_lateral);
    }
    f->slip = KODE_IS_POINT(f->kode) ? 0.0 : hypot(f->right_lateral, f->reverse);
}

/**
 * 按 Kode 解释文件第 7、8 列
 *
 * @param[out]  f             保存解释结果的有限断层结构体
 * @param[in]   rake_format   是否使用表头标识的 rake/net slip 格式
 * @param[in]   where         错误位置描述
 */
static void set_fault_components(FINITE_FAULT *f, bool rake_format, const char *where)
{
    f->rake_format = rake_format;
    if(rake_format && f->kode != KODE_RTLAT_REVERSE){
        GRTRaiseError("%s: Coulomb rake/net slip format only supports Kode=100.", where);
    }

    switch(f->kode){
        case KODE_RTLAT_REVERSE:
            if(rake_format){
                f->right_lateral = -f->value2 * cos(f->value1 * DEG1);
                f->reverse = f->value2 * sin(f->value1 * DEG1);
            } else {
                f->right_lateral = f->value1;
                f->reverse = f->value2;
            }
            break;
        case KODE_RTLAT_TENSILE:
            f->right_lateral = f->value1;
            f->tensile = f->value2;
            break;
        case KODE_TENSILE_REVERSE:
            f->tensile = f->value1;
            f->reverse = f->value2;
            break;
        case KODE_POINT_DC:
            f->right_lateral = f->value1;
            f->reverse = f->value2;
            break;
        case KODE_POINT_TENSILE_INFLATE:
            f->tensile = f->value1;
            f->inflate = f->value2;
            break;
        default:
            GRTRaiseError("%s: unsupported Coulomb Kode=%u, expected %u, %u, %u, %u or %u.", where,
                f->kode, KODE_RTLAT_REVERSE, KODE_RTLAT_TENSILE, KODE_TENSILE_REVERSE,
                KODE_POINT_DC, KODE_POINT_TENSILE_INFLATE);
    }
}


FINITE_FAULT *grt_finite_fault_load_coulomb(const char *path, size_t *nfault)
{
    if(path == NULL || nfault == NULL){
        GRTRaiseError("path/nfault is NULL.");
    }

    GRTCheckFileExist(path);
    FILE *fp = GRTCheckOpenFile(path, "r");

    char *line = NULL;
    size_t nlen = 0;
    if(grt_getline(&line, &nlen, fp) <= 0){
        fclose(fp);
        GRT_SAFE_FREE_PTR(line);
        GRTRaiseError("read Coulomb fault header of %s failed.", path);
    }

    // Coulomb 表格有 11 个数值列；首行第一个 token 是 ID 列的 # 标记
    char header[COULOMB_HEADER_MAX_TOKENS][COULOMB_HEADER_TOKEN_SIZE] = {{0}};
    int nheader = sscanf(line,
        "%31s %31s %31s %31s %31s %31s %31s %31s %31s %31s %31s %31s",
        header[0], header[1], header[2], header[3], header[4], header[5],
        header[6], header[7], header[8], header[9], header[10], header[11]);
    if(nheader < COULOMB_HEADER_MIN_TOKENS || strcmp(header[0], "#") != 0){
        fclose(fp);
        GRT_SAFE_FREE_PTR(line);
        GRTRaiseError("invalid Coulomb fault header in %s: expected # plus 10 field labels for 11 data columns.", path);
    }

    // 仅第 7 个数据列的完整 token "rake" 表示 rake/net-slip 格式
    bool rake_format = strcmp(header[6], "rake") == 0;
    size_t path_length = strlen(path);
    bool suffix_rake = path_length >= 4 && strcmp(path + path_length - 4, ".inr") == 0;
    if(suffix_rake && !rake_format){
        GRTRaiseWarning(
            "Coulomb file \"%s\" has .inr suffix but the seventh header column is \"%s\", "
            "not the exact token \"rake\"; using the header-defined component format.",
            path, header[6]);
    } else if(!suffix_rake && rake_format){
        GRTRaiseWarning(
            "Coulomb file \"%s\" has the exact token \"rake\" in the seventh header column "
            "but no .inr suffix; using the header-defined rake/net-slip format.", path);
    }

    // 第二行是与 11 个数据列对应的占位行
    if(grt_getline(&line, &nlen, fp) <= 0){
        fclose(fp);
        GRT_SAFE_FREE_PTR(line);
        GRTRaiseError("read Coulomb fault placeholder header of %s failed.", path);
    }

    FINITE_FAULT *faults = NULL;
    size_t n = 0;
    size_t line_number = 2;
    while(grt_getline(&line, &nlen, fp) != -1){
        ++line_number;

        real_t dum1, kode_value;
        real_t east_begin, north_begin, east_end, north_end;
        real_t value1, value2, dip, top, bot;
        int nscan = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &dum1, &east_begin, &north_begin, &east_end, &north_end,
            &kode_value, &value1, &value2, &dip, &top, &bot);
        if(nscan != 11){
            fclose(fp);
            GRT_SAFE_FREE_PTR(line);
            GRT_SAFE_FREE_PTR(faults);
            GRTRaiseError("parse Coulomb fault data at line %zu of %s failed.", line_number, path);
        }

        if(fabs(kode_value - round(kode_value)) > 1e-8 || kode_value < 0.0){
            fclose(fp);
            GRT_SAFE_FREE_PTR(line);
            GRT_SAFE_FREE_PTR(faults);
            GRTRaiseError("invalid Coulomb Kode at line %zu of %s.", line_number, path);
        }

        faults = (FINITE_FAULT *)realloc(faults, sizeof(FINITE_FAULT) * (n + 1));
        FINITE_FAULT *f = faults + n;
        memset(f, 0, sizeof(*f));

        f->east_begin = east_begin;
        f->north_begin = north_begin;
        f->east_end = east_end;
        f->north_end = north_end;
        f->kode = (unsigned int)round(kode_value);
        f->value1 = value1;
        f->value2 = value2;
        f->dip = dip;
        f->top = top;
        f->bot = bot;

        char where[256];
        snprintf(where, sizeof(where), "in %s line %zu", path, line_number);
        set_fault_components(f, rake_format, where);
        check_fault_geometry(f, where);

        set_fault_derived(f);
        n++;
    }

    GRT_SAFE_FREE_PTR(line);
    fclose(fp);

    if(n == 0){
        GRTRaiseError("no fault in %s.", path);
    }

    *nfault = n;
    return faults;
}


FINITE_FAULT *grt_finite_fault_from_option(const char *option, size_t *nfault, real_t *dL, real_t *dW)
{
    if((option == NULL) || (nfault == NULL) || (dL == NULL) || (dW == NULL)){
        GRTRaiseError("finite fault option is incomplete.");
    }

    char *option_copy = strdup(option);
    // 将文件路径和可选的剖分尺寸拆分到独立字符串中
    char *path = strtok(option_copy, "+");
    char *token = strtok(NULL, "+");
    if((path == NULL) || (*path == '\0') || (strtok(NULL, "+") != NULL)){
        GRT_SAFE_FREE_PTR(option_copy);
        GRTRaiseError("Error in finite fault option. expected <fault>[+i<dL>/<dW>]. Use \"-h\" for help.");
    }

    *dL = 0.0;
    *dW = 0.0;
    if(token != NULL){
        // 解析 +i<dL>/<dW>，未提供时保留不剖分标记
        char extra;
        if((token[0] != 'i') || (sscanf(token + 1, "%lf/%lf%c", dL, dW, &extra) != 2)){
            GRT_SAFE_FREE_PTR(option_copy);
            GRTRaiseError("Error in finite fault option. expected +i<dL>/<dW>. Use \"-h\" for help.");
        }
        if((*dL <= 0.0) || (*dW <= 0.0)){
            GRT_SAFE_FREE_PTR(option_copy);
            GRTRaiseError("Error in finite fault option. dL and dW must be positive. Use \"-h\" for help.");
        }
    }

    // 统一由有限断层模块读取文件，调用方只保留数组和剖分尺寸
    FINITE_FAULT *faults = grt_finite_fault_load_coulomb(path, nfault);
    GRT_SAFE_FREE_PTR(option_copy);
    return faults;
}


void grt_finite_fault_free(FINITE_FAULT *faults)
{
    free(faults);
}


void grt_finite_fault_subdiv(
    const FINITE_FAULT *fault, real_t dL, real_t dW,
    real_t *W, real_t *L, size_t *nW, size_t *nL)
{
    if((dL <= 0.0) != (dW <= 0.0)){
        GRTRaiseError("dL and dW must both be positive or both be omitted.");
    }
    check_fault_geometry(fault, "finite fault");

    *W = (fault->bot - fault->top) / sin(DEG1 * fault->dip);
    *L = hypot(fault->east_end - fault->east_begin, fault->north_end - fault->north_begin);
    if(*W <= 0.0 || *L <= 0.0){
        GRTRaiseError("fault along-dip/along-strike length must be positive (W=%.6g, L=%.6g).", *W, *L);
    }
    if(dL <= 0.0){
        *nW = 1;
        *nL = 1;
    } else {
        *nW = GRT_MAX(1, (size_t)ceil(*W / dW));
        *nL = GRT_MAX(1, (size_t)ceil(*L / dL));
    }
}


void grt_finite_fault_subfault(
    const FINITE_FAULT *fault,
    real_t dL, real_t dW, real_t W, real_t L,
    size_t iW, size_t iL,
    FINITE_SUBFAULT *sub)
{
    if((dL <= 0.0) != (dW <= 0.0)){
        GRTRaiseError("dL and dW must both be positive or both be omitted.");
    }
    if(dL <= 0.0){
        // 不剖分时，当前有限断层本身就是唯一的子断层
        dL = L;
        dW = W;
    }

    // 末块可短于 dW/dL，中心取该块中点：i*d + size/2，而不是 (i+0.5)*size
    real_t width = GRT_MIN(dW, W - iW * dW);
    real_t length = GRT_MIN(dL, L - iL * dL);
    if(width <= 0.0 || length <= 0.0){
        GRTRaiseError("nonpositive subfault size at (iW=%zu, iL=%zu).", iW, iL);
    }
    real_t w = iW * dW + 0.5 * width;
    real_t l = iL * dL + 0.5 * length;

    real_t sind = sin(DEG1 * fault->dip);
    real_t cosd = cos(DEG1 * fault->dip);
    real_t sins = sin(DEG1 * fault->strike);
    real_t coss = cos(DEG1 * fault->strike);

    real_t hproj = w * cosd;
    real_t east0 = fault->east_begin + l * sins;
    real_t north0 = fault->north_begin + l * coss;

    sub->width = width;
    sub->length = length;
    sub->depsrc = fault->top + w * sind;
    sub->east = east0 + hproj * coss;
    sub->north = north0 - hproj * sins;
    // slip(m) * width(km) * length(km) → cm^3
    sub->potency = fault->slip * width * length * 1e12;
}
