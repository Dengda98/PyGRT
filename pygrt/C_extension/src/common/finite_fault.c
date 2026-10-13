/**
 * @file   finite_fault.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * sCoulomb (simplified Coulomb) 格式有限断层：读入、衍生量与几何剖分
 *
 */

#include <ctype.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include "grt/common/finite_fault.h"
#include "grt/common/util.h"
#include "grt/common/model.h"
#include "grt/common/radiation.h"
#include "grt/dynamic/signals.h"

#define COULOMB_HEADER_MIN_TOKENS 11  ///< 表头至少包含 # 和十个字段标签
#define COULOMB_HEADER_MAX_TOKENS 12  ///< 允许 dip angle 使用两个 token
#define COULOMB_HEADER_TOKEN_SIZE 32  ///< 表头字段缓冲长度

/**
 * 检查有限断层的走向长度、倾角和深度范围
 *
 * @param[in]  f      有限断层结构体
 * @param[in]  where  错误位置描述
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
 * 由顶边向量建立有限断层的走向
 *
 * @param[in,out]  f  有限断层结构体
 */
static void set_fault_strike(FINITE_FAULT *f)
{
    f->strike = 1.0 / DEG1 * atan2(f->east_end - f->east_begin, f->north_end - f->north_begin);
    if(f->strike < 0.0) {
        f->strike += 360.0;
    }
}

/**
 * 按 Kode 解释文件第 7、8 列
 *
 * @param[out]  f            保存解释结果的有限断层结构体
 * @param[in]   rake_format  是否使用表头标识的 rake/net slip 格式
 * @param[in]   where        错误位置描述
 */
static void set_fault_components(FINITE_FAULT *f, bool rake_format, const char *where)
{
    if(rake_format && f->kode != KODE_RTLAT_REVERSE){
        GRTRaiseError("%s: sCoulomb rake/net slip format only supports Kode=100.", where);
    }

    if(rake_format && fabs(f->value1) > 180.0) {
        GRTRaiseError("%s: rake must be in [-180, 180].", where);
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
            GRTRaiseError("%s: unsupported sCoulomb Kode=%u, expected %u, %u, %u, %u or %u.", where,
                f->kode, KODE_RTLAT_REVERSE, KODE_RTLAT_TENSILE, KODE_TENSILE_REVERSE,
                KODE_POINT_DC, KODE_POINT_TENSILE_INFLATE);
    }
}


/**
 * 将剪切源项统一为非负强度和实际滑动方向，并检查滑动角范围
 *
 * @param[in,out]  term   基本源型分量，非剪切源保持带符号强度
 * @param[in]      where  错误位置描述
 */
static void normalize_shear_term(FINITE_SOURCE_TERM *term, const char *where)
{
    if(term->type != GRT_SYN_DC) {
        return;
    }

    // 负剪切强度改为正值，同时将滑动方向反转 180 度，角度仍位于 [-180, 180]
    real_t *rake = &term->mechanism[2];
    if(term->scale < 0.0) {
        term->scale = -term->scale;
        *rake += *rake > 0.0 ? -180.0 : 180.0;
    }
    bool undefined = term->scale == 0.0 && *rake == GRT_FINITE_FAULT_UNDEFINED_RAKE;
    if(!undefined && fabs(*rake) > 180.0) {
        GRTRaiseError("%s: shear-source rake must be in [-180, 180].", where);
    }
}


/**
 * 按 sCoulomb Kode 建立断层共享的源项或接收机制，矩形源按单位面积保存强度
 *
 * @param[in,out]  fault        已建立走向和原始滑动分量的断层
 * @param[in]      is_source    是否为震源断层
 * @param[in]      rake_format  是否使用 rake/net slip 格式
 * @param[in]      where        错误位置描述
 */
static void set_fault_terms(FINITE_FAULT *fault, bool is_source, bool rake_format, const char *where)
{
    real_t factor = KODE_IS_POINT(fault->kode) ? 1e6 : 1e12;
    GRT_SYN_TYPE types[2] = {GRT_SYN_DC, GRT_SYN_DC};
    real_t strengths[2] = {0}, rakes[2] = {0};
    int ncomponents = 1;

    // 只在此处解释 Kode，动态解和静态解使用同一套源型及符号约定
    switch(fault->kode) {
        case KODE_RTLAT_REVERSE:
        case KODE_POINT_DC:
            strengths[0] = hypot(fault->right_lateral, fault->reverse);
            if(strengths[0] == 0.0) {
                // 无滑动分量时，rake 格式仍可提供显式方向，分量格式则标记为未定义
                rakes[0] = rake_format ? fault->value1 : GRT_FINITE_FAULT_UNDEFINED_RAKE;
            } else {
                rakes[0] = atan2(fault->reverse, -fault->right_lateral) / DEG1;
            }
            break;
        case KODE_RTLAT_TENSILE:
            strengths[0] = fault->right_lateral;
            rakes[0] = 180.0;
            types[1] = GRT_SYN_TS;
            strengths[1] = fault->tensile;
            ncomponents = 2;
            break;
        case KODE_TENSILE_REVERSE:
            types[0] = GRT_SYN_TS;
            strengths[0] = fault->tensile;
            strengths[1] = fault->reverse;
            rakes[1] = 90.0;
            ncomponents = 2;
            break;
        case KODE_POINT_TENSILE_INFLATE:
            types[0] = GRT_SYN_TS;
            strengths[0] = fault->tensile;
            types[1] = GRT_SYN_EX;
            strengths[1] = fault->inflate;
            ncomponents = 2;
            break;
        default:
            GRTRaiseError("Unsupported sCoulomb Kode=%u.", fault->kode);
    }

    // 接收断层只使用显式 rake 定义投影方向，不随 net slip 的符号或大小改变
    if(!is_source && rake_format) {
        rakes[0] = fault->value1;
    }

    fault->nterms = 0;
    for(int i = 0; i < ncomponents; ++i) {
        FINITE_SOURCE_TERM term = {.type = types[i], .scale = is_source ? strengths[i] * factor : 0.0, .with_mu = is_source,
                                   .mechanism = {fault->strike, fault->dip, rakes[i]}};
        normalize_shear_term(&term, where);

        // 接收断层和零滑动的纯剪切源保留 terms[0] 的机制，但不计入参与合成的非零源项
        if(!is_source || strengths[i] == 0.0) {
            if(ncomponents == 1 && term.type == GRT_SYN_DC) {
                fault->terms[0] = term;
            }
            continue;
        }
        fault->terms[fault->nterms++] = term;
        for(int im = 0; im < GRT_SRC_M_NUM; ++im) {
            fault->required_gf_sources[im] |= grt_source_has_component(term.type, im);
        }
    }
}


/**
 * 将无延迟的源时间函数阻尼到格林函数使用的复频率时间域
 * @param[in,out]  stf    源时间函数
 * @param[in]      npts   样本数
 * @param[in]      decay  每个采样点的阻尼因子
 */
static void damp_stf(real_t *stf, int npts, real_t decay)
{
    real_t damping = 1;
    for(int n = 0; n < npts; ++n) {
        stf[n] *= damping;
        damping *= decay;
    }
}


FINITE_FAULT *grt_finite_fault_load_coulomb(const char *path, size_t *nfault, bool is_source, real_t dt, real_t stf_decay, const char *stf_option)
{
    if(path == NULL || nfault == NULL){
        GRTRaiseError("path/nfault is NULL.");
    }

    GRTCheckFileExist(path);
    FILE *fp = GRTCheckOpenFile(path, "r");

    char *line = NULL;
    size_t nlen = 0;
    if(grt_getline(&line, &nlen, fp) <= 0){
        GRTRaiseError("read sCoulomb fault header of %s failed.", path);
    }

    // sCoulomb 表格有 11 个数值列；第一列为占位符，首行对应 # 标记
    char header[COULOMB_HEADER_MAX_TOKENS][COULOMB_HEADER_TOKEN_SIZE] = {{0}};
    int nheader = sscanf(line,
        "%31s %31s %31s %31s %31s %31s %31s %31s %31s %31s %31s %31s",
        header[0], header[1], header[2], header[3], header[4], header[5],
        header[6], header[7], header[8], header[9], header[10], header[11]);
    if(nheader < COULOMB_HEADER_MIN_TOKENS || strcmp(header[0], "#") != 0){
        GRTRaiseError("invalid sCoulomb fault header in %s: expected # plus 10 field labels for 11 data columns.", path);
    }

    // 仅第 7 个数据列的完整 token "rake" 表示 rake/net-slip 格式
    bool rake_format = strcmp(header[6], "rake") == 0;
    size_t path_length = strlen(path);
    bool suffix_rake = path_length >= 4 && strcmp(path + path_length - 4, ".inr") == 0;
    if(suffix_rake && !rake_format){
        GRTRaiseWarning(
            "sCoulomb file \"%s\" has .inr suffix but the seventh header column is \"%s\", "
            "not the exact token \"rake\"; using the header-defined component format.",
            path, header[6]);
    } else if(!suffix_rake && rake_format){
        GRTRaiseWarning(
            "sCoulomb file \"%s\" has the exact token \"rake\" in the seventh header column "
            "but no .inr suffix; using the header-defined rake/net-slip format.", path);
    }

    // 第二行是与 11 个数据列对应的占位行
    if(grt_getline(&line, &nlen, fp) <= 0){
        GRTRaiseError("read sCoulomb fault placeholder header of %s failed.", path);
    }

    // 全局选项只解析一次，dt<0 仅保存原始延迟，dt>0 再生成 STF 并阻尼
    real_t *global_stf = NULL, global_delay = 0;
    int global_npts = 0;
    if(is_source && dt != 0 && stf_option) {
        global_stf = grt_time_function_from_option(stf_option, dt, &global_npts, &global_delay);
        if(dt > 0) {
            damp_stf(global_stf, global_npts, stf_decay);
        }
    }

    // 每行对应一个共享破裂过程的断层，全局 STF 优先于行末设置
    FINITE_FAULT *faults = NULL;
    size_t n = 0;
    size_t nempty = 0;
    size_t line_number = 2;
    while(grt_getline(&line, &nlen, fp) != -1){
        ++line_number;

        if(grt_is_empty_line(line)){
            ++nempty;
            continue;
        }
        if(grt_is_comment_line(line)) {
            continue;
        }

        real_t ID, kode_value;  // 第一列为占位符，约定填 1，不参与计算
        real_t east_begin, north_begin, east_end, north_end;
        real_t value1, value2, dip, top, bot;
        int nscan = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &ID, &east_begin, &north_begin, &east_end, &north_end,
            &kode_value, &value1, &value2, &dip, &top, &bot);
        if(nscan != 11){
            GRTRaiseError("parse sCoulomb fault data at line %zu of %s failed.", line_number, path);
        }

        if(is_source && kode_value != KODE_RTLAT_REVERSE && kode_value != KODE_RTLAT_TENSILE && kode_value != KODE_TENSILE_REVERSE &&
           kode_value != KODE_POINT_DC && kode_value != KODE_POINT_TENSILE_INFLATE){
            GRTRaiseError("invalid sCoulomb Kode at line %zu of %s.", line_number, path);
        }
        // 接收断层仅支持 Kode=100，滑动角由表头选择的格式确定
        if(!is_source && kode_value != KODE_RTLAT_REVERSE) {
            GRTRaiseError("Receiver faults require Kode=100 (line %zu of %s).", line_number, path);
        }

        faults = GRT_SAFE_REALLOC(faults, sizeof(FINITE_FAULT) * (n + 1));
        FINITE_FAULT *f = faults + n;
        memset(f, 0, sizeof(*f));

        // 保留文件中的几何和原始滑动量，实际源项及接收机制统一由 terms 管理
        f->east_begin = east_begin;
        f->north_begin = north_begin;
        f->east_end = east_end;
        f->north_end = north_end;
        f->kode = (unsigned int)kode_value;
        f->value1 = value1;
        f->value2 = value2;
        f->dip = dip;
        f->top = top;
        f->bot = bot;

        char where[256];
        snprintf(where, sizeof(where), "in %s line %zu", path, line_number);
        set_fault_components(f, rake_format, where);
        check_fault_geometry(f, where);
        set_fault_strike(f);
        set_fault_terms(f, is_source, rake_format, where);

        // dt=0 忽略 STF，其他模式统一解析延迟，全局选项优先于行末设置
        if(is_source && dt != 0) {
            char *tail = NULL;
            if(!stf_option) {
                grt_trim_whitespace(line);
                tail = line + strlen(line);
                while(tail > line && !isspace((unsigned char)tail[-1])) {
                    --tail;
                }
            }
            bool row_time = tail && strncmp(tail, "-D", 2) == 0;
            // 只有动态合成要求所有记录统一提供或省略 STF，仅提取延迟时允许混用
            if(dt > 0 && !stf_option && n > 0 && faults[0].stf_explicit != row_time) {
                GRTRaiseError("sCoulomb time functions must be present on every row or none (%s:%zu).", path, line_number);
            }

            f->stf_explicit = stf_option != NULL || row_time;
            if(stf_option) {
                f->stf_npts = global_npts;
                f->stf_delay = global_delay;
                if(dt > 0) {
                    f->stfd = GRT_SAFE_MALLOC(global_npts * sizeof(*f->stfd));
                    memcpy(f->stfd, global_stf, global_npts * sizeof(*f->stfd));
                }
            } else {
                // 仅动态震源解析自定义波形的相对路径，dt<0 不读取波形文件
                char *option = NULL;
                const char *slash = strrchr(path, '/');
                if(dt > 0 && row_time && strncmp(tail, "-D0/", 4) == 0 && tail[4] != '/' &&
                   !(isalpha((unsigned char)tail[4]) && tail[5] == ':') && slash) {
                    GRT_SAFE_ASPRINTF(&option, "-D0/%.*s/%s", (int)(slash - path), path, tail + 4);
                }
                f->stfd = grt_time_function_from_option(option ? option : (row_time ? tail : NULL), dt, &f->stf_npts, &f->stf_delay);
                GRT_SAFE_FREE_PTR(option);
                if(dt > 0) {
                    damp_stf(f->stfd, f->stf_npts, stf_decay);
                }
            }
        }

        n++;
    }

    GRT_SAFE_FREE_PTR(global_stf);
    GRT_SAFE_FREE_PTR(line);
    fclose(fp);

    if(nempty > 0){
        GRTRaiseWarning("skip %zu empty lines.", nempty);
    }

    if(n == 0){
        GRTRaiseError("no fault in %s.", path);
    }

    *nfault = n;
    return faults;
}


FINITE_FAULT *grt_finite_fault_from_option(const char *option, size_t *nfault, bool *has_i, real_t *dL, real_t *dW,
                                           bool is_source, real_t dt, real_t stf_decay, const char *stf_option)
{
    if((option == NULL) || (nfault == NULL) || (has_i == NULL) || (dL == NULL) || (dW == NULL)){
        GRTRaiseError("finite fault option is incomplete.");
    }

    char *option_copy = strdup(option);

    // 将文件路径和可选的剖分尺寸拆分到独立字符串中
    char *path = strtok(option_copy, "+");
    char *token = strtok(NULL, "+");
    if((path == NULL) || (*path == '\0') || (strtok(NULL, "+") != NULL)){
        GRTRaiseError("Error in finite fault option. expected <fault>[+i<dL>/<dW>]. Use \"-h\" for help.");
    }

    *dL = 0.0;
    *dW = 0.0;
    *has_i = token != NULL;
    if(*has_i){
        char extra;
        if((token[0] != 'i') || (sscanf(token + 1, "%lf/%lf%c", dL, dW, &extra) != 2)){
            GRTRaiseError("Error in finite fault option. expected +i<dL>/<dW>. Use \"-h\" for help.");
        }
        if(!((*dL == 0.0 && *dW == 0.0) || (*dL > 0.0 && *dW > 0.0))){
            GRTRaiseError("Error in finite fault option. dL and dW must both be zero or both be positive. Use \"-h\" for help.");
        }
    }

    FINITE_FAULT *faults = grt_finite_fault_load_coulomb(path, nfault, is_source, dt, stf_decay, stf_option);
    // 所有模块的 -U 共用几何剖分，未指定尺寸时只取整条断层的中心
    if(!is_source) {
        for(size_t i = 0; i < *nfault; ++i) {
            grt_finite_fault_subdiv(&faults[i], *dL, *dW, 0, NULL);
        }
    }
    GRT_SAFE_FREE_PTR(option_copy);
    return faults;
}


void grt_finite_fault_free(size_t nfault, FINITE_FAULT *faults)
{
    if(faults) {
        for(size_t i = 0; i < nfault; ++i) {
            GRT_SAFE_FREE_PTR(faults[i].east);
            GRT_SAFE_FREE_PTR(faults[i].north);
            GRT_SAFE_FREE_PTR(faults[i].depth);
            GRT_SAFE_FREE_PTR(faults[i].width);
            GRT_SAFE_FREE_PTR(faults[i].length);
            GRT_SAFE_FREE_PTR(faults[i].scalar_moments);
            GRT_SAFE_FREE_PTR(faults[i].stfd);
        }
    }
    GRT_SAFE_FREE_PTR(faults);
}


/**
 * 按子断层实际中心深度查询源介质，并保存完整等效张量的标量矩
 *
 * @param[in,out]  fault   已填充当前子断层几何的断层
 * @param[in]      isub    当前子断层索引
 * @param[in]      nlayer  原始模型层数
 * @param[in]      modarr  原始模型矩阵
 */
static void set_subfault_moment(FINITE_FAULT *fault, size_t isub, size_t nlayer, const real_t (*modarr)[GRT_MODARR_NCOL])
{
    real_t vp, vs, rho;
    grt_modarr_medium_at_depth(nlayer, modarr, fault->depth[isub], &vp, &vs, &rho);
    real_t area = KODE_IS_FINITE(fault->kode) ? fault->width[isub] * fault->length[isub] : 1;
    real_t tensor[6] = {0};

    // 混合源先叠加为完整张量，再求标量矩，不能把各源型的标量矩直接相加
    for(int t = 0; t < fault->nterms; ++t) {
        const FINITE_SOURCE_TERM *term = &fault->terms[t];
        real_t part[6];
        real_t scale = term->scale * area * (term->with_mu ? vs * vs * rho * 1e10 : 1);
        grt_source_moment_tensor(term->type, scale, vs > 0 ? vp / vs : 0, term->mechanism, part);
        for(int c = 0; c < 6; ++c) {
            tensor[c] += part[c];
        }
    }
    fault->scalar_moments[isub] = sqrt(0.5 * (tensor[0] * tensor[0] + tensor[3] * tensor[3] + tensor[5] * tensor[5]) +
                                       tensor[1] * tensor[1] + tensor[2] * tensor[2] + tensor[4] * tensor[4]);
}


void grt_finite_fault_subdiv(FINITE_FAULT *fault, real_t dL, real_t dW, size_t nlayer, const real_t (*modarr)[GRT_MODARR_NCOL])
{
    // 普通点源保持单点，其他 sCoulomb 源按断层面尺寸确定剖分数量
    real_t W = 0, L = 0;
    fault->nW = fault->nL = 1;
    if(fault->kode != KODE_POINT_GENERIC) {
        W = (fault->bot - fault->top) / sin(DEG1 * fault->dip);
        L = hypot(fault->east_end - fault->east_begin, fault->north_end - fault->north_begin);
        if(dL <= 0.0) {
            dL = L;
            dW = W;
        } else {
            fault->nW = GRT_MAX(1, (size_t)ceil(W / dW));
            fault->nL = GRT_MAX(1, (size_t)ceil(L / dL));
        }
    }

    // 一次分配全部子源几何，只为显式指定 STF 的动态有限源保存标量矩
    size_t nsub = fault->nW * fault->nL;
    fault->east           = GRT_SAFE_CALLOC(nsub, sizeof(*fault->east));
    fault->north          = GRT_SAFE_CALLOC(nsub, sizeof(*fault->north));
    fault->depth          = GRT_SAFE_CALLOC(nsub, sizeof(*fault->depth));
    fault->width          = GRT_SAFE_CALLOC(nsub, sizeof(*fault->width));
    fault->length         = GRT_SAFE_CALLOC(nsub, sizeof(*fault->length));
    fault->scalar_moments = modarr && fault->stf_explicit && fault->kode != KODE_POINT_GENERIC ?
                           GRT_SAFE_CALLOC(nsub, sizeof(*fault->scalar_moments)) : NULL;

    // 普通点源直接使用给定坐标，不需要断层面方向和尺寸
    if(fault->kode == KODE_POINT_GENERIC) {
        fault->east[0]  = fault->east_begin;
        fault->north[0] = fault->north_begin;
        fault->depth[0] = fault->top;
        if(fault->depth[0] < 0.0) {
            GRTRaiseError("Source point depth must be nonnegative (%.6g km).", fault->depth[0]);
        }
        return;
    }

    // 同一断层的方向余弦只计算一次，末块使用剩余尺寸及其实际中心
    real_t sind = sin(DEG1 * fault->dip);
    real_t cosd = cos(DEG1 * fault->dip);
    real_t sins = sin(DEG1 * fault->strike);
    real_t coss = cos(DEG1 * fault->strike);
    for(size_t iW = 0; iW < fault->nW; ++iW) {
        real_t width = GRT_MIN(dW, W - iW * dW);
        real_t w = iW * dW + 0.5 * width;
        for(size_t iL = 0; iL < fault->nL; ++iL) {
            real_t length = GRT_MIN(dL, L - iL * dL);
            real_t l = iL * dL + 0.5 * length;
            if(width <= 0.0 || length <= 0.0) {
                GRTRaiseError("Nonpositive subfault size at (iW=%zu, iL=%zu).", iW, iL);
            }
            size_t isub = iW * fault->nL + iL;
            fault->width[isub]  = width;
            fault->length[isub] = length;
            fault->depth[isub]  = fault->top + w * sind;
            if(fault->depth[isub] < 0.0) {
                GRTRaiseError("Negative subfault point depth %.6g km at (iW=%zu, iL=%zu).", fault->depth[isub], iW, iL);
            }
            fault->east[isub]   = fault->east_begin + l * sins + w * cosd * coss;
            fault->north[isub]  = fault->north_begin + l * coss - w * cosd * sins;

            // 按子源实际中心处的介质计算标量矩，用于保存有限源总矩率
            if(fault->scalar_moments) {
                set_subfault_moment(fault, isub, nlayer, modarr);
            }
        }
    }
}

FINITE_FAULT *grt_finite_fault_from_point(real_t depth, GRT_SYN_TYPE type, real_t scale, bool with_mu,
                                        const real_t *mechanism, real_t dt, real_t stf_decay, const char *stf_option)
{
    // 将普通点源封装为单个源项，以复用有限源剖分、合成和保存流程
    FINITE_FAULT *fault = GRT_SAFE_CALLOC(1, sizeof(*fault));
    fault->kode = KODE_POINT_GENERIC;
    fault->top = fault->bot = depth;
    fault->nterms = 1;
    fault->terms[0] = (FINITE_SOURCE_TERM){.type = type, .scale = scale, .with_mu = with_mu};
    memcpy(fault->terms[0].mechanism, mechanism, sizeof(fault->terms[0].mechanism));
    normalize_shear_term(&fault->terms[0], "point source");
    for(int im = 0; im < GRT_SRC_M_NUM; ++im) {
        fault->required_gf_sources[im] = grt_source_has_component(type, im);
    }
    // dt<0 仅提取延迟，dt>0 生成并阻尼时间函数，未指定时使用脉冲
    if(dt != 0) {
        fault->stf_explicit = stf_option != NULL;
        fault->stfd = grt_time_function_from_option(stf_option, dt, &fault->stf_npts, &fault->stf_delay);
        if(dt > 0) {
            damp_stf(fault->stfd, fault->stf_npts, stf_decay);
        }
    }
    return fault;
}
