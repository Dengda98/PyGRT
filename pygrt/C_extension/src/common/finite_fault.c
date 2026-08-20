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

/** 检查 dip ∈ (0, 90] 且 bot > top */
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


void grt_finite_fault_set_derived(FINITE_FAULT *f)
{
    f->strike = 1.0 / DEG1 * atan2(f->east_end - f->east_begin, f->north_end - f->north_begin);
    if(f->kode == KODE_POINT_TENSILE_INFLATE){
        f->rake = 0.0;
        f->slip = 0.0;
    } else {
        f->rake = 1.0 / DEG1 * atan2(f->reverse, -f->right_lateral);
        f->slip = KODE_IS_POINT(f->kode) ? 0.0 : hypot(f->right_lateral, f->reverse);
    }
}

/** 按 Kode 解释文件第 7、8 列 */
static void set_fault_components(FINITE_FAULT *f, bool rake_format, const char *where)
{
    f->rake_format = rake_format;
    if(rake_format && f->kode != KODE_RTLAT_REVERSE){
        GRTRaiseError("%s: .inr rake/net slip format only supports Kode=100.", where);
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
    // 跳过两行表头（列名 + 占位行）
    for(int ih = 0; ih < 2; ++ih){
        if(grt_getline(&line, &nlen, fp) <= 0){
            fclose(fp);
            GRT_SAFE_FREE_PTR(line);
            GRTRaiseError("read header of %s failed.", path);
        }
    }

    FINITE_FAULT *faults = NULL;
    size_t n = 0;
    size_t line_number = 2;
    size_t path_length = strlen(path);
    bool rake_format = path_length >= 4 && strcmp(path + path_length - 4, ".inr") == 0;
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

        grt_finite_fault_set_derived(f);
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


void grt_finite_fault_free(FINITE_FAULT *faults)
{
    free(faults);
}


void grt_finite_fault_subdiv(
    const FINITE_FAULT *fault, real_t dL, real_t dW,
    real_t *W, real_t *L, size_t *nW, size_t *nL)
{
    if(dL <= 0.0 || dW <= 0.0){
        GRTRaiseError("dL and dW must be positive.");
    }
    check_fault_geometry(fault, "finite fault");

    *W = (fault->bot - fault->top) / sin(DEG1 * fault->dip);
    *L = hypot(fault->east_end - fault->east_begin, fault->north_end - fault->north_begin);
    if(*W <= 0.0 || *L <= 0.0){
        GRTRaiseError("fault along-dip/along-strike length must be positive (W=%.6g, L=%.6g).", *W, *L);
    }
    *nW = GRT_MAX(1, (size_t)ceil(*W / dW));
    *nL = GRT_MAX(1, (size_t)ceil(*L / dL));
}


void grt_finite_fault_subfault(
    const FINITE_FAULT *fault,
    real_t dL, real_t dW, real_t W, real_t L,
    size_t iW, size_t iL,
    FINITE_SUBFAULT *sub)
{
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
