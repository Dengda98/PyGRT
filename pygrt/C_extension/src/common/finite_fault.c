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
    f->rake = 1.0 / DEG1 * atan2(f->reverse, - f->right_lateral);
    f->slip = hypot(f->right_lateral, f->reverse); // m
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
    while(grt_getline(&line, &nlen, fp) != -1){
        faults = (FINITE_FAULT *)realloc(faults, sizeof(FINITE_FAULT) * (n + 1));
        FINITE_FAULT *f = faults + n;
        memset(f, 0, sizeof(*f));

        real_t dum1, dum2;
        int nscan = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &dum1, &f->east_begin, &f->north_begin, &f->east_end, &f->north_end,
            &dum2, &f->right_lateral, &f->reverse, &f->dip, &f->top, &f->bot);
        if(nscan != 11){
            fclose(fp);
            GRT_SAFE_FREE_PTR(line);
            GRT_SAFE_FREE_PTR(faults);
            GRTRaiseError("parse line %zu of %s failed.", n + 3, path);
        }

        char where[256];
        snprintf(where, sizeof(where), "in %s line %zu", path, n + 3);
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
