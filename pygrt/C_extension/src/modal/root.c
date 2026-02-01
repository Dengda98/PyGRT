/**
 * @file   root.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-05
 * 
 *     总函数：寻找Rayleigh波和Love波的久期函数零点
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/common/search.h"
#include "grt/common/progressbar.h"

#include "grt/modal/root.h"
#include "grt/modal/saroot.h"

static size_t grt_get_approx_maxnroots(GRT_MODEL1D *mod1d, const real_t freq, const real_t cmax)
{
    real_t res = 0;
    grt_mod1d_xa_xb(mod1d, PI2*freq/cmax);

    for(size_t iy=0; iy < mod1d->n - 1; ++iy){
        res += fabs(cimag(mod1d->xa[iy] + ((mod1d->Vb[iy] > 0.0) ? mod1d->xb[iy] : 0.0))) * mod1d->Thk[iy];
        // printf("xa=" GRT_CMPLX_FMT ", xb=" GRT_CMPLX_FMT "\n", GRT_CMPLX_SPLIT(xa), GRT_CMPLX_SPLIT(xb));
    }
    res *= 2.0*freq;
    res /= cmax;
    return GRT_MAX((size_t)res, 1);
}


static void get_dc_min_max(GRT_MODEL1D *mod1d, const real_t freq, const real_t cmin, const real_t cmax, real_t *pt_dcmin, real_t *pt_dcmax)
{
    size_t approx_n = grt_get_approx_maxnroots(mod1d, freq, cmax);
    *pt_dcmax = (cmax - cmin) / approx_n;
    *pt_dcmin = *pt_dcmax * 1e-4;
    // printf("# dcmin=%e, dcmax=%e, n=%d, cmax=%f\n", *pt_dcmin, *pt_dcmax, approx_n, cmax);
}

/** 将对单个频率的处理包装成函数，方便被其它函数并行调用 */
static void get_secular_roots_single_freq(GRT_MODEL1D *mod1d, EIGENV_METHOD *eigmet, EIGENV *eigv)
{
    get_dc_min_max(mod1d, mod1d->omega/PI2, eigmet->cmin, eigmet->cmax, &eigv->dcmin, &eigv->dcmax);

    int secRaylType = 0; // 使用哪种Rayleigh波的secular function
    size_t iref = 0; // 使用哪一层的secular function

    // 对于首层
    iref = 0;
    if(eigmet->wtype == GRT_DISPERSION_RAYL){
        // Rayleigh 波
        // 搜索“基阶”进行补充
        secRaylType = 0;
        grt_sa_secular_roots(mod1d, eigmet, eigv, secRaylType, iref);
    } 
    else if(eigmet->wtype == GRT_DISPERSION_LOVE){
        // Love 波 跳过液体层
        if(! mod1d->isLiquid[iref]){
            secRaylType = 1;
            grt_sa_secular_roots(mod1d, eigmet, eigv, secRaylType, iref);
        }
    }
    else {
        GRTRaiseError("Wrong execution.");
    }

    // 检查低速层
    for(iref=1; iref<mod1d->n-1; ++iref){
        // 跳过两侧固体的递增 Vp 和 递增 Vs
        if( ! mod1d->isLiquid[iref-1] && ! mod1d->isLiquid[iref] && mod1d->Vb[iref-1] <= mod1d->Vb[iref] && mod1d->Va[iref-1] <= mod1d->Va[iref]){
            continue;
        }
        // 跳过两侧液体的递增 Vp
        if( mod1d->isLiquid[iref-1] && mod1d->isLiquid[iref] && mod1d->Va[iref-1] <= mod1d->Va[iref]){
            continue;
        }
        // 对于 Love 波跳过液体层
        if(eigmet->wtype == GRT_DISPERSION_LOVE && mod1d->isLiquid[iref]){
            continue;
        }
        secRaylType = 1;
        grt_sa_secular_roots(mod1d, eigmet, eigv, secRaylType, iref);
    }
}


void grt_get_secular_roots(GRT_MODEL1D *mod1d, EIGENV_METHOD *eigmet, const bool print_progressbar)
{
    real_t *freqs = eigmet->freqs;
    size_t nf = eigmet->nf;

    mod1d->omgref = PI2*freqs[nf-1];
    
    // 进度条变量 
    int progress=0;

#define __CALL_SECULAR_ROOT(mod1d, iw) \
    mod1d->omega = PI2*freqs[iw]; \
    get_secular_roots_single_freq(mod1d, eigmet, eigmet->eigv+iw);

    if(eigmet->print_sec){
        __CALL_SECULAR_ROOT(mod1d, 0);
        goto FINISH;
    }

#ifdef _OPENMP
    // 如果外部已在并行区域中，此处串行
    if (omp_in_parallel()) {
        for(size_t iw = 0; iw < nf; ++iw){
            __CALL_SECULAR_ROOT(mod1d, iw);
            // 记录进度条变量 
            progress++;
            if(print_progressbar) grt_printprogressBar("Computing Dispersion Curves: ", progress*100/nf);
        }
    }
    // 否则正常使用并行
    else{
        #pragma omp parallel for schedule(guided) default(shared) 
        for(size_t iw = 0; iw < nf; ++iw){
            GRT_MODEL1D *local_mod1d = grt_copy_mod1d(mod1d);
            __CALL_SECULAR_ROOT(local_mod1d, iw);

            // 记录进度条变量 
            #pragma omp critical
            {   
                progress++;
                if(print_progressbar) grt_printprogressBar("Computing Dispersion Curves: ", progress*100/nf);
            }
            grt_free_mod1d(local_mod1d);
        }
    }

#else
    // 串行
    for(size_t iw = 0; iw < nf; ++iw){
        __CALL_SECULAR_ROOT(mod1d, iw);
        // 记录进度条变量 
        progress++;
        if(print_progressbar) grt_printprogressBar("Computing Dispersion Curves: ", progress*100/nf);
    }
#endif

FINISH:
}