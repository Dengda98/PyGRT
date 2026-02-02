/**
 * @file   modal_def.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-01
 * 
 *                   
 */

#include "grt/modal/modal_def.h"
#include "grt/common/search.h"

void grt_free_eigenv(EIGENV *eigv)
{
    GRT_SAFE_FREE_PTR(eigv->c_roots);
    GRT_SAFE_FREE_PTR(eigv->c_roots_iref);
    GRT_SAFE_FREE_PTR(eigv->g_roots);
}

void grt_free_eigenv_info(EIGENV_INFO *eigmet)
{
    GRT_SAFE_FREE_PTR(eigmet->freqs);
    for(size_t i = 0; i < eigmet->nf; ++i){
        grt_free_eigenv(eigmet->eigv + i);
    }
    GRT_SAFE_FREE_PTR(eigmet->eigv);
    GRT_SAFE_FREE_PTR(eigmet);
}


// 从全部频散信息中，筛选出部分信息
void grt_filter_eigenfn_info(
    const size_t nf, const real_t *freqs, const bool def_freq_range, 
    const size_t nmode, const size_t *modes, 
    EIGENV_INFO *eigmet, EIGENFN_INFO *eigfnmet)
{
    // 先申请全部内存，然后再缩
    eigfnmet->nf = 0;
    eigfnmet->freqs = (real_t *)calloc(eigmet->nf, sizeof(real_t));
    eigfnmet->nmode = 0;
    eigfnmet->modes = (size_t *)calloc(eigmet->nmode, sizeof(size_t));
    eigfnmet->eigv = (EIGENV *)calloc(eigmet->nf, sizeof(EIGENV));
    eigfnmet->wtype = eigmet->wtype;

    const real_t *in_freqs = eigmet->freqs;

    // 筛选阶次
    for(size_t i = 0; i < eigmet->nmode; ++i){
        bool match = false;
        if(modes == NULL){
            match = true;
        }
        else{
            match = (grt_findElement_size_t(modes, nmode, i) >= 0);
        }
        
        if(match){
            eigfnmet->modes[eigfnmet->nmode++] = i;
        }
    }
    eigfnmet->modes = (size_t *)realloc(eigfnmet->modes, sizeof(size_t)*eigfnmet->nmode);

    // 筛选频率
    for(size_t iw = 0; iw < eigmet->nf; ++iw){
        bool match = false;
        
        if(freqs == NULL){
            match = true;
        }
        else if(def_freq_range){
            match = (in_freqs[iw] >= freqs[0] && in_freqs[iw] <= freqs[1]);
        }
        else {
            size_t idx = grt_findClosest_real_t(freqs, nf, in_freqs[iw]);
            match = (fabs(freqs[idx] - in_freqs[iw]) < 1e-8);
        }
        
        if(! match) continue;

        eigfnmet->freqs[eigfnmet->nf] = in_freqs[iw];

        // 写入频散
        EIGENV *eigvtmp = &eigfnmet->eigv[eigfnmet->nf];
        eigvtmp->c_roots = (real_t *)calloc(eigmet->nmode, sizeof(real_t));
        eigvtmp->g_roots  = (real_t *)calloc(eigmet->nmode, sizeof(real_t));
        eigvtmp->c_roots_iref = (size_t *)calloc(eigmet->nmode, sizeof(size_t));
        eigvtmp->n = 0;

        EIGENV *in_eigv = &eigmet->eigv[iw];
        for(size_t i = 0; i < eigfnmet->nmode; ++i){
            size_t imode = eigfnmet->modes[i];
            if(imode >= in_eigv->n) break;

            eigvtmp->c_roots[eigvtmp->n] = in_eigv->c_roots[imode];
            eigvtmp->c_roots_iref[eigvtmp->n] = in_eigv->c_roots_iref[imode];

            eigvtmp->n++;
        }
        eigvtmp->c_roots = (real_t *)realloc(eigvtmp->c_roots, sizeof(real_t)*eigvtmp->n);
        eigvtmp->g_roots  = (real_t *)realloc(eigvtmp->g_roots, sizeof(real_t)*eigvtmp->n);
        eigvtmp->c_roots_iref = (size_t *)realloc(eigvtmp->c_roots_iref, sizeof(size_t)*eigvtmp->n);

        eigfnmet->nf++;
    }
    eigfnmet->freqs = (real_t *)realloc(eigfnmet->freqs, sizeof(real_t)*eigfnmet->nf);
}