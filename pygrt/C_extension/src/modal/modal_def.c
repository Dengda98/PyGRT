/**
 * @file   modal_def.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-01
 * 
 *                   
 */

#include "grt/modal/modal_def.h"

void grt_free_eigenv(EIGENV *eigv)
{
    GRT_SAFE_FREE_PTR(eigv->c_roots);
    GRT_SAFE_FREE_PTR(eigv->c_roots_iref);
}

void grt_free_eigenv_method(EIGENV_METHOD *eigmet)
{
    GRT_SAFE_FREE_PTR(eigmet->freqs);
    for(size_t i = 0; i < eigmet->nf; ++i){
        grt_free_eigenv(eigmet->eigv + i);
    }
    GRT_SAFE_FREE_PTR(eigmet->eigv);
    GRT_SAFE_FREE_PTR(eigmet);
}


