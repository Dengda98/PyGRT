/**
 * @file   static_util.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-06
 * 
 * 一些关于静态解的辅助函数
 * 
 */

#include "grt/static/static_util.h"
#include "grt/integral/kernel.h"


real_t grt_predict_static_kmax(MODEL1D *mod1d, real_t kmax_ref, size_t *Ncount)
{
    if(mod1d->depsrc == mod1d->deprcv){
        return kmax_ref;
    }

    MODEL1D_STATE *mstat = grt_init_mod1d_state(mod1d);
    grt_update_mod1d_state_omega(mstat, 1.0, true);

    // 利用振幅来估计一个合适的积分上限
    real_t Fmax = 0.0, F = 0.0, kmax = 0.0;
    cplxChnlGrid QWV = {0}, QWVz = {0};
    kmax = kmax_ref * 1e-3;
    size_t ncount = 0;
    while(kmax < kmax_ref){
        F = 0.0;
        grt_static_kernel(mstat, kmax, QWV, true, QWVz);
        GRT_LOOP_ChnlGrid(im, c){
            F += sqrt(kmax) * (fabs(QWV[im][c]) + fabs(QWVz[im][c]));
        }
        Fmax = GRT_MAX(F, Fmax);

        if(F < 1e-3 * Fmax) break;

        kmax *= 1.2;
        ncount++;
    }
    grt_free_mod1d_state(mstat);

    kmax = GRT_MIN(kmax, kmax_ref);

    if(Ncount != NULL)  *Ncount = ncount;

    return kmax;
}