/**
 * @file   static_util.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-06
 * 
 * 一些关于静态解的辅助函数
 * 
 */

#include <string.h>

#include "grt/static/static_util.h"
#include "grt/integral/kernel.h"

/** 衰减到 0 */
static real_t _decay_zero(MODEL1D *mod1d, real_t kmax_ref, size_t *Ncount)
{
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

/** 衰减到常数 */
static real_t _decay_constant(MODEL1D *mod1d, real_t kmax_ref, size_t *Ncount)
{
    MODEL1D_STATE *mstat = grt_init_mod1d_state(mod1d);
    grt_update_mod1d_state_omega(mstat, 1.0, true);

    // 利用振幅来估计一个合适的积分上限
    real_t dF = 0.0, F = 0.0, kmax = 0.0;
    cplxChnlGrid QWV = {0}, QWVz = {0};
    cplxChnlGrid QWV_const = {0}, QWVz_const = {0};
    kmax = kmax_ref * 1e-3;
    size_t ncount = 0;
    while(kmax < kmax_ref){
        dF = F = 0.0;
        grt_static_kernel(mstat, kmax, QWV, true, QWVz);

        // 对于收敛到常数的讨论，有必要细化到不同震源核函数的收敛形态
        GRT_LOOP_ChnlGrid(im, c){
            // 将 QWV, QWVz 乘上系数
            if(im == GRT_SRC_M_VF_INDEX || im == GRT_SRC_M_HF_INDEX){
                QWV[im][c] *= kmax;
            } else {
                QWVz[im][c] /= kmax;
            }

            F += fabs(QWV[im][c]) + fabs(QWVz[im][c]);
            dF += fabs(QWV[im][c] - QWV_const[im][c]) + fabs(QWVz[im][c] - QWVz_const[im][c]);
        }
        memcpy(QWV_const, QWV, sizeof(cplxChnlGrid));
        memcpy(QWVz_const, QWVz, sizeof(cplxChnlGrid));

        if(dF < 1e-3 * F) break;

        kmax *= 1.2;
        ncount++;
    }
    grt_free_mod1d_state(mstat);

    kmax = GRT_MIN(kmax, kmax_ref);

    if(Ncount != NULL)  *Ncount = ncount;

    return kmax;
}

real_t grt_predict_static_kmax(MODEL1D *mod1d, real_t kmax_ref, size_t *Ncount)
{
    if(mod1d->depsrc == mod1d->deprcv){
        return _decay_constant(mod1d, kmax_ref, Ncount);
    } else {
        return _decay_zero(mod1d, kmax_ref, Ncount);
    }

}