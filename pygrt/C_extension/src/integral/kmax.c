/**
 * @file   kmax.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-06
 * 
 * 波数积分辅助函数
 * 
 */

#include <string.h>

#include "grt/integral/kmax.h"
#include "grt/integral/kernel.h"

/** 衰减到 0 */
static real_t _decay_zero(
    MODEL1D_STATE *mstat, GRT_KernelFunc kerfunc,
    real_t kmax_init, real_t kmax_ref, size_t *Ncount)
{
    // 利用振幅来估计一个合适的积分上限
    real_t Fmax = 0.0, F = 0.0, kmax = 0.0;
    cplxChnlGrid QWV = {0}, QWVz = {0};
    kmax = kmax_init;
    size_t ncount = 0;
    size_t nval = 0;
    while(kmax < kmax_ref){
        F = 0.0;
        kerfunc(mstat, kmax, QWV, true, QWVz);
        GRT_LOOP_ChnlGrid(im, c){
            F += sqrt(kmax) * (fabs(QWV[im][c]) + fabs(QWVz[im][c]));
        }
        Fmax = GRT_MAX(F, Fmax);

        if(F < 1e-3 * Fmax){
            nval++;
        } else {
            nval = 0;
        }

        if(nval >= 3) break;

        kmax *= 1.2;
        ncount++;
    }
    kmax = GRT_MIN(kmax, kmax_ref);

    if(Ncount != NULL)  *Ncount = ncount;

    return kmax;
}

/** 衰减到常数 */
static real_t _decay_constant(
    MODEL1D_STATE *mstat, GRT_KernelFunc kerfunc,
    real_t kmax_init, real_t kmax_ref, size_t *Ncount)
{
    // 利用振幅来估计一个合适的积分上限
    real_t dF = 0.0, F = 0.0, kmax = 0.0;
    cplxChnlGrid QWV = {0}, QWVz = {0};
    cplxChnlGrid QWV_const = {0}, QWVz_const = {0};
    kmax = kmax_init;
    size_t ncount = 0;
    size_t nval = 0;
    while(kmax < kmax_ref){
        dF = F = 0.0;
        kerfunc(mstat, kmax, QWV, true, QWVz);

        // 对于收敛到常数的讨论，有必要细化到不同震源核函数的收敛形态
        GRT_LOOP_ChnlGrid(im, c){
            // 将 QWV, QWVz 乘上系数
            if(GRT_SRC_M_INDEX_IS_FORCE(im)){
                QWV[im][c] *= kmax;
            } else {
                QWVz[im][c] /= kmax;
            }

            F += fabs(QWV[im][c]) + fabs(QWVz[im][c]);
            dF += fabs(QWV[im][c] - QWV_const[im][c]) + fabs(QWVz[im][c] - QWVz_const[im][c]);
        }
        memcpy(QWV_const, QWV, sizeof(cplxChnlGrid));
        memcpy(QWVz_const, QWVz, sizeof(cplxChnlGrid));

        if(dF < 1e-3 * F){
            nval++;
        } else {
            nval = 0;
        }

        if(nval >= 3)  break;

        kmax *= 1.2;
        ncount++;

        // printf("ncount = %zu, F = %.3e, dF = %.3e, kmax = %.3e, kref = %.3e\n", ncount, F, dF, kmax, kmax_ref);
    }
    kmax = GRT_MIN(kmax, kmax_ref);

    if(Ncount != NULL)  *Ncount = ncount;

    return kmax;
}

real_t grt_predict_kmax(
    MODEL1D_STATE *mstat, GRT_KernelFunc kerfunc,
    real_t kmax_init, real_t kmax_ref, size_t *Ncount)
{
    if(mstat->mod1d->depsrc == mstat->mod1d->deprcv){
        return _decay_constant(mstat, kerfunc, kmax_init, kmax_ref, Ncount);
    } else {
        return _decay_zero(mstat, kerfunc, kmax_init, kmax_ref, Ncount);
    }
}
