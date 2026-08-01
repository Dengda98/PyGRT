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

/**
 * kmax 搜索参数
 *
 * 步长采用 log(k) 上等间距（自适应几何因子），使静态（k_init << k_ref）
 * 与动态高频（k_init 接近 k_ref）都能在约 GRT_KMAX_SEARCH_N 步内扫完区间。
 * 同时限制绝对步长不超过 kmax_ref / GRT_KMAX_SEARCH_N，避免后期几何
 * 递进产生过大的 dk。
 */
#define GRT_KMAX_SEARCH_N          40      ///< 覆盖 [k_init, k_ref] 的目标步数
#define GRT_KMAX_FACTOR_MIN        1.01    ///< 几何因子下限，保证每步有推进
#define GRT_KMAX_FACTOR_MAX        1.50    ///< 几何因子上限，避免步长过大
#define GRT_KMAX_NVAL_NEED         3       ///< 连续满足收敛准则的次数
#define GRT_KMAX_MAX_STEP_DIVISOR  GRT_KMAX_SEARCH_N
#define GRT_KMAX_DECAY_ZERO_TOL    1e-2    ///< 衰减到 0：F < tol * Fmax
#define GRT_KMAX_DECAY_CONST_TOL   5e-2    ///< 衰减到常数：dF < tol * Fmax


/** 根据搜索区间计算自适应几何因子 */
static real_t _kmax_step_factor(real_t kmax_init, real_t kmax_ref)
{
    real_t ratio, factor;
    if(kmax_init <= 0.0 || kmax_ref <= kmax_init){
        return GRT_KMAX_FACTOR_MIN;
    }
    ratio = kmax_ref / kmax_init;
    factor = pow(ratio, 1.0 / (real_t)GRT_KMAX_SEARCH_N);
    factor = GRT_MIN(GRT_MAX(factor, GRT_KMAX_FACTOR_MIN), GRT_KMAX_FACTOR_MAX);
    return factor;
}

/** 将搜索前期的最低判定波数限制在 [kmax_init, kmax_ref]。 */
static real_t _kmax_low(real_t kmax_init, real_t kmax_low, real_t kmax_ref)
{
    return GRT_MIN(GRT_MAX(kmax_low, kmax_init), kmax_ref);
}


/** 衰减到 0 */
static real_t _decay_zero(
    MODEL1D_STATE *mstat, GRT_KernelFunc kerfunc,
    real_t kmax_init, real_t kmax_low, real_t kmax_ref, size_t *Ncount)
{
    real_t Fmax = 0.0, F = 0.0, kmax = 0.0;
    real_t factor = 0.0, dk_max = 0.0, dk = 0.0, k_low_eff = 0.0;
    cplxChnlGrid QWV = {0}, QWVz = {0};
    size_t ncount = 0;
    size_t nval = 0;

    if(kmax_init >= kmax_ref){
        if(Ncount != NULL) *Ncount = 0;
        return kmax_ref;
    }

    factor = _kmax_step_factor(kmax_init, kmax_ref);
    dk_max = kmax_ref / GRT_KMAX_MAX_STEP_DIVISOR;
    k_low_eff = _kmax_low(kmax_init, kmax_low, kmax_ref);
    kmax = kmax_init;

    while(kmax < kmax_ref){
        F = 0.0;
        kerfunc(mstat, kmax, QWV, true, QWVz);
        GRT_LOOP_ChnlGrid(im, c){
            F += sqrt(kmax) * (fabs(QWV[im][c]) + fabs(QWVz[im][c]));
        }
        Fmax = GRT_MAX(F, Fmax);

        if(kmax >= k_low_eff && F < GRT_KMAX_DECAY_ZERO_TOL * Fmax){
            nval++;
        } else {
            nval = 0;
        }

        if(nval >= GRT_KMAX_NVAL_NEED) break;

        dk = GRT_MIN((factor - 1.0) * kmax, dk_max);
        kmax += dk;
        ncount++;
    }
    kmax = GRT_MIN(kmax, kmax_ref);

    if(Ncount != NULL)  *Ncount = ncount;

    return kmax;
}

/** 衰减到常数 */
static real_t _decay_constant(
    MODEL1D_STATE *mstat, GRT_KernelFunc kerfunc,
    real_t kmax_init, real_t kmax_low, real_t kmax_ref, size_t *Ncount)
{
    real_t dF = 0.0, F = 0.0, Fmax = 0.0, kmax = 0.0;
    real_t factor = 0.0, dk_max = 0.0, dk = 0.0, k_low_eff = 0.0;
    cplxChnlGrid QWV = {0}, QWVz = {0};
    cplxChnlGrid QWV_const = {0}, QWVz_const = {0};
    size_t ncount = 0;
    size_t nval = 0;

    if(kmax_init >= kmax_ref){
        if(Ncount != NULL) *Ncount = 0;
        return kmax_ref;
    }

    factor = _kmax_step_factor(kmax_init, kmax_ref);
    dk_max = kmax_ref / GRT_KMAX_MAX_STEP_DIVISOR;
    k_low_eff = _kmax_low(kmax_init, kmax_low, kmax_ref);
    kmax = kmax_init;

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
        Fmax = GRT_MAX(F, Fmax);
        memcpy(QWV_const, QWV, sizeof(cplxChnlGrid));
        memcpy(QWVz_const, QWVz, sizeof(cplxChnlGrid));

        // dF 是当前核函数与前一采样点（收敛值估计）之间的差值。
        if(kmax >= k_low_eff && dF < GRT_KMAX_DECAY_CONST_TOL * (Fmax + 1e-30)){
            nval++;
        } else {
            nval = 0;
        }

        if(nval >= GRT_KMAX_NVAL_NEED)  break;

        dk = GRT_MIN((factor - 1.0) * kmax, dk_max);
        kmax += dk;
        ncount++;
    }
    kmax = GRT_MIN(kmax, kmax_ref);

    if(Ncount != NULL)  *Ncount = ncount;

    return kmax;
}

real_t grt_predict_kmax(
    MODEL1D_STATE *mstat, GRT_KernelFunc kerfunc,
    real_t kmax_init, real_t kmax_low, real_t kmax_ref, size_t *Ncount)
{
    if(mstat->mod1d->depsrc == mstat->mod1d->deprcv){
        return _decay_constant(mstat, kerfunc, kmax_init, kmax_low, kmax_ref, Ncount);
    } else {
        return _decay_zero(mstat, kerfunc, kmax_init, kmax_low, kmax_ref, Ncount);
    }
}
