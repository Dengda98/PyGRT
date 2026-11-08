/**
 * @file   eigenfn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-07
 * 
 * 使用 广义反射透射系数矩阵 将某层的垂直波函数传播到整个模型的物理层上下两侧，
 * 再计算任意深度的本征函数
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grt/dynamic/layer.h"

#include "grt/common/matrix.h"
#include "grt/common/model.h"
#include "grt/common/const.h"

#include "grt/modal/secular.h"
#include "grt/modal/eigenfn.h"
#include "grt/modal/energy.h"
#include "grt/modal/modal_util.h"
#include "grt/modal/root.h"

/**
 * 已知 z_{iy-1}+ 的垂直波函数，计算 z_{iy}+ 和 z_{iy}- 的垂直波函数 P-SV
 * 
 * @param[in]      mstat                模型结构体指针
 * @param[in]      iy                   待求解层位
 * @param[in]      M_RL                 广义矩阵 RD_RL
 * @param[in]      potRaylLove1         z_{iy-1}+ 的垂直波函数
 * @param[out]     potRaylLove2         z_{iy}+ 的垂直波函数
 * @param[out]     potRaylLove2_Up      z_{iy}- 的垂直波函数
 * 
 */
static void grt_potential_propagate_UpDown_Rayl(
    const MODEL1D_STATE *mstat, const size_t iy, const RT_MATRIX *M_RL,
    const cplx_t potRaylLove1[GRT_RAYL_DIM], cplx_t potRaylLove2[GRT_RAYL_DIM], cplx_t potRaylLove2_Up[GRT_RAYL_DIM])
{
    RT_MATRIX *M = &(RT_MATRIX){};
    grt_reset_RT_matrix_PSV(M);

    grt_RT_matrix_PSV(mstat, iy, M);

    real_t thk0 = mstat->mod1d->Thk[iy-1];
    cplx_t xa0 = mstat->xa[iy-1];
    cplx_t xb0 = mstat->xb[iy-1];

    cplx_t exa, exb;
    cplx_t E[2][2] = GRT_INIT_ZERO_2x2_MATRIX; 
    
    real_t eigenK = creal(mstat->k);
    exa = exp(- eigenK*thk0*xa0);
    exb = exp(- eigenK*thk0*xb0);
    E[0][0] = exa;
    E[1][1] = exb;
    
    // 先延时得到 \phi+_{j-}
    grt_cmat2x1_mul(E, potRaylLove1 + 2, potRaylLove2_Up + 2);

    // \phi+_{j+}
    cplx_t tmp[2][2];
    grt_cmat2x2_mul(M->RU, M_RL->RD, tmp);
    grt_cmat2x2_one_sub(tmp);
    grt_cmat2x2_inv(tmp, tmp);
    grt_cmat2x2_mul(tmp, M->TD, tmp);
    grt_cmat2x1_mul(tmp,   potRaylLove2_Up + 2, potRaylLove2 + 2);
    // \phi-_{j+}
    grt_cmat2x1_mul(M_RL->RD, potRaylLove2 + 2, potRaylLove2);

    // \phi-_{j-}
    cplx_t tmp1[2];
    grt_cmat2x1_mul(M->RD, potRaylLove2_Up + 2, tmp1);
    potRaylLove2_Up[0] = tmp1[0];
    potRaylLove2_Up[1] = tmp1[1];
    grt_cmat2x1_mul(M->TU, potRaylLove2, tmp1);
    potRaylLove2_Up[0] += tmp1[0];
    potRaylLove2_Up[1] += tmp1[1];

}


/** 已知 z_{iy-1}+ 的垂直波函数，计算 z_{iy}+ 和 z_{iy}- 的垂直波函数 SH */
static void grt_potential_propagate_UpDown_Love(
    const MODEL1D_STATE *mstat, const size_t iy, const RT_MATRIX *M_RL,
    const cplx_t potRaylLove1[GRT_LOVE_DIM], cplx_t potRaylLove2[GRT_LOVE_DIM], cplx_t potRaylLove2_Up[GRT_LOVE_DIM])
{
    RT_MATRIX *M = &(RT_MATRIX){};
    grt_reset_RT_matrix_SH(M);

    grt_RT_matrix_SH(mstat, iy, M);

    real_t eigenK = creal(mstat->k);
    real_t thk0 = mstat->mod1d->Thk[iy-1];
    cplx_t xb0 = mstat->xb[iy-1];
    bool *isLiquid = mstat->mod1d->isLiquid;

    cplx_t exb = exp(- eigenK*thk0*xb0);
        
    // 先延时得到 \phi+_{j-}
    // 如果上层为液体
    if(isLiquid[iy-1]){
        potRaylLove2_Up[0] = potRaylLove2_Up[1] = 0.0;
    } else {
        potRaylLove2_Up[1] = potRaylLove1[1] * exb;
    }

    // \phi+_{j+}, \phi-_{j+}, 
    // 如果下层为液体
    if(isLiquid[iy]){
        potRaylLove2[0] = potRaylLove2[1] = 0.0;
    } else {
        potRaylLove2[1] = M->TDL*potRaylLove2_Up[1] / (1.0 - M->RUL*M_RL->RDL);
        potRaylLove2[0] = M_RL->RDL * potRaylLove2[1];
    }

    // \phi-_{j-}
    // 如果上层为液体
    if(isLiquid[iy-1]){
        potRaylLove2_Up[0] = potRaylLove2_Up[1] = 0.0;
    } else {
        potRaylLove2_Up[0] = M->RDL * potRaylLove2_Up[1] + M->TUL * potRaylLove2[0];
    }
    
}

/**
 * 已知 z_{iy+1}+ 的垂直波函数，计算 z_{iy}+ 和 z_{iy+1}- 的垂直波函数 P-SV
 * 
 * @param[in]      mstat                模型结构体指针
 * @param[in]      iy                   待求解层位
 * @param[in]      M_FR                 广义矩阵 RU_FR
 * @param[out]     potRaylLove1         z_{iy}+ 的垂直波函数
 * @param[in]      potRaylLove2         z_{iy+1}+ 的垂直波函数
 * @param[out]     potRaylLove2_Up      z_{iy+1}- 的垂直波函数
 * 
 */
static void grt_potential_propagate_DownUp_Rayl(
    const MODEL1D_STATE *mstat, const size_t iy, const RT_MATRIX *M_FR,
    cplx_t potRaylLove1[GRT_RAYL_DIM], const cplx_t potRaylLove2[GRT_RAYL_DIM], cplx_t potRaylLove2_Up[GRT_RAYL_DIM])
{
    RT_MATRIX *M = &(RT_MATRIX){};
    grt_reset_RT_matrix_PSV(M);

    grt_RT_matrix_PSV(mstat, iy+1, M);
    
    real_t thk0 = mstat->mod1d->Thk[iy];
    cplx_t xa0 = mstat->xa[iy];
    cplx_t xb0 = mstat->xb[iy];

    cplx_t exa, exb;
    cplx_t E[2][2] = GRT_INIT_ZERO_2x2_MATRIX; 
    
    real_t eigenK = creal(mstat->k);
    exa = exp(- eigenK*thk0*xa0);
    exb = exp(- eigenK*thk0*xb0);
    E[0][0] = exa;
    E[1][1] = exb;

    cplx_t RU_FR_delay[2][2] = {0};
    memcpy(RU_FR_delay, M_FR->RU, sizeof(cplx_t)*4);
    grt_cmat2x2_mul(RU_FR_delay, E, RU_FR_delay);
    grt_cmat2x2_mul(E, RU_FR_delay, RU_FR_delay);

    // \phi-_{j+1}-
    // 这里计算 (I - RD^(j) * RU_FR_delay) ^ (-1) 不适用于液体-固体界面，
    // 因为从公式推导可得到，液体-固体界面对应的 |I - RD^(j) * RU_FR_delay| == 0， 其逆不存在
    if(potRaylLove2_Up[0] == 0.0){
        cplx_t tmp[2][2];
        grt_cmat2x2_mul(M->RD, RU_FR_delay, tmp);
        grt_cmat2x2_one_sub(tmp);
        grt_cmat2x2_inv(tmp, tmp);
        grt_cmat2x2_mul(tmp, M->TU, tmp);
        grt_cmat2x1_mul(tmp, potRaylLove2, potRaylLove2_Up);
        // \phi+_{j+1}-
        grt_cmat2x1_mul(RU_FR_delay, potRaylLove2_Up, potRaylLove2_Up + 2);
    }

    // \phi-_{j}+
    grt_cmat2x1_mul(E, potRaylLove2_Up, potRaylLove1);
    // \phi+_{j}+
    grt_cmat2x1_mul(M_FR->RU, potRaylLove1, potRaylLove1 + 2);

}

/** 已知 z_{iy+1}+ 的垂直波函数，计算 z_{iy}+ 和 z_{iy+1}- 的垂直波函数 SH */
static void grt_potential_propagate_DownUp_Love(
    const MODEL1D_STATE *mstat, const size_t iy, const RT_MATRIX *M_FR,
    cplx_t potRaylLove1[GRT_LOVE_DIM], const cplx_t potRaylLove2[GRT_LOVE_DIM], cplx_t potRaylLove2_Up[GRT_LOVE_DIM])
{
    // 如果上层为液体层
    if(mstat->mod1d->isLiquid[iy]){
        potRaylLove1[0]    = potRaylLove1[1]    = 0.0;
        potRaylLove2_Up[0] = potRaylLove2_Up[1] = 0.0;
        return;  // 直接返回
    }

    RT_MATRIX *M = &(RT_MATRIX){};
    grt_reset_RT_matrix_SH(M);

    grt_RT_matrix_SH(mstat, iy+1, M);

    real_t eigenK = creal(mstat->k);
    real_t thk0 = mstat->mod1d->Thk[iy];
    cplx_t xb0 = mstat->xb[iy];

    cplx_t exb = exp(- eigenK*thk0*xb0);
    cplx_t ex2b = exb * exb;
    cplx_t RUL_FR_delay = M_FR->RUL * ex2b;
    // \phi-_{j+1}-
    potRaylLove2_Up[0] = M->TUL * potRaylLove2[0] / (1.0 - M->RDL*RUL_FR_delay);
    // \phi+_{j+1}-
    potRaylLove2_Up[1] = RUL_FR_delay * potRaylLove2_Up[0];
    
    // \phi-_{j}+
    potRaylLove1[0] = potRaylLove2_Up[0] * exb;
    // \phi+_{j}+
    potRaylLove1[1] = potRaylLove1[0] * M_FR->RUL;
}


void grt_get_mod_potential_Up_Down_Rayl(
    MODEL1D_STATE *mstat, const size_t iref, const cplx_t potRaylLove[GRT_RAYL_DIM], const cplx_t potRaylLoveUp[GRT_RAYL_DIM], 
    cplx_t (*mod_potRaylLove_Down)[GRT_RAYL_DIM], cplx_t (*mod_potRaylLove_Up)[GRT_RAYL_DIM])
{
    size_t nlay = mstat->mod1d->n;

    // 先置零
    memset(mod_potRaylLove_Down, 0, sizeof(cplx_t)*GRT_RAYL_DIM*nlay);
    memset(mod_potRaylLove_Up,   0, sizeof(cplx_t)*GRT_RAYL_DIM*nlay);

    // 赋值 iref 层的垂直波函数
    memcpy(&mod_potRaylLove_Down[iref][0], potRaylLove, sizeof(cplx_t)*GRT_RAYL_DIM);
    memcpy(&mod_potRaylLove_Up[iref][0], potRaylLoveUp, sizeof(cplx_t)*GRT_RAYL_DIM);

    // 计算 RD_RL, RU_FR 矩阵
    RT_MATRIX *Mall_FR = GRT_SAFE_CALLOC(nlay, sizeof(RT_MATRIX));
    RT_MATRIX *Mall_RL = GRT_SAFE_CALLOC(nlay, sizeof(RT_MATRIX));
    grt_GRT_matrix_allLayer_Rayl(mstat, creal(mstat->k), Mall_RL, Mall_FR);

    // 向下传播
    for(size_t iy = iref+1; iy < nlay; ++iy){
        grt_potential_propagate_UpDown_Rayl(
            mstat, iy, &Mall_RL[iy], 
            mod_potRaylLove_Down[iy-1], mod_potRaylLove_Down[iy], mod_potRaylLove_Up[iy]);
    }

    // 向上传播
    for(size_t iy = iref; iy-- >0;){
        grt_potential_propagate_DownUp_Rayl(
            mstat, iy, &Mall_FR[iy], 
            mod_potRaylLove_Down[iy], mod_potRaylLove_Down[iy+1], mod_potRaylLove_Up[iy+1]);
    }

    GRT_SAFE_FREE_PTR(Mall_FR);
    GRT_SAFE_FREE_PTR(Mall_RL);
}


void grt_get_mod_potential_Up_Down_Love(
    MODEL1D_STATE *mstat, const size_t iref, const cplx_t potRaylLove[GRT_LOVE_DIM],
    cplx_t (*mod_potRaylLove_Down)[GRT_LOVE_DIM], cplx_t (*mod_potRaylLove_Up)[GRT_LOVE_DIM])
{
    size_t nlay = mstat->mod1d->n;

    // 先置零
    memset(mod_potRaylLove_Down, 0, sizeof(cplx_t)*GRT_LOVE_DIM*nlay);
    memset(mod_potRaylLove_Up,   0, sizeof(cplx_t)*GRT_LOVE_DIM*nlay);

    // 赋值 iref 层的垂直波函数
    memcpy(&mod_potRaylLove_Down[iref][0], potRaylLove, sizeof(cplx_t)*GRT_LOVE_DIM);

    // 计算 RDL_RL, RUL_FR 矩阵
    RT_MATRIX *Mall_FR = GRT_SAFE_CALLOC(nlay, sizeof(RT_MATRIX));
    RT_MATRIX *Mall_RL = GRT_SAFE_CALLOC(nlay, sizeof(RT_MATRIX));
    grt_GRT_matrix_allLayer_Love(mstat, creal(mstat->k), Mall_RL, Mall_FR);

    // 向下传播
    for(size_t iy = iref + 1; iy < nlay; ++iy){
        grt_potential_propagate_UpDown_Love(
            mstat, iy, &Mall_RL[iy], 
            mod_potRaylLove_Down[iy-1], mod_potRaylLove_Down[iy], mod_potRaylLove_Up[iy]);
    }

    // 向上传播
    for(size_t iy = iref; iy-- > 0;){
        grt_potential_propagate_DownUp_Love(
            mstat, iy, &Mall_FR[iy], 
            mod_potRaylLove_Down[iy], mod_potRaylLove_Down[iy+1], mod_potRaylLove_Up[iy+1]);
    }

    GRT_SAFE_FREE_PTR(Mall_FR);
    GRT_SAFE_FREE_PTR(Mall_RL);
}



void grt_get_mod_potential_Up_Down(
    MODEL1D_STATE *mstat, const DISPER_TYPE wtype, const size_t ncols, const size_t iref, 
    const cplx_t potRaylLove[ncols], const cplx_t potRaylLoveUp[ncols], 
    cplx_t (*mod_potRaylLove_Down)[ncols], cplx_t (*mod_potRaylLove_Up)[ncols])
{
    if(wtype == GRT_DISPERSION_RAYL && ncols == GRT_RAYL_DIM){
        grt_get_mod_potential_Up_Down_Rayl(mstat, iref, potRaylLove, potRaylLoveUp, mod_potRaylLove_Down, mod_potRaylLove_Up);
    } 
    else if(wtype == GRT_DISPERSION_LOVE && ncols == GRT_LOVE_DIM) {
        grt_get_mod_potential_Up_Down_Love(mstat, iref, potRaylLove, mod_potRaylLove_Down, mod_potRaylLove_Up);
    }
    else {
        GRTRaiseError("Wrong execution.");
    }
}


void grt_get_eigenfn_single_depth_Rayl(
    const MODEL1D_STATE *mstat, const cplx_t (*mod_potRaylLove_Down)[GRT_RAYL_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_RAYL_DIM],
    cplx_t T0[GRT_RAYL_DIM][GRT_RAYL_DIM], const bool reuseT, const real_t zsamp, const size_t ziref, cplx_t eigenfn[4])
{
    real_t dz = zsamp - mstat->mod1d->Dep[ziref];
    real_t thk = mstat->mod1d->Thk[ziref];
    cplx_t xa=0.0, xb=0.0;
    memset(eigenfn, 0, sizeof(cplx_t)*4);

    real_t eigenK = creal(mstat->k);
    xa = mstat->xa[ziref];
    xb = mstat->xb[ziref];

    if(! reuseT){
        grt_get_layer_D(mstat, ziref, false, 1, T0);
    }
    
    // 计算该点的垂直波函数，再得到本征函数
    // 复制 z_{j+1}- 侧的上传波场
    memcpy(eigenfn, &mod_potRaylLove_Up[ziref+1][0], sizeof(cplx_t)*2);
    // 复制 z_{j}+ 侧的下传波场
    memcpy(eigenfn + 2, &mod_potRaylLove_Down[ziref][0] + 2, sizeof(cplx_t)*2);

    cplx_t E[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    E[0][0] = exp(- eigenK * xa * dz);
    E[1][1] = exp(- eigenK * xb * dz);

    // 向下传播 \phi+
    grt_cmat2x1_mul(E, eigenfn + 2, eigenfn + 2);
    // 向上传播 \phi-
    E[0][0] = exp(- eigenK * xa * (thk - dz));
    E[1][1] = exp(- eigenK * xb * (thk - dz));
    grt_cmat2x1_mul(E, eigenfn, eigenfn);

    // 本征函数
    grt_cmatmx1_mul(GRT_RAYL_DIM, GRT_RAYL_DIM, T0, eigenfn, eigenfn);  

}


void grt_get_eigenfn_single_depth_Love(
    const MODEL1D_STATE *mstat, const cplx_t (*mod_potRaylLove_Down)[GRT_LOVE_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_LOVE_DIM],
    cplx_t T0[GRT_LOVE_DIM][GRT_LOVE_DIM], const bool reuseT, const real_t zsamp, const size_t ziref, cplx_t eigenfn[4])
{
    real_t dz = zsamp - mstat->mod1d->Dep[ziref];
    real_t thk = mstat->mod1d->Thk[ziref];
    cplx_t xb=0.0;
    memset(eigenfn, 0, sizeof(cplx_t)*4);

    real_t eigenK = creal(mstat->k);
    xb = mstat->xb[ziref];

    // 直接跳过液体层
    if(mstat->mod1d->isLiquid[ziref])  return;

    if(! reuseT){
        grt_get_layer_T(mstat, ziref, false, T0);
    }
    
    // 计算该点的垂直波函数，再得到本征函数
    // 复制 z_{j+1}- 侧的上传波场
    eigenfn[0] = mod_potRaylLove_Up[ziref+1][0];
    // 复制 z_{j}+ 侧的下传波场
    eigenfn[1] = mod_potRaylLove_Down[ziref][1];

    cplx_t E = exp(- eigenK * xb * dz);

    // 向下传播 \phi+
    eigenfn[1] *= E;
    // 向上传播 \phi-
    E = exp(- eigenK * xb * (thk - dz));
    eigenfn[0] *= E;

    // 本征函数
    grt_cmat2x1_mul(T0, eigenfn, eigenfn);  

}


void grt_get_eigenfn_depths_Rayl(
    const MODEL1D_STATE *mstat, const cplx_t (*mod_potRaylLove_Down)[GRT_RAYL_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_RAYL_DIM],
    const EIGENFN_INFO *eigfnmet, EIGENFN *eigfn)
{
    // 逐个循环处理采样点
    size_t last_ziref = -1;
    cplx_t T0[GRT_RAYL_DIM][GRT_RAYL_DIM];

    for(size_t iz = 0; iz < eigfnmet->nz; ++iz){
        real_t zsamp = eigfnmet->zs[iz];

        // 判断层位
        size_t ziref = eigfnmet->z_irefs[iz];
        
        // 利用了 zsamp 有序的性质计算本征函数
        grt_get_eigenfn_single_depth_Rayl(
            mstat, mod_potRaylLove_Down, mod_potRaylLove_Up, 
            T0, (last_ziref==ziref), zsamp, ziref, eigfn->fn[iz]);
        
        last_ziref = ziref;
    }
}


void grt_get_eigenfn_depths_Love(
    const MODEL1D_STATE *mstat, const cplx_t (*mod_potRaylLove_Down)[GRT_LOVE_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_LOVE_DIM],
    const EIGENFN_INFO *eigfnmet, EIGENFN *eigfn)
{
    // 逐个循环处理采样点
    size_t last_ziref = -1;
    cplx_t T0[GRT_LOVE_DIM][GRT_LOVE_DIM];

    for(size_t iz = 0; iz < eigfnmet->nz; ++iz){
        real_t zsamp = eigfnmet->zs[iz];

        // 判断层位
        size_t ziref = eigfnmet->z_irefs[iz];
        
        // 利用了 zsamp 有序的性质计算本征函数
        grt_get_eigenfn_single_depth_Love(
            mstat, mod_potRaylLove_Down, mod_potRaylLove_Up, 
            T0, (last_ziref==ziref), zsamp, ziref, eigfn->fn[iz]);
        
        last_ziref = ziref;
    }
}


void grt_get_eigenfn_depths(
    const MODEL1D_STATE *mstat, const DISPER_TYPE wtype, const size_t ncols, 
    const cplx_t (*mod_potRaylLove_Down)[ncols], const cplx_t (*mod_potRaylLove_Up)[ncols],
    const EIGENFN_INFO *eigfnmet, EIGENFN *eigfn)
{
    if(wtype == GRT_DISPERSION_RAYL && ncols == GRT_RAYL_DIM){
        grt_get_eigenfn_depths_Rayl(mstat, mod_potRaylLove_Down, mod_potRaylLove_Up, eigfnmet, eigfn);
    } 
    else if(wtype == GRT_DISPERSION_LOVE && ncols == GRT_LOVE_DIM) {
        grt_get_eigenfn_depths_Love(mstat, mod_potRaylLove_Down, mod_potRaylLove_Up, eigfnmet, eigfn);
    }
    else {
        GRTRaiseError("Wrong execution.");
    }
}


/**
 * 计算一个频率下的本征函数、能量积分和群速度，工作缓冲区由当前线程复用
 *
 * @param[in,out]  state                当前线程的模型状态
 * @param[in,out]  result               本征函数结果
 * @param[in]      iw                   频率索引
 * @param[in]      ncols                垂直波函数维数
 * @param[in]      energy               是否计算能量积分和群速度
 * @param[in,out]  potential_down       各物理层下行波工作缓冲区
 * @param[in,out]  potential_up         各物理层上行波工作缓冲区
 */
static void compute_eigenfn_frequency(MODEL1D_STATE *state, EIGENFN_INFO *result, size_t iw, int ncols, bool energy,
    cplx_t (*potential_down)[ncols], cplx_t (*potential_up)[ncols])
{
    real_t omega = PI2 * result->freqs[iw];
    EIGENV *roots = &result->eigv[iw];
    grt_update_mod1d_state_omega(state, omega, true);
    for(size_t ic = 0; ic < roots->n; ++ic){
        real_t phase = roots->c_roots[ic];
        real_t k = omega / phase;
        state->c_phase = phase;
        EIGENFN *mode = &result->eigfn[iw][ic];
        cplx_t potential[ncols], potential_top[ncols];
        cplx_t secular = 0.0;
        memset(potential, 0, sizeof(potential));
        memset(potential_top, 0, sizeof(potential_top));

        // 从搜根时的参考层传播到全部物理界面
        grt_secular_function_potential(state, phase, roots->c_roots_iref[ic], result->wtype, &secular, potential, potential_top);
        grt_get_mod_potential_Up_Down(state, result->wtype, ncols, roots->c_roots_iref[ic],
            potential, potential_top, potential_down, potential_up);
        if(mode->fn != NULL){
            grt_get_eigenfn_depths(state, result->wtype, ncols, potential_down, potential_up, result, mode);
        }
        if(energy){
            grt_energy_integrals(state, result->wtype, ncols, potential_down, potential_up, result, mode);
            cplx_t numerator = mode->egyint[1];
            if(result->wtype == GRT_DISPERSION_RAYL){
                numerator -= mode->egyint[2] / (2.0 * k);
            }
            roots->u_roots[ic] = creal(numerator / (phase * mode->egyint[0]));
        }
    }
}


/**
 * 偏移频率后重新计算相速度核，通过差分计算群速度核
 *
 * @param[in]      mod1d                原分层模型，副本用于微移频率计算
 * @param[in,out]  result               原频率的相速度核和群速度，写入群速度核
 * @param[in]      kernel_dz            与原结果一致的敏感核分段间距
 * @param[in]      shift                正频率偏移量，Hz
 */
static void compute_shifted_group_sensitivity(const MODEL1D *mod1d, EIGENFN_INFO *result, real_t kernel_dz, real_t shift)
{
    real_t *freqs = GRT_SAFE_CALLOC(result->nf, sizeof(real_t));
    for(size_t iw = 0; iw < result->nf; ++iw){
        freqs[iw] = result->freqs[iw] + shift;
    }
    size_t nmode = result->nmode > 0 ? result->modes[result->nmode - 1] + 1 : 0;
    EIGENV_INFO *roots = grt_init_eigenv_info(result->nf, freqs, result->wtype, nmode);
    GRT_SAFE_FREE_PTR(freqs);

    // 搜根会更新模型的参考频率，使用副本避免影响调用方或外层并行任务
    MODEL1D *shifted_model = grt_copy_mod1d(mod1d);
    grt_compute_eigenv(shifted_model, roots, false);
    EIGENFN_INFO *shifted = GRT_SAFE_CALLOC(1, sizeof(*shifted));
    grt_filter_eigenfn_info(0, NULL, false, result->nmode, result->modes, roots, shifted);
    EIGENFN_OPTIONS options = {.phase_kernel = true, .kernel_dz = kernel_dz};
    grt_compute_eigenfn(shifted_model, shifted, &options);
    grt_group_sensitivity_difference(result, shifted);

    grt_free_eigenfn_info(shifted);
    grt_free_eigenv_info(roots);
    grt_free_mod1d(shifted_model);
}


void grt_compute_eigenfn(MODEL1D *mod1d, EIGENFN_INFO *result, const EIGENFN_OPTIONS *options)
{
    if(options->kernel_dz < 0.0 || (options->nz > 0 && !options->zs)){
        GRTRaiseError("Invalid eigenfunction depth sampling.");
    }
    if(options->group_kernel_shift < 0.0){
        GRTRaiseError("Group-kernel frequency shift must be nonnegative.");
    }
    bool phase_kernel = options->phase_kernel || options->group_kernel;
    bool energy = options->energy || phase_kernel;
    if(!options->eigenfunctions && !energy){
        return;
    }
    const size_t NLAY = mod1d->n;

    // 本征函数默认在物理层界面采样，也可直接传入任意升序深度
    if(options->eigenfunctions){
        result->nz = options->nz > 0 ? options->nz : NLAY;
        result->zs = GRT_SAFE_CALLOC(result->nz, sizeof(real_t));
        result->z_irefs = GRT_SAFE_CALLOC(result->nz, sizeof(size_t));
        memcpy(result->zs, options->nz > 0 ? options->zs : mod1d->Dep, sizeof(real_t) * result->nz);
        size_t iref = 0;
        for(size_t iz = 0; iz < result->nz; ++iz){
            if(result->zs[iz] < 0.0 || (iz > 0 && result->zs[iz] <= result->zs[iz - 1])){
                GRTRaiseError("Eigenfunction depths must be nonnegative and strictly ascending.");
            }
            while(iref + 1 < NLAY && mod1d->Dep[iref + 1] <= result->zs[iz]){
                ++iref;
            }
            result->z_irefs[iz] = iref;
        }
    }

    // 能量积分的默认分段为原物理层，指定 dz 时仅细分有限层
    if(energy){
        if(options->kernel_dz == 0.0){
            result->cpar_nz = NLAY;
            result->cpar_zs = GRT_SAFE_CALLOC(NLAY, sizeof(real_t));
            result->cpar_z_irefs = GRT_SAFE_CALLOC(NLAY, sizeof(size_t));
            memcpy(result->cpar_zs, mod1d->Dep, sizeof(real_t) * NLAY);
            for(size_t iz = 0; iz < NLAY; ++iz){
                result->cpar_z_irefs[iz] = iz;
            }
        }
        else {
            for(size_t iy = 0; iy < NLAY; ++iy){
                real_t thk = mod1d->Thk[iy];
                real_t dz = options->kernel_dz;
                for(real_t h = 0.0; h < thk; h += dz){
                    if(fabs(thk - h) < dz * 1e-5){
                        break;
                    }
                    size_t iz = result->cpar_nz++;
                    result->cpar_zs = GRT_SAFE_REALLOC(result->cpar_zs, sizeof(real_t) * result->cpar_nz);
                    result->cpar_z_irefs = GRT_SAFE_REALLOC(result->cpar_z_irefs, sizeof(size_t) * result->cpar_nz);
                    result->cpar_zs[iz] = mod1d->Dep[iy] + h;
                    result->cpar_z_irefs[iz] = iy;
                    if(iy == NLAY - 1){
                        break;
                    }
                }
            }
        }
    }

    result->eigfn = GRT_SAFE_CALLOC(result->nf, sizeof(EIGENFN *));
    for(size_t iw = 0; iw < result->nf; ++iw){
        EIGENV *roots = &result->eigv[iw];
        result->eigfn[iw] = GRT_SAFE_CALLOC(roots->n, sizeof(EIGENFN));
        for(size_t ic = 0; ic < roots->n; ++ic){
            EIGENFN *mode = &result->eigfn[iw][ic];
            mode->eigenC = roots->c_roots[ic];
            if(options->eigenfunctions){
                mode->fn = GRT_SAFE_CALLOC(result->nz, sizeof(*mode->fn));
            }
            if(phase_kernel){
                mode->csens = GRT_SAFE_CALLOC(result->cpar_nz, sizeof(*mode->csens));
            }
            if(options->group_kernel){
                mode->usens = GRT_SAFE_CALLOC(result->cpar_nz, sizeof(*mode->usens));
            }
        }
    }

    // 外层已并行时，本次调用仅使用当前线程；否则按频率并行并复用每线程缓冲区
    int ncols = result->wtype == GRT_DISPERSION_RAYL ? GRT_RAYL_DIM : GRT_LOVE_DIM;
    #pragma omp parallel if(!omp_in_parallel())
    {
        MODEL1D_STATE *state = grt_init_mod1d_state(mod1d);
        cplx_t (*potential_down)[ncols] = GRT_SAFE_CALLOC(NLAY, sizeof(*potential_down));
        cplx_t (*potential_up)[ncols] = GRT_SAFE_CALLOC(NLAY + 1, sizeof(*potential_up));
        #pragma omp for schedule(guided)
        for(size_t iw = 0; iw < result->nf; ++iw){
            compute_eigenfn_frequency(state, result, iw, ncols, energy, potential_down, potential_up);
        }
        GRT_SAFE_FREE_PTR(potential_down);
        GRT_SAFE_FREE_PTR(potential_up);
        grt_free_mod1d_state(state);
    }
    if(options->group_kernel){
        if(options->group_kernel_shift > 0.0){
            compute_shifted_group_sensitivity(mod1d, result, options->kernel_dz, options->group_kernel_shift);
        }
        else {
            grt_group_sensitivity(result);
        }
    }
}
