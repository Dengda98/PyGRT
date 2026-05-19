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
    
    real_t eigenK = mstat->k;
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

    real_t eigenK = mstat->k;
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
    
    real_t eigenK = mstat->k;
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

    real_t eigenK = mstat->k;
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
    RT_MATRIX *Mall_FR = (RT_MATRIX *)calloc(nlay, sizeof(RT_MATRIX));
    RT_MATRIX *Mall_RL = (RT_MATRIX *)calloc(nlay, sizeof(RT_MATRIX));
    grt_GRT_matrix_allLayer_Rayl(mstat, mstat->k, Mall_RL, Mall_FR);

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
    RT_MATRIX *Mall_FR = (RT_MATRIX *)calloc(nlay, sizeof(RT_MATRIX));
    RT_MATRIX *Mall_RL = (RT_MATRIX *)calloc(nlay, sizeof(RT_MATRIX));
    grt_GRT_matrix_allLayer_Love(mstat, mstat->k, Mall_RL, Mall_FR);

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

    real_t eigenK = mstat->k;
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

    real_t eigenK = mstat->k;
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

