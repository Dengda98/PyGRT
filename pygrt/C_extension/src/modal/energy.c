/**
 * @file   energy.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-07
 * 
 * 根据每层界面垂直波函数，计算能量积分和群速度
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
#include "grt/modal/energy.h"

/** 
 * 求解特定形式的积分, 用于求解 Rayleigh 能量积分的理论解
 * \int_0^H x^T * diag[E1(z), E2(z)] * M4x4 * diag[E1(z), E2(z)] * y dz，
 * 其中 E1 = diag[e^{-a(H-h)}, e^{-b(H-h)}], E2 = diag[e^{-ah}, e^{-bh}], 
 * 
 * @param[in]     H            层厚
 * @param[in]     xa           P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]     xb           S波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]     k            本征波数
 * @param[in]     M            提前构造好的 M 矩阵
 * @param[in]     potRaylLove  垂直波函数,用于构造 x,y 矢量
 * 
 * @return     积分理论解
 * 
 */
static cplx_t _solve_typical_integral_Rayl(
    const real_t H, const cplx_t xa, const cplx_t xb, 
    const real_t k, const cplx_t M[GRT_RAYL_DIM][GRT_RAYL_DIM], const cplx_t potRaylLove[GRT_RAYL_DIM])
{
    // 将维度降至四项 2x2 ，然后计算各项积分的解析解
    // 辅助变量
    cplx_t c11, c12, c21, c22, x1, x2, y1, y2;

    // 厚度 H 为负认作正无穷，由于exa=exb=0.0，res3和res4结果不受影响（均为0）
    cplx_t exa = (H>=0.0) ? exp(- k*H*xa) : 0.0;
    cplx_t exb = (H>=0.0) ? exp(- k*H*xb) : 0.0;
    cplx_t ex2a = exa*exa;
    cplx_t ex2b = exb*exb;
    cplx_t exab = exa*exb;

    cplx_t a = k*xa;
    cplx_t b = k*xb;

    // 1, phi(-)^T * E1 * M11 * E1 * phi(-)
    c11 = M[0][0];    c12 = M[0][1];    
    c21 = M[1][0];    c22 = M[1][1];   
    x1 = (potRaylLove[0]);
    x2 = (potRaylLove[1]);
    y1 = potRaylLove[0];
    y2 = potRaylLove[1];
    cplx_t res1 = 
           0.5*c11*x1*y1/a * (1.0 - ex2a) + 0.5*c22*x2*y2/b * (1.0 - ex2b)
         + (c21*x2*y1 + c12*x1*y2)/(a + b) * (1.0 - exab);

    // 2, phi(+)^T * E2 * M22 * E2 * phi(+), 与 1 表达式一致
    c11 = M[2][2];    c12 = M[2][3];    
    c21 = M[3][2];    c22 = M[3][3];   
    x1 = (potRaylLove[2]);
    x2 = (potRaylLove[3]);
    y1 = potRaylLove[2];
    y2 = potRaylLove[3];
    cplx_t res2 = 
           0.5*c11*x1*y1/a * (1.0 - ex2a) + 0.5*c22*x2*y2/b * (1.0 - ex2b)
         + (c21*x2*y1 + c12*x1*y2)/(a + b) * (1.0 - exab);

    // 3, phi(-)^T * E1 * M12 * E2 * phi(+)
    c11 = M[0][2];    c12 = M[0][3];    
    c21 = M[1][2];    c22 = M[1][3];   
    x1 = (potRaylLove[0]);
    x2 = (potRaylLove[1]);
    y1 = potRaylLove[2];
    y2 = potRaylLove[3];
    cplx_t res3 = 
          H*(c11*x1*y1*exa + c22*x2*y2*exb) + 1.0/(a-b) * (exb - exa) * (c12*x1*y2 + c21*x2*y1);
    
    // 4, phi(+)^T * E2 * M21 * E1 * phi(-), 与 3 表达式一致
    c11 = M[2][0];    c12 = M[2][1];    
    c21 = M[3][0];    c22 = M[3][1];   
    x1 = (potRaylLove[2]);
    x2 = (potRaylLove[3]);
    y1 = potRaylLove[0];
    y2 = potRaylLove[1];
    cplx_t res4 = 
          H*(c11*x1*y1*exa + c22*x2*y2*exb) + 1.0/(a-b) * (exb - exa) * (c12*x1*y2 + c21*x2*y1);
    
    cplx_t res = res1 + res2 + res3 + res4;
    return res;
}


/** 
 * 求解特定形式的积分, 用于求解 Love 能量积分的理论解
 * \int_0^H x^T * E(z) * M2x2 * E(z) * x dz
 * 其中 E = diag[e^{-b(H-h)}, e^{-bh}]
 * 
 * @param[in]     H            层厚
 * @param[in]     xb           S波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]     k            本征波数
 * @param[in]     M            提前构造好的 M 矩阵
 * @param[in]     potRaylLove  垂直波函数,用于提取 x 矢量
 * 
 * @return     积分理论解
 * 
 */
static cplx_t _solve_typical_integral_Love(
    const real_t H, const cplx_t xb, 
    const real_t k, const cplx_t M[GRT_LOVE_DIM][GRT_LOVE_DIM], const cplx_t potRaylLove[GRT_LOVE_DIM])
{
    // 辅助变量
    cplx_t c11, c12, c21, c22, x1, x2, y1, y2;

    // 厚度 H 为负认作正无穷，由于exa=exb=0.0，res3和res4结果不受影响（均为0）
    cplx_t exb = (H>=0.0) ? exp(- k*H*xb) : 0.0;
    cplx_t ex2b = exb*exb;

    cplx_t b = k*xb;

    c11 = M[0][0];    c12 = M[0][1];    
    c21 = M[1][0];    c22 = M[1][1];   
    x1 = potRaylLove[0];
    x2 = potRaylLove[1];
    y1 = potRaylLove[0];
    y2 = potRaylLove[1];

    cplx_t res = H*(c12*x1*y2 + c21*x2*y1)*exb + 0.5*c11*x1*y1/b*(1.0 - ex2b) + 0.5*c22*x2*y2/b*(1.0 - ex2b);

    return res;
}



/** 
 * 给定某层的物性参数，计算单层的 Rayleigh 波能量积分的 6 个子项,
 * 其被积函数分别为
 *     1. (q)^2
 *     2. (w)^2
 *     3. q*(dw/dz)
 *     4. w*(dq/dz)
 *     5. (dq/dz)^2
 *     6. (dw/dz)^2
 * 
 * @param[in]     mod1d        模型结构体指针
 * @param[in]     iy           层位索引
 * @param[in]     H            层厚
 * @param[in]     potRaylLove  垂直波函数
 * @param[out]    sub_egyint   6 个积分子项结果
 * 
 */
static void energy_integrals_single_layer_Rayl(
    const GRT_MODEL1D *mod1d, const size_t iy, const real_t H, const cplx_t potRaylLove[GRT_RAYL_DIM], cplx_t sub_egyint[6])
{
    cplx_t M[GRT_RAYL_DIM][GRT_RAYL_DIM] = {0};
    cplx_t D11[2][2] = {0}, D12[2][2] = {0};
    cplx_t D11_uiz[2][2] = {0}, D12_uiz[2][2] = {0};
    cplx_t Diag[2][2] = {0};
    cplx_t D2x4[2][4] = {0}, D4x2[4][2] = {0};
    cplx_t D2x4_uiz[2][4] = {0}, D4x2_uiz[4][2] = {0};
    cplx_t tmp4x2[4][2] = {0};

    // 构造 M 矩阵
    grt_get_layer_D11(mod1d, iy, D11);
    grt_get_layer_D12(mod1d, iy, D12);
    grt_get_layer_D11_uiz(mod1d, iy, D11_uiz);
    grt_get_layer_D12_uiz(mod1d, iy, D12_uiz);

    grt_cmatmxn_block_assign(2, 4, D2x4, 0, 0, 2, 2, D11);
    grt_cmatmxn_block_assign(2, 4, D2x4, 0, 2, 2, 2, D12);
    grt_cmatmxn_transpose(2, 4, D2x4, D4x2);

    grt_cmatmxn_block_assign(2, 4, D2x4_uiz, 0, 0, 2, 2, D11_uiz);
    grt_cmatmxn_block_assign(2, 4, D2x4_uiz, 0, 2, 2, 2, D12_uiz);
    grt_cmatmxn_transpose(2, 4, D2x4_uiz, D4x2_uiz);
    
    real_t k  = mod1d->k;
    cplx_t xa = mod1d->xa[iy];
    cplx_t xb = mod1d->xb[iy];

    // 1. (q)^2
    memset(Diag, 0, sizeof(cplx_t)*4);
    Diag[0][0] = 1.0;
    grt_cmatmxn_mul(4, 2, 2, D4x2, Diag, tmp4x2);
    grt_cmatmxn_mul(4, 2, 4, tmp4x2, D2x4, M);
    sub_egyint[0] = _solve_typical_integral_Rayl(H, xa, xb, k, M, potRaylLove);

    // 2. (w)^2
    memset(Diag, 0, sizeof(cplx_t)*4);
    Diag[1][1] = 1.0;
    grt_cmatmxn_mul(4, 2, 2, D4x2, Diag, tmp4x2);
    grt_cmatmxn_mul(4, 2, 4, tmp4x2, D2x4, M);
    sub_egyint[1] = _solve_typical_integral_Rayl(H, xa, xb, k, M, potRaylLove);

    // 3. q*(dw/dz)
    memset(Diag, 0, sizeof(cplx_t)*4);
    Diag[0][1] = 1.0;
    grt_cmatmxn_mul(4, 2, 2, D4x2, Diag, tmp4x2);
    grt_cmatmxn_mul(4, 2, 4, tmp4x2, D2x4_uiz, M);
    sub_egyint[2] = _solve_typical_integral_Rayl(H, xa, xb, k, M, potRaylLove);

    // 4. w*(dq/dz)
    memset(Diag, 0, sizeof(cplx_t)*4);
    Diag[1][0] = 1.0;
    grt_cmatmxn_mul(4, 2, 2, D4x2, Diag, tmp4x2);
    grt_cmatmxn_mul(4, 2, 4, tmp4x2, D2x4_uiz, M);
    sub_egyint[3] = _solve_typical_integral_Rayl(H, xa, xb, k, M, potRaylLove);

    // 5. (dq/dz)^2
    memset(Diag, 0, sizeof(cplx_t)*4);
    Diag[0][0] = 1.0;
    grt_cmatmxn_mul(4, 2, 2, D4x2_uiz, Diag, tmp4x2);
    grt_cmatmxn_mul(4, 2, 4, tmp4x2, D2x4_uiz, M);
    sub_egyint[4] = _solve_typical_integral_Rayl(H, xa, xb, k, M, potRaylLove);

    // 6. (dw/dz)^2
    memset(Diag, 0, sizeof(cplx_t)*4);
    Diag[1][1] = 1.0;
    grt_cmatmxn_mul(4, 2, 2, D4x2_uiz, Diag, tmp4x2);
    grt_cmatmxn_mul(4, 2, 4, tmp4x2, D2x4_uiz, M);
    sub_egyint[5] = _solve_typical_integral_Rayl(H, xa, xb, k, M, potRaylLove);

}


/** 
 * 给定某层的物性参数，计算单层的 Love 波能量积分的 2 个子项,
 * 其被积函数分别为
 *     1. (v)^2
 *     2. (dv/dz)^2
 * 
 * @param[in]     mod1d        模型结构体指针
 * @param[in]     iy           层位索引
 * @param[in]     H            层厚
 * @param[in]     potRaylLove  垂直波函数
 * @param[out]    sub_egyint   2 个积分子项结果
 * 
 */
static void energy_integrals_single_layer_Love(
    const GRT_MODEL1D *mod1d, const size_t iy, const real_t H, const cplx_t potRaylLove[GRT_LOVE_DIM], cplx_t sub_egyint[2])
{
    cplx_t M[2][2] = {0};

    real_t k  = mod1d->k;
    cplx_t xb = mod1d->xb[iy];

    // 1. (v)^2
    memset(M, 0, sizeof(cplx_t)*4);
    M[0][0] = M[0][1] = M[1][0] = M[1][1] = 1.0;
    sub_egyint[0] = GRT_SQUARE(k) * _solve_typical_integral_Love(H, xb, k, M, potRaylLove);

    // 2. (dv/dz)^2
    memset(M, 0, sizeof(cplx_t)*4);
    M[0][0] = M[1][1] = 1;
    M[0][1] = M[1][0] = -1;
    sub_egyint[1] = GRT_SQUARE(k*k*xb) * _solve_typical_integral_Love(H, xb, k, M, potRaylLove);
}

/**
 * 计算某层的 Love 波相速度敏感核表达式中的分子项
 * 
 * @param[in]     rho           密度
 * @param[in]     omega         圆频率
 * @param[in]     k             本征波数
 * @param[in]     mu            拉梅常数 mu
 * @param[in]     sub_egyint    2 个积分子项结果
 * @param[out]    K             敏感核
 * 
 */
static void phase_sensitivity_numerator_single_layer_Love(
    const real_t rho, const real_t omega, const real_t k, const cplx_t mu, const cplx_t sub_egyint[2], cplx_t K[GRT_SNSTVTY_MAX])
{
    cplx_t K_mu, K_rho;
    
    K_mu = mu * (k*k * sub_egyint[0] + sub_egyint[1]);
    K_rho = - rho * GRT_SQUARE(omega) * sub_egyint[0];

    K[0] = 0.0;
    // K_beta
    K[1] = 2.0 * K_mu;
    // K_rho'
    K[2] = K_mu + K_rho;
}

/**
 * 计算某层的 Rayleigh 波相速度敏感核表达式中的分子项
 * 
 * @param[in]     rho           密度
 * @param[in]     omega         圆频率
 * @param[in]     k             本征波数
 * @param[in]     lambda        拉梅常数 lambda
 * @param[in]     mu            拉梅常数 mu
 * @param[in]     sub_egyint    6 个积分子项结果
 * @param[out]    K             敏感核
 * 
 */
static void phase_sensitivity_numerator_single_layer_Rayl(
    const real_t rho, const real_t omega, const real_t k, const cplx_t lambda, const cplx_t mu, const cplx_t sub_egyint[6], cplx_t K[GRT_SNSTVTY_MAX])
{
    cplx_t K_lam, K_mu, K_rho;
    real_t kk = k*k;

    K_lam = lambda * (kk*sub_egyint[0] - 2.0*k*sub_egyint[2] + sub_egyint[5]);
    K_mu = mu * ( kk*sub_egyint[1] + 2.0*k*sub_egyint[3] + sub_egyint[4] + 2.0*kk*sub_egyint[0] + 2.0*sub_egyint[5] );
    K_rho = - rho * GRT_SQUARE(omega) * (sub_egyint[0] + sub_egyint[1]);

    // // K_alpha
    K[0] = 2.0 * (lambda + 2.0*mu)/lambda * K_lam;
    // K_beta
    K[1] = (2.0*K_mu - 4.0*mu/lambda*K_lam);
    // K_rho'
    K[2] = K_lam + K_mu + K_rho;

    // K_alpha
    // K[0] = K_lam;
    // // K_beta
    // K[1] = K_mu;
    // // K_rho'
    // K[2] = K_rho;
}



void grt_energy_integrals_Rayl(
    const GRT_MODEL1D *mod1d, const cplx_t (*mod_potRaylLove_Down)[GRT_RAYL_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_RAYL_DIM], 
    const EIGENFN_INFO *eigfnmet, EIGENFN *eigfn)
{
    // 先置零
    memset(eigfn->egyint, 0, sizeof(cplx_t)*GRT_EGYINTS_MAX);

    cplx_t potRaylLove[GRT_RAYL_DIM], potRaylLove_bak[GRT_RAYL_DIM] = {0};
    cplx_t xa=0.0, xb=0.0, mu=0.0, lambda=0.0;
    real_t thk=0.0, rho=0.0;

    // 将每层的能量积分结果累加
    size_t ziref, last_ziref;
    ziref = last_ziref = -1;
    real_t dz=0.0, h=0.0;

    real_t eigenK = mod1d->k;
    size_t cpar_nz = eigfnmet->cpar_nz;

    // 如果没有定义cpar系列数组，则表示使用模型数组，
    // 此时要求 cpar_nz == mod1d->n
    if((eigfnmet->cpar_zs==NULL || eigfnmet->cpar_z_irefs==NULL) && eigfnmet->cpar_nz != mod1d->n){
        GRTRaiseError("cpar_xx variables error.\n");
    }

    for(size_t iz = 0; iz < cpar_nz; ++iz){
        ziref = (eigfnmet->cpar_z_irefs!=NULL)? eigfnmet->cpar_z_irefs[iz] : iz;

        // 应用深度有序的性质
        if(ziref != last_ziref){
            xa = mod1d->xa[ziref];
            xb = mod1d->xb[ziref];
            lambda = mod1d->lambda[ziref];
            mu = mod1d->mu[ziref];
            thk = mod1d->Thk[ziref];  
            rho = mod1d->Rho[ziref];

            // 构造垂直波函数
            memcpy(potRaylLove_bak, &mod_potRaylLove_Up[ziref+1][0], sizeof(cplx_t)*2);
            memcpy(potRaylLove_bak+2, &mod_potRaylLove_Down[ziref][0]+2, sizeof(cplx_t)*2);
        }

        if(eigfnmet->cpar_zs!=NULL){
            dz = (iz < cpar_nz-1)? eigfnmet->cpar_zs[iz+1] - eigfnmet->cpar_zs[iz] : -1.0;   // 这里负厚度用于表示正无穷
            h = eigfnmet->cpar_zs[iz] - mod1d->Dep[ziref];
        } else {
            // 模型物理层
            dz = (iz < cpar_nz-1)? mod1d->Thk[ziref] : -1.0;
            h = 0.0;
        }

        // 为薄层重新定义上行和下行垂直波函数
        cplx_t E[GRT_RAYL_DIM][GRT_RAYL_DIM] = {0};
        E[0][0] = (iz < cpar_nz-1)? exp( - eigenK * xa * (thk - h - dz)) : 0.0;
        E[1][1] = (iz < cpar_nz-1)? exp( - eigenK * xb * (thk - h - dz)) : 0.0;
        E[2][2] = exp( - eigenK * xa * h);
        E[3][3] = exp( - eigenK * xb * h);
        grt_cmatmx1_mul(GRT_RAYL_DIM, GRT_RAYL_DIM, E, potRaylLove_bak, potRaylLove);

        // 计算能量积分所需的子项, 这里负厚度用于表示正无穷
        cplx_t sub_egyint[6];
        energy_integrals_single_layer_Rayl(mod1d, ziref, dz, potRaylLove, sub_egyint);

        // I1
        eigfn->egyint[0] += 0.5 * rho * (sub_egyint[0] + sub_egyint[1]);
        // I2
        eigfn->egyint[1] += 0.5 * ((lambda+2.0*mu)*sub_egyint[0] + mu*sub_egyint[1]);
        // I3
        eigfn->egyint[2] += lambda*sub_egyint[2] - mu*sub_egyint[3];
        // I4
        eigfn->egyint[3] += 0.5 * (mu*sub_egyint[4] + (lambda+2.0*mu)*sub_egyint[5]);

        // 计算相速度敏感核
        if(eigfn->csens!=NULL){
            phase_sensitivity_numerator_single_layer_Rayl(rho, creal(mod1d->omega), eigenK, lambda, mu, sub_egyint, eigfn->csens[iz]);
            // 去掉系数
            eigfn->csens[iz][0] /= mod1d->Va[ziref]/mod1d->c_phase;
            eigfn->csens[iz][1] /= mod1d->Vb[ziref]/mod1d->c_phase;
            eigfn->csens[iz][2] /= rho/mod1d->c_phase;
        }

        last_ziref = ziref;
    }

    // 修正敏感核系数
    if(eigfn->csens!=NULL){
        cplx_t dom = 4.0*GRT_SQUARE(eigenK)*eigfn->egyint[1] - 2.0*eigenK*eigfn->egyint[2];
        if(dom == 0.0){
            GRTRaiseError("Wrong execution.");
        }
        for(size_t iz = 0; iz < cpar_nz; ++iz){
            for(int j = 0; j < GRT_SNSTVTY_MAX; ++j){
                eigfn->csens[iz][j] /= dom;
            }
        }
    }

}



void grt_energy_integrals_Love(
    const GRT_MODEL1D *mod1d, const cplx_t (*mod_potRaylLove_Down)[GRT_LOVE_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_LOVE_DIM], 
    const EIGENFN_INFO *eigfnmet, EIGENFN *eigfn)
{
    // 先置零
    memset(eigfn->egyint, 0, sizeof(cplx_t)*GRT_EGYINTS_MAX);

    cplx_t potRaylLove[2], potRaylLove_bak[2] = {0};
    cplx_t xb=0.0, mu=0.0;
    real_t thk=0.0, rho=0.0;
    
    // 将每层的能量积分结果累加
    size_t ziref, last_ziref;
    ziref = last_ziref = -1;
    real_t dz=0.0, h=0.0;

    real_t eigenK = mod1d->k;
    size_t cpar_nz = eigfnmet->cpar_nz;

    // 如果没有定义cpar系列数组，则表示使用模型数组，
    // 此时要求 cpar_nz == mod1d->n
    if((eigfnmet->cpar_zs==NULL || eigfnmet->cpar_z_irefs==NULL) && eigfnmet->cpar_nz != mod1d->n){
        GRTRaiseError("cpar_xx variables error.\n");
    }

    for(size_t iz=0; iz<cpar_nz; ++iz){
        ziref = (eigfnmet->cpar_z_irefs!=NULL)? eigfnmet->cpar_z_irefs[iz] : iz;

        // 应用深度有序的性质
        if(ziref != last_ziref){
            xb = mod1d->xb[ziref];
            mu = mod1d->mu[ziref];
            thk = mod1d->Thk[ziref];  
            rho = mod1d->Rho[ziref];

            // 构造垂直波函数
            potRaylLove_bak[0] = mod_potRaylLove_Up[ziref+1][0];
            potRaylLove_bak[1] = mod_potRaylLove_Down[ziref][1];
        }

        // 跳过液体层， Love 波无位移应力，不会为能量积分和敏感核有贡献
        if(mod1d->isLiquid[ziref])  continue;

        if(eigfnmet->cpar_zs!=NULL){
            dz = (iz < cpar_nz-1)? eigfnmet->cpar_zs[iz+1] - eigfnmet->cpar_zs[iz] : -1.0;   // 这里负厚度用于表示正无穷
            h = eigfnmet->cpar_zs[iz] - mod1d->Dep[ziref];
        } else {
            // 模型物理层
            dz = (iz < cpar_nz-1)? mod1d->Thk[ziref] : -1.0;
            h = 0.0;
        }
        

        // 为薄层重新定义上行和下行垂直波函数
        cplx_t E[2][2] = {0};
        E[0][0] = (iz < cpar_nz-1)? exp( - eigenK * xb * (thk - h - dz)) : 0.0;
        E[1][1] = exp( - eigenK * xb * h);
        grt_cmat2x1_mul(E, potRaylLove_bak, potRaylLove);

        // 计算能量积分所需的子项
        cplx_t sub_egyint[2] = {0};
        energy_integrals_single_layer_Love(mod1d, ziref, dz, potRaylLove, sub_egyint);

        // I1
        eigfn->egyint[0] += 0.5 * rho * sub_egyint[0];

        // I2
        eigfn->egyint[1] += 0.5 * mu * sub_egyint[0];

        // 位移和位移偏导的交叉项为0
        eigfn->egyint[2] += 0.0;

        // I3
        eigfn->egyint[3] += 0.5 * mu * sub_egyint[1];

        // 计算相速度敏感核
        if(eigfn->csens!=NULL){
            phase_sensitivity_numerator_single_layer_Love(rho, creal(mod1d->omega), eigenK, mu, sub_egyint, eigfn->csens[iz]);
            // 去掉系数
            eigfn->csens[iz][1] /= mod1d->Vb[ziref]/mod1d->c_phase;
            eigfn->csens[iz][2] /= rho/mod1d->c_phase;
        }

    }
    
    // 修正敏感核系数
    if(eigfn->csens!=NULL){
        cplx_t dom = 4.0*GRT_SQUARE(eigenK)*eigfn->egyint[1];
        if(dom == 0.0){
            GRTRaiseError("Wrong execution.");
        }
        for(size_t iz = 0; iz < cpar_nz; ++iz){
            for(int j = 0; j < GRT_SNSTVTY_MAX; ++j){
                eigfn->csens[iz][j] /= dom;
            }
        }
    }

}



void grt_energy_integrals(
    const GRT_MODEL1D *mod1d, const DISPER_TYPE wtype, const size_t ncols, 
    const cplx_t (*mod_potRaylLove_Down)[ncols], const cplx_t (*mod_potRaylLove_Up)[ncols], 
    const EIGENFN_INFO *eigfnmet, EIGENFN *eigfn)
{
    if(wtype == GRT_DISPERSION_RAYL && ncols == GRT_RAYL_DIM){
        grt_energy_integrals_Rayl(mod1d, mod_potRaylLove_Down, mod_potRaylLove_Up, eigfnmet, eigfn);
    }
    else if(wtype == GRT_DISPERSION_LOVE && ncols == GRT_LOVE_DIM) {
        grt_energy_integrals_Love(mod1d, mod_potRaylLove_Down, mod_potRaylLove_Up, eigfnmet, eigfn);
    }
    else {
        GRTRaiseError("Wrong execution.");
    }
}


