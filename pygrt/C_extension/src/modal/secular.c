/**
 * @file   secular.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-07
 * 
 *     计算Rayleigh波和Love波的久期函数
 * 
 */

#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>

#include "grt/dynamic/layer.h"
// #include "grt/dynamic/propagate.h"

#include "grt/common/RT_matrix.h"
#include "grt/common/model.h"
#include "grt/common/matrix.h"
#include "grt/common/search.h"

#include "grt/modal/secular.h"


#define __DEBUG_SECULAR__  1


// void grt_get_mod_all_RD_RU_Rayl(
//     const GRT_MODEL1D *mod1d, const cplx_t omega, const real_t k,
//     cplx_t RD_RL0[mod1d->n][2][2], cplx_t RU_FR0[mod1d->n][2][2], int *stats)
// {
//     // 先置零
//     if(RD_RL0!=NULL) memset(RD_RL0, 0, sizeof(cplx_t)*mod1d->n*4); 
//     if(RU_FR0!=NULL) memset(RU_FR0, 0, sizeof(cplx_t)*mod1d->n*4); 

// #define X(F) cplx_t (*F)[2][2] = (cplx_t (*)[2][2])calloc(mod1d->n, sizeof(cplx_t[2][2]));
//     /* RD_RL */ \
//     X(RD_RLs)  X(RU_RLs)  X(TD_RLs)  X(TU_RLs) \
//     /* FR (实际先计算ZR，再递推到FR) */ \
//     X(RD_FRs)  X(RU_FRs)  X(TD_FRs)  X(TU_FRs) 
// #undef X

//     // 初始化透射矩阵
//     for(size_t iy=0; iy<mod1d->n; ++iy){
//     #define X(F) F[iy][0][0] = F[iy][1][1] = 1.0;
//         X(TD_RLs)  X(TU_RLs)  X(TD_FRs)  X(TU_FRs)  
//     #undef X
//     }

//     // 自由表面的反射系数
//     cplx_t R_tilt[2][2]; // SH波在自由表面的反射系数为1，不必定义变量
    
//     // 定义物理层内的反射透射系数矩阵，相对于界面上的系数矩阵增加了时间延迟因子
//     cplx_t RD[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
//     cplx_t RU[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
//     cplx_t TD[2][2] = GRT_INIT_IDENTITY_2x2_MATRIX;
//     cplx_t TU[2][2] = GRT_INIT_IDENTITY_2x2_MATRIX;

//     // 模型参数
//     // 后缀0，1分别代表上层和下层
//     real_t thk0, thk1, Rho0, Rho1;
//     cplx_t mu0, mu1;
//     cplx_t xa0=0.0, xb0=0.0, xa1=0.0, xb1=0.0;
//     cplx_t cbcb0=0.0, cbcb1=0.0;
//     cplx_t caca1=0.0;
//     cplx_t top_xa=0.0, top_xb=0.0;
//     cplx_t top_cbcb=0.0;

//     // 相速度
//     cplx_t c_phase = omega / k;

//     // 循环每一层
//     for(size_t iy=0; iy<mod1d->n; ++iy){
//         // 赋值上层 
//         thk0 = thk1;
//         Rho0 = Rho1;
//         mu0 = mu1;
//         xa0 = xa1;
//         xb0 = xb1;
//         cbcb0 = cbcb1;

//         // 更新模型参数
//         thk1 = mod1d->Thk[iy];
//         Rho1 = mod1d->Rho[iy];
//         mu1 = mod1d->mu[iy];
//         grt_get_mod1d_xa_xb(mod1d, iy, c_phase, &caca1, &xa1, &cbcb1, &xb1);

//         if(0==iy){
//             top_xa = xa1;
//             top_xb = xb1;
//             top_cbcb = cbcb1;
//             continue;
//         }

//         // 计算该层的反射透射系数
//         // 对第iy层的系数矩阵赋值，加入时间延迟因子(第iy-1界面与第iy界面之间)
//         grt_RT_matrix_PSV(
//             Rho0, xa0, xb0, cbcb0, mu0, 
//             Rho1, xa1, xb1, cbcb1, mu1, 
//             omega, k, 
//             RD, RU, TD, TU, stats);
//         grt_delay_RT_matrix_PSV(xa0, xb0, thk0, k, RD, RU, TD, TU);
//         if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
        
//         // FR
//         if(iy == 1){
//             GRT_RT_PSV_ASSIGN(FRs[iy]);
//         }
//         else{
//             grt_recursion_RT_matrix_PSV(
//                 RD_FRs[iy-1], RU_FRs[iy-1], TD_FRs[iy-1], TU_FRs[iy-1],
//                 RD, RU, TD, TU,
//                 RD_FRs[iy], RU_FRs[iy], TD_FRs[iy], TU_FRs[iy], stats);
//         }
        
//         // RL
//         GRT_RT_PSV_ASSIGN(RLs[iy-1]);
//         for(size_t my=0; my<iy-1; ++my){
//             grt_recursion_RT_matrix_PSV(
//                 RD_RLs[my], RU_RLs[my], TD_RLs[my], TU_RLs[my],
//                 RD, RU, TD, TU,
//                 RD_RLs[my], RU_RLs[my], TD_RLs[my], TU_RLs[my], stats);
//             if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
//         }
//     } // END for loop 
//     //===================================================================================

//     // 递推 RU_FR
//     grt_topfree_RU_PSV(top_xa, top_xb, top_cbcb, k, R_tilt, stats);
//     if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
//     for(size_t iy=0; iy<mod1d->n; ++iy){
//         grt_recursion_RU_PSV(
//             R_tilt, RD_FRs[iy], RU_FRs[iy], TD_FRs[iy],
//             TU_FRs[iy],
//             RU_FRs[iy], NULL, stats);
//         if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
//     }
    
//     // 赋值
//     memcpy(RD_RL0, RD_RLs, sizeof(cplx_t)*mod1d->n*4);
//     memcpy(RU_FR0, RU_FRs, sizeof(cplx_t)*mod1d->n*4);

//     BEFORE_RETURN:

// #define X(F) GRT_SAFE_FREE_PTR(F);
//     X(RD_RLs)  X(RU_RLs)  X(TD_RLs)  X(TU_RLs) \
//     X(RD_FRs)  X(RU_FRs)  X(TD_FRs)  X(TU_FRs) 
// #undef X

//     return;
// }


// void grt_get_mod_all_RD_RU_Love(
//     const GRT_MODEL1D *mod1d, cplx_t omega, real_t k,
//     cplx_t RDL_RL0[mod1d->n], cplx_t RUL_FR0[mod1d->n], int *stats)
// {
//     // 先置零
//     if(RDL_RL0!=NULL) memset(RDL_RL0, 0, sizeof(cplx_t)*mod1d->n); 
//     if(RUL_FR0!=NULL) memset(RUL_FR0, 0, sizeof(cplx_t)*mod1d->n); 

// #define X(F) cplx_t *F = (cplx_t *)calloc(mod1d->n, sizeof(cplx_t));
//     /* RDL_RL */ \
//     X(RDL_RLs)  X(RUL_RLs)  X(TDL_RLs)  X(TUL_RLs) \
//     /* FR (实际先计算ZR，再递推到FR) */ \
//     X(RDL_FRs)  X(RUL_FRs)  X(TDL_FRs)  X(TUL_FRs) 
// #undef X

//     // 初始化透射矩阵
//     for(size_t iy=0; iy<mod1d->n; ++iy){
//     #define X(F) F[iy] = 1.0;
//         X(TDL_RLs)  X(TUL_RLs)  X(TDL_FRs)  X(TUL_FRs)  
//     #undef X
//     }

//     // 定义物理层内的反射透射系数矩阵，相对于界面上的系数矩阵增加了时间延迟因子
//     cplx_t RDL = 0.0;
//     cplx_t RUL = 0.0;
//     cplx_t TDL = 1.0;
//     cplx_t TUL = 1.0;

//     // 模型参数
//     // 后缀0，1分别代表上层和下层
//     real_t thk0, thk1;
//     cplx_t mu0, mu1;
//     cplx_t cbcb1;
//     cplx_t xb0, xb1;

//     // 相速度
//     cplx_t c_phase = omega / k;

//     for(size_t iy=0; iy<mod1d->n; ++iy){
//         // 赋值上层 
//         thk0 = thk1;
//         mu0 = mu1;
//         xb0 = xb1;

//         // 更新模型参数
//         thk1 = mod1d->Thk[iy];
//         mu1 = mod1d->mu[iy];
//         grt_get_mod1d_xa_xb(mod1d, iy, c_phase, NULL, NULL, &cbcb1, &xb1);

//         if(0==iy){
//             continue;
//         }

//         // 计算该层的反射透射系数
//         // 对第iy层的系数矩阵赋值，加入时间延迟因子(第iy-1界面与第iy界面之间)
//         grt_RT_matrix_SH(
//             xb0, mu0, 
//             xb1, mu1, 
//             omega, k, 
//             &RDL, &RUL, &TDL, &TUL);
//         grt_delay_RT_matrix_SH(xb0, thk0, k, &RDL, &RUL, &TDL, &TUL);

//         // FR
//         if(iy == 1){
//             GRT_RT_SH_ASSIGN(FRs[iy]);
//         }
//         else{
//             grt_recursion_RT_matrix_SH(
//                 RDL_FRs[iy-1], RUL_FRs[iy-1], TDL_FRs[iy-1], TUL_FRs[iy-1],
//                 RDL, RUL, TDL, TUL,
//                 &RDL_FRs[iy], &RUL_FRs[iy], &TDL_FRs[iy], &TUL_FRs[iy], stats);
//         }

//         // RL
//         GRT_RT_SH_ASSIGN(RLs[iy-1]);
//         for(size_t my=0; my<iy-1; ++my){
//             grt_recursion_RT_matrix_SH(
//                 RDL_RLs[my], RUL_RLs[my], TDL_RLs[my], TUL_RLs[my],
//                 RDL, RUL, TDL, TUL,
//                 &RDL_RLs[my], &RUL_RLs[my], &TDL_RLs[my], &TUL_RLs[my], stats);
//             if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
//         }

//     } // END for loop 
//     //===================================================================================

//     // 递推 RU_FR
//     for(size_t iy=0; iy<mod1d->n; ++iy){
//         grt_recursion_RU_SH(
//             1.0, RDL_FRs[iy], RUL_FRs[iy], TDL_FRs[iy],
//             TUL_FRs[iy],
//             &RUL_FRs[iy], NULL, stats);
//         if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
//     }

//     // 赋值
//     memcpy(RDL_RL0, RDL_RLs, sizeof(cplx_t)*mod1d->n);
//     memcpy(RUL_FR0, RUL_FRs, sizeof(cplx_t)*mod1d->n);

//     BEFORE_RETURN:

// #define X(F) GRT_SAFE_FREE_PTR(F);
//     X(RDL_RLs)  X(RUL_RLs)  X(TDL_RLs)  X(TUL_RLs) \
//     X(RDL_FRs)  X(RUL_FRs)  X(TDL_FRs)  X(TUL_FRs) 
// #undef X

//     return;
// }


void grt_kernel_Rayl(GRT_MODEL1D *mod1d, const real_t k, const size_t iref)
{
    mod1d->k = k;
    grt_mod1d_xa_xb(mod1d, k);

    // 广义层的反射透射系数
    RT_MATRIX *M_BL = &mod1d->M_BL;
    RT_MATRIX *M_RS = &mod1d->M_RS;
    RT_MATRIX *M_FA = &mod1d->M_FA;
    // 顶底界面
    RT_MATRIX *M_top = &mod1d->M_top;
    RT_MATRIX *M_bot = &mod1d->M_bot;

    grt_reset_RT_matrix_PSV(M_BL);
    grt_reset_RT_matrix_PSV(M_RS);
    grt_reset_RT_matrix_PSV(M_FA);
    grt_reset_RT_matrix_PSV(M_top);
    grt_reset_RT_matrix_PSV(M_bot);

    // 从顶到底，计算RD_RL
    for(size_t iy = 1; iy < mod1d->n-1; ++iy)
    {
        // 定义物理层内的反射透射系数矩阵
        RT_MATRIX *M = &(RT_MATRIX){};
        grt_reset_RT_matrix_PSV(M);

        // 对第iy层的系数矩阵赋值
        grt_RT_matrix_PSV(mod1d, iy, M);
        if(M->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

        // TEMP!! 可以再划分为 PSV 和 SH
        // 加入时间延迟因子(第iy-1界面与第iy界面之间)
        grt_delay_RT_matrix_PSV(mod1d, iy, M);

        // FR
        if(iy <= iref){
            if(iy == 1){
                memcpy(M_FA, M, sizeof(*M));
            } else {
                grt_recursion_RT_matrix_PSV(M_FA, M, M_FA);  // 写入原地址
                if(M_FA->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
            }
        }
        // RL
        else {
            if(iy == iref+1){
                memcpy(M_BL, M, sizeof(*M));
            } else {
                grt_recursion_RT_matrix_PSV(M_BL, M, M_BL);  // 写入原地址
                if(M_BL->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
            } 
        }
    } // END for loop 
    //===================================================================================

    // 递推 RU_FR
    grt_topbound_RU_PSV(mod1d);
    if(M_top->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
    grt_recursion_RU_PSV(M_top, M_FA, M_FA);
    if(M_FA->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

    // 递推 RD_RL
    grt_botbound_RD_PSV(mod1d);
    if(M_bot->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
    grt_recursion_RD_PSV(M_BL, M_bot, M_BL);
    if(M_BL->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

    BEFORE_RETURN:
    return;
}


void grt_kernel_Love(GRT_MODEL1D *mod1d, const real_t k, const size_t iref)
{
    mod1d->k = k;
    grt_mod1d_xa_xb(mod1d, k);

    // 广义层的反射透射系数
    RT_MATRIX *M_BL = &mod1d->M_BL;
    RT_MATRIX *M_RS = &mod1d->M_RS;
    RT_MATRIX *M_FA = &mod1d->M_FA;
    // 顶底界面
    RT_MATRIX *M_top = &mod1d->M_top;
    RT_MATRIX *M_bot = &mod1d->M_bot;

    grt_reset_RT_matrix_SH(M_BL);
    grt_reset_RT_matrix_SH(M_RS);
    grt_reset_RT_matrix_SH(M_FA);
    grt_reset_RT_matrix_SH(M_top);
    grt_reset_RT_matrix_SH(M_bot);

    // 从顶到底，计算RD_RL
    for(size_t iy = 1; iy < mod1d->n-1; ++iy)
    {
        // 定义物理层内的反射透射系数矩阵
        RT_MATRIX *M = &(RT_MATRIX){};
        grt_reset_RT_matrix_SH(M);

        // 对第iy层的系数矩阵赋值
        grt_RT_matrix_SH(mod1d, iy, M);
        
        // fprintf(stderr, "RDL="GRT_CMPLX_FMT"\n", GRT_CMPLX_SPLIT(M->RDL));
        // fprintf(stderr, "RUL="GRT_CMPLX_FMT"\n", GRT_CMPLX_SPLIT(M->RUL));
        // fprintf(stderr, "TDL="GRT_CMPLX_FMT"\n", GRT_CMPLX_SPLIT(M->TDL));
        // fprintf(stderr, "TUL="GRT_CMPLX_FMT"\n", GRT_CMPLX_SPLIT(M->TUL));
        if(M->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

        // TEMP!! 可以再划分为 PSV 和 SH
        // 加入时间延迟因子(第iy-1界面与第iy界面之间)
        grt_delay_RT_matrix_SH(mod1d, iy, M);

        // FR
        if(iy <= iref){
            if(iy == 1){
                memcpy(M_FA, M, sizeof(*M));
            } else {
                grt_recursion_RT_matrix_SH(M_FA, M, M_FA);  // 写入原地址
                if(M_FA->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
            }
        }
        // RL
        else {
            if(iy == iref+1){
                memcpy(M_BL, M, sizeof(*M));
            } else {
                grt_recursion_RT_matrix_SH(M_BL, M, M_BL);  // 写入原地址
                if(M_BL->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
            } 
        }
        
        // fprintf(stderr, "iy=%zu, RDL=%f\n", iy, creal(mod1d->M_BL.RDL));
    } // END for loop 
    //===================================================================================

    // 递推 RU_FR
    grt_topbound_RU_SH(mod1d);
    if(M_top->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
    grt_recursion_RU_SH(M_top, M_FA, M_FA);
    if(M_FA->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

    // 递推 RD_RL
    grt_botbound_RD_SH(mod1d);
    if(M_bot->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
    grt_recursion_RD_SH(M_BL, M_bot, M_BL);
    if(M_BL->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

    BEFORE_RETURN:
    return;
}



void grt_secular_function_potential_Rayl(
    GRT_MODEL1D *mod1d, const real_t cphase,
    const int secRaylType, const size_t iref, cplx_t *psec, cplx_t *ppot)
{
    // iref > 0 时， secRaylType必须等于1，否则报错
    if(iref > 0 && secRaylType != 1){
        fprintf(stderr, "iref > 0 && secRaylType != 1\n");
        exit(EXIT_FAILURE);
    }

    real_t k = mod1d->omega/cphase;
    grt_kernel_Rayl(mod1d, k, iref);

    cplx_t tmp1[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    cplx_t tmp2[2][2] = GRT_INIT_ZERO_2x2_MATRIX;

    if(secRaylType==0) {
        // 对于首层，改为 D21*RD_RL + D22
        // 计算顶层的D22 
        cplx_t D22[2][2], D21[2][2];
        cplx_t top_xa=0.0, top_xb=0.0;
        cplx_t top_caca=0.0, top_cbcb=0.0;

        top_xa = mod1d->xa[0];
        top_xb = mod1d->xb[0];
        top_caca = mod1d->caca[0];
        top_cbcb = mod1d->cbcb[0];

        // 使用归一化的 D21, D22
        D21[0][0] = 1.0 - 0.5*top_cbcb;   D21[0][1] = top_xb;
        D21[1][0] = top_xa;               D21[1][1] = 1.0 - 0.5*top_cbcb;

        D22[0][0] = 1.0 - 0.5*top_cbcb;   D22[0][1] = - top_xb;
        D22[1][0] = - top_xa;             D22[1][1] = 1.0 - 0.5*top_cbcb;

        // grt_get_layer_D21(top_xa, top_xb, top_kbkb, top_mu, omega, k, D21);
        // grt_get_layer_D22(top_xa, top_xb, top_kbkb, top_mu, omega, k, D22);
        
        // 先适当缩放
        grt_cmat2x2_k(mod1d->M_BL.RD, 1.0/(k*k), tmp1);

        grt_cmat2x2_mul(D21, tmp1, tmp1);
        // 逆缩放
        grt_cmat2x2_k(tmp1, (k*k), tmp1);

        grt_cmat2x2_add(D22, tmp1, tmp1);
        // 适当乘上一个系数，将久期函数的值域置于合适的范围（大致在[-1,1]），方便后续找零点判断幅值
        grt_cmat2x2_k(tmp1, 1.0/(2.0*sqrt(top_caca*top_cbcb)), tmp1);
    }
    else {
        // grt_cmat2x2_mul(RU_FR, RD_RL, tmp1);
        // grt_cmat2x2_one_sub(tmp1);

        // 先适当缩放
        grt_cmat2x2_k(mod1d->M_BL.RD, 1.0/(k*k), tmp1);
        grt_cmat2x2_k(mod1d->M_FA.RU, 1.0/(k*k), tmp2);

        grt_cmat2x2_mul(tmp2, tmp1, tmp1);
        // 逆缩放
        grt_cmat2x2_k(tmp1, (k*k*k*k), tmp1);
        grt_cmat2x2_one_sub(tmp1);
    }
    *psec = tmp1[0][0]*tmp1[1][1] - tmp1[0][1]*tmp1[1][0];

    // TEMP!!!
    // 这里先临时加个处理，不使用if-else调整代码
    // 如果 iref 位于液体层，则 RD、RU 仅在 [0][0] 位置有值，久期函数的计算要改变
    if(mod1d->isLiquid[iref]){
        // 不管 secRaylType
        *psec = 1.0 - mod1d->M_FA.RU[0][0]*mod1d->M_BL.RD[0][0];
    }

    // 返回对应的垂直波函数
    if(ppot != NULL){
        // 假设一个比例
        ppot[2] = - tmp1[0][1] / sqrt( GRT_SQUARE(fabs(tmp1[0][0])) + GRT_SQUARE(fabs(tmp1[0][1])) ); 
        ppot[3] = + tmp1[0][0] / sqrt( GRT_SQUARE(fabs(tmp1[0][0])) + GRT_SQUARE(fabs(tmp1[0][1])) );
        grt_cmat2x1_mul(mod1d->M_BL.RD, ppot+2, ppot);

        // TEMP!!!
        // 这里先临时加个处理，不使用if-else调整代码
        // 如果 iref 位于液体层，则 RD、RU 仅在 [0][0] 位置有值，久期函数的计算要改变
        if(mod1d->isLiquid[iref]){
            // 不管 secRaylType
            ppot[2] = 1.0;
            ppot[3] = 0.0;
            ppot[0] = mod1d->M_BL.RD[0][0];
            ppot[1] = 0.0;
        }
    }

    return;
}


void grt_secular_function_potential_Love(
    GRT_MODEL1D *mod1d, const real_t cphase,
    const size_t iref, cplx_t *psec, cplx_t *ppot)
{
    real_t k = mod1d->omega/cphase;
    grt_kernel_Love(mod1d, k, iref);
    
    // 1 - RUL_FR*RDL_RL
    // 对于首层，1 - RDL_RL
    // fprintf(stderr, "%f, %f, %f\n", creal(mod1d->c_phase), creal(mod1d->M_FA.RUL), creal(mod1d->M_BL.RDL));
    *psec = 1.0 - mod1d->M_FA.RUL*mod1d->M_BL.RDL;

    // TEMP!!!
    // 这里先临时加个处理，不使用if-else调整代码
    // 对于液体层不应该运行该函数
    if(mod1d->isLiquid[iref]){
        GRTRaiseError("Wrong execution.\n");
    }

    // 返回对应的垂直波函数
    if(ppot != NULL){
        // 假设一个值
        ppot[1] = 1.0;
        ppot[0] = mod1d->M_BL.RDL;
    }

    return;
}

void grt_secular_function_potential(
    GRT_MODEL1D *mod1d, const real_t cphase,
    const int secRaylType, const size_t iref, const DISPER_TYPE wtype, cplx_t *psec, cplx_t *ppot)
{
    if(wtype == GRT_DISPERSION_RAYL){
        grt_secular_function_potential_Rayl(mod1d, cphase, secRaylType, iref, psec, ppot);
    } 
    else if(wtype == GRT_DISPERSION_LOVE){
        grt_secular_function_potential_Love(mod1d, cphase, iref, psec, ppot);
    }
    else {
        GRTRaiseError("Wrong execution.");
    }
}

void grt_secular_function(
    GRT_MODEL1D *mod1d, const real_t cphase,
    const int secRaylType, const size_t iref, const DISPER_TYPE wtype, cplx_t *psec)
{
    grt_secular_function_potential(mod1d, cphase, secRaylType, iref, wtype, psec, NULL);
}





size_t grt_goldensection_search_root(
    GRT_MODEL1D *mod1d, const real_t c1, const real_t c2, 
    const int secRaylType, const size_t iref, const DISPER_TYPE wtype,
    const real_t rtol, cplx_t *root_sec, real_t *root_c)
{
    real_t c[2] = {c1, c2};
    real_t x[2];
    x[0] = c[1] - GOLDEN_RATIO*(c[1] - c[0]);
    x[1] = c[0] + GOLDEN_RATIO*(c[1] - c[0]);
    cplx_t sec[2];
    real_t secABS[2];

    for(int i=0; i<2; ++i){
        grt_secular_function(mod1d, x[i], secRaylType, iref, wtype, &sec[i]);
        secABS[i] = fabs(sec[i]);
    }
    size_t is=0;

    #if __DEBUG_SECULAR__ == 1
    fprintf(stderr, "# f=%f, root_sec2=%.9e%+.9ej|%.9e%+.9ej, c1=%.16e, c2=%.16e, dc=%.16e, root_c=%.9e, %d\n", 
            creal(mod1d->omega)/PI2, creal(sec[0]), cimag(sec[0]), creal(sec[1]), cimag(sec[1]), x[0], x[1], x[1]-x[0], 
            x[0], is);
    fflush(stderr);
    #endif
    
    do{
        if(secABS[1] > secABS[0]){
            c[1] = x[1];
            sec[1] = sec[0];
            secABS[1] = secABS[0];
            x[1] = x[0];
            x[0] = c[1] - GOLDEN_RATIO*(c[1] - c[0]);
            grt_secular_function(mod1d, x[0], secRaylType, iref, wtype, &sec[0]);
            secABS[0] = fabs(sec[0]);
        } else {
            c[0] = x[0];
            sec[0] = sec[1];
            secABS[0] = secABS[1];
            x[0] = x[1];
            x[1] = c[0] + GOLDEN_RATIO*(c[1] - c[0]);
            grt_secular_function(mod1d, x[1], secRaylType, iref, wtype, &sec[1]);
            secABS[1] = fabs(sec[1]);
        }

        // 接近达到浮点数精度极限
        if(c[1] - c[0] < 1e-15)  break;

        ++is;
        #if __DEBUG_SECULAR__ == 1
        fprintf(stderr, "# f=%f, root_sec2=%.9e%+.9ej|%.9e%+.9ej, secABS=%.9e, c1=%.16e, c2=%.16e, dc=%.19e, root_c=%.9e, %d\n", 
            creal(mod1d->omega)/PI2, creal(sec[0]), cimag(sec[0]), creal(sec[1]), cimag(sec[1]), secABS[0], x[0], x[1], x[1]-x[0], 
            x[0], is);
        fflush(stderr);
        #endif
        
    } while(is < 800);

    #if __DEBUG_SECULAR__ == 1
    fprintf(stderr, "# f=%f, root_sec2=%.9e%+.9ej|%.9e%+.9ej, dc=%e, root_c=%.9e, %d\n", 
            creal(mod1d->omega)/PI2, creal(sec[0]), cimag(sec[0]), creal(sec[1]), cimag(sec[1]), x[1]-x[0], 
            x[0], is);
    fflush(stderr);
    #endif

    // 是否为零点
    if(secABS[0] < rtol){
        *root_sec = sec[0];
        *root_c = x[0];
        #if __DEBUG_SECULAR__ == 1
        fprintf(stderr, "#      f=%f, root_sec=%.9e%+.9ej, root_c=%.9e, %d\n", creal(mod1d->omega)/PI2, creal(sec[0]), cimag(sec[0]), x[0], is);
        fflush(stderr);
        #endif
    } else {
        *root_sec = -999.0;
        *root_c = -1.0;
    }

    return is;
}



// /** 
//  * （仅备份）二分区间法确定久期函数零点 
//  * 
//  * @param[in]      mod1d        模型结构体指针
//  * @param[in]      w            圆频率
//  * @param[in]      kmid         区间中点波数
//  * @param[in]      dk0          对称区间长度一半
//  * @param[in]      secmid       区间中点久期函数值
//  * @param[in]      secmidABS    区间中点久期函数绝对值
//  * @param[in]      secRaylType  Rayl久期函数类型
//  * @param[in]      iref         久期函数层位
//  * @param[in]      isRayl       Rayleigh or Love
//  * @param[out]     root_sec     零点的久期函数值
//  * @param[out]     root_k       零点的波数
//  * 
//  */
// size_t grt_bisearch_root(
//     GRT_MODEL1D *mod1d, const real_t kmid, const real_t dk0,
//     const cplx_t secmid, const real_t secmidABS, const int secRaylType, const size_t iref, const bool isRayl,
//     const real_t rtol, cplx_t *root_sec, real_t *root_k)
// {
//     real_t k[3], dk;
//     cplx_t sec[3];
//     real_t secABS[3];
//     bool recompute_sec[3] = {false, false, false};
//     // MYINT dot3stat = 0;
//     int inv_stats=GRT_INVERSE_SUCCESS;

//     dk = dk0;
//     k[1] = kmid;
//     sec[1] = secmid;
//     secABS[1] = secmidABS;
//     dk *= 0.5;
//     for(int i=0; i<=2; i+=2){
//         k[i] = k[1] + (i-1)*dk;
//         recompute_sec[i] = true;
//     }

//     size_t is=0;
//     // char stats_str[2000]={'\0'};
//     do{
//         // 重新计算三点值
//         for(int i=0; i<3; ++i){
//             if(!recompute_sec[i])  continue;

//             grt_secular_function(mod1d, w, k[i], secRaylType, iref, isRayl, &sec[i], &inv_stats);
//             secABS[i] = fabs(sec[i]);
//             recompute_sec[i] = false;
//         }


//         // 是否为零点，设计一个特别小的tol（tol < zero_tol）用于提前退出循环
//         // 由于浮点数计算误差，零点的幅值不会无限逼近0，此处的tol是理想化误差，某些情况下的零点够不着
//         if(secABS[1] < GRT_MIN(1e-8, GRT_ZEROROOT_TOL))  break;
        
//         // 判断修改方向，定义新的三点
//         // 中间点为极小值点
//         if(secABS[0] >= secABS[1] && secABS[1] <= secABS[2]){
//             dk *= 0.5;
//             // 更新两边
//             for(int i=0; i<=2; i+=2){
//                 k[i] = k[1] + (i-1)*dk;
//                 recompute_sec[i] = true;
//             }
//             // dot3stat = 0;
//             ++is;
//         }
//         // 左侧为极小值点
//         else if(secABS[0] < secABS[1] && secABS[1] < secABS[2]){ 
//             for(int i=2; i>0; --i){
//                 k[i] = k[i-1];
//                 sec[i] = sec[i-1];
//                 secABS[i] = secABS[i-1];
//             }
//             k[0] = k[1] - dk;
//             recompute_sec[0] = true;
//             // dot3stat = 1;
//         }
//         // 右侧为极小值点
//         else if(secABS[0] > secABS[1] && secABS[1] > secABS[2]){ 
//             for(int i=0; i<2; ++i){
//                 k[i] = k[i+1];
//                 sec[i] = sec[i+1];
//                 secABS[i] = secABS[i+1];
//             }
//             k[2] = k[1] + dk;
//             recompute_sec[2] = true;
//             // dot3stat = 2;
//         }
//         else {
//             break;
//         }

//         // sprintf(stats_str, "%s %d", stats_str, dot3stat);
        
//     } while(is < 50);

//     #if __DEBUG_SECULAR__ == 1
//     fprintf(stderr, "# f=%f, root_sec3=%e%+ej|%e%+ej|%e%+ej, dk=%e, root_c=%f, %d\n", 
//             w/PI2, creal(sec[0]), cimag(sec[0]), creal(sec[1]), cimag(sec[1]), creal(sec[2]), cimag(sec[2]), k[2]-k[0],
//             w/k[1], is);
//     fflush(stderr);
//     #endif

//     // 是否为零点
//     if(secABS[1] < GRT_ZEROROOT_TOL){
//         *root_sec = sec[1];
//         *root_k = k[1];
//         #if __DEBUG_SECULAR__ == 1
//         fprintf(stderr, "#      f=%f, root_sec=%e%+ej, root_c=%f, %d\n", w/PI2, creal(sec[1]), cimag(sec[1]), w/k[1], is);
//         fflush(stderr);
//         #endif
//     } else {
//         *root_sec = -999.0;
//         *root_k = -1.0;
//     }

//     return is;
// }



// /** 备份，固定间隔打印久期函数 */
// int get_secular(int argc, char **argv){
//     char *command = argv[0];

//     char *s_modelpath = argv[1];
//     real_t freq = atof(argv[2]);
//     GRT_PYMODEL1D *pymod;
//     // 读入模型文件（暂设置为不支持液体层）
//     if((pymod = read_pymod_from_file(command, s_modelpath, -1, -1, true)) ==NULL){
//         exit(EXIT_FAILURE);
//     }
//     // pymod1d -> mod1d
//     GRT_MODEL1D *mod1d = init_mod1d(pymod->n);
//     get_mod1d(pymod, mod1d);


//     // Vs最大最小速度
//     real_t vbmax = pymod->Vb[findMinMax_real_t(pymod->Vb, pymod->n, true)];
//     real_t vbmin = pymod->Vb[findMinMax_real_t(pymod->Vb, pymod->n, false)];

//     // print_mod1d(mod1d);

//     cplx_t secRayl, secLove;
//     MYINT secRaylType = 1;
//     MYINT iref = 4;

//     real_t omega, k;
//     omega = PI2*freq;

//     update_omega(mod1d, omega);
//     // print_mod1d(mod1d);
//     // getchar();

//     MYINT inv_stats=GRT_INVERSE_SUCCESS;

//     for(real_t c=vbmin*0.9; c<vbmax*1.1; c+=0.0001){
//     // for(real_t c=3.4; c<3.9; c+=1e-5){
//         k = omega/c;
//         secular_function(mod1d, omega, k, secRaylType, iref, &secRayl, &secLove, &inv_stats);
//         printf("%e %e %e %e %e \n", c, creal(secRayl), cimag(secRayl), creal(secLove), cimag(secLove));
//     }

//     exit(EXIT_SUCCESS);
// }