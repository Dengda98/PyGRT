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


#define __DEBUG_SECULAR__  0


void grt_GRT_matrix_allLayer_Rayl(GRT_MODEL1D *mod1d, const real_t k, RT_MATRIX *Mall_RL, RT_MATRIX *Mall_FR)
{
    mod1d->k = k;
    grt_mod1d_xa_xb(mod1d, k);

    // 顶底界面
    RT_MATRIX *M_top = &mod1d->M_top;
    RT_MATRIX *M_bot = &mod1d->M_bot;

    for(size_t iy = 0; iy < mod1d->n; ++iy){
        grt_reset_RT_matrix_PSV(&Mall_RL[iy]);
        grt_reset_RT_matrix_PSV(&Mall_FR[iy]);
    }

    // 循环每一层
    for(size_t iy = 1; iy < mod1d->n-1; ++iy)
    {
        // 定义物理层内的反射透射系数矩阵
        RT_MATRIX *M = &(RT_MATRIX){};
        grt_reset_RT_matrix_PSV(M);

        // 对第iy层的系数矩阵赋值
        grt_RT_matrix_PSV(mod1d, iy, M);
        
        // 加入时间延迟因子(第iy-1界面与第iy界面之间)
        grt_delay_RT_matrix_PSV(mod1d, iy, M);
        
        // FR
        if(iy == 1){
            memcpy(&Mall_FR[iy], M, sizeof(*M));
        }
        else{
            grt_recursion_RT_matrix_PSV(&Mall_FR[iy-1], M, &Mall_FR[iy]);
        }
        
        // RL
        memcpy(&Mall_RL[iy-1], M, sizeof(*M));
        for(size_t my = 0; my < iy-1; ++my){
            grt_recursion_RT_matrix_PSV(&Mall_RL[my], M, &Mall_RL[my]);
        }
    } // END for loop 
    //===================================================================================

    // 递推 RU_FR
    grt_topbound_RU_PSV(mod1d);
    for(size_t iy = 0; iy < mod1d->n; ++iy){
        grt_recursion_RU_PSV(M_top, &Mall_FR[iy], &Mall_FR[iy]);
    }

    // ---------- 只算到 n-1 层 --------------
    // 递推 RU_FR
    grt_topbound_RU_PSV(mod1d);
    for(size_t iy = 0; iy < mod1d->n-1; ++iy){
        grt_recursion_RU_PSV(M_top, &Mall_FR[iy], &Mall_FR[iy]);
    }

    // 递推 RD_RL
    grt_botbound_RD_PSV(mod1d);
    for(size_t iy = 0; iy < mod1d->n-1; ++iy){
        grt_recursion_RD_PSV(&Mall_RL[iy], M_bot, &Mall_RL[iy]);
    }
    
}


void grt_GRT_matrix_allLayer_Love(GRT_MODEL1D *mod1d, const real_t k, RT_MATRIX *Mall_RL, RT_MATRIX *Mall_FR)
{
    mod1d->k = k;
    grt_mod1d_xa_xb(mod1d, k);

    // 顶底界面
    RT_MATRIX *M_top = &mod1d->M_top;
    RT_MATRIX *M_bot = &mod1d->M_bot;

    for(size_t iy = 0; iy < mod1d->n; ++iy){
        grt_reset_RT_matrix_SH(&Mall_RL[iy]);
        grt_reset_RT_matrix_SH(&Mall_FR[iy]);
    }

    // 循环每一层
    for(size_t iy = 1; iy < mod1d->n-1; ++iy)
    {
        // 定义物理层内的反射透射系数矩阵
        RT_MATRIX *M = &(RT_MATRIX){};
        grt_reset_RT_matrix_SH(M);

        // 对第iy层的系数矩阵赋值
        grt_RT_matrix_SH(mod1d, iy, M);
        
        // 加入时间延迟因子(第iy-1界面与第iy界面之间)
        grt_delay_RT_matrix_SH(mod1d, iy, M);
        
        // FR
        if(iy == 1){
            memcpy(&Mall_FR[iy], M, sizeof(*M));
        }
        else{
            grt_recursion_RT_matrix_SH(&Mall_FR[iy-1], M, &Mall_FR[iy]);
        }
        
        // RL
        memcpy(&Mall_RL[iy-1], M, sizeof(*M));
        for(size_t my = 0; my < iy-1; ++my){
            grt_recursion_RT_matrix_SH(&Mall_RL[my], M, &Mall_RL[my]);
        }
    } // END for loop 
    //===================================================================================

    // ---------- 只算到 n-1 层 --------------
    // 递推 RU_FR
    grt_topbound_RU_SH(mod1d);
    for(size_t iy = 0; iy < mod1d->n-1; ++iy){
        grt_recursion_RU_SH(M_top, &Mall_FR[iy], &Mall_FR[iy]);
    }

    // 递推 RD_RL
    grt_botbound_RD_SH(mod1d);
    for(size_t iy = 0; iy < mod1d->n-1; ++iy){
        grt_recursion_RD_SH(&Mall_RL[iy], M_bot, &Mall_RL[iy]);
    }
}

void grt_GRT_matrix_Rayl(GRT_MODEL1D *mod1d, const real_t k, const size_t iref)
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


void grt_GRT_matrix_Love(GRT_MODEL1D *mod1d, const real_t k, const size_t iref)
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
    grt_GRT_matrix_Rayl(mod1d, k, iref);

    cplx_t tmp1[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    cplx_t tmp2[2][2] = GRT_INIT_ZERO_2x2_MATRIX;

    if(secRaylType==0) {
        // !! WARNING
        // 这里这种改写仅适用于顶层为自由表面边界条件的情况
        // 例如对于刚性条件就需要使用 D11 和 D12 矩阵

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
    grt_GRT_matrix_Love(mod1d, k, iref);
    
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