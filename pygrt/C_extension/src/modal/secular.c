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

#include "grt/common/RT_matrix.h"
#include "grt/common/model.h"
#include "grt/common/matrix.h"
#include "grt/common/search.h"

#include "grt/modal/secular.h"


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
        // 界面两侧性质不同（固-液），计算广义矩阵时跳过该层，在组成久期函数时会单独处理该界面
        if(iy == iref && (mod1d->isLiquid[iref-1] ^ mod1d->isLiquid[iref]) == 1){
            continue;
        }

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

real_t grt_secular_function_cbegin(const GRT_MODEL1D *mod1d, const size_t iref, DISPER_TYPE wtype)
{
    real_t cbegin = 0.0;

    if(wtype == GRT_DISPERSION_RAYL){
        // 逻辑应与 grt_secular_function_potential_Rayl 函数相互呼应
        // 界面两侧均为固体
        if(! mod1d->isLiquid[iref] && (iref==0 || (! mod1d->isLiquid[iref-1]))){
            if(iref > 0)  cbegin = mod1d->Vb[iref];
        }
        // 界面两侧均为液体
        else if(mod1d->isLiquid[iref] && (iref==0 || (mod1d->isLiquid[iref-1]))){
            cbegin = mod1d->Va[iref];
        }
    }
    else if(wtype == GRT_DISPERSION_LOVE){
        cbegin = mod1d->Vb[iref];
    }
    
    return cbegin;
}

void grt_secular_function_potential_Rayl(
    GRT_MODEL1D *mod1d, const real_t cphase, const size_t iref, cplx_t *psec, cplx_t ppot[GRT_RAYL_DIM], cplx_t ppotUp[GRT_RAYL_DIM])
{
    real_t k = creal(mod1d->omega)/cphase;
    grt_GRT_matrix_Rayl(mod1d, k, iref);
    
    real_t cref = 0.0;

    // 界面两侧均为固体
    if(! mod1d->isLiquid[iref] && (iref==0 || (! mod1d->isLiquid[iref-1]))){
        cplx_t Det[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
        cref = mod1d->Vb[0];
        if(iref==0 && cphase < cref) {
            // !! WARNING
            // 这里这种改写仅适用于顶层为自由表面边界条件的情况
            // 例如对于刚性条件就需要使用 D11 和 D12 矩阵

            // 对于首层，改为 D21*RD_RL + D22
            // 计算顶层的D22 
            cplx_t D22[2][2], D21[2][2];
            cplx_t top_xa=0.0, top_xb=0.0;
            cplx_t top_cbcb=0.0;

            top_xa = mod1d->xa[0];
            top_xb = mod1d->xb[0];
            top_cbcb = mod1d->cbcb[0];

            // 使用归一化的 D21, D22
            D21[0][0] = 1.0 - 0.5*top_cbcb;   D21[0][1] = top_xb;
            D21[1][0] = top_xa;               D21[1][1] = 1.0 - 0.5*top_cbcb;

            D22[0][0] = 1.0 - 0.5*top_cbcb;   D22[0][1] = - top_xb;
            D22[1][0] = - top_xa;             D22[1][1] = 1.0 - 0.5*top_cbcb;

            grt_cmat2x2_mul(D21, mod1d->M_BL.RD, Det);
            grt_cmat2x2_add(D22, Det, Det);
            *psec = Det[0][0]*Det[1][1] - Det[0][1]*Det[1][0];
            // 适当控制久期函数幅值
            *psec /= 2.0;
        }
        else {
            grt_cmat2x2_mul(mod1d->M_FA.RU, mod1d->M_BL.RD, Det);
            grt_cmat2x2_one_sub(Det);
            *psec = Det[0][0]*Det[1][1] - Det[0][1]*Det[1][0];
        }

        // 返回对应的垂直波函数
        if(ppot != NULL){
            // 假设一个比例
            ppot[2] = - Det[0][1] / sqrt( GRT_SQUARE(fabs(Det[0][0])) + GRT_SQUARE(fabs(Det[0][1])) ); 
            ppot[3] = + Det[0][0] / sqrt( GRT_SQUARE(fabs(Det[0][0])) + GRT_SQUARE(fabs(Det[0][1])) );
            grt_cmat2x1_mul(mod1d->M_BL.RD, ppot+2, ppot);
        }
    }
    // 界面两侧均为液体
    else if(mod1d->isLiquid[iref] && (iref==0 || (mod1d->isLiquid[iref-1]))){
        *psec = 1.0 - mod1d->M_FA.RU[0][0]*mod1d->M_BL.RD[0][0];
        // 返回对应的垂直波函数
        if(ppot != NULL){
            ppot[2] = 1.0;
            ppot[3] = 0.0;
            ppot[0] = mod1d->M_BL.RD[0][0];
            ppot[1] = 0.0;
        }
    }
    // 界面两侧为液固/固液，此时 iref >= 1
    // Qin * diag{RU, RD} - Qout
    else {
        // TODO
        // 这一块儿的公式有一些变体比较有意思，固体-固体界面也适用，有空放到在线文档里
        cplx_t Qin[3][3] = {0}, Qout[3][3] = {0};
        grt_Q_mat3x3_ls_PSV(mod1d, iref, Qin, Qout);
        cplx_t R3[3][3] = {0};
        cplx_t RU_FA_delay[2][2] = {0};
        if(mod1d->isLiquid[iref-1]){
            cref = GRT_MIN(mod1d->Va[iref-1], mod1d->Vb[iref]);
            cplx_t ex2a = exp( - mod1d->k * mod1d->Thk[iref-1] * mod1d->xa[iref-1]);
            ex2a *= ex2a;
            R3[0][0] = RU_FA_delay[0][0] = ex2a * mod1d->M_FA.RU[0][0];
            R3[1][1] = mod1d->M_BL.RD[0][0];
            R3[1][2] = mod1d->M_BL.RD[0][1];
            R3[2][1] = mod1d->M_BL.RD[1][0];
            R3[2][2] = mod1d->M_BL.RD[1][1];
        } else {
            cref = GRT_MIN(mod1d->Vb[iref-1], mod1d->Va[iref]);
            real_t thk = mod1d->Thk[iref-1];
            cplx_t xa1 = mod1d->xa[iref-1];
            cplx_t xb1 = mod1d->xb[iref-1];
            real_t k = mod1d->k;
            
            cplx_t exa, exb, ex2a, ex2b, exab;
            exa = exp(- k*thk*xa1);
            exb = exp(- k*thk*xb1);
            ex2a = exa * exa;
            ex2b = exb * exb;
            exab = exa * exb;

            R3[0][0] = mod1d->M_BL.RD[0][0];
            R3[1][1] = RU_FA_delay[0][0] = ex2a * mod1d->M_FA.RU[0][0];
            R3[1][2] = RU_FA_delay[0][1] = exab * mod1d->M_FA.RU[0][1];
            R3[2][1] = RU_FA_delay[1][0] = exab * mod1d->M_FA.RU[1][0];
            R3[2][2] = RU_FA_delay[1][1] = ex2b * mod1d->M_FA.RU[1][1];
        }

        if(cphase < cref){
            grt_cmatmxn_mul(3, 3, 3, Qin, R3, R3);
            grt_cmatmxn_sub(3, 3, R3, Qout, R3);
            // 公式推导可证明，由于 Qin 和 Qout 在 [2][0] 恒为 0 
            // 再结合原始的分块矩阵 R3 的特点，此时计算结果也在 [2][0] 位置恒为 0
        } else {
            cplx_t RT[3][3] = {0};
            grt_RT_mat3x3_ls_PSV(mod1d, iref, RT);
            grt_cmatmxn_mul(3, 3, 3, RT, R3, R3);
            for(int i = 0; i < 3; ++i){
                for(int j = 0; j < 3; ++j){
                    if(i == j){
                        R3[i][j] = 1.0 - R3[i][j];
                    } else {
                        R3[i][j] = - R3[i][j];
                    }
                }
            }
        }

        // 3x3 矩阵行列式
        *psec =  R3[0][0] * (R3[1][1] * R3[2][2] - R3[1][2] * R3[2][1])
               - R3[0][1] * (R3[1][0] * R3[2][2] - R3[1][2] * R3[2][0])
               + R3[0][2] * (R3[1][0] * R3[2][1] - R3[1][1] * R3[2][0]);

        // 适当调整幅值
        if(cphase < cref && cphase > 1.0){
            *psec /= 4.0*cphase*cphase;
        }

        // 返回对应的垂直波函数
        if(ppot != NULL && ppotUp != NULL){
            // 做叉积得到基础解
            cplx_t tmp[3];
            tmp[0] = R3[1][1] * R3[2][2] - R3[1][2] * R3[2][1];
            tmp[1] = R3[1][2] * R3[2][0] - R3[1][0] * R3[2][2];
            tmp[2] = R3[1][0] * R3[2][1] - R3[1][1] * R3[2][0];
            real_t amp = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
            tmp[0] /= amp;
            tmp[1] /= amp;
            tmp[2] /= amp;

            if(mod1d->isLiquid[iref-1]){
                ppot[2] = tmp[1]; 
                ppot[3] = tmp[2];
                grt_cmat2x1_mul(mod1d->M_BL.RD, ppot + 2, ppot);
                ppotUp[0] = tmp[0];
                ppotUp[1] = 0.0;
                ppotUp[2] = RU_FA_delay[0][0] * tmp[0];
                ppotUp[3] = 0.0;
            } else {
                ppot[2] = tmp[0];
                ppot[3] = 0.0;
                ppot[0] = mod1d->M_BL.RD[0][0] * tmp[0];
                ppot[1] = 0.0;
                ppotUp[0] = tmp[1];
                ppotUp[1] = tmp[2];
                grt_cmat2x1_mul(RU_FA_delay, ppotUp, ppotUp + 2);
            }
        }
    }

    return;
}


void grt_secular_function_potential_Love(
    GRT_MODEL1D *mod1d, const real_t cphase, const size_t iref, cplx_t *psec, cplx_t ppot[GRT_LOVE_DIM])
{
    real_t k = creal(mod1d->omega)/cphase;
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
    GRT_MODEL1D *mod1d, const real_t cphase, const size_t iref, const DISPER_TYPE wtype, cplx_t *psec, cplx_t *ppot, cplx_t *ppotUp)
{
    mod1d->neval++;
    if(wtype == GRT_DISPERSION_RAYL){
        grt_secular_function_potential_Rayl(mod1d, cphase, iref, psec, ppot, ppotUp);
    } 
    else if(wtype == GRT_DISPERSION_LOVE){
        grt_secular_function_potential_Love(mod1d, cphase, iref, psec, ppot);
    }
    else {
        GRTRaiseError("Wrong execution.");
    }
}

void grt_secular_function(
    GRT_MODEL1D *mod1d, const real_t cphase, const size_t iref, const DISPER_TYPE wtype, cplx_t *psec)
{
    grt_secular_function_potential(mod1d, cphase, iref, wtype, psec, NULL, NULL);
}
