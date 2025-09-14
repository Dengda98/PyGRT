/**
 * @file   static_propagate.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 * 以下代码实现的是 静态广义反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *              透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 *
 */


#include <stdio.h>
#include <complex.h>

#include "grt/static/static_propagate.h"
#include "grt/static/static_layer.h"
#include "grt/static/static_source.h"
#include "grt/common/recursion.h"
#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/common/matrix.h"


void grt_static_kernel(
    const GRT_MODEL1D *mod1d, MYCOMPLEX omega, MYREAL k, MYCOMPLEX QWV[GRT_SRC_M_NUM][GRT_QWV_NUM],
    bool calc_uiz, MYCOMPLEX QWV_uiz[GRT_SRC_M_NUM][GRT_QWV_NUM], MYINT *stats)
{   
    // 初始化qwv为0
    for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
        for(MYINT j=0; j<GRT_QWV_NUM; ++j){
            QWV[i][j] = 0.0;
            if(calc_uiz)  QWV_uiz[i][j] = 0.0;
        }
    }


    bool ircvup = mod1d->ircvup;
    MYINT isrc = mod1d->isrc; // 震源所在虚拟层位, isrc>=1
    MYINT ircv = mod1d->ircv; // 接收点所在虚拟层位, ircv>=1, ircv != isrc
    MYINT imin, imax; // 相对浅层深层层位
    imin = GRT_MIN(mod1d->isrc, mod1d->ircv);
    imax = GRT_MAX(mod1d->isrc, mod1d->ircv);

    // 初始化广义反射透射系数矩阵
    // BL
    MYCOMPLEX RD_BL[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL_BL = 0.0;
    MYCOMPLEX RU_BL[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL_BL = 0.0;
    MYCOMPLEX TD_BL[2][2] = GRT_INIT_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TDL_BL = 1.0;
    MYCOMPLEX TU_BL[2][2] = GRT_INIT_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TUL_BL = 1.0;
    // AL
    MYCOMPLEX RD_AL[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL_AL = 0.0;
    // RS
    MYCOMPLEX RD_RS[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL_RS = 0.0;
    MYCOMPLEX RU_RS[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL_RS = 0.0;
    MYCOMPLEX TD_RS[2][2] = GRT_INIT_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TDL_RS = 1.0;
    MYCOMPLEX TU_RS[2][2] = GRT_INIT_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TUL_RS = 1.0;
    // FA (实际先计算ZA，再递推到FA)
    MYCOMPLEX RD_FA[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL_FA = 0.0;
    MYCOMPLEX RU_FA[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL_FA = 0.0;
    MYCOMPLEX TD_FA[2][2] = GRT_INIT_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TDL_FA = 1.0;
    MYCOMPLEX TU_FA[2][2] = GRT_INIT_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TUL_FA = 1.0;
    // FB (实际先计算ZB，再递推到FB)
    MYCOMPLEX RU_FB[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL_FB = 0.0;

    // 定义物理层内的反射透射系数矩阵，相对于界面上的系数矩阵增加了时间延迟因子
    MYCOMPLEX RD[2][2], RDL, TD[2][2], TDL;
    MYCOMPLEX RU[2][2], RUL, TU[2][2], TUL;

    // 自由表面的反射系数
    MYCOMPLEX R_tilt[2][2] = GRT_INIT_ZERO_2x2_MATRIX; // SH波在自由表面的反射系数为1，不必定义变量

    // 接收点处的接收矩阵
    MYCOMPLEX R_EV[2][2], R_EVL;
    
    // 接收点处的接收矩阵(转为位移导数ui_z的(B_m, C_m, P_m)系分量)
    MYCOMPLEX uiz_R_EV[2][2], uiz_R_EVL;

    // 模型参数
    // 后缀0，1分别代表上层和下层
    MYREAL thk0, thk1;
    MYCOMPLEX mu0, mu1;
    MYCOMPLEX delta0, delta1;
    MYCOMPLEX top_delta = 0.0;
    MYCOMPLEX src_delta = 0.0;
    MYCOMPLEX rcv_delta = 0.0;
    

    // 从顶到底进行矩阵递推, 公式(5.5.3)
    for(MYINT iy=0; iy<mod1d->n; ++iy){ // 因为n>=3, 故一定会进入该循环
        // 赋值上层 
        thk0 = thk1;
        mu0 = mu1;
        delta0 = delta1;

        // 更新模型参数
        thk1 = mod1d->Thk[iy];
        mu1 = mod1d->mu[iy];
        delta1 = mod1d->delta[iy];

        if(0==iy){
            top_delta = delta1;
            continue;
        }

        // 确定上下层的物性参数
        if(ircv==iy){
            rcv_delta = delta1;
        } else if(isrc==iy){
            src_delta = delta1;
        }

        // 这里和动态解情况不同，即使是震源层、接收层这种虚拟层位也需要显式计算R/T矩阵
        // 对第iy层的系数矩阵赋值
        grt_static_RT_matrix_PSV(
            delta0, mu0,
            delta1, mu1,
            thk0, k, // 使用iy-1层的厚度
            RD, RU, TD, TU);
        grt_static_RT_matrix_SH(
            mu0, mu1,
            thk0, k, // 使用iy-1层的厚度
            &RDL, &RUL, &TDL, &TUL);
        // 加入时间延迟因子(第iy-1界面与第iy界面之间)
        grt_static_delay_RT_matrix_PSV(thk0, k, RD, RU, TD, TU);
        grt_static_delay_RT_matrix_SH(thk0, k, &RDL, &RUL, &TDL, &TUL);

        // FA
        if(iy <= imin){
            if(iy == 1){ // 初始化FA
                GRT_RT_PSV_ASSIGN(FA);
                GRT_RT_SH_ASSIGN(FA);
            } else { // 递推FA
                grt_recursion_RT_matrix(
                    RD_FA, RDL_FA, RU_FA, RUL_FA, 
                    TD_FA, TDL_FA, TU_FA, TUL_FA,
                    RD, RDL, RU, RUL, 
                    TD, TDL, TU, TUL,
                    RD_FA, &RDL_FA, RU_FA, &RUL_FA, 
                    TD_FA, &TDL_FA, TU_FA, &TUL_FA, stats);  
            }
        }
        // RS
        else if(iy <= imax){
            if(iy == imin+1){// 初始化RS
                GRT_RT_PSV_ASSIGN(RS);
                GRT_RT_SH_ASSIGN(RS);
            } else { // 递推RS
                grt_recursion_RT_matrix(
                    RD_RS, RDL_RS, RU_RS, RUL_RS, 
                    TD_RS, TDL_RS, TU_RS, TUL_RS,
                    RD, RDL, RU, RUL, 
                    TD, TDL, TU, TUL,
                    RD_RS, &RDL_RS, RU_RS, &RUL_RS, 
                    TD_RS, &TDL_RS, TU_RS, &TUL_RS, stats);  // 写入原地址
            }
        } 
        // BL
        else {
            if(iy == imax+1){// 初始化BL
                GRT_RT_PSV_ASSIGN(BL);
                GRT_RT_SH_ASSIGN(BL);
            } else { // 递推BL
                // 只有 RD 矩阵最终会被使用到
                grt_recursion_RT_matrix(
                    RD_BL, RDL_BL, RU_BL, RUL_BL, 
                    TD_BL, TDL_BL, TU_BL, TUL_BL,
                    RD, RDL, RU, RUL, 
                    TD, TDL, TU, TUL,
                    RD_BL, &RDL_BL, RU_BL, &RUL_BL, 
                    TD_BL, &TDL_BL, TU_BL, &TUL_BL, stats);  // 写入原地址
            }

        } // END if

    } // END for loop
    //===================================================================================


    // 计算震源系数
    MYCOMPLEX src_coef_PSV[GRT_SRC_M_NUM][GRT_QWV_NUM-1][2] = {0};
    MYCOMPLEX src_coef_SH[GRT_SRC_M_NUM][2] = {0};
    grt_static_source_coef_PSV(src_delta, k, src_coef_PSV);
    grt_static_source_coef_SH(k, src_coef_SH);

    // 临时中转矩阵 (temperary)
    MYCOMPLEX tmpR2[2][2], tmp2x2[2][2], tmpRL, tmp2x2_uiz[2][2], tmpRL_uiz;
    MYCOMPLEX inv_2x2T[2][2], invT;

    // 递推RU_FA
    grt_static_topfree_RU_PSV(top_delta, R_tilt);
    grt_recursion_RU(
        R_tilt, 1.0, 
        RD_FA, RDL_FA,
        RU_FA, RUL_FA, 
        TD_FA, TDL_FA,
        TU_FA, TUL_FA,
        RU_FA, &RUL_FA, NULL, NULL, stats);

    // 根据震源和台站相对位置，计算最终的系数
    if(ircvup){ // A接收  B震源
        // 计算R_EV
        grt_static_wave2qwv_REV_PSV(ircvup, RU_FA, R_EV);
        grt_static_wave2qwv_REV_SH(RUL_FA, &R_EVL);

        // 递推RU_FS
        grt_recursion_RU(
            RU_FA, RUL_FA, // 已从ZR变为FR，加入了自由表面的效应
            RD_RS, RDL_RS,
            RU_RS, RUL_RS, 
            TD_RS, TDL_RS,
            TU_RS, TUL_RS,
            RU_FB, &RUL_FB, inv_2x2T, &invT, stats);
        
        // 公式(5.7.12-14)
        grt_cmat2x2_mul(RD_BL, RU_FB, tmpR2);
        grt_cmat2x2_one_sub(tmpR2);
        grt_cmat2x2_inv(tmpR2, tmpR2, stats);// (I - xx)^-1
        grt_cmat2x2_mul(inv_2x2T, tmpR2, tmp2x2);

        if(calc_uiz) grt_cmat2x2_assign(tmp2x2, tmp2x2_uiz); // 为后续计算空间导数备份

        grt_cmat2x2_mul(R_EV, tmp2x2, tmp2x2);
        tmpRL = R_EVL * invT  / (1.0 - RDL_BL * RUL_FB);

        for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
            grt_construct_qwv(ircvup, tmp2x2, tmpRL, RD_BL, RDL_BL, src_coef_PSV[i], src_coef_SH[i], QWV[i]);
        }
        

        if(calc_uiz){
            grt_static_wave2qwv_z_REV_PSV(rcv_delta, ircvup, k, RU_FA, uiz_R_EV);
            grt_static_wave2qwv_z_REV_SH(ircvup, k, RUL_FA, &uiz_R_EVL);
            grt_cmat2x2_mul(uiz_R_EV, tmp2x2_uiz, tmp2x2_uiz);
            tmpRL_uiz = tmpRL / R_EVL * uiz_R_EVL;
            
            for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
                grt_construct_qwv(ircvup, tmp2x2_uiz, tmpRL_uiz, RD_BL, RDL_BL, src_coef_PSV[i], src_coef_SH[i], QWV_uiz[i]);
            }    
        }
    }
    else { // A震源  B接收

        // 计算R_EV
        grt_static_wave2qwv_REV_PSV(ircvup, RD_BL, R_EV);    
        grt_static_wave2qwv_REV_SH(RDL_BL, &R_EVL);    

        // 递推RD_SL
        grt_recursion_RD(
            RD_RS, RDL_RS,
            RU_RS, RUL_RS,
            TD_RS, TDL_RS,
            TU_RS, TUL_RS,
            RD_BL, RDL_BL,
            RD_AL, &RDL_AL, inv_2x2T, &invT, stats);
        
        // 公式(5.7.26-27)
        grt_cmat2x2_mul(RU_FA, RD_AL, tmpR2);
        grt_cmat2x2_one_sub(tmpR2);
        grt_cmat2x2_inv(tmpR2, tmpR2, stats);// (I - xx)^-1
        grt_cmat2x2_mul(inv_2x2T, tmpR2, tmp2x2);
        
        if(calc_uiz) grt_cmat2x2_assign(tmp2x2, tmp2x2_uiz); // 为后续计算空间导数备份

        grt_cmat2x2_mul(R_EV, tmp2x2, tmp2x2);
        tmpRL = R_EVL * invT / (1.0 - RUL_FA * RDL_AL);
        
        for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
            grt_construct_qwv(ircvup, tmp2x2, tmpRL, RU_FA, RUL_FA, src_coef_PSV[i], src_coef_SH[i], QWV[i]);
        }

        if(calc_uiz){
            grt_static_wave2qwv_z_REV_PSV(rcv_delta, ircvup, k, RD_BL, uiz_R_EV);    
            grt_static_wave2qwv_z_REV_SH(ircvup, k, RDL_BL, &uiz_R_EVL);    
            grt_cmat2x2_mul(uiz_R_EV, tmp2x2_uiz, tmp2x2_uiz);
            tmpRL_uiz = tmpRL / R_EVL * uiz_R_EVL;
            
            for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
                grt_construct_qwv(ircvup, tmp2x2_uiz, tmpRL_uiz, RU_FA, RUL_FA, src_coef_PSV[i], src_coef_SH[i], QWV_uiz[i]);
            }
        }
    } // END if
}
