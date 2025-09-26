/**
 * @file   propagate.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码通过递推公式实现 广义反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient 
 *               matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699
 * 
 */

#include <stdio.h>
#include <complex.h>
#include <string.h>

#include "grt/dynamic/propagate.h"
#include "grt/dynamic/layer.h"
#include "grt/dynamic/source.h"
#include "grt/common/recursion.h"
#include "grt/common/model.h"
#include "grt/common/matrix.h"
#include "grt/common/prtdbg.h"


void grt_kernel(
    const GRT_MODEL1D *mod1d, MYCOMPLEX omega, MYREAL k, MYCOMPLEX QWV[GRT_SRC_M_NUM][GRT_QWV_NUM],
    bool calc_uiz,
    MYCOMPLEX QWV_uiz[GRT_SRC_M_NUM][GRT_QWV_NUM], MYINT *stats)
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
    // bool ircvup = true;
    // MYINT isrc = 2;
    // MYINT ircv = 1;
    // MYINT imin=1, imax=2;
    

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
    MYCOMPLEX RD[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL = 0.0;
    MYCOMPLEX RU[2][2] = GRT_INIT_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL = 0.0;
    MYCOMPLEX TD[2][2] = GRT_INIT_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TDL = 1.0;
    MYCOMPLEX TU[2][2] = GRT_INIT_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TUL = 1.0;

    // 自由表面的反射系数
    MYCOMPLEX R_tilt[2][2] = GRT_INIT_ZERO_2x2_MATRIX; // SH波在自由表面的反射系数为1，不必定义变量

    // 接收点处的接收矩阵(转为位移u的(B_m, C_m, P_m)系分量)
    MYCOMPLEX R_EV[2][2], R_EVL;

    // 接收点处的接收矩阵(转为位移导数ui_z的(B_m, C_m, P_m)系分量)
    MYCOMPLEX uiz_R_EV[2][2], uiz_R_EVL;

    // 模型参数
    // 后缀0，1分别代表上层和下层
    MYREAL thk0, thk1, Rho0, Rho1;
    MYCOMPLEX mu0, mu1;
    MYCOMPLEX xa0=0.0, xb0=0.0, xa1=0.0, xb1=0.0;
    MYCOMPLEX top_xa=0.0, top_xb=0.0;
    MYCOMPLEX rcv_xa=0.0, rcv_xb=0.0;
    MYCOMPLEX src_xa=0.0, src_xb=0.0;
    MYCOMPLEX cbcb0=0.0, cbcb1=0.0;
    MYCOMPLEX caca1=0.0;
    MYCOMPLEX src_caca=0.0, src_cbcb=0.0;
    MYCOMPLEX top_cbcb=0.0;

    // 相速度
    MYCOMPLEX c_phase = omega / k;

    // 从顶到底进行矩阵递推, 公式(5.5.3)
    for(MYINT iy=0; iy<mod1d->n; ++iy){ // 因为n>=3, 故一定会进入该循环

        // 赋值上层 
        thk0 = thk1;
        Rho0 = Rho1;
        mu0 = mu1;
        xa0 = xa1;
        xb0 = xb1;
        cbcb0 = cbcb1;

        // 更新模型参数
        thk1 = mod1d->Thk[iy];
        Rho1 = mod1d->Rho[iy];
        mu1 = mod1d->mu[iy];
        grt_get_mod1d_xa_xb(mod1d, iy, c_phase, &caca1, &xa1, &cbcb1, &xb1);

        if(0==iy){
            top_xa = xa1;
            top_xb = xb1;
            top_cbcb = cbcb1;
            continue;
        }

        // 确定上下层的物性参数
        if(ircv==iy){
            rcv_xa = xa1;
            rcv_xb = xb1;
        } else if(isrc==iy){
            src_xa = xa1;
            src_xb = xb1;
            src_caca = caca1;
            src_cbcb = cbcb1;
        } else {
            // 对第iy层的系数矩阵赋值，加入时间延迟因子(第iy-1界面与第iy界面之间)
            grt_RT_matrix_PSV(
                Rho0, xa0, xb0, cbcb0, mu0, 
                Rho1, xa1, xb1, cbcb1, mu1, 
                omega, k, 
                RD, RU, TD, TU, stats);
            grt_RT_matrix_SH(
                xb0, mu0, 
                xb1, mu1, 
                omega, k, 
                &RDL, &RUL, &TDL, &TUL);
            if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
            grt_delay_RT_matrix_PSV(xa0, xb0, thk0, k, RD, RU, TD, TU);
            grt_delay_RT_matrix_SH(xb0, thk0, k, &RDL, &RUL, &TDL, &TUL);
        }

#if Print_GRTCOEF == 1
        // TEST-------------------------------------------------------------
        // fprintf(stderr, "k=%f. iy=%d\n", k, iy);
        // fprintf(stderr, "RD\n");
        // cmatmxn_print(2, 2, RD);
        // fprintf(stderr, "RDL="GRT_CMPLX_FMT"\n", creal(RDL), cimag(RDL));
        // fprintf(stderr, "RU\n");
        // cmatmxn_print(2, 2, RU);
        // fprintf(stderr, "RUL="GRT_CMPLX_FMT"\n", creal(RUL), cimag(RUL));
        // fprintf(stderr, "TD\n");
        // cmatmxn_print(2, 2, TD);
        // fprintf(stderr, "TDL="GRT_CMPLX_FMT"\n", creal(TDL), cimag(TDL));
        // fprintf(stderr, "TU\n");
        // cmatmxn_print(2, 2, TU);
        // fprintf(stderr, "TUL="GRT_CMPLX_FMT"\n", creal(TUL), cimag(TUL));
        // if(creal(omega)==PI2*15e-4 && iy==5){
        // fprintf(stderr, GRT_REAL_FMT, k);
        // for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(RD[i][j]), cimag(RD[i][j]));
        // for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(RU[i][j]), cimag(RU[i][j]));
        // for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(TD[i][j]), cimag(TD[i][j]));
        // for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(TU[i][j]), cimag(TU[i][j]));
        // fprintf(stderr, "\n");
        // }
        // TEST-------------------------------------------------------------
#endif
        // FA
        if(iy < imin){ 
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
                if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
            }
        } 
        else if(iy==imin){ // 虚拟层位，可对递推公式简化
            grt_recursion_RT_matrix_virtual(
                xa0, xb0, thk0, k,
                RU_FA, &RUL_FA, 
                TD_FA, &TDL_FA, TU_FA, &TUL_FA);
        }
        // RS
        else if(iy < imax){
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
                if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
            }
        } 
        else if(iy==imax){ // 虚拟层位，可对递推公式简化
            grt_recursion_RT_matrix_virtual(
                xa0, xb0, thk0, k,
                RU_RS, &RUL_RS, 
                TD_RS, &TDL_RS, TU_RS, &TUL_RS);
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
                if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
            }
        } // END if


    } // END for loop 
    //===================================================================================

    // return;


    // 计算震源系数
    MYCOMPLEX src_coef_PSV[GRT_SRC_M_NUM][GRT_QWV_NUM-1][2] = {0};
    MYCOMPLEX src_coef_SH[GRT_SRC_M_NUM][2] = {0};
    grt_source_coef_PSV(src_xa, src_xb, src_caca, src_cbcb, k, src_coef_PSV);
    grt_source_coef_SH(src_xb, src_cbcb, k, src_coef_SH);

    // 临时中转矩阵 (temperary)
    MYCOMPLEX tmpR2[2][2], tmp2x2[2][2], tmpRL, tmp2x2_uiz[2][2], tmpRL2;
    MYCOMPLEX inv_2x2T[2][2], invT;

    // 递推RU_FA
    grt_topfree_RU_PSV(top_xa, top_xb, top_cbcb, k, R_tilt, stats);
    if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
    grt_recursion_RU(
        R_tilt, 1.0, 
        RD_FA, RDL_FA,
        RU_FA, RUL_FA, 
        TD_FA, TDL_FA,
        TU_FA, TUL_FA,
        RU_FA, &RUL_FA, NULL, NULL, stats);
    if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
    
    // 根据震源和台站相对位置，计算最终的系数
    if(ircvup){ // A接收  B震源

        // 计算R_EV
        grt_wave2qwv_REV_PSV(rcv_xa, rcv_xb, ircvup, k, RU_FA, R_EV);
        grt_wave2qwv_REV_SH(rcv_xb, k, RUL_FA, &R_EVL);

        // 递推RU_FS
        grt_recursion_RU(
            RU_FA, RUL_FA, // 已从ZR变为FR，加入了自由表面的效应
            RD_RS, RDL_RS,
            RU_RS, RUL_RS, 
            TD_RS, TDL_RS,
            TU_RS, TUL_RS,
            RU_FB, &RUL_FB, inv_2x2T, &invT, stats);
        if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
        
#if Print_GRTCOEF == 1
        // TEST-------------------------------------------------------------
        if(creal(omega)==PI2*0.1){
        fprintf(stderr, GRT_REAL_FMT, k);
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(RD_BL[i][j]), cimag(RD_BL[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(RU_BL[i][j]), cimag(RU_BL[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(TD_BL[i][j]), cimag(TD_BL[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(TU_BL[i][j]), cimag(TU_BL[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(RD_RS[i][j]), cimag(RD_RS[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(RU_RS[i][j]), cimag(RU_RS[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(TD_RS[i][j]), cimag(TD_RS[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(TU_RS[i][j]), cimag(TU_RS[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(RD_FA[i][j]), cimag(RD_FA[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(RU_FA[i][j]), cimag(RU_FA[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(TD_FA[i][j]), cimag(TD_FA[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(TU_FA[i][j]), cimag(TU_FA[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(RU_FB[i][j]), cimag(RU_FB[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(R_tilt[i][j]), cimag(R_tilt[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(R_EV[i][j]), cimag(R_EV[i][j]));
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(inv_2x2T[i][j]), cimag(inv_2x2T[i][j]));
        cmat2x2_mul(RD_BL, RU_FB, tmpR2);
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(tmpR2[i][j]), cimag(tmpR2[i][j]));
        cmat2x2_one_sub(tmpR2);
        cmat2x2_inv(tmpR2, tmpR2, stats);// (I - xx)^-1
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(tmpR2[i][j]), cimag(tmpR2[i][j]));
        cmat2x2_mul(inv_2x2T, tmpR2, tmp2x2);
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(tmp2x2[i][j]), cimag(tmp2x2[i][j]));
        cmat2x2_mul(R_EV, tmp2x2, tmp2x2);
        for(int i=0; i<2; ++i)  for(int j=0; j<2; ++j)  fprintf(stderr, GRT_CMPLX_FMT, creal(tmp2x2[i][j]), cimag(tmp2x2[i][j]));
        fprintf(stderr, "\n");
        }
        // TEST-------------------------------------------------------------
#endif

        // 公式(5.7.12-14)
        grt_cmat2x2_mul(RD_BL, RU_FB, tmpR2);
        grt_cmat2x2_one_sub(tmpR2);
        grt_cmat2x2_inv(tmpR2, tmpR2, stats);// (I - xx)^-1
        if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
        grt_cmat2x2_mul(inv_2x2T, tmpR2, tmp2x2);

        if(calc_uiz) grt_cmat2x2_assign(tmp2x2, tmp2x2_uiz); // 为后续计算空间导数备份

        grt_cmat2x2_mul(R_EV, tmp2x2, tmp2x2);
        tmpRL = invT  / (1.0 - RDL_BL * RUL_FB);
        tmpRL2 = R_EVL * tmpRL;

        for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
            grt_construct_qwv(ircvup, tmp2x2, tmpRL2, RD_BL, RDL_BL, src_coef_PSV[i], src_coef_SH[i], QWV[i]);
        }


        if(calc_uiz){
            grt_wave2qwv_z_REV_PSV(rcv_xa, rcv_xb, ircvup, k, RU_FA, uiz_R_EV);
            grt_wave2qwv_z_REV_SH(rcv_xb, ircvup, k, RUL_FA, &uiz_R_EVL);
            grt_cmat2x2_mul(uiz_R_EV, tmp2x2_uiz, tmp2x2_uiz);
            tmpRL2 = uiz_R_EVL * tmpRL;

            for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
                grt_construct_qwv(ircvup, tmp2x2_uiz, tmpRL2, RD_BL, RDL_BL, src_coef_PSV[i], src_coef_SH[i], QWV_uiz[i]);
            }    
        }
    } 
    else { // A震源  B接收

        // 计算R_EV
        grt_wave2qwv_REV_PSV(rcv_xa, rcv_xb, ircvup, k, RD_BL, R_EV);    
        grt_wave2qwv_REV_SH(rcv_xb, k, RDL_BL, &R_EVL);    

        // 递推RD_SL
        grt_recursion_RD(
            RD_RS, RDL_RS,
            RU_RS, RUL_RS,
            TD_RS, TDL_RS,
            TU_RS, TUL_RS,
            RD_BL, RDL_BL,
            RD_AL, &RDL_AL, inv_2x2T, &invT, stats);
        if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
        
        // 公式(5.7.26-27)
        grt_cmat2x2_mul(RU_FA, RD_AL, tmpR2);
        grt_cmat2x2_one_sub(tmpR2);
        grt_cmat2x2_inv(tmpR2, tmpR2, stats);// (I - xx)^-1
        if(*stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
        grt_cmat2x2_mul(inv_2x2T, tmpR2, tmp2x2);

        if(calc_uiz) grt_cmat2x2_assign(tmp2x2, tmp2x2_uiz); // 为后续计算空间导数备份

        grt_cmat2x2_mul(R_EV, tmp2x2, tmp2x2);
        tmpRL = invT / (1.0 - RUL_FA * RDL_AL);
        tmpRL2 = R_EVL * tmpRL;

        for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
            grt_construct_qwv(ircvup, tmp2x2, tmpRL2, RU_FA, RUL_FA, src_coef_PSV[i], src_coef_SH[i], QWV[i]);
        }


        if(calc_uiz){
            grt_wave2qwv_z_REV_PSV(rcv_xa, rcv_xb, ircvup, k, RD_BL, uiz_R_EV);    
            grt_wave2qwv_z_REV_SH(rcv_xb, ircvup, k, RDL_BL, &uiz_R_EVL);    
            grt_cmat2x2_mul(uiz_R_EV, tmp2x2_uiz, tmp2x2_uiz);
            tmpRL2 = uiz_R_EVL * tmpRL;
            
            for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
                grt_construct_qwv(ircvup, tmp2x2_uiz, tmpRL2, RU_FA, RUL_FA, src_coef_PSV[i], src_coef_SH[i], QWV_uiz[i]);
            }
        }

    } // END if



    BEFORE_RETURN:

    // 对一些特殊情况的修正
    // 当震源和场点均位于地表时，可理论验证DS分量恒为0，这里直接赋0以避免后续的精度干扰
    if(mod1d->depsrc == 0.0 && mod1d->deprcv == 0.0)
    {
        for(MYINT c=0; c<GRT_QWV_NUM; ++c){
            QWV[GRT_SRC_M_DS_INDEX][c] = 0.0;
            if(calc_uiz)  QWV_uiz[GRT_SRC_M_DS_INDEX][c] = 0.0;
        }
    }

}






void grt_construct_qwv(
    bool ircvup, 
    const MYCOMPLEX R1[2][2], MYCOMPLEX RL1, 
    const MYCOMPLEX R2[2][2], MYCOMPLEX RL2, 
    const MYCOMPLEX coef_PSV[GRT_QWV_NUM-1][2], const MYCOMPLEX coef_SH[2], 
    MYCOMPLEX qwv[GRT_QWV_NUM])
{
    MYCOMPLEX qw0[2], qw1[2], v0;
    MYCOMPLEX coefD[2] = {coef_PSV[0][0], coef_PSV[1][0]};
    MYCOMPLEX coefU[2] = {coef_PSV[0][1], coef_PSV[1][1]};
    if(ircvup){
        grt_cmat2x1_mul(R2, coefD, qw0);
        qw0[0] += coefU[0]; qw0[1] += coefU[1]; 
        v0 = RL1 * (RL2*coef_SH[0] + coef_SH[1]);
    } else {
        grt_cmat2x1_mul(R2, coefU, qw0);
        qw0[0] += coefD[0]; qw0[1] += coefD[1]; 
        v0 = RL1 * (coef_SH[0] + RL2*coef_SH[1]);
    }
    grt_cmat2x1_mul(R1, qw0, qw1);

    qwv[0] = qw1[0];
    qwv[1] = qw1[1];
    qwv[2] = v0;
}