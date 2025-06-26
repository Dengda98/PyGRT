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

#include "dynamic/propagate.h"
#include "dynamic/layer.h"
#include "dynamic/source.h"
#include "common/recursion.h"
#include "common/model.h"
#include "common/matrix.h"
#include "common/prtdbg.h"

#define CMAT_ASSIGN_SPLIT 0  // 2x2的小矩阵赋值合并为1个循环，程序速度提升微小




void kernel(
    const MODEL1D *mod1d, MYCOMPLEX omega, MYREAL k, MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM],
    bool calc_uiz,
    MYCOMPLEX QWV_uiz[SRC_M_NUM][QWV_NUM], MYINT *stats)
{
    // 初始化qwv为0
    for(MYINT i=0; i<SRC_M_NUM; ++i){
        for(MYINT j=0; j<QWV_NUM; ++j){
            QWV[i][j] = CZERO;
            if(calc_uiz)  QWV_uiz[i][j] = CZERO;
        }
    }

    bool ircvup = mod1d->ircvup;
    MYINT isrc = mod1d->isrc; // 震源所在虚拟层位, isrc>=1
    MYINT ircv = mod1d->ircv; // 接收点所在虚拟层位, ircv>=1, ircv != isrc
    MYINT imin, imax; // 相对浅层深层层位
    imin = mod1d->imin;
    imax = mod1d->imax;
    // bool ircvup = true;
    // MYINT isrc = 2;
    // MYINT ircv = 1;
    // MYINT imin=1, imax=2;
    

    // 初始化广义反射透射系数矩阵
    // BL
    MYCOMPLEX RD_BL[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL_BL = CZERO;
    MYCOMPLEX RU_BL[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL_BL = CZERO;
    MYCOMPLEX TD_BL[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TDL_BL = CONE;
    MYCOMPLEX TU_BL[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TUL_BL = CONE;
    // AL
    MYCOMPLEX RD_AL[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL_AL = CZERO;
    // RS
    MYCOMPLEX RD_RS[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL_RS = CZERO;
    MYCOMPLEX RU_RS[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL_RS = CZERO;
    MYCOMPLEX TD_RS[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TDL_RS = CONE;
    MYCOMPLEX TU_RS[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TUL_RS = CONE;
    // FA (实际先计算ZA，再递推到FA)
    MYCOMPLEX RD_FA[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL_FA = CZERO;
    MYCOMPLEX RU_FA[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL_FA = CZERO;
    MYCOMPLEX TD_FA[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TDL_FA = CONE;
    MYCOMPLEX TU_FA[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TUL_FA = CONE;
    // FB (实际先计算ZB，再递推到FB)
    MYCOMPLEX RU_FB[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL_FB = CZERO;

    // 抽象指针 
    // BL
    MYCOMPLEX *const pRDL_BL = &RDL_BL;
    MYCOMPLEX *const pRUL_BL = &RUL_BL;
    MYCOMPLEX *const pTDL_BL = &TDL_BL;
    MYCOMPLEX *const pTUL_BL = &TUL_BL;
    // AL
    MYCOMPLEX *const pRDL_AL = &RDL_AL;
    // RS
    MYCOMPLEX *const pRDL_RS = &RDL_RS;
    MYCOMPLEX *const pRUL_RS = &RUL_RS;
    MYCOMPLEX *const pTDL_RS = &TDL_RS;
    MYCOMPLEX *const pTUL_RS = &TUL_RS;
    // FA
    MYCOMPLEX *const pRDL_FA = &RDL_FA;
    MYCOMPLEX *const pRUL_FA = &RUL_FA;
    MYCOMPLEX *const pTDL_FA = &TDL_FA;
    MYCOMPLEX *const pTUL_FA = &TUL_FA;
    // FB 
    MYCOMPLEX *const pRUL_FB = &RUL_FB;

    
    // 定义物理层内的反射透射系数矩阵，相对于界面上的系数矩阵增加了时间延迟因子
    MYCOMPLEX RD[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL = CZERO;
    MYCOMPLEX RU[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL = CZERO;
    MYCOMPLEX TD[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TDL = CONE;
    MYCOMPLEX TU[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TUL = CONE;
    MYCOMPLEX *const pRDL = &RDL;
    MYCOMPLEX *const pTDL = &TDL;
    MYCOMPLEX *const pRUL = &RUL;
    MYCOMPLEX *const pTUL = &TUL;


    // 自由表面的反射系数
    MYCOMPLEX R_tilt[2][2] = INIT_C_ZERO_2x2_MATRIX; // SH波在自由表面的反射系数为1，不必定义变量

    // 接收点处的接收矩阵(转为位移u的(B_m, C_m, P_m)系分量)
    MYCOMPLEX R_EV[2][2], R_EVL;
    MYCOMPLEX *const pR_EVL = &R_EVL;

    // 接收点处的接收矩阵(转为位移导数ui_z的(B_m, C_m, P_m)系分量)
    MYCOMPLEX uiz_R_EV[2][2], uiz_R_EVL;
    MYCOMPLEX *const puiz_R_EVL = &uiz_R_EVL;
    

    // 模型参数
    // 后缀0，1分别代表上层和下层
    LAYER *lay = NULL;
    MYREAL mod1d_thk0, mod1d_thk1, mod1d_Rho0, mod1d_Rho1;
    MYCOMPLEX mod1d_mu0, mod1d_mu1;
    MYCOMPLEX mod1d_kaka1, mod1d_kbkb0, mod1d_kbkb1;
    MYCOMPLEX mod1d_xa0, mod1d_xb0, mod1d_xa1, mod1d_xb1;
    MYCOMPLEX top_xa=RZERO, top_xb=RZERO, top_kbkb=RZERO;
    MYCOMPLEX rcv_xa=RZERO, rcv_xb=RZERO;
    MYCOMPLEX src_xa=RZERO, src_xb=RZERO, src_kaka=RZERO, src_kbkb=RZERO;


    // 从顶到底进行矩阵递推, 公式(5.5.3)
    for(MYINT iy=0; iy<mod1d->n; ++iy){ // 因为n>=3, 故一定会进入该循环
        lay = mod1d->lays + iy;

        // 赋值上层 
        mod1d_thk0 = mod1d_thk1;
        mod1d_Rho0 = mod1d_Rho1;
        mod1d_mu0 = mod1d_mu1;
        mod1d_kbkb0 = mod1d_kbkb1;
        mod1d_xa0 = mod1d_xa1;
        mod1d_xb0 = mod1d_xb1;

        // 更新模型参数
        mod1d_thk1 = lay->thk;
        mod1d_Rho1 = lay->Rho;
        mod1d_mu1 = lay->mu;
        mod1d_kaka1 = lay->kaka;
        mod1d_kbkb1 = lay->kbkb;
        mod1d_xa1 = sqrt(RONE - mod1d_kaka1/(k*k));
        mod1d_xb1 = sqrt(RONE - mod1d_kbkb1/(k*k));
        // mod1d_xa1 = sqrt(k*k - mod1d_kaka1)/k;
        // mod1d_xb1 = sqrt(k*k - mod1d_kbkb1)/k;

        if(0==iy){
            top_xa = mod1d_xa1;
            top_xb = mod1d_xb1;
            top_kbkb = mod1d_kbkb1;
            continue;
        }

        // 确定上下层的物性参数
        if(ircv==iy){
            rcv_xa = mod1d_xa1;
            rcv_xb = mod1d_xb1;
        } else if(isrc==iy){
            src_xa = mod1d_xa1;
            src_xb = mod1d_xb1;
            src_kaka = mod1d_kaka1;
            src_kbkb = mod1d_kbkb1;
        } else {
            // 对第iy层的系数矩阵赋值，加入时间延迟因子(第iy-1界面与第iy界面之间)
            calc_RT_2x2(
                mod1d_Rho0, mod1d_xa0, mod1d_xb0, mod1d_kbkb0, mod1d_mu0, 
                mod1d_Rho1, mod1d_xa1, mod1d_xb1, mod1d_kbkb1, mod1d_mu1, 
                mod1d_thk0, // 使用iy-1层的厚度
                omega, k, 
                RD, pRDL, RU, pRUL, 
                TD, pTDL, TU, pTUL, stats);
            if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN;
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
#if CMAT_ASSIGN_SPLIT == 1
                cmat2x2_assign(RD, RD_FA);  RDL_FA = RDL;
                cmat2x2_assign(RU, RU_FA);  RUL_FA = RUL;
                cmat2x2_assign(TD, TD_FA);  TDL_FA = TDL;
                cmat2x2_assign(TU, TU_FA);  TUL_FA = TUL;
#else 
                for(MYINT kk=0; kk<2; ++kk){
                    for(MYINT pp=0; pp<2; ++pp){
                        RD_FA[kk][pp] = RD[kk][pp];
                        RU_FA[kk][pp] = RU[kk][pp];
                        TD_FA[kk][pp] = TD[kk][pp];
                        TU_FA[kk][pp] = TU[kk][pp];
                    }
                }
                RDL_FA = RDL;
                RUL_FA = RUL;
                TDL_FA = TDL;
                TUL_FA = TUL;
#endif
            } else { // 递推FA

                recursion_RT_2x2(
                    RD_FA, RDL_FA, RU_FA, RUL_FA, 
                    TD_FA, TDL_FA, TU_FA, TUL_FA,
                    RD, RDL, RU, RUL, 
                    TD, TDL, TU, TUL,
                    RD_FA, pRDL_FA, RU_FA, pRUL_FA, 
                    TD_FA, pTDL_FA, TU_FA, pTUL_FA, stats);  
                if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN;
            }
        } 
        else if(iy==imin){ // 虚拟层位，可对递推公式简化
            recursion_RT_2x2_imaginary(
                mod1d_xa0, mod1d_xb0, mod1d_thk0, k,
                RU_FA, pRUL_FA, 
                TD_FA, pTDL_FA, TU_FA, pTUL_FA);
        }
        // RS
        else if(iy < imax){
            if(iy == imin+1){// 初始化RS
#if CMAT_ASSIGN_SPLIT == 1
                cmat2x2_assign(RD, RD_RS);  RDL_RS = RDL;
                cmat2x2_assign(RU, RU_RS);  RUL_RS = RUL;
                cmat2x2_assign(TD, TD_RS);  TDL_RS = TDL;
                cmat2x2_assign(TU, TU_RS);  TUL_RS = TUL;
#else
                for(MYINT kk=0; kk<2; ++kk){
                    for(MYINT pp=0; pp<2; ++pp){
                        RD_RS[kk][pp] = RD[kk][pp];
                        RU_RS[kk][pp] = RU[kk][pp];
                        TD_RS[kk][pp] = TD[kk][pp];
                        TU_RS[kk][pp] = TU[kk][pp];
                    }
                }
                RDL_RS = RDL;
                RUL_RS = RUL;
                TDL_RS = TDL;
                TUL_RS = TUL;
#endif
            } else { // 递推RS
                recursion_RT_2x2(
                    RD_RS, RDL_RS, RU_RS, RUL_RS, 
                    TD_RS, TDL_RS, TU_RS, TUL_RS,
                    RD, RDL, RU, RUL, 
                    TD, TDL, TU, TUL,
                    RD_RS, pRDL_RS, RU_RS, pRUL_RS, 
                    TD_RS, pTDL_RS, TU_RS, pTUL_RS, stats);  // 写入原地址
                if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN;
            }
        } 
        else if(iy==imax){ // 虚拟层位，可对递推公式简化
            recursion_RT_2x2_imaginary(
                mod1d_xa0, mod1d_xb0, mod1d_thk0, k,
                RU_RS, pRUL_RS, 
                TD_RS, pTDL_RS, TU_RS, pTUL_RS);
        }
        // BL
        else {
            if(iy == imax+1){// 初始化BL
#if CMAT_ASSIGN_SPLIT == 1
                cmat2x2_assign(RD, RD_BL);  RDL_BL = RDL;
                cmat2x2_assign(RU, RU_BL);  RUL_BL = RUL;
                cmat2x2_assign(TD, TD_BL);  TDL_BL = TDL;
                cmat2x2_assign(TU, TU_BL);  TUL_BL = TUL;
#else 
                for(MYINT kk=0; kk<2; ++kk){
                    for(MYINT pp=0; pp<2; ++pp){
                        RD_BL[kk][pp] = RD[kk][pp];
                        RU_BL[kk][pp] = RU[kk][pp];
                        TD_BL[kk][pp] = TD[kk][pp];
                        TU_BL[kk][pp] = TU[kk][pp];
                    }
                }
                RDL_BL = RDL;
                RUL_BL = RUL;
                TDL_BL = TDL;
                TUL_BL = TUL;
#endif
            } else { // 递推BL

                // 这个IF纯粹是为了优化，因为不论是SL还是RL，只有RD矩阵最终会被使用到
                // 这里最终只把RD矩阵的值记录下来，其它的舍去，以减少部分运算
                if(iy < mod1d->n - 1){
                    recursion_RT_2x2(
                        RD_BL, RDL_BL, RU_BL, RUL_BL, 
                        TD_BL, TDL_BL, TU_BL, TUL_BL,
                        RD, RDL, RU, RUL, 
                        TD, TDL, TU, TUL,
                        RD_BL, pRDL_BL, RU_BL, pRUL_BL, 
                        TD_BL, pTDL_BL, TU_BL, pTUL_BL, stats);  // 写入原地址
                } else {
                    recursion_RT_2x2(
                        RD_BL, RDL_BL, RU_BL, RUL_BL, 
                        TD_BL, TDL_BL, TU_BL, TUL_BL,
                        RD, RDL, RU, RUL, 
                        TD, TDL, TU, TUL,
                        RD_BL, pRDL_BL, NULL, NULL, 
                        NULL, NULL, NULL, NULL, stats);  // 写入原地址
                }
                if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN;
            }
        } // END if


    } // END for loop 
    //===================================================================================

    // return;


    // 计算震源系数
    MYCOMPLEX src_coef[SRC_M_NUM][QWV_NUM][2] = {0};
    source_coef(src_xa, src_xb, src_kaka, src_kbkb, k, src_coef);

    // 临时中转矩阵 (temperary)
    MYCOMPLEX tmpR2[2][2], tmp2x2[2][2], tmpRL, tmp2x2_uiz[2][2], tmpRL2;
    MYCOMPLEX inv_2x2T[2][2], invT;

    // 递推RU_FA
    calc_R_tilt(top_xa, top_xb, top_kbkb, k, R_tilt, stats);
    if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN;
    recursion_RU(
        R_tilt, RONE, 
        RD_FA, RDL_FA,
        RU_FA, RUL_FA, 
        TD_FA, TDL_FA,
        TU_FA, TUL_FA,
        RU_FA, pRUL_FA, NULL, NULL, stats);
    if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN;
    
    // 根据震源和台站相对位置，计算最终的系数
    if(ircvup){ // A接收  B震源

        // 计算R_EV
        calc_R_EV(rcv_xa, rcv_xb, ircvup, k, RU_FA, RUL_FA, R_EV, pR_EVL);

        // 递推RU_FS
        recursion_RU(
            RU_FA, RUL_FA, // 已从ZR变为FR，加入了自由表面的效应
            RD_RS, RDL_RS,
            RU_RS, RUL_RS, 
            TD_RS, TDL_RS,
            TU_RS, TUL_RS,
            RU_FB, pRUL_FB, inv_2x2T, &invT, stats);
        if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN;
        
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
        cmat2x2_mul(RD_BL, RU_FB, tmpR2);
        cmat2x2_one_sub(tmpR2);
        cmat2x2_inv(tmpR2, tmpR2, stats);// (I - xx)^-1
        if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN;
        cmat2x2_mul(inv_2x2T, tmpR2, tmp2x2);

        if(calc_uiz) cmat2x2_assign(tmp2x2, tmp2x2_uiz); // 为后续计算空间导数备份

        cmat2x2_mul(R_EV, tmp2x2, tmp2x2);
        tmpRL = invT  / (RONE - RDL_BL * RUL_FB);
        tmpRL2 = R_EVL * tmpRL;

        for(MYINT i=0; i<SRC_M_NUM; ++i){
            get_qwv(ircvup, tmp2x2, tmpRL2, RD_BL, RDL_BL, src_coef[i], QWV[i]);
        }


        if(calc_uiz){
            calc_uiz_R_EV(rcv_xa, rcv_xb, ircvup, k, RU_FA, RUL_FA, uiz_R_EV, puiz_R_EVL);
            cmat2x2_mul(uiz_R_EV, tmp2x2_uiz, tmp2x2_uiz);
            tmpRL2 = uiz_R_EVL * tmpRL;

            for(MYINT i=0; i<SRC_M_NUM; ++i){
                get_qwv(ircvup, tmp2x2_uiz, tmpRL2, RD_BL, RDL_BL, src_coef[i], QWV_uiz[i]);
            }    
        }
    } 
    else { // A震源  B接收

        // 计算R_EV
        calc_R_EV(rcv_xa, rcv_xb, ircvup, k, RD_BL, RDL_BL, R_EV, pR_EVL);    

        // 递推RD_SL
        recursion_RD(
            RD_RS, RDL_RS,
            RU_RS, RUL_RS,
            TD_RS, TDL_RS,
            TU_RS, TUL_RS,
            RD_BL, RDL_BL,
            RD_AL, pRDL_AL, inv_2x2T, &invT, stats);
        if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN;
        
        // 公式(5.7.26-27)
        cmat2x2_mul(RU_FA, RD_AL, tmpR2);
        cmat2x2_one_sub(tmpR2);
        cmat2x2_inv(tmpR2, tmpR2, stats);// (I - xx)^-1
        if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN;
        cmat2x2_mul(inv_2x2T, tmpR2, tmp2x2);

        if(calc_uiz) cmat2x2_assign(tmp2x2, tmp2x2_uiz); // 为后续计算空间导数备份

        cmat2x2_mul(R_EV, tmp2x2, tmp2x2);
        tmpRL = invT / (RONE - RUL_FA * RDL_AL);
        tmpRL2 = R_EVL * tmpRL;

        for(MYINT i=0; i<SRC_M_NUM; ++i){
            get_qwv(ircvup, tmp2x2, tmpRL2, RU_FA, RUL_FA, src_coef[i], QWV[i]);
        }


        if(calc_uiz){
            calc_uiz_R_EV(rcv_xa, rcv_xb, ircvup, k, RD_BL, RDL_BL, uiz_R_EV, puiz_R_EVL);    
            cmat2x2_mul(uiz_R_EV, tmp2x2_uiz, tmp2x2_uiz);
            tmpRL2 = uiz_R_EVL * tmpRL;
            
            for(MYINT i=0; i<SRC_M_NUM; ++i){
                get_qwv(ircvup, tmp2x2_uiz, tmpRL2, RU_FA, RUL_FA, src_coef[i], QWV_uiz[i]);
            }
        }

    } // END if



    BEFORE_RETURN:

    // 对一些特殊情况的修正
    // 当震源和场点均位于地表时，可理论验证DS分量恒为0，这里直接赋0以避免后续的精度干扰
    if(mod1d->lays[isrc].dep == RZERO && mod1d->lays[ircv].dep == RZERO)
    {
        for(MYINT c=0; c<QWV_NUM; ++c){
            QWV[SRC_M_DS_INDEX][c] = CZERO;
            if(calc_uiz)  QWV_uiz[SRC_M_DS_INDEX][c] = CZERO;
        }
    }

}

