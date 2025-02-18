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

#include "propagate.h"
#include "bessel.h"
#include "model.h"
#include "matrix.h"
#include "layer.h"
#include "source.h"
#include "prtdbg.h"

#define CMAT_ASSIGN_SPLIT 0  // 2x2的小矩阵赋值合并为1个循环，程序速度提升微小




void kernel(
    const MODEL1D *mod1d, MYCOMPLEX omega, MYREAL k,
    MYCOMPLEX EXP_qwv[3][3], MYCOMPLEX VF_qwv[3][3], MYCOMPLEX HF_qwv[3][3], MYCOMPLEX DC_qwv[3][3])
{
    // 初始化qwv为0
    for(MYINT i=0; i<3; ++i){
        for(MYINT j=0; j<3; ++j){
            if(EXP_qwv!=NULL) EXP_qwv[i][j] = RZERO;
            if(VF_qwv!=NULL)  VF_qwv[i][j] = RZERO;
            if(HF_qwv!=NULL)  HF_qwv[i][j] = RZERO;
            if(DC_qwv!=NULL)  DC_qwv[i][j] = RZERO;
        }
    }

    bool ircvup = mod1d->ircvup;
    MYINT isrc = mod1d->isrc; // 震源所在虚拟层位, isrc>=1
    MYINT ircv = mod1d->ircv; // 接收点所在虚拟层位, ircv>=1, ircv != isrc
    MYINT imin, imax; // 相对浅层深层层位
    imin = mod1d->imin;
    imax = mod1d->imax;
    

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
    MYCOMPLEX RD[2][2], RDL, TD[2][2], TDL;
    MYCOMPLEX RU[2][2], RUL, TU[2][2], TUL;
    MYCOMPLEX *const pRDL = &RDL;
    MYCOMPLEX *const pTDL = &TDL;
    MYCOMPLEX *const pRUL = &RUL;
    MYCOMPLEX *const pTUL = &TUL;


    // 自由表面的反射系数
    MYCOMPLEX R_tilt[2][2] = INIT_C_ZERO_2x2_MATRIX; // SH波在自由表面的反射系数为1，不必定义变量

    // 接收点处的接收矩阵
    MYCOMPLEX R_EV[2][2], R_EVL;
    MYCOMPLEX *const pR_EVL = &R_EVL;
    

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
        mod1d_xa1 = CSQRT(RONE - mod1d_kaka1/(k*k));
        mod1d_xb1 = CSQRT(RONE - mod1d_kbkb1/(k*k));

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
                TD, pTDL, TU, pTUL);
        }

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
                    TD_FA, pTDL_FA, TU_FA, pTUL_FA);  
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
                    TD_RS, pTDL_RS, TU_RS, pTUL_RS);  // 写入原地址
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
                        TD_BL, pTDL_BL, TU_BL, pTUL_BL);  // 写入原地址
                } else {
                    recursion_RT_2x2(
                        RD_BL, RDL_BL, RU_BL, RUL_BL, 
                        TD_BL, TDL_BL, TU_BL, TUL_BL,
                        RD, RDL, RU, RUL, 
                        TD, TDL, TU, TUL,
                        RD_BL, pRDL_BL, NULL, NULL, 
                        NULL, NULL, NULL, NULL);  // 写入原地址
                }
                
            }
        } // END if


    } // END for loop 
    //===================================================================================

    // return;


    // 计算震源系数
    MYCOMPLEX EXP[3][3][2], VF[3][3][2], HF[3][3][2], DC[3][3][2];
    MYCOMPLEX (*pEXP)[3][2] = (EXP_qwv!=NULL)? EXP : NULL; 
    MYCOMPLEX (*pVF)[3][2]  = (VF_qwv!=NULL)?  VF  : NULL; 
    MYCOMPLEX (*pHF)[3][2]  = (HF_qwv!=NULL)?  HF  : NULL; 
    MYCOMPLEX (*pDC)[3][2]  = (DC_qwv!=NULL)?  DC  : NULL; 
    for(MYINT i=0; i<3; ++i){
        for(MYINT j=0; j<3; ++j){
            for(MYINT p=0; p<2; ++p){
                EXP[i][j][p] = VF[i][j][p] = HF[i][j][p] = DC[i][j][p] = RZERO;
            }
        }
    }
    source_coef(src_xa, src_xb, src_kaka, src_kbkb, omega, k, pEXP, pVF, pHF, pDC);

    // 临时中转矩阵 (temperary)
    MYCOMPLEX tmpR1[2][2], tmpR2[2][2], tmp2x2[2][2], tmpRL;
    MYCOMPLEX inv_2x2T[2][2], invT;

    // 递推RU_FA
    calc_R_tilt(top_xa, top_xb, top_kbkb, k, R_tilt);
    recursion_RU(
        R_tilt, RONE, 
        RD_FA, RDL_FA,
        RU_FA, RUL_FA, 
        TD_FA, TDL_FA,
        TU_FA, TUL_FA,
        RU_FA, pRUL_FA, NULL, NULL);

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
            RU_FB, pRUL_FB, inv_2x2T, &invT);
        
        // 公式(5.7.12-14)
        cmat2x2_mul(R_EV, inv_2x2T, tmpR1);
        cmat2x2_mul(RD_BL, RU_FB, tmpR2);
        cmat2x2_one_sub(tmpR2);
        cmat2x2_inv(tmpR2, tmpR2);// (I - xx)^-1
        cmat2x2_mul(tmpR1, tmpR2, tmp2x2);
        
        tmpRL = R_EVL * invT  / (RONE - RDL_BL * RUL_FB);
        for(MYINT m=0; m<3; ++m){
            if(0==m){
                // 爆炸源
                if(EXP_qwv!=NULL) get_qwv(ircvup, tmp2x2, tmpRL, RD_BL, RDL_BL, EXP[m], EXP_qwv[m]);
                // 垂直力源
                if(VF_qwv!=NULL)  get_qwv(ircvup, tmp2x2, tmpRL, RD_BL, RDL_BL, VF[m], VF_qwv[m]);
            }
            
            // 水平力源
            if(1==m && HF_qwv!=NULL) get_qwv(ircvup, tmp2x2, tmpRL, RD_BL, RDL_BL, HF[m], HF_qwv[m]);

            // 剪切位错
            if(DC_qwv!=NULL)  get_qwv(ircvup, tmp2x2, tmpRL, RD_BL, RDL_BL, DC[m], DC_qwv[m]);
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
            RD_AL, pRDL_AL, inv_2x2T, &invT);
        
        // 公式(5.7.26-27)
        cmat2x2_mul(R_EV, inv_2x2T, tmpR1);
        cmat2x2_mul(RU_FA, RD_AL, tmpR2);
        cmat2x2_one_sub(tmpR2);
        cmat2x2_inv(tmpR2, tmpR2);// (I - xx)^-1
        cmat2x2_mul(tmpR1, tmpR2, tmp2x2);
        tmpRL = R_EVL * invT / (RONE - RUL_FA * RDL_AL);
        for(MYINT m=0; m<3; ++m){
            if(0==m){
                // 爆炸源
                if(EXP_qwv!=NULL) get_qwv(ircvup, tmp2x2, tmpRL, RU_FA, RUL_FA, EXP[m], EXP_qwv[m]);
                // 垂直力源
                if(VF_qwv!=NULL)  get_qwv(ircvup, tmp2x2, tmpRL, RU_FA, RUL_FA, VF[m], VF_qwv[m]);
            }
            
            // 水平力源
            if(1==m && HF_qwv!=NULL) get_qwv(ircvup, tmp2x2, tmpRL, RU_FA, RUL_FA, HF[m], HF_qwv[m]);

            // 剪切位错
            if(DC_qwv!=NULL)  get_qwv(ircvup, tmp2x2, tmpRL, RU_FA, RUL_FA, DC[m], DC_qwv[m]);

        }
    } // END if


}



void int_Pk(
    MYREAL k, MYREAL r, 
    // F(ki,w)， 第一个维度3代表阶数m=0,1,2，第二个维度3代表三类系数qm,wm,vm 
    const MYCOMPLEX EXP_qwv[3][3], const MYCOMPLEX VF_qwv[3][3], 
    const MYCOMPLEX HF_qwv[3][3],  const MYCOMPLEX DC_qwv[3][3], 
    // F(ki,w)Jm(ki*r)ki，维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型
    MYCOMPLEX EXP_J[3][4], MYCOMPLEX VF_J[3][4], 
    MYCOMPLEX HF_J[3][4],  MYCOMPLEX DC_J[3][4] )
{
    MYREAL bj0k, bj1k, bj2k;
    MYREAL kr = k*r;

    bessel012(kr, &bj0k, &bj1k, &bj2k); 

    bj0k = bj0k*k;
    bj1k = bj1k*k;
    bj2k = bj2k*k;

    
    if(EXP_qwv!=NULL){
    // 公式(5.6.22), 将公式分解为F(k,w)Jm(kr)k的形式
    // m=0 爆炸源
    EXP_J[0][0] = - EXP_qwv[0][0]*bj1k;
    EXP_J[0][2] =   EXP_qwv[0][1]*bj0k;
    }

    if(VF_qwv!=NULL){
    // m=0 垂直力源
    VF_J[0][0] = - VF_qwv[0][0]*bj1k;
    VF_J[0][2] =   VF_qwv[0][1]*bj0k;
    }

    if(HF_qwv!=NULL){
    // m=1 水平力源
    HF_J[1][0]  =   HF_qwv[1][0]*bj0k;         // q1*J0*k
    HF_J[1][1]  = - (HF_qwv[1][0] + HF_qwv[1][2])*bj1k/(kr);    // - (q1+v1)*J1*k/kr
    HF_J[1][2]  =   HF_qwv[1][1]*bj1k;         // w1*J1*k
    HF_J[1][3]  = - HF_qwv[1][2]*bj0k;         // -v1*J0*k
    }

    if(DC_qwv!=NULL){
    // m=0 双力偶源
    DC_J[0][0] = - DC_qwv[0][0]*bj1k;
    DC_J[0][2] =   DC_qwv[0][1]*bj0k;

    // m=1 双力偶源
    DC_J[1][0]  =   DC_qwv[1][0]*bj0k;         // q1*J0*k
    DC_J[1][1]  = - (DC_qwv[1][0] + DC_qwv[1][2])*bj1k/(kr);    // - (q1+v1)*J1*k/kr
    DC_J[1][2]  =   DC_qwv[1][1]*bj1k;         // w1*J1*k
    DC_J[1][3]  = - DC_qwv[1][2]*bj0k;         // -v1*J0*k

    // m=2 双力偶源
    DC_J[2][0]  =   DC_qwv[2][0]*bj1k;         // q2*J1*k
    DC_J[2][1]  = - RTWO*(DC_qwv[2][0] + DC_qwv[2][2])*bj2k/(kr);    // - (q2+v2)*J2*k/kr
    DC_J[2][2]  =   DC_qwv[2][1]*bj2k;         // w2*J2*k
    DC_J[2][3]  = - DC_qwv[2][2]*bj1k;         // -v2*J1*k
    }
}



void merge_Pk(
    // F(ki,w)Jm(ki*r)ki，维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型
    const MYCOMPLEX sum_EXP_J[3][4], const MYCOMPLEX sum_VF_J[3][4], 
    const MYCOMPLEX sum_HF_J[3][4],  const MYCOMPLEX sum_DC_J[3][4], 
    // 累积求和，维度2代表Z、R分量，维度3代表Z、R、T分量 
    MYCOMPLEX tol_EXP[2], MYCOMPLEX tol_VF[2], MYCOMPLEX tol_HF[3],
    MYCOMPLEX tol_DD[2],  MYCOMPLEX tol_DS[3], MYCOMPLEX tol_SS[3])
{   
    if(sum_EXP_J!=NULL){
    tol_EXP[0] = sum_EXP_J[0][2];
    tol_EXP[1] = sum_EXP_J[0][0];
    }

    if(sum_VF_J!=NULL){
    tol_VF[0] = sum_VF_J[0][2];
    tol_VF[1] = sum_VF_J[0][0];
    }

    if(sum_HF_J!=NULL){
    tol_HF[0] = sum_HF_J[1][2];
    tol_HF[1] = sum_HF_J[1][0] + sum_HF_J[1][1];
    tol_HF[2] = - sum_HF_J[1][1] + sum_HF_J[1][3];
    }

    if(sum_DC_J!=NULL){
    tol_DD[0] = sum_DC_J[0][2];
    tol_DD[1] = sum_DC_J[0][0];
    
    tol_DS[0] = sum_DC_J[1][2];
    tol_DS[1] = sum_DC_J[1][0] + sum_DC_J[1][1];
    tol_DS[2] = - sum_DC_J[1][1] + sum_DC_J[1][3];

    tol_SS[0] = sum_DC_J[2][2];
    tol_SS[1] = sum_DC_J[2][0] + sum_DC_J[2][1];
    tol_SS[2] = - sum_DC_J[2][1] + sum_DC_J[2][3];
    }
}




void recursion_RD(
    const MYCOMPLEX RD1[2][2], MYCOMPLEX RDL1, const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT)
{
    MYCOMPLEX tmp1[2][2], tmp2[2][2], inv1;

    // RD, RDL
    cmat2x2_mul(RU1, RD2, tmp1);
    cmat2x2_one_sub(tmp1);
    cmat2x2_inv(tmp1, tmp1);
    cmat2x2_mul(tmp1, TD1, tmp2);
    if(inv_2x2T!=NULL) cmat2x2_assign(tmp2, inv_2x2T);

    cmat2x2_mul(RD2, tmp2, tmp1);
    cmat2x2_mul(TU1, tmp1, tmp2);
    cmat2x2_add(RD1, tmp2, RD);
    inv1 = RONE / (RONE - RUL1*RDL2) * TDL1;
    *RDL = RDL1 + TUL1*RDL2*inv1;
    if(invT!=NULL)  *invT = inv1;
}


void recursion_TD(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, 
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, 
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT)
{
    MYCOMPLEX tmp1[2][2], tmp2[2][2], inv1;

    // TD, TDL
    cmat2x2_mul(RU1, RD2, tmp2);
    cmat2x2_one_sub(tmp2);
    cmat2x2_inv(tmp2, tmp1);
    cmat2x2_mul(tmp1, TD1, tmp2);
    if(inv_2x2T!=NULL)  cmat2x2_assign(tmp2, inv_2x2T);
    cmat2x2_mul(TD2, tmp2, TD);
    
    inv1 = RONE / (RONE - RUL1*RDL2) * TDL1;
    *TDL = TDL2 * inv1;

    if(invT!=NULL) *invT = inv1;
}


void recursion_RU(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, const MYCOMPLEX RU2[2][2], MYCOMPLEX RUL2,
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT)
{
    MYCOMPLEX tmp1[2][2], tmp2[2][2], inv1;

    // RU, RUL
    cmat2x2_mul(RD2, RU1, tmp2);
    cmat2x2_one_sub(tmp2);
    cmat2x2_inv(tmp2, tmp1);
    cmat2x2_mul(tmp1, TU2, tmp2);
    if(inv_2x2T!=NULL)  cmat2x2_assign(tmp2, inv_2x2T);

    cmat2x2_mul(RU1, tmp2, tmp1); 
    cmat2x2_mul(TD2, tmp1, tmp2);
    cmat2x2_add(RU2, tmp2, RU);
    inv1 = RONE / (RONE - RUL1*RDL2) * TUL2;
    *RUL = RUL2 + TDL2*RUL1*inv1; 

    if(invT!=NULL)  *invT = inv1;
}


void recursion_TU(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2,
    const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT)
{
    MYCOMPLEX tmp1[2][2], tmp2[2][2], inv1;

    // TU, TUL
    cmat2x2_mul(RD2, RU1, tmp2);
    cmat2x2_one_sub(tmp2);
    cmat2x2_inv(tmp2, tmp1);
    cmat2x2_mul(tmp1, TU2, tmp2);
    if(inv_2x2T!=NULL) cmat2x2_assign(tmp2, inv_2x2T);
    cmat2x2_mul(TU1, tmp2, TU);
    
    inv1 = RONE / (RONE - RUL1*RDL2) * TUL2;
    *TUL = TUL1 * inv1;

    if(invT!=NULL)  *invT = inv1;

}




void recursion_RT_2x2(
    const MYCOMPLEX RD1[2][2], MYCOMPLEX RDL1, const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, const MYCOMPLEX RU2[2][2], MYCOMPLEX RUL2,
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL,
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL)
{

    // 临时矩阵
    MYCOMPLEX tmp1[2][2], tmp2[2][2];
    MYCOMPLEX inv0, inv1T;

    inv0 = RONE / (RONE - RUL1*RDL2);
    // return;

    // Rayleigh RD,TD
    if( RD!=NULL || TD!=NULL ){
        cmat2x2_mul(RU1, RD2, tmp1);
        cmat2x2_one_sub(tmp1);
        cmat2x2_inv(tmp1, tmp1);
        cmat2x2_mul(tmp1, TD1, tmp2);

        // TD
        if(TD!=NULL){
            cmat2x2_mul(TD2, tmp2, TD); // 相同的逆阵，节省计算量
        }

        // RD
        if(RD!=NULL){
            cmat2x2_mul(RD2, tmp2, tmp1);
            cmat2x2_mul(TU1, tmp1, tmp2);
            cmat2x2_add(RD1, tmp2, RD);
        }
    }

    // Rayleigh RU,TU
    if( RU!=NULL || TU!=NULL ){
        cmat2x2_mul(RD2, RU1, tmp1);
        cmat2x2_one_sub(tmp1);
        cmat2x2_inv(tmp1, tmp1);
        cmat2x2_mul(tmp1, TU2, tmp2);

        // TU
        if(TU!=NULL){
            cmat2x2_mul(TU1, tmp2, TU);
        }

        // RU
        if(RU!=NULL){
            cmat2x2_mul(RU1, tmp2, tmp1);
            cmat2x2_mul(TD2, tmp1, tmp2);
            cmat2x2_add(RU2, tmp2, RU);
        }
    }


    // Love RDL,TDL
    if(RDL!=NULL || TDL!=NULL){
        inv1T = inv0 * TDL1;
        // TDL
        if(TDL!=NULL){
            *TDL = TDL2 * inv1T;
        }
        // RDL
        if(RDL!=NULL){
            *RDL = RDL1 + TUL1*RDL2*inv1T;
        }
    }

    // Love RUL,TUL
    if(RUL!=NULL || TUL!=NULL){
        inv1T = inv0 * TUL2;
        // TUL
        if(TUL!=NULL){
            *TUL = TUL1 * inv1T;
        }

        // RUL
        if(RUL!=NULL){
            *RUL = RUL2 + TDL2*RUL1 *inv1T; 
        }
    }

}


void recursion_RT_2x2_imaginary(
    MYCOMPLEX xa1, MYCOMPLEX xb1, MYREAL thk, MYREAL k, // 使用上层的厚度
    MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL)
{
    MYCOMPLEX exa, exb, exab, ex2a, ex2b; 
    exa = CEXP(-k*thk*xa1);
    exb = CEXP(-k*thk*xb1);

    exab = exa * exb;
    ex2a = exa * exa;
    ex2b = exb * exb;


    // 虚拟层位不是介质物理间断面
    RU[0][0] *= ex2a;    RU[0][1] *= exab;  
    RU[1][0] *= exab;    RU[1][1] *= ex2b;  
    
    TD[0][0] *= exa;     TD[0][1] *= exa; 
    TD[1][0] *= exb;     TD[1][1] *= exb;

    TU[0][0] *= exa;     TU[0][1] *= exb; 
    TU[1][0] *= exa;     TU[1][1] *= exb;

    *RUL *= ex2b;
    *TDL *= exb;
    *TUL *= exb;
}




void get_qwv(
    bool ircvup, 
    const MYCOMPLEX R1[2][2], MYCOMPLEX RL1, 
    const MYCOMPLEX R2[2][2], MYCOMPLEX RL2, 
    const MYCOMPLEX coef[3][2], MYCOMPLEX qwv[3])
{
    MYCOMPLEX qw0[2], qw1[2], v0;
    MYCOMPLEX coefD[2] = {coef[0][0], coef[1][0]};
    MYCOMPLEX coefU[2] = {coef[0][1], coef[1][1]};
    if(ircvup){
        cmat2x1_mul(R2, coefD, qw0);
        qw0[0] += coefU[0]; qw0[1] += coefU[1]; 
        v0 = RL1 * (RL2*coef[2][0] + coef[2][1]);
    } else {
        cmat2x1_mul(R2, coefU, qw0);
        qw0[0] += coefD[0]; qw0[1] += coefD[1]; 
        v0 = RL1 * (coef[2][0] + RL2*coef[2][1]);
    }
    cmat2x1_mul(R1, qw0, qw1);

    qwv[0] = qw1[0];
    qwv[1] = qw1[1];
    qwv[2] = v0;
}

