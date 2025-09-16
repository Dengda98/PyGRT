/**
 * @file   layer.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 P-SV 波和 SH 波的反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient 
 *               matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>

#include "grt/dynamic/layer.h"
#include "grt/common/model.h"
#include "grt/common/prtdbg.h"
#include "grt/common/matrix.h"

#include "grt/common/checkerror.h"

void grt_topfree_RU_PSV(MYCOMPLEX xa0, MYCOMPLEX xb0, MYCOMPLEX cbcb0, MYREAL k, MYCOMPLEX R_tilt[2][2], MYINT *stats)
{
    if(cbcb0 != 0.0){
        // 固体表面
        // 公式(5.3.10-14)
        MYCOMPLEX Delta = 0.0;
        // MYREAL kk = k*k; 
        MYCOMPLEX cbcb02 = 0.25*cbcb0*cbcb0;

        // 对公式(5.3.10-14)进行重新整理，对浮点数友好一些
        Delta = -1.0 + xa0*xb0 + cbcb0 - cbcb02;
        if(Delta == 0.0){
            *stats = GRT_INVERSE_FAILURE;
            return;
        }
        R_tilt[0][0] = (1.0 + xa0*xb0 - cbcb0 + cbcb02) / Delta;
        R_tilt[0][1] = 2.0 * xb0 * (1.0 - 0.5*cbcb0) / Delta;
        R_tilt[1][0] = 2.0 * xa0 * (1.0 - 0.5*cbcb0) / Delta;
        R_tilt[1][1] = R_tilt[0][0];
    }
    else {
        // 液体表面
        R_tilt[0][0] = -1.0;
        R_tilt[1][1] = R_tilt[0][1] = R_tilt[1][0] = 0.0;
    }
}


void grt_wave2qwv_REV_PSV(
    MYCOMPLEX xa_rcv, MYCOMPLEX xb_rcv, bool ircvup, MYREAL k, 
    const MYCOMPLEX R[2][2], MYCOMPLEX R_EV[2][2])
{
    MYCOMPLEX D11[2][2], D12[2][2];
    if(xb_rcv != 1.0){
        // 位于固体层
        // 公式(5.2.19)
        D11[0][0] = k;         D11[0][1] = k*xb_rcv;
        D11[1][0] = k*xa_rcv;  D11[1][1] = k;
        D12[0][0] = k;         D12[0][1] = -k*xb_rcv;
        D12[1][0] = -k*xa_rcv; D12[1][1] = k;
    } else {
        // 位于液体层
        D11[0][0] = k;         D11[0][1] = 0.0;
        D11[1][0] = k*xa_rcv;  D11[1][1] = 0.0;
        D12[0][0] = k;         D12[0][1] = 0.0;
        D12[1][0] = -k*xa_rcv; D12[1][1] = 0.0;
    }

    // 公式(5.7.7,25)
    if(ircvup){// 震源更深
        grt_cmat2x2_mul(D12, R, R_EV);
        grt_cmat2x2_add(D11, R_EV, R_EV);
    } else { // 接收点更深
        grt_cmat2x2_mul(D11, R, R_EV);
        grt_cmat2x2_add(D12, R_EV, R_EV);
    }
}


void grt_wave2qwv_REV_SH(MYCOMPLEX xb_rcv, MYREAL k, MYCOMPLEX RL, MYCOMPLEX *R_EVL)
{
    if(xb_rcv != 1.0){
        // 位于固体层
        // 公式(5.2.19)
        *R_EVL = (1.0 + (RL))*k;
    } else {
        // 位于液体层
        *R_EVL = 0.0;
    }
}


void grt_wave2qwv_z_REV_PSV(
    MYCOMPLEX xa_rcv, MYCOMPLEX xb_rcv, bool ircvup,
    MYREAL k, 
    const MYCOMPLEX R[2][2], MYCOMPLEX R_EV[2][2])
{
    // 将垂直波函数转为ui,z在(B_m, P_m, C_m)系下的分量
    // 新推导的公式
    MYCOMPLEX ak = k*k*xa_rcv;
    MYCOMPLEX bk = k*k*xb_rcv;
    MYCOMPLEX bb = xb_rcv*bk;
    MYCOMPLEX aa = xa_rcv*ak;
    MYCOMPLEX D11[2][2] = {{ak, bb}, {aa, bk}};
    MYCOMPLEX D12[2][2] = {{-ak, bb}, {aa, -bk}};

    // 位于液体层
    if(xb_rcv == 1.0){
        D11[0][1] = D11[1][1] = D12[0][1] = D12[1][1] = 0.0;
    }

    // 公式(5.7.7,25)
    if(ircvup){// 震源更深
        grt_cmat2x2_mul(D12, R, R_EV);
        grt_cmat2x2_add(D11, R_EV, R_EV);
    } else { // 接收点更深
        grt_cmat2x2_mul(D11, R, R_EV);
        grt_cmat2x2_add(D12, R_EV, R_EV);
    }
}    


void grt_wave2qwv_z_REV_SH(MYCOMPLEX xb_rcv, bool ircvup, MYREAL k, MYCOMPLEX RL, MYCOMPLEX *R_EVL)
{
    // 将垂直波函数转为ui,z在(B_m, P_m, C_m)系下的分量
    // 新推导的公式
    MYCOMPLEX bk = k*k*xb_rcv;

    if(xb_rcv != 1.0){
        // 位于固体层
        if(ircvup){// 震源更深
            *R_EVL = (1.0 - (RL))*bk;
        } else { // 接收点更深
            *R_EVL = (RL - 1.0)*bk;
        }
    } else {
        // 位于液体层
        *R_EVL = 0.0;
    }
}    


void grt_RT_matrix_ll_PSV(
    MYREAL Rho1, MYCOMPLEX xa1,
    MYREAL Rho2, MYCOMPLEX xa2,
    MYCOMPLEX omega, MYREAL k,
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2], 
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2], MYINT *stats)
{
    MYCOMPLEX A = xa1*Rho2 + xa2*Rho1;
    if(A==0.0){
        *stats = GRT_INVERSE_FAILURE;
        return;
    }

    RD[0][0] = (xa1*Rho2 - xa2*Rho1)/A;  
    RD[0][1] = RD[1][0] = RD[1][1] = 0.0;
    
    RU[0][0] = (xa2*Rho1 - xa1*Rho2)/A;
    RU[0][1] = RU[1][0] = RU[1][1] = 0.0;

    TD[0][0] = 2.0*xa1*Rho1/A;
    TD[0][1] = TD[1][0] = TD[1][1] = 0.0;

    TU[0][0] = 2.0*xa2*Rho2/A;
    TU[0][1] = TU[1][0] = TU[1][1] = 0.0;

}

void grt_RT_matrix_ll_SH(MYCOMPLEX *RDL, MYCOMPLEX *RUL, MYCOMPLEX *TDL, MYCOMPLEX *TUL)
{
    *RDL = 0.0;
    *RUL = 0.0;
    *TDL = 0.0;
    *TUL = 0.0;
}



void grt_RT_matrix_ls_PSV(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX cbcb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX cbcb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k,
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2], 
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2], MYINT *stats)
{
    // 后缀1表示上层的液体的物理参数，后缀2表示下层的固体的物理参数
    // 若mu2==0, 则下层为液体，参数需相互交换 

    // 讨论液-固 or 固-液
    bool isfluidUp = (mu1 == 0.0);  // 上层是否为液体
    MYINT sgn = 1;
    if(isfluidUp && mu2 == 0.0){
        GRTRaiseError("Error: fluid-fluid interface is not allowed in function %s\n", __func__);
    }

    // 使用指针
    MYCOMPLEX (*pRD)[2], (*pRU)[2];
    MYCOMPLEX (*pTD)[2], (*pTU)[2];
    if(isfluidUp){
        pRD = RD; pRU = RU; 
        pTD = TD; pTU = TU; 
    } else {
        pRD = RU; pRU = RD;
        pTD = TU; pTU = TD;
        GRT_SWAP(MYREAL, Rho1, Rho2);
        GRT_SWAP(MYCOMPLEX, xa1, xa2);
        GRT_SWAP(MYCOMPLEX, xb1, xb2);
        GRT_SWAP(MYCOMPLEX, cbcb1, cbcb2);
        GRT_SWAP(MYCOMPLEX, mu1, mu2);
        sgn = -1;
    }

    
    // 定义一些中间变量来简化运算和书写
    MYCOMPLEX lamka1k = Rho1*GRT_SQUARE(omega/k);
    MYCOMPLEX kb2k = cbcb2;
    MYCOMPLEX Og2k = 1.0 - 0.5*kb2k;
    MYCOMPLEX Og2k2 = Og2k*Og2k;
    MYCOMPLEX A = 2.0*Og2k2*xa1*mu2 + 0.5*lamka1k*kb2k*xa2 - 2.0*mu2*xa1*xa2*xb2;
    MYCOMPLEX B = 2.0*Og2k2*xa1*mu2 - 0.5*lamka1k*kb2k*xa2 + 2.0*mu2*xa1*xa2*xb2;
    MYCOMPLEX C = 2.0*Og2k2*xa1*mu2 + 0.5*lamka1k*kb2k*xa2 + 2.0*mu2*xa1*xa2*xb2;
    MYCOMPLEX D = 2.0*Og2k2*xa1*mu2 - 0.5*lamka1k*kb2k*xa2 - 2.0*mu2*xa1*xa2*xb2;

    if(A == 0.0){
        *stats = GRT_INVERSE_FAILURE;
        return;
    }
    
    // 按液体层在上层处理
    pRD[0][0] = D/A; 
    pRD[0][1] = pRD[1][0] = pRD[1][1] = 0.0;

    pRU[0][0] = - B/A;
    pRU[0][1] = - 4.0*Og2k*xa1*xb2*mu2/A * sgn;
    pRU[1][0] = pRU[0][1]/xb2 * xa2;
    pRU[1][1] = - C/A;

    pTD[0][0] = - 2.0*Og2k*xa1*lamka1k/A;      pTD[0][1] = 0.0;
    pTD[1][0] = pTD[0][0]/Og2k*xa2 * sgn;       pTD[1][1] = 0.0;

    pTU[0][0] = - 2.0*Og2k*xa2*mu2*kb2k/A;     pTU[0][1] = pTU[0][0]/Og2k*xb2 * sgn;
    pTU[1][0] = pTU[1][1] = 0.0;

}


void grt_RT_matrix_ls_SH(
    MYCOMPLEX xb1, MYCOMPLEX mu1, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k,
    MYCOMPLEX *RDL, MYCOMPLEX *RUL, 
    MYCOMPLEX *TDL, MYCOMPLEX *TUL)
{
    // 后缀1表示上层的液体的物理参数，后缀2表示下层的固体的物理参数
    // 若mu2==0, 则下层为液体，参数需相互交换 

    // 讨论液-固 or 固-液
    bool isfluidUp = (mu1 == 0.0);  // 上层是否为液体
    if(isfluidUp && mu2 == 0.0){
        GRTRaiseError("Error: fluid-fluid interface is not allowed in function %s\n", __func__);
    }

    // 使用指针
    MYCOMPLEX *pRDL, *pRUL;
    MYCOMPLEX *pTDL, *pTUL;
    if(isfluidUp){
        pRDL = RDL; pRUL = RUL;
        pTDL = TDL; pTUL = TUL;
    } else {
        pRDL = RUL; pRUL = RDL;
        pTDL = TUL; pTUL = TDL;
    }

    *pRDL = 0.0;
    *pRUL = 1.0;
    *pTDL = 0.0;
    *pTUL = 0.0;
}



void grt_RT_matrix_ss_PSV(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX cbcb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX cbcb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k,
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2],
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2], MYINT *stats)
{
    MYCOMPLEX tmp;

    // 定义一些中间变量来简化运算和书写
    // MYREAL kk = k*k;
    MYCOMPLEX dmu = mu1/mu2 - 1.0; // mu1 - mu2; 分子分母同除mu2
    MYCOMPLEX dmu2 = dmu*dmu;

    MYCOMPLEX mu1cbcb1 = mu1/mu2*cbcb1;// mu1*kb1_k2;
    MYCOMPLEX mu2cbcb2 = cbcb2; // mu2*kb2_k2;

    MYREAL rho12 = Rho1 / Rho2;
    MYREAL rho21 = Rho2 / Rho1;

    // 从原公式上，分母包含5项，但前四项会随着k的增大迅速超过最后一项
    // 最后一项要小前几项10余个数量级，但计算结果还是保持在最后一项的量级，
    // 这种情况会受到浮点数的有效位数的限制，64bit的双精度double大概就是15-16位，
    // 故会发生严重精度损失的情况。目前只在实部上观察到这个现象，虚部基本都在相近量级(或许是相对不明显)
    // 
    // 以下对公式重新整理，提出k的高阶项，以避免上述问题
    MYCOMPLEX Delta;
    Delta =   dmu2*(1.0-xa1*xb1)*(1.0-xa2*xb2) + mu1cbcb1*dmu*(rho21*(1.0-xa1*xb1) - (1.0-xa2*xb2)) 
            + 0.25*mu1cbcb1*mu2cbcb2*(rho12*(1.0-xa2*xb2) + rho21*(1.0-xa1*xb1) - 2.0 - (xa1*xb2+xa2*xb1));

    if( Delta == 0.0 ){
        // printf("# zero Delta_inv=%e+%eJ\n", creal(Delta_inv), cimag(Delta_inv));
        *stats = GRT_INVERSE_FAILURE;
        return;
    } 

    // REFELCTION
    //------------------ RD -----------------------------------
    // rpp+
    RD[0][0] = ( - dmu2*(1.0+xa1*xb1)*(1.0-xa2*xb2) - mu1cbcb1*dmu*(rho21*(1.0+xa1*xb1) - (1.0-xa2*xb2))
                    - 0.25*mu1cbcb1*mu2cbcb2*(rho12*(1.0-xa2*xb2) + rho21*(1.0+xa1*xb1) - 2.0 + (xa1*xb2-xa2*xb1))) / Delta;
    // rsp+
    RD[0][1] = ( - dmu2*(1.0-xa2*xb2) + 0.5*mu1cbcb1*dmu*((1.0-xa2*xb2) - 2.0*rho21) 
                    + 0.25*mu1cbcb1*mu2cbcb2*(1.0-rho21)) / Delta * (-2.0*xb1);
    // rps+
    RD[1][0] = RD[0][1]*(xa1/xb1);
    // rss+
    RD[1][1] = ( - dmu2*(1.0+xa1*xb1)*(1.0-xa2*xb2) - mu1cbcb1*dmu*(rho21*(1.0+xa1*xb1) - (1.0-xa2*xb2))
                    - 0.25*mu1cbcb1*mu2cbcb2*(rho12*(1.0-xa2*xb2) + rho21*(1.0+xa1*xb1) - 2.0 - (xa1*xb2-xa2*xb1))) / Delta;
    //------------------ RU -----------------------------------
    // rpp-
    RU[0][0] = ( - dmu2*(1.0-xa1*xb1)*(1.0+xa2*xb2) - mu1cbcb1*dmu*(rho21*(1.0-xa1*xb1) - (1.0+xa2*xb2))
                    - 0.25*mu1cbcb1*mu2cbcb2*(rho12*(1.0+xa2*xb2) + rho21*(1.0-xa1*xb1) - 2.0 - (xa1*xb2-xa2*xb1))) / Delta;
    // rsp-
    RU[0][1] = ( - dmu2*(1.0-xa1*xb1) - 0.5*mu1cbcb1*dmu*(rho21*(1.0-xa1*xb1) - 2.0)
                    + 0.25*mu1cbcb1*mu2cbcb2*(1.0-rho12)) / Delta * (2.0*xb2);
    // rps-
    RU[1][0] = RU[0][1]*(xa2/xb2);
    // rss-
    RU[1][1] = ( - dmu2*(1.0-xa1*xb1)*(1.0+xa2*xb2) - mu1cbcb1*dmu*(rho21*(1.0-xa1*xb1) - (1.0+xa2*xb2))
                    - 0.25*mu1cbcb1*mu2cbcb2*(rho12*(1.0+xa2*xb2) + rho21*(1.0-xa1*xb1) - 2.0 + (xa1*xb2-xa2*xb1))) / Delta;

    // REFRACTION
    tmp = mu1cbcb1*xa1*(dmu*(xb2-xb1) - 0.5*mu1cbcb1*(rho21*xb1+xb2)) / Delta;
    TD[0][0] = tmp;     TU[0][0] = (rho21*xa2/xa1) * tmp;
    tmp = mu1cbcb1*xb1*(dmu*(1.0-xa1*xb2) - 0.5*mu1cbcb1*(1.0-rho21)) / Delta;
    TD[0][1] = tmp;     TU[1][0] = (rho21*xa2/xb1) * tmp;
    tmp = mu1cbcb1*xa1*(dmu*(1.0-xa2*xb1) - 0.5*mu1cbcb1*(1.0-rho21)) / Delta;
    TD[1][0] = tmp;     TU[0][1] = (rho21*xb2/xa1) * tmp;
    tmp = mu1cbcb1*xb1*(dmu*(xa2-xa1) - 0.5*mu1cbcb1*(rho21*xa1+xa2)) / Delta;
    TD[1][1] = tmp;     TU[1][1] = (rho21*xb2/xb1) * tmp;
}


void grt_RT_matrix_ss_SH(
    MYCOMPLEX xb1, MYCOMPLEX mu1, 
    MYCOMPLEX xb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX *RDL, MYCOMPLEX *RUL, 
    MYCOMPLEX *TDL, MYCOMPLEX *TUL)
{
    MYCOMPLEX tmp;

    // REFELCTION
    *RUL = (mu2*xb2 - mu1*xb1) / (mu2*xb2 + mu1*xb1) ;
    *RDL = - (*RUL);

    // REFRACTION
    tmp = 2.0 / (mu2*xb2 + mu1*xb1);
    *TDL = mu1*xb1 * tmp;
    *TUL = mu2*xb2 * tmp;
}



void grt_RT_matrix_PSV(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX cbcb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX cbcb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2],
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2], MYINT *stats)
{
    // 根据界面两侧的具体情况选择函数
    if(mu1 != 0.0 && mu2 != 0.0){
        grt_RT_matrix_ss_PSV(
            Rho1, xa1, xb1, cbcb1, mu1, 
            Rho2, xa2, xb2, cbcb2, mu2, 
            omega, k, 
            RD, RU, TD, TU, stats);
    }
    else if(mu1 == 0.0 && mu2 == 0.0){
        grt_RT_matrix_ll_PSV(
            Rho1, xa1,
            Rho2, xa2,
            omega, k, 
            RD, RU, TD, TU, stats);
    }
    else{
        grt_RT_matrix_ls_PSV(
            Rho1, xa1, xb1, cbcb1, mu1, 
            Rho2, xa2, xb2, cbcb2, mu2, 
            omega, k, 
            RD, RU, TD, TU, stats);
    }
}


void grt_RT_matrix_SH(
    MYCOMPLEX xb1, MYCOMPLEX mu1, 
    MYCOMPLEX xb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX *RDL, MYCOMPLEX *RUL, 
    MYCOMPLEX *TDL, MYCOMPLEX *TUL)
{
    // 根据界面两侧的具体情况选择函数
    if(mu1 != 0.0 && mu2 != 0.0){
        grt_RT_matrix_ss_SH(
            xb1, mu1, 
            xb2, mu2, 
            omega, k, 
            RDL, RUL, TDL, TUL);
    }
    else if(mu1 == 0.0 && mu2 == 0.0){
        grt_RT_matrix_ll_SH(
            RDL, RUL, TDL, TUL);
    }
    else{
        grt_RT_matrix_ls_SH(
            xb1, mu1, mu2, 
            omega, k, 
            RDL, RUL, TDL, TUL);
    }
}


void grt_delay_RT_matrix_PSV(
    MYCOMPLEX xa1, MYCOMPLEX xb1, 
    MYREAL thk, MYREAL k,
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2], 
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2])
{
    MYCOMPLEX exa, exb, ex2a, ex2b, exab;
    exa = exp(- k*thk*xa1);
    exb = exp(- k*thk*xb1);
    ex2a = exa * exa;
    ex2b = exb * exb;
    exab = exa * exb;

    RD[0][0] *= ex2a;   RD[0][1] *= exab;
    RD[1][0] *= exab;   RD[1][1] *= ex2b;

    TD[0][0] *= exa;    TD[0][1] *= exb;
    TD[1][0] *= exa;    TD[1][1] *= exb;

    TU[0][0] *= exa;    TU[0][1] *= exa;
    TU[1][0] *= exb;    TU[1][1] *= exb;
}


void grt_delay_RT_matrix_SH(
    MYCOMPLEX xb1, 
    MYREAL thk, MYREAL k,
    MYCOMPLEX *RDL, MYCOMPLEX *RUL, 
    MYCOMPLEX *TDL, MYCOMPLEX *TUL)
{
    MYCOMPLEX exb, ex2b;
    exb = exp(- k*thk*xb1);
    ex2b = exb * exb;

    *RDL *= ex2b;
    *TDL *= exb;
    *TUL *= exb;
}


void grt_get_layer_D(
    MYCOMPLEX xa, MYCOMPLEX xb, MYCOMPLEX kbkb, MYCOMPLEX mu,
    MYCOMPLEX omega, MYREAL rho, MYREAL k, MYCOMPLEX D[4][4], bool inverse, MYINT liquid_invtype)
{
    // 第iy层物理量
    MYCOMPLEX Omg;
    if(xb != 1.0){
        Omg = k*k - 0.5*kbkb;
        if( ! inverse ){
            D[0][0] = k;               D[0][1] = k*xb;             D[0][2] = k;               D[0][3] = -k*xb;     
            D[1][0] = k*xa;            D[1][1] = k;                D[1][2] = -k*xa;           D[1][3] = k;   
            D[2][0] = 2*mu*Omg;        D[2][1] = 2*k*mu*k*xb;      D[2][2] = 2*mu*Omg;        D[2][3] = -2*k*mu*k*xb;   
            D[3][0] = 2*k*mu*k*xa;     D[3][1] = 2*mu*Omg;         D[3][2] = -2*k*mu*k*xa;    D[3][3] = 2*mu*Omg;   
        } else {
            D[0][0] = -2*k*mu*k*xa*k*xb;  D[0][1] = 2*mu*Omg*k*xb;       D[0][2] = k*xa*k*xb;            D[0][3] = -k*k*xb;     
            D[1][0] = 2*mu*Omg*k*xa;      D[1][1] = -2*k*mu*k*xa*k*xb;   D[1][2] = -k*k*xa;              D[1][3] = k*xa*k*xb;   
            D[2][0] = -2*k*mu*k*xa*k*xb;  D[2][1] = -2*mu*Omg*k*xb;      D[2][2] = k*xa*k*xb;            D[2][3] = k*k*xb;   
            D[3][0] = -2*mu*Omg*k*xa;     D[3][1] = -2*k*mu*k*xa*k*xb;   D[3][2] = k*k*xa;               D[3][3] = k*xa*k*xb;
            for(MYINT i=0; i<4; ++i){
                for(MYINT j=0; j<4; ++j){
                    D[i][j] /= - 2*mu*kbkb*k*xa*k*xb;
                }
            }
        }
    } else {
        Omg = rho * GRT_SQUARE(omega) / k;
        if(liquid_invtype == 1){
            if( ! inverse ){
                D[0][0] = k;            D[0][1] = 0.0;      D[0][2] = k;               D[0][3] = 0.0;     
                D[1][0] = k*xa;         D[1][1] = 0.0;      D[1][2] = -k*xa;           D[1][3] = 0.0;   
                D[2][0] = - k*Omg;      D[2][1] = 0.0;      D[2][2] = - k*Omg;         D[2][3] = 0.0;   
                D[3][0] = 0.0;          D[3][1] = 0.0;      D[3][2] = 0.0;             D[3][3] = 0.0;   
            } else {
                D[0][0] = xa;                 D[0][1] = 1.0 + Omg*Omg;       D[0][2] = - xa * Omg;           D[0][3] = 0.0;     
                D[1][0] = 0.0;                D[1][1] = 0.0;                 D[1][2] = 0.0;                  D[1][3] = 0.0;   
                D[2][0] = xa;                 D[2][1] = - (1.0 + Omg*Omg);   D[2][2] = - xa * Omg;           D[2][3] = 0.0;   
                D[3][0] = 0.0;                D[3][1] = 0.0;                 D[3][2] = 0.0;                  D[3][3] = 0.0;
                for(MYINT i=0; i<4; ++i){
                    for(MYINT j=0; j<4; ++j){
                        D[i][j] /= 2*k*xa*(1.0 + Omg*Omg);
                    }
                }
            }
        } 
        else if(liquid_invtype == 2){
            // 此处液体层内的 D 和 D^{-1} 只考虑了 w, \sigma_R 两项，这是由边界条件决定的
            if( ! inverse ){
                D[0][0] = 0.0;          D[0][1] = 0.0;      D[0][2] = 0.0;             D[0][3] = 0.0;     
                D[1][0] = k*xa;         D[1][1] = 0.0;      D[1][2] = - k*xa;          D[1][3] = 0.0;   
                D[2][0] = - k*Omg;      D[2][1] = 0.0;      D[2][2] = - k*Omg;         D[2][3] = 0.0;   
                D[3][0] = 0.0;          D[3][1] = 0.0;      D[3][2] = 0.0;             D[3][3] = 0.0;   
            } else {
                D[0][0] = 0.0;         D[0][1] = 0.5/(k*xa);       D[0][2] = - 0.5 / (k*Omg);      D[0][3] = 0.0;     
                D[1][0] = 0.0;         D[1][1] = 0.0;              D[1][2] = 0.0;                  D[1][3] = 0.0;   
                D[2][0] = 0.0;         D[2][1] = - 0.5/(k*xa);     D[2][2] = - 0.5 / (k*Omg);      D[2][3] = 0.0;   
                D[3][0] = 0.0;         D[3][1] = 0.0;              D[3][2] = 0.0;                  D[3][3] = 0.0;
            }
        }
        else {
            GRTRaiseError("Wrong execution of function %s.", __func__);
        }
    }
    
}

void grt_get_layer_D11(
    MYCOMPLEX xa, MYCOMPLEX xb, MYREAL k, MYCOMPLEX D[2][2])
{
    // 第iy层物理量
    if(xb != 1.0){
        D[0][0] = k;        D[0][1] = k*xb;
        D[1][0] = k*xa;     D[1][1] = k;   
    } else {
        D[0][0] = k;        D[0][1] = 0.0;
        D[1][0] = k*xa;     D[1][1] = 0.0;   
    }
    
}

void grt_get_layer_D12(
    MYCOMPLEX xa, MYCOMPLEX xb, MYREAL k, MYCOMPLEX D[2][2])
{
    // 第iy层物理量
    if(xb != 1.0){
        D[0][0] = k;        D[0][1] = -k*xb;
        D[1][0] = -k*xa;    D[1][1] = k;   
    } else {
        D[0][0] = k;        D[0][1] = 0.0;
        D[1][0] = -k*xa;    D[1][1] = 0.0;   
    }
    
}

void grt_get_layer_D11_uiz(
    MYCOMPLEX xa, MYCOMPLEX xb, MYREAL k, MYCOMPLEX D[2][2])
{
    // 第iy层物理量
    MYCOMPLEX a = k*xa;
    MYCOMPLEX b = k*xb;

    if(xb != 1.0){
        D[0][0] = a*k;     D[0][1] = b*b;
        D[1][0] = a*a;     D[1][1] = b*k;   
    } else {
        D[0][0] = a*k;     D[0][1] = 0.0;
        D[1][0] = a*a;     D[1][1] = 0.0;   
    }
}

void grt_get_layer_D12_uiz(
    MYCOMPLEX xa, MYCOMPLEX xb, MYREAL k, MYCOMPLEX D[2][2])
{
    // 第iy层物理量
    MYCOMPLEX a = k*xa;
    MYCOMPLEX b = k*xb;

    if(xb != 1.0){
        D[0][0] = - a*k;     D[0][1] = b*b;
        D[1][0] = a*a;       D[1][1] = - b*k;   
    } else {
        D[0][0] = - a*k;     D[0][1] = 0.0;
        D[1][0] = a*a;       D[1][1] = 0.0;   
    }
}

void grt_get_layer_D21(
    MYCOMPLEX xa, MYCOMPLEX xb, MYCOMPLEX kbkb, MYCOMPLEX mu,
    MYCOMPLEX omega, MYREAL rho, MYREAL k, MYCOMPLEX D[2][2])
{
    // 第iy层物理量
    MYCOMPLEX Omg;
    if(xb != 1.0){
        Omg = k*k - 0.5*kbkb;
        D[0][0] = 2*mu*Omg;        D[0][1] = 2*k*mu*k*xb;
        D[1][0] = 2*k*mu*k*xa;     D[1][1] = 2*mu*Omg;   
    } else {
        D[0][0] = - rho * GRT_SQUARE(omega);        D[0][1] = 0.0;
        D[1][0] = 0.0;                              D[1][1] = 0.0;   
    }
    
}

void grt_get_layer_D22(
    MYCOMPLEX xa, MYCOMPLEX xb, MYCOMPLEX kbkb, MYCOMPLEX mu,
    MYCOMPLEX omega, MYREAL rho, MYREAL k, MYCOMPLEX D[2][2])
{
    // 第iy层物理量
    MYCOMPLEX Omg;
    if(xb != 1.0){
        Omg = k*k - 0.5*kbkb;
        D[0][0] = 2*mu*Omg;        D[0][1] = -2*k*mu*k*xb;
        D[1][0] = -2*k*mu*k*xa;    D[1][1] = 2*mu*Omg;   
    } else {
        D[0][0] = - rho * GRT_SQUARE(omega);        D[0][1] = 0.0;
        D[1][0] = 0.0;                              D[1][1] = 0.0;   
    }
}

void grt_get_layer_T(
    MYCOMPLEX xb, MYCOMPLEX mu,
    MYCOMPLEX omega, MYREAL k, MYCOMPLEX T[2][2], bool inverse)
{
    // 液体层不应该使用该函数
    if(xb == 1.0){
        GRTRaiseError("Wrong execution of function %s.", __func__);
    }

    if( ! inverse ){
        T[0][0] = k;              T[0][1] = k;
        T[1][0] = mu*k*k*xb;      T[1][1] = - mu*k*k*xb;
    } else{
        T[0][0] = mu*k*xb;      T[0][1] = 1;
        T[1][0] = mu*k*xb;      T[1][1] = - 1;
        for(MYINT i=0; i<2; ++i){
            for(MYINT j=0; j<2; ++j){
                T[i][j] *= 1/(2*mu*k*k*xb);
            }
        }
    }
}

void grt_get_layer_E_Love(MYCOMPLEX xb1, MYREAL thk, MYREAL k, MYCOMPLEX E[2][2], bool inverse)
{
    MYCOMPLEX exb = exp(k*thk*xb1); 

    memset(E, 0, sizeof(MYCOMPLEX) * 4);
    if(! inverse){
        E[0][0] = exb;
        E[1][1] = 1.0/exb;
    } else {
        E[0][0] = 1.0/exb;
        E[1][1] = exb;
    }
    
}

void grt_get_layer_E_Rayl(
    MYCOMPLEX xa1, MYCOMPLEX xb1, MYREAL thk, MYREAL k, MYCOMPLEX E[4][4], bool inverse)
{
    MYCOMPLEX exa, exb; 

    exa = exp(k*thk*xa1);
    exb = exp(k*thk*xb1);

    memset(E, 0, sizeof(MYCOMPLEX) * 16);

    if( ! inverse){
        E[0][0] = exa;
        E[1][1] = exb;
        E[2][2] = 1.0/exa;
        E[3][3] = 1.0/exb;
    } else {
        E[0][0] = 1.0/exa;
        E[1][1] = 1.0/exb;
        E[2][2] = exa;
        E[3][3] = exb;
    }
}

void grt_RT_matrix_from_4x4(
    MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX kbkb1, MYCOMPLEX mu1, MYREAL rho1, 
    MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX kbkb2, MYCOMPLEX mu2, MYREAL rho2,
    MYCOMPLEX omega, MYREAL thk,
    MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats)
{
    MYCOMPLEX D1_inv[4][4], D2[4][4], Q[4][4];

    grt_get_layer_D(xa1, xb1, kbkb1, mu1, omega, rho1, k, D1_inv, true, 2);
    grt_get_layer_D(xa2, xb2, kbkb2, mu2, omega, rho2, k, D2,    false, 2);

    grt_cmatmxn_mul(4, 4, 4, D1_inv, D2, Q);

    MYCOMPLEX exa, exb; 

    exa = exp(-k*thk*xa1);
    exb = exp(-k*thk*xb1);

    MYCOMPLEX E[4][4] = {0};
    E[0][0] = exa;
    E[1][1] = exb;
    E[2][2] = 1/exa;
    E[3][3] = 1/exb;
    grt_cmatmxn_mul(4, 4, 4, E, Q, Q);

    // 对Q矩阵划分子矩阵 
    MYCOMPLEX Q11[2][2], Q12[2][2], Q21[2][2], Q22[2][2];
    grt_cmatmxn_block(4, 4, Q, 0, 0, 2, 2, Q11);
    grt_cmatmxn_block(4, 4, Q, 0, 2, 2, 2, Q12);
    grt_cmatmxn_block(4, 4, Q, 2, 0, 2, 2, Q21);
    grt_cmatmxn_block(4, 4, Q, 2, 2, 2, 2, Q22);

    // 计算反射透射系数 
    // TD
    grt_cmat2x2_inv(Q22, TD, stats);
    // RD
    grt_cmat2x2_mul(Q12, TD, RD); 
    // RU
    grt_cmat2x2_mul(TD, Q21, RU);
    grt_cmat2x2_k(RU, -1, RU);
    // TU
    grt_cmat2x2_mul(Q12, RU, TU);
    grt_cmat2x2_add(Q11, TU, TU);

    *RDL = (mu1*xb1 - mu2*xb2) / (mu1*xb1 + mu2*xb2) * exa*exa;
    *RUL = - (*RDL);
    *TDL = 2.0*mu1*xb1/(mu1*xb1 + mu2*xb2) * exb;
    *TUL = 2.0*mu2*xb2/(mu1*xb1 + mu2*xb2) * exb;

    
}