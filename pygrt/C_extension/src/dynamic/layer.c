/**
 * @file   layer.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient 
 *               matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699
 * 
 */

#include <stdio.h>
#include <complex.h>

#include "dynamic/layer.h"
#include "common/model.h"
#include "common/prtdbg.h"
#include "common/matrix.h"



void calc_R_tilt(
    MYCOMPLEX xa0, MYCOMPLEX xb0, MYCOMPLEX kbkb0, MYREAL k, MYCOMPLEX R_tilt[2][2], MYINT *stats)
{   

    // // 公式(5.3.10-14)
    MYCOMPLEX Delta = RZERO;
    MYREAL kk = k*k; 
    MYCOMPLEX kbkb_k2inv = kbkb0/kk;
    MYCOMPLEX kbkb_k4inv = RQUART*kbkb_k2inv*kbkb_k2inv;

    // 对公式(5.3.10-14)进行重新整理，对浮点数友好一些
    Delta = -RONE + xa0*xb0 + kbkb_k2inv - kbkb_k4inv;
    if(Delta == CZERO){
        *stats = INVERSE_FAILURE;
        return;
    }
    R_tilt[0][0] = (RONE + xa0*xb0 - kbkb_k2inv + kbkb_k4inv) / Delta;
    R_tilt[0][1] = RTWO * xb0 * (RONE - RHALF*kbkb_k2inv) / Delta;
    R_tilt[1][0] = RTWO * xa0 * (RONE - RHALF*kbkb_k2inv) / Delta;
    R_tilt[1][1] = R_tilt[0][0];
}



void calc_R_EV(
    MYCOMPLEX xa_rcv, MYCOMPLEX xb_rcv, bool ircvup,
    MYREAL k, 
    const MYCOMPLEX R[2][2], MYCOMPLEX RL, 
    MYCOMPLEX R_EV[2][2], MYCOMPLEX *R_EVL)
{
    // 公式(5.2.19)
    MYCOMPLEX D11[2][2] = {{k, k*xb_rcv}, {k*xa_rcv, k}};
    MYCOMPLEX D12[2][2] = {{k, -k*xb_rcv}, {-k*xa_rcv, k}};

    // 公式(5.7.7,25)
    if(ircvup){// 震源更深
        cmat2x2_mul(D12, R, R_EV);
        cmat2x2_add(D11, R_EV, R_EV);
    } else { // 接收点更深
        cmat2x2_mul(D11, R, R_EV);
        cmat2x2_add(D12, R_EV, R_EV);
    }
    *R_EVL = (RONE + (RL))*k;

}


void calc_uiz_R_EV(
    MYCOMPLEX xa_rcv, MYCOMPLEX xb_rcv, bool ircvup,
    MYREAL k, 
    const MYCOMPLEX R[2][2], MYCOMPLEX RL, 
    MYCOMPLEX R_EV[2][2], MYCOMPLEX *R_EVL)
{
    // 将势函数转为ui,z在(B_m, P_m, C_m)系下的分量
    // 新推导的公式
    MYCOMPLEX ak = k*k*xa_rcv;
    MYCOMPLEX bk = k*k*xb_rcv;
    MYCOMPLEX bb = xb_rcv*bk;
    MYCOMPLEX aa = xa_rcv*ak;
    MYCOMPLEX D11[2][2] = {{ak, bb}, {aa, bk}};
    MYCOMPLEX D12[2][2] = {{-ak, bb}, {aa, -bk}};

    // 公式(5.7.7,25)
    if(ircvup){// 震源更深
        cmat2x2_mul(D12, R, R_EV);
        cmat2x2_add(D11, R_EV, R_EV);
        *R_EVL = (RONE - (RL))*bk;
    } else { // 接收点更深
        cmat2x2_mul(D11, R, R_EV);
        cmat2x2_add(D12, R_EV, R_EV);
        *R_EVL = (RL - RONE)*bk;
    }
    
}    
    


void calc_RT_2x2(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX kbkb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX kbkb2, MYCOMPLEX mu2, 
    MYREAL thk, // 使用上层的厚度
    MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats)
{
    
    MYCOMPLEX exa, exb, exab, ex2a, ex2b; 
    MYCOMPLEX tmp;

    exa = exp(-k*thk*xa1);
    exb = exp(-k*thk*xb1);

    exab = exa * exb;
    ex2a = exa * exa;
    ex2b = exb * exb;


    bool computeRayl = true;
    bool computeLove = true;
    if(RD==NULL || RU==NULL || TD==NULL || TU==NULL) computeRayl=false;
    if(RDL==NULL || RUL==NULL || TDL==NULL || TUL==NULL) computeLove=false;
    

    // 定义一些中间变量来简化运算和书写
    MYREAL kk = k*k;
    MYCOMPLEX dmu = mu1/mu2 - RONE; // mu1 - mu2; 分子分母同除mu2
    MYCOMPLEX dmu2 = dmu*dmu;

    MYCOMPLEX mu1kb1_k2 = mu1/mu2*kbkb1/kk;// mu1*kb1_k2;
    MYCOMPLEX mu2kb2_k2 = kbkb2/kk; // mu2*kb2_k2;

    MYREAL rho12 = Rho1 / Rho2;
    MYREAL rho21 = Rho2 / Rho1;

    // 从原公式上，分母包含5项，但前四项会随着k的增大迅速超过最后一项
    // 最后一项要小前几项10余个数量级，但计算结果还是保持在最后一项的量级，
    // 这种情况会受到浮点数的有效位数的限制，64bit的双精度double大概就是15-16位，
    // 故会发生严重精度损失的情况。目前只在实部上观察到这个现象，虚部基本都在相近量级(或许是相对不明显)
    // 
    // 以下对公式重新整理，提出k的高阶项，以避免上述问题
    MYCOMPLEX Delta;
    Delta =   dmu2*(RONE-xa1*xb1)*(RONE-xa2*xb2) + mu1kb1_k2*dmu*(rho21*(RONE-xa1*xb1) - (RONE-xa2*xb2)) 
            + RQUART*mu1kb1_k2*mu2kb2_k2*(rho12*(RONE-xa2*xb2) + rho21*(RONE-xa1*xb1) - RTWO - (xa1*xb2+xa2*xb1));

    if( Delta == CZERO ){
        // printf("# zero Delta_inv=%e+%eJ\n", creal(Delta_inv), cimag(Delta_inv));
        *stats = INVERSE_FAILURE;
        return;
    } 

    // REFELCTION
    if(computeRayl){
        //------------------ RD -----------------------------------
        // rpp+
        RD[0][0] = ( - dmu2*(RONE+xa1*xb1)*(RONE-xa2*xb2) - mu1kb1_k2*dmu*(rho21*(RONE+xa1*xb1) - (RONE-xa2*xb2))
                     - RQUART*mu1kb1_k2*mu2kb2_k2*(rho12*(RONE-xa2*xb2) + rho21*(RONE+xa1*xb1) - RTWO + (xa1*xb2-xa2*xb1))) / Delta * ex2a;
        // rsp+
        RD[0][1] = ( - dmu2*(RONE-xa2*xb2) + RHALF*mu1kb1_k2*dmu*((RONE-xa2*xb2) - RTWO*rho21) 
                     + RQUART*mu1kb1_k2*mu2kb2_k2*(RONE-rho21)) / Delta * (-RTWO*xb1) * exab;
        // rps+
        RD[1][0] = RD[0][1]*(xa1/xb1);
        // rss+
        RD[1][1] = ( - dmu2*(RONE+xa1*xb1)*(RONE-xa2*xb2) - mu1kb1_k2*dmu*(rho21*(RONE+xa1*xb1) - (RONE-xa2*xb2))
                     - RQUART*mu1kb1_k2*mu2kb2_k2*(rho12*(RONE-xa2*xb2) + rho21*(RONE+xa1*xb1) - RTWO - (xa1*xb2-xa2*xb1))) / Delta * ex2b;
        //------------------ RU -----------------------------------
        // rpp-
        RU[0][0] = ( - dmu2*(RONE-xa1*xb1)*(RONE+xa2*xb2) - mu1kb1_k2*dmu*(rho21*(RONE-xa1*xb1) - (RONE+xa2*xb2))
                     - RQUART*mu1kb1_k2*mu2kb2_k2*(rho12*(RONE+xa2*xb2) + rho21*(RONE-xa1*xb1) - RTWO - (xa1*xb2-xa2*xb1))) / Delta;
        // rsp-
        RU[0][1] = ( - dmu2*(RONE-xa1*xb1) - RHALF*mu1kb1_k2*dmu*(rho21*(RONE-xa1*xb1) - RTWO)
                     + RQUART*mu1kb1_k2*mu2kb2_k2*(RONE-rho12)) / Delta * (RTWO*xb2);
        // rps-
        RU[1][0] = RU[0][1]*(xa2/xb2);
        // rss-
        RU[1][1] = ( - dmu2*(RONE-xa1*xb1)*(RONE+xa2*xb2) - mu1kb1_k2*dmu*(rho21*(RONE-xa1*xb1) - (RONE+xa2*xb2))
                     - RQUART*mu1kb1_k2*mu2kb2_k2*(rho12*(RONE+xa2*xb2) + rho21*(RONE-xa1*xb1) - RTWO + (xa1*xb2-xa2*xb1))) / Delta;
    }
    if(computeLove){
        *RUL = (mu2*xb2 - mu1*xb1) / (mu2*xb2 + mu1*xb1) ;
        *RDL = - (*RUL) * ex2b;
    }
    


    // REFRACTION
    if(computeRayl){
        tmp = mu1kb1_k2*xa1*(dmu*(xb2-xb1) - RHALF*mu1kb1_k2*(rho21*xb1+xb2)) / Delta * exa;
        TD[0][0] = tmp;     TU[0][0] = (rho21*xa2/xa1) * tmp;
        tmp = mu1kb1_k2*xb1*(dmu*(RONE-xa1*xb2) - RHALF*mu1kb1_k2*(RONE-rho21)) / Delta * exb;
        TD[0][1] = tmp;     TU[1][0] = (rho21*xa2/xb1) * tmp;
        tmp = mu1kb1_k2*xa1*(dmu*(RONE-xa2*xb1) - RHALF*mu1kb1_k2*(RONE-rho21)) / Delta * exa;
        TD[1][0] = tmp;     TU[0][1] = (rho21*xb2/xa1) * tmp;
        tmp = mu1kb1_k2*xb1*(dmu*(xa2-xa1) - RHALF*mu1kb1_k2*(rho21*xa1+xa2)) / Delta * exb;
        TD[1][1] = tmp;     TU[1][1] = (rho21*xb2/xb1) * tmp;
    }
    if(computeLove){
        tmp = RTWO / (mu2*xb2 + mu1*xb1)  * exb;
        *TDL = mu1*xb1 * tmp;
        *TUL = mu2*xb2 * tmp;
    }
  
}



void get_layer_D(
    MYCOMPLEX xa, MYCOMPLEX xb, MYCOMPLEX kbkb, MYCOMPLEX mu,
    MYCOMPLEX omega, MYREAL k, MYCOMPLEX D[4][4], bool inverse)
{
    // 第iy层物理量
    MYCOMPLEX Omg;
    Omg = k*k - RHALF*kbkb;

    if( ! inverse ){
        D[0][0] = k;            D[0][1] = k*xb;             D[0][2] = k;            D[0][3] = -k*xb;     
        D[1][0] = k*xa;            D[1][1] = k;             D[1][2] = -k*xa;           D[1][3] = k;   
        D[2][0] = 2*mu*Omg;     D[2][1] = 2*k*mu*k*xb;      D[2][2] = 2*mu*Omg;     D[2][3] = -2*k*mu*k*xb;   
        D[3][0] = 2*k*mu*k*xa;     D[3][1] = 2*mu*Omg;      D[3][2] = -2*k*mu*k*xa;    D[3][3] = 2*mu*Omg;   
    } else {
        D[0][0] = -2*k*mu*k*xa*k*xb;  D[0][1] = 2*mu*Omg*k*xb;    D[0][2] = k*xa*k*xb;            D[0][3] = -k*k*xb;     
        D[1][0] = 2*mu*Omg*k*xa;   D[1][1] = -2*k*mu*k*xa*k*xb;   D[1][2] = -k*k*xa;           D[1][3] = k*xa*k*xb;   
        D[2][0] = -2*k*mu*k*xa*k*xb;  D[2][1] = -2*mu*Omg*k*xb;   D[2][2] = k*xa*k*xb;            D[2][3] = k*k*xb;   
        D[3][0] = -2*mu*Omg*k*xa;  D[3][1] = -2*k*mu*k*xa*k*xb;   D[3][2] = k*k*xa;            D[3][3] = k*xa*k*xb;
        for(MYINT i=0; i<4; ++i){
            for(MYINT j=0; j<4; ++j){
                D[i][j] *= (-1/(2*mu*kbkb*k*xa*k*xb));
            }
        }
    }
}




void calc_RT_2x2_from_4x4(
    MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX kbkb1, MYCOMPLEX mu1, 
    MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX kbkb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL thk,
    MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats)
{

    MYCOMPLEX D1_inv[4][4], D2[4][4], Q[4][4];

    get_layer_D(xa1, xb1, kbkb1, mu1, omega, k, D1_inv, true);
    get_layer_D(xa2, xb2, kbkb2, mu2, omega, k, D2,    false);

    cmatmxn_mul(4, 4, 4, D1_inv, D2, Q);

    MYCOMPLEX exa, exb; 

    exa = exp(-k*thk*xa1);
    exb = exp(-k*thk*xb1);

    MYCOMPLEX E[4][4] = {0};
    E[0][0] = exa;
    E[1][1] = exb;
    E[2][2] = 1/exa;
    E[3][3] = 1/exb;
    cmatmxn_mul(4, 4, 4, E, Q, Q);

    // fprintf(stderr, "Q\n");
    // cmatmxn_print(4, 4, Q);

    // 对Q矩阵划分子矩阵 
    MYCOMPLEX Q11[2][2], Q12[2][2], Q21[2][2], Q22[2][2];
    cmatmxn_block(4, 4, Q, 0, 0, 2, 2, Q11);
    cmatmxn_block(4, 4, Q, 0, 2, 2, 2, Q12);
    cmatmxn_block(4, 4, Q, 2, 0, 2, 2, Q21);
    cmatmxn_block(4, 4, Q, 2, 2, 2, 2, Q22);

    // 计算反射透射系数 
    // TD
    cmat2x2_inv(Q22, TD, stats);
    // RD
    cmat2x2_mul(Q12, TD, RD); 
    // RU
    cmat2x2_mul(TD, Q21, RU);
    cmat2x2_k(RU, -1, RU);
    // TU
    cmat2x2_mul(Q12, RU, TU);
    cmat2x2_add(Q11, TU, TU);

    *RDL = (mu1*xb1 - mu2*xb2) / (mu1*xb1 + mu2*xb2) * exa*exa;
    *RUL = - (*RDL);
    *TDL = RTWO*mu1*xb1/(mu1*xb1 + mu2*xb2) * exb;
    *TUL = RTWO*mu2*xb2/(mu1*xb1 + mu2*xb2) * exb;

    
}