/**
 * @file   static_layer.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 * 以下代码实现的是 静态反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *              透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 *
 */

#include <stdio.h>
#include <complex.h>
#include <stdbool.h>

#include "static/static_layer.h"
#include "common/model.h"
#include "common/matrix.h"

void calc_static_R_tilt(MYCOMPLEX delta1, MYCOMPLEX R_tilt[2][2]){
    // 公式(6.3.12)
    R_tilt[0][0] = R_tilt[1][1] = CZERO;
    R_tilt[0][1] = -delta1;
    R_tilt[1][0] = -RONE/delta1;
}

void calc_static_R_EV(
    bool ircvup,
    const MYCOMPLEX R[2][2], MYCOMPLEX RL, 
    MYCOMPLEX R_EV[2][2], MYCOMPLEX *R_EVL)
{
    MYCOMPLEX D11[2][2] = {{RONE, -RONE}, {RONE, RONE}};
    MYCOMPLEX D12[2][2] = {{RONE, -RONE}, {-RONE, -RONE}};

    // 公式(6.3.35,37)
    if(ircvup){// 震源更深
        cmat2x2_mul(D12, R, R_EV);
        cmat2x2_add(D11, R_EV, R_EV);
    } else { // 接收点更深
        cmat2x2_mul(D11, R, R_EV);
        cmat2x2_add(D12, R_EV, R_EV);
    }
    *R_EVL = (RONE + (RL));
}

void calc_static_uiz_R_EV(
    MYCOMPLEX delta1, bool ircvup, MYREAL k, 
    const MYCOMPLEX R[2][2], MYCOMPLEX RL, 
    MYCOMPLEX R_EV[2][2], MYCOMPLEX *R_EVL)
{
    // 新推导公式
    MYCOMPLEX kd2 = RTWO*k*delta1;
    MYCOMPLEX D11[2][2] = {{k, -k-kd2}, {k, k-kd2}};
    MYCOMPLEX D12[2][2] = {{-k, k+kd2}, {k, k-kd2}};
    if(ircvup){// 震源更深
        cmat2x2_mul(D12, R, R_EV);
        cmat2x2_add(D11, R_EV, R_EV);
        *R_EVL = (RONE - (RL))*k;
    } else { // 接收点更深
        cmat2x2_mul(D11, R, R_EV);
        cmat2x2_add(D12, R_EV, R_EV);
        *R_EVL = (RL - RONE)*k;
    }
}


void calc_static_RT_2x2(
    MYCOMPLEX delta1, MYCOMPLEX mu1, 
    MYCOMPLEX delta2, MYCOMPLEX mu2, 
    MYREAL thk, MYREAL k,
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL)
{
    // 公式(6.3.18)

    MYCOMPLEX ex, ex2; 

    ex = exp(-k*thk);
    ex2 = ex*ex;

    MYCOMPLEX dmu = mu1 - mu2;
    MYCOMPLEX amu = mu1 + mu2;
    MYCOMPLEX A112 = mu1*delta1 + mu2;
    MYCOMPLEX A221 = mu2*delta2 + mu1;
    MYCOMPLEX B = mu1*delta1 - mu2*delta2;
    MYCOMPLEX del11 = delta1*delta1;
    MYREAL k2 = k*k;
    MYREAL thk2 = thk*thk;

    // REFELCTION
    //------------------ RD -----------------------------------
    RD[0][0] = -RTWO*delta1*k*thk*dmu/A112 * ex2;
    RD[0][1] = - ( RFOUR*del11*k2*thk2*A221*dmu + A112*B ) / (A221*A112) * ex2;
    RD[1][0] = - dmu/A112 * ex2;
    RD[1][1] = RD[0][0];
    //------------------ RU -----------------------------------
    RU[0][0] = RZERO;
    RU[0][1] = B/A112;
    RU[1][0] = dmu/A221;
    RU[1][1] = RZERO;

    *RDL = dmu/amu * ex2;
    *RUL = - dmu/amu;

    // Transmission
    //------------------ TD -----------------------------------
    TD[0][0] = mu1*(RONE+delta1)/(A112) * ex;
    TD[0][1] = RTWO*mu1*delta1*k*thk*(RONE+delta1)/(A112) * ex;
    TD[1][0] = RZERO;
    TD[1][1] = TD[0][0]*A112/A221;
    //------------------ TU -----------------------------------
    TU[0][0] = mu2*(RONE+delta2)/A221 * ex;
    TU[0][1] = RTWO*delta1*k*thk*mu2*(RONE+delta2)/A112 * ex;
    TU[1][0] = RZERO;
    TU[1][1] = TU[0][0]*A221/A112;

    *TDL = RTWO*mu1/amu * ex;
    *TUL = (*TDL)*mu2/mu1;

    // printf("delta1=%.5e%+.5ej, delta2=%.5e%+.5ej, mu1=%.5e%+.5ej, mu2=%.5e%+.5ej, thk=%e, k=%e\n", 
    //         creal(delta1),cimag(delta1),creal(delta2),cimag(delta2),creal(mu1),cimag(mu1),creal(mu2),cimag(mu2),
    //         thk, k);
    // cmatmxn_print(2, 2, RD);
    // cmatmxn_print(2, 2, RU);
    // cmatmxn_print(2, 2, TD);
    // cmatmxn_print(2, 2, TU);
    // printf("-----------------------------\n");
}