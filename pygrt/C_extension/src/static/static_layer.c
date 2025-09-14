/**
 * @file   static_layer.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 * 以下代码实现的是 P-SV 波和 SH 波的静态反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *              透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 *
 */

#include <stdio.h>
#include <complex.h>
#include <stdbool.h>

#include "grt/static/static_layer.h"
#include "grt/common/model.h"
#include "grt/common/matrix.h"

void grt_static_topfree_RU_PSV(MYCOMPLEX delta1, MYCOMPLEX R_tilt[2][2]){
    // 公式(6.3.12)
    R_tilt[0][0] = R_tilt[1][1] = 0.0;
    R_tilt[0][1] = -delta1;
    R_tilt[1][0] = -1.0/delta1;
}

void grt_static_wave2qwv_REV_PSV(bool ircvup, const MYCOMPLEX R[2][2], MYCOMPLEX R_EV[2][2])
{
    MYCOMPLEX D11[2][2] = {{1.0, -1.0}, {1.0, 1.0}};
    MYCOMPLEX D12[2][2] = {{1.0, -1.0}, {-1.0, -1.0}};

    // 公式(6.3.35,37)
    if(ircvup){// 震源更深
        grt_cmat2x2_mul(D12, R, R_EV);
        grt_cmat2x2_add(D11, R_EV, R_EV);
    } else { // 接收点更深
        grt_cmat2x2_mul(D11, R, R_EV);
        grt_cmat2x2_add(D12, R_EV, R_EV);
    }
}

void grt_static_wave2qwv_REV_SH(MYCOMPLEX RL, MYCOMPLEX *R_EVL)
{
    *R_EVL = (1.0 + (RL));
}

void grt_static_wave2qwv_z_REV_PSV(
    MYCOMPLEX delta1, bool ircvup, MYREAL k, 
    const MYCOMPLEX R[2][2], MYCOMPLEX R_EV[2][2])
{
    // 新推导公式
    MYCOMPLEX kd2 = 2.0*k*delta1;
    MYCOMPLEX D11[2][2] = {{k, -k-kd2}, {k, k-kd2}};
    MYCOMPLEX D12[2][2] = {{-k, k+kd2}, {k, k-kd2}};
    if(ircvup){// 震源更深
        grt_cmat2x2_mul(D12, R, R_EV);
        grt_cmat2x2_add(D11, R_EV, R_EV);
    } else { // 接收点更深
        grt_cmat2x2_mul(D11, R, R_EV);
        grt_cmat2x2_add(D12, R_EV, R_EV);
    }
}

void grt_static_wave2qwv_z_REV_SH(bool ircvup, MYREAL k, MYCOMPLEX RL, MYCOMPLEX *R_EVL)
{
    // 新推导公式
    if(ircvup){// 震源更深
        *R_EVL = (1.0 - (RL))*k;
    } else { // 接收点更深
        *R_EVL = (RL - 1.0)*k;
    }
}


void grt_static_RT_matrix_PSV(
    MYCOMPLEX delta1, MYCOMPLEX mu1, 
    MYCOMPLEX delta2, MYCOMPLEX mu2, 
    MYREAL thk, MYREAL k,
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2], MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2])
{
    // 公式(6.3.18)
    MYCOMPLEX dmu = mu1 - mu2;
    MYCOMPLEX A112 = mu1*delta1 + mu2;
    MYCOMPLEX A221 = mu2*delta2 + mu1;
    MYCOMPLEX B = mu1*delta1 - mu2*delta2;
    MYCOMPLEX del11 = delta1*delta1;
    MYREAL k2 = k*k;
    MYREAL thk2 = thk*thk;

    // Reflection
    //------------------ RD -----------------------------------
    RD[0][0] = -2.0*delta1*k*thk*dmu/A112;
    RD[0][1] = - ( 4.0*del11*k2*thk2*A221*dmu + A112*B ) / (A221*A112);
    RD[1][0] = - dmu/A112;
    RD[1][1] = RD[0][0];
    //------------------ RU -----------------------------------
    RU[0][0] = 0.0;
    RU[0][1] = B/A112;
    RU[1][0] = dmu/A221;
    RU[1][1] = 0.0;

    // Transmission
    //------------------ TD -----------------------------------
    TD[0][0] = mu1*(1.0+delta1)/(A112);
    TD[0][1] = 2.0*mu1*delta1*k*thk*(1.0+delta1)/(A112);
    TD[1][0] = 0.0;
    TD[1][1] = TD[0][0]*A112/A221;
    //------------------ TU -----------------------------------
    TU[0][0] = mu2*(1.0+delta2)/A221;
    TU[0][1] = 2.0*delta1*k*thk*mu2*(1.0+delta2)/A112;
    TU[1][0] = 0.0;
    TU[1][1] = TU[0][0]*A221/A112;
}


void grt_static_RT_matrix_SH(
    MYCOMPLEX mu1, MYCOMPLEX mu2, 
    MYREAL thk, MYREAL k,
    MYCOMPLEX *RDL, MYCOMPLEX *RUL, MYCOMPLEX *TDL, MYCOMPLEX *TUL)
{
    // 公式(6.3.18)
    MYCOMPLEX dmu = mu1 - mu2;
    MYCOMPLEX amu = mu1 + mu2;

    // Reflection
    *RDL = dmu/amu;
    *RUL = - dmu/amu;

    // Transmission
    *TDL = 2.0*mu1/amu;
    *TUL = (*TDL)*mu2/mu1;
}


void grt_static_delay_RT_matrix_PSV(
    MYREAL thk, MYREAL k,
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2], 
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2])
{
    MYCOMPLEX ex, ex2;
    ex = exp(- k*thk);
    ex2 = ex * ex;

    RD[0][0] *= ex2;   RD[0][1] *= ex2;
    RD[1][0] *= ex2;   RD[1][1] *= ex2;

    TD[0][0] *= ex;    TD[0][1] *= ex;
    TD[1][0] *= ex;    TD[1][1] *= ex;

    TU[0][0] *= ex;    TU[0][1] *= ex;
    TU[1][0] *= ex;    TU[1][1] *= ex;
}


void grt_static_delay_RT_matrix_SH(
    MYREAL thk, MYREAL k,
    MYCOMPLEX *RDL, MYCOMPLEX *RUL, 
    MYCOMPLEX *TDL, MYCOMPLEX *TUL)
{
    MYCOMPLEX ex, ex2;
    ex = exp(- k*thk);
    ex2 = ex * ex;

    *RDL *= ex2;
    *TDL *= ex;
    *TUL *= ex;
}