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

/* 定义用于提取相邻两层物性参数的宏 */
#define MODEL_2LAYS_ATTRIB(T, N) \
    T N##1 = mod1d->N[iy-1];\
    T N##2 = mod1d->N[iy];\


void grt_static_topfree_RU(GRT_MODEL1D *mod1d)
{
    cplx_t delta1 = mod1d->delta[0];
    // 公式(6.3.12)
    mod1d->M_top.RU[0][0] = mod1d->M_top.RU[1][1] = 0.0;
    mod1d->M_top.RU[0][1] = -delta1;
    mod1d->M_top.RU[1][0] = -1.0/delta1;
    mod1d->M_top.RUL = 1.0;
}

void grt_static_wave2qwv_REV_PSV(GRT_MODEL1D *mod1d)
{
    cplx_t D11[2][2] = {{1.0, -1.0}, {1.0, 1.0}};
    cplx_t D12[2][2] = {{1.0, -1.0}, {-1.0, -1.0}};

    // 公式(6.3.35,37)
    if(mod1d->ircvup){// 震源更深
        grt_cmat2x2_mul(D12, mod1d->M_FA.RU, mod1d->R_EV);
        grt_cmat2x2_add(D11, mod1d->R_EV, mod1d->R_EV);
    } else { // 接收点更深
        grt_cmat2x2_mul(D11, mod1d->M_BL.RD, mod1d->R_EV);
        grt_cmat2x2_add(D12, mod1d->R_EV, mod1d->R_EV);
    }
}

void grt_static_wave2qwv_REV_SH(GRT_MODEL1D *mod1d)
{
    if(mod1d->ircvup){// 震源更深
        mod1d->R_EVL = 1.0 + mod1d->M_FA.RUL;
    } else {
        mod1d->R_EVL = 1.0 + mod1d->M_BL.RDL;
    }
}

void grt_static_wave2qwv_z_REV_PSV(GRT_MODEL1D *mod1d)
{
    real_t k = mod1d->k;
    size_t ircv = mod1d->ircv;
    cplx_t delta1 = mod1d->delta[ircv];

    // 新推导公式
    cplx_t kd2 = 2.0*k*delta1;
    cplx_t D11[2][2] = {{k, -k-kd2}, {k, k-kd2}};
    cplx_t D12[2][2] = {{-k, k+kd2}, {k, k-kd2}};
    if(mod1d->ircvup){// 震源更深
        grt_cmat2x2_mul(D12, mod1d->M_FA.RU, mod1d->uiz_R_EV);
        grt_cmat2x2_add(D11, mod1d->uiz_R_EV, mod1d->uiz_R_EV);
    } else { // 接收点更深
        grt_cmat2x2_mul(D11, mod1d->M_BL.RD, mod1d->uiz_R_EV);
        grt_cmat2x2_add(D12, mod1d->uiz_R_EV, mod1d->uiz_R_EV);
    }
}

void grt_static_wave2qwv_z_REV_SH(GRT_MODEL1D *mod1d)
{
    real_t k = mod1d->k;
    // 新推导公式
    if(mod1d->ircvup){// 震源更深
        mod1d->uiz_R_EVL = (1.0 - mod1d->M_FA.RUL)*k;
    } else { // 接收点更深
        mod1d->uiz_R_EVL = (mod1d->M_BL.RDL - 1.0)*k;
    }
}


void grt_static_RT_matrix_PSV(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M)
{
    MODEL_2LAYS_ATTRIB(cplx_t, mu);
    MODEL_2LAYS_ATTRIB(cplx_t, delta);
    real_t thk = mod1d->Thk[iy-1];
    real_t k = mod1d->k;

    // 公式(6.3.18)
    cplx_t dmu = mu1 - mu2;
    cplx_t A112 = mu1*delta1 + mu2;
    cplx_t A221 = mu2*delta2 + mu1;
    cplx_t B = mu1*delta1 - mu2*delta2;
    cplx_t del11 = delta1*delta1;
    real_t k2 = k*k;
    real_t thk2 = thk*thk;

    // Reflection
    //------------------ RD -----------------------------------
    M->RD[0][0] = -2.0*delta1*k*thk*dmu/A112;
    M->RD[0][1] = - ( 4.0*del11*k2*thk2*A221*dmu + A112*B ) / (A221*A112);
    M->RD[1][0] = - dmu/A112;
    M->RD[1][1] = M->RD[0][0];
    //------------------ RU -----------------------------------
    M->RU[0][0] = 0.0;
    M->RU[0][1] = B/A112;
    M->RU[1][0] = dmu/A221;
    M->RU[1][1] = 0.0;

    // Transmission
    //------------------ TD -----------------------------------
    M->TD[0][0] = mu1*(1.0+delta1)/(A112);
    M->TD[0][1] = 2.0*mu1*delta1*k*thk*(1.0+delta1)/(A112);
    M->TD[1][0] = 0.0;
    M->TD[1][1] = M->TD[0][0]*A112/A221;
    //------------------ TU -----------------------------------
    M->TU[0][0] = mu2*(1.0+delta2)/A221;
    M->TU[0][1] = 2.0*delta1*k*thk*mu2*(1.0+delta2)/A112;
    M->TU[1][0] = 0.0;
    M->TU[1][1] = M->TU[0][0]*A221/A112;
}


void grt_static_RT_matrix_SH(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M)
{
    MODEL_2LAYS_ATTRIB(cplx_t, mu);
    
    // 公式(6.3.18)
    cplx_t dmu = mu1 - mu2;
    cplx_t amu = mu1 + mu2;

    // Reflection
    M->RDL = dmu/amu;
    M->RUL = - dmu/amu;

    // Transmission
    M->TDL = 2.0*mu1/amu;
    M->TUL = (M->TDL)*mu2/mu1;
}


void grt_static_delay_RT_matrix(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M)
{
    real_t thk = mod1d->Thk[iy-1];
    real_t k = mod1d->k;
    
    cplx_t ex, ex2;
    ex = exp(- k*thk);
    ex2 = ex * ex;

    M->RD[0][0] *= ex2;   M->RD[0][1] *= ex2;
    M->RD[1][0] *= ex2;   M->RD[1][1] *= ex2;

    M->TD[0][0] *= ex;    M->TD[0][1] *= ex;
    M->TD[1][0] *= ex;    M->TD[1][1] *= ex;

    M->TU[0][0] *= ex;    M->TU[0][1] *= ex;
    M->TU[1][0] *= ex;    M->TU[1][1] *= ex;

    M->RDL *= ex2;
    M->TDL *= ex;
    M->TUL *= ex;
}