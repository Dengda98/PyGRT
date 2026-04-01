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
#include <string.h>

#include "grt/static/static_layer.h"
#include "grt/common/model.h"
#include "grt/common/matrix.h"

/* 定义用于提取相邻两层物性参数的宏 */
#define MODEL_2LAYS_ATTRIB(MM, T, N) \
    T N##1 = MM->N[iy-1];\
    T N##2 = MM->N[iy];\


static int __freeBound_R_PSV(real_t d, real_t k, cplx_t delta1, cplx_t M[2][2])
{
    if(d <= 0.0){
        // 公式(6.3.12)
        M[0][0] = M[1][1] = 0.0;
        M[0][1] = -delta1;
        M[1][0] = -1.0/delta1;
    } else {
        M[0][0] = M[1][1] = -2.0*d*k;
        M[0][1] = delta1 * (-4.0*d*d*k*k - 1.0);
        M[1][0] = -1.0/delta1;
    }
    return GRT_INVERSE_SUCCESS;
}

static int __freeBound_R_SH(cplx_t *ML)
{
    *ML = 1.0;
    return GRT_INVERSE_SUCCESS;
}

static int __rigidBound_R_PSV(real_t d, real_t k, cplx_t delta1, cplx_t M[2][2])
{
    if(d <= 0.0){
        // 公式(6.3.12)
        M[0][0] = M[1][1] = 0.0;
        M[0][1] = M[1][0] = 1.0;
    } else {
        M[0][0] = M[1][1] = 2.0*delta1*d*k;
        M[0][1] = (4.0*GRT_SQUARE(delta1*d*k) + 1.0);
        M[1][0] = 1.0;
    }
    return GRT_INVERSE_SUCCESS;
}

static int __rigidBound_R_SH(cplx_t *ML)
{
    *ML = -1.0;
    return GRT_INVERSE_SUCCESS;
}

void grt_static_topbound_RU_PSV(MODEL1D_STATE *mstat)
{
    MODEL1D *mod1d = mstat->mod1d;
    cplx_t delta = mstat->delta[0];
    real_t k = mstat->k;
    if(mod1d->topbound == GRT_BOUND_FREE){
        mstat->M_top.stats = __freeBound_R_PSV(0.0, k, delta, mstat->M_top.RU);
    }
    else if(mod1d->topbound == GRT_BOUND_RIGID){
        mstat->M_top.stats = __rigidBound_R_PSV(0.0, k, delta, mstat->M_top.RU);
    }
    else if(mod1d->topbound == GRT_BOUND_HALFSPACE){
        memset(mstat->M_top.RU, 0, sizeof(cplx_t)*4);
        mstat->M_top.stats = GRT_INVERSE_SUCCESS;
    }
    else{
        GRTRaiseError("Wrong execution.");
    }
    // RU 不需要时延
}

void grt_static_topbound_RU_SH(MODEL1D_STATE *mstat)
{
    MODEL1D *mod1d = mstat->mod1d;
    if(mod1d->topbound == GRT_BOUND_FREE){
        mstat->M_top.stats = __freeBound_R_SH(&mstat->M_top.RUL);
    }
    else if(mod1d->topbound == GRT_BOUND_RIGID){
        mstat->M_top.stats = __rigidBound_R_SH(&mstat->M_top.RUL);
    }
    else if(mod1d->topbound == GRT_BOUND_HALFSPACE){
        mstat->M_top.RUL = 0.0;
        mstat->M_top.stats = GRT_INVERSE_SUCCESS;
    }
    else{
        GRTRaiseError("Wrong execution.");
    }
    // RU 不需要时延
}

void grt_static_botbound_RD_PSV(MODEL1D_STATE *mstat)
{
    MODEL1D *mod1d = mstat->mod1d;
    size_t nlay = mod1d->n;
    cplx_t delta = mstat->delta[nlay-2];
    real_t thk = mod1d->Thk[nlay-2];
    real_t k = mstat->k;
    if(mod1d->botbound == GRT_BOUND_FREE){
        mstat->M_bot.stats = __freeBound_R_PSV(thk, k, delta, mstat->M_bot.RD);
    }
    else if(mod1d->botbound == GRT_BOUND_RIGID){
        mstat->M_bot.stats = __rigidBound_R_PSV(thk, k, delta, mstat->M_bot.RD);
    }
    else if(mod1d->botbound == GRT_BOUND_HALFSPACE){
        grt_static_RT_matrix_PSV(mstat, nlay-1, &mstat->M_bot);
    }
    else{
        GRTRaiseError("Wrong execution.");
    }

    // 时延 RD
    cplx_t ex, ex2;
    ex = exp(- k*thk);
    ex2 = ex * ex;

    mstat->M_bot.RD[0][0] *= ex2;   mstat->M_bot.RD[0][1] *= ex2;
    mstat->M_bot.RD[1][0] *= ex2;   mstat->M_bot.RD[1][1] *= ex2;
}

void grt_static_botbound_RD_SH(MODEL1D_STATE *mstat)
{
    MODEL1D *mod1d = mstat->mod1d;
    size_t nlay = mod1d->n;
    real_t thk = mod1d->Thk[nlay-2];
    real_t k = mstat->k;
    if(mod1d->botbound == GRT_BOUND_FREE){
        mstat->M_bot.stats = __freeBound_R_SH(&mstat->M_bot.RDL);
    }
    else if(mod1d->botbound == GRT_BOUND_RIGID){
        mstat->M_bot.stats = __rigidBound_R_SH(&mstat->M_bot.RDL);
    }
    else if(mod1d->botbound == GRT_BOUND_HALFSPACE){
        grt_static_RT_matrix_SH(mstat, nlay-1, &mstat->M_bot);
    }
    else{
        GRTRaiseError("Wrong execution.");
    }

    // 时延 RD
    cplx_t ex, ex2;
    ex = exp(- k*thk);
    ex2 = ex * ex;

    mstat->M_bot.RDL *= ex2;
}

void grt_static_wave2qwv_REV_PSV(MODEL1D_STATE *mstat)
{
    cplx_t D11[2][2] = {{1.0, -1.0}, {1.0, 1.0}};
    cplx_t D12[2][2] = {{1.0, -1.0}, {-1.0, -1.0}};

    // 公式(6.3.35,37)
    if(mstat->mod1d->ircvup){// 震源更深
        grt_cmat2x2_mul(D12, mstat->M_FA.RU, mstat->R_EV);
        grt_cmat2x2_add(D11, mstat->R_EV, mstat->R_EV);
    } else { // 接收点更深
        grt_cmat2x2_mul(D11, mstat->M_BL.RD, mstat->R_EV);
        grt_cmat2x2_add(D12, mstat->R_EV, mstat->R_EV);
    }
}

void grt_static_wave2qwv_REV_SH(MODEL1D_STATE *mstat)
{
    if(mstat->mod1d->ircvup){// 震源更深
        mstat->R_EVL = 1.0 + mstat->M_FA.RUL;
    } else {
        mstat->R_EVL = 1.0 + mstat->M_BL.RDL;
    }
}

void grt_static_wave2qwv_z_REV_PSV(MODEL1D_STATE *mstat)
{
    MODEL1D *mod1d = mstat->mod1d;
    real_t k = mstat->k;
    size_t ircv = mod1d->ircv;
    cplx_t delta1 = mstat->delta[ircv];

    // 新推导公式
    cplx_t kd2 = 2.0*k*delta1;
    cplx_t D11[2][2] = {{k, -k-kd2}, {k, k-kd2}};
    cplx_t D12[2][2] = {{-k, k+kd2}, {k, k-kd2}};
    if(mod1d->ircvup){// 震源更深
        grt_cmat2x2_mul(D12, mstat->M_FA.RU, mstat->uiz_R_EV);
        grt_cmat2x2_add(D11, mstat->uiz_R_EV, mstat->uiz_R_EV);
    } else { // 接收点更深
        grt_cmat2x2_mul(D11, mstat->M_BL.RD, mstat->uiz_R_EV);
        grt_cmat2x2_add(D12, mstat->uiz_R_EV, mstat->uiz_R_EV);
    }
}

void grt_static_wave2qwv_z_REV_SH(MODEL1D_STATE *mstat)
{
    real_t k = mstat->k;
    // 新推导公式
    if(mstat->mod1d->ircvup){// 震源更深
        mstat->uiz_R_EVL = (1.0 - mstat->M_FA.RUL)*k;
    } else { // 接收点更深
        mstat->uiz_R_EVL = (mstat->M_BL.RDL - 1.0)*k;
    }
}


void grt_static_RT_matrix_PSV(const MODEL1D_STATE *mstat, const size_t iy, RT_MATRIX *M)
{
    MODEL_2LAYS_ATTRIB(mstat, cplx_t, mu);
    MODEL_2LAYS_ATTRIB(mstat, cplx_t, delta);
    real_t thk = mstat->mod1d->Thk[iy-1];
    real_t k = mstat->k;

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


void grt_static_RT_matrix_SH(const MODEL1D_STATE *mstat, const size_t iy, RT_MATRIX *M)
{
    MODEL_2LAYS_ATTRIB(mstat, cplx_t, mu);
    
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


void grt_static_delay_RT_matrix(const MODEL1D_STATE *mstat, const size_t iy, RT_MATRIX *M)
{
    real_t thk = mstat->mod1d->Thk[iy-1];
    real_t k = mstat->k;
    
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