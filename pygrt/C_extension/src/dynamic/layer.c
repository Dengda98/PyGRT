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

/* 定义用于提取相邻两层物性参数的宏 */
#define MODEL_2LAYS_ATTRIB(T, N) \
    T N##1 = mod1d->N[iy-1];\
    T N##2 = mod1d->N[iy];\


void grt_topfree_RU(GRT_MODEL1D *mod1d)
{
    cplx_t cbcb0 = mod1d->cbcb[0];
    if(cbcb0 != 0.0){
        // 固体表面
        // 公式(5.3.10-14)
        cplx_t Delta = 0.0;
        cplx_t cbcb02 = 0.25*cbcb0*cbcb0;

        cplx_t xa0 = mod1d->xa[0];
        cplx_t xb0 = mod1d->xb[0];

        // 对公式(5.3.10-14)进行重新整理，对浮点数友好一些
        Delta = -1.0 + xa0*xb0 + cbcb0 - cbcb02;
        if(Delta == 0.0){
            mod1d->M_top.stats = GRT_INVERSE_FAILURE;
            return;
        }
        Delta = 1.0 / Delta;
        mod1d->M_top.RU[0][0] = (1.0 + xa0*xb0 - cbcb0 + cbcb02) * Delta;
        mod1d->M_top.RU[0][1] = 2.0 * xb0 * (1.0 - 0.5*cbcb0) * Delta;
        mod1d->M_top.RU[1][0] = 2.0 * xa0 * (1.0 - 0.5*cbcb0) * Delta;
        mod1d->M_top.RU[1][1] = mod1d->M_top.RU[0][0];
        mod1d->M_top.RUL = 1.0;
    }
    else {
        // 液体表面
        mod1d->M_top.RU[0][0] = -1.0;
        mod1d->M_top.RU[1][1] = mod1d->M_top.RU[0][1] = mod1d->M_top.RU[1][0] = 0.0;
        mod1d->M_top.RUL = 0.0;
    }
}


void grt_wave2qwv_REV_PSV(GRT_MODEL1D *mod1d)
{
    size_t ircv = mod1d->ircv;
    cplx_t xa = mod1d->xa[ircv];
    cplx_t xb = mod1d->xb[ircv];
    real_t k = mod1d->k;

    cplx_t D11[2][2], D12[2][2];
    if(xb != 1.0){
        // 位于固体层
        // 公式(5.2.19)
        D11[0][0] = k;         D11[0][1] = k*xb;
        D11[1][0] = k*xa;      D11[1][1] = k;
        D12[0][0] = k;         D12[0][1] = -k*xb;
        D12[1][0] = -k*xa;     D12[1][1] = k;
    } else {
        // 位于液体层
        D11[0][0] = k;         D11[0][1] = 0.0;
        D11[1][0] = k*xa;      D11[1][1] = 0.0;
        D12[0][0] = k;         D12[0][1] = 0.0;
        D12[1][0] = -k*xa;     D12[1][1] = 0.0;
    }

    // 公式(5.7.7,25)
    if(mod1d->ircvup){// 震源更深
        grt_cmat2x2_mul(D12, mod1d->M_FA.RU, mod1d->R_EV);
        grt_cmat2x2_add(D11, mod1d->R_EV, mod1d->R_EV);
    } else { // 接收点更深
        grt_cmat2x2_mul(D11, mod1d->M_BL.RD, mod1d->R_EV);
        grt_cmat2x2_add(D12, mod1d->R_EV, mod1d->R_EV);
    }
}


void grt_wave2qwv_REV_SH(GRT_MODEL1D *mod1d)
{
    size_t ircv = mod1d->ircv;
    cplx_t xb = mod1d->xb[ircv];
    real_t k = mod1d->k;

    if(xb != 1.0){
        // 位于固体层
        // 公式(5.2.19)
        if(mod1d->ircvup){// 震源更深
            mod1d->R_EVL = (1.0 + mod1d->M_FA.RUL)*k;
        } else {
            mod1d->R_EVL = (1.0 + mod1d->M_BL.RDL)*k;
        }
    } else {
        // 位于液体层
        mod1d->R_EVL = 0.0;
    }
}


void grt_wave2qwv_z_REV_PSV(GRT_MODEL1D *mod1d)
{
    size_t ircv = mod1d->ircv;
    cplx_t xa = mod1d->xa[ircv];
    cplx_t xb = mod1d->xb[ircv];
    real_t k = mod1d->k;


    // 将垂直波函数转为ui,z在(B_m, P_m, C_m)系下的分量
    // 新推导的公式
    cplx_t ak = k*k*xa;
    cplx_t bk = k*k*xb;
    cplx_t bb = xb*bk;
    cplx_t aa = xa*ak;
    cplx_t D11[2][2] = {{ak, bb}, {aa, bk}};
    cplx_t D12[2][2] = {{-ak, bb}, {aa, -bk}};

    // 位于液体层
    if(xb == 1.0){
        D11[0][1] = D11[1][1] = D12[0][1] = D12[1][1] = 0.0;
    }

    // 公式(5.7.7,25)
    if(mod1d->ircvup){// 震源更深
        grt_cmat2x2_mul(D12, mod1d->M_FA.RU, mod1d->uiz_R_EV);
        grt_cmat2x2_add(D11, mod1d->uiz_R_EV, mod1d->uiz_R_EV);
    } else { // 接收点更深
        grt_cmat2x2_mul(D11, mod1d->M_BL.RD, mod1d->uiz_R_EV);
        grt_cmat2x2_add(D12, mod1d->uiz_R_EV, mod1d->uiz_R_EV);
    }
}    


void grt_wave2qwv_z_REV_SH(GRT_MODEL1D *mod1d)
{
    size_t ircv = mod1d->ircv;
    cplx_t xb = mod1d->xb[ircv];
    real_t k = mod1d->k;
    
    // 将垂直波函数转为ui,z在(B_m, P_m, C_m)系下的分量
    // 新推导的公式
    cplx_t bk = k*k*xb;

    if(xb != 1.0){
        // 位于固体层
        if(mod1d->ircvup){// 震源更深
            mod1d->uiz_R_EVL = (1.0 - mod1d->M_FA.RUL)*bk;
        } else { // 接收点更深
            mod1d->uiz_R_EVL = (mod1d->M_BL.RDL - 1.0)*bk;
        }
    } else {
        // 位于液体层
        mod1d->uiz_R_EVL = 0.0;
    }
}    


void grt_RT_matrix_ll_PSV(const GRT_MODEL1D *mod1d, size_t iy, RT_MATRIX *M)
{
    MODEL_2LAYS_ATTRIB(cplx_t, xa);
    MODEL_2LAYS_ATTRIB(real_t, Rho);

    cplx_t A = xa1*Rho2 + xa2*Rho1;
    if(A==0.0){
        M->stats = GRT_INVERSE_FAILURE;
        return;
    }

    M->RD[0][0] = (xa1*Rho2 - xa2*Rho1)/A;  
    M->RD[0][1] = M->RD[1][0] = M->RD[1][1] = 0.0;
    
    M->RU[0][0] = (xa2*Rho1 - xa1*Rho2)/A;
    M->RU[0][1] = M->RU[1][0] = M->RU[1][1] = 0.0;

    M->TD[0][0] = 2.0*xa1*Rho1/A;
    M->TD[0][1] = M->TD[1][0] = M->TD[1][1] = 0.0;

    M->TU[0][0] = 2.0*xa2*Rho2/A;
    M->TU[0][1] = M->TU[1][0] = M->TU[1][1] = 0.0;

}

void grt_RT_matrix_ll_SH(RT_MATRIX *M)
{
    M->RDL = 0.0;
    M->RUL = 0.0;
    M->TDL = 0.0;
    M->TUL = 0.0;
}



void grt_RT_matrix_ls_PSV(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M)
{
    MODEL_2LAYS_ATTRIB(cplx_t, xa);
    MODEL_2LAYS_ATTRIB(cplx_t, xb);
    MODEL_2LAYS_ATTRIB(cplx_t, mu);
    MODEL_2LAYS_ATTRIB(cplx_t, cbcb);
    MODEL_2LAYS_ATTRIB(real_t, Rho);

    // 后缀1表示上层的液体的物理参数，后缀2表示下层的固体的物理参数
    // 若mu2==0, 则下层为液体，参数需相互交换 

    // 讨论液-固 or 固-液
    bool isfluidUp = (mu1 == 0.0);  // 上层是否为液体
    int sgn = 1;
    if(isfluidUp && mu2 == 0.0){
        GRTRaiseError("fluid-fluid interface is not allowed.\n");
    }

    // 使用指针
    cplx_t (*pRD)[2], (*pRU)[2];
    cplx_t (*pTD)[2], (*pTU)[2];
    if(isfluidUp){
        pRD = M->RD; pRU = M->RU; 
        pTD = M->TD; pTU = M->TU; 
    } else {
        pRD = M->RU; pRU = M->RD;
        pTD = M->TU; pTU = M->TD;
        GRT_SWAP(real_t, Rho1, Rho2);
        GRT_SWAP(cplx_t, xa1, xa2);
        GRT_SWAP(cplx_t, xb1, xb2);
        GRT_SWAP(cplx_t, cbcb1, cbcb2);
        GRT_SWAP(cplx_t, mu1, mu2);
        sgn = -1;
    }

    
    // 定义一些中间变量来简化运算和书写
    cplx_t lamka1k = Rho1*GRT_SQUARE(mod1d->c_phase);
    cplx_t kb2k = cbcb2;
    cplx_t Og2k = 1.0 - 0.5*kb2k;
    cplx_t Og2k2 = Og2k*Og2k;
    cplx_t A = 2.0*Og2k2*xa1*mu2 + 0.5*lamka1k*kb2k*xa2 - 2.0*mu2*xa1*xa2*xb2;
    cplx_t B = 2.0*Og2k2*xa1*mu2 - 0.5*lamka1k*kb2k*xa2 + 2.0*mu2*xa1*xa2*xb2;
    cplx_t C = 2.0*Og2k2*xa1*mu2 + 0.5*lamka1k*kb2k*xa2 + 2.0*mu2*xa1*xa2*xb2;
    cplx_t D = 2.0*Og2k2*xa1*mu2 - 0.5*lamka1k*kb2k*xa2 - 2.0*mu2*xa1*xa2*xb2;

    if(A == 0.0){
        M->stats = GRT_INVERSE_FAILURE;
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


void grt_RT_matrix_ls_SH(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M)
{
    // TEMPORARY!!!!!!
    // 之后不再使用mu或其它变量来判断是否为液体，直接定义一个新的数组
    MODEL_2LAYS_ATTRIB(cplx_t, mu);
    
    // 后缀1表示上层的液体的物理参数，后缀2表示下层的固体的物理参数
    // 若mu2==0, 则下层为液体，参数需相互交换 

    // 讨论液-固 or 固-液
    bool isfluidUp = (mu1 == 0.0);  // 上层是否为液体
    if(isfluidUp && mu2 == 0.0){
        GRTRaiseError("fluid-fluid interface is not allowed.\n");
    }

    // 使用指针
    cplx_t *pRDL, *pRUL;
    cplx_t *pTDL, *pTUL;
    if(isfluidUp){
        pRDL = &M->RDL; pRUL = &M->RUL;
        pTDL = &M->TDL; pTUL = &M->TUL;
    } else {
        pRDL = &M->RUL; pRUL = &M->RDL;
        pTDL = &M->TUL; pTUL = &M->TDL;
    }

    *pRDL = 0.0;
    *pRUL = 1.0;
    *pTDL = 0.0;
    *pTUL = 0.0;
}



void grt_RT_matrix_ss_PSV(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M)
{
    MODEL_2LAYS_ATTRIB(cplx_t, xa);
    MODEL_2LAYS_ATTRIB(cplx_t, xb);
    MODEL_2LAYS_ATTRIB(cplx_t, mu);
    MODEL_2LAYS_ATTRIB(cplx_t, cbcb);
    MODEL_2LAYS_ATTRIB(real_t, Rho);

    // 定义一些中间变量来简化运算和书写
    // real_t kk = k*k;
    cplx_t rmu = mu1/mu2;
    cplx_t dmu = rmu - 1.0; // mu1 - mu2; 分子分母同除mu2
    cplx_t dmu2 = dmu*dmu;

    cplx_t mu1cbcb1 = rmu*cbcb1;// mu1*kb1_k2;
    cplx_t mu2cbcb2 = cbcb2; // mu2*kb2_k2;

    real_t rho12 = Rho1 / Rho2;
    real_t rho21 = Rho2 / Rho1;

    // 从原公式上，分母包含5项，但前四项会随着k的增大迅速超过最后一项
    // 最后一项要小前几项10余个数量级，但计算结果还是保持在最后一项的量级，
    // 这种情况会受到浮点数的有效位数的限制，64bit的双精度double大概就是15-16位，
    // 故会发生严重精度损失的情况。目前只在实部上观察到这个现象，虚部基本都在相近量级(或许是相对不明显)
    // 
    // 以下对公式重新整理，提出k的高阶项，以避免上述问题
    cplx_t Delta;
    Delta =   dmu2*(1.0-xa1*xb1)*(1.0-xa2*xb2) + mu1cbcb1*dmu*(rho21*(1.0-xa1*xb1) - (1.0-xa2*xb2)) 
            + 0.25*mu1cbcb1*mu2cbcb2*(rho12*(1.0-xa2*xb2) + rho21*(1.0-xa1*xb1) - 2.0 - (xa1*xb2+xa2*xb1));

    if( Delta == 0.0 ){
        // printf("# zero Delta_inv=%e+%eJ\n", creal(Delta_inv), cimag(Delta_inv));
        M->stats = GRT_INVERSE_FAILURE;
        return;
    } 

    Delta = 1.0 / Delta;

    cplx_t tmp;

    // REFELCTION
    //------------------ RD -----------------------------------
    // rpp+
    M->RD[0][0] = ( - dmu2*(1.0+xa1*xb1)*(1.0-xa2*xb2) - mu1cbcb1*dmu*(rho21*(1.0+xa1*xb1) - (1.0-xa2*xb2))
                    - 0.25*mu1cbcb1*mu2cbcb2*(rho12*(1.0-xa2*xb2) + rho21*(1.0+xa1*xb1) - 2.0 + (xa1*xb2-xa2*xb1))) * Delta;
    // rsp+
    tmp = ( - dmu2*(1.0-xa2*xb2) + 0.5*mu1cbcb1*dmu*((1.0-xa2*xb2) - 2.0*rho21) 
                    + 0.25*mu1cbcb1*mu2cbcb2*(1.0-rho21)) * Delta * (-2.0);
    M->RD[0][1] = xb1*tmp;
    // rps+
    M->RD[1][0] = xa1*tmp;
    // rss+
    M->RD[1][1] = ( - dmu2*(1.0+xa1*xb1)*(1.0-xa2*xb2) - mu1cbcb1*dmu*(rho21*(1.0+xa1*xb1) - (1.0-xa2*xb2))
                    - 0.25*mu1cbcb1*mu2cbcb2*(rho12*(1.0-xa2*xb2) + rho21*(1.0+xa1*xb1) - 2.0 - (xa1*xb2-xa2*xb1))) * Delta;
    //------------------ RU -----------------------------------
    // rpp-
    M->RU[0][0] = ( - dmu2*(1.0-xa1*xb1)*(1.0+xa2*xb2) - mu1cbcb1*dmu*(rho21*(1.0-xa1*xb1) - (1.0+xa2*xb2))
                    - 0.25*mu1cbcb1*mu2cbcb2*(rho12*(1.0+xa2*xb2) + rho21*(1.0-xa1*xb1) - 2.0 - (xa1*xb2-xa2*xb1))) * Delta;
    // rsp-
    tmp = ( - dmu2*(1.0-xa1*xb1) - 0.5*mu1cbcb1*dmu*(rho21*(1.0-xa1*xb1) - 2.0)
                    + 0.25*mu1cbcb1*mu2cbcb2*(1.0-rho12)) * Delta * (2.0);
    M->RU[0][1] = xb2*tmp;
    // rps-
    M->RU[1][0] = xa2*tmp;
    // rss-
    M->RU[1][1] = ( - dmu2*(1.0-xa1*xb1)*(1.0+xa2*xb2) - mu1cbcb1*dmu*(rho21*(1.0-xa1*xb1) - (1.0+xa2*xb2))
                    - 0.25*mu1cbcb1*mu2cbcb2*(rho12*(1.0+xa2*xb2) + rho21*(1.0-xa1*xb1) - 2.0 + (xa1*xb2-xa2*xb1))) * Delta;

    // REFRACTION
    tmp = mu1cbcb1*(dmu*(xb2-xb1) - 0.5*mu1cbcb1*(rho21*xb1+xb2)) * Delta;
    M->TD[0][0] = xa1*tmp;     M->TU[0][0] = (rho21*xa2) * tmp;
    tmp = mu1cbcb1*(dmu*(1.0-xa1*xb2) - 0.5*mu1cbcb1*(1.0-rho21)) * Delta;
    M->TD[0][1] = xb1*tmp;     M->TU[1][0] = (rho21*xa2) * tmp;
    tmp = mu1cbcb1*(dmu*(1.0-xa2*xb1) - 0.5*mu1cbcb1*(1.0-rho21)) * Delta;
    M->TD[1][0] = xa1*tmp;     M->TU[0][1] = (rho21*xb2) * tmp;
    tmp = mu1cbcb1*(dmu*(xa2-xa1) - 0.5*mu1cbcb1*(rho21*xa1+xa2)) * Delta;
    M->TD[1][1] = xb1*tmp;     M->TU[1][1] = (rho21*xb2) * tmp;
}


void grt_RT_matrix_ss_SH(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M)
{
    MODEL_2LAYS_ATTRIB(cplx_t, xb);
    MODEL_2LAYS_ATTRIB(cplx_t, mu);

    // REFELCTION
    M->RUL = (mu2*xb2 - mu1*xb1) / (mu2*xb2 + mu1*xb1) ;
    M->RDL = - (M->RUL);

    // REFRACTION
    cplx_t tmp;
    tmp = 2.0 / (mu2*xb2 + mu1*xb1);
    M->TDL = mu1*xb1 * tmp;
    M->TUL = mu2*xb2 * tmp;
}



void grt_RT_matrix_PSV(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M)
{
    // TEMPORARY!!!!!!
    // 之后不再使用mu或其它变量来判断是否为液体，直接定义一个新的数组
    MODEL_2LAYS_ATTRIB(cplx_t, mu);
    
    // 根据界面两侧的具体情况选择函数
    if(mu1 != 0.0 && mu2 != 0.0){
        grt_RT_matrix_ss_PSV(mod1d, iy, M);
    }
    else if(mu1 == 0.0 && mu2 == 0.0){
        grt_RT_matrix_ll_PSV(mod1d, iy, M);
    }
    else{
        grt_RT_matrix_ls_PSV(mod1d, iy, M);
    }
}


void grt_RT_matrix_SH(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M)
{
    // TEMPORARY!!!!!!
    // 之后不再使用mu或其它变量来判断是否为液体，直接定义一个新的数组
    MODEL_2LAYS_ATTRIB(cplx_t, mu);
    
    // 根据界面两侧的具体情况选择函数
    if(mu1 != 0.0 && mu2 != 0.0){
        grt_RT_matrix_ss_SH(mod1d, iy, M);
    }
    else if(mu1 == 0.0 && mu2 == 0.0){
        grt_RT_matrix_ll_SH(M);
    }
    else{
        grt_RT_matrix_ls_SH(mod1d, iy, M);
    }
}


void grt_delay_RT_matrix(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M)
{
    real_t thk = mod1d->Thk[iy-1];
    cplx_t xa1 = mod1d->xa[iy-1];
    cplx_t xb1 = mod1d->xb[iy-1];
    real_t k = mod1d->k;

    cplx_t exa, exb, ex2a, ex2b, exab;
    exa = exp(- k*thk*xa1);
    exb = exp(- k*thk*xb1);
    // exa = 0.9;
    // exb = 0.9;
    ex2a = exa * exa;
    ex2b = exb * exb;
    exab = exa * exb;

    M->RD[0][0] *= ex2a;   M->RD[0][1] *= exab;
    M->RD[1][0] *= exab;   M->RD[1][1] *= ex2b;

    M->TD[0][0] *= exa;    M->TD[0][1] *= exb;
    M->TD[1][0] *= exa;    M->TD[1][1] *= exb;

    M->TU[0][0] *= exa;    M->TU[0][1] *= exa;
    M->TU[1][0] *= exb;    M->TU[1][1] *= exb;

    M->RDL *= ex2b;
    M->TDL *= exb;
    M->TUL *= exb;
}



void grt_delay_GRT_matrix(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M)
{
    real_t thk = mod1d->Thk[iy-1];
    cplx_t xa1 = mod1d->xa[iy-1];
    cplx_t xb1 = mod1d->xb[iy-1];
    real_t k = mod1d->k;
    
    cplx_t exa, exb, ex2a, ex2b, exab;
    exa = exp(- k*thk*xa1);
    exb = exp(- k*thk*xb1);
    ex2a = exa * exa;
    ex2b = exb * exb;
    exab = exa * exb;

    M->RU[0][0] *= ex2a;    M->RU[0][1] *= exab;  
    M->RU[1][0] *= exab;    M->RU[1][1] *= ex2b;  
    
    M->TD[0][0] *= exa;     M->TD[0][1] *= exa; 
    M->TD[1][0] *= exb;     M->TD[1][1] *= exb;

    M->TU[0][0] *= exa;     M->TU[0][1] *= exb; 
    M->TU[1][0] *= exa;     M->TU[1][1] *= exb;

    M->RUL *= ex2b;
    M->TDL *= exb;
    M->TUL *= exb;
}


void grt_get_layer_D(
    cplx_t xa, cplx_t xb, cplx_t kbkb, cplx_t mu,
    cplx_t omega, real_t rho, real_t k, cplx_t D[4][4], bool inverse, int liquid_invtype)
{
    // 第iy层物理量
    cplx_t Omg;
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
            for(int i=0; i<4; ++i){
                for(int j=0; j<4; ++j){
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
                for(int i=0; i<4; ++i){
                    for(int j=0; j<4; ++j){
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
            GRTRaiseError("Wrong execution.");
        }
    }
    
}

void grt_get_layer_D11(
    cplx_t xa, cplx_t xb, real_t k, cplx_t D[2][2])
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
    cplx_t xa, cplx_t xb, real_t k, cplx_t D[2][2])
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
    cplx_t xa, cplx_t xb, real_t k, cplx_t D[2][2])
{
    // 第iy层物理量
    cplx_t a = k*xa;
    cplx_t b = k*xb;

    if(xb != 1.0){
        D[0][0] = a*k;     D[0][1] = b*b;
        D[1][0] = a*a;     D[1][1] = b*k;   
    } else {
        D[0][0] = a*k;     D[0][1] = 0.0;
        D[1][0] = a*a;     D[1][1] = 0.0;   
    }
}

void grt_get_layer_D12_uiz(
    cplx_t xa, cplx_t xb, real_t k, cplx_t D[2][2])
{
    // 第iy层物理量
    cplx_t a = k*xa;
    cplx_t b = k*xb;

    if(xb != 1.0){
        D[0][0] = - a*k;     D[0][1] = b*b;
        D[1][0] = a*a;       D[1][1] = - b*k;   
    } else {
        D[0][0] = - a*k;     D[0][1] = 0.0;
        D[1][0] = a*a;       D[1][1] = 0.0;   
    }
}

void grt_get_layer_D21(
    cplx_t xa, cplx_t xb, cplx_t kbkb, cplx_t mu,
    cplx_t omega, real_t rho, real_t k, cplx_t D[2][2])
{
    // 第iy层物理量
    cplx_t Omg;
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
    cplx_t xa, cplx_t xb, cplx_t kbkb, cplx_t mu,
    cplx_t omega, real_t rho, real_t k, cplx_t D[2][2])
{
    // 第iy层物理量
    cplx_t Omg;
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
    cplx_t xb, cplx_t mu,
    cplx_t omega, real_t k, cplx_t T[2][2], bool inverse)
{
    // 液体层不应该使用该函数
    if(xb == 1.0){
        GRTRaiseError("Wrong execution.");
    }

    if( ! inverse ){
        T[0][0] = k;              T[0][1] = k;
        T[1][0] = mu*k*k*xb;      T[1][1] = - mu*k*k*xb;
    } else{
        T[0][0] = mu*k*xb;      T[0][1] = 1;
        T[1][0] = mu*k*xb;      T[1][1] = - 1;
        for(int i=0; i<2; ++i){
            for(int j=0; j<2; ++j){
                T[i][j] *= 1/(2*mu*k*k*xb);
            }
        }
    }
}

void grt_get_layer_E_Love(cplx_t xb1, real_t thk, real_t k, cplx_t E[2][2], bool inverse)
{
    cplx_t exb = exp(k*thk*xb1); 

    memset(E, 0, sizeof(cplx_t) * 4);
    if(! inverse){
        E[0][0] = exb;
        E[1][1] = 1.0/exb;
    } else {
        E[0][0] = 1.0/exb;
        E[1][1] = exb;
    }
    
}

void grt_get_layer_E_Rayl(
    cplx_t xa1, cplx_t xb1, real_t thk, real_t k, cplx_t E[4][4], bool inverse)
{
    cplx_t exa, exb; 

    exa = exp(k*thk*xa1);
    exb = exp(k*thk*xb1);

    memset(E, 0, sizeof(cplx_t) * 16);

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
    cplx_t xa1, cplx_t xb1, cplx_t kbkb1, cplx_t mu1, real_t rho1, 
    cplx_t xa2, cplx_t xb2, cplx_t kbkb2, cplx_t mu2, real_t rho2,
    cplx_t omega, real_t thk,
    real_t k, 
    cplx_t RD[2][2], cplx_t *RDL, cplx_t RU[2][2], cplx_t *RUL, 
    cplx_t TD[2][2], cplx_t *TDL, cplx_t TU[2][2], cplx_t *TUL, int *stats)
{
    cplx_t D1_inv[4][4], D2[4][4], Q[4][4];

    grt_get_layer_D(xa1, xb1, kbkb1, mu1, omega, rho1, k, D1_inv, true, 2);
    grt_get_layer_D(xa2, xb2, kbkb2, mu2, omega, rho2, k, D2,    false, 2);

    grt_cmatmxn_mul(4, 4, 4, D1_inv, D2, Q);

    cplx_t exa, exb; 

    exa = exp(-k*thk*xa1);
    exb = exp(-k*thk*xb1);

    cplx_t E[4][4] = {0};
    E[0][0] = exa;
    E[1][1] = exb;
    E[2][2] = 1/exa;
    E[3][3] = 1/exb;
    grt_cmatmxn_mul(4, 4, 4, E, Q, Q);

    // 对Q矩阵划分子矩阵 
    cplx_t Q11[2][2], Q12[2][2], Q21[2][2], Q22[2][2];
    grt_cmatmxn_block(4, 4, Q, 0, 0, 2, 2, Q11);
    grt_cmatmxn_block(4, 4, Q, 0, 2, 2, 2, Q12);
    grt_cmatmxn_block(4, 4, Q, 2, 0, 2, 2, Q21);
    grt_cmatmxn_block(4, 4, Q, 2, 2, 2, 2, Q22);

    // 计算反射透射系数 
    // TD
    grt_cmat2x2_inv(Q22, TD);
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