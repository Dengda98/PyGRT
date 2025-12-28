/**
 * @file   k_integ.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 *     将被积函数的逐点值累加成积分值
 *                   
 */


#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "grt/common/const.h"
#include "grt/common/bessel.h"
#include "grt/integral/quadratic.h"
#include "grt/integral/k_integ.h"


K_INTEG * grt_init_K_INTEG(const bool calc_upar, const size_t nr)
{
    K_INTEG *K = (K_INTEG *)calloc(1, sizeof(K_INTEG));
    K->calc_upar = calc_upar;
    K->nr = nr;
    K->sumJ  = (cplxIntegGrid *)calloc(nr, sizeof(cplxIntegGrid));
    K->sumJz = K->sumJr = NULL;
    if(calc_upar){
        K->sumJz = (cplxIntegGrid *)calloc(nr, sizeof(cplxIntegGrid));
        K->sumJr = (cplxIntegGrid *)calloc(nr, sizeof(cplxIntegGrid));
    }
    return K;
}

K_INTEG * grt_copy_K_INTEG(const K_INTEG *K)
{
    K_INTEG *K0 = (K_INTEG *)calloc(1, sizeof(K_INTEG));
    // 先直接赋值，实现浅拷贝
    *K0 = *K;

    // 对指针部分再重新申请内存并赋值，实现深拷贝
    size_t nr = K->nr;
    K0->sumJ  = (cplxIntegGrid *)calloc(nr, sizeof(cplxIntegGrid));
    memcpy(K0->sumJ, K->sumJ, sizeof(cplxIntegGrid)*nr);
    K0->sumJz = K0->sumJr = NULL;
    if(K->calc_upar){
        K0->sumJz = (cplxIntegGrid *)calloc(nr, sizeof(cplxIntegGrid));
        memcpy(K0->sumJz, K->sumJz, sizeof(cplxIntegGrid)*nr);
        K0->sumJr = (cplxIntegGrid *)calloc(nr, sizeof(cplxIntegGrid));
        memcpy(K0->sumJr, K->sumJr, sizeof(cplxIntegGrid)*nr);
    }
    return K0;
}

void grt_free_K_INTEG(K_INTEG *K)
{
    GRT_SAFE_FREE_PTR(K->sumJ);
    GRT_SAFE_FREE_PTR(K->sumJz);
    GRT_SAFE_FREE_PTR(K->sumJr);
    GRT_SAFE_FREE_PTR(K);
}


void grt_int_Pk(real_t k, real_t r, const cplxChnlGrid QWV, bool calc_uir, cplxIntegGrid SUM)
{
    real_t bjmk[GRT_MORDER_MAX+1] = {0};
    real_t kr = k*r;
    real_t kr_inv = 1.0/kr;
    real_t kcoef = k;

    real_t Jmcoef[GRT_MORDER_MAX+1] = {0};

    grt_bessel012(kr, &bjmk[0], &bjmk[1], &bjmk[2]); 
    if(calc_uir){
        real_t bjmk0[GRT_MORDER_MAX+1] = {0};
        for(int i=0; i<=GRT_MORDER_MAX; ++i)  bjmk0[i] = bjmk[i];

        grt_besselp012(kr, &bjmk[0], &bjmk[1], &bjmk[2]); 
        kcoef = k*k;

        for(int i=1; i<=GRT_MORDER_MAX; ++i)  Jmcoef[i] = kr_inv * (-kr_inv * bjmk0[i] + bjmk[i]);
    } 
    else {
        for(int i=1; i<=GRT_MORDER_MAX; ++i)  Jmcoef[i] = bjmk[i]*kr_inv;
    }

    for(int i=1; i<=GRT_MORDER_MAX; ++i)  Jmcoef[i] *= kcoef;
    for(int i=0; i<=GRT_MORDER_MAX; ++i)  bjmk[i] *= kcoef;


    // 公式(5.6.22), 将公式分解为F(k,w)Jm(kr)k的形式
    for(int i=0; i<GRT_SRC_M_NUM; ++i){
        int modr = GRT_SRC_M_ORDERS[i];  // 对应m阶数
        if(modr == 0){
            SUM[i][0] = - QWV[i][0] * bjmk[1];   // - q0*J1*k
            SUM[i][2] =   QWV[i][1] * bjmk[0];   //   w0*J0*k
        }
        else{
            SUM[i][0]  =   QWV[i][0] * bjmk[modr-1];         // qm*Jm-1*k
            SUM[i][1]  = - modr*(QWV[i][0] + QWV[i][2]) * Jmcoef[modr];    // - m*(qm+vm)*Jm*k/kr
            SUM[i][2]  =   QWV[i][1] * bjmk[modr];           // wm*Jm*k
            SUM[i][3]  = - QWV[i][2] * bjmk[modr-1];         // -vm*Jm-1*k
        }
    }

}


void grt_int_Pk_filon(real_t k, real_t r, bool iscos, const cplxChnlGrid QWV, bool calc_uir, cplxIntegGrid SUM)
{
    real_t phi0 = 0.0;
    if(! iscos)  phi0 = - HALFPI;  // 在cos函数中添加的相位差，用于计算sin函数

    real_t kr = k*r;
    real_t kcoef = sqrt(k);
    real_t bjmk[GRT_MORDER_MAX+1] = {0};

    if(calc_uir){
        kcoef *= k;
        // 使用bessel递推公式 Jm'(x) = m/x * Jm(x) - J_{m+1}(x)
        // 考虑大震中距，忽略第一项，再使用bessel渐近公式
        bjmk[0] = - cos(kr - THREEQUARTERPI - phi0);
        bjmk[1] = - cos(kr - FIVEQUARTERPI - phi0);
        bjmk[2] = - cos(kr - SEVENQUARTERPI - phi0);
    } else {
        bjmk[0] = cos(kr - QUARTERPI - phi0);
        bjmk[1] = cos(kr - THREEQUARTERPI - phi0);
        bjmk[2] = cos(kr - FIVEQUARTERPI - phi0);
    }

    for(int i=0; i<=GRT_MORDER_MAX; ++i)  bjmk[i] *= kcoef;
    
    // 公式(5.6.22), 将公式分解为F(k,w)Jm(kr)k的形式，忽略近场项
    for(int i=0; i<GRT_SRC_M_NUM; ++i){
        int modr = GRT_SRC_M_ORDERS[i];  // 对应m阶数
        if(modr == 0){
            SUM[i][0] = - QWV[i][0] * bjmk[1];   // - q0*J1*k
            SUM[i][2] =   QWV[i][1] * bjmk[0];   //   w0*J0*k
        }
        else{
            SUM[i][0]  =   QWV[i][0] * bjmk[modr-1];         // qm*Jm-1*k
            SUM[i][1]  =   0.0;    // - m*(qm+vm)*Jm*k/kr
            SUM[i][2]  =   QWV[i][1] * bjmk[modr];           // wm*Jm*k
            SUM[i][3]  = - QWV[i][2] * bjmk[modr-1];         // -vm*Jm-1*k
        }
    }
}


/**
 * 计算 (a*k^2 + b*k + c) * cos(kr - (2m+1)/4) 在[k1, k2]上的积分，
 * 如果kodr0==1，则说明对r取偏导，需计算(a*k^3 + b*k^2 + c*k) * cos(kr - (2m+1)/4) 在[k1, k2]上的积分，
 */
static cplx_t interg_quad_cos(
    cplx_t a, cplx_t b, cplx_t c, int modr, real_t r, real_t k1, real_t k2, int kodr0)
{
    real_t s1, s2, c1, c2;
    real_t k1r, k2r, phi;
    k1r = k1*r;
    k2r = k2*r;
    phi = (2.0*modr+1.0)/4.0 * PI;
    s1 = sin(k1r-phi);   c1 = cos(k1r-phi);
    s2 = sin(k2r-phi);   c2 = cos(k2r-phi);

    cplx_t res;
    if(kodr0==0){
        res = + 2.0*a*(s1 - s2) / (r*r*r)
              + ( -a*k1*k1*s1 + a*k2*k2*s2 - b*k1*s1 + b*k2*s2 - c*s1 + c*s2 ) / r
              + ( -2.0*a*k1*c1 + 2.0*a*k2*c2 - b*c1 + b*c2 ) / (r*r);
    }
    else if(kodr0==1){
        real_t rr = r*r;
        real_t rrr = rr*r;
        real_t rrrr = rrr*r;
        real_t kk1 = k1*k1;
        real_t kk2 = k2*k2;
        real_t kkk1 = kk1*k1;
        real_t kkk2 = kk2*k2;
        res = + 6.0*a*(c1 - c2) / rrrr
              + ( -a*kkk1*s1 + a*kkk2*s2 -b*kk1*s1 + b*kk2*s2 - c*k1*s1 + c*k2*s2 ) / r
              + ( -3.0*a*kk1*c1 + 3.0*a*kk2*c2 - 2.0*b*k1*c1 + 2.0*b*k2*c2 - c*c1 + c*c2 ) / rr
              + 2.0*( 3.0*a*k1*s1 - 3.0*a*k2*s2 + b*s1 - b*s2 ) / rrr;
    }
    else{
        GRTRaiseError("WRONG kodr0 in interg_quad_cos().\n");
    }

    return res;                
}



void grt_int_Pk_sa_filon(const real_t k3[3], real_t r, const cplxChnlGrid QWV3[3], bool calc_uir, cplxIntegGrid SUM)
{
    // 使用bessel递推公式 Jm'(x) = m/x * Jm(x) - J_{m+1}(x)
    // 考虑大震中距，忽略第一项，再使用bessel渐近公式
    int modr0 = (calc_uir)? 1 : 0;
    int kodr0 = (calc_uir)? 1 : 0;
    int sgn = (calc_uir)? -1 : 1;

    // 对sqrt(k)*F(k,w)进行二次曲线拟合，再计算 (a*k^2 + b*k + c) * cos(kr - (2m+1)/4) 的积分
    // 拟合二次函数的参数
    cplxChnlGrid quad_a={0};
    cplxChnlGrid quad_b={0};
    cplxChnlGrid quad_c={0};
    GRT_LOOP_ChnlGrid(im, c){
        int modr = GRT_SRC_M_ORDERS[im];
        if(modr==0 && GRT_QWV_CODES[c] == 'v')  continue;
        cplx_t F3[3];
        for(int d=0; d<3; ++d)  F3[d] = QWV3[d][im][c] * sqrt(k3[d]) * sgn;

        // 拟合参数
        grt_quad_term(k3, F3, &quad_a[im][c], &quad_b[im][c], &quad_c[im][c]);
    }

    real_t kmin = k3[0], kmax = k3[2];
    // 公式(5.6.22), 将公式分解为F(k,w)Jm(kr)k的形式
    for(int im=0; im<GRT_SRC_M_NUM; ++im){
        int modr = GRT_SRC_M_ORDERS[im];  // 对应m阶数
        if(modr == 0){
            SUM[im][0] = - interg_quad_cos(quad_a[im][0],quad_b[im][0],quad_c[im][0],1+modr0,r,kmin,kmax,kodr0);   // - q0*J1*k
            SUM[im][2] =   interg_quad_cos(quad_a[im][1],quad_b[im][1],quad_c[im][1],0+modr0,r,kmin,kmax,kodr0);   //   w0*J0*k
        }
        else{
            SUM[im][0]  =   interg_quad_cos(quad_a[im][0],quad_b[im][0],quad_c[im][0],modr-1+modr0,r,kmin,kmax,kodr0);         // qm*Jm-1*k
            SUM[im][1]  =   0.0;    // - m*(qm+vm)*Jm*k/kr
            SUM[im][2]  =   interg_quad_cos(quad_a[im][1],quad_b[im][1],quad_c[im][1],modr+modr0,r,kmin,kmax,kodr0);           // wm*Jm*k
            SUM[im][3]  = - interg_quad_cos(quad_a[im][2],quad_b[im][2],quad_c[im][2],modr-1+modr0,r,kmin,kmax,kodr0);;         // -vm*Jm-1*k
        }
    }
}



void grt_merge_Pk(const cplxIntegGrid sumJ, cplxChnlGrid tol)
{   
    for(int i=0; i<GRT_SRC_M_NUM; ++i){
        int modr = GRT_SRC_M_ORDERS[i];
        if(modr==0){
            tol[i][0] = sumJ[i][2];
            tol[i][1] = sumJ[i][0];
        }
        else {
            tol[i][0] = sumJ[i][2];
            tol[i][1] = sumJ[i][0] + sumJ[i][1];
            tol[i][2] = - sumJ[i][1] + sumJ[i][3];
        }
    }
}
