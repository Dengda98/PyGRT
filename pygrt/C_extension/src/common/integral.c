/**
 * @file   integral.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 *     将被积函数的逐点值累加成积分值
 *                   
 */


#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "common/integral.h"
#include "common/const.h"
#include "common/bessel.h"
#include "common/quadratic.h"



void int_Pk(
    MYREAL k, MYREAL r, 
    // F(ki,w)， 第一个维度表示不同震源，不同阶数，第二个维度3代表三类系数qm,wm,vm 
    const MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM],
    // F(ki,w)Jm(ki*r)ki，第一个维度表示不同震源，不同阶数，第二个维度代表4种类型的F(k,w)Jm(kr)k的类型
    bool calc_uir,
    MYCOMPLEX SUM[SRC_M_NUM][INTEG_NUM])
{
    MYREAL bjmk[MORDER_MAX+1] = {0};
    MYREAL kr = k*r;
    MYREAL kr_inv = RONE/kr;
    MYREAL kcoef = k;

    MYREAL Jmcoef[MORDER_MAX+1] = {0};

    bessel012(kr, &bjmk[0], &bjmk[1], &bjmk[2]); 
    if(calc_uir){
        MYREAL bjmk0[MORDER_MAX+1] = {0};
        for(MYINT i=0; i<=MORDER_MAX; ++i)  bjmk0[i] = bjmk[i];

        besselp012(kr, &bjmk[0], &bjmk[1], &bjmk[2]); 
        kcoef = k*k;

        for(MYINT i=1; i<=MORDER_MAX; ++i)  Jmcoef[i] = kr_inv * (-kr_inv * bjmk0[i] + bjmk[i]);
    } 
    else {
        for(MYINT i=1; i<=MORDER_MAX; ++i)  Jmcoef[i] = bjmk[i]*kr_inv;
    }

    for(MYINT i=1; i<=MORDER_MAX; ++i)  Jmcoef[i] *= kcoef;
    for(MYINT i=0; i<=MORDER_MAX; ++i)  bjmk[i] *= kcoef;


    // 公式(5.6.22), 将公式分解为F(k,w)Jm(kr)k的形式
    for(MYINT i=0; i<SRC_M_NUM; ++i){
        MYINT modr = SRC_M_ORDERS[i];  // 对应m阶数
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


void int_Pk_filon(
    MYREAL k, MYREAL r, bool iscos,
    const MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM],
    bool calc_uir,
    MYCOMPLEX SUM[SRC_M_NUM][INTEG_NUM])
{
    MYREAL phi0 = 0.0;
    if(! iscos)  phi0 = - HALFPI;  // 在cos函数中添加的相位差，用于计算sin函数

    MYREAL kr = k*r;
    MYREAL kcoef = sqrt(k);
    MYREAL bjmk[MORDER_MAX+1] = {0};

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

    for(MYINT i=0; i<=MORDER_MAX; ++i)  bjmk[i] *= kcoef;
    
    // 公式(5.6.22), 将公式分解为F(k,w)Jm(kr)k的形式，忽略近场项
    for(MYINT i=0; i<SRC_M_NUM; ++i){
        MYINT modr = SRC_M_ORDERS[i];  // 对应m阶数
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
static MYCOMPLEX interg_quad_cos(
    MYCOMPLEX a, MYCOMPLEX b, MYCOMPLEX c, MYINT modr, MYREAL r, MYREAL k1, MYREAL k2, MYINT kodr0)
{
    MYREAL s1, s2, c1, c2;
    MYREAL k1r, k2r, phi;
    k1r = k1*r;
    k2r = k2*r;
    phi = (2.0*modr+1.0)/4.0 * PI;
    s1 = sin(k1r-phi);   c1 = cos(k1r-phi);
    s2 = sin(k2r-phi);   c2 = cos(k2r-phi);

    MYCOMPLEX res;
    if(kodr0==0){
        res = + 2.0*a*(s1 - s2) / (r*r*r)
              + ( -a*k1*k1*s1 + a*k2*k2*s2 - b*k1*s1 + b*k2*s2 - c*s1 + c*s2 ) / r
              + ( -2.0*a*k1*c1 + 2.0*a*k2*c2 - b*c1 + b*c2 ) / (r*r);
    }
    else if(kodr0==1){
        MYREAL rr = r*r;
        MYREAL rrr = rr*r;
        MYREAL rrrr = rrr*r;
        MYREAL kk1 = k1*k1;
        MYREAL kk2 = k2*k2;
        MYREAL kkk1 = kk1*k1;
        MYREAL kkk2 = kk2*k2;
        res = + 6.0*a*(c1 - c2) / rrrr
              + ( -a*kkk1*s1 + a*kkk2*s2 -b*kk1*s1 + b*kk2*s2 - c*k1*s1 + c*k2*s2 ) / r
              + ( -3.0*a*kk1*c1 + 3.0*a*kk2*c2 - 2.0*b*k1*c1 + 2.0*b*k2*c2 - c*c1 + c*c2 ) / rr
              + 2.0*( 3.0*a*k1*s1 - 3.0*a*k2*s2 + b*s1 - b*s2 ) / rrr;
    }
    else{
        fprintf(stderr, "WRONG kodr0 in interg_quad_cos().\n");
        exit(EXIT_FAILURE);
    }

    return res;                
}



void int_Pk_sa_filon(
    const MYREAL k3[3], MYREAL r, 
    const MYCOMPLEX QWV3[3][SRC_M_NUM][QWV_NUM],
    bool calc_uir,
    MYCOMPLEX SUM[SRC_M_NUM][INTEG_NUM])
{
    // 使用bessel递推公式 Jm'(x) = m/x * Jm(x) - J_{m+1}(x)
    // 考虑大震中距，忽略第一项，再使用bessel渐近公式
    MYINT modr0 = (calc_uir)? 1 : 0;
    MYINT kodr0 = (calc_uir)? 1 : 0;
    MYINT sgn = (calc_uir)? -1 : 1;

    // 对sqrt(k)*F(k,w)进行二次曲线拟合，再计算 (a*k^2 + b*k + c) * cos(kr - (2m+1)/4) 的积分
    // 拟合二次函数的参数
    MYCOMPLEX quad_a[SRC_M_NUM][QWV_NUM]={0};
    MYCOMPLEX quad_b[SRC_M_NUM][QWV_NUM]={0};
    MYCOMPLEX quad_c[SRC_M_NUM][QWV_NUM]={0};
    for(MYINT im=0; im<SRC_M_NUM; ++im){
        MYINT modr = SRC_M_ORDERS[im];
        for(MYINT c=0; c<QWV_NUM; ++c){
            if(modr==0 && qwvchs[c] == 'v')  continue;

            MYCOMPLEX F3[3];
            for(MYINT d=0; d<3; ++d)  F3[d] = QWV3[d][im][c] * sqrt(k3[d]) * sgn;

            // 拟合参数
            quad_term(k3, F3, &quad_a[im][c], &quad_b[im][c], &quad_c[im][c]);
        }
    }

    MYREAL kmin = k3[0], kmax = k3[2];
    // 公式(5.6.22), 将公式分解为F(k,w)Jm(kr)k的形式
    for(MYINT im=0; im<SRC_M_NUM; ++im){
        MYINT modr = SRC_M_ORDERS[im];  // 对应m阶数
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



void merge_Pk(
    // F(ki,w)Jm(ki*r)ki，
    const MYCOMPLEX sum_J[SRC_M_NUM][INTEG_NUM], 
    // 累积求和，Z、R、T分量 
    MYCOMPLEX tol[SRC_M_NUM][CHANNEL_NUM])
{   
    for(MYINT i=0; i<SRC_M_NUM; ++i){
        MYINT modr = SRC_M_ORDERS[i];
        if(modr==0){
            tol[i][0] = sum_J[i][2];
            tol[i][1] = sum_J[i][0];
        }
        else {
            tol[i][0] = sum_J[i][2];
            tol[i][1] = sum_J[i][0] + sum_J[i][1];
            tol[i][2] = - sum_J[i][1] + sum_J[i][3];
        }
    }
}
