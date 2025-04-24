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

#include "common/integral.h"
#include "common/const.h"
#include "common/bessel.h"



void int_Pk(
    MYREAL k, MYREAL r, 
    // F(ki,w)， 第一个维度3代表阶数m=0,1,2，第二个维度3代表三类系数qm,wm,vm 
    const MYCOMPLEX EXP_qwv[3][3], const MYCOMPLEX VF_qwv[3][3], 
    const MYCOMPLEX HF_qwv[3][3],  const MYCOMPLEX DC_qwv[3][3], 
    // F(ki,w)Jm(ki*r)ki，维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型
    bool calc_uir,
    MYCOMPLEX EXP_J[3][4], MYCOMPLEX VF_J[3][4], 
    MYCOMPLEX HF_J[3][4],  MYCOMPLEX DC_J[3][4])
{
    MYREAL bj0k, bj1k, bj2k;
    MYREAL kr = k*r;
    MYREAL kr_inv = RONE/kr;
    MYREAL kcoef = k;

    MYREAL J1coef, J2coef;

    bessel012(kr, &bj0k, &bj1k, &bj2k); 
    if(calc_uir){
        MYREAL j1, j2;
        j1 = bj1k;
        j2 = bj2k;
        besselp012(kr, &bj0k, &bj1k, &bj2k); 
        kcoef = k*k;

        J1coef = kr_inv * (-kr_inv * j1 + bj1k);
        J2coef = kr_inv * (-kr_inv * j2 + bj2k);
    } else {
        J1coef = bj1k*kr_inv;
        J2coef = bj2k*kr_inv;
    }

    J1coef *= kcoef;
    J2coef *= kcoef;

    bj0k *= kcoef;
    bj1k *= kcoef;
    bj2k *= kcoef;

    
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
    HF_J[1][1]  = - (HF_qwv[1][0] + HF_qwv[1][2])*J1coef;    // - (q1+v1)*J1*k/kr
    HF_J[1][2]  =   HF_qwv[1][1]*bj1k;         // w1*J1*k
    HF_J[1][3]  = - HF_qwv[1][2]*bj0k;         // -v1*J0*k
    }

    if(DC_qwv!=NULL){
    // m=0 剪切源
    DC_J[0][0] = - DC_qwv[0][0]*bj1k;
    DC_J[0][2] =   DC_qwv[0][1]*bj0k;

    // m=1 剪切源
    DC_J[1][0]  =   DC_qwv[1][0]*bj0k;         // q1*J0*k
    DC_J[1][1]  = - (DC_qwv[1][0] + DC_qwv[1][2])*J1coef;    // - (q1+v1)*J1*k/kr
    DC_J[1][2]  =   DC_qwv[1][1]*bj1k;         // w1*J1*k
    DC_J[1][3]  = - DC_qwv[1][2]*bj0k;         // -v1*J0*k

    // m=2 剪切源
    DC_J[2][0]  =   DC_qwv[2][0]*bj1k;         // q2*J1*k
    DC_J[2][1]  = - RTWO*(DC_qwv[2][0] + DC_qwv[2][2])*J2coef;    // - (q2+v2)*J2*k/kr
    DC_J[2][2]  =   DC_qwv[2][1]*bj2k;         // w2*J2*k
    DC_J[2][3]  = - DC_qwv[2][2]*bj1k;         // -v2*J1*k
    }
}


void int_Pk_filon(
    MYREAL k, MYREAL r, bool iscos,
    const MYCOMPLEX EXP_qwv[3][3], const MYCOMPLEX VF_qwv[3][3], 
    const MYCOMPLEX HF_qwv[3][3],  const MYCOMPLEX DC_qwv[3][3], 
    bool calc_uir,
    MYCOMPLEX EXP_J[3][4], MYCOMPLEX VF_J[3][4], 
    MYCOMPLEX HF_J[3][4],  MYCOMPLEX DC_J[3][4] )
{
    MYREAL phi0 = 0.0;
    if(! iscos)  phi0 = - HALFPI;  // 在cos函数中添加的相位差，用于计算sin函数

    MYREAL kr = k*r;
    MYREAL kr_inv = RONE/kr;
    MYREAL kcoef = SQRT(k);
    MYCOMPLEX bj0k, bj1k, bj2k;

    MYCOMPLEX J1coef, J2coef;

    if(calc_uir){
        kcoef *= k;
        // 使用bessel递推公式 Jm'(x) = m/x * Jm(x) - J_{m+1}(x)
        // 考虑大震中距，忽略第一项，再使用bessel渐近公式
        bj0k = - COS(kr - THREEQUARTERPI - phi0);
        bj1k = - COS(kr - FIVEQUARTERPI - phi0);
        bj2k = - COS(kr - SEVENQUARTERPI - phi0);
    } else {
        bj0k = COS(kr - QUARTERPI - phi0);
        bj1k = COS(kr - THREEQUARTERPI - phi0);
        bj2k = COS(kr - FIVEQUARTERPI - phi0);
    }
    J1coef = bj1k*kr_inv;
    J2coef = bj2k*kr_inv;

    J1coef *= kcoef;
    J2coef *= kcoef;

    bj0k *= kcoef;
    bj1k *= kcoef;
    bj2k *= kcoef;

    
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
    HF_J[1][1]  = - (HF_qwv[1][0] + HF_qwv[1][2])*J1coef;    // - (q1+v1)*J1*k/kr
    HF_J[1][2]  =   HF_qwv[1][1]*bj1k;         // w1*J1*k
    HF_J[1][3]  = - HF_qwv[1][2]*bj0k;         // -v1*J0*k
    }

    if(DC_qwv!=NULL){
    // m=0 剪切源
    DC_J[0][0] = - DC_qwv[0][0]*bj1k;
    DC_J[0][2] =   DC_qwv[0][1]*bj0k;

    // m=1 剪切源
    DC_J[1][0]  =   DC_qwv[1][0]*bj0k;         // q1*J0*k
    DC_J[1][1]  = - (DC_qwv[1][0] + DC_qwv[1][2])*J1coef;    // - (q1+v1)*J1*k/kr
    DC_J[1][2]  =   DC_qwv[1][1]*bj1k;         // w1*J1*k
    DC_J[1][3]  = - DC_qwv[1][2]*bj0k;         // -v1*J0*k

    // m=2 剪切源
    DC_J[2][0]  =   DC_qwv[2][0]*bj1k;         // q2*J1*k
    DC_J[2][1]  = - RTWO*(DC_qwv[2][0] + DC_qwv[2][2])*J2coef;    // - (q2+v2)*J2*k/kr
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
