/**
 * @file   recursion.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 * 以下代码通过递推公式计算两层的广义反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "common/recursion.h"
#include "common/const.h"
#include "common/matrix.h"



void recursion_RD(
    const MYCOMPLEX RD1[2][2], MYCOMPLEX RDL1, const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats)
{
    MYCOMPLEX tmp1[2][2], tmp2[2][2], inv1;

    // RD, RDL
    cmat2x2_mul(RU1, RD2, tmp1);
    cmat2x2_one_sub(tmp1);
    cmat2x2_inv(tmp1, tmp1, stats);  if(*stats==INVERSE_FAILURE)  return;
    cmat2x2_mul(tmp1, TD1, tmp2);
    if(inv_2x2T!=NULL) cmat2x2_assign(tmp2, inv_2x2T);

    cmat2x2_mul(RD2, tmp2, tmp1);
    cmat2x2_mul(TU1, tmp1, tmp2);
    cmat2x2_add(RD1, tmp2, RD);
    inv1 = RONE / (RONE - RUL1*RDL2) * TDL1;
    *RDL = RDL1 + TUL1*RDL2*inv1;
    if(invT!=NULL)  *invT = inv1;
}


void recursion_TD(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, 
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, 
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats)
{
    MYCOMPLEX tmp1[2][2], tmp2[2][2], inv1;

    // TD, TDL
    cmat2x2_mul(RU1, RD2, tmp2);
    cmat2x2_one_sub(tmp2);
    cmat2x2_inv(tmp2, tmp1, stats);  if(*stats==INVERSE_FAILURE)  return;
    cmat2x2_mul(tmp1, TD1, tmp2);
    if(inv_2x2T!=NULL)  cmat2x2_assign(tmp2, inv_2x2T);
    cmat2x2_mul(TD2, tmp2, TD);
    
    inv1 = RONE / (RONE - RUL1*RDL2) * TDL1;
    *TDL = TDL2 * inv1;

    if(invT!=NULL) *invT = inv1;
}


void recursion_RU(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, const MYCOMPLEX RU2[2][2], MYCOMPLEX RUL2,
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats)
{
    MYCOMPLEX tmp1[2][2], tmp2[2][2], inv1;

    // RU, RUL
    cmat2x2_mul(RD2, RU1, tmp2);
    cmat2x2_one_sub(tmp2);
    cmat2x2_inv(tmp2, tmp1, stats);  if(*stats==INVERSE_FAILURE)  return;
    cmat2x2_mul(tmp1, TU2, tmp2);
    if(inv_2x2T!=NULL)  cmat2x2_assign(tmp2, inv_2x2T);

    cmat2x2_mul(RU1, tmp2, tmp1); 
    cmat2x2_mul(TD2, tmp1, tmp2);
    cmat2x2_add(RU2, tmp2, RU);
    inv1 = RONE / (RONE - RUL1*RDL2) * TUL2;
    *RUL = RUL2 + TDL2*RUL1*inv1; 

    if(invT!=NULL)  *invT = inv1;
}


void recursion_TU(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2,
    const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats)
{
    MYCOMPLEX tmp1[2][2], tmp2[2][2], inv1;

    // TU, TUL
    cmat2x2_mul(RD2, RU1, tmp2);
    cmat2x2_one_sub(tmp2);
    cmat2x2_inv(tmp2, tmp1, stats);  if(*stats==INVERSE_FAILURE)  return;
    cmat2x2_mul(tmp1, TU2, tmp2);
    if(inv_2x2T!=NULL) cmat2x2_assign(tmp2, inv_2x2T);
    cmat2x2_mul(TU1, tmp2, TU);
    
    inv1 = RONE / (RONE - RUL1*RDL2) * TUL2;
    *TUL = TUL1 * inv1;

    if(invT!=NULL)  *invT = inv1;

}




void recursion_RT_2x2(
    const MYCOMPLEX RD1[2][2], MYCOMPLEX RDL1, const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, const MYCOMPLEX RU2[2][2], MYCOMPLEX RUL2,
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL,
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats)
{

    // 临时矩阵
    MYCOMPLEX tmp1[2][2], tmp2[2][2];
    MYCOMPLEX inv0, inv1T;

    inv0 = RONE / (RONE - RUL1*RDL2);
    // return;

    // Rayleigh RD,TD
    if( RD!=NULL || TD!=NULL ){
        cmat2x2_mul(RU1, RD2, tmp1);
        cmat2x2_one_sub(tmp1);
        cmat2x2_inv(tmp1, tmp1, stats);  if(*stats==INVERSE_FAILURE)  return;
        cmat2x2_mul(tmp1, TD1, tmp2);

        // TD
        if(TD!=NULL){
            cmat2x2_mul(TD2, tmp2, TD); // 相同的逆阵，节省计算量
        }

        // RD
        if(RD!=NULL){
            cmat2x2_mul(RD2, tmp2, tmp1);
            cmat2x2_mul(TU1, tmp1, tmp2);
            cmat2x2_add(RD1, tmp2, RD);
        }
    }

    // Rayleigh RU,TU
    if( RU!=NULL || TU!=NULL ){
        cmat2x2_mul(RD2, RU1, tmp1);
        cmat2x2_one_sub(tmp1);
        cmat2x2_inv(tmp1, tmp1, stats);  if(*stats==INVERSE_FAILURE)  return;
        cmat2x2_mul(tmp1, TU2, tmp2);

        // TU
        if(TU!=NULL){
            cmat2x2_mul(TU1, tmp2, TU);
        }

        // RU
        if(RU!=NULL){
            cmat2x2_mul(RU1, tmp2, tmp1);
            cmat2x2_mul(TD2, tmp1, tmp2);
            cmat2x2_add(RU2, tmp2, RU);
        }
    }


    // Love RDL,TDL
    if(RDL!=NULL || TDL!=NULL){
        inv1T = inv0 * TDL1;
        // TDL
        if(TDL!=NULL){
            *TDL = TDL2 * inv1T;
        }
        // RDL
        if(RDL!=NULL){
            *RDL = RDL1 + TUL1*RDL2*inv1T;
        }
    }

    // Love RUL,TUL
    if(RUL!=NULL || TUL!=NULL){
        inv1T = inv0 * TUL2;
        // TUL
        if(TUL!=NULL){
            *TUL = TUL1 * inv1T;
        }

        // RUL
        if(RUL!=NULL){
            *RUL = RUL2 + TDL2*RUL1 *inv1T; 
        }
    }

}


void recursion_RT_2x2_imaginary(
    MYCOMPLEX xa1, MYCOMPLEX xb1, MYREAL thk, MYREAL k, // 使用上层的厚度
    MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL)
{
    MYCOMPLEX exa, exb, exab, ex2a, ex2b; 
    exa = exp(-k*thk*xa1);
    exb = exp(-k*thk*xb1);

    exab = exa * exb;
    ex2a = exa * exa;
    ex2b = exb * exb;


    // 虚拟层位不是介质物理间断面
    RU[0][0] *= ex2a;    RU[0][1] *= exab;  
    RU[1][0] *= exab;    RU[1][1] *= ex2b;  
    
    TD[0][0] *= exa;     TD[0][1] *= exa; 
    TD[1][0] *= exb;     TD[1][1] *= exb;

    TU[0][0] *= exa;     TU[0][1] *= exb; 
    TU[1][0] *= exa;     TU[1][1] *= exb;

    *RUL *= ex2b;
    *TDL *= exb;
    *TUL *= exb;
}




void get_qwv(
    bool ircvup, 
    const MYCOMPLEX R1[2][2], MYCOMPLEX RL1, 
    const MYCOMPLEX R2[2][2], MYCOMPLEX RL2, 
    const MYCOMPLEX coef[QWV_NUM][2], MYCOMPLEX qwv[QWV_NUM])
{
    MYCOMPLEX qw0[2], qw1[2], v0;
    MYCOMPLEX coefD[2] = {coef[0][0], coef[1][0]};
    MYCOMPLEX coefU[2] = {coef[0][1], coef[1][1]};
    if(ircvup){
        cmat2x1_mul(R2, coefD, qw0);
        qw0[0] += coefU[0]; qw0[1] += coefU[1]; 
        v0 = RL1 * (RL2*coef[2][0] + coef[2][1]);
    } else {
        cmat2x1_mul(R2, coefU, qw0);
        qw0[0] += coefD[0]; qw0[1] += coefD[1]; 
        v0 = RL1 * (coef[2][0] + RL2*coef[2][1]);
    }
    cmat2x1_mul(R1, qw0, qw1);

    qwv[0] = qw1[0];
    qwv[1] = qw1[1];
    qwv[2] = v0;
}

