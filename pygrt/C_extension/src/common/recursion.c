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

#include "grt/common/recursion.h"
#include "grt/common/const.h"
#include "grt/common/matrix.h"


void grt_recursion_RD(
    const MYCOMPLEX RD1[2][2], MYCOMPLEX RDL1, const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats)
{
    grt_recursion_RD_PSV(
        RD1, RU1, TD1, TU1,
        RD2,
        RD, inv_2x2T, stats);
    grt_recursion_RD_SH(
        RDL1, RUL1, TDL1, TUL1,
        RDL2,
        RDL, invT, stats);
}

void grt_recursion_RD_PSV(
    const MYCOMPLEX RD1[2][2], const MYCOMPLEX RU1[2][2],
    const MYCOMPLEX TD1[2][2], const MYCOMPLEX TU1[2][2],
    const MYCOMPLEX RD2[2][2],
    MYCOMPLEX RD[2][2], MYCOMPLEX inv_2x2T[2][2], MYINT *stats)
{
    MYCOMPLEX tmp1[2][2], tmp2[2][2];

    // RD, RDL
    grt_cmat2x2_mul(RU1, RD2, tmp1);
    grt_cmat2x2_one_sub(tmp1);
    grt_cmat2x2_inv(tmp1, tmp1, stats);  if(*stats==GRT_INVERSE_FAILURE)  return;
    grt_cmat2x2_mul(tmp1, TD1, tmp2);
    if(inv_2x2T!=NULL) grt_cmat2x2_assign(tmp2, inv_2x2T);

    grt_cmat2x2_mul(RD2, tmp2, tmp1);
    grt_cmat2x2_mul(TU1, tmp1, tmp2);
    grt_cmat2x2_add(RD1, tmp2, RD);
}

void grt_recursion_RD_SH(
    MYCOMPLEX RDL1, MYCOMPLEX RUL1,
    MYCOMPLEX TDL1, MYCOMPLEX TUL1,
    MYCOMPLEX RDL2, MYCOMPLEX *RDL, MYCOMPLEX *invT, MYINT *stats)
{
    MYCOMPLEX inv1;

    inv1 = 1.0 - RUL1*RDL2;
    if(inv1 == 0.0){
        *stats=GRT_INVERSE_FAILURE;
        return;
    }
    inv1 = 1.0 / inv1 * TDL1;
    *RDL = RDL1 + TUL1*RDL2*inv1;
    if(invT!=NULL)  *invT = inv1;
}


void grt_recursion_TD(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, 
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, 
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats)
{
    grt_recursion_TD_PSV(
        RU1, TD1, RD2, TD2,
        TD, inv_2x2T, stats);
    grt_recursion_TD_SH(
        RUL1, TDL1, RDL2, TDL2,
        TDL, invT, stats);
}

void grt_recursion_TD_PSV(
    const MYCOMPLEX RU1[2][2], const MYCOMPLEX TD1[2][2],
    const MYCOMPLEX RD2[2][2], const MYCOMPLEX TD2[2][2],
    MYCOMPLEX TD[2][2], MYCOMPLEX inv_2x2T[2][2], MYINT *stats)
{
    MYCOMPLEX tmp1[2][2], tmp2[2][2];

    // TD, TDL
    grt_cmat2x2_mul(RU1, RD2, tmp2);
    grt_cmat2x2_one_sub(tmp2);
    grt_cmat2x2_inv(tmp2, tmp1, stats);  if(*stats==GRT_INVERSE_FAILURE)  return;
    grt_cmat2x2_mul(tmp1, TD1, tmp2);
    if(inv_2x2T!=NULL)  grt_cmat2x2_assign(tmp2, inv_2x2T);
    grt_cmat2x2_mul(TD2, tmp2, TD);
}

void grt_recursion_TD_SH(
    MYCOMPLEX RUL1, MYCOMPLEX TDL1, 
    MYCOMPLEX RDL2, MYCOMPLEX TDL2, 
    MYCOMPLEX *TDL, MYCOMPLEX *invT, MYINT *stats)
{
    MYCOMPLEX inv1;

    inv1 = 1.0 - RUL1*RDL2;
    if(inv1 == 0.0){
        *stats=GRT_INVERSE_FAILURE;
        return;
    }
    inv1 = 1.0 / inv1 * TDL1;
    *TDL = TDL2 * inv1;
    if(invT!=NULL) *invT = inv1;
}

void grt_recursion_RU(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, const MYCOMPLEX RU2[2][2], MYCOMPLEX RUL2,
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats)
{
    grt_recursion_RU_PSV(
        RU1, RD2, RU2, TD2,
        TU2,
        RU, inv_2x2T, stats);
    grt_recursion_RU_SH(
        RUL1, RDL2, RUL2, TDL2,
        TUL2,
        RUL, invT, stats);
}

void grt_recursion_RU_PSV(
    const MYCOMPLEX RU1[2][2], const MYCOMPLEX RD2[2][2], 
    const MYCOMPLEX RU2[2][2], const MYCOMPLEX TD2[2][2], 
    const MYCOMPLEX TU2[2][2], 
    MYCOMPLEX RU[2][2], MYCOMPLEX inv_2x2T[2][2], MYINT *stats)
{
    MYCOMPLEX tmp1[2][2], tmp2[2][2];

    // RU, RUL
    grt_cmat2x2_mul(RD2, RU1, tmp2);
    grt_cmat2x2_one_sub(tmp2);
    grt_cmat2x2_inv(tmp2, tmp1, stats);  if(*stats==GRT_INVERSE_FAILURE)  return;
    grt_cmat2x2_mul(tmp1, TU2, tmp2);
    if(inv_2x2T!=NULL)  grt_cmat2x2_assign(tmp2, inv_2x2T);

    grt_cmat2x2_mul(RU1, tmp2, tmp1); 
    grt_cmat2x2_mul(TD2, tmp1, tmp2);
    grt_cmat2x2_add(RU2, tmp2, RU);
}


void grt_recursion_RU_SH(
    MYCOMPLEX RUL1, MYCOMPLEX RDL2,
    MYCOMPLEX RUL2, MYCOMPLEX TDL2, 
    MYCOMPLEX TUL2,
    MYCOMPLEX *RUL, MYCOMPLEX *invT, MYINT *stats)
{
    MYCOMPLEX inv1;

    inv1 = 1.0 - RUL1*RDL2;
    if(inv1 == 0.0){
        *stats=GRT_INVERSE_FAILURE;
        return;
    }
    inv1 = 1.0 / inv1 * TUL2;
    *RUL = RUL2 + TDL2*RUL1*inv1; 
    if(invT!=NULL)  *invT = inv1;
}


void grt_recursion_TU(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2,
    const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats)
{
    grt_recursion_TU_PSV(
        RU1, TU1, RD2, TU2,
        TU, inv_2x2T, stats);
    grt_recursion_TU_SH(
        RUL1, TUL1, RDL2, TUL2,
        TUL, invT, stats);
}


void grt_recursion_TU_PSV(
    const MYCOMPLEX RU1[2][2], const MYCOMPLEX TU1[2][2],
    const MYCOMPLEX RD2[2][2], const MYCOMPLEX TU2[2][2],
    MYCOMPLEX TU[2][2], MYCOMPLEX inv_2x2T[2][2], MYINT *stats)
{
    MYCOMPLEX tmp1[2][2], tmp2[2][2];

    // TU, TUL
    grt_cmat2x2_mul(RD2, RU1, tmp2);
    grt_cmat2x2_one_sub(tmp2);
    grt_cmat2x2_inv(tmp2, tmp1, stats);  if(*stats==GRT_INVERSE_FAILURE)  return;
    grt_cmat2x2_mul(tmp1, TU2, tmp2);
    if(inv_2x2T!=NULL) grt_cmat2x2_assign(tmp2, inv_2x2T);
    grt_cmat2x2_mul(TU1, tmp2, TU);
}



void grt_recursion_TU_SH(
    MYCOMPLEX RUL1, MYCOMPLEX TUL1,
    MYCOMPLEX RDL2, MYCOMPLEX TUL2,
    MYCOMPLEX *TUL, MYCOMPLEX *invT, MYINT *stats)
{
    MYCOMPLEX inv1;

    inv1 = 1.0 - RUL1*RDL2;
    if(inv1 == 0.0){
        *stats=GRT_INVERSE_FAILURE;
        return;
    }
    inv1 = 1.0 / inv1 * TUL2;
    *TUL = TUL1 * inv1;
    if(invT!=NULL)  *invT = inv1;
}


void grt_recursion_RT_matrix(
    const MYCOMPLEX RD1[2][2], MYCOMPLEX RDL1, const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, const MYCOMPLEX RU2[2][2], MYCOMPLEX RUL2,
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL,
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats)
{
    grt_recursion_RT_matrix_PSV(
        RD1, RU1, TD1, TU1,
        RD2, RU2, TD2, TU2,
        RD, RU, TD, TU, stats);
    grt_recursion_RT_matrix_SH(
        RDL1, RUL1, TDL1, TUL1,
        RDL2, RUL2, TDL2, TUL2,
        RDL, RUL, TDL, TUL, stats);
}


void grt_recursion_RT_matrix_PSV(
    const MYCOMPLEX RD1[2][2], const MYCOMPLEX RU1[2][2],
    const MYCOMPLEX TD1[2][2], const MYCOMPLEX TU1[2][2],
    const MYCOMPLEX RD2[2][2], const MYCOMPLEX RU2[2][2],
    const MYCOMPLEX TD2[2][2], const MYCOMPLEX TU2[2][2],
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2],
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2], MYINT *stats)
{
    // 临时矩阵
    MYCOMPLEX tmp1[2][2], tmp2[2][2];

    grt_cmat2x2_mul(RU1, RD2, tmp1);
    grt_cmat2x2_one_sub(tmp1);
    grt_cmat2x2_inv(tmp1, tmp1, stats);  if(*stats==GRT_INVERSE_FAILURE)  return;
    grt_cmat2x2_mul(tmp1, TD1, tmp2);

    // TD
    grt_cmat2x2_mul(TD2, tmp2, TD); // 相同的逆阵，节省计算量

    // RD
    grt_cmat2x2_mul(RD2, tmp2, tmp1);
    grt_cmat2x2_mul(TU1, tmp1, tmp2);
    grt_cmat2x2_add(RD1, tmp2, RD);

    grt_cmat2x2_mul(RD2, RU1, tmp1);
    grt_cmat2x2_one_sub(tmp1);
    grt_cmat2x2_inv(tmp1, tmp1, stats);  if(*stats==GRT_INVERSE_FAILURE)  return;
    grt_cmat2x2_mul(tmp1, TU2, tmp2);

    // TU
    grt_cmat2x2_mul(TU1, tmp2, TU);

    // RU
    grt_cmat2x2_mul(RU1, tmp2, tmp1);
    grt_cmat2x2_mul(TD2, tmp1, tmp2);
    grt_cmat2x2_add(RU2, tmp2, RU);
}


void grt_recursion_RT_matrix_SH(
    MYCOMPLEX RDL1, MYCOMPLEX RUL1,
    MYCOMPLEX TDL1, MYCOMPLEX TUL1,
    MYCOMPLEX RDL2, MYCOMPLEX RUL2,
    MYCOMPLEX TDL2, MYCOMPLEX TUL2,
    MYCOMPLEX *RDL, MYCOMPLEX *RUL,
    MYCOMPLEX *TDL, MYCOMPLEX *TUL, MYINT *stats)
{
    // 临时
    MYCOMPLEX inv0, inv1T;

    inv0 = 1.0 - RUL1*RDL2;
    if(inv0 == 0.0){
        *stats=GRT_INVERSE_FAILURE;
        return;
    }
    inv0 = 1.0 / inv0;

    inv1T = inv0 * TDL1;
    // TDL
    *TDL = TDL2 * inv1T;
    // RDL
    *RDL = RDL1 + TUL1*RDL2*inv1T;

    inv1T = inv0 * TUL2;
    // TUL
    *TUL = TUL1 * inv1T;

    // RUL
    *RUL = RUL2 + TDL2*RUL1 *inv1T; 
}


void grt_recursion_RT_matrix_virtual(
    MYCOMPLEX xa1, MYCOMPLEX xb1, MYREAL thk, MYREAL k, // 使用上层的厚度
    MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL)
{
    grt_recursion_RT_matrix_PSV_virtual(
        xa1, xb1, thk, k,
        RU, TD, TU);
    grt_recursion_RT_matrix_SH_virtual(
        xb1, thk, k,
        RUL, TDL, TUL);
}


void grt_recursion_RT_matrix_PSV_virtual(
    MYCOMPLEX xa1, MYCOMPLEX xb1, MYREAL thk, MYREAL k, // 使用上层的厚度
    MYCOMPLEX RU[2][2], MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2])
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
}


void grt_recursion_RT_matrix_SH_virtual(
    MYCOMPLEX xb1, MYREAL thk, MYREAL k, // 使用上层的厚度
    MYCOMPLEX *RUL, MYCOMPLEX *TDL, MYCOMPLEX *TUL)
{
    MYCOMPLEX exb, ex2b; 
    exb = exp(-k*thk*xb1);
    ex2b = exb * exb;

    // 虚拟层位不是介质物理间断面
    *RUL *= ex2b;
    *TDL *= exb;
    *TUL *= exb;
}



