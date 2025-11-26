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

#include "grt/common/const.h"
#include "grt/common/matrix.h"
#include "grt/common/RT_matrix.h"


void grt_recursion_RD(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    grt_recursion_RD_PSV(M1, M2, M);
    grt_recursion_RD_SH(M1, M2, M);
}

void grt_recursion_RD_PSV(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    cplx_t tmp1[2][2], tmp2[2][2];

    // RD, RDL
    grt_cmat2x2_mul(M1->RU, M2->RD, tmp1);
    grt_cmat2x2_one_sub(tmp1);
    if((M->stats = grt_cmat2x2_inv(tmp1, tmp1)) == GRT_INVERSE_FAILURE) return;
    grt_cmat2x2_mul(tmp1, M1->TD, tmp2);
    grt_cmat2x2_assign(tmp2, M->invT);

    grt_cmat2x2_mul(M2->RD, tmp2, tmp1);
    grt_cmat2x2_mul(M1->TU, tmp1, tmp2);
    grt_cmat2x2_add(M1->RD, tmp2, M->RD);
}

void grt_recursion_RD_SH(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    cplx_t inv1;

    inv1 = 1.0 - M1->RUL * M2->RDL;
    if(inv1 == 0.0){
        M->stats=GRT_INVERSE_FAILURE;
        return;
    }
    inv1 = 1.0 / inv1 * M1->TDL;
    M->RDL = M1->RDL + M1->TUL * M2->RDL * inv1;
    M->invTL = inv1;
}


void grt_recursion_TD(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    grt_recursion_TD_PSV(M1, M2, M);
    grt_recursion_TD_SH(M1, M2, M);
}

void grt_recursion_TD_PSV(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    cplx_t tmp1[2][2], tmp2[2][2];

    // TD, TDL
    grt_cmat2x2_mul(M1->RU, M2->RD, tmp2);
    grt_cmat2x2_one_sub(tmp2);
    if((M->stats = grt_cmat2x2_inv(tmp2, tmp1)) == GRT_INVERSE_FAILURE) return;
    grt_cmat2x2_mul(tmp1, M1->TD, tmp2);
    grt_cmat2x2_assign(tmp2, M->invT);
    grt_cmat2x2_mul(M2->TD, tmp2, M->TD);
}

void grt_recursion_TD_SH(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    cplx_t inv1;

    inv1 = 1.0 - M1->RUL * M2->RDL;
    if(inv1 == 0.0){
        M->stats=GRT_INVERSE_FAILURE;
        return;
    }
    inv1 = 1.0 / inv1 * M1->TDL;
    M->TDL = M2->TDL * inv1;
    M->invTL = inv1;
}

void grt_recursion_RU(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    grt_recursion_RU_PSV(M1, M2, M);
    grt_recursion_RU_SH(M1, M2, M);
}

void grt_recursion_RU_PSV(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    cplx_t tmp1[2][2], tmp2[2][2];

    // RU, RUL
    grt_cmat2x2_mul(M2->RD, M1->RU, tmp2);
    grt_cmat2x2_one_sub(tmp2);
    if((M->stats = grt_cmat2x2_inv(tmp2, tmp1)) == GRT_INVERSE_FAILURE) return;
    grt_cmat2x2_mul(tmp1, M2->TU, tmp2);
    grt_cmat2x2_assign(tmp2, M->invT);

    grt_cmat2x2_mul(M1->RU, tmp2, tmp1); 
    grt_cmat2x2_mul(M2->TD, tmp1, tmp2);
    grt_cmat2x2_add(M2->RU, tmp2, M->RU);
}


void grt_recursion_RU_SH(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    cplx_t inv1;

    inv1 = 1.0 - M1->RUL * M2->RDL;
    if(inv1 == 0.0){
        M->stats=GRT_INVERSE_FAILURE;
        return;
    }
    inv1 = 1.0 / inv1 * M2->TUL;
    M->RUL = M2->RUL + M2->TDL * M1->RUL * inv1; 
    M->invTL = inv1;
}


void grt_recursion_TU(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    grt_recursion_TU_PSV(M1, M2, M);
    grt_recursion_TU_SH(M1, M2, M);
}


void grt_recursion_TU_PSV(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    cplx_t tmp1[2][2], tmp2[2][2];

    // TU, TUL
    grt_cmat2x2_mul(M2->RD, M1->RU, tmp2);
    grt_cmat2x2_one_sub(tmp2);
    if((M->stats = grt_cmat2x2_inv(tmp2, tmp1)) == GRT_INVERSE_FAILURE) return;
    grt_cmat2x2_mul(tmp1, M2->TU, tmp2);
    grt_cmat2x2_assign(tmp2, M->invT);
    grt_cmat2x2_mul(M1->TU, tmp2, M->TU);
}



void grt_recursion_TU_SH(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    cplx_t inv1;

    inv1 = 1.0 - M1->RUL * M2->RDL;
    if(inv1 == 0.0){
        M->stats=GRT_INVERSE_FAILURE;
        return;
    }
    inv1 = 1.0 / inv1 * M2->TUL;
    M->TUL = M1->TUL * inv1;
    M->invTL = inv1;
}


void grt_recursion_RT_matrix(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    grt_recursion_RT_matrix_PSV(M1, M2, M);
    grt_recursion_RT_matrix_SH(M1, M2, M);
}


void grt_recursion_RT_matrix_PSV(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    // 临时矩阵
    cplx_t tmp1[2][2], tmp2[2][2];

    grt_cmat2x2_mul(M1->RU, M2->RD, tmp1);
    grt_cmat2x2_one_sub(tmp1);
    if((M->stats = grt_cmat2x2_inv(tmp1, tmp1)) == GRT_INVERSE_FAILURE) return;
    grt_cmat2x2_mul(tmp1, M1->TD, tmp2);

    // TD
    grt_cmat2x2_mul(M2->TD, tmp2, M->TD); // 相同的逆阵，节省计算量

    // RD
    grt_cmat2x2_mul(M2->RD, tmp2, tmp1);
    grt_cmat2x2_mul(M1->TU, tmp1, tmp2);
    grt_cmat2x2_add(M1->RD, tmp2, M->RD);

    grt_cmat2x2_mul(M2->RD, M1->RU, tmp1);
    grt_cmat2x2_one_sub(tmp1);
    if((M->stats = grt_cmat2x2_inv(tmp1, tmp1)) == GRT_INVERSE_FAILURE) return;
    grt_cmat2x2_mul(tmp1, M2->TU, tmp2);

    // TU
    grt_cmat2x2_mul(M1->TU, tmp2, M->TU);

    // RU
    grt_cmat2x2_mul(M1->RU, tmp2, tmp1);
    grt_cmat2x2_mul(M2->TD, tmp1, tmp2);
    grt_cmat2x2_add(M2->RU, tmp2, M->RU);
}


void grt_recursion_RT_matrix_SH(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M)
{
    // 临时
    cplx_t inv0, inv1T;

    inv0 = 1.0 - M1->RUL * M2->RDL;
    if(inv0 == 0.0){
        M->stats=GRT_INVERSE_FAILURE;
        return;
    }
    inv0 = 1.0 / inv0;

    inv1T = inv0 * M1->TDL;
    // TDL
    M->TDL = M2->TDL * inv1T;
    // RDL
    M->RDL = M1->RDL + M1->TUL * M2->RDL * inv1T;

    inv1T = inv0 * M2->TUL;
    // TUL
    M->TUL = M1->TUL * inv1T;

    // RUL
    M->RUL = M2->RUL + M2->TDL * M1->RUL * inv1T; 
}



// ============================================================================
//                               Debug
// ============================================================================

void grt_print_RT_matrix(const RT_MATRIX *M)
{
    fprintf(stderr, "RD\n");
    grt_cmatmxn_print(2, 2, M->RD);
    fprintf(stderr, "RDL="GRT_CMPLX_FMT"\n", GRT_CMPLX_SPLIT(M->RDL));
    fprintf(stderr, "RU\n");
    grt_cmatmxn_print(2, 2, M->RU);
    fprintf(stderr, "RUL="GRT_CMPLX_FMT"\n", GRT_CMPLX_SPLIT(M->RUL));
    fprintf(stderr, "TD\n");
    grt_cmatmxn_print(2, 2, M->TD);
    fprintf(stderr, "TDL="GRT_CMPLX_FMT"\n", GRT_CMPLX_SPLIT(M->TDL));
    fprintf(stderr, "TU\n");
    grt_cmatmxn_print(2, 2, M->TU);
    fprintf(stderr, "TUL="GRT_CMPLX_FMT"\n", GRT_CMPLX_SPLIT(M->TUL));
}