/**
 * @file   RTmatrix.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 *     R/T 系数矩阵，以及通过递推公式计算两层的广义 R/T 系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *
 */

#pragma once

#include "grt/common/const.h"

/* 使用结构体管理上行和下行的 R/T 矩阵 */
typedef struct {
    cplx_t RD[2][2];     ///< P-SV 下传反射系数矩阵
    cplx_t RU[2][2];     ///< P-SV 上传反射系数矩阵
    cplx_t TD[2][2];     ///< P-SV 下传透射系数矩阵
    cplx_t TU[2][2];     ///< P-SV 上传透射系数矩阵
    cplx_t RDL;          ///< SH 下传反射系数
    cplx_t RUL;          ///< SH 上传反射系数
    cplx_t TDL;          ///< SH 下传透射系数
    cplx_t TUL;          ///< SH 上传透射系数

    /* 一些辅助变量 */
    cplx_t invT[2][2];
    cplx_t invTL;

    /* 状态变量 */
    int stats;   ///< 是否有除零错误，非0为异常值
} RT_MATRIX;


/** 初始化 R/T 矩阵 */
#define grt_init_RT_matrix(M) \
    RT_MATRIX *M = &(RT_MATRIX){\
        .RD = GRT_INIT_ZERO_2x2_MATRIX,\
        .RU = GRT_INIT_ZERO_2x2_MATRIX,\
        .TD = GRT_INIT_IDENTITY_2x2_MATRIX,\
        .TU = GRT_INIT_IDENTITY_2x2_MATRIX,\
        .RDL = 0.0,\
        .RUL = 0.0,\
        .TDL = 1.0,\
        .TUL = 1.0,\
        .invT = GRT_INIT_ZERO_2x2_MATRIX,\
        .invTL = 0.0,\
        .stats = GRT_INVERSE_SUCCESS \
    };\


/** 
 * 合并 recursion_RD_PSV(SH) ，仅计算RD/RDL
 * 
 * @param[in]      M1            上层 R/T 系数
 * @param[in]      M2            下层 R/T 系数
 * @param[out]     M2            合并 R/T 系数
 * 
 */
void grt_recursion_RD(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);

/** 根据公式(5.5.3(1))进行递推 P-SV ，仅计算RD */
void grt_recursion_RD_PSV(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);

/** 根据公式(5.5.3(1))进行递推 SH ，仅计算RDL */
void grt_recursion_RD_SH(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);


/** 合并 recursion_TD_PSV(SH)，仅计算TD/TDL */
void grt_recursion_TD(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);

/** 根据公式(5.5.3(2))进行递推 P-SV ，仅计算TD */
void grt_recursion_TD_PSV(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);

/** 根据公式(5.5.3(2))进行递推 SH ，仅计算TDL */
void grt_recursion_TD_SH(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);


/** 合并 recursion_RU_PSV(SH)，仅计算RU/RUL */
void grt_recursion_RU(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);

/** 根据公式(5.5.3(3))进行递推 P-SV ，仅计算RU */
void grt_recursion_RU_PSV(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);

/** 根据公式(5.5.3(3))进行递推 SH ，仅计算RUL */
void grt_recursion_RU_SH(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);


/** 合并 recursion_TU_PSV(SH)，仅计算TU/TUL  */
void grt_recursion_TU(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);

/** 根据公式(5.5.3(4))进行递推 P-SV ，仅计算TU */
void grt_recursion_TU_PSV(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);

/** 根据公式(5.5.3(4))进行递推 SH ，仅计算TUL */
void grt_recursion_TU_SH(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);


/** 合并 recursion_RT_matrix_PSV(SH) */
void grt_recursion_RT_matrix(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);

/**
 * 根据公式(5.5.3)进行递推 P-SV ，相当于对应四个函数合并，
 * 内部使用了共有变量防止重复计算
 */
void grt_recursion_RT_matrix_PSV(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);

/**
 * 根据公式(5.5.3)进行递推 SH ，相当于对应四个函数合并，
 * 内部使用了共有变量防止重复计算
 */
void grt_recursion_RT_matrix_SH(const RT_MATRIX *M1, const RT_MATRIX *M2, RT_MATRIX *M);


/** 打印 R/T 矩阵 */
void grt_print_RT_matrix(const RT_MATRIX *M);