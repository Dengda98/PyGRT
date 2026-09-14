/**
 * @file   lamb_poly.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-09
 *
 *    Lamb 问题闭合解中使用的两类多项式
 */

#pragma once

#include "grt/common/const.h"


#define LAMB_POLY_X_SIZE 24
#define LAMB_POLY_TBAR_SIZE 8
#define LAMB_X_POLY_EPS 1e-28


/** 仅含变量 x 的一元复系数多项式 */
typedef struct {
    cplx_t c[LAMB_POLY_X_SIZE]; ///< c[m] 是 x^m 的系数
    int degree;                 ///< x 的最高次数
} LAMB_X_POLY;


/** 按 tbar^r*x^m 保存系数的二维实系数多项式 */
typedef struct {
    real_t c[LAMB_POLY_TBAR_SIZE][LAMB_POLY_X_SIZE]; ///< c[r][m] 是 tbar^r*x^m 的系数
    int degree;                                      ///< x 的最高次数
    int tbar_degree;                                 ///< tbar 的最高次数
} LAMB_TBAR_X_POLY;


/** 构造一个常数 x 多项式 */
LAMB_X_POLY grt_lamb_x_poly_const(const cplx_t value);

/** 构造积分变量 x */
LAMB_X_POLY grt_lamb_x_poly_x(void);

/** 根据非零系数重新计算 x 多项式次数 */
void grt_lamb_x_poly_trim(LAMB_X_POLY *poly);

/** 相乘两个 x 多项式 */
LAMB_X_POLY grt_lamb_x_poly_mul(const LAMB_X_POLY a, const LAMB_X_POLY b);

/** 构造一个二次因子 */
LAMB_X_POLY grt_lamb_x_poly_factor(const cplx_t root, const bool plus);

/** 计算 x 多项式的值 */
cplx_t grt_lamb_x_poly_eval(const LAMB_X_POLY *poly, const cplx_t x);

/** 对两个 x 多项式做带余除法 */
void grt_lamb_x_poly_divide(
    const LAMB_X_POLY numerator, const LAMB_X_POLY denominator,
    LAMB_X_POLY *quotient, LAMB_X_POLY *remainder);


/** 清空一个 tbar-x 多项式 */
void grt_lamb_tbar_x_poly_zero(LAMB_TBAR_X_POLY *poly);

/** 构造一个常数 tbar-x 多项式 */
void grt_lamb_tbar_x_poly_const(LAMB_TBAR_X_POLY *poly, const real_t value);

/** 构造积分变量 x */
void grt_lamb_tbar_x_poly_x(LAMB_TBAR_X_POLY *poly);

/** 构造归一化时间变量 tbar */
void grt_lamb_tbar_x_poly_tbar(LAMB_TBAR_X_POLY *poly);

/** 根据非零系数重新计算 tbar-x 多项式次数 */
void grt_lamb_tbar_x_poly_trim(LAMB_TBAR_X_POLY *poly);

/** 相加两个 tbar-x 多项式 */
void grt_lamb_tbar_x_poly_add(const LAMB_TBAR_X_POLY *left, const LAMB_TBAR_X_POLY *right,
                              LAMB_TBAR_X_POLY *result);

/** 将缩放后的 tbar-x 多项式加到另一个多项式 */
void grt_lamb_tbar_x_poly_add_scaled(const LAMB_TBAR_X_POLY *base, const real_t scale,
                                     const LAMB_TBAR_X_POLY *term, LAMB_TBAR_X_POLY *result);

/** 缩放一个 tbar-x 多项式 */
void grt_lamb_tbar_x_poly_scale(const LAMB_TBAR_X_POLY *source, const real_t scale,
                                LAMB_TBAR_X_POLY *result);

/** 通过二维卷积相乘两个 tbar-x 多项式 */
void grt_lamb_tbar_x_poly_mul(const LAMB_TBAR_X_POLY *left, const LAMB_TBAR_X_POLY *right,
                              LAMB_TBAR_X_POLY *result);

/** 计算两个 tbar-x 多项式的乘积并缩放 */
void grt_lamb_tbar_x_poly_product2_scaled(const LAMB_TBAR_X_POLY *first, const LAMB_TBAR_X_POLY *second,
                                          const real_t scale, LAMB_TBAR_X_POLY *result);

/** 计算三个 tbar-x 多项式的乘积 */
void grt_lamb_tbar_x_poly_product3(const LAMB_TBAR_X_POLY *first, const LAMB_TBAR_X_POLY *second,
                                   const LAMB_TBAR_X_POLY *third, LAMB_TBAR_X_POLY *result);

/** 计算三个 tbar-x 多项式的乘积并缩放 */
void grt_lamb_tbar_x_poly_product3_scaled(const LAMB_TBAR_X_POLY *first, const LAMB_TBAR_X_POLY *second,
                                          const LAMB_TBAR_X_POLY *third, const real_t scale,
                                          LAMB_TBAR_X_POLY *result);

/** 计算四个 tbar-x 多项式的乘积并缩放 */
void grt_lamb_tbar_x_poly_product4_scaled(const LAMB_TBAR_X_POLY *first, const LAMB_TBAR_X_POLY *second,
                                          const LAMB_TBAR_X_POLY *third, const LAMB_TBAR_X_POLY *fourth,
                                          const real_t scale, LAMB_TBAR_X_POLY *result);

/** 计算五个 tbar-x 多项式的乘积并缩放 */
void grt_lamb_tbar_x_poly_product5_scaled(const LAMB_TBAR_X_POLY *first, const LAMB_TBAR_X_POLY *second,
                                          const LAMB_TBAR_X_POLY *third, const LAMB_TBAR_X_POLY *fourth,
                                          const LAMB_TBAR_X_POLY *fifth, const real_t scale,
                                          LAMB_TBAR_X_POLY *result);

/** 计算六个 tbar-x 多项式的乘积并缩放 */
void grt_lamb_tbar_x_poly_product6_scaled(const LAMB_TBAR_X_POLY *first, const LAMB_TBAR_X_POLY *second,
                                          const LAMB_TBAR_X_POLY *third, const LAMB_TBAR_X_POLY *fourth,
                                          const LAMB_TBAR_X_POLY *fifth, const LAMB_TBAR_X_POLY *sixth,
                                          const real_t scale, LAMB_TBAR_X_POLY *result);
