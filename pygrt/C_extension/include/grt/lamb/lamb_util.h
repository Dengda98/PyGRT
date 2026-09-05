/**
 * @file   lamb_util.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    一些使用广义闭合解求解 Lamb 问题过程中可能用到的辅助函数
 */

#pragma once

#include <stdio.h>

#include "grt/common/const.h"

#define LAMB_POLY_SIZE 17
#define LAMB_POLY_EPS 1e-28

/* 仅用于数值稳定性警告，不改变 Lamb 问题的数学定义域 */
#define LAMB_NU_WARNING_MARGIN 1e-3
#define LAMB_SURFACE_DEPTH_WARNING_RATIO 1e-3

/** Lamb 模块共用的 x 多项式 */
typedef struct {
    cplx_t c[LAMB_POLY_SIZE];
    int degree;
} LAMB_POLY;

LAMB_POLY grt_lamb_poly_const(const cplx_t value);

LAMB_POLY grt_lamb_poly_x(void);

void grt_lamb_poly_trim(LAMB_POLY *poly);

LAMB_POLY grt_lamb_poly_mul(const LAMB_POLY a, const LAMB_POLY b);

LAMB_POLY grt_lamb_poly_factor(const cplx_t root, const bool plus);

cplx_t grt_lamb_poly_eval(const LAMB_POLY *poly, const cplx_t x);

void grt_lamb_poly_divide(
    const LAMB_POLY numerator, const LAMB_POLY denominator,
    LAMB_POLY *quotient, LAMB_POLY *remainder);

/** 判断复数的虚部是否可以视为零 */
bool grt_lamb_is_real(const cplx_t value);

/** 对接近零的实数取平方根 */
real_t grt_lamb_positive_sqrt(const real_t value, const char *name);

/** 将椭圆积分参数限制在有效范围内 */
real_t grt_lamb_clamp_elliptic_parameter(const real_t value, const real_t tolerance, const char *name);

/** 求解一元三次方程的三个复根 */
void grt_lamb_cubic_roots(const real_t a, const real_t b, const real_t c, cplx_t roots[3]);

/** 计算以实数为自变量的时间多项式 */
cplx_t grt_lamb_eval_time_coeff(const cplx_t *coefficient, const int degree, const real_t t);

/** 用三个点的二次插值计算一阶导数 */
real_t grt_lamb_derivative_three_points(
    const real_t x0, const real_t x1, const real_t x2,
    const real_t f0, const real_t f1, const real_t f2, const real_t x);

/** 对 Lamb 模块的时间积分项求时间导数 */
void grt_lamb_differentiate_Fk(
    const real_t *ts, const int nt, const real_t (*Fk)[3][3][3],
    real_t (*dG)[3][3][3]);

/** 解析 Lamb 模块的空间导数输出路径 */
void grt_lamb_parse_derivative_paths(const char *argument, char **source_path, char **receiver_path);

/** 输出 Lamb 模块的 Green 函数序列 */
void grt_lamb_print_green_series(FILE *fp, const real_t *ts, const int nt, const real_t (*G)[3][3]);

/** 输出 Lamb 模块的一类空间导数序列 */
void grt_lamb_print_derivative_series(FILE *fp, const real_t *ts, const int nt, const real_t (*dG)[3][3][3], const bool source);

/**
 * 求解如下一元三次形式的 Rayleigh 方程的根,  其中 \f$ \nu \f$ 为泊松比
 * \f[
 *       x^3 - \dfrac{2\nu^2 + 1}{2(1 - \nu)} x^2
 *     + \dfrac{4\nu^3 - 4\nu^2 + 4\nu - 1}{4(1 - \nu)^2} x 
 *     - \dfrac{\nu^4}{8(1-\nu)^3} = 0
 * \f]
 * 
 * 
 * @param[in]      nu    泊松比， (0, 0.5)
 * @param[out]     y3    三个根，其中 y3[2] 为正根
 */
void grt_rayleigh1_roots(real_t nu, cplx_t y3[3]);


/**
 * 求解如下一元三次形式的 Rayleigh 方程的根,  其中 \f$ m=\dfrac{1}{2}\dfrac{1-2\nu}{1-\nu}, \nu \f$ 为泊松比
 * \f[
 *       x^3 + \dfrac{2m - 3}{2(1 - m)} x^2
 *     + \dfrac{1}{2(1-m)} x
 *     - \dfrac{1}{16(1-m)} = 0
 * \f]
 * 
 * 
 * @param[in]      m     系数 m
 * @param[out]     y3    三个根，其中 y3[2] 为正根
 */
void grt_rayleigh2_roots(real_t m, cplx_t y3[3]);

/**
 * 做如下多项式求值， \f$ \sum_{m=0}^n C_{2m+o} y^m \f$
 * 
 * @param[in]    C       数组 C
 * @param[in]    n       最高幂次 n
 * @param[in]    y       自变量 y
 * @param[in]    o       偏移量
 * 
 * @return    多项式结果
 * 
 */
cplx_t grt_evalpoly2(const cplx_t *C, const int n, const cplx_t y, const int offset);
