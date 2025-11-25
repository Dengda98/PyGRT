/**
 * @file   lamb_util.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    一些使用广义闭合解求解 Lamb 问题过程中可能用到的辅助函数
 */

#pragma once

#include "grt/common/const.h"

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