/**
 * @file   elliptic.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    使用渐近解和 Carlson 算法计算三类完全椭圆积分，参考：
 *    
 *        1. Abramowitz, Milton, and Irene A. Stegun. 1964. Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables. Vol. 55. US Government printing office.
 *        2. Carlson, B. C. 1995. “Numerical Computation of Real or Complex Elliptic Integrals.” Numerical Algorithms 10 (1): 13–26. https://doi.org/10.1007/BF02198293.
 *        3. 张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社
 */

#pragma once

#include "grt/common/const.h"

/**
 * 使用 Abramowitz and Stegun (1964, P591-592) 提供的渐近表达式计算第一类完全椭圆积分， \f$ m \in (0,1) \f$
 * \f[
 *  K(m) = \int_0^{\pi/2} \dfrac{1}{\sqrt{1 - m \sin^2\theta}}\ d\theta
 * \f]  
 * 
 * @param[in]     m     参数 m
 * @return    积分值
 * 
 */
real_t grt_ellipticK(const real_t m);


/**
 * 使用 Abramowitz and Stegun (1964, P591-592) 提供的渐近表达式计算第二类完全椭圆积分， \f$ m \in (0,1) \f$
 * \f[
 *  E(m) = \int_0^{\pi/2} \sqrt{1 - m \sin^2\theta}\ d\theta
 * \f]  
 * 
 * @param[in]     m     参数 m
 * @return    积分值
 * 
 */
real_t grt_ellipticE(const real_t m);


/**
 * 使用 Carlson (1995) 算法计算第三类完全椭圆积分， \f$ m \in (0,1) \f$
 * \f[
 *  {\it \Pi}(m) = \int_0^{\pi/2} \dfrac{1}{(1-nx^2)\sqrt{1-x^2}\sqrt{1 - m \sin^2\theta}}\ d\theta
 * \f]  
 * 
 * @param[in]     n     参数 n
 * @param[in]     m     参数 m
 * @return    积分值
 * 
 */
cplx_t grt_ellipticPi(const cplx_t n, const real_t m);
