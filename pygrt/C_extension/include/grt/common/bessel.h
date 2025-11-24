/**
 * @file   bessel.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 *                   
 */

#pragma once

#include "grt/common/const.h"

/**
 * 计算Bessel函数 \f$ J_m(x), m=0,1,2 \f$ 
 * 
 * @param[in]   x          自变量 
 * @param[out]  bj0        \f$ J_0(x) \f$
 * @param[out]  bj1        \f$ J_1(x) \f$
 * @param[out]  bj2        \f$ J_2(x) \f$
 * 
 */
void grt_bessel012(real_t x, real_t *bj0, real_t *bj1, real_t *bj2);


/**
 * 计算Bessel函数的一阶导数 \f$ J_m^{'}(x), m=0,1,2 \f$ 
 * 
 * @param[in]       x          自变量 
 * @param[in,out]   bj0        传入 \f$ J_0(x) \f$, 返回\f$ J_0^{'}(x) \f$
 * @param[in,out]   bj1        传入 \f$ J_1(x) \f$, 返回\f$ J_1^{'}(x) \f$
 * @param[in,out]   bj2        传入 \f$ J_2(x) \f$, 返回\f$ J_2^{'}(x) \f$
 * 
 */
void grt_besselp012(real_t x, real_t *bj0, real_t *bj1, real_t *bj2);