/**
 * @file   bessel.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 *                   
 */

#pragma once

#include "common/const.h"

/**
 * 计算Bessel函数 \f$ J_m(x), m=0,1,2 \f$ 
 * 
 * @param x          自变量 
 * @param bj0  (out) \f$ J_0(x) \f$
 * @param bj1  (out) \f$ J_1(x) \f$
 * @param bj2  (out) \f$ J_2(x) \f$
 * 
 */
void bessel012(MYREAL x, MYREAL *bj0, MYREAL *bj1, MYREAL *bj2);


/**
 * 计算Bessel函数的一阶导数 \f$ J_m^{'}(x), m=0,1,2 \f$ 
 * 
 * @param x          自变量 
 * @param bj0  (inout) 传入 \f$ J_0(x) \f$, 返回\f$ J_0^{'}(x) \f$
 * @param bj1  (inout) 传入 \f$ J_1(x) \f$, 返回\f$ J_1^{'}(x) \f$
 * @param bj2  (inout) 传入 \f$ J_2(x) \f$, 返回\f$ J_2^{'}(x) \f$
 * 
 */
void besselp012(MYREAL x, MYREAL *bj0, MYREAL *bj1, MYREAL *bj2);