/**
 * @file   attenuation.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 
 */


#pragma once 

#include "common/const.h"

/**
 *  品质因子Q 对 波速的影响, Futterman causal Q, 参考Aki&Richards, 2009, Chapter 5.5.1
 * \f[
 *  c(\omega) = c(2\pi)\times (1 + \frac{1}{\pi Q} \log(\frac{\omega}{2\pi}) + \frac{i}{2Q}) 
 * \f] 
 * 其中虚数部分的正负号和书中不同，是因为书中使用的傅里叶变换的e指数符号和我们通常使用的相反。
 * 
 * @param[in]    Qinv     1/Q
 * @param[in]    omega    复数频率\f$ \tilde{\omega} =\omega - i\zeta \f$ 
 * 
 * @return atncoef 系数因子，作用在 \f$ k=\omega / c(\omega)\f$的计算
 */
MYCOMPLEX attenuation_law(MYREAL Qinv, MYCOMPLEX omega);

/**
 * attenuation_law函数在python中被调用的版本，长度2的数组分别表示复数的实部和虚部
 */
void py_attenuation_law(MYREAL Qinv, MYREAL omg[2], MYREAL atte[2]);