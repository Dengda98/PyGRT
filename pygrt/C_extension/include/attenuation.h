/**
 * @file   attenuation.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 
 */


#pragma once 

#include "const.h"

/**
 *  品质因子Q 对 波速的影响, Futterman causal Q, 参考Aki&Richards, 2009, Chapter 5.5.1
 * \f[
 *  c(\omega) = c(2\pi)\times (1 + \frac{1}{\pi Q} \log(\frac{\omega}{2\pi}) + \frac{i}{2Q}) 
 * \f] 
 * 
 * @param Qinv     (in) 1/Q
 * @param omega    (in)复数频率\f$ \tilde{\omega} =\omega - i\zeta \f$ 
 * @return atncoef 系数因子，作用在 \f$ k=\omega / c(\omega)\f$的计算
 */
MYCOMPLEX attenuation_law(MYREAL Qinv, MYCOMPLEX omega);