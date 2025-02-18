/**
 * @file   source.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 震源系数————爆炸源，垂直力源，水平力源，剪切源， 参考：
 * 
 *          1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.    
 *                 
 */

#pragma once 

#include "const.h"
#include "model.h"


/**
 * 根据公式(4.6.6)，(4.6.15)，(4.6.21,26)，(4.8.34)计算不同震源的震源系数，
 * 数组形状[3][3][2]，代表在[i][j][p]时表示m=i阶时的
 * P(j=0),SV(j=1),SH(j=2)的震源系数(分别可记为q,w,v)，且分为下行波(p=0)和上行波(p=1). 
 * 
 * @param     src_xa      (in)震源层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param     src_xb      (in)震源层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param     src_kaka   (in)震源层的P波水平波数的平方 \f$ k_a^2=(\frac{\omega}{V_a})^2 \f$
 * @param     src_kbkb   (in)震源层的S波水平波数的平方 \f$ k_b^2=(\frac{\omega}{V_b})^2 \f$
 * @param  omega     (in)复数频率
 * @param  k         (in)波数
 * @param  EXP[3][3][2]    (out)爆炸源的震源系数，下同
 * @param  VF[3][3][2]     (out)垂直力源
 * @param  HF[3][3][2]     (out)水平力源
 * @param  DC[3][3][2]     (out)双力偶源
 * 
 */
void source_coef(
    MYCOMPLEX src_xa, MYCOMPLEX src_xb, MYCOMPLEX src_kaka, MYCOMPLEX src_kbkb, 
    MYCOMPLEX omega, MYREAL k,
    MYCOMPLEX EXP[3][3][2], MYCOMPLEX VF[3][3][2], MYCOMPLEX HF[3][3][2], MYCOMPLEX DC[3][3][2]);
