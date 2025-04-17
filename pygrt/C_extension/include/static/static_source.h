/**
 * @file   static_source.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 * 以下代码实现的是 静态震源系数————爆炸源，垂直力源，水平力源，剪切源， 参考：
 *             1. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *                透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 *
 */
#pragma once

#include "common/const.h"
 
/**
 * 计算不同震源的静态震源系数，文献/书中仅提供剪切源的震源系数，其它震源系数重新推导
 * 
 * 数组形状[3][3][2]，代表在[i][j][p]时表示m=i阶时的
 * P(j=0),SV(j=1),SH(j=2)的震源系数(分别可记为q,w,v)，且分为下行波(p=0)和上行波(p=1). 
 * 
 * @param  delta           (in)震源层的\f$ \Delta \f$
 * @param  k               (in)波数
 * @param  EXP[3][3][2]    (out)爆炸源的震源系数，下同
 * @param  VF[3][3][2]     (out)垂直力源
 * @param  HF[3][3][2]     (out)水平力源
 * @param  DC[3][3][2]     (out)剪切源
 */
void static_source_coef(
    MYCOMPLEX delta, MYREAL k,
    MYCOMPLEX EXP[3][3][2], MYCOMPLEX VF[3][3][2], MYCOMPLEX HF[3][3][2], MYCOMPLEX DC[3][3][2]);
 