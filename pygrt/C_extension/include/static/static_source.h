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
 * 数组形状代表在[i][j][p]时表示i类震源的
 * P(j=0),SV(j=1),SH(j=2)的震源系数(分别对应q,w,v)，且分为下行波(p=0)和上行波(p=1). 
 * 
 * @param[in]     delta    震源层的\f$ \Delta \f$
 * @param[in]     k        波数
 * @param[out]    coef     震源系数 \f$ P_m, SV_m, SH_m \f$
 */
void static_source_coef(
    MYCOMPLEX delta, MYREAL k, MYCOMPLEX coef[SRC_M_NUM][QWV_NUM][2]);
 