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

#include "grt/common/const.h"
#include "grt/common/model.h"
 
/**
 * 计算不同震源的静态震源系数，文献/书中仅提供剪切源的震源系数，其它震源系数重新推导
 * 
 * 数组形状代表在[i][j][p]时表示i类震源的
 * P(j=0),SV(j=1)的震源系数(分别对应q,w)，且分为下行波(p=0)和上行波(p=1). 
 * 
 * @param[in]     mod1d    模型结构体指针
 * @param[out]    coef     震源系数 \f$ P_m, SV_m, SH_m \f$
 */
void grt_static_source_coef_PSV(const GRT_MODEL1D *mod1d, cplx_t coef[GRT_SRC_M_NUM][GRT_QWV_NUM-1][2]);

/* SH 波的静态震源系数 */
void grt_static_source_coef_SH(const GRT_MODEL1D *mod1d, cplx_t coef[GRT_SRC_M_NUM][2]);