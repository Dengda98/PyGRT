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

#include "grt/common/const.h"
#include "grt/common/model.h"


/**
 * 根据公式(4.6.6)，(4.6.15)，(4.6.21,26)，(4.8.34)计算不同震源不同阶数的震源系数，
 * 数组形状代表在[i][j][p]时表示i类震源的
 * P(j=0),SV(j=1)的震源系数(分别对应q,w)，且分为下行波(p=0)和上行波(p=1). 
 * 
 * @param[in,out]     mod1d        模型结构体指针，结果保存在 src_coef
 * 
 */
void grt_source_coef(GRT_MODEL1D *mod1d);

/* P-SV 波的震源系数  */
void grt_source_coef_PSV(GRT_MODEL1D *mod1d);

/* SH 波的震源系数  */
void grt_source_coef_SH(GRT_MODEL1D *mod1d);
