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
 * 数组形状代表在[i][j]时表示i类震源的P(j=0),SV(j=1),SH(j=2)的震源系数(分别对应q,w,v). 
 * 
 * @param[in]     mod1d        模型结构体指针
 * @param[out]    src_coefD    下行震源系数
 * @param[out]    src_coefU    上行震源系数
 * 
 */
void grt_source_coef(const GRT_MODEL1D *mod1d, cplxChnlGrid src_coefD, cplxChnlGrid src_coefU);

/* P-SV 波的震源系数  */
void grt_source_coef_PSV(const GRT_MODEL1D *mod1d, cplxChnlGrid src_coefD, cplxChnlGrid src_coefU);

/* SH 波的震源系数  */
void grt_source_coef_SH(const GRT_MODEL1D *mod1d, cplxChnlGrid src_coefD, cplxChnlGrid src_coefU);
