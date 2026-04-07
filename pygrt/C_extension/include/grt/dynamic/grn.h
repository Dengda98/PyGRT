/**
 * @file   grn.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是利用多种波数积分类方法计算理论地震图，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient 
 *               matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699
 *   
 *                 
 */

#pragma once 

#include "grt/common/model.h"
#include "grt/dynamic/grnspec.h"
#include "grt/integral/integ_method.h"

/**
 * 积分计算Z, R, T三个分量格林函数的频谱的核心函数（被Python调用）  
 * 
 * @param[in,out]   mod1d            `MODEL1D` 结构体指针 
 * @param[in,out]   Kmet             波数积分相关参数的结构体指针
 * @param[in,out]   grn              计算的格林函数频谱
 * @param[in]       print_progressbar        是否打印进度条
 * 
 */ 
void grt_integ_grn_spec(MODEL1D *mod1d, K_INTEG_METHOD *Kmet, GRNSPEC *grn, const bool print_progressbar);




