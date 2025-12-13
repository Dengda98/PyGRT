/**
 * @file   static_grn.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 * 以下代码实现的是 广义反射透射系数矩阵+离散波数法 计算静态格林函数，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *              透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 * 
 */



#pragma once


#include "grt/common/model.h"
#include "grt/integral/integ_method.h"


/**
 * 积分计算Z, R, T三个分量静态格林函数的核心函数
 * 
 * @param[in,out]      mod1d            `GRT_MODEL1D` 结构体指针 
 * @param[in]      nr               震中距数量
 * @param[in]      rs               震中距数组 
 * @param[in,out]   Kmet            波数积分相关参数的结构体指针
 * @param[in]       calc_upar         是否计算位移u的空间导数
 * @param[out]      grn               浮点数数组，不同震源不同阶数的静态格林函数的Z、R、T分量
 * @param[out]      grn_uiz           浮点数数组，不同震源不同阶数的ui_z的Z、R、T分量
 * @param[out]      grn_uir           浮点数数组，不同震源不同阶数的ui_r的Z、R、T分量
 * 
 * @param[in]       statsstr           积分过程输出目录
 * 
 */
void grt_integ_static_grn(
    GRT_MODEL1D *mod1d, size_t nr, real_t *rs, K_INTEG_METHOD *Kmet,
    bool calc_upar, 
    realChnlGrid grn[nr],
    realChnlGrid grn_uiz[nr],
    realChnlGrid grn_uir[nr],
    const char *statsstr);