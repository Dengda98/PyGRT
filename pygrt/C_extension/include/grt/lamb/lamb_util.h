/**
 * @file   lamb_util.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    一些使用广义闭合解求解 Lamb 问题过程中可能用到的辅助函数
 */

#pragma once

#include "grt/common/const.h"

/**
 * 求解一元三次形式的 Rayleigh 方程的根,  \f$ x^3 + ax^2 + bx + c = 0 \f$
 * 
 * @param[in]      a     系数 a
 * @param[in]      b     系数 b
 * @param[in]      c     系数 c
 * @param[out]     y3    三个根，其中 y3[2] 为正根
 */
void grt_roots3(const MYREAL a, const MYREAL b, const MYREAL c, MYCOMPLEX y3[3]);

/**
 * 做如下多项式求值， \f$ \sum_{m=0}^n C_{2m+o} y^m \f$
 * 
 * @param[in]    C       数组 C
 * @param[in]    n       最高幂次 n
 * @param[in]    y       自变量 y
 * @param[in]    o       偏移量
 * 
 * @return    多项式结果
 * 
 */
MYCOMPLEX grt_evalpoly2(const MYCOMPLEX *C, const int n, const MYCOMPLEX y, const int offset);