/**
 * @file   static_propagate.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 * 以下代码实现的是 静态广义反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *              透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 *
 */

#pragma once 

#include "common/const.h"
#include "common/model.h"


/**
 * 静态kernel函数根据(5.5.3)式递推计算静态广义反射透射矩阵。递推公式适用于动态和静态情况。
 * 函数参数与动态kernel函数保持一致，具体说明详见`dynamic/propagate.h`。
 * 
 */
void static_kernel(
    const MODEL1D *mod1d, MYCOMPLEX omega, MYREAL k, MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM],
    bool calc_uiz, MYCOMPLEX QWV_uiz[SRC_M_NUM][QWV_NUM], MYINT *stats);