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
 * 此处omega未使用，传入0即可
 */
void static_kernel(
    const MODEL1D *mod1d, MYCOMPLEX omega, MYREAL k,
    MYCOMPLEX EXP_qwv[3][3], MYCOMPLEX VF_qwv[3][3], MYCOMPLEX HF_qwv[3][3], MYCOMPLEX DC_qwv[3][3],
    bool calc_uiz,
    MYCOMPLEX EXP_uiz_qwv[3][3], MYCOMPLEX VF_uiz_qwv[3][3], MYCOMPLEX HF_uiz_qwv[3][3], MYCOMPLEX DC_uiz_qwv[3][3]);