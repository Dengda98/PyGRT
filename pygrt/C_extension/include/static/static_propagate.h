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

void static_kernel(
    const MODEL1D *mod1d, MYREAL k,
    MYCOMPLEX EXP_qwv[3][3], MYCOMPLEX VF_qwv[3][3], MYCOMPLEX HF_qwv[3][3], MYCOMPLEX DC_qwv[3][3]);