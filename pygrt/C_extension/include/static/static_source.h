/**
 * @file   static_source.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 * 以下代码实现的是 静态震源系数————剪切源， 参考：
 *             1. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *                透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 *
 */
#pragma once

#include "common/const.h"
 
 
void static_source_coef(
    MYCOMPLEX delta, 
    MYCOMPLEX EXP[3][3][2], MYCOMPLEX VF[3][3][2], MYCOMPLEX HF[3][3][2], MYCOMPLEX DC[3][3][2]);
 