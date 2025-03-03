/**
 * @file   static_layer.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 * 以下代码实现的是 静态反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *              透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 *
 */

#pragma once

#include <stdbool.h>

#include "const.h"


void calc_static_R_tilt(MYCOMPLEX delta1, MYCOMPLEX R_tilt[2][2]);


void calc_static_R_EV(
    bool ircvup,
    const MYCOMPLEX R[2][2], MYCOMPLEX RL, 
    MYCOMPLEX R_EV[2][2], MYCOMPLEX *R_EVL);


void calc_static_RT_2x2(
    MYCOMPLEX delta1, MYCOMPLEX mu1, 
    MYCOMPLEX delta2, MYCOMPLEX mu2, 
    MYREAL thk, MYREAL k,
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL);    