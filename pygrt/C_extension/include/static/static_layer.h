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

#include "common/const.h"

/**
 * 计算自由表面的静态反射系数，公式(6.3.12)
 * 
 * @param     delta1           (in)表层的 \f$ \Delta = \frac{\lambda + \mu}{\lambda + 3\mu} \f$
 * @param     R_tilt[2][2]     (out)P-SV系数矩阵，SH系数为1
 * 
 */
void calc_static_R_tilt(MYCOMPLEX delta1, MYCOMPLEX R_tilt[2][2]);


/**
 * 计算接收点位置的静态接收矩阵，将波场转为位移，公式(6.3.35,37)
 * 
 * @param     ircvup      (in)接收点是否浅于震源层
 * @param     k           (in)波数
 * @param     R[2][2]     (in)P-SV波场
 * @param     RL          (in)SH波场
 * @param     R_EV[2][2]  (out)P-SV接收函数矩阵
 * @param     R_EVL       (out)SH接收函数值
 * 
 */
void calc_static_R_EV(
    bool ircvup,
    const MYCOMPLEX R[2][2], MYCOMPLEX RL, 
    MYCOMPLEX R_EV[2][2], MYCOMPLEX *R_EVL);


/**
 * 计算接收点位置的ui_z的静态接收矩阵，即将波场转为ui_z。
 * 公式本质是推导ui_z关于q_m, w_m, v_m的连接矩阵（就是应力推导过程的一部分）
 * 
 * @param     delta1      (in)接收层的 \f$ \Delta \f$
 * @param     ircvup      (in)接收点是否浅于震源层
 * @param     k           (in)波数
 * @param     R[2][2]     (in)P-SV波场
 * @param     RL          (in)SH波场
 * @param     R_EV[2][2]  (out)P-SV接收函数矩阵
 * @param     R_EVL       (out)SH接收函数值
 * 
 */
void calc_static_uiz_R_EV(
    MYCOMPLEX delta1, bool ircvup, MYREAL k, 
    const MYCOMPLEX R[2][2], MYCOMPLEX RL, 
    MYCOMPLEX R_EV[2][2], MYCOMPLEX *R_EVL);


/**
 * 计算界面的静态反射系数RD/RDL/RU/RUL, 静态透射系数TD/TDL/TU/TUL, 包括时间延迟因子，
 * 后缀L表示SH波的系数, 其余表示P-SV波的系数, 根据公式(6.3.18)  
 * 
 * @param      delta1    (in)上层的 \f$ \Delta \f$
 * @param      mu1       (in)上层的剪切模量
 * @param      delta2    (in)下层的 \f$ \Delta \f$
 * @param      mu2       (in)下层的剪切模量
 * @param      thk       (in)上层层厚
 * @param      k         (in)波数
 * @param      RD[2][2]  (out)P-SV 下传反射系数矩阵
 * @param      RDL       (out)SH 下传反射系数
 * @param      RU[2][2]  (out)P-SV 上传反射系数矩阵
 * @param      RUL       (out)SH 上传反射系数
 * @param      TD[2][2]  (out)P-SV 下传透射系数矩阵
 * @param      TDL       (out)SH 下传透射系数
 * @param      TU[2][2]  (out)P-SV 上传透射系数矩阵
 * @param      TUL       (out)SH 上传透射系数
 * 
 */
void calc_static_RT_2x2(
    MYCOMPLEX delta1, MYCOMPLEX mu1, 
    MYCOMPLEX delta2, MYCOMPLEX mu2, 
    MYREAL thk, MYREAL k,
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL);    