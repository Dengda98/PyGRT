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
 * @param[in]     delta1            表层的 \f$ \Delta = \frac{\lambda + \mu}{\lambda + 3\mu} \f$
 * @param[out]    R_tilt            P-SV系数矩阵，SH系数为1
 * 
 */
void calc_static_R_tilt(MYCOMPLEX delta1, MYCOMPLEX R_tilt[2][2]);


/**
 * 计算接收点位置的静态接收矩阵，将波场转为位移，公式(6.3.35,37)
 * 
 * @param[in]      ircvup          接收点是否浅于震源层
 * @param[in]      k               波数
 * @param[in]      R               P-SV波场
 * @param[in]      RL              SH波场
 * @param[out]     R_EV            P-SV接收函数矩阵
 * @param[out]     R_EVL           SH接收函数值
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
 * @param[in]      delta1          接收层的 \f$ \Delta \f$
 * @param[in]      ircvup          接收点是否浅于震源层
 * @param[in]      k               波数
 * @param[in]      R               P-SV波场
 * @param[in]      RL              SH波场
 * @param[out]     R_EV            P-SV接收函数矩阵
 * @param[out]     R_EVL           SH接收函数值
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
 * @param[in]       delta1        上层的 \f$ \Delta \f$
 * @param[in]       mu1           上层的剪切模量
 * @param[in]       delta2        下层的 \f$ \Delta \f$
 * @param[in]       mu2           下层的剪切模量
 * @param[in]       thk           上层层厚
 * @param[in]       k             波数
 * @param[out]      RD            P-SV 下传反射系数矩阵
 * @param[out]      RDL           SH 下传反射系数
 * @param[out]      RU            P-SV 上传反射系数矩阵
 * @param[out]      RUL           SH 上传反射系数
 * @param[out]      TD            P-SV 下传透射系数矩阵
 * @param[out]      TDL           SH 下传透射系数
 * @param[out]      TU            P-SV 上传透射系数矩阵
 * @param[out]      TUL           SH 上传透射系数
 * 
 */
void calc_static_RT_2x2(
    MYCOMPLEX delta1, MYCOMPLEX mu1, 
    MYCOMPLEX delta2, MYCOMPLEX mu2, 
    MYREAL thk, MYREAL k,
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL);    