/**
 * @file   layer.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient 
 *               matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699
 *            
 */

#pragma once

#include "model.h"
#include "const.h"

/**
 * 计算自由表面的反射系数，公式(5.3.10-14) 
 * 
 * @param     xa0        (in)表层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param     xb0        (in)表层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param     kbkb0     (in)表层的S波水平波数的平方 \f$ k_b^2=(\frac{\omega}{V_b})^2 \f$
 * @param     k         (in)波数
 * @param     R_tilt[2][2]     (out)P-SV系数矩阵，SH系数为1
 * 
 */
void calc_R_tilt(
    MYCOMPLEX xa0, MYCOMPLEX xb0, MYCOMPLEX kbkb0, MYREAL k, MYCOMPLEX R_tilt[2][2]);


/**
 * 计算接收点位置的接收矩阵，将波场转为位移，公式(5.2.19) + (5.7.7,25)
 * 
 * @param     xa_rcv      (in)接受层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param     xb_rcv      (in)接受层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param     ircvup      (in)接收点是否浅于震源层
 * @param     k           (in)波数
 * @param     R[2][2]     (in)P-SV波场
 * @param     RL          (in)SH波场
 * @param     R_EV[2][2]  (out)P-SV接收函数矩阵
 * @param     R_EVL       (out)SH接收函数值
 * 
 */
void calc_R_EV(
    MYCOMPLEX xa_rcv, MYCOMPLEX xb_rcv, bool ircvup,
    MYREAL k, 
    const MYCOMPLEX R[2][2], MYCOMPLEX RL, 
    MYCOMPLEX R_EV[2][2], MYCOMPLEX *R_EVL);


/**
 * 计算界面的反射系数RD/RDL/RU/RUL, 透射系数TD/TDL/TU/TUL, 包括时间延迟因子，
 * 后缀L表示SH波的系数, 其余表示P-SV波的系数, 根据公式(5.4.14)和(5.4.31)计算系数   
 * 
 * @note   对公式(5.4.14)进行了重新整理。原公式各项之间的数量级差别过大，浮点数计算损失精度严重。
 * 
 * @param      Rho1      (in)上层的密度
 * @param      xa1       (in)上层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param      xb1       (in)上层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param      kbkb1     (in)上层的S波水平波数的平方 \f$ k_b^2=(\frac{\omega}{V_b})^2 \f$
 * @param      mu1       (in)上层的剪切模量
 * @param      Rho2      (in)下层的密度
 * @param      xa2       (in)下层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param      xb2       (in)下层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param      kbkb2     (in)下层的S波水平波数的平方 \f$ k_b^2=(\frac{\omega}{V_b})^2 \f$
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
void calc_RT_2x2(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX kbkb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX kbkb2, MYCOMPLEX mu2, 
    MYREAL thk,
    MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL);



/**
 * 【未使用】
 * 被calc_RT_2x2_from_4x4函数调用，生成该层的连接P-SV应力位移矢量与垂直波函数的D矩阵(或其逆矩阵)，
 * 见公式(5.2.19-20)
 * 
 * @param      xa        (in)P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param      xb        (in)S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param      kbkb      (in)S波水平波数的平方 \f$ k_b^2=(\frac{\omega}{V_b})^2 \f$
 * @param      mu        (in)剪切模量
 * @param      omega     (in)复数频率  \f$ \tilde{\omega} =\omega - i\omega_I \f$ 
 * @param      k         (in)波数
 * 
 * @param      D[4][4]   (out) D矩阵(或其逆矩阵)
 * @param      inverse   (in)  是否生成逆矩阵
 * 
 */
void get_layer_D(
    MYCOMPLEX xa, MYCOMPLEX xb, MYCOMPLEX kbkb, MYCOMPLEX mu, 
    MYCOMPLEX omega, MYREAL k, MYCOMPLEX D[4][4], bool inverse);



/**
 *  和calc_RT_2x2函数解决相同问题（但未包含时间延迟因子），但没有使用显式推导的公式，而是直接做矩阵运算，
 * 暂未加入相位延迟矩阵
 * 函数接口也和 calc_RT_2x2函数 类似
 */
void calc_RT_2x2_from_4x4(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX kbkb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX kbkb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL);
