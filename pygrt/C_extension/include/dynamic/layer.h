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

#include "common/model.h"
#include "common/const.h"

/**
 * 计算自由表面的反射系数，公式(5.3.10-14) 
 * 
 * @param[in]     xa0            表层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]     xb0            表层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]     kbkb0          表层的S波水平波数的平方 \f$ k_b^2=(\frac{\omega}{V_b})^2 \f$
 * @param[in]     k              波数
 * @param[out]    R_tilt         P-SV系数矩阵，SH系数为1
 * @param[out]    stats          状态代码，是否有除零错误，非0为异常值
 * 
 */
void calc_R_tilt(
    MYCOMPLEX xa0, MYCOMPLEX xb0, MYCOMPLEX kbkb0, MYREAL k, MYCOMPLEX R_tilt[2][2], MYINT *stats);


/**
 * 计算接收点位置的接收矩阵，将波场转为位移，公式(5.2.19) + (5.7.7,25)
 * 
 * @param[in]     xa_rcv          接受层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]     xb_rcv          接受层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]     ircvup          接收点是否浅于震源层
 * @param[in]     k               波数
 * @param[in]     R               P-SV波场
 * @param[in]     RL              SH波场
 * @param[out]    R_EV            P-SV接收函数矩阵
 * @param[out]    R_EVL           SH接收函数值
 * 
 */
void calc_R_EV(
    MYCOMPLEX xa_rcv, MYCOMPLEX xb_rcv, bool ircvup,
    MYREAL k, 
    const MYCOMPLEX R[2][2], MYCOMPLEX RL, 
    MYCOMPLEX R_EV[2][2], MYCOMPLEX *R_EVL);


/**
 * 计算接收点位置的ui_z的接收矩阵，即将波场转为ui_z。
 * 公式本质是推导ui_z关于q_m, w_m, v_m的连接矩阵（就是应力推导过程的一部分）
 * 
 * @param[in]     xa_rcv          接受层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]     xb_rcv          接受层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]     ircvup          接收点是否浅于震源层
 * @param[in]     k               波数
 * @param[in]     R               P-SV波场
 * @param[in]     RL              SH波场
 * @param[out]    R_EV            P-SV接收函数矩阵
 * @param[out]    R_EVL           SH接收函数值
 * 
 */
void calc_uiz_R_EV(
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
 * @param[in]      Rho1          上层的密度
 * @param[in]      xa1           上层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]      xb1           上层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]      kbkb1         上层的S波水平波数的平方 \f$ k_b^2=(\frac{\omega}{V_b})^2 \f$
 * @param[in]      mu1           上层的剪切模量
 * @param[in]      Rho2          下层的密度
 * @param[in]      xa2           下层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]      xb2           下层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]      kbkb2         下层的S波水平波数的平方 \f$ k_b^2=(\frac{\omega}{V_b})^2 \f$
 * @param[in]      mu2           下层的剪切模量
 * @param[in]      thk           上层层厚
 * @param[in]      omega         角频率
 * @param[in]      k             波数
 * @param[out]     RD            P-SV 下传反射系数矩阵
 * @param[out]     RDL           SH 下传反射系数
 * @param[out]     RU            P-SV 上传反射系数矩阵
 * @param[out]     RUL           SH 上传反射系数
 * @param[out]     TD            P-SV 下传透射系数矩阵
 * @param[out]     TDL           SH 下传透射系数
 * @param[out]     TU            P-SV 上传透射系数矩阵
 * @param[out]     TUL           SH 上传透射系数
 * @param[out]     stats         状态代码，是否有除零错误，非0为异常值
 * 
 */
void calc_RT_2x2(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX kbkb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX kbkb2, MYCOMPLEX mu2, 
    MYREAL thk, MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats);

/** 固-固 界面，函数参数与 calc_RT_2x2 函数一致 */
void calc_RT_ss_2x2(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX kbkb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX kbkb2, MYCOMPLEX mu2, 
    MYREAL thk, MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats);

/** 液-液 界面，函数参数与 calc_RT_2x2 函数一致 */
void calc_RT_ll_2x2(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX kbkb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX kbkb2, MYCOMPLEX mu2, 
    MYREAL thk, MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats);

/** 液-固 界面，函数参数与 calc_RT_2x2 函数一致 */
void calc_RT_ls_2x2(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX kbkb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX kbkb2, MYCOMPLEX mu2, 
    MYREAL thk, MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats);

/**
 * 【未使用，仅用于代码测试】
 * 被calc_RT_2x2_from_4x4函数调用，生成该层的连接P-SV应力位移矢量与垂直波函数的D矩阵(或其逆矩阵)，
 * 见公式(5.2.19-20)
 * 
 * @param[in]      xa            P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]      xb            S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]      kbkb          S波水平波数的平方 \f$ k_b^2=(\frac{\omega}{V_b})^2 \f$
 * @param[in]      mu            剪切模量
 * @param[in]      k             波数
 * 
 * @param[out]      D                D矩阵(或其逆矩阵)
 * @param[in]       inverse          是否生成逆矩阵
 * 
 */
void get_layer_D(
    MYCOMPLEX xa, MYCOMPLEX xb, MYCOMPLEX kbkb, MYCOMPLEX mu, 
    MYCOMPLEX omega, MYREAL k, MYCOMPLEX D[4][4], bool inverse);



/**
 *  【未使用，仅用于代码测试】
 *  和calc_RT_2x2函数解决相同问题，但没有使用显式推导的公式，而是直接做矩阵运算，
 *  函数接口也和 calc_RT_2x2函数 类似
 */
void calc_RT_2x2_from_4x4(
    MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX kbkb1, MYCOMPLEX mu1, 
    MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX kbkb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL thk,
    MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats);
