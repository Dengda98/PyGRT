/**
 * @file   layer.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 P-SV 波和 SH 波的反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient 
 *               matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699
 *            
 */

#pragma once

#include "grt/common/model.h"
#include "grt/common/const.h"

/**
 * 计算自由表面的 P-SV 波反射系数，公式(5.3.10-14) 
 * 
 * @param[in]     xa0            表层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]     xb0            表层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]     cbcb0          相速度与表层的S波速度比值的平方 \f$ c_b^2=(\frac{c}{V_b})^2 \f$
 * @param[in]     k              波数
 * @param[out]    R_tilt         P-SV系数矩阵，SH系数为1
 * @param[out]    stats          状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_topfree_RU_PSV(MYCOMPLEX xa0, MYCOMPLEX xb0, MYCOMPLEX cbcb0, MYREAL k, MYCOMPLEX R_tilt[2][2], MYINT *stats);


/**
 * 计算接收点位置的 P-SV 波接收矩阵，将波场转为位移，公式(5.2.19) + (5.7.7,25)
 * 
 * @param[in]     xa_rcv          接受层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]     xb_rcv          接受层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]     ircvup          接收点是否浅于震源层
 * @param[in]     k               波数
 * @param[in]     R               P-SV波场
 * @param[out]    R_EV            P-SV接收函数矩阵
 * 
 */
void grt_wave2qwv_REV_PSV(
    MYCOMPLEX xa_rcv, MYCOMPLEX xb_rcv, bool ircvup,
    MYREAL k, 
    const MYCOMPLEX R[2][2], MYCOMPLEX R_EV[2][2]);

/**
 * 计算接收点位置的 SH 波接收矩阵，将波场转为位移，公式(5.2.19) + (5.7.7,25)
 * 
 * @param[in]     xb_rcv          接受层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]     k               波数
 * @param[in]     RL              SH波场
 * @param[out]    R_EVL           SH接收函数值
 * 
 */
void grt_wave2qwv_REV_SH(
    MYCOMPLEX xb_rcv,
    MYREAL k, 
    MYCOMPLEX RL, MYCOMPLEX *R_EVL);


/**
 * 计算接收点位置的ui_z的 P-SV 波接收矩阵，即将波场转为ui_z。
 * 公式本质是推导ui_z关于q_m, w_m, v_m的连接矩阵（就是应力推导过程的一部分）
 * 
 * @param[in]     xa_rcv          接受层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]     xb_rcv          接受层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]     ircvup          接收点是否浅于震源层
 * @param[in]     k               波数
 * @param[in]     R               P-SV波场
 * @param[out]    R_EV            P-SV接收函数矩阵
 * 
 */
void grt_wave2qwv_z_REV_PSV(
    MYCOMPLEX xa_rcv, MYCOMPLEX xb_rcv, bool ircvup,
    MYREAL k, 
    const MYCOMPLEX R[2][2], MYCOMPLEX R_EV[2][2]);


/**
 * 计算接收点位置的ui_z的 SH 波接收矩阵，即将波场转为ui_z。
 * 公式本质是推导ui_z关于q_m, w_m, v_m的连接矩阵（就是应力推导过程的一部分）
 * 
 * @param[in]     xb_rcv          接受层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]     ircvup          接收点是否浅于震源层
 * @param[in]     k               波数
 * @param[in]     RL              SH波场
 * @param[out]    R_EVL           SH接收函数值
 * 
 */
void grt_wave2qwv_z_REV_SH(
    MYCOMPLEX xb_rcv, bool ircvup,
    MYREAL k, 
    MYCOMPLEX RL, MYCOMPLEX *R_EVL);


/**
 * 计算界面的 P-SV 波反射透射系数 RD/RU/TD/TU,
 * 根据公式(5.4.14)计算系数   
 * 
 * @note   对公式(5.4.14)进行了重新整理。原公式各项之间的数量级差别过大，浮点数计算损失精度严重。
 * 
 * @param[in]      Rho1          上层的密度
 * @param[in]      xa1           上层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]      xb1           上层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]      cbcb1         相速度与上层的S波速度比值的平方 \f$ c_b^2=(\frac{c}{V_b})^2 \f$
 * @param[in]      mu1           上层的剪切模量
 * @param[in]      Rho2          下层的密度
 * @param[in]      xa2           下层的P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]      xb2           下层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]      cbcb2         相速度与下层的S波速度比值的平方 \f$ c_b^2=(\frac{c}{V_b})^2 \f$
 * @param[in]      mu2           下层的剪切模量
 * @param[in]      omega         角频率
 * @param[in]      k             波数
 * @param[out]     RD            P-SV 下传反射系数矩阵
 * @param[out]     RU            P-SV 上传反射系数矩阵
 * @param[out]     TD            P-SV 下传透射系数矩阵
 * @param[out]     TU            P-SV 上传透射系数矩阵
 * @param[out]     stats         状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_RT_matrix_PSV(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX cbcb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX cbcb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2],
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2], MYINT *stats);


/**
 * 计算界面的 SH 波反射透射系数 RDL/RUL/TDL/TUL,
 * 根据公式(5.4.31)计算系数   
 * 
 * @param[in]      xb1           上层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]      mu1           上层的剪切模量
 * @param[in]      xb2           下层的S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]      mu2           下层的剪切模量
 * @param[in]      omega         角频率
 * @param[in]      k             波数
 * @param[out]     RDL           SH 下传反射系数
 * @param[out]     RUL           SH 上传反射系数
 * @param[out]     TDL           SH 下传透射系数
 * @param[out]     TUL           SH 上传透射系数
 * @param[out]     stats         状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_RT_matrix_SH(
    MYCOMPLEX xb1, MYCOMPLEX mu1, 
    MYCOMPLEX xb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX *RDL, MYCOMPLEX *RUL, 
    MYCOMPLEX *TDL, MYCOMPLEX *TUL);

/** 液-液 界面，函数参数见 calc_RT_PSV 函数 */
void grt_RT_matrix_ll_PSV(
    MYREAL Rho1, MYCOMPLEX xa1,
    MYREAL Rho2, MYCOMPLEX xa2,
    MYCOMPLEX omega, MYREAL k,
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2], 
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2], MYINT *stats);

/** 液-液 界面，函数参数见 calc_RT_SH 函数 */
void grt_RT_matrix_ll_SH(
    MYCOMPLEX *RDL, MYCOMPLEX *RUL, 
    MYCOMPLEX *TDL, MYCOMPLEX *TUL);

/** 液-固 界面，函数参数见 calc_RT_PSV 函数 */
void grt_RT_matrix_ls_PSV(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX cbcb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX cbcb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2], 
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2], MYINT *stats);

/** 液-固 界面，函数参数见 calc_RT_SH 函数 */
void grt_RT_matrix_ls_SH(
    MYCOMPLEX xb1, MYCOMPLEX mu1, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k,
    MYCOMPLEX *RDL, MYCOMPLEX *RUL, 
    MYCOMPLEX *TDL, MYCOMPLEX *TUL);

/** 固-固 界面，函数参数见 calc_RT_PSV 函数 */
void grt_RT_matrix_ss_PSV(
    MYREAL Rho1, MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX cbcb1, MYCOMPLEX mu1, 
    MYREAL Rho2, MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX cbcb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2],
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2], MYINT *stats);

/** 固-固 界面，函数参数见 calc_RT_SH 函数 */
void grt_RT_matrix_ss_SH(
    MYCOMPLEX xb1, MYCOMPLEX mu1, 
    MYCOMPLEX xb2, MYCOMPLEX mu2, 
    MYCOMPLEX omega, MYREAL k, 
    MYCOMPLEX *RDL, MYCOMPLEX *RUL, 
    MYCOMPLEX *TDL, MYCOMPLEX *TUL);

/**
 * 为 P-SV 波的 R/T 矩阵添加时间延迟因子
 * 
 * @param[in]      xa1            P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]      xb1            S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]      thk            厚度
 * @param[in]      k              波数
 * @param[in,out]     RD          P-SV 下传反射系数矩阵
 * @param[in,out]     RU          P-SV 上传反射系数矩阵
 * @param[in,out]     TD          P-SV 下传透射系数矩阵
 * @param[in,out]     TU          P-SV 上传透射系数矩阵     
 * 
 */
void grt_delay_RT_matrix_PSV(
    MYCOMPLEX xa1, MYCOMPLEX xb1, 
    MYREAL thk, MYREAL k,
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2], 
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2]);

/**
 * 为 SH 波的 R/T 矩阵添加时间延迟因子
 * 
 * @param[in]      xb1            S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]      thk            厚度
 * @param[in]      k              波数
 * @param[in,out]     RDL         SH 下传反射系数矩阵
 * @param[in,out]     RUL         SH 上传反射系数矩阵
 * @param[in,out]     TDL         SH 下传透射系数矩阵
 * @param[in,out]     TUL         SH 上传透射系数矩阵     
 * 
 */
void grt_delay_RT_matrix_SH(
    MYCOMPLEX xb1, 
    MYREAL thk, MYREAL k,
    MYCOMPLEX *RDL, MYCOMPLEX *RUL, 
    MYCOMPLEX *TDL, MYCOMPLEX *TUL);


/**
 * 计算该层的连接 P-SV 应力位移矢量与垂直波函数的D矩阵(或其逆矩阵)，
 * 见公式(5.2.19-20)
 * 
 * @param[in]      xa            P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]      xb            S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]      kbkb          S波水平波数的平方 \f$ k_b^2=(\frac{\omega}{V_b})^2 \f$
 * @param[in]      mu            剪切模量
 * @param[in]      omega         角频率
 * @param[in]      rho           密度
 * @param[in]      k             波数
 * @param[out]     D             D矩阵(或其逆矩阵)
 * @param[in]      inverse       是否生成逆矩阵
 * @param[in]      liquid_invtype   对于液体层，矩阵会有很多零，至少第二列、第四列和第四行均为零；
 *                                  剩余部分根据所选类型进行讨论：
 *                                  [1] 其余6项保留， \f$ 2\mu\Omega \f$ 退化为 \f$ - \rho \omega^2 \f$ ;
 *                                  [2] 在 [1] 基础上第一行也置零，这用于满足液体层的边界条件；
 *                                  对应逆矩阵使用伪逆。 
 * 
 */
void grt_get_layer_D(
    MYCOMPLEX xa, MYCOMPLEX xb, MYCOMPLEX kbkb, MYCOMPLEX mu, 
    MYCOMPLEX omega, MYREAL rho, MYREAL k, MYCOMPLEX D[4][4], bool inverse, MYINT liquid_invtype);

/** 子矩阵 D11，函数参数见 get_layer_D 函数 */
void grt_get_layer_D11(
    MYCOMPLEX xa, MYCOMPLEX xb, MYREAL k, MYCOMPLEX D[2][2]);

/** 子矩阵 D12，函数参数见 get_layer_D 函数 */
void grt_get_layer_D12(
    MYCOMPLEX xa, MYCOMPLEX xb, MYREAL k, MYCOMPLEX D[2][2]);

/** 子矩阵 D21，函数参数见 get_layer_D 函数 */
void grt_get_layer_D21(
    MYCOMPLEX xa, MYCOMPLEX xb, MYCOMPLEX kbkb, MYCOMPLEX mu,
    MYCOMPLEX omega, MYREAL rho, MYREAL k, MYCOMPLEX D[2][2]);

/** 子矩阵 D22，函数参数见 get_layer_D 函数 */
void grt_get_layer_D22(
    MYCOMPLEX xa, MYCOMPLEX xb, MYCOMPLEX kbkb, MYCOMPLEX mu,
    MYCOMPLEX omega, MYREAL rho, MYREAL k, MYCOMPLEX D[2][2]);

/** 子矩阵 D11_uiz，后缀uiz表示连接位移对z的偏导和垂直波函数，函数参数见 get_layer_D 函数 */
void grt_get_layer_D11_uiz(
    MYCOMPLEX xa, MYCOMPLEX xb, MYREAL k, MYCOMPLEX D[2][2]);

/** 子矩阵 D12_uiz，函数参数见 get_layer_D 函数 */
void grt_get_layer_D12_uiz(
    MYCOMPLEX xa, MYCOMPLEX xb, MYREAL k, MYCOMPLEX D[2][2]);


/**
 * 计算该层的连接 SH 应力位移矢量与垂直波函数的 T 矩阵(或其逆矩阵)，
 * 见公式(5.2.21-22)
 * 
 * @param[in]      xb            S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]      mu            剪切模量
 * @param[in]      omega         角频率
 * @param[in]      k             波数
 * @param[out]     T             T矩阵(或其逆矩阵)
 * @param[in]      inverse       是否生成逆矩阵
 * 
 */
void grt_get_layer_T(
    MYCOMPLEX xb, MYCOMPLEX mu,
    MYCOMPLEX omega, MYREAL k, MYCOMPLEX T[2][2], bool inverse);

/** 计算 P-SV 型垂直波函数的时间延迟矩阵，公式(5.2.27) */
void grt_get_layer_E_Rayl(MYCOMPLEX xa1, MYCOMPLEX xb1, MYREAL thk, MYREAL k, MYCOMPLEX E[4][4], bool inverse);

/** 计算 SH 型垂直波函数的时间延迟矩阵，公式(5.2.28) */
void grt_get_layer_E_Love(MYCOMPLEX xb1, MYREAL thk, MYREAL k, MYCOMPLEX E[2][2], bool inverse);



/**
 *  【未维护，未使用，仅用于内部代码测试】
 *  和 calc_RT_PSV(SH) 函数解决相同问题，但没有使用显式推导的公式，而是直接做矩阵运算，
 *  函数接口也类似
 */
void grt_RT_matrix_from_4x4(
    MYCOMPLEX xa1, MYCOMPLEX xb1, MYCOMPLEX kbkb1, MYCOMPLEX mu1, MYREAL rho1, 
    MYCOMPLEX xa2, MYCOMPLEX xb2, MYCOMPLEX kbkb2, MYCOMPLEX mu2, MYREAL rho2,
    MYCOMPLEX omega, MYREAL thk,
    MYREAL k, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats);