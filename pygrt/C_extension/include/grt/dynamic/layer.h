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
#include "grt/common/RT_matrix.h"

/**
 * 计算自由表面的反射系数，公式(5.3.10-14) 
 * 
 * 
 * @param[in,out]     mod1d          模型结构体指针，结果保存在 Mtop
 * 
 */
void grt_topfree_RU(GRT_MODEL1D *mod1d);


/**
 * 计算接收点位置的 P-SV 波接收矩阵，将波场转为位移，公式(5.2.19) + (5.7.7,25)
 * 
 * @param[in,out]     mod1d           模型结构体指针，结果保存在 R_EV
 * 
 */
void grt_wave2qwv_REV_PSV(GRT_MODEL1D *mod1d);

/**
 * 计算接收点位置的 SH 波接收矩阵，将波场转为位移，公式(5.2.19) + (5.7.7,25)
 * 
 * @param[in,out]     mod1d           模型结构体指针，结果保存在 R_EVL
 * 
 */
void grt_wave2qwv_REV_SH(GRT_MODEL1D *mod1d);


/**
 * 计算接收点位置的ui_z的 P-SV 波接收矩阵，即将波场转为ui_z。
 * 公式本质是推导ui_z关于q_m, w_m, v_m的连接矩阵（就是应力推导过程的一部分）
 * 
 * @param[in,out]     mod1d           模型结构体指针，结果保存在 uiz_R_EV
 * 
 */
void grt_wave2qwv_z_REV_PSV(GRT_MODEL1D *mod1d);


/**
 * 计算接收点位置的ui_z的 SH 波接收矩阵，即将波场转为ui_z。
 * 公式本质是推导ui_z关于q_m, w_m, v_m的连接矩阵（就是应力推导过程的一部分）
 * 
 * @param[in,out]     mod1d           模型结构体指针，结果保存在 uiz_R_EVL
 * 
 */
void grt_wave2qwv_z_REV_SH(GRT_MODEL1D *mod1d);


/**
 * 计算界面的 P-SV 波反射透射系数 RD/RU/TD/TU,
 * 根据公式(5.4.14)计算系数   
 * 
 * @note   对公式(5.4.14)进行了重新整理。原公式各项之间的数量级差别过大，浮点数计算损失精度严重。
 * 
 * @param[in]      mod1d         模型结构体指针
 * @param[in]      iy            层位索引
 * @param[out]     M             R/T矩阵
 * 
 */
void grt_RT_matrix_PSV(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M);


/**
 * 计算界面的 SH 波反射透射系数 RDL/RUL/TDL/TUL,
 * 根据公式(5.4.31)计算系数   
 * 
 * @param[in]      mod1d         模型结构体指针
 * @param[in]      iy            层位索引
 * @param[out]     M             R/T矩阵
 * 
 */
void grt_RT_matrix_SH(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M);

/** 液-液 界面 */
void grt_RT_matrix_ll_PSV(const GRT_MODEL1D *mod1d, size_t iy, RT_MATRIX *M);

/** 液-液 界面 */
void grt_RT_matrix_ll_SH(RT_MATRIX *M);

/** 液-固 界面 */
void grt_RT_matrix_ls_PSV(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M);

/** 液-固 界面 */
void grt_RT_matrix_ls_SH(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M);

/** 固-固 界面 */
void grt_RT_matrix_ss_PSV(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M);

/** 固-固 界面 */
void grt_RT_matrix_ss_SH(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M);

/**
 * 为 R/T 矩阵添加时间延迟因子
 * 
 * @param[in]      mod1d         模型结构体指针
 * @param[in]      iy            层位索引
 * @param[out]     M             R/T矩阵    
 * 
 */
void grt_delay_RT_matrix(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M);

/**
 * 为虚拟层的广义 R/T 矩阵添加时间延迟因子
 * 
 * @param[in]      mod1d         模型结构体指针
 * @param[in]      iy            层位索引
 * @param[out]     M             R/T矩阵    
 * 
 */
void grt_delay_GRT_matrix(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M);

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
    cplx_t xa, cplx_t xb, cplx_t kbkb, cplx_t mu, 
    cplx_t omega, real_t rho, real_t k, cplx_t D[4][4], bool inverse, int liquid_invtype);

/** 子矩阵 D11，函数参数见 get_layer_D 函数 */
void grt_get_layer_D11(
    cplx_t xa, cplx_t xb, real_t k, cplx_t D[2][2]);

/** 子矩阵 D12，函数参数见 get_layer_D 函数 */
void grt_get_layer_D12(
    cplx_t xa, cplx_t xb, real_t k, cplx_t D[2][2]);

/** 子矩阵 D21，函数参数见 get_layer_D 函数 */
void grt_get_layer_D21(
    cplx_t xa, cplx_t xb, cplx_t kbkb, cplx_t mu,
    cplx_t omega, real_t rho, real_t k, cplx_t D[2][2]);

/** 子矩阵 D22，函数参数见 get_layer_D 函数 */
void grt_get_layer_D22(
    cplx_t xa, cplx_t xb, cplx_t kbkb, cplx_t mu,
    cplx_t omega, real_t rho, real_t k, cplx_t D[2][2]);

/** 子矩阵 D11_uiz，后缀uiz表示连接位移对z的偏导和垂直波函数，函数参数见 get_layer_D 函数 */
void grt_get_layer_D11_uiz(
    cplx_t xa, cplx_t xb, real_t k, cplx_t D[2][2]);

/** 子矩阵 D12_uiz，函数参数见 get_layer_D 函数 */
void grt_get_layer_D12_uiz(
    cplx_t xa, cplx_t xb, real_t k, cplx_t D[2][2]);


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
    cplx_t xb, cplx_t mu,
    cplx_t omega, real_t k, cplx_t T[2][2], bool inverse);

/** 计算 P-SV 型垂直波函数的时间延迟矩阵，公式(5.2.27) */
void grt_get_layer_E_Rayl(cplx_t xa1, cplx_t xb1, real_t thk, real_t k, cplx_t E[4][4], bool inverse);

/** 计算 SH 型垂直波函数的时间延迟矩阵，公式(5.2.28) */
void grt_get_layer_E_Love(cplx_t xb1, real_t thk, real_t k, cplx_t E[2][2], bool inverse);



/**
 *  【未维护，未使用，仅用于内部代码测试】
 *  和 calc_RT_PSV(SH) 函数解决相同问题，但没有使用显式推导的公式，而是直接做矩阵运算，
 *  函数接口也类似
 */
void grt_RT_matrix_from_4x4(
    cplx_t xa1, cplx_t xb1, cplx_t kbkb1, cplx_t mu1, real_t rho1, 
    cplx_t xa2, cplx_t xb2, cplx_t kbkb2, cplx_t mu2, real_t rho2,
    cplx_t omega, real_t thk,
    real_t k, 
    cplx_t RD[2][2], cplx_t *RDL, cplx_t RU[2][2], cplx_t *RUL, 
    cplx_t TD[2][2], cplx_t *TDL, cplx_t TU[2][2], cplx_t *TUL, int *stats);