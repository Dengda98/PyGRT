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

#include "grt/common/const.h"
#include "grt/common/RT_matrix.h"
#include "grt/common/model.h"

/**
 * 计算自由表面的静态反射系数，公式(6.3.12)
 * 
 * @param[in]     mod1d             模型结构体指针
 * @param[out]    R_tilt            R/T矩阵，仅填充RU
 * 
 */
void grt_static_topfree_RU(const GRT_MODEL1D *mod1d, RT_MATRIX *M);


/**
 * 计算接收点位置的 P-SV 静态接收矩阵，将波场转为位移，公式(6.3.35)
 * 
 * @param[in]      mod1d           模型结构体指针
 * @param[in]      R               P-SV波场
 * @param[out]     R_EV            P-SV接收函数矩阵
 * 
 */
void grt_static_wave2qwv_REV_PSV(const GRT_MODEL1D *mod1d, const cplx_t R[2][2], cplx_t R_EV[2][2]);

/**
 * 计算接收点位置的 SH 静态接收矩阵，将波场转为位移，公式(6.3.37)
 * 
 * @param[in]      mod1d           模型结构体指针
 * @param[in]      RL              SH波场
 * @param[out]     R_EVL           SH接收函数值
 * 
 */
void grt_static_wave2qwv_REV_SH(const GRT_MODEL1D *mod1d, cplx_t RL, cplx_t *R_EVL);

/**
 * 计算接收点位置的ui_z的 P-SV 静态接收矩阵，即将波场转为ui_z。
 * 公式本质是推导ui_z关于q_m, w_m, v_m的连接矩阵（就是应力推导过程的一部分）
 * 
 * @param[in]      mod1d           模型结构体指针
 * @param[in]      R               P-SV波场
 * @param[out]     R_EV            P-SV接收函数矩阵
 * 
 */
void grt_static_wave2qwv_z_REV_PSV(const GRT_MODEL1D *mod1d, const cplx_t R[2][2], cplx_t R_EV[2][2]);

/**
 * 计算接收点位置的ui_z的 SH 静态接收矩阵，即将波场转为ui_z。
 * 公式本质是推导ui_z关于q_m, w_m, v_m的连接矩阵（就是应力推导过程的一部分）
 * 
 * @param[in]      mod1d           模型结构体指针
 * @param[in]      RL              SH波场
 * @param[out]     R_EVL           SH接收函数值
 * 
 */
void grt_static_wave2qwv_z_REV_SH(const GRT_MODEL1D *mod1d, cplx_t RL, cplx_t *R_EVL);


/**
 * 计算界面的 P-SV 波静态反射透射系数RD/RU/TD/TU
 * 根据公式(6.3.18)  
 * 
 * @param[in]      mod1d         模型结构体指针
 * @param[in]      iy            层位索引
 * @param[out]     M             R/T矩阵
 * 
 */
void grt_static_RT_matrix_PSV(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M);

/**
 * 计算界面的 SH 波静态反射透射系数RDL/RUL/TDL/TUL
 * 根据公式(6.3.18)  
 * 
 * @param[in]      mod1d         模型结构体指针
 * @param[in]      iy            层位索引
 * @param[out]     M             R/T矩阵
 * 
 */
void grt_static_RT_matrix_SH(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M);


/**
 * 为静态 R/T 矩阵添加时间延迟因子
 * 
 * @param[in]      mod1d         模型结构体指针
 * @param[in]      iy            层位索引
 * @param[out]     M             R/T矩阵
 * 
 */
void grt_static_delay_RT_matrix(const GRT_MODEL1D *mod1d, const size_t iy, RT_MATRIX *M);