/**
 * @file   lamb_util.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    一些使用广义闭合解求解 Lamb 问题过程中可能用到的辅助函数
 */

#pragma once

#include <stddef.h>
#include <stdio.h>

#include "grt/common/const.h"

/* 仅用于数值稳定性警告，不改变 Lamb 问题的数学定义域 */
#define LAMB_NU_WARNING_MARGIN 1e-3
#define LAMB_SURFACE_DEPTH_WARNING_RATIO 1e-3

/** 三类 Lamb 解共用的震相掩码 */
typedef enum {
    GRT_LAMB_PHASE_P   = 1u << 0,  ///< 直达 P 波
    GRT_LAMB_PHASE_S   = 1u << 1,  ///< 直达 S 波
    GRT_LAMB_PHASE_R   = 1u << 2,  ///< Rayleigh 波
    GRT_LAMB_PHASE_PP  = 1u << 3,  ///< 自由表面反射 PP 波
    GRT_LAMB_PHASE_SS  = 1u << 4,  ///< 自由表面反射 SS 波
    GRT_LAMB_PHASE_PS  = 1u << 5,  ///< PS 转换波
    GRT_LAMB_PHASE_SP  = 1u << 6,  ///< SP 转换波
    GRT_LAMB_PHASE_SPS = 1u << 7,  ///< sPs 滑行波
} GRT_LAMB_PHASE;

#define GRT_LAMB1_PHASES    (GRT_LAMB_PHASE_P | GRT_LAMB_PHASE_S | GRT_LAMB_PHASE_R)
#define GRT_LAMB2_PHASES    (GRT_LAMB_PHASE_P | GRT_LAMB_PHASE_S | GRT_LAMB_PHASE_PS | GRT_LAMB_PHASE_SP)
#define GRT_LAMB3_PHASES    (GRT_LAMB2_PHASES | GRT_LAMB_PHASE_PP | GRT_LAMB_PHASE_SS | GRT_LAMB_PHASE_SPS)
#define GRT_LAMB_ALL_PHASES (GRT_LAMB1_PHASES | GRT_LAMB3_PHASES)

/**
 * 在参数准备阶段解析震相列表，求解器仅接收解析后的掩码
 *
 * @param[in]  phase_list  以逗号分隔的震相名称，NULL 表示选择全部可用震相
 * @param[in]  supported   当前求解器支持的震相掩码
 * @return 选中的震相掩码，无效名称及重复名称给出警告并忽略，结果为零时输出全零波形
 */
unsigned int grt_lamb_parse_phase_list(const char *phase_list, unsigned int supported);

/** 判断复数的虚部是否可以视为零 */
bool grt_lamb_is_real(const cplx_t value);

/** 对接近零的实数取平方根 */
real_t grt_lamb_positive_sqrt(const real_t value, const char *name);

/** 将椭圆积分参数限制在有效范围内 */
real_t grt_lamb_clamp_elliptic_parameter(const real_t value, const real_t tolerance, const char *name);

/** 求解一元三次方程的三个复根 */
void grt_lamb_cubic_roots(const real_t a, const real_t b, const real_t c, cplx_t roots[3]);

/** 计算以实数为自变量的时间多项式 */
cplx_t grt_lamb_eval_time_coeff(const cplx_t *coefficient, const int degree, const real_t t);

/** 用三个点的二次插值计算一阶导数 */
real_t grt_lamb_derivative_three_points(
    const real_t x0, const real_t x1, const real_t x2,
    const real_t f0, const real_t f1, const real_t f2, const real_t x);

/** 对 Lamb 模块的时间积分项求时间导数 */
void grt_lamb_differentiate_Fk(
    const real_t *ts, const int nt, const real_t (*Fk)[3][3][3],
    real_t (*dG)[3][3][3]);

/** 对 Fkk[time][k][k'][i][j] 连续求两次时间导数 */
void grt_lamb_differentiate_Fkk(
    const real_t *ts, const int nt, const real_t (*Fkk)[3][3][3][3],
    real_t (*dG)[3][3][3][3]);

/** 解析 Lamb 模块的空间导数输出路径 */
void grt_lamb_parse_derivative_paths(const char *argument, char **source_path, char **receiver_path);

/** 解析 Lamb 模块的源点、接收点和混合二阶导数输出路径 */
void grt_lamb_parse_derivative_paths_with_mixed(
    const char *argument, char **source_path, char **receiver_path, char **mixed_path);

/** 输出 Lamb 模块的 Green 函数序列 */
void grt_lamb_print_green_series(FILE *fp, const real_t *ts, const int nt, const real_t (*G)[3][3]);

/** 输出 Lamb 模块的一类空间导数序列 */
void grt_lamb_print_derivative_series(FILE *fp, const real_t *ts, const int nt, const real_t (*dG)[3][3][3], const bool source);

/** 输出混合二阶空间导数序列，dG[time][k][k'][i][j] */
void grt_lamb_print_mixed_derivative_series(
    FILE *fp, const real_t *ts, const int nt, const real_t (*dG)[3][3][3][3]);

/**
 * 求解如下一元三次形式的 Rayleigh 方程的根,  其中 \f$ \nu \f$ 为泊松比
 * \f[
 *       x^3 - \dfrac{2\nu^2 + 1}{2(1 - \nu)} x^2
 *     + \dfrac{4\nu^3 - 4\nu^2 + 4\nu - 1}{4(1 - \nu)^2} x 
 *     - \dfrac{\nu^4}{8(1-\nu)^3} = 0
 * \f]
 * 
 * 
 * @param[in]      nu    泊松比， (0, 0.5)
 * @param[out]     y3    三个根，其中 y3[2] 为正根
 */
void grt_rayleigh1_roots(real_t nu, cplx_t y3[3]);


/**
 * 求解如下一元三次形式的 Rayleigh 方程的根,  其中 \f$ m=\dfrac{1}{2}\dfrac{1-2\nu}{1-\nu}, \nu \f$ 为泊松比
 * \f[
 *       x^3 + \dfrac{2m - 3}{2(1 - m)} x^2
 *     + \dfrac{1}{2(1-m)} x
 *     - \dfrac{1}{16(1-m)} = 0
 * \f]
 * 
 * 
 * @param[in]      m     系数 m
 * @param[out]     y3    三个根，其中 y3[2] 为正根
 */
void grt_rayleigh2_roots(real_t m, cplx_t y3[3]);

/**
 * 做如下多项式求值， \f$ \sum_{m=0}^n C_{2m+o} y^m \f$
 * 
 * @param[in]    C       数组 C
 * @param[in]    n       最高幂次 n
 * @param[in]    y       自变量 y
 * @param[in]    o       偏移量
 * 
 * @return    多项式结果
 * 
 */
cplx_t grt_evalpoly2(const cplx_t *C, const int n, const cplx_t y, const int offset);
