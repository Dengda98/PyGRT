/**
 * @file   lamb3.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-09
 *
 *    使用广义闭合解求解第三类 Lamb 问题，参考：
 *
 *        张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社
 */

#pragma once

#include "grt/common/const.h"


/** 第三类 Lamb 解中可单独选择的震相项 */
typedef enum {
    GRT_LAMB3_PHASE_P  = 1u << 0,  ///< 直达 P 波
    GRT_LAMB3_PHASE_S  = 1u << 1,  ///< 直达 S 波
    GRT_LAMB3_PHASE_PP = 1u << 2,  ///< 自由表面反射 PP 波
    GRT_LAMB3_PHASE_SS = 1u << 3,  ///< 自由表面反射 SS 波
    GRT_LAMB3_PHASE_PS = 1u << 4,  ///< 自由表面反射 PS 转换波
    GRT_LAMB3_PHASE_SP = 1u << 5,  ///< 自由表面反射 SP 转换波
    GRT_LAMB3_PHASE_SPS = 1u << 6,  ///< sPs 滑行波
} GRT_LAMB3_PHASE;


/**
 * 计算第三类 Lamb 问题的无量纲震相到时
 *
 * @param[in]    nu             泊松比，(0, 0.5)
 * @param[in]    R              源点和接收点之间的水平距离
 * @param[in]    depsrc         源点深度
 * @param[in]    deprcv         接收点深度
 * @param[out]   tP             直达 P 波到时
 * @param[out]   tPP            反射 PP 波到时
 * @param[out]   tSS            反射 SS 波到时
 * @param[out]   tPS            PS 转换波到时
 * @param[out]   tSP            SP 转换波到时
 * @param[out]   t_sPs          sPs 波到时，不存在时为负数
 */
void grt_compute_lamb3_travt(
    const real_t nu, const real_t R, const real_t depsrc, const real_t deprcv,
    real_t *tP, real_t *tPP, real_t *tSS, real_t *tPS, real_t *tSP, real_t *t_sPs);


/**
 * 使用广义闭合解求解第三类 Lamb 问题
 *
 * @param[in]    nu              泊松比，(0, 0.5)，距任一边界小于 1e-3 时给出警告
 * @param[in]    ts              无量纲时间序列 tbar=t/(r/beta)=beta*t/r
 * @param[in]    nt              时间序列点数
 * @param[in]    R               源点和接收点之间的水平距离，应为正数
 * @param[in]    depsrc          源点深度 x3'，必须大于零
 * @param[in]    deprcv          接收点深度 x3，必须大于零
 * @param[in]    azimuth         水平距离方位角，单位度，[0, 360]
 * @param[in]    phase_list      以逗号分隔的震相名称，NULL 表示选择全部震相
 * @param[out]   G               无量纲阶跃力位移，G[time][i][j]
 * @param[out]   dG_source       无量纲源点导数，dG_source[time][k'][i][j]
 * @param[out]   dG_receiver     无量纲接收点导数，dG_receiver[time][k][i][j]
 * @param[out]   dG_mixed        无量纲混合二阶导数，dG_mixed[time][k][k'][i][j]
 *
 * 当四个输出指针同时为 NULL 时，仅将时间和 G 输出到标准输出
 * 支持的震相为 P、S、PP、SS、PS、SP 和 sPs，SS 与 sPs 可以分别选择
 * 其中 k 为接收点方向，k' 为源点方向。使用 phase_list 后如果没有有效震相，输出全零波形并给出警告
 */
void grt_solve_lamb3(
    const real_t nu, const real_t *ts, const int nt,
    const real_t R, const real_t depsrc, const real_t deprcv,
    const real_t azimuth, const char *phase_list, real_t (*G)[3][3], real_t (*dG_source)[3][3][3],
    real_t (*dG_receiver)[3][3][3], real_t (*dG_mixed)[3][3][3][3]);
