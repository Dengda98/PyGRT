/**
 * @file   lamb2.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 *    使用广义闭合解求解第二类 Lamb 问题，参考：
 *
 *        张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社
 *
 *    本实现直接使用预先推导的多项式系数，避免运行时构造多项式
 */

#pragma once

#include "grt/common/const.h"

/**
 * 使用广义闭合解求解第二类 Lamb 问题
 *
 * @param[in]    nu             泊松比，(0, 0.5)，距任一边界小于 1e-3 时给出警告
 * @param[in]    ts             无量纲时间序列 tbar=t/(r/beta)=beta*t/r
 * @param[in]    nt             时间序列点数
 * @param[in]    R              源点到接收点的水平距离，必须为正数
 * @param[in]    source_depth   源点深度，与 receiver_depth 恰好一个大于零
 * @param[in]    receiver_depth 接收点深度，与 source_depth 恰好一个大于零
 * @param[in]    azimuth        方位角，单位度，[0, 360]
 * @param[out]   G              无量纲阶跃力位移，G[time][i][j]
 * @param[out]   dG_source      无量纲源点导数，dG_source[time][k'][i][j]
 * @param[out]   dG_receiver    无量纲接收点导数，dG_receiver[time][k][i][j]
 *
 * 当三个输出指针同时为 NULL 时，仅将时间和 G 输出到标准输出。
 */
void grt_solve_lamb2(
    const real_t nu, const real_t *ts, const int nt,
    const real_t R, const real_t source_depth, const real_t receiver_depth, const real_t azimuth,
    real_t (*G)[3][3], real_t (*dG_source)[3][3][3], real_t (*dG_receiver)[3][3][3]);
