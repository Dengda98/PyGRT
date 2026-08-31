/**
 * @file   lamb2.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 *    使用广义闭合解求解第二类 Lamb 问题，参考：
 *
 *        张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社
 */

#pragma once

#include "grt/common/const.h"


/**
 * 使用广义闭合解求解第二类 Lamb 问题
 *
 * 这里的观测点位于地表，震源位于地下。theta 是从震源指向观测点的
 * 射线与竖直向上方向的夹角，phi 是该射线在水平面内的方位角。两个角
 * 均以度为单位传入。
 *
 * 输出的 Green 函数和空间导数均为无量纲量。dG_source 按照
 * dG_source[k'][i][j] 的顺序存储源点坐标导数，dG_receiver 按照
 * dG_receiver[k][i][j] 的顺序存储接收点坐标导数。按照式 (7.4.1.4)，
 * 空间导数中的时间导数使用时间序列上的数值差分计算。若物理阶跃力 Green
 * 函数为 G^H，则 Gbar = pi^2 mu r G^H，且两类无量纲空间导数均为
 * pi^2 mu r^2 乘以相应的物理空间导数。
 *
 * @param[in]    nu       泊松比，(0, 0.5)
 * @param[in]    ts       归一化时间序列
 * @param[in]    nt       时间序列点数
 * @param[in]    theta    震源射线角，单位度，(0, 90)
 * @param[in]    azimuth  方位角，单位度，[0, 360]
 * @param[out]   G              无量纲阶跃力位移，G[time][i][j]
 * @param[out]   dG_source      无量纲源点导数，dG_source[time][k'][i][j]
 * @param[out]   dG_receiver    无量纲接收点导数，dG_receiver[time][k][i][j]
 *
 * 当三个输出指针同时为 NULL 时，结果输出到标准输出。
 */
void grt_solve_lamb2(
    const real_t nu, const real_t *ts, const int nt,
    const real_t theta, const real_t azimuth,
    real_t (*G)[3][3], real_t (*dG_source)[3][3][3],
    real_t (*dG_receiver)[3][3][3]);
