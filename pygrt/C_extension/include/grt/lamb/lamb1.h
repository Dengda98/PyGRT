/**
 * @file   lamb1.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    使用广义闭合解求解第一类 Lamb 问题，参考：
 * 
 *        张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社
 */

#pragma once

#include "grt/common/const.h"


/**
 * 计算第一类 Lamb 问题的无量纲震相到时
 *
 * @param[in]    nu        泊松比，(0, 0.5)
 * @param[out]   tP        直达 P 波到时
 * @param[out]   tR        Rayleigh 波到时
 */
void grt_compute_lamb1_travt(const real_t nu, real_t *tP, real_t *tR);


/**
 * 使用广义闭合解求解第一类 Lamb 问题
 * 
 * @param[in]    nu           泊松比， (0, 0.5)
 * @param[in]    ts           无量纲时间序列 tbar=t/(r/beta)=beta*t/r
 * @param[in]    nt           时间序列点数
 * @param[in]    azimuth      初始震源到台站的方位角，单位度，[0, 360]
 * @param[in]    cbar         沿 x1 正方向的无量纲运动源速度 c/beta，0 表示固定源
 * @param[out]   u            记录结果的指针，如果为NULL则输出到标准输出
 *
 * cbar 非零时采用第 9 章的亚 Rayleigh 速度竖向运动点载荷闭合解
 * 台站不能位于 x1 轴上，且 cbar 必须小于 vR/beta
 * 这时仅 u[time][i][2] 为非零运动源结果
 */
void grt_solve_lamb1(
    const real_t nu, const real_t *ts, const int nt, const real_t azimuth, const real_t cbar, real_t (*u)[3][3]);
