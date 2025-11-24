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
#include "grt/lamb/lamb_util.h"


/**
 * 使用广义闭合解求解第一类 Lamb 问题
 * 
 * @param[in]    nu        泊松比， (0, 0.5)
 * @param[in]    ts        归一化时间序列
 * @param[in]    nt        时间序列点数
 * @param[in]    azimuth   方位角，单位度
 * @param[out]   u         记录结果的指针，如果为NULL则输出到标准输出
 * 
 */
void grt_solve_lamb1(
    const real_t nu, const real_t *ts, const int nt, const real_t azimuth, real_t (*u)[3][3]);