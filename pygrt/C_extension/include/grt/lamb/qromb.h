/**
 * @file   qromb.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-09
 *
 *    Romberg numerical integration used by the Lamb-problem fallbacks
 */

#pragma once

#include "grt/common/const.h"

typedef real_t (*GRT_LAMB_INTEGRAND)(real_t x, void *userdata);

/**
 * 用梯形积分逐次分半并进行 Richardson 外推，计算有限区间积分
 *
 * @param[in]    func       被积函数
 * @param[in]    a          积分下限
 * @param[in]    b          积分上限
 * @param[in]    eps        相对收敛容限
 * @param[in]    userdata   传递给被积函数的用户数据
 *
 * @return       积分结果
 */
real_t grt_lamb_qromb(
    GRT_LAMB_INTEGRAND func, real_t a, real_t b, real_t eps, void *userdata);
