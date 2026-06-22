/**
 * @file   static_util.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-06
 * 
 * 一些关于静态解的辅助函数
 * 
 */

#pragma once

#include "grt/common/model.h"

/**
 * 根据核函数的振幅，估计静态解波数积分中合适的 kmax
 * 
 * @param[in]      mod1d        模型指针
 * @param[in]      kmax_ref     最大上限
 * @param[out]     Ncount       估计过程中计算静态核函数的次数
 * @return     估计的 kmax
 */
real_t grt_predict_static_kmax(MODEL1D *mod1d, real_t kmax_ref, size_t *Ncount);
