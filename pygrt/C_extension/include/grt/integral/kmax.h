/**
 * @file   kmax.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-06
 * 
 * 波数积分辅助函数
 * 
 */

#pragma once

#include "grt/integral/kernel.h"

/**
 * 根据核函数的振幅，估计波数积分中合适的 kmax。
 *
 * 调用方须先按对应的静态或动态频率更新 mstat。对于静态解，通常从
 * kmax_ref * 1e-3 开始；对于动态解，通常从 omega / vmin 开始。
 * 
 * @param[in,out]  mstat        已设置频率的模型状态
 * @param[in]      kerfunc      待检查的核函数
 * @param[in]      kmax_init    扫描的初始波数
 * @param[in]      kmax_ref     最大上限
 * @param[out]     Ncount       估计过程中计算核函数的次数，可为 NULL
 * @return     估计的 kmax
 */
real_t grt_predict_kmax(
    MODEL1D_STATE *mstat, GRT_KernelFunc kerfunc,
    real_t kmax_init, real_t kmax_ref, size_t *Ncount);
