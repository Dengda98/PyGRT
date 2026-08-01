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
 * 调用方须先按对应的静态或动态频率更新 mstat，通常从 dk 开始
 *
 * 搜索在 log(k) 上等间距推进（自适应几何因子），目标约 40 步覆盖 [kmax_init, kmax_ref]
 * 在 kmax_low 之前不判断收敛
 * 同深度使用当前核函数与前一采样点的差值 dF 相对于搜索过程中的最大振幅 Fmax 判断逼近常数；
 * 异深度用振幅相对峰值衰减判断逼近 0
 * 
 * @param[in,out]  mstat        已设置频率的模型状态
 * @param[in]      kerfunc      待检查的核函数
 * @param[in]      kmax_init    扫描的初始波数
 * @param[in]      kmax_low     允许判断收敛的最低波数
 * @param[in]      kmax_ref     最大上限
 * @param[out]     Ncount       估计过程中计算核函数的次数，可为 NULL
 * @return     估计的 kmax
 */
real_t grt_predict_kmax(
    MODEL1D_STATE *mstat, GRT_KernelFunc kerfunc,
    real_t kmax_init, real_t kmax_low, real_t kmax_ref, size_t *Ncount);
