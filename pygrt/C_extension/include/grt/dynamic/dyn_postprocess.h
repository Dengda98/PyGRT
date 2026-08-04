/**
 * @file   dyn_postprocess.h
 * @brief  动态位移偏导后处理（应力/应变/旋转）
 * @date   2026-07
 */

#pragma once

#include "grt/common/const.h"

/**
 * 由动态位移偏导合成应变张量。
 * 数组布局：u[分量][采样点]、upar[偏导方向][分量][采样点]、
 * res[第二分量][第一分量][采样点]。
 */
void grt_compute_strain(
    size_t npts, float dist, float *const u[GRT_CHANNEL_NUM],
    float *const upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    float *const res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM], bool rot2ZNE);

/** 由动态位移偏导合成旋转张量。数组布局同 grt_compute_strain()。 */
void grt_compute_rotation(
    size_t npts, float dist, float *const u[GRT_CHANNEL_NUM],
    float *const upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    float *const res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM], bool rot2ZNE);

/**
 * 在频域由动态位移偏导合成应力张量。
 * 介质参数约定与命令行 stress 子模块的 SAC 头段一致。
 */
void grt_compute_stress(
    size_t npts, float dt, float dist, float va, float vb, float rho,
    float Qainv, float Qbinv, float *const u[GRT_CHANNEL_NUM],
    float *const upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    float *const res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM], bool rot2ZNE);
