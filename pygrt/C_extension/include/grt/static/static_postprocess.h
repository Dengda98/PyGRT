/**
 * @file   static_postprocess.h
 * @brief  静态位移偏导后处理（应力/应变/旋转）
 * @date   2026-07
 */

#pragma once

#include "grt/common/const.h"

/**
 * 由静态位移偏导合成对称应力张量。
 *
 * 数组布局：u[分量][点]、upar[偏导方向][分量][点]、
 * res[第二分量][第一分量][点]。仅写入 res 的上三角分量。
 *
 * ZRT 联络项：r≠0 用 u/r；r=0 改用 ∂_r u（upar[1][*]），
 * 与 syn 中轴点处 (1/r)∂_θ 有限部分配套。
 */
void grt_static_compute_stress(
    size_t nnorth, size_t neast, const real_t *norths, const real_t *easts,
    real_t *const u[GRT_CHANNEL_NUM],
    real_t *const upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    real_t *const res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    bool rot2ZNE, real_t mu, real_t lam);

/**
 * 由静态位移偏导合成对称应变张量。
 * 数组布局同 grt_static_compute_stress()。
 */
void grt_static_compute_strain(
    size_t nnorth, size_t neast, const real_t *norths, const real_t *easts,
    real_t *const u[GRT_CHANNEL_NUM],
    real_t *const upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    real_t *const res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM], bool rot2ZNE);

/**
 * 由静态位移偏导合成反对称旋转张量。
 * 数组布局同 grt_static_compute_stress()。仅写入 res 的非对角分量。
 */
void grt_static_compute_rotation(
    size_t nnorth, size_t neast, const real_t *norths, const real_t *easts,
    real_t *const u[GRT_CHANNEL_NUM],
    real_t *const upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    real_t *const res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM], bool rot2ZNE);
