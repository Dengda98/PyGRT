/**
 * @file   finite_fault.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * Coulomb 格式有限断层：读入、衍生量与几何剖分
 *
 */

#pragma once

#include "grt/common/const.h"

/** Coulomb 程序格式的有限断层（及衍生量） */
typedef struct {
    real_t east_begin;
    real_t north_begin;
    real_t east_end;
    real_t north_end;
    real_t right_lateral;   ///< m
    real_t reverse;         ///< m
    real_t dip;             ///< degree
    real_t top;             ///< km
    real_t bot;             ///< km

    // 衍生量，由 grt_finite_fault_set_derived 填充
    real_t strike;          ///< degree
    real_t rake;            ///< degree
    real_t slip;            ///< m, hypot(right_lateral, reverse)
} FINITE_FAULT;

/** 剖分后的单个子断层（点源几何） */
typedef struct {
    real_t east;      ///< 中心 east (km)
    real_t north;     ///< 中心 north (km)
    real_t depsrc;    ///< 中心深度 (km)
    real_t width;     ///< 沿倾向边长 (km)
    real_t length;    ///< 沿走向边长 (km)
    real_t potency;   ///< 矩势 (cm^3) = slip(m)*width*length*1e12
} FINITE_SUBFAULT;

/**
 * 由 Coulomb 原始字段计算 strike / rake / slip
 *
 * @param[in,out]  f   有限断层结构体
 */
void grt_finite_fault_set_derived(FINITE_FAULT *f);

/**
 * 读取 Coulomb 格式有限断层文件
 *
 * 跳过前两行表头，逐行解析 11 列中的断层字段，并对每条调用
 * grt_finite_fault_set_derived。要求 dip ∈ (0, 90] 且 bot > top
 * 调用方负责 grt_finite_fault_free
 *
 * @param[in]   path     文件路径
 * @param[out]  nfault   读入的断层段数
 * @return      新分配的 FINITE_FAULT 数组，失败则报错退出
 */
FINITE_FAULT *grt_finite_fault_load_coulomb(const char *path, size_t *nfault);

/**
 * 释放 grt_finite_fault_load_coulomb 返回的数组
 *
 * @param[in,out]  faults   断层数组，可为 NULL
 */
void grt_finite_fault_free(FINITE_FAULT *faults);

/**
 * 由断层几何与 dL/dW 得到倾向/走向长度及剖分数
 *
 * 要求 dip ∈ (0, 90]、bot > top，且 dL/dW > 0
 *
 * @param[in]   fault  有限断层（需已有 dip / 端点坐标）
 * @param[in]   dL     沿走向剖分间隔 (km)
 * @param[in]   dW     沿倾向剖分间隔 (km)
 * @param[out]  W      沿倾向长度 (km)
 * @param[out]  L      沿走向长度 (km)
 * @param[out]  nW     倾向方向子断层数
 * @param[out]  nL     走向方向子断层数
 */
void grt_finite_fault_subdiv(
    const FINITE_FAULT *fault, real_t dL, real_t dW,
    real_t *W, real_t *L, size_t *nW, size_t *nL);

/**
 * 计算 (iW, iL) 子断层中心几何与矩势
 *
 * 要求 fault 已 set_derived，且 W/L/nW/nL 由 grt_finite_fault_subdiv 给出
 * 末块可短于 dL/dW，中心取该块中点：i*d + size/2
 *
 * @param[in]   fault  有限断层
 * @param[in]   dL     沿走向剖分间隔 (km)
 * @param[in]   dW     沿倾向剖分间隔 (km)
 * @param[in]   W      沿倾向总长 (km)
 * @param[in]   L      沿走向总长 (km)
 * @param[in]   iW     倾向方向子断层索引 [0, nW)
 * @param[in]   iL     走向方向子断层索引 [0, nL)
 * @param[out]  sub    子断层几何
 */
void grt_finite_fault_subfault(
    const FINITE_FAULT *fault,
    real_t dL, real_t dW, real_t W, real_t L,
    size_t iW, size_t iL,
    FINITE_SUBFAULT *sub);
