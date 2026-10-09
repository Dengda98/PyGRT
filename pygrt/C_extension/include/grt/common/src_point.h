/**
 * @file   src_point.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-10-03
 *
 *    动态解和静态解共用的震源点展开与基本源型分量判断
 */

#pragma once

#include "grt/common/finite_fault.h"

/** 单个震源点，断层由调用方持有，须比点数组存活更久 */
typedef struct {
    real_t north;              ///< 北向坐标，km
    real_t east;               ///< 东向坐标，km
    real_t depth;              ///< 深度，km
    const FINITE_FAULT *fault;  ///< 所属断层，借用指针
    size_t isub;               ///< 剖分索引，沿走向变化最快
} SRC_POINT;

/**
 * 将已剖分的断层展开为源点数组，点源 Kode 和通用单点源只产生一个中心点
 * @param[in]   nfault  断层数量
 * @param[in]   faults  已调用 grt_finite_fault_subdiv 的断层数组
 * @param[out]  npts    源点数量
 * @return      新分配的源点数组，调用方用 free 释放
 */
SRC_POINT *grt_src_points_from_faults(size_t nfault, const FINITE_FAULT *faults, size_t *npts);

