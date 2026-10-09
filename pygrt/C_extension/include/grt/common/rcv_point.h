/**
 * @file   rcv_point.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * 动态解和静态解共用的接收点列表
 * 网格 (-X/-Y 或延用库坐标) 与任意点文件 (-Q) 均展开为点列
 *
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "grt/common/const.h"
#include "grt/common/finite_fault.h"

/** 单个接收点，有限断层由调用方持有，须比点数组存活更久 */
typedef struct {
    real_t north;               ///< 北向坐标，km
    real_t east;                ///< 东向坐标，km
    real_t depth;               ///< 深度，km

    bool has_mechanism;         ///< 是否包含逐点机制，有限断层机制从 fault 获取
    real_t strike;              ///< 逐点走向，degree
    real_t dip;                 ///< 逐点倾角，degree
    real_t rake;                ///< 逐点滑动角，degree

    const FINITE_FAULT *fault;  ///< 所属有限断层，普通接收点为 NULL，借用指针
    size_t isub;                ///< 剖分索引，沿走向变化最快
} RCV_POINT;

/**
 * 由极坐标创建单个接收点，调用方用 free 释放
 * @param[in]  dist     震中距，km
 * @param[in]  azimuth  方位角，degree
 * @param[in]  depth    接收深度，km
 */
RCV_POINT *grt_rcv_points_from_polar(real_t dist, real_t azimuth, real_t depth);

/**
 * 将规则网格展开为接收点数组，沿 east 方向变化最快，调用方用 free 释放
 * @param[in]  nnorth  north 方向点数
 * @param[in]  norths  north 坐标，km
 * @param[in]  neast   east 方向点数
 * @param[in]  easts   east 坐标，km
 * @param[in]  depth   接收深度，km
 */
RCV_POINT *grt_rcv_points_from_grid(size_t nnorth, const real_t *norths, size_t neast, const real_t *easts, real_t depth);

/**
 * 从 ASCII 文件读取接收点，每行 north east depth [strike dip rake]，调用方用 free 释放
 * @param[in]   path  文件路径
 * @param[out]  npts  接收点数量
 */
RCV_POINT *grt_rcv_points_from_file(const char *path, size_t *npts);

/**
 * 将已剖分的有限接收断层展开为中心点数组，调用方用 free 释放
 * @param[in]   nfault  断层数量
 * @param[in]   faults  已调用 grt_finite_fault_subdiv 的断层数组
 * @param[out]  npts    接收点数量
 */
RCV_POINT *grt_rcv_points_from_faults(size_t nfault, const FINITE_FAULT *faults, size_t *npts);
