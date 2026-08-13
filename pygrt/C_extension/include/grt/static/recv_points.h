/**
 * @file   recv_points.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * 静态 syn / 后处理用的接收点列表
 * 网格 (-X/-Y 或延用库坐标) 与任意点文件 (-Q) 均展开为点列
 *
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "grt/common/const.h"

/** layout 字符串：写入 nc 全局属性，读端据此分支 */
#define GRT_RECV_LAYOUT_GRID   "grid"
#define GRT_RECV_LAYOUT_POINTS "points"

/**
 * 接收点坐标列表（三个分量各自连续存放）
 *
 * norths/easts/depths 长度均为 npts，单位 km
 * is_grid 为真时可由 nnorth/neast 还原二维网格，ipt = ieast + inorth*neast
 */
typedef struct {
    size_t npts;
    real_t *norths;
    real_t *easts;
    real_t *depths;

    bool is_grid;
    size_t nnorth;
    size_t neast;
} GRT_RECV_POINTS;

/**
 * 由 north/east 轴与单一深度展开为点列（is_grid=true）
 *
 * @param[in]   nnorth   north 方向点数
 * @param[in]   norths   north 坐标 (km)
 * @param[in]   neast    east 方向点数
 * @param[in]   easts    east 坐标 (km)
 * @param[in]   depth    接收深度 (km)
 * @return      新分配的 GRT_RECV_POINTS*，调用方负责 grt_recv_points_free
 */
GRT_RECV_POINTS *grt_recv_points_from_grid(
    size_t nnorth, const real_t *norths,
    size_t neast,  const real_t *easts,
    real_t depth);

/**
 * 从 ASCII 文件读任意接收点（is_grid=false）
 *
 * 每行 north east depth (km)，# 开头为注释
 *
 * @param[in]   path   文件路径
 * @return      新分配的 GRT_RECV_POINTS*
 */
GRT_RECV_POINTS *grt_recv_points_from_file(const char *path);

/**
 * 释放 GRT_RECV_POINTS（含 norths/easts/depths）
 *
 * @param[in,out]  pts   可为 NULL
 */
void grt_recv_points_free(GRT_RECV_POINTS *pts);

/**
 * 从已打开的 nc 判断是否为 points 布局
 *
 * 优先读全局属性 layout；若无属性则看是否存在 point 维
 *
 * @param[in]   ncid   已打开的 nc id
 * @return      true 表示 points，false 表示 grid
 */
bool grt_recv_nc_is_points(int ncid);
