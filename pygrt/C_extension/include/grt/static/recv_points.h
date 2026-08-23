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
#include "grt/common/finite_fault.h"

/** layout 字符串：写入 nc 全局属性，读端据此分支 */
#define GRT_RECV_LAYOUT_GRID   "grid"
#define GRT_RECV_LAYOUT_POINTS "points"

/** NetCDF 接收点布局 */
typedef enum {
    GRT_RECV_NC_LAYOUT_GRID = 0,           ///< 规则网格布局
    GRT_RECV_NC_LAYOUT_POINTS              ///< 一维接收点布局
} GRT_RECV_NC_LAYOUT;

/** 接收点坐标及可选的有限接收断层信息 */
typedef struct {
    size_t npts;            ///< 接收点总数
    real_t *norths;         ///< 北向坐标数组 (km)
    real_t *easts;          ///< 东向坐标数组 (km)
    real_t *depths;         ///< 深度坐标数组 (km)

    bool is_grid;           ///< 是否由二维规则网格展开
    size_t nnorth;          ///< north 方向点数
    size_t neast;           ///< east 方向点数

    bool has_geometry;      ///< 任意接收点是否包含逐点几何
    real_t *strikes;        ///< 任意接收点的走向 (degree)
    real_t *dips;           ///< 任意接收点的倾角 (degree)
    real_t *rakes;          ///< 任意接收点的滑动角 (degree)

    bool is_fault;          ///< 是否由有限接收断层剖分得到
    size_t nfault;          ///< 有限接收断层数量
    size_t *nsubs;          ///< 每条有限接收断层的子断层数量
    size_t *offsets;        ///< 每条有限接收断层在 point 数组中的排他性结束索引
    real_t *fstrikes;       ///< 每条有限接收断层的走向 (degree)
    real_t *fdips;          ///< 每条有限接收断层的倾角 (degree)
    real_t *frakes;         ///< 每条有限接收断层的滑动角 (degree)
    size_t *stksizes;       ///< 每条有限接收断层沿走向的子断层数量
    size_t *dipsizes;       ///< 每条有限接收断层沿倾向的子断层数量
} GRT_RECV_POINTS;

/** 已打开 NetCDF 文件中的接收坐标及维度信息 */
typedef struct {
    GRT_RECV_NC_LAYOUT layout;     ///< 接收点布局类型
    size_t npts;                   ///< 展平后的接收点总数
    real_t *norths;                ///< 展平后的北向坐标数组 (km)
    real_t *easts;                 ///< 展平后的东向坐标数组 (km)

    int dimids[2];                 ///< 接收坐标对应的 NetCDF 维度 ID
} GRT_RECV_NC_INFO;

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
 * 每行可以是 north east depth (km)，也可以在其后增加
 * strike dip rake (degree)，# 开头为注释
 *
 * @param[in]   path   文件路径
 * @return      新分配的 GRT_RECV_POINTS*
 */
GRT_RECV_POINTS *grt_recv_points_from_file(const char *path);

/**
 * 释放 GRT_RECV_POINTS（含坐标和可选接收断层几何）
 *
 * @param[in,out]  pts   可为 NULL
 */
void grt_recv_points_free(GRT_RECV_POINTS *pts);

/**
 * 从有限断层数组生成接收点列表
 *
 * dL/dW 均为非正值时不剖分，每条有限断层只生成一个中心点
 * dL/dW 均为正值时按沿走向/沿倾向尺寸剖分
 * 生成的点数组同时保存每条有限断层的索引和几何信息
 *
 * @param[in]   nfault   有限断层数量
 * @param[in]   faults   已读取并建立衍生量的有限断层数组
 * @param[in]   dL       沿走向子断层尺寸 (km)
 * @param[in]   dW       沿倾向子断层尺寸 (km)
 * @return      新分配的 GRT_RECV_POINTS*，调用方负责 grt_recv_points_free
 */
GRT_RECV_POINTS *grt_recv_points_from_faults(
    size_t nfault, const FINITE_FAULT *faults, real_t dL, real_t dW);

/**
 * 获取已打开 NetCDF 文件的接收布局
 *
 * @param[in]  ncid   已打开的 NetCDF 文件 ID
 * @return            接收布局类型
 */
GRT_RECV_NC_LAYOUT grt_recv_nc_get_layout(int ncid);

/**
 * 读取已打开 NetCDF 文件中的接收坐标布局
 *
 * @param[in]   ncid   已打开的 NetCDF 文件 ID
 * @param[out]  info   接收坐标和维度信息
 */
void grt_recv_nc_info_load(int ncid, GRT_RECV_NC_INFO *info);

/**
 * 释放 GRT_RECV_NC_INFO
 *
 * @param[in,out]  info   接收坐标和维度信息
 */
void grt_recv_nc_info_free(GRT_RECV_NC_INFO *info);

/**
 * 从已打开的 nc 判断是否为 points 布局
 *
 * 读取全局属性 layout，输入文件必须显式保存该属性
 *
 * @param[in]   ncid   已打开的 nc id
 * @return      true 表示 points，false 表示 grid
 */
bool grt_recv_nc_is_points(int ncid);
