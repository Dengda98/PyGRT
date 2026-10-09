/**
 * @file   static_nc.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * 静态解 NetCDF 接收布局读入和结果输出
 *
 */

#pragma once

#include "grt/common/rcv_point.h"

/** layout 字符串：写入 nc 全局属性，读端据此分支 */
#define GRT_RCV_LAYOUT_GRID   "grid"    ///< 静态 NetCDF 规则网格布局属性
#define GRT_RCV_LAYOUT_POINTS "points"  ///< 静态 NetCDF 任意接收点布局属性
#define GRT_RCV_LAYOUT_FAULTS "faults"  ///< 静态 NetCDF 有限接收断层布局属性

/** NetCDF 接收点布局 */
typedef enum {
    GRT_RCV_NC_LAYOUT_GRID = 0,  ///< 规则网格布局
    GRT_RCV_NC_LAYOUT_POINTS,    ///< 任意接收点布局
    GRT_RCV_NC_LAYOUT_FAULTS     ///< 有限接收断层布局
} GRT_RCV_NC_LAYOUT;

/** 静态 NetCDF 接收布局，写出时借用点列和断层，读入时按需分配坐标数组 */
typedef struct {
    GRT_RCV_NC_LAYOUT layout;    ///< 接收点布局类型
    size_t npts;                 ///< 展平后的接收点总数
    const RCV_POINT *rcvs;       ///< 写出时借用的点列，grid 沿 east 方向变化最快

    // grid 布局专用
    size_t nnorth;               ///< 网格 north 方向点数
    size_t neast;                ///< 网格 east 方向点数

    // faults 布局专用
    size_t nfault;               ///< 有限接收断层数量
    const FINITE_FAULT *faults;  ///< 写出时借用的断层数组，点列按断层及剖分索引排列
    int nfault_dimid;            ///< 读入文件中的 nfault 维度 ID，仅对该文件有效

    // 从 NetCDF 读入的维度信息，ID 仅对该文件有效
    int ndims;                   ///< 结果变量维度数，grid 为 2，points 和 faults 为 1
    int dimids[2];               ///< 结果变量维度 ID，grid 按 north/east 排列

    // 按需读入的展平坐标，由本结构体管理
    real_t *norths;              ///< 北向坐标数组 (km)，长度为 npts
    real_t *easts;               ///< 东向坐标数组 (km)，长度为 npts
} RCV_NC_INFO;

/**
 * 获取已打开 NetCDF 文件的接收布局
 *
 * @param[in]  ncid   已打开的 NetCDF 文件 ID
 * @return            接收布局类型
 */
GRT_RCV_NC_LAYOUT grt_rcv_nc_get_layout(int ncid);

/**
 * 读取已打开 NetCDF 文件中的接收布局和维度，不分配坐标数组
 * info 会被重新初始化，调用前不能持有尚未释放的坐标数组
 *
 * @param[in]   ncid   已打开的 NetCDF 文件 ID
 * @param[out]  info   接收布局和维度信息
 */
void grt_rcv_nc_info_load(int ncid, RCV_NC_INFO *info);

/**
 * 按布局读取并展平接收坐标，调用前需已读取同一文件的布局和维度
 *
 * @param[in]      ncid   已打开的 NetCDF 文件 ID
 * @param[in,out]  info   接收布局及坐标，已有坐标数组会被替换
 */
void grt_rcv_nc_info_load_coordinates(int ncid, RCV_NC_INFO *info);

/**
 * 释放自行分配的坐标数组并清空接收信息，不释放借用的点列和断层
 *
 * @param[in,out]  info   接收布局及坐标
 */
void grt_rcv_nc_info_free(RCV_NC_INFO *info);

/**
 * 向处于定义模式的文件写入公共接收布局、介质和静态位移结果
 * 调用方负责创建文件、写入模块专属属性及关闭文件，本函数结束定义模式后写入数据
 * grid 的统一接收深度和介质保存为全局属性，points 和 faults 的深度和介质保存为逐点变量
 *
 * @param[in]  ncid       已创建且处于定义模式的 NetCDF 文件 ID
 * @param[in]  receivers  接收布局
 * @param[in]  nlayer     原始模型层数，均匀半空间为 1
 * @param[in]  modarr     原始模型矩阵，每行 Thk/Va/Vb/Rho/Qa/Qb
 * @param[in]  rot2ZNE    是否使用 ZNE 分量，否则为 ZRT
 * @param[in]  calc_upar  是否记录位移偏导
 * @param[in]  syn        位移数组，按接收点及分量排列
 * @param[in]  syn_upar   位移偏导数组，按接收点、求导方向及分量排列，不记录时可为 NULL
 */
void grt_static_nc_write(int ncid, const RCV_NC_INFO *receivers, size_t nlayer, const real_t (*modarr)[GRT_MODARR_NCOL],
                         bool rot2ZNE, bool calc_upar, const real_t (*syn)[GRT_CHANNEL_NUM],
                         const real_t (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM]);
