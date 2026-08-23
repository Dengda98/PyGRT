/**
 * @file   static_output.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * 静态位移模块的公共 NetCDF 输出
 *
 */

#pragma once

#include "grt/static/recv_points.h"

/**
 * 查询指定深度处的接收介质参数
 *
 * @param[in]   context   介质查询上下文
 * @param[in]   depth     接收点深度 (km)
 * @param[out]  va        P 波速度 (km/s)
 * @param[out]  vb        S 波速度 (km/s)
 * @param[out]  rho       密度 (g/cm^3)
 */
typedef void (*GRT_STATIC_MEDIUM_FUNC)(
    void *context, real_t depth, real_t *va, real_t *vb, real_t *rho);

/** 静态位移输出描述 */
typedef struct {
    const char *path;                            ///< 输出 NetCDF 文件路径
    const char *channels;                        ///< 输出分量编码
    const char *compute_type;                    ///< 震源类型属性
    const char *coordinate;                      ///< 坐标系属性
    bool calc_upar;                              ///< 是否包含位移偏导
    bool rot2ZNE;                                ///< 是否已旋转到 ZNE 分量
    bool has_depsrc;                             ///< 是否写入震源深度属性
    real_t depsrc;                               ///< 点源深度 (km)
    bool has_elastic_params;                     ///< 是否写入均匀介质弹性参数
    real_t alpha;                                ///< 均匀半空间 alpha 参数
    real_t lambda;                               ///< 均匀半空间 lambda 参数
    real_t mu;                                   ///< 均匀半空间 mu 参数
    const GRT_RECV_POINTS *recv;                 ///< 规则网格、任意点或有限接收断层点列表
    GRT_STATIC_MEDIUM_FUNC get_medium;           ///< 接收介质查询回调函数
    void *medium_context;                        ///< 接收介质查询回调上下文
    const real_t (*syn)[GRT_CHANNEL_NUM];        ///< 位移数组
    const real_t (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
                                                 ///< 位移偏导数组
} GRT_STATIC_NC_OUTPUT;

/**
 * 按公共 NetCDF 布局写出静态位移结果
 *
 * @param[in]  output   静态位移输出描述
 */
void grt_static_save_nc(const GRT_STATIC_NC_OUTPUT *output);
