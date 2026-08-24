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

#define KODE_RTLAT_REVERSE          100  ///< 矩形右旋/逆冲断层
#define KODE_RTLAT_TENSILE          200  ///< 矩形右旋/张裂断层
#define KODE_TENSILE_REVERSE        300  ///< 矩形张裂/逆冲断层
#define KODE_POINT_DC               400  ///< 点双力偶源
#define KODE_POINT_TENSILE_INFLATE  500  ///< 点张裂/膨胀源
#define KODE_IS_FINITE(kode)        ((kode) == KODE_RTLAT_REVERSE || (kode) == KODE_RTLAT_TENSILE || (kode) == KODE_TENSILE_REVERSE)
#define KODE_IS_POINT(kode)         ((kode) == KODE_POINT_DC || (kode) == KODE_POINT_TENSILE_INFLATE)

/** 滑动方向未定义时使用的 rake 值，超出合法角度范围 */
#define GRT_FINITE_FAULT_UNDEFINED_RAKE (-999.0)

/**
 * Coulomb 程序格式的有限断层（及衍生量）
 *
 * Coulomb 表格的断层数据区严格包含 11 个数值列，依次为
 * ID、X-start、Y-start、X-fin、Y-fin、Kode、value1、value2、dip、top、bot
 * 其中 X/Y 在 PyGRT 中分别对应 east/north
 *
 * 文件前两行是 Coulomb 表头：第一行首个 token 必须为 #，其后给出 X-start、Y-start、X-fin、Y-fin、
 * Kode、value1、value2、dip、top、bot 共 10 个字段标签，第二行给出 11 个占位字段
 * Coulomb 常见表头中的 "dip angle" 可写成两个空白分隔的 token，但数据区仍为 11 列
 * 第一行第 7 个数据列的标签精确为 "rake" 时，value1/value2 按 rake (degree)/net slip (m) 解释
 * 文件名后缀不参与格式选择，仅用于发现与表头标识不一致的情况
 *
 * Kode=100：矩形断层，value1=右旋滑动 (m)，value2=逆冲滑动 (m)
 * Kode=200：矩形断层，value1=右旋滑动 (m)，value2=张裂开度 (m)
 * Kode=300：矩形断层，value1=张裂开度 (m)，value2=逆冲滑动 (m)
 * Kode=400：点双力偶，value1=右旋滑动 potency (m^3)，value2=逆冲滑动 potency (m^3)
 * Kode=500：点张裂/膨胀源，value1=张裂 potency (m^3)，value2=膨胀 potency (m^3)
 *
 * value1/value2 保留文件第 7、8 列的值
 * right_lateral/reverse/tensile/inflate 是按 Kode 解释后的量
 */
typedef struct {
    real_t east_begin;       ///< 上边界起点 east 坐标 (km)
    real_t north_begin;      ///< 上边界起点 north 坐标 (km)
    real_t east_end;         ///< 上边界终点 east 坐标 (km)
    real_t north_end;        ///< 上边界终点 north 坐标 (km)

    unsigned int kode;      ///< Coulomb 标识符，只允许 100、200、300、400、500
    bool rake_format;       ///< 是否按表头标识的 rake/net slip 格式解释
    real_t value1;          ///< 文件第 7 列，物理意义由表头标识与 kode 共同决定
    real_t value2;          ///< 文件第 8 列，物理意义由表头标识与 kode 共同决定

    real_t right_lateral;   ///< Kode 100/200 为 m，Kode 400 为 m^3
    real_t reverse;         ///< Kode 100/300 为 m，Kode 400 为 m^3
    real_t tensile;         ///< Kode 200/300 为 m，Kode 500 为 m^3
    real_t inflate;         ///< Kode 500 为 m^3
    real_t dip;             ///< degree
    real_t top;             ///< km
    real_t bot;             ///< km

    // 衍生量，由有限断层读取函数直接填充
    real_t strike;          ///< 走向 (degree)
    real_t rake;            ///< 滑动角 (degree)，未定义时为 GRT_FINITE_FAULT_UNDEFINED_RAKE
    real_t slip;            ///< Kode 100/200/300 的剪切滑动合量，单位 m
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
 * 读取 Coulomb 格式有限断层文件
 *
 * 第一行首个 token 为 #，其后给出 10 个字段标签；第二行为 11 个占位字段
 * 读取前两行 Coulomb 表头，逐行读取 11 列断层数据并建立衍生量
 * 第 7 个数据列的标签精确为 rake 时按 rake/net slip 格式解释
 * 要求 dip ∈ (0, 90]、bot > top 且沿走向长度 > 0
 * 调用方负责 grt_finite_fault_free
 *
 * @param[in]   path     文件路径
 * @param[out]  nfault   读入的断层段数
 * @return      新分配的 FINITE_FAULT 数组，失败则报错退出
 */
FINITE_FAULT *grt_finite_fault_load_coulomb(const char *path, size_t *nfault);

/**
 * 解析有限断层选项并读取 Coulomb 格式有限断层文件
 *
 * 选项格式为 <fault>[+i<dL>/<dW>]，不提供 +i 时将 dL/dW 置为非正值
 * 返回的断层数组已经建立衍生量，调用方负责 grt_finite_fault_free
 *
 * @param[in]   option   有限断层选项值，不含 -C 或 -R 选项字符
 * @param[out]  nfault   读入的断层段数
 * @param[out]  dL       沿走向剖分尺寸 (km)，未指定时为非正值
 * @param[out]  dW       沿倾向剖分尺寸 (km)，未指定时为非正值
 * @return      新分配的 FINITE_FAULT 数组
 */
FINITE_FAULT *grt_finite_fault_from_option(const char *option, size_t *nfault, real_t *dL, real_t *dW);

/**
 * 释放 grt_finite_fault_load_coulomb 返回的数组
 *
 * @param[in,out]  faults   断层数组，可为 NULL
 */
void grt_finite_fault_free(FINITE_FAULT *faults);

/**
 * 由断层几何与 dL/dW 得到倾向/走向长度及剖分数
 *
 * 要求 dip ∈ (0, 90]、bot > top
 * dL/dW 均为非正值时不剖分，返回整个有限断层及 1×1 的剖分数量
 * dL/dW 必须同时为正值或同时为非正值
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
 * 要求 fault 已由有限断层读取函数建立衍生量，且 W/L/nW/nL
 * 由 grt_finite_fault_subdiv 给出
 * dL/dW 均为非正值时使用 W/L 作为单个子断层尺寸
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
