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

/**
 * Coulomb 程序格式的有限断层（及衍生量）
 *
 * 原始断层数据区的 11 列依次为
 * ID、X-start、Y-start、X-fin、Y-fin、Kode、slip1、slip2、dip、top、bot
 * 其中 X/Y 在 PyGRT 中分别对应 east/north
 *
 * .inp 文件的 slip1/slip2 列含义为右旋/逆冲或 Kode 指定的物理量
 * .inr 文件仅允许 Kode=100，slip1/slip2 列含义为 rake (degree)/net slip (m)
 *
 * Kode=100：矩形断层，slip1=右旋滑动 (m)，slip2=逆冲滑动 (m)
 * Kode=200：矩形断层，slip1=右旋滑动 (m)，slip2=张裂开度 (m)
 * Kode=300：矩形断层，slip1=张裂开度 (m)，slip2=逆冲滑动 (m)
 * Kode=400：点双力偶，slip1=走向滑动 potency (m^3)，slip2=倾向滑动 potency (m^3)
 * Kode=500：点张裂/膨胀源，slip1=张裂 potency (m^3)，slip2=膨胀 potency (m^3)
 *
 * value1/value2 保留文件第 7、8 列的值
 * right_lateral/reverse/tensile/inflate 是按 Kode 解释后的量
 */
typedef struct {
    real_t east_begin;
    real_t north_begin;
    real_t east_end;
    real_t north_end;

    unsigned int kode;      ///< Coulomb 标识符，只允许 100、200、300、400、500
    bool rake_format;       ///< 是否按 rake/net slip 格式解释
    real_t value1;          ///< 文件第 7 列，物理意义由 kode 决定
    real_t value2;          ///< 文件第 8 列，物理意义由 kode 决定

    real_t right_lateral;   ///< Kode 100/200 为 m，Kode 400 为 m^3
    real_t reverse;         ///< Kode 100/300 为 m，Kode 400 为 m^3
    real_t tensile;         ///< Kode 200/300 为 m，Kode 500 为 m^3
    real_t inflate;         ///< Kode 500 为 m^3
    real_t dip;             ///< degree
    real_t top;             ///< km
    real_t bot;             ///< km

    // 衍生量，由 grt_finite_fault_set_derived 填充
    real_t strike;          ///< degree
    real_t rake;            ///< degree
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
 * 由 Coulomb 原始字段计算 strike / rake / slip
 *
 * @param[in,out]  f   有限断层结构体
 */
void grt_finite_fault_set_derived(FINITE_FAULT *f);

/**
 * 读取 Coulomb 格式有限断层文件
 *
 * 跳过前两行表头，逐行读取断层数据，并对每条调用 grt_finite_fault_set_derived
 * 支持 .inp 与 .inr，要求 dip ∈ (0, 90]、bot > top 且沿走向长度 > 0
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
