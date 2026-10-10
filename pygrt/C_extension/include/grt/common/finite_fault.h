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

#define KODE_POINT_GENERIC           0   ///< 通用单点源，仅用于内部构造，不属于 Coulomb 输入
#define KODE_RTLAT_REVERSE          100  ///< 矩形右旋/逆冲断层
#define KODE_RTLAT_TENSILE          200  ///< 矩形右旋/张裂断层
#define KODE_TENSILE_REVERSE        300  ///< 矩形张裂/逆冲断层
#define KODE_POINT_DC               400  ///< 点双力偶源
#define KODE_POINT_TENSILE_INFLATE  500  ///< 点张裂/膨胀源
#define KODE_IS_FINITE(kode)        ((kode) == KODE_RTLAT_REVERSE || (kode) == KODE_RTLAT_TENSILE || (kode) == KODE_TENSILE_REVERSE)  ///< 矩形断层
#define KODE_IS_POINT(kode)         ((kode) == KODE_POINT_GENERIC || (kode) == KODE_POINT_DC || (kode) == KODE_POINT_TENSILE_INFLATE)  ///< 单点源

#define GRT_FINITE_FAULT_UNDEFINED_RAKE (-999.0)  ///< 滑动方向未定义的标记

/** 断层的一个基本源型分量，矩形断层强度按单位面积保存 */
typedef struct {
    GRT_SYN_TYPE type;                    ///< 基本源型
    real_t scale;                         ///< 点源为地震矩、力或矩势，矩形断层为单位面积矩势 (cm^3/km^2)
    bool with_mu;                         ///< 是否乘源点剪切模量
    real_t mechanism[GRT_MECHANISM_NUM];  ///< 机制参数
} FINITE_SOURCE_TERM;

/**
 * Coulomb 程序格式的有限断层（及衍生量）
 *
 * Coulomb 表格的断层数据区严格包含 11 个数值列，依次为
 * 占位符、X-start、Y-start、X-fin、Y-fin、Kode、value1、value2、dip、top、bot
 * 第一列为占位符，约定全部填 1，不参与计算
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
 * 接收断层仅支持 Kode=100，滑动角由 rake 格式或右旋、逆冲分量格式确定
 */
typedef struct {
    real_t east_begin;               ///< 上边界起点 east 坐标 (km)
    real_t north_begin;              ///< 上边界起点 north 坐标 (km)
    real_t east_end;                 ///< 上边界终点 east 坐标 (km)
    real_t north_end;                ///< 上边界终点 north 坐标 (km)

    unsigned int kode;               ///< 震源为 Coulomb Kode 或通用单点源 0，接收断层仅支持 100
    real_t value1;                   ///< 文件第 7 列，物理意义由表头标识与 kode 共同决定
    real_t value2;                   ///< 文件第 8 列，物理意义由表头标识与 kode 共同决定

    real_t right_lateral;            ///< Kode 100/200 为 m，Kode 400 为 m^3
    real_t reverse;                  ///< Kode 100/300 为 m，Kode 400 为 m^3
    real_t tensile;                  ///< Kode 200/300 为 m，Kode 500 为 m^3
    real_t inflate;                  ///< Kode 500 为 m^3
    real_t dip;                      ///< degree
    real_t top;                      ///< km
    real_t bot;                      ///< km

    // 衍生量，由有限断层读取函数直接填充
    real_t strike;                   ///< 走向 (degree)
    real_t rake;                     ///< 滑动角 (degree)，未定义时为 GRT_FINITE_FAULT_UNDEFINED_RAKE

    // 剖分结果，由模块确定尺寸后一次性填充，索引 iW*nL+iL 沿走向变化最快
    size_t nL;                       ///< 沿走向剖分数
    size_t nW;                       ///< 沿倾向剖分数
    real_t *east;                    ///< 子断层中心 east 坐标 (km)，长度为 nW*nL
    real_t *north;                   ///< 子断层中心 north 坐标 (km)，长度为 nW*nL
    real_t *depth;                   ///< 子断层中心深度 (km)，长度为 nW*nL
    real_t *width;                   ///< 子断层沿倾向边长 (km)，长度为 nW*nL
    real_t *length;                  ///< 子断层沿走向边长 (km)，长度为 nW*nL
    real_t *scalar_moments;          ///< 子断层完整等效张量的标量矩 (dyne-cm)，长度为 nW*nL，仅动态有限源保存矩率时计算

    // 断层共享的震源描述，STF 在动态震源读入时生成
    bool stf_explicit;                          ///< 是否显式指定行内或全局时间函数
    real_t *stfd;                               ///< 无延迟、已阻尼的源时间函数，静态解及接收断层为 NULL
    int stf_npts;                               ///< 源时间函数样本数
    real_t stf_delay;                           ///< 采样对齐后的显式破裂延迟，s
    FINITE_SOURCE_TERM terms[2];                ///< 基本源型分量，同一断层的子源共用
    int nterms;                                 ///< 非零基本源型分量数，0、1 或 2
    bool required_gf_sources[GRT_SRC_M_NUM];    ///< 震源所需基本源型，接收断层不计算
} FINITE_FAULT;

/**
 * 读取 Coulomb 格式有限断层文件
 *
 * 第一行首个 token 为 #，其后给出 10 个字段标签；第二行为 11 个占位字段
 * 读取前两行 Coulomb 表头，逐行读取 11 列断层数据并建立衍生量
 * 每行末尾支持给定完整 -D 参数以指定各断层的破裂过程
 * 动态震源在读入阶段生成无延迟、已阻尼的 STF 并保存采样对齐后的延迟，未指定时使用脉冲
 * 静态解和接收断层忽略行末的时间函数，接收断层不建立基本源型分量
 * 指定全局 STF 时忽略所有行末内容，各断层使用相同的时间函数和延迟
 * 第 7 个数据列的标签精确为 rake 时按 rake/net slip 格式解释
 * 接收断层仅支持 Kode=100，忽略滑动量大小；rake 格式直接保留第 7 列滑动角，与第 8 列无关
 * 接收断层的非 rake 格式以第 7、8 列作为右旋、逆冲分量确定滑动方向
 * 要求 dip ∈ (0, 90]、bot > top >= 0 且沿走向长度 > 0
 * 调用方负责 grt_finite_fault_free
 *
 * @param[in]   path        文件路径
 * @param[out]  nfault      读入的断层段数
 * @param[in]   is_source   是否为震源断层，接收断层只读取几何和接收机制
 * @param[in]   dt          动态震源采样间隔，s；静态解及接收断层传 0，不生成 STF
 * @param[in]   stf_decay   每个采样点的阻尼因子 exp(-wI*dt)，Lamb 解和静态解传 1
 * @param[in]   stf_option  完整的全局 -D 选项，NULL 时使用行末设置或默认脉冲
 * @return      新分配的 FINITE_FAULT 数组，失败则报错退出
 */
FINITE_FAULT *grt_finite_fault_load_coulomb(const char *path, size_t *nfault, bool is_source, real_t dt, real_t stf_decay, const char *stf_option);

/**
 * 解析有限断层选项并读取 Coulomb 格式有限断层文件
 *
 * 选项格式为 <fault>[+i<dL>/<dW>]，dL/dW 必须同时为零或同时为正值
 * 不提供 +i 时 has_i 为 false，dL/dW 均为零；+i0/0 表示不再剖分
 * 返回的断层数组已经建立衍生量，调用方负责 grt_finite_fault_free
 * 接收断层在此处完成几何剖分，未提供 +i 时每条断层只取中心点，与模块和格林函数库无关
 *
 * @param[in]   option      有限断层选项值，不含 -C 或 -U 选项字符
 * @param[out]  nfault      读入的断层段数
 * @param[out]  has_i       是否显式指定 +i
 * @param[out]  dL          沿走向剖分尺寸 (km)，未指定时为零
 * @param[out]  dW          沿倾向剖分尺寸 (km)，未指定时为零
 * @param[in]   is_source   是否为震源断层，接收断层只读取几何和接收机制
 * @param[in]   dt          动态震源采样间隔，s；静态解及接收断层传 0，不生成 STF
 * @param[in]   stf_decay   每个采样点的阻尼因子 exp(-wI*dt)，Lamb 解和静态解传 1
 * @param[in]   stf_option  完整的全局 -D 选项，NULL 时使用行末设置或默认脉冲
 * @return      新分配的 FINITE_FAULT 数组
 */
FINITE_FAULT *grt_finite_fault_from_option(const char *option, size_t *nfault, bool *has_i, real_t *dL, real_t *dW,
                                           bool is_source, real_t dt, real_t stf_decay, const char *stf_option);

/**
 * 按断层数量释放各断层的剖分结果、时间函数及数组本身
 *
 * @param[in]      nfault  断层数量
 * @param[in,out]  faults  断层数组，可为 NULL
 */
void grt_finite_fault_free(size_t nfault, FINITE_FAULT *faults);

/**
 * 一次性计算并保存全部子断层的中心坐标、长宽及可选的动态标量矩
 *
 * dL/dW 必须同时为正值或同时为非正值，均为非正值时使用整个断层
 * 末块可短于 dL/dW，中心取实际尺寸的中点
 * 作为震源的 Kode=400/500 应传入非正的 dL/dW，只生成一个断层面中心点
 * 通用单点源忽略剖分尺寸，使用原始坐标，长宽均为零
 * 动态有限源显式指定 STF 且提供模型时，按子断层实际中心深度查询介质并计算标量矩，界面取上层
 * 标量矩仅由模型、震源机制和剖分结果确定，不随格林函数库查询方式或时间函数改变
 * 其他情况只剖分几何，scalar_moments 为 NULL
 * 每条断层只初始化一次，结果由 grt_finite_fault_free 释放
 *
 * @param[in,out]  fault   有限断层，需已有走向等衍生量且尚未剖分
 * @param[in]      dL      沿走向剖分间隔 (km)
 * @param[in]      dW      沿倾向剖分间隔 (km)
 * @param[in]      nlayer  原始模型层数，不提供模型时为 0
 * @param[in]      modarr  原始模型矩阵，每行 Thk/Va/Vb/Rho/Qa/Qb，不提供模型时为 NULL
 */
void grt_finite_fault_subdiv(FINITE_FAULT *fault, real_t dL, real_t dW, size_t nlayer, const real_t (*modarr)[GRT_MODARR_NCOL]);

/**
 * 将通用单点源构造成一个断层，水平位置为原点，kode=0
 * @param[in]  depth       震源深度，km
 * @param[in]  type        基本源型
 * @param[in]  scale       地震矩、力或矩势
 * @param[in]  with_mu     是否乘源点剪切模量
 * @param[in]  mechanism   机制参数
 * @param[in]  dt          动态震源采样间隔，s；静态解传 0，不生成 STF
 * @param[in]  stf_decay   每个采样点的阻尼因子 exp(-wI*dt)，Lamb 解和静态解传 1
 * @param[in]  stf_option  完整的 -D 选项，NULL 时使用默认脉冲
 * @return     新分配的单点断层，调用方负责 grt_finite_fault_free
 */
FINITE_FAULT *grt_finite_fault_from_point(real_t depth, GRT_SYN_TYPE type, real_t scale, bool with_mu,
                                        const real_t *mechanism, real_t dt, real_t stf_decay, const char *stf_option);
