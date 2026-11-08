/**
 * @file   modal_def.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-01
 * 
 *                   
 */

#pragma once

#include <stdbool.h>

#include "grt/common/const.h"

#define GRT_EIGENV_SATOL  1e-3   ///< 自适应搜根的默认收敛容差
#define GRT_EIGENV_CGAP   1e-7   ///< 默认相速度根间隔，km/s
#define GRT_EIGENV_RTOL   1e-4   ///< 默认久期函数振幅阈值
#define GRT_EIGENV_VGAP   1e-7   ///< 根与模型速度的默认最小间隔，km/s

/** 某频率的频散值 */
typedef struct {
    real_t *c_roots;       ///< 本征值
    real_t *u_roots;
    size_t *c_roots_iref;  ///< 每个本征值来自哪层的久期函数
    size_t n;              ///< 每个频率下找到的零点数量
} EIGENV;

/** 与频散计算方法相关的参数 */
typedef struct {
    size_t nf;
    real_t *freqs;
    size_t nmode;
    size_t *modes;
    DISPER_TYPE wtype;
    bool print_sec;  ///< 仅打印单一频率的久期函数，不搜根

    real_t cmin;   ///< 搜根范围
    real_t cmax;
    bool manual_crange;  ///< 是否手动给定了相速度范围
    size_t iref;
    bool manual_iref;   ///< 是否手动给定了某层的久期函数
    
    // 两个搜根方法的参数
    real_t satol;  ///< 自适应方法的精度
    real_t uniform_dc;  ///< 固定间隔方法的搜索间隔

    real_t cgap;   ///< 搜根时限制不同阶之间相速度最小间隔
    real_t rtol;   ///< 极小值处判断为零点的久期函数幅值的阈值
    real_t vgap;   ///< 零点与模型速度值之间的最小间隔

    size_t neval;  ///< 计算久期函数次数

    EIGENV *eigv;  ///< 不同频率的频散
} EIGENV_INFO;

/** 某频率的某个本征值对应的本征函数和频散敏感核 */
typedef struct {
    real_t eigenC;

    cplx_t (*fn)[4];  // [nz][4]
    
    cplx_t egyint[GRT_EGYINTS_MAX];  ///< 能量积分
    cplx_t (*csens)[GRT_SNSTVTY_MAX];  // 相速度敏感核     [cpar_nz][3]
    cplx_t (*usens)[GRT_SNSTVTY_MAX];  // 群速度速度敏感核 [cpar_nz][3]
} EIGENFN;


typedef struct {
    size_t nf;
    real_t *freqs;
    size_t nmode;
    size_t *modes;
    DISPER_TYPE wtype;

    size_t nz;
    real_t *zs;
    size_t *z_irefs;

    size_t cpar_nz;
    real_t *cpar_zs;
    size_t *cpar_z_irefs;

    EIGENV *eigv;     ///< （筛选后）不同频率的频散
    EIGENFN **eigfn;  // [nf][cnum]
} EIGENFN_INFO;

/** 本征函数、能量积分和敏感核的计算选项，零初始化表示不请求额外计算 */
typedef struct {
    bool eigenfunctions;         ///< 是否计算深度采样的本征函数
    bool energy;                 ///< 是否计算能量积分和群速度
    bool phase_kernel;           ///< 是否计算相速度核，同时计算能量积分和群速度
    bool group_kernel;           ///< 是否计算群速度核，同时计算相速度核
    real_t group_kernel_shift;   ///< 群速度核的正频率偏移量，Hz，0 表示使用相邻频点差分
    real_t kernel_dz;            ///< 敏感核分段间距，0 表示直接使用物理层
    size_t nz;                   ///< 本征函数采样深度数量，0 表示使用物理层界面
    const real_t *zs;            ///< 本征函数采样深度 [nz]，严格升序，km，仅借用
} EIGENFN_OPTIONS;

/**
 * 初始化相速度计算结果，复制频率并设置默认搜根参数
 *
 * @param[in]      nf                   频率数量，至少为 1
 * @param[in]      freqs                正频率 [nf]，严格升序或严格降序，Hz
 * @param[in]      wtype                Rayleigh 或 Love 波
 * @param[in]      nmode                搜索的阶数数量，1 表示基阶，0 表示尽可能搜索所有阶
 * @return   新分配的结果，使用 grt_free_eigenv_info 释放
 */
EIGENV_INFO *grt_init_eigenv_info(size_t nf, const real_t *freqs, DISPER_TYPE wtype, size_t nmode);

/**
 * 复制全部频散结果，初始化尚未计算本征函数的结果对象
 *
 * @param[in]      eigmet               grt_compute_eigenv 或 grt_read_dispersion 得到的结果
 * @return   独立拥有频率和频散数组的结果，使用 grt_free_eigenfn_info 释放
 */
EIGENFN_INFO *grt_init_eigenfn_info(const EIGENV_INFO *eigmet);

/** 释放 eigv 内部内存 */
void grt_free_eigenv(EIGENV *eigv);

/** 释放 eigmet 内部及其本身内存 */
void grt_free_eigenv_info(EIGENV_INFO *eigmet);

/** 释放 eigfnmet 内部及其本身内存 */
void grt_free_eigenfn_info(EIGENFN_INFO *eigfnmet);

/** 从全部频散信息中，筛选出部分信息 */
void grt_filter_eigenfn_info(
    const size_t nf, const real_t *freqs, const bool def_freq_range, 
    const size_t nmode, const size_t *modes, 
    const EIGENV_INFO *eigmet, EIGENFN_INFO *eigfnmet);
