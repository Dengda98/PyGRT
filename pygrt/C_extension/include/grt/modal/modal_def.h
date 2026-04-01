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

/** 某频率的频散值 */
typedef struct {
    real_t *c_roots;       ///< 本征值， 两个维度分别为频率和根数
    real_t *g_roots;
    size_t *c_roots_iref;  ///< 每个本征值来自哪层的久期函数
    size_t n;              ///< 每个频率下找到的零点数量
} EIGENV;

/** 与频散计算方法相关的参数 */
typedef struct {
    size_t nf;
    real_t *freqs;
    size_t nmode;    ///< 最高阶数，例如最高为2阶， 则 nmode=3 (基阶，1阶，2阶)
    DISPER_TYPE wtype;
    bool print_sec;  ///< 仅打印单一频率的久期函数，不搜根

    real_t cmin;   ///< 搜根范围
    real_t cmax;
    
    // 两个搜根方法的参数
    real_t satol;  ///< 自适应方法的精度
    real_t uniform_dc;  ///< 固定间隔方法的搜索间隔

    real_t cgap;   ///< 搜根时限制不同阶之间相速度最小间隔
    real_t rtol;   ///< 极小值处判断为零点的久期函数幅值的阈值
    real_t vgap;   ///< 零点与模型速度值之间的最小间隔

    size_t neval;

    EIGENV *eigv;  ///< 不同频率的频散
} EIGENV_INFO;

/** 某频率的某个本征值对应的本征函数和频散敏感核 */
typedef struct {
    real_t eigenC;

    cplx_t (*fn)[4];  // [nz][4]
    
    cplx_t egyint[GRT_EGYINTS_MAX];  ///< 能量积分
    cplx_t (*csens)[GRT_SNSTVTY_MAX];  // 相速度敏感核     [cpar_nz][3]
    cplx_t (*gsens)[GRT_SNSTVTY_MAX];  // 群速度速度敏感核 [cpar_nz][3]
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

void grt_free_eigenv(EIGENV *eigv);

void grt_free_eigenv_info(EIGENV_INFO *eigmet);

void grt_filter_eigenfn_info(
    const size_t nf, const real_t *freqs, const bool def_freq_range, 
    const size_t nmode, const size_t *modes, 
    EIGENV_INFO *eigmet, EIGENFN_INFO *eigfnmet);