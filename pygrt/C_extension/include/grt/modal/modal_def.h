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

/** 用于存储频散值的结构体 */
typedef struct {
    real_t dcmin;   ///< 设定搜根间隔的范围
    real_t dcmax;
    real_t *c_roots;       ///< 本征值， 两个维度分别为频率和根数
    size_t *c_roots_iref;  ///< 每个本征值来自哪层的久期函数
    size_t n;              ///< 每个频率下找到的零点数量
} EIGENV;

/** 与频散计算方法相关的参数 */
typedef struct {
    bool io_period;
    size_t nf;
    real_t *freqs;
    size_t nmode;    ///< 最高阶数，例如最高为2阶， 则 nmode=3 (基阶，1阶，2阶)
    DISPER_TYPE wtype;
    bool print_sec;  ///< 仅打印单一频率的久期函数，不搜根

    real_t cmin;   ///< 搜根范围，若为负数则由程序自行确定
    real_t cmax;
    
    real_t satol;  ///< 自适应方法的精度
    real_t cgap;   ///< 搜根时限制不同阶之间相速度最小间隔
    real_t rtol;   ///< 极小值处判断为零点的久期函数幅值的阈值
    real_t vgap;   ///< 零点与模型速度值之间的最小间隔

    EIGENV *eigv;  ///< 不同频率的频散
} EIGENV_METHOD;

void grt_free_eigenv(EIGENV *eigv);

void grt_free_eigenv_method(EIGENV_METHOD *eigmet);