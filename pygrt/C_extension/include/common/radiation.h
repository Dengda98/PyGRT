/**
 * @file   radiation.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-06
 * 
 *    计算不同震源的辐射因子
 * 
 */

#pragma once

#include "common/const.h"

#define GRT_SYN_COMPUTE_EX 0   ///< 计算爆炸源
#define GRT_SYN_COMPUTE_SF 1   ///< 计算单力源
#define GRT_SYN_COMPUTE_DC 2   ///< 计算剪切源
#define GRT_SYN_COMPUTE_MT 3   ///< 计算矩张量源


/**
 * 设置每个震源的方向因子
 * 
 * @param[out]      srcRadi          方向因子，[3]表示ZRT三分量，[6]表示6个震源(EX,VF,HF,DD,DS,SS)
 * @param[in]       computeType      要计算的震源类型，使用宏定义
 * @param[in]       par_theta        方向因子中是否对theta(az)求导
 * @param[in]       M0               放大系数，对于剪切源、爆炸源、张量震源，M0是标量地震矩；对于单力源，M0是放大系数
 * @param[in]       coef             放大系数，用于位移空间导数的计算
 * @param[in]       azrad            弧度制的方位角
 * @param[in]       mchn             震源机制参数，
 *                                   对于单力源，mchn={fn, fe, fz}，
 *                                   对于剪切源，mchn={strike, dip, rake}，
 *                                   对于张量源，mchn={Mxx, Mxy, Mxz, Myy, Myz, Mzz}
 */
void set_source_radiation(
    double srcRadi[SRC_M_NUM][CHANNEL_NUM], const int computeType, const bool par_theta,
    const double M0, const double coef, const double azrad, const double mchn[MECHANISM_NUM]
);