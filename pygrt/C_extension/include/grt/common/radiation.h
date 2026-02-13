/**
 * @file   radiation.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-06
 * 
 *    计算不同震源的辐射因子
 * 
 */

#pragma once

#include <stdbool.h>

#include "grt/common/const.h"

/**
 * 设置每个震源的方向因子
 * 
 * @param[out]      srcRadi          方向因子，[3]表示ZRT三分量，[6]表示6个震源(EX,VF,HF,DD,DS,SS)
 * @param[in]       computeType      要计算的震源类型
 * @param[in]       par_theta        方向因子中是否对theta(az)求导
 * @param[in]       M0               放大系数，对于剪切源、爆炸源、张量震源，M0是标量地震矩；对于单力源，M0是放大系数
 * @param[in]       coef             放大系数，用于位移空间导数的计算
 * @param[in]       VpVs_ratio       震源层的 Vp/Vs 比值，用于张裂源
 * @param[in]       azrad            弧度制的方位角
 * @param[in]       mchn             震源机制参数，
 *                                   对于单力源，mchn={fn, fe, fz}，
 *                                   对于剪切源，mchn={strike, dip, rake}，
 *                                   对于张量源，mchn={Mxx, Mxy, Mxz, Myy, Myz, Mzz}
 */
void grt_set_source_radiation(
    realChnlGrid srcRadi, const GRT_SYN_TYPE computeType, const bool par_theta,
    const real_t M0, const real_t coef, const real_t VpVs_ratio, const real_t azrad, const real_t mchn[GRT_MECHANISM_NUM]
);