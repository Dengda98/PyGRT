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
 * @param[out]  srcRadi      方向因子，[3]表示ZRT三分量，[6]表示6个震源(EX,VF,HF,DD,DS,SS)
 * @param[in]   source_type  要计算的震源类型
 * @param[in]   par_theta    方向因子中是否对theta(az)求导
 * @param[in]   scale        放大系数，对于剪切源、爆炸源、张量震源，scale 是标量地震矩；对于单力源，scale 是力
 * @param[in]   coef         附加放大系数，用于位移空间导数计算或量纲换算
 * @param[in]   VpVs_ratio   震源层的 Vp/Vs 比值，用于张裂源
 * @param[in]   azrad        弧度制的方位角
 * @param[in]   mchn         震源机制参数，
 *                                   对于单力源，mchn={fn, fe, fz}，
 *                                   对于剪切源，mchn={strike, dip, rake}，
 *                                   对于张量源，mchn={Mxx, Mxy, Mxz, Myy, Myz, Mzz}
 */
void grt_set_source_radiation(
    realChnlGrid srcRadi, const GRT_SYN_TYPE source_type, const bool par_theta,
    const real_t scale, const real_t coef, const real_t VpVs_ratio, const real_t azrad, const real_t mchn[GRT_MECHANISM_NUM]
);
/**
 * 按辐射系数的源强约定构造 NED 矩张量
 * @param[in]   type        源类型，不用于单力源
 * @param[in]   scale       源强，dyne-cm
 * @param[in]   VpVs_ratio  源点 Vp/Vs
 * @param[in]   mchn        震源机制
 * @param[out]  tensor      Mxx、Mxy、Mxz、Myy、Myz、Mzz
 */
void grt_source_moment_tensor(GRT_SYN_TYPE type, real_t scale, real_t VpVs_ratio, const real_t mchn[GRT_MECHANISM_NUM], real_t tensor[6]);


/**
 * 判断震源类型是否使用某个基本格林函数源型
 * @param[in]  source_type  震源类型
 * @param[in]  im           基本格林函数源型索引
 */
bool grt_source_has_component(GRT_SYN_TYPE source_type, int im);

/**
 * 检查液体震源介质只允许不乘剪切模量的爆炸源
 * @param[in]  source_type  震源类型
 * @param[in]  with_mu      是否按 -Su 乘源点剪切模量
 * @param[in]  vs           源点 S 波速度，km/s
 */
void grt_check_source_medium(GRT_SYN_TYPE source_type, bool with_mu, real_t vs);
