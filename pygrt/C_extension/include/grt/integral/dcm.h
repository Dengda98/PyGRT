/**
 * @file   dcm.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-12
 * 
 *     DCM 的校正系数
 *                   
 */

#pragma once

#include "grt/common/const.h"
#include "grt/integral/k_integ.h"

/**
 * 对 DCM 结果进行校正
 * 
 * @param[in]        nr      震中距数
 * @param[in]        rs      震中距数组
 * @param[in,out]    Kint    K_INTEG 结构体指针，对应着积分结果
 * @param[in]        keep_nearfield   是否要作用于近场项(p==1项)
 */
void grt_dcm_correction(size_t nr, real_t *rs, K_INTEG *Kint, bool keep_nearfield);