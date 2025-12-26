/**
 * @file   syn.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-12
 * 
 *    将合成地震图部分单独用一个源文件和若干函数来管理
 * 
 */

#pragma once

#include "grt/common/const.h"
#include "grt/common/sacio2.h"

/**
 * 根据已有系数（方向因子）合成理论地震图
 * 
 * @param[in]      srcRadi    不同震源不同分量的系数
 * @param[in]      computeType   要计算的震源类型，使用宏定义
 * @param[in]      dirpath    格林函数所在目录路径
 * @param[in]      prefix     格林函数文件名前缀
 * @param[out]     synsac     三分量 SACTRACE
 */
void grt_syn(const realChnlGrid srcRadi, const int computeType, const char *dirpath, const char *prefix, SACTRACE *synsac[GRT_CHANNEL_NUM]);