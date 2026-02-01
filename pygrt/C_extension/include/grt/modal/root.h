/**
 * @file   root.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-05
 * 
 *     总函数：寻找Rayleigh波和Love波的久期函数零点
 * 
 */

#pragma once 

#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/modal/modal_def.h"

/**
 * 寻找 Rayleigh 波和 Love 波的久期函数零点的总函数
 * 
 * @param[inout]   mod1d        模型结构体指针
 * @param[in]      nfreq        频点数
 * @param[in]      freqs        频率值
 * @param[in]      nmode        最大零点数（<=0 则找全部零点）
 * @param[in]      isRayl       Rayleigh or Love
 * @param[in]      tol          自适应采样精度
 * @param[in]      print_progressbar    是否打印进度条
 * @param[in]      secfreq_idx  >=0 则打印对应频点索引值的久期函数
 * @param[in]      cmin         相速度搜索下界
 * @param[in]      cmax         相速度搜索上界
 * @param[out]     c_roots      不同频率的零点本征值
 * @param[out]     c_roots_iref    不同频率的零点所使用的久期函数层位
 * @param[out]     c_roots_num  不同频率的零点数量
 * 
 */
void grt_get_secular_roots(GRT_MODEL1D *mod1d, EIGENV_METHOD *eigmet, const bool print_progressbar);
// void grt_get_secular_roots(
//     GRT_MODEL1D *mod1d, EIGENV_METHOD *eigmet,
//     const bool print_progressbar, const size_t secfreq_idx);