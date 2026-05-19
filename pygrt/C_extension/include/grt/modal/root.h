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
 * @param[in,out]   mod1d        模型结构体指针
 * @param[in,out]   eigmet       频散数据结构体指针
 * @param[in]       print_progressbar    是否打印进度条
 */
void grt_get_secular_roots(MODEL1D *mod1d, EIGENV_INFO *eigmet, const bool print_progressbar);