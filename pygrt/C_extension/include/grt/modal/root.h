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

/**
 * 按模型自动确定搜索范围并计算频散，支持严格升序或严格降序频率
 * 使用 grt_init_eigenv_info 初始化结果，可在调用前覆盖默认搜根参数
 * 每个结果对象仅计算一次，完成后 nmode 和 modes 记录实际找到的阶数
 * 外层已并行时内部按频率串行，否则使用 OpenMP 并行
 *
 * @param[in,out]  mod1d                分层模型，要求自由表面和底部半空间
 * @param[in,out]  eigmet               搜根选项及频散结果，保持输入频率顺序
 * @param[in]      print_progressbar    是否显示进度
 */
void grt_compute_eigenv(MODEL1D *mod1d, EIGENV_INFO *eigmet, bool print_progressbar);
