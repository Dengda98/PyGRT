/**
 * @file   saroot.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-05
 * 
 *     使用自适应采样寻找Rayleigh波和Love波的久期函数零点
 * 
 */


#pragma once 

#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/modal/modal_def.h"


/**
 * 使用自适应采样寻找久期函数零点
 * 
 * @param[in]      mod1d        模型结构体指针
 * @param[in]      V_sort       升序的模型速度
 * @param[in]      w            圆频率
 * @param[in]      c10          相速度搜索下界
 * @param[in]      tol          自适应收敛精度
 * @param[in]      c20          相速度搜索上界
 * @param[in]      isRayl       Rayleigh or Love
 * @param[in]      secRaylType  Rayl久期函数类型
 * @param[in]      iref         久期函数层位
 * @param[in]      print_sec    是否打印久期函数
 * @param[in]      nmode        最大零点个数, 如果 nmode<=0 则找全部零点
 * @param[out]     pt_c_roots   存储零点的本征值
 * @param[out]     pt_c_roots_iref   存储所使用的久期函数层位
 * @param[out]     Nroot        最终零点数量
 * 
 */
void grt_sa_secular_roots(GRT_MODEL1D *mod1d, EIGENV_METHOD *eigmet, EIGENV *eigv, const int secRaylType, const size_t iref);