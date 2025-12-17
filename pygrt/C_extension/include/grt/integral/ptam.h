/**
 * @file   ptam.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 峰谷平均法 ，参考：
 * 
 *         1. 张海明. 2021. 地震学中的Lamb问题（上）. 科学出版社
 *         2. Zhang, H. M., Chen, X. F., & Chang, S. (2003). 
 *               An efficient numerical method for computing synthetic seismograms 
 *               for a layered half-space with sources and receivers at close or same depths. 
 *               Seismic motion, lithospheric structures, earthquake and volcanic sources: 
 *               The Keiiti Aki volume, 467-486.
 *                   
 */

#pragma once 

#include <stdio.h>

#include "grt/common/model.h"
#include "grt/integral/kernel.h"
#include "grt/integral/k_integ.h"


/**
 * 峰谷平均法 Peak-Trough Averaging Method，最后收敛的积分结果以三维数组的形式返回，
 * 
 * @param[in,out] mod1d         `GRT_MODEL1D` 结构体指针
 * @param[in]     k0            先前的积分已经进行到了波数k0
 * @param[in]     predk         先前的积分使用的积分间隔dk，因为峰谷平均法使用的
 *                              积分间隔会和之前的不一致，这里传入该系数以做预先调整
 * @param[in]     nr            震中距数量
 * @param[in]     rs            震中距数组  
 * 
 * @param[in,out]  K            用于存储积分的结构体
 * 
 * @param[out]    ptam_fstatsnr      峰谷平均法过程文件指针数组
 * @param[in]     kerfunc            计算核函数的函数指针
 * 
 * 
 */
void grt_PTA_method(
    GRT_MODEL1D *mod1d, real_t k0, real_t predk,
    size_t nr, real_t *rs, K_INTEG *K, FILE *ptam_fstatsnr[nr][2], GRT_KernelFunc kerfunc);

