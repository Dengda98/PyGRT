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
#include "grt/common/kernel.h"



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
 * @param[out]    sum_J0        积分值
 * 
 * @param[in]     calc_upar       是否计算位移u的空间导数
 * @param[out]    sum_uiz_J0      uiz的积分值
 * @param[out]    sum_uir_J0      uir的积分值
 * 
 * @param[out]    ptam_fstatsnr      峰谷平均法过程文件指针数组
 * @param[in]     kerfunc            计算核函数的函数指针
 * 
 * 
 */
void grt_PTA_method(
    GRT_MODEL1D *mod1d, real_t k0, real_t predk,
    size_t nr, real_t *rs,
    cplx_t sum_J0[nr][GRT_SRC_M_NUM][GRT_INTEG_NUM],
    bool calc_upar,
    cplx_t sum_uiz_J0[nr][GRT_SRC_M_NUM][GRT_INTEG_NUM],
    cplx_t sum_uir_J0[nr][GRT_SRC_M_NUM][GRT_INTEG_NUM],
    FILE *ptam_fstatsnr[nr][2], GRT_KernelFunc kerfunc);





/**
 * 观察连续3个点的函数值的实部变化，判断是波峰(1)还是波谷(-1), 并计算对应值。
 * 
 * @param[in]     idx1        阶数索引
 * @param[in]     idx2        积分类型索引 
 * @param[in]     arr         存有连续三个点的函数值的数组 
 * @param[in]     k           三个点的起始波数
 * @param[in]     dk          三个点的波数间隔，这样使用k和dk定义了三个点的位置
 * @param[out]    pk          估计的波峰或波谷处的波数
 * @param[out]    value       估计的波峰或波谷处的函数值
 * 
 * @return    波峰(1)，波谷(-1)，其它(0)
 *  
 */
int grt_cplx_peak_or_trough(
    int idx1, int idx2, const cplx_t arr[GRT_PTAM_WINDOW_SIZE][GRT_SRC_M_NUM][GRT_INTEG_NUM], 
    real_t k, real_t dk, real_t *pk, cplx_t *value);


/**
 * 递归式地计算缩减序列的值，
 * \f[
 * M_i = 0.5\times (M_i + M_{i+1})
 * \f]
 * 
 * @param[in]         n1          数组长度 
 * @param[in,out]     arr         振荡的数组，最终收敛值在第一个，arr[0] 
 * 
 */
void grt_cplx_shrink(size_t n1, cplx_t *arr);