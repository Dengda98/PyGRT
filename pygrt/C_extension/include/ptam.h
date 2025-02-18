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

#include "model.h"


#define PTAM_MAX_PEAK_TROUGH 36  ///< 最后统计波峰波谷的目标数量



/**
 * 峰谷平均法 Peak-Trough Averaging Method，最后收敛的积分结果以三维数组的形式返回，
 * 形状为[nr][3][4], 分别代表震中距、阶数(m=0,1,2) 和4种积分类型(p=0,1,2,3)  
 * 
 * @param    mod1d     (in)`MODEL1D` 结构体指针
 * @param    k0        (in)先前的积分已经进行到了波数k0
 * @param    predk     (in)先前的积分使用的积分间隔dk，因为峰谷平均法使用的
 *                     积分间隔会和之前的不一致，这里传入该系数以做预先调整
 * @param    rmin      (in)最小震中距 
 * @param    rmax      (in)最大震中距 
 * @param    omega     (in)复数频率 
 * @param    nr        (in)震中距数量
 * @param    rs        (in)震中距数组  
 * @param    sum_EXP_J0[nr][3][4]   (out)爆炸源
 * @param    sum_VF_J0[nr][3][4]    (out)垂直力源
 * @param    sum_HF_J0[nr][3][4]    (out)水平力源
 * @param    sum_DC_J0[nr][3][4]    (out)双力偶源
 * @param    fstats[nr]          (out)波数积分过程文件指针
 * @param    ptam_fstats[nr]     (out)峰谷平均法过程文件指针
 * 
 * 
 */
void PTA_method(
    const MODEL1D *mod1d, MYREAL k0, MYREAL predk, MYREAL rmin, MYREAL rmax, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs, 
    MYCOMPLEX sum_EXP_J0[nr][3][4], MYCOMPLEX sum_VF_J0[nr][3][4],  
    MYCOMPLEX sum_HF_J0[nr][3][4],  MYCOMPLEX sum_DC_J0[nr][3][4],  
    FILE *(fstats[nr]), FILE *(ptam_fstats[nr]));





/**
 * 观察连续3个点的函数值的实部变化，判断是波峰(1)还是波谷(-1), 并计算对应值。
 * 其中存储函数值的数组形状为[3][3][4], 分别代表
 * 连续3个点、阶数(m=0,1,2) 和4种积分类型(p=0,1,2,3)  
 * 
 * @param     idx1    (in)阶数索引
 * @param     idx2    (in)积分类型索引 
 * @param     arr[3][3][4]   (in)存有连续三个点的函数值的数组 
 * @param     k       (in)三个点的起始波数
 * @param     dk      (in)三个点的波数间隔，这样使用k和dk定义了三个点的位置
 * @param     pk      (out)估计的波峰或波谷处的波数
 * @param     value   (out)估计的波峰或波谷处的函数值
 * 
 * @return    波峰(1)，波谷(-1)，其它(0)
 *  
 */
MYINT cplx_peak_or_trough(
    MYINT idx1, MYINT idx2, const MYCOMPLEX arr[3][3][4], 
    MYREAL k, MYREAL dk, MYREAL *pk, MYCOMPLEX *value);


/**
 * 递归式地计算缩减序列的值，
 * \f[
 * M_i = 0.5\times (M_i + M_{i+1})
 * \f]
 * 
 * @param     n1      (in)数组长度 
 * @param     arr     (inout)振荡的数组，最终收敛值在第一个，arr[0] 
 * 
 */
void cplx_shrink(MYINT n1, MYCOMPLEX *arr);