/**
 * @file   grt.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 广义反射透射系数矩阵+离散波数法 计算理论地震图，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient 
 *               matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699
 *   
 *                 
 */

#pragma once 

#include "model.h"

/**
 * 设置OpenMP多线程数
 * 
 * @param   num_threads    (in)线程数
 */
void set_num_threads(int num_threads);


/**
 * 积分计算Z, R, T三个分量格林函数的频谱的核心函数（被C函数调用）
 * 
 * @param      pymod1d      (in)`PYMODEL1D` 结构体指针 
 * @param      nf1          (in)开始计算频谱的频率索引值, 总范围在[nf1, nf2]
 * @param      nf2          (in)结束计算频谱的频率索引值, 总范围在[nf1, nf2]
 * @param      nf           (in)所有频点个数
 * @param      freqs        (in)所有频点的频率值（包括未计算的）
 * @param      nr           (in)震中距数量
 * @param      rs           (in)震中距数组 
 * @param      wI           (in)虚频率, \f$ \tilde{\omega} =\omega - i\wI \f$ 
 * @param      vmin_ref     (in)参考最小速度，用于定义波数积分的上限
 * @param      keps         (in)波数积分的收敛条件，要求在某震中距下所有格林函数都收敛，为负数代表不提前判断收敛，按照波数积分上限进行积分 
 * @param      ampk         (in)影响波数k积分上限的系数，见下方
 * @param      iwk0         (in)k0是否取随频率变化的线性关系，即 \f$ k_{max} = k_{max} * f/f_{max} \f$ 
 * @param      k0           (in)波数积分的上限 \f$ \tilde{k_{max}}=\sqrt{(k_{0}*\pi/hs)^2 + (ampk*w/vmin_{ref})^2} \f$ ，k循环必须退出, hs=max(震源和台站深度差,1.0) 
 * @param      Length       (in)波数k积分间隔 \f$ dk=2\pi/(fabs(L)*r_{max}) \f$ , 如果为负数，则使用线性Filon积分
 * @param      print_progressbar    (in)是否打印进度条
 * @param      EXPgrn[nr][2]      (out)复数数组，爆炸源的Z、R分量频谱结果
 * @param      VFgrn[nr][2]       (out)复数数组，垂直力源的Z、R分量频谱结果
 * @param      HFgrn[nr][3]       (out)复数数组，水平力源的Z、R、T分量频谱结果
 * @param      DDgrn[nr][2]       (out)复数数组，45度倾滑的Z、R分量频谱结果
 * @param      DSgrn[nr][3]       (out)复数数组，90度倾滑的Z、R、T分量频谱结果
 * @param      SSgrn[nr][3]       (out)复数数组，90度走滑的Z、R、T分量频谱结果
 * 
 * @param      statsstr          (in) 积分结果输出目录
 * @param      nstatsidxs        (in) 输出积分结果的特定频点的个数
 * @param      statsidxs         (in) 特定频点的索引值
 * 
 */ 
void integ_grn_spec_in_C(
    PYMODEL1D *pymod1d, MYINT nf1, MYINT nf2, MYINT nf, MYREAL *freqs,  
    MYINT nr, MYREAL *rs, MYREAL wI, 
    MYREAL vmin_ref, MYREAL keps, MYREAL ampk, bool iwk0, MYREAL k0, MYREAL Length,       
    bool print_progressbar, 

    // 返回值，维度2代表Z、R分量，维度3代表Z、R、T分量
    MYCOMPLEX *EXPcplx[nr][2], // EXZ, EXR 的实部和虚部
    MYCOMPLEX *VFcplx[nr][2],  // VFZ, VFR 的实部和虚部
    MYCOMPLEX *HFcplx[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    MYCOMPLEX *DDcplx[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    MYCOMPLEX *DScplx[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    MYCOMPLEX *SScplx[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]

    const char *statsstr, // 积分结果输出
    MYINT  nstatsidxs, // 仅输出特定频点
    MYINT *statsidxs
);


/**
 * 积分计算Z, R, T三个分量格林函数的频谱的核心函数（被Python调用）  
 * 
 * @param      pymod1d      (in)`PYMODEL1D` 结构体指针 
 * @param      nf1          (in)开始计算频谱的频率索引值, 总范围在[nf1, nf2]
 * @param      nf2          (in)结束计算频谱的频率索引值, 总范围在[nf1, nf2]
 * @param      nf           (in)所有频点个数
 * @param      freqs        (in)所有频点的频率值（包括未计算的）
 * @param      nr           (in)震中距数量
 * @param      rs           (in)震中距数组 
 * @param      wI           (in)虚频率, \f$ \tilde{\omega} =\omega - i\wI \f$ 
 * @param      vmin_ref     (in)参考最小速度，用于定义波数积分的上限
 * @param      keps         (in)波数积分的收敛条件，要求在某震中距下所有格林函数都收敛，为负数代表不提前判断收敛，按照波数积分上限进行积分 
 * @param      ampk         (in)影响波数k积分上限的系数，见下方
 * @param      iwk0         (in)k0是否取随频率变化的线性关系，即 \f$ k_{max} = k_{max} * f/f_{max} \f$ 
 * @param      k0           (in)波数积分的上限 \f$ \tilde{k_{max}}=\sqrt{(k_{0}*\pi/hs)^2 + (ampk*w/vmin_{ref})^2} \f$ ，k循环必须退出, hs=max(震源和台站深度差,1.0) 
 * @param      Length       (in)波数k积分间隔 \f$ dk=2\pi/(fabs(L)*r_{max}) \f$ , 如果为负数，则使用线性Filon积分
 * @param      print_progressbar    (in)是否打印进度条
 * @param      EXPgrn[nr][2]      (out)`GRN` 结构体指针，爆炸源的Z、R分量频谱结果
 * @param      VFgrn[nr][2]       (out)`GRN` 结构体指针，垂直力源的Z、R分量频谱结果
 * @param      HFgrn[nr][3]       (out)`GRN` 结构体指针，水平力源的Z、R、T分量频谱结果
 * @param      DDgrn[nr][2]       (out)`GRN` 结构体指针，45度倾滑的Z、R分量频谱结果
 * @param      DSgrn[nr][3]       (out)`GRN` 结构体指针，90度倾滑的Z、R、T分量频谱结果
 * @param      SSgrn[nr][3]       (out)`GRN` 结构体指针，90度走滑的Z、R、T分量频谱结果
 * 
 * @param      statsstr          (in) 积分结果输出目录
 * @param      nstatsidxs        (in) 输出积分结果的特定频点的个数
 * @param      statsidxs         (in) 特定频点的索引值
 * 
 */ 
void integ_grn_spec(
    PYMODEL1D *pymod1d, MYINT nf1, MYINT nf2, MYINT nf, MYREAL *freqs,  
    MYINT nr, MYREAL *rs, MYREAL wI, 
    MYREAL vmin_ref, MYREAL keps, MYREAL ampk, bool iwk0, MYREAL k0, MYREAL Length,       
    bool print_progressbar, 

    // 返回值，维度2代表Z、R分量，维度3代表Z、R、T分量
    GRN *EXPgrn[nr][2], // EXZ, EXR 的实部和虚部
    GRN *VFgrn[nr][2],  // VFZ, VFR 的实部和虚部
    GRN *HFgrn[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    GRN *DDgrn[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    GRN *DSgrn[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    GRN *SSgrn[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]

    const char *statsstr, // 积分结果输出
    MYINT  nstatsidxs, // 仅输出特定频点
    MYINT *statsidxs
);




