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

#include "common/model.h"

/**
 * 设置OpenMP多线程数
 * 
 * @param[in]   num_threads        线程数
 */
void set_num_threads(int num_threads);



/**
 * 积分计算Z, R, T三个分量格林函数的频谱的核心函数（被Python调用）  
 * 
 * @param[in]      pymod1d          `PYMODEL1D` 结构体指针 
 * @param[in]      nf1              开始计算频谱的频率索引值, 总范围在[nf1, nf2]
 * @param[in]      nf2              结束计算频谱的频率索引值, 总范围在[nf1, nf2]
 * @param[in]      freqs            所有频点的频率值（包括未计算的）
 * @param[in]      nr               震中距数量
 * @param[in]      rs               震中距数组 
 * @param[in]      wI               虚频率, \f$ \tilde{\omega} =\omega - i \omega_I  \f$ 
 * @param[in]      vmin_ref         参考最小速度，用于定义波数积分的上限
 * @param[in]      keps             波数积分的收敛条件，要求在某震中距下所有格林函数都收敛，为负数代表不提前判断收敛，按照波数积分上限进行积分 
 * @param[in]      ampk             影响波数k积分上限的系数，见下方
 * @param[in]      k0               波数积分的上限 \f$ \tilde{k_{max}}=\sqrt{(k_{0}*\pi/hs)^2 + (ampk*w/vmin_{ref})^2} \f$ ，k循环必须退出, hs=max(震源和台站深度差,1.0) 
 * @param[in]      Length           波数k积分间隔 \f$ dk=2\pi/(fabs(L)*r_{max}) \f$ 
 * @param[in]      filonLength      Filon积分间隔
 * @param[in]      safilonTol       自适应Filon积分采样精度
 * @param[in]      filonCut         波数积分和Filon积分的分割点
 * @param[in]      print_progressbar        是否打印进度条
 * 
 * @param[out]      grn               不同震源不同阶数的格林函数的Z、R、T分量频谱结果
 * 
 * @param[in]       calc_upar         是否计算位移u的空间导数
 * @param[out]      grn_uiz           不同震源不同阶数的ui_z的Z、R、T分量频谱结果
 * @param[out]      grn_uir           不同震源不同阶数的ui_r的Z、R、T分量频谱结果
 * 
 * @param[in]       statsstr               积分过程输出目录
 * @param[in]       nstatsidxs             输出积分过程的特定频点的个数
 * @param[in]       statsidxs              特定频点的索引值
 * 
 */ 
void integ_grn_spec(
    PYMODEL1D *pymod1d, MYINT nf1, MYINT nf2, MYREAL *freqs,  
    MYINT nr, MYREAL *rs, MYREAL wI, 
    MYREAL vmin_ref, MYREAL keps, MYREAL ampk, MYREAL k0, MYREAL Length, MYREAL filonLength, MYREAL safilonTol, MYREAL filonCut,             
    bool print_progressbar, 

    // 返回值，代表Z、R、T分量
    MYCOMPLEX *grn[nr][SRC_M_NUM][CHANNEL_NUM],

    bool calc_upar,
    MYCOMPLEX *grn_uiz[nr][SRC_M_NUM][CHANNEL_NUM],
    MYCOMPLEX *grn_uir[nr][SRC_M_NUM][CHANNEL_NUM],

    const char *statsstr, // 积分过程输出
    MYINT  nstatsidxs, // 仅输出特定频点
    MYINT *statsidxs
);




