/**
 * @file   grnspec.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-04
 * 
 * 将格林函数频谱用一个结构体包裹
 *   
 *                 
 */

#pragma once

#include "grt/common/const.h"
#include "grt/common/myfftw.h"
#include "grt/common/sacio.h"

/** 不同震中距、不同震源、不同分量的格林函数频谱 */
typedef struct {
    size_t nf;            ///< 总频点数
    real_t *freqs;        ///< 频率数组 freqs[nf]
    size_t nf1;
    size_t nf2;           ///< 待计算的频率索引 
    size_t nr;            ///< 震中距数量
    real_t *rs;           ///< 震中距数组 
    real_t wI;            ///< 虚频率, \f$ \tilde{\omega} =\omega - i \omega_I  \f$ 
    bool keepAllFreq;     ///< 是否计算所有频点，不论频率多低
    bool calc_upar;       ///< 是否计算位移u的空间导数
    
    pcplxChnlGrid *u;     ///< 不同震源不同阶数的格林函数的Z、R、T分量频谱结果
    pcplxChnlGrid *uiz;   ///< 不同震源不同阶数的格林函数的Z、R、T分量对z偏导的频谱结果
    pcplxChnlGrid *uir;   ///< 不同震源不同阶数的格林函数的Z、R、T分量对r偏导的频谱结果

    char *statsstr;       ///< 积分结果输出路径
    size_t  nstatsidxs;   ///< 仅输出特定频点的对应数量
    size_t *statsidxs;    ///< 对应频点的频率索引 statsidxs[nstatsidxs]
} GRNSPEC;

/** 申请 u, uiz, uir 的内存 */
void grt_grnspec_allocate_u(GRNSPEC *grn);

/** 释放 u, uiz, uir 的内存 */
void grt_grnspec_free_u(GRNSPEC *grn);

/**
 * 将频谱 u, uiz, uir 变换到时域后以 SAC 格式保存到本地
 * 
 * @param[in]    grn          格林函数频谱结构体
 * @param[in]    travtPS      不同震源距的初至P、S到时
 * @param[in]    begintimes   不同震中距的波形时移
 * @param[in]    outputdirs   不同震中距的保存目录
 * @param[in]    fh           控制反傅里叶变换的结构体
 * @param[in,out]   sac       SACTRACE 原型，在头段变量中记录了基本信息
 * @param[in]   validChnls    要保存的分量，例如全波解为 "ZRT"， Rayleigh面波解为 "ZR", Love面波解为 "T"
 * @param[in]   skipImagComps 跳过虚频率的补偿
 * @param[in]    saveEX       保存爆炸源结果
 * @param[in]    saveVF       保存垂直力源结果
 * @param[in]    saveHF       保存水平力源结果
 * @param[in]    saveDC       保存水平力源结果
 * 
 */
void grt_grnspec_write_sac(
    const GRNSPEC *grn, const real_t (*travtPS)[2], const real_t *begintimes, char **outputdirs, GRT_FFTW_HOLDER *fh, SACTRACE *sac, 
    const char *validChnls, const bool skipImagComps, const bool saveEX, const bool saveVF, const bool saveHF, const bool saveDC);
