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
 * 对频谱执行 IFFT，按统一格式将时域波形保存为 SAC 文件
 * 全波解和面波解共用此接口，refhead 保存当前源深和接收深度的介质头段
 * @param[in]      refhead        当前源台深度组合的 SAC 参考头段
 * @param[in]      grn            当前深度组合的格林函数频谱
 * @param[in]      travtPS        各震中距的初至 P、S 到时
 * @param[in]      begintimes     各震中距的波形起点，s
 * @param[in]      outputdirs     各震中距的 SAC 输出目录
 * @param[in,out]  fh             反傅里叶变换缓冲
 * @param[in]      validChnls     保存的分量，全波解为 ZRT，Rayleigh 为 ZR，Love 为 T
 * @param[in]      skipImagComps  是否跳过虚频率补偿
 * @param[in]      saveEX         是否保存爆炸源
 * @param[in]      saveVF         是否保存垂直力源
 * @param[in]      saveHF         是否保存水平力源
 * @param[in]      saveDC         是否保存剪切源
 */
void grt_grnspec_save_waveforms(
    const SACHEAD *refhead, const GRNSPEC *grn, const real_t (*travtPS)[2], const real_t *begintimes,
    char *const *outputdirs, FFTW_HOLDER *fh,
    const char *validChnls, bool skipImagComps, bool saveEX, bool saveVF, bool saveHF, bool saveDC);
