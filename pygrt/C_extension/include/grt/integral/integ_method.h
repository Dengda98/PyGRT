/**
 * @file   integ_method.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-12
 * 
 * 使用结构体来管理与波数积分相关的参数和方法
 * 
 */

#pragma once

#include <stdio.h>
#include <stdbool.h>

#include "grt/common/const.h"
#include "grt/integral/k_integ.h"
#include "grt/integral/kernel.h"


// 描述不同波数积分方法的结构体
typedef struct {
    real_t k0;      ///< 波数积分的上限 \f$ \tilde{k_{max}}=\sqrt{(k_{0}*\pi/hs)^2 + (ampk*w/vmin_{ref})^2} \f$ ，k循环必须退出, hs=max(震源和台站深度差,1.0) 
    real_t ampk;    ///< 影响波数k积分上限的系数
    real_t keps;    ///< 波数积分的收敛条件，要求在某震中距下所有格林函数都收敛，为负数代表不提前判断收敛，按照波数积分上限进行积分 
    real_t vmin;    ///< 参考最小速度，用于定义波数积分的上限

    real_t kcut;    ///< 波数积分和Filon积分的分割点

    real_t kmax;    ///< 全局波数最大值，程序运行中会随频率变动

    real_t dk;      ///< DWM 的波数积分间隔

    bool applyFIM;  ///< 是否使用 FIM
    real_t filondk; ///< FIM 的波数积分间隔

    bool applySAFIM;  ///< 是否使用 SAFIM
    real_t sa_tol;    ///< SAFIM 的收敛极限
    
    // 积分显式收敛方法
    bool applyDCM;    ///< 是否使用 DCM
    bool applyPTAM;   ///< 是否使用 PTAM
    
    FILE *fstats;  ///< 保存核函数的文件指针
    FILE *(*ptam_fstatsnr)[2];  ///< 保存 PTAM 中核函数以及波峰波谷的文件指针
    
} K_INTEG_METHOD;


/**
 * 初始化 K_INTEG_METHOD 结构体中的文件指针
 * 
 * @param[in]          nr           震中距数量
 * @param[in]          rs           震中距数组 
 * @param[in]          statsstr     积分过程输出目录
 * @param[in]          suffix       文件名后缀
 * @param[in,out]      Kmet         K_INTEG_METHOD 结构体指针
 * 
 */
void grt_KMET_init_fstats(
    const size_t nr, const real_t *rs, 
    const char *statsstr, const char *suffix, K_INTEG_METHOD *Kmet);

/**
 * 关闭 K_INTEG_METHOD 结构体中的文件指针，并释放内存
 * 
 * @param[in]          nr           震中距数量
 * @param[in,out]      Kmet         K_INTEG_METHOD 结构体指针
 */
void grt_KMET_destroy_fstats(const size_t nr, K_INTEG_METHOD *Kmet);




/**
 * 发起波数积分的总函数
 * 
 * @param[in,out]      mod1d        `GRT_MODEL1D` 结构体指针
 * @param[in]          nr           震中距数量
 * @param[in]          rs           震中距数组 
 * @param[in,out]      Kmet         K_INTEG_METHOD 结构体指针
 * @param[in]          calc_upar    是否计算位移u的空间导数
 * @param[in]          kerfunc      计算核函数的函数指针
 * 
 * @return    K_INTEG 结构体指针，对应着积分结果
 * 
 */
K_INTEG * grt_wavenumber_integral(
    GRT_MODEL1D *mod1d, size_t nr, real_t *rs, K_INTEG_METHOD *Kmet, bool calc_upar, GRT_KernelFunc kerfunc);
