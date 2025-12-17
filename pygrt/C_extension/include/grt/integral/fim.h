/**
 * @file   fim.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是基于线性插值的Filon积分，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.   
 *         2. 纪晨, 姚振兴. 1995. 区域地震范围的宽频带理论地震图算法研究. 地球物理学报. 38(4)    
 *               
 */

#pragma once 

#include <stdio.h>

#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/integral/kernel.h"
#include "grt/integral/k_integ.h"


/**
 * 基于线性插值的Filon积分(5.9.6-11), 在大震中距下对Bessel函数取零阶近似，得
 * \f[
 * J_m(x) \approx \sqrt{\frac{2}{\pi x}} \cos(x - \frac{m \pi}{2} - \frac{\pi}{4})
 * \f]
 * 其中\f$x=kr\f$.
 * 
 * 
 * @param[in,out]  mod1d         `GRT_MODEL1D` 结构体指针
 * @param[in]      k0            前一部分的波数积分结束点k值
 * @param[in]      dk0           前一部分的波数积分间隔
 * @param[in]      filondk       filon积分间隔
 * @param[in]      kmax          波数积分的上限
 * @param[in]      keps          波数积分的收敛条件，要求在某震中距下所有格林函数都收敛
 * @param[in]      nr            震中距数量
 * @param[in]      rs            震中距数组
 *
 * @param[in,out]  K             用于存储积分的结构体
 * 
 * @param[out]    fstats         文件指针，保存不同k值的格林函数积分核函数
 * @param[in]     kerfunc        计算核函数的函数指针
 * 
 * @return  k        积分截至时的波数
 */
real_t grt_linear_filon_integ(
    GRT_MODEL1D *mod1d, real_t k0, real_t dk0, real_t dk, real_t kmax, real_t keps,
    size_t nr, real_t *rs, K_INTEG *K, FILE *fstats, GRT_KernelFunc kerfunc);


