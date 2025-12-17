/**
 * @file   safim.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-4-27
 * 
 * 以下代码实现的是自适应Filon积分，参考：
 * 
 *         Chen, X., and H. Zhang (2001). An Efficient Method for Computing Green’s Functions for a 
 *         Layered Half-Space at Large Epicentral Distances, Bulletin of the Seismological Society of America 91, 
 *         no. 4, 858–869, doi: 10.1785/0120000113.
 *               
 */

#pragma once 

#include <stdio.h>

#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/integral/kernel.h"
#include "grt/integral/k_integ.h"



/**
 * 自适应Filon积分, 在大震中距下对Bessel函数取零阶近似，得
 * \f[
 * J_m(x) \approx \sqrt{\frac{2}{\pi x}} \cos(x - \frac{m \pi}{2} - \frac{\pi}{4})
 * \f]
 * 其中\f$x=kr\f$.
 * 
 * 
 * @param[in,out]  mod1d         `GRT_MODEL1D` 结构体指针
 * @param[in]      k0            前一部分的波数积分结束点k值
 * @param[in]      dk0           前一部分的波数积分间隔
 * @param[in]      tol           自适应Filon积分的采样精度
 * @param[in]      kmax          波数积分的上限
 * @param[in]      kref          将k区间整体分为[dk0, kref]和[kref, kmax]，后一段使用更宽松的拟合规则
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
real_t grt_sa_filon_integ(
    GRT_MODEL1D *mod1d, real_t k0, real_t dk0, real_t tol, real_t kmax, real_t kref,
    size_t nr, real_t *rs, K_INTEG *K, FILE *fstats, GRT_KernelFunc kerfunc);


