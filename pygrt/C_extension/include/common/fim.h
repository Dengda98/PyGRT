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

#include "common/const.h"
#include "common/model.h"
#include "common/kernel.h"



/**
 * 基于线性插值的Filon积分(5.9.6-11), 在大震中距下对Bessel函数取零阶近似，得
 * \f[
 * J_m(x) \approx \sqrt{\frac{2}{\pi x}} \cos(x - \frac{m \pi}{2} - \frac{\pi}{4})
 * \f]
 * 其中\f$x=kr\f$.
 * 
 * 
 * @param[in]      mod1d         `MODEL1D` 结构体指针
 * @param[in]      k0            前一部分的波数积分结束点k值
 * @param[in]      dk0           前一部分的波数积分间隔
 * @param[in]      filondk       filon积分间隔
 * @param[in]      kmax          波数积分的上限
 * @param[in]      keps          波数积分的收敛条件，要求在某震中距下所有格林函数都收敛
 * @param[in]      omega         复数频率
 * @param[in]      nr            震中距数量
 * @param[in]      rs            震中距数组
 *
 * @param[out]    sum_J          积分值
 * 
 * @param[in]     calc_upar      是否计算位移u的空间导数
 * @param[out]    sum_uiz_J      uiz的积分值
 * @param[out]    sum_uir_J      uir的积分值
 * 
 * @param[out]    fstats         文件指针，保存不同k值的格林函数积分核函数
 * @param[in]     kerfunc        计算核函数的函数指针
 * @param[out]    stats          状态代码，是否有除零错误，非0为异常值
 * 
 * @return  k        积分截至时的波数
 */
MYREAL linear_filon_integ(
    const MODEL1D *mod1d, MYREAL k0, MYREAL dk0, MYREAL filondk, MYREAL kmax, MYREAL keps, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_J[nr][SRC_M_NUM][INTEG_NUM],
    bool calc_upar,
    MYCOMPLEX sum_uiz_J[nr][SRC_M_NUM][INTEG_NUM],
    MYCOMPLEX sum_uir_J[nr][SRC_M_NUM][INTEG_NUM],
    FILE *fstats, KernelFunc kerfunc, MYINT *stats);


