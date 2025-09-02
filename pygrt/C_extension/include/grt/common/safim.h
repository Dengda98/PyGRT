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

#include "common/const.h"
#include "common/model.h"
#include "common/kernel.h"



/**
 * 自适应Filon积分, 在大震中距下对Bessel函数取零阶近似，得
 * \f[
 * J_m(x) \approx \sqrt{\frac{2}{\pi x}} \cos(x - \frac{m \pi}{2} - \frac{\pi}{4})
 * \f]
 * 其中\f$x=kr\f$.
 * 
 * 
 * @param[in]      mod1d         `MODEL1D` 结构体指针
 * @param[in]      vmin          最小速度，用于将k区间整体分为两段，在自适应过程中第二段使用更宽松的拟合规则
 * @param[in]      k0            前一部分的波数积分结束点k值
 * @param[in]      dk0           前一部分的波数积分间隔
 * @param[in]      tol           自适应Filon积分的采样精度
 * @param[in]      kmax          波数积分的上限
 * @param[in]      omega         复数频率
 * @param[in]      nr            震中距数量
 * @param[in]      rs            震中距数组
 *
 * @param[out]    sum_J0         积分值
 * 
 * @param[in]     calc_upar      是否计算位移u的空间导数
 * @param[out]    sum_uiz_J0     uiz的积分值
 * @param[out]    sum_uir_J0     uir的积分值
 * 
 * @param[out]    fstats         文件指针，保存不同k值的格林函数积分核函数
 * @param[in]     kerfunc        计算核函数的函数指针
 * @param[out]    stats          状态代码，是否有除零错误，非0为异常值
 * 
 * @return  k        积分截至时的波数
 */
MYREAL sa_filon_integ(
    const MODEL1D *mod1d, MYREAL vmin, MYREAL k0, MYREAL dk0, MYREAL tol, MYREAL kmax, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_J0[nr][SRC_M_NUM][INTEG_NUM],
    bool calc_upar,
    MYCOMPLEX sum_uiz_J0[nr][SRC_M_NUM][INTEG_NUM],
    MYCOMPLEX sum_uir_J0[nr][SRC_M_NUM][INTEG_NUM],
    FILE *fstats, KernelFunc kerfunc, MYINT *stats);


