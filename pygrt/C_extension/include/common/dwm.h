/**
 * @file   dwm.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 使用离散波数法求积分，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient 
 *               matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699
 * 
 */

#pragma once 

#include <stdio.h>

#include "common/model.h"
#include "common/kernel.h"


/**
 * 传统的离散波数积分，结果以三维数组的形式返回，形状分别代表震中距、不同震源不同阶数
 * 和4种积分类型(p=0,1,2,3)
 * 
 * @param[in]     mod1d         `MODEL1D` 结构体指针
 * @param[in]     dk            波数积分间隔
 * @param[in]     kmax          波数积分的上限
 * @param[in]     keps          波数积分的收敛条件，要求在某震中距下所有格林函数都收敛，为负数代表不提前判断收敛，按照波数积分上限进行积分
 * @param[in]     omega         复数频率
 * @param[in]     nr            震中距数量
 * @param[in]     rs            震中距数组
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
MYREAL discrete_integ(
    const MODEL1D *mod1d, MYREAL dk, MYREAL kmax, MYREAL keps, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_J[nr][SRC_M_NUM][INTEG_NUM],
    bool calc_upar,
    MYCOMPLEX sum_uiz_J[nr][SRC_M_NUM][INTEG_NUM],
    MYCOMPLEX sum_uir_J[nr][SRC_M_NUM][INTEG_NUM],
    FILE *fstats, KernelFunc kerfunc, MYINT *stats);
