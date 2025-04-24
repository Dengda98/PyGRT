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
 * 其中\f$x=kr\f$. 结果以三维数组的形式返回，形状为[nr][3][4], 分别代表震中距、阶数(m=0,1,2)
 * 和4种积分类型(p=0,1,2,3)
 * 
 * 
 * @param  mod1d     (in)`MODEL1D` 结构体指针
 * @param  k0        (in)前一部分的波数积分结束点k值
 * @param  dk0       (in)前一部分的波数积分间隔
 * @param  filondk   (in)filon积分间隔
 * @param  kmax      (in)波数积分的上限
 * @param  keps      (in)波数积分的收敛条件，要求在某震中距下所有格林函数都收敛
 * @param  omega     (in)复数频率
 * @param  nr        (in)震中距数量
 * @param  rs        (in)震中距数组
 *
 * @param  sum_EXP_J[nr][3][4]  (out)爆炸源
 * @param  sum_VF_J[nr][3][4]   (out)垂直力源
 * @param  sum_HF_J[nr][3][4]   (out)水平力源
 * @param  sum_DC_J[nr][3][4]   (out)剪切源
 * 
 * @param  calc_upar       (in)是否计算位移u的空间导数
 * @param  sum_EXP_uiz_J[nr][3][4]  (out)爆炸源产生的ui_z(位移u对坐标z的偏导)，下同
 * @param  sum_VF_uiz_J[nr][3][4]   (out)垂直力源
 * @param  sum_HF_uiz_J[nr][3][4]   (out)水平力源
 * @param  sum_DC_uiz_J[nr][3][4]   (out)剪切源
 * @param  sum_EXP_uir_J[nr][3][4]  (out)爆炸源产生的ui_r(位移u对坐标r的偏导)，下同
 * @param  sum_VF_uir_J[nr][3][4]   (out)垂直力源
 * @param  sum_HF_uir_J[nr][3][4]   (out)水平力源
 * @param  sum_DC_uir_J[nr][3][4]   (out)剪切源
 * 
 * @param  fstats               (out)文件指针，保存不同k值的格林函数积分核函数
 * @param  kerfunc              (in)计算核函数的函数指针
 * 
 * @return  k        积分截至时的波数
 */
MYREAL linear_filon_integ(
    const MODEL1D *mod1d, MYREAL k0, MYREAL dk0, MYREAL filondk, MYREAL kmax, MYREAL keps, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_EXP_J[nr][3][4], MYCOMPLEX sum_VF_J[nr][3][4],  
    MYCOMPLEX sum_HF_J[nr][3][4],  MYCOMPLEX sum_DC_J[nr][3][4],  
    bool calc_upar,
    MYCOMPLEX sum_EXP_uiz_J[nr][3][4], MYCOMPLEX sum_VF_uiz_J[nr][3][4],  
    MYCOMPLEX sum_HF_uiz_J[nr][3][4],  MYCOMPLEX sum_DC_uiz_J[nr][3][4],  
    MYCOMPLEX sum_EXP_uir_J[nr][3][4], MYCOMPLEX sum_VF_uir_J[nr][3][4],  
    MYCOMPLEX sum_HF_uir_J[nr][3][4],  MYCOMPLEX sum_DC_uir_J[nr][3][4],  
    FILE *fstats, KernelFunc kerfunc);


