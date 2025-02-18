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

#include "const.h"
#include "model.h"



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
 * @param  dk        (in)波数积分间隔
 * @param  kmax      (in)波数积分的上限
 * @param  keps      (in)波数积分的收敛条件，要求在某震中距下所有格林函数都收敛
 * @param  omega     (in)复数频率
 * @param  nr        (in)震中距数量
 * @param  rs        (in)震中距数组
 *
 * @param  sum_EXP_J[nr][3][4]  (out)爆炸源
 * @param  sum_VF_J[nr][3][4]   (out)垂直力源
 * @param  sum_HF_J[nr][3][4]   (out)水平力源
 * @param  sum_DC_J[nr][3][4]   (out)双力偶源
 * @param  fstats[nr]           (out)不同震中距的格林函数积分过程文件
 * 
 * @return  k        积分截至时的波数
 */
MYREAL linear_filon_integ(
    const MODEL1D *mod1d, MYREAL dk, MYREAL kmax, MYREAL keps, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_EXP_J[nr][3][4], MYCOMPLEX sum_VF_J[nr][3][4],  
    MYCOMPLEX sum_HF_J[nr][3][4],  MYCOMPLEX sum_DC_J[nr][3][4],  
    FILE *(fstats[nr]));



/**
 *  和int_Pk函数类似，不过是计算核函数和渐近Bessel函数的乘积，其中涉及两种数组形状：
 *    + [3][3]. 存储的是核函数，第一个维度3代表阶数m=0,1,2，第二个维度3代表三类系数qm,wm,vm  
 *    + [3][4]. 存储的是该dk区间内的积分值，维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型
 * 
 * 
 * @param     k     (in)波数  
 * @param     r     (in)震中距 
 * @param    EXP_qwv[3][3]    (in)爆炸源核函数
 * @param    VF_qwv[3][3]     (in)垂直力源核函数
 * @param    HF_qwv[3][3]     (in)水平力源核函数
 * @param    DC_qwv[3][3]     (in)双力偶源核函数 
 * @param    EXP_J[3][4]      (out)爆炸源，该dk区间内的积分值，下同
 * @param    VF_J[3][4]       (out)垂直力源
 * @param    HF_J[3][4]       (out)水平力源
 * @param    DC_J[3][4]       (out)双力偶源
 *  
 */
void int_Pk_filon(
    MYREAL k, MYREAL r, 
    const MYCOMPLEX EXP_qwv[3][3], const MYCOMPLEX VF_qwv[3][3], 
    const MYCOMPLEX HF_qwv[3][3],  const MYCOMPLEX DC_qwv[3][3], 
    MYCOMPLEX EXP_J[3][4], MYCOMPLEX VF_J[3][4], 
    MYCOMPLEX HF_J[3][4],  MYCOMPLEX DC_J[3][4] );

