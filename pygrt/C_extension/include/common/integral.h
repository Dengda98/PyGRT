/**
 * @file   integral.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 *     将被积函数的逐点值累加成积分值
 *                   
 */

#pragma once 

#include "common/const.h"

/**
 * 计算核函数和Bessel函数的乘积，相当于计算了一个小积分区间内的值。参数中涉及两种数组形状：
 *    + [3][3]. 存储的是核函数，第一个维度3代表阶数m=0,1,2，第二个维度3代表三类系数qm,wm,vm  
 *    + [3][4]. 存储的是该dk区间内的积分值，维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型
 * 
 * 
 * @param     k     (in)波数
 * @param     r     (in)震中距 
 * @param    EXP_qwv[3][3]    (in)爆炸源核函数
 * @param    VF_qwv[3][3]     (in)垂直力源核函数
 * @param    HF_qwv[3][3]     (in)水平力源核函数
 * @param    DC_qwv[3][3]     (in)剪切源核函数
 * @param    calc_uir         (in)是否计算ui_r（位移u对坐标r的偏导）
 * @param    EXP_J[3][4]      (out)爆炸源，该dk区间内的积分值，下同
 * @param    VF_J[3][4]       (out)垂直力源
 * @param    HF_J[3][4]       (out)水平力源
 * @param    DC_J[3][4]       (out)剪切源
 * 
 */
void int_Pk(
    MYREAL k, MYREAL r, 
    const MYCOMPLEX EXP_qwv[3][3], const MYCOMPLEX VF_qwv[3][3], 
    const MYCOMPLEX HF_qwv[3][3],  const MYCOMPLEX DC_qwv[3][3], 
    bool calc_uir,
    MYCOMPLEX EXP_J[3][4], MYCOMPLEX VF_J[3][4], 
    MYCOMPLEX HF_J[3][4],  MYCOMPLEX DC_J[3][4] );




/**
 * 将最终计算好的多个积分值，按照公式(5.6.22)组装成3分量。数组形状[3][4]，\
 * 存储的是最终的积分值，维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型
 * 
 * @param    sum_EXP_J[3][4]      (in)爆炸源，最终的积分值，下同
 * @param    sum_VF_J[3][4]       (in)垂直力源
 * @param    sum_HF_J[3][4]       (in)水平力源
 * @param    sum_DC_J[3][4]       (in)剪切源
 * @param    tol_EXP[2]           (out)爆炸源的Z、R分量频谱结果
 * @param    tol_VF[2]            (out)垂直力源的Z、R分量频谱结果
 * @param    tol_HF[3]            (out)水平力源的Z、R、T分量频谱结果
 * @param    tol_DD[2]            (out)45度倾滑的Z、R分量频谱结果
 * @param    tol_DS[3]            (out)90度倾滑的Z、R、T分量频谱结果
 * @param    tol_SS[3]            (out)90度走滑的Z、R、T分量频谱结果
 */
void merge_Pk(
    const MYCOMPLEX sum_EXP_J[3][4], const MYCOMPLEX sum_VF_J[3][4], 
    const MYCOMPLEX sum_HF_J[3][4],  const MYCOMPLEX sum_DC_J[3][4], 
    MYCOMPLEX tol_EXP[2], MYCOMPLEX tol_VF[2], MYCOMPLEX tol_HF[3],
    MYCOMPLEX tol_DD[2],  MYCOMPLEX tol_DS[3], MYCOMPLEX tol_SS[3]);



/**
 *  和int_Pk函数类似，不过是计算核函数和渐近Bessel函数的乘积 sqrt(k) * F(k,w) * cos ，其中涉及两种数组形状：
 *    + [3][3]. 存储的是核函数，第一个维度3代表阶数m=0,1,2，第二个维度3代表三类系数qm,wm,vm  
 *    + [3][4]. 存储的是该dk区间内的积分值，维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型
 * 
 * 
 * @param     k     (in)波数  
 * @param     r     (in)震中距 
 * @param     iscos           (in)计算sin函数
 * @param    EXP_qwv[3][3]    (in)爆炸源核函数
 * @param    VF_qwv[3][3]     (in)垂直力源核函数
 * @param    HF_qwv[3][3]     (in)水平力源核函数
 * @param    DC_qwv[3][3]     (in)剪切源核函数 
 * @param    calc_uir         (in)是否计算ui_r（位移u对坐标r的偏导）
 * @param    EXP_J[3][4]      (out)爆炸源，该dk区间内的积分值，下同
 * @param    VF_J[3][4]       (out)垂直力源
 * @param    HF_J[3][4]       (out)水平力源
 * @param    DC_J[3][4]       (out)剪切源
 *  
 */
void int_Pk_filon(
    MYREAL k, MYREAL r, bool iscos,
    const MYCOMPLEX EXP_qwv[3][3], const MYCOMPLEX VF_qwv[3][3], 
    const MYCOMPLEX HF_qwv[3][3],  const MYCOMPLEX DC_qwv[3][3], 
    bool calc_uir, 
    MYCOMPLEX EXP_J[3][4], MYCOMPLEX VF_J[3][4], 
    MYCOMPLEX HF_J[3][4],  MYCOMPLEX DC_J[3][4] );
