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
 *    + QWV. 存储的是核函数，第一个维度不同震源，不同阶数，第二个维度3代表三类系数qm,wm,vm  
 *    + SUM. 存储的是该dk区间内的积分值，第一个维度不同震源，不同阶数，维度4代表4种类型的F(k,w)Jm(kr)k的类型
 * 
 * 
 * @param[in]     k              波数
 * @param[in]     r              震中距 
 * @param[in]     QWV            不同震源，不同阶数的核函数 \f$ q_m, w_m, v_m \f$
 * @param[in]     calc_uir       是否计算ui_r（位移u对坐标r的偏导）
 * @param[out]    SUM            该dk区间内的积分值
 * 
 */
void int_Pk(
    MYREAL k, MYREAL r, 
    const MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM],
    bool calc_uir,
    MYCOMPLEX SUM[SRC_M_NUM][INTEG_NUM]);




/**
 * 将最终计算好的多个积分值，按照公式(5.6.22)组装成3分量。
 * 
 * @param[in]     sum_J           积分结果
 * @param[out]    tol             Z、R、T分量结果
 */
void merge_Pk(
    const MYCOMPLEX sum_J[SRC_M_NUM][INTEG_NUM], MYCOMPLEX tol[SRC_M_NUM][CHANNEL_NUM]);



/**
 *  和int_Pk函数类似，不过是计算核函数和渐近Bessel函数的乘积 sqrt(k) * F(k,w) * cos ，其中涉及两种数组形状：
 *    + QWV. 存储的是核函数，第一个维度不同震源，不同阶数，第二个维度3代表三类系数qm,wm,vm  
 *    + SUM. 存储的是该dk区间内的积分值，第一个维度不同震源，不同阶数，维度4代表4种类型的F(k,w)Jm(kr)k的类型
 * 
 * 
 * @param[in]     k              波数
 * @param[in]     r              震中距 
 * @param[in]     iscos          是否使用cos函数，否则使用sin函数
 * @param[in]     QWV            不同震源，不同阶数的核函数 \f$ q_m, w_m, v_m \f$
 * @param[in]     calc_uir       是否计算ui_r（位移u对坐标r的偏导）
 * @param[out]    SUM            该dk区间内的积分值
 *  
 */
void int_Pk_filon(
    MYREAL k, MYREAL r, bool iscos,
    const MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM],
    bool calc_uir,
    MYCOMPLEX SUM[SRC_M_NUM][INTEG_NUM]);


/**
 * 对sqrt(k)*F(k,w)进行二次曲线拟合，再计算 (a*k^2 + b*k + c) * cos(kr - (2m+1)/4) 的积分，其中涉及两种数组形状：
 *    + QWV. 存储的是核函数，第一个维度不同震源，不同阶数，第二个维度3代表三类系数qm,wm,vm  
 *    + SUM. 存储的是该三点区间内的积分值，第一个维度不同震源，不同阶数，维度4代表4种类型的F(k,w)Jm(kr)k的类型
 * 
 * @param[in]     k3            三点等距波数
 * @param[in]     r             震中距 
 * @param[in]     QWV3          k3对应的不同震源，不同阶数的核函数 \f$ q_m, w_m, v_m \f$
 * @param[in]     calc_uir      是否计算ui_r（位移u对坐标r的偏导）
 * @param[out]    SUM           该三点区间内的积分值
 * 
 */
void int_Pk_sa_filon(
    const MYREAL k3[3], MYREAL r, 
    const MYCOMPLEX QWV3[3][SRC_M_NUM][QWV_NUM],
    bool calc_uir,
    MYCOMPLEX SUM[SRC_M_NUM][INTEG_NUM]);