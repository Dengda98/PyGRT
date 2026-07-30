/**
 * @file   k_integ.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 *     将被积函数的逐点值累加成积分值
 *                   
 */

#pragma once 

#include <stdbool.h>

#include "grt/common/const.h"

typedef struct {
    bool calc_upar;

    // 不同震源不同阶数的核函数 F(k, w) 
    cplxChnlGrid QWV;
    cplxChnlGrid QWVz;
    cplxChnlGrid QWV_raw;  ///< 不受 DCM 影响，仅为计算出的核函数
    cplxChnlGrid QWVz_raw;

    // 最大波数时的核函数，用于 DCM
    bool applyDCM;
    real_t kmax;
    cplxChnlGrid QWV_kmax;
    cplxChnlGrid QWVz_kmax;

    size_t nr;
    cplxIntegGrid SUM;     ///< 被积函数，用于临时存储
    cplxIntegGrid *sumJ;   ///< 积分值
    cplxIntegGrid *sumJz;  ///< z偏导的积分值
    cplxIntegGrid *sumJr;  ///< r偏导的积分值
} K_INTEG;

/**
 * 初始化 K_INTEG
 * 
 * @param[in]    calc_upar     是否计算空间导数，这将决定是否申请对应变量的内存
 * @param[in]    nr            震中距个数
 * 
 * @return      K_INTEG 结构体指针
 */
K_INTEG * grt_init_K_INTEG(const bool calc_upar, const size_t nr);

/**
 * 复制 K_INTEG 
 * 
 * @param[in]   K     源 K_INTEG 结构体指针
 * @return      复制的K_INTEG 结构体指针
 */
K_INTEG * grt_copy_K_INTEG(const K_INTEG *K);

/**
 * 释放 K_INTEG 指针内存
 * 
 * @param[in,out]   K     K_INTEG 结构体指针
 */
void grt_free_K_INTEG(K_INTEG *K);

/** 
 * 波数积分过程对 DCM 对原积分的分解
 * 
 * @param[in]     kmax         DCM 校正所用的最大波数
 * @param[in]     k            当前积分波数
 * @param[in,out] QWV          当前波数处的核函数；函数返回时已减去 DCM 分解项
 * @param[in]     QWV_kmax     最大波数处的核函数，用于构造 DCM 分解项
 * @param[in]     calc_upar    是否同时校正核函数的 z 方向空间导数
 * @param[in,out] QWVz         当前波数处核函数的 z 方向导数；仅当 calc_upar 为 true 时修改
 * @param[in]     QWVz_kmax    最大波数处核函数的 z 方向导数，用于构造 DCM 分解项
 */
void grt_int_dcm_revision(real_t kmax, real_t k, cplxChnlGrid QWV, cplxChnlGrid QWV_kmax, bool calc_upar, cplxChnlGrid QWVz, cplxChnlGrid QWVz_kmax);

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
void grt_int_Pk(real_t k, real_t r, const cplxChnlGrid QWV, bool calc_uir, cplxIntegGrid SUM);




/**
 * 将最终计算好的多个积分值，按照公式(5.6.22)组装成3分量。
 * 
 * @param[in]     sumJ            积分结果
 * @param[out]    tol             Z、R、T分量结果
 */
void grt_merge_Pk(const cplxIntegGrid sumJ, cplxChnlGrid tol);



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
void grt_int_Pk_filon(real_t k, real_t r, bool iscos, const cplxChnlGrid QWV, bool calc_uir, cplxIntegGrid SUM);


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
void grt_int_Pk_sa_filon(const real_t k3[3], real_t r, const cplxChnlGrid QWV3[3], bool calc_uir, cplxIntegGrid SUM);