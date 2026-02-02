/**
 * @file   secular.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-07
 * 
 *     计算Rayleigh波和Love波的久期函数
 * 
 */

#pragma once 

#include "grt/common/const.h"
#include "grt/common/model.h"

// /** 
//  * 计算指定层位 iref 处的两个 P-SV 广义矩阵： RD_RL, RU_FR 
//  * 
//  * @param[in]      mod1d        模型结构体指针
//  * @param[in]      omega        圆频率
//  * @param[in]      k            本征波数
//  * @param[in]      iref         当前层位
//  * @param[out]     RD_RL0       RD_RL 矩阵
//  * @param[out]     RU_FR0       RU_FR 矩阵
//  * @param[out]     stats        状态代码，是否有除零错误，非0为异常值
//  */
// void grt_get_mod_iref_RD_RU_Rayl(
//     const GRT_MODEL1D *mod1d, const cplx_t omega, const real_t k, const size_t iref, 
//     cplx_t RD_RL0[2][2], cplx_t RU_FR0[2][2], int *stats);

// /** 计算所有层位的两个  P-SV 广义矩阵： RD_RL, RU_FR  */
// void grt_get_mod_all_RD_RU_Rayl(
//     const GRT_MODEL1D *mod1d, const cplx_t omega, const real_t k,
//     cplx_t RD_RL0[mod1d->n][2][2], cplx_t RU_FR0[mod1d->n][2][2], int *stats);

// /** 
//  * 计算指定层位 iref 处的两个 SH 广义矩阵(标量)： RDL_RL, RUL_FR 
//  * 
//  * @param[in]      mod1d        模型结构体指针
//  * @param[in]      omega        圆频率
//  * @param[in]      k            本征波数
//  * @param[in]      iref         当前层位
//  * @param[out]     RDL_RL0      RDL_RL 值
//  * @param[out]     RUL_FR0      RUL_FR 值
//  * @param[out]     stats        状态代码，是否有除零错误，非0为异常值
//  */
void grt_GRT_matrix_Rayl(GRT_MODEL1D *mod1d, const real_t k, const size_t iref);

/** 计算所有层位的两个 SH 广义矩阵(标量)： RDL_RL, RUL_FR   */
void grt_GRT_matrix_Love(GRT_MODEL1D *mod1d, const real_t k, const size_t iref);

/**
 * 计算 Rayleigh 波的久期函数
 * 
 * @param[in]      mod1d        模型结构体指针
 * @param[in]      cphase       本征相速度
 * @param[in]      secRaylType  Rayl久期函数类型，首层可为0，即 (D21*RD_RL + D22)；其余为1
 * @param[in]      iref         当前层位
 * @param[out]     psec         久期函数值
 * @param[out]     ppot         对应垂直波函数，不需要则设置为NULL
 */
void grt_secular_function_potential_Rayl(
    GRT_MODEL1D *mod1d, const real_t cphase,
    const int secRaylType, const size_t iref, cplx_t *psec, cplx_t *ppot);

/**
 * 计算 Love 波的久期函数
 * 
 * @param[in]      mod1d        模型结构体指针
 * @param[in]      cphase       本征相速度
 * @param[in]      iref         当前层位
 * @param[out]     psec         久期函数值
 * @param[out]     ppot         对应垂直波函数，不需要则设置为NULL
 */
void grt_secular_function_potential_Love(
    GRT_MODEL1D *mod1d, real_t cphase,
    const size_t iref, cplx_t *psec, cplx_t *ppot);

/** 合并 secular_function_potential_Rayl(Love) */
void grt_secular_function_potential(
    GRT_MODEL1D *mod1d, real_t cphase,
    const int secRaylType, const size_t iref, const DISPER_TYPE wtype, cplx_t *psec, cplx_t *ppot);

/** 只求久期函数值，即在 secular_function_potential 中设置 ppot=NULL */
void grt_secular_function(
    GRT_MODEL1D *mod1d, real_t cphase,
    const int secRaylType, const size_t iref, const DISPER_TYPE wtype, cplx_t *psec);

/** 
 * 黄金分割法确定久期函数零点 
 * 
 * @param[in]      mod1d        模型结构体指针
 * @param[in]      w            圆频率
 * @param[in]      c1           左区间相速度
 * @param[in]      c2           右区间相速度
 * @param[in]      secRaylType  Rayl久期函数类型
 * @param[in]      iref         久期函数层位
 * @param[in]      wtype        频散类型
 * @param[in]      rtol         久期函数零点附近的幅值阈值
 * @param[out]     root_sec     零点的久期函数值
 * @param[out]     root_c       零点的相速度
 * 
 */
size_t grt_goldensection_search_root(
    GRT_MODEL1D *mod1d,  const real_t c1, const real_t c2, 
    const int secRaylType, const size_t iref, const DISPER_TYPE wtype,
    const real_t rtol, cplx_t *root_sec, real_t *root_c);