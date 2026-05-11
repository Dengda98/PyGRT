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

/** 
 * 计算所有层位的两个 P-SV 广义矩阵： RD_RL, RU_FR 
 * 
 * @param[in]      mstat        模型结构体指针
 * @param[in]      k            本征波数
 * @param[out]     Mall_RL      所有层的 RD_RL 矩阵
 * @param[out]     Mall_FR      所有层的 RU_FR 矩阵
 */
void grt_GRT_matrix_allLayer_Rayl(MODEL1D_STATE *mstat, const real_t k, RT_MATRIX *Mall_RL, RT_MATRIX *Mall_FR);

/** 计算所有层位的两个  P-SV 广义矩阵： RD_RL, RU_FR  */
void grt_GRT_matrix_allLayer_Love(MODEL1D_STATE *mstat, const real_t k, RT_MATRIX *Mall_RL, RT_MATRIX *Mall_FR);

/** 
 * 计算指定层位 iref 处的两个 SH 广义矩阵(标量)： RDL_RL, RUL_FR 
 * 
 * @param[in]      mstat        模型结构体指针
 * @param[in]      k            本征波数
 * @param[in]      iref         当前层位
 */
void grt_GRT_matrix_Rayl(MODEL1D_STATE *mstat, const real_t k, const size_t iref);

/** 计算所有层位的两个 SH 广义矩阵(标量)： RDL_RL, RUL_FR   */
void grt_GRT_matrix_Love(MODEL1D_STATE *mstat, const real_t k, const size_t iref);



/**
 * 以 iref 层为参考层定义久期函数时，相速度搜索的一个比较合适的起始位置
 * 
 * @param[in]     mstat        模型结构体指针
 * @param[in]     iref         当前层位
 * @param[in]     wtype        频散类型
 * @return        相速度搜索起始点
 */
real_t grt_secular_function_cbegin(const MODEL1D_STATE *mstat, const size_t iref, DISPER_TYPE wtype);


/**
 * 计算 Rayleigh 波的久期函数
 * 
 * @param[in]      mstat        模型结构体指针
 * @param[in]      cphase       本征相速度
 * @param[in]      iref         当前层位
 * @param[out]     psec         久期函数值
 * @param[out]     ppot         对应z+处垂直波函数，不需要则设置为NULL
 * @param[out]     ppotUp       对应z-处垂直波函数，不需要则设置为NULL【这个仅在液固界面才会用到】
 */
void grt_secular_function_potential_Rayl(
    MODEL1D_STATE *mstat, const real_t cphase, const size_t iref, cplx_t *psec, cplx_t ppot[GRT_RAYL_DIM], cplx_t ppotUp[GRT_RAYL_DIM]);

/**
 * 计算 Love 波的久期函数
 * 
 * @param[in]      mstat        模型结构体指针
 * @param[in]      cphase       本征相速度
 * @param[in]      iref         当前层位
 * @param[out]     psec         久期函数值
 * @param[out]     ppot         对应垂直波函数，不需要则设置为NULL
 */
void grt_secular_function_potential_Love(
    MODEL1D_STATE *mstat, real_t cphase, const size_t iref, cplx_t *psec, cplx_t ppot[GRT_LOVE_DIM]);

/** 合并 secular_function_potential_Rayl(Love) */
void grt_secular_function_potential(
    MODEL1D_STATE *mstat, real_t cphase, const size_t iref, const DISPER_TYPE wtype, cplx_t *psec, cplx_t *ppot, cplx_t *ppotUp);

/** 只求久期函数值，即在 secular_function_potential 中设置 ppot=NULL, ppotUp=NULL */
void grt_secular_function(
    MODEL1D_STATE *mstat, real_t cphase, const size_t iref, const DISPER_TYPE wtype, cplx_t *psec);
