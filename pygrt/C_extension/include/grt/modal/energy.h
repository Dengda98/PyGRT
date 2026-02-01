/**
 * @file   energy.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-07
 * 
 * 根据每层界面垂直波函数，计算能量积分和群速度
 * 
 */

#pragma once

#include "grt/common/model.h"
#include "grt/common/const.h"


/**
 * 计算 Rayleigh 波的能量积分和敏感核
 * 
 * @param[in]      mod1d         模型结构体指针
 * @param[in]      omega         圆频率
 * @param[in]      eigenK        本征波数
 * @param[in]      mod_potRaylLove_Down 每层 z_j+ 的垂直波函数
 * @param[in]      mod_potRaylLove_Up   每层 z_j- 的垂直波函数
 * @param[out]     EgyInt        能量积分
 * @param[in]      cpar_nz       计算敏感核划分的薄层数量
 * @param[in]      cpar_zs       薄层深度
 * @param[in]      cpar_z_irefs  薄层所在层位
 * @param[out]     phase_K       敏感核
 * 
 */
void grt_energy_integrals_Rayl(
    const GRT_MODEL1D *mod1d, const real_t omega, const real_t eigenK,
    const cplx_t (*mod_potRaylLove_Down)[GRT_RAYL_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_RAYL_DIM], cplx_t *EgyInt, 
    const size_t cpar_nz, const real_t *cpar_zs, const size_t *cpar_z_irefs, cplx_t phase_K[cpar_nz][GRT_SNSTVTY_MAX]);


/**
 * 计算 Love 波的能量积分和敏感核
 * 
 * @param[in]      mod1d         模型结构体指针
 * @param[in]      omega         圆频率
 * @param[in]      eigenK        本征波数
 * @param[in]      mod_potRaylLove_Down 每层 z_j+ 的垂直波函数
 * @param[in]      mod_potRaylLove_Up   每层 z_j- 的垂直波函数
 * @param[out]     EgyInt        能量积分
 * @param[in]      cpar_nz       计算敏感核划分的薄层数量
 * @param[in]      cpar_zs       薄层深度
 * @param[in]      cpar_z_irefs  薄层所在层位
 * @param[out]     phase_K       敏感核
 * 
 */
void grt_energy_integrals_Love(
    const GRT_MODEL1D *mod1d, const real_t omega, const real_t eigenK,
    const cplx_t (*mod_potRaylLove_Down)[GRT_LOVE_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_LOVE_DIM], cplx_t *EgyInt, 
    const size_t cpar_nz, const real_t *cpar_zs, const size_t *cpar_z_irefs, cplx_t phase_K[cpar_nz][GRT_SNSTVTY_MAX]);


/** 计算能量积分, 合并 energy_integrals_Rayl(Love) */
void grt_energy_integrals(
    const GRT_MODEL1D *mod1d, const real_t omega, const real_t eigenK, const size_t ncols, 
    const cplx_t (*mod_potRaylLove_Down)[ncols], const cplx_t (*mod_potRaylLove_Up)[ncols], 
    cplx_t *EgyInt, const size_t cpar_nz, const real_t *cpar_zs, const size_t *cpar_z_irefs, 
    cplx_t phase_K[cpar_nz][GRT_SNSTVTY_MAX], const bool isRayl);


/** 
 * 根据能量积分计算 Rayleight 波和 Love 波的群速度 
 * 
 * @param[in]       omega       圆频率
 * @param[in]       eigenK      本征波数
 * @param[in]       EgyInt      能量积分
 * @param[in]       isRayl      Rayleigh or Love
 * 
 * @return     群速度
 * 
 */
real_t grt_get_group_velocity(const real_t omega, const real_t eigenK, const cplx_t EgyInt[GRT_EGYINTS_MAX], const bool isRayl);