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
#include "grt/modal/modal_def.h"


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
    const GRT_MODEL1D *mod1d, const cplx_t (*mod_potRaylLove_Down)[GRT_RAYL_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_RAYL_DIM], 
    const EIGENFN_INFO *eigfnmet, EIGENFN *eigfn);


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
    const GRT_MODEL1D *mod1d, const cplx_t (*mod_potRaylLove_Down)[GRT_LOVE_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_LOVE_DIM], 
    const EIGENFN_INFO *eigfnmet, EIGENFN *eigfn);


/** 计算能量积分, 合并 energy_integrals_Rayl(Love) */
void grt_energy_integrals(
    const GRT_MODEL1D *mod1d, const DISPER_TYPE wtype, const size_t ncols, 
    const cplx_t (*mod_potRaylLove_Down)[ncols], const cplx_t (*mod_potRaylLove_Up)[ncols], 
    const EIGENFN_INFO *eigfnmet, EIGENFN *eigfn);

