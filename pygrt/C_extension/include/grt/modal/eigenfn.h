/**
 * @file   eigenfn.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-07
 * 
 * 使用 广义反射透射系数矩阵 将某层的垂直波函数传播到整个模型的物理层上下两侧，
 * 再计算任意深度的本征函数
 * 
 */

#pragma once

#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/modal/modal_def.h"



/**
 * 已知某个界面 z_j+ 上的垂直波函数，计算模型每个分界面 z_j+ 和 z_j- 的垂直波函数 P-SV
 * 
 * @param[in]      mod1d                模型结构体指针
 * @param[in]      omega                圆频率
 * @param[in]      eigenK               本征波数
 * @param[out]     mod_potRaylLove_Down 每层 z_j+ 的垂直波函数
 * @param[out]     mod_potRaylLove_Up   每层 z_j- 的垂直波函数
 * @param[in]      iy                   已知垂直波函数的层位
 * @param[in]      potRaylLove          已知垂直波函数
 * @param[out]     stats                状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_get_mod_potential_Up_Down_Rayl(
    GRT_MODEL1D *mod1d, const size_t iref, const cplx_t potRaylLove[GRT_RAYL_DIM],
    cplx_t (*mod_potRaylLove_Down)[GRT_RAYL_DIM], cplx_t (*mod_potRaylLove_Up)[GRT_RAYL_DIM]);

/** 已知某个界面 z_j+ 上的垂直波函数，计算模型每个分界面 z_j+ 和 z_j- 的垂直波函数 SH, 参数见 get_mod_potential_Up_Down_Rayl */
void grt_get_mod_potential_Up_Down_Love(
    GRT_MODEL1D *mod1d, const size_t iref, const cplx_t potRaylLove[GRT_LOVE_DIM],
    cplx_t (*mod_potRaylLove_Down)[GRT_LOVE_DIM], cplx_t (*mod_potRaylLove_Up)[GRT_LOVE_DIM]);

/** 合并 get_mod_potential_Up_Down_Rayl(Love) */
void grt_get_mod_potential_Up_Down(
    GRT_MODEL1D *mod1d, const DISPER_TYPE wtype, const size_t ncols, const size_t iref, const cplx_t potRaylLove[ncols], 
    cplx_t (*mod_potRaylLove_Down)[ncols], cplx_t (*mod_potRaylLove_Up)[ncols]);


/**
 * 计算多个深度处的本征函数 P-SV
 * 
 * @param[in]      mod1d                模型结构体指针
 * @param[in]      omega                圆频率
 * @param[in]      eigenK               本征波数
 * @param[in]      mod_potRaylLove_Down 每层 z_j+ 的垂直波函数
 * @param[in]      mod_potRaylLove_Up   每层 z_j- 的垂直波函数
 * @param[in]      nz                   沿深度采样数量
 * @param[in]      zsamps               有序深度数组
 * @param[in]      z_irefs              每个采样点所在层位
 * @param[out]     eigenfns             每个采样点本征函数
 * 
 */
void grt_get_eigenfn_depths_Rayl(
    const GRT_MODEL1D *mod1d, const cplx_t (*mod_potRaylLove_Down)[GRT_RAYL_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_RAYL_DIM],
    const EIGENFN_INFO *eigfnmet, EIGENFN *eigfn);

/** 计算多个深度处的本征函数 SH, 参数见 get_eigenfn_depths_Rayl */
void grt_get_eigenfn_depths_Love(
    const GRT_MODEL1D *mod1d, const cplx_t (*mod_potRaylLove_Down)[GRT_LOVE_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_LOVE_DIM],
    const EIGENFN_INFO *eigfnmet, EIGENFN *eigfn);

/** 合并 get_eigenfn_depths_Rayl(Love)  */
void grt_get_eigenfn_depths(
    const GRT_MODEL1D *mod1d, const DISPER_TYPE wtype, const size_t ncols, 
    const cplx_t (*mod_potRaylLove_Down)[ncols], const cplx_t (*mod_potRaylLove_Up)[ncols],
    const EIGENFN_INFO *eigfnmet, EIGENFN *eigfn);