/**
 * @file   modal_util.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 *     频散相关的辅助函数
 *                   
 */

#pragma once

#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/modal/modal_def.h"


/** 
 * 输出相速度频散结果
 * 
 * @param[in]    filepath         输出路径
 * @param[in]    full_command     计算频散的完整命令
 * @param[in]    modelname        模型文件名（不含路径）
 * @param[in]    mod1d            模型，使用其中未插入虚拟层的原始 ``modarr``；写入 NC 时仅保存前四列
 * @param[in]    eigmet           频散数据结构体指针
 */
void grt_output_cdisp(
    const char *filepath, const char *full_command, const char *modelname,
    const MODEL1D *mod1d, EIGENV_INFO *eigmet);

/** 
 * 输出群速度频散结果 
 * 
 * @param[in]    filepath         输出路径
 * @param[in]    eigfnmet         本征函数数据结构体指针
 */
void grt_output_udisp(const char *filepath, EIGENFN_INFO *eigfnmet);

/** 
 * 输出本征函数结果
 * 
 * @param[in]    filepath         输出路径
 * @param[in]    ncols            4 for Rayleigh, 2 for Love
 * @param[in]    eigfnmet         本征函数数据结构体指针
 */
void grt_output_eigenfns(const char *filepath, const int ncols, EIGENFN_INFO *eigfnmet);

/** 
 * 输出能量积分结果
 * 
 * @param[in]    filepath         输出路径
 * @param[in]    eigfnmet         本征函数数据结构体指针
 */
void grt_output_energy_integrals(const char *filepath, EIGENFN_INFO *eigfnmet);

/** 
 * 输出相速度/群速度敏感核
 * 
 * @param[in]    filepath         输出路径
 * @param[in]    char_uc          'c' for phase velocity, 'u' for group velocity
 * @param[in]    eigfnmet         本征函数数据结构体指针
 * 
 */
void grt_output_sensitivity(const char *filepath, const char *char_uc, EIGENFN_INFO *eigfnmet);

/** 
 * 根据计算好的相速度、群速度和相速度敏感核，计算群速度敏感核
 * 
 * 根据群速度与相速度的关系： U = c^2 / ( c - w*dc/dw )
 * 推导得到敏感核（其中 X 表示 alpha/beta/rho ）
 *     dU/dX = U/c * (2 - U/c) * dc/dX + (U/c)^2 * w * d(dc/dX)/dw
 * 后一项存在对相速度敏感核在频率上做差分，
 * 因此输入的频率间隔大小会直接影响群速度敏感核精度
 * 
 * @param[in,out]    eigfnmet      本征函数数据结构体指针
 */
void grt_group_sensitivity(EIGENFN_INFO *eigfnmet);



/** 
 * 读取相/群速度频散结果
 * 
 * @param[in]     filepath         读入频散文件路径
 * @param[out]    eigmet           频散数据结构体指针
 * @param[out]    pt_modelname     模型文件名，可为 NULL
 *                                 群速度文件不含此项，非 NULL 时置为 NULL
 * @param[out]    pt_mod1d         由内嵌模型构建的 `MODEL1D`，可为 NULL；非 NULL 时由调用方释放
 *                                 群速度文件不含内嵌模型，此时若非 NULL 则报错
 */
void grt_read_dispersion(
    const char *filepath, EIGENV_INFO *eigmet, char **pt_modelname, MODEL1D **pt_mod1d);
