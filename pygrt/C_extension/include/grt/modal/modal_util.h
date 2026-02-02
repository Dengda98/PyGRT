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


/** 输出相速度频散结果 */
void grt_output_cdisp(const char *filepath, const char *full_command, const char *modelpath, EIGENV_INFO *eigmet);

/** 输出群速度频散结果 */
void grt_output_gdisp(const char *filepath, EIGENFN_INFO *eigfnmet);

/** 读取频散文件 */
void grt_read_mod1d_cdisp(const char *s_filepath, EIGENV_INFO *eigmet, const real_t depsrc, const real_t deprcv, GRT_MODEL1D *mod1d);

/* 输出本征函数结果 */
void grt_output_eigenfns(const char *filepath, const int ncols, EIGENFN_INFO *eigfnmet);

/* 输出能量积分结果 */
void grt_output_energy_integrals(const char *filepath, EIGENFN_INFO *eigfnmet);

/* 输出敏感核 */
void grt_output_sensitivity(const char *filepath, const char *UC, EIGENFN_INFO *eigfnmet);

/** 根据计算好的相速度、群速度和相速度敏感核，计算群速度敏感核 */
void grt_group_sensitivity(EIGENFN_INFO *eigfnmet);