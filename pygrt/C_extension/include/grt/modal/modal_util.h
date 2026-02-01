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
void grt_output_cdisp(const char *filepath, const char *full_command, const char *modelpath, EIGENV_METHOD *eigmet);

// /** 输出群速度频散结果 */
// void grt_output_gdisp(
//     const char *filepath, 
//     const size_t nfreq, const real_t *freqs, const bool io_period, const bool isRayl,
//     const size_t *modes,
//     const real_t *const cdisp[nfreq], const real_t *const gdisp[nfreq], const size_t cdisp_num[nfreq]);

// /** 读取频散文件 */
// void grt_read_mod1d_cdisp(
//     const char *command, const char *s_filepath, char **pt_modelpath, const bool def_freq_period, 
//     size_t *pt_nfreq, real_t **pt_freqs, bool *pt_io_period, bool *pt_readin_period,
//     size_t *pt_nmode, size_t **pt_modes, 
//     real_t ***pt_cdisp, size_t ***pt_cdisp_iref, size_t **pt_cdisp_num,
//     bool *pt_isRayl, GRT_MODEL1D **pt_mod1d, real_t depsrc, real_t deprcv);

// /* 输出本征函数结果 */
// void grt_output_eigenfns(
//     const char *filepath, 
//     const size_t nfreq, const real_t *freqs, const bool io_period, const bool isRayl,
//     const size_t *modes, const size_t ncols, 
//     const real_t *const cdisp[nfreq], const size_t cdisp_num[nfreq], 
//     const size_t nz, const real_t zsamps[nz],
//     const cplx_t (*const *eigenfns)[nz][ncols]);

// /* 输出能量积分结果 */
// void grt_output_energy_integrals(
//     const char *filepath, 
//     const size_t nfreq, const real_t *freqs, const bool io_period, const bool isRayl,
//     const size_t *modes, 
//     const real_t *const cdisp[nfreq], const size_t cdisp_num[nfreq], 
//     const cplx_t (* const *EgyInt)[GRT_EGYINTS_MAX]);

// /* 输出敏感核 */
// void grt_output_sensitivity(
//     const char *filepath, 
//     const size_t nfreq, const real_t *freqs, const bool io_period, const bool isRayl,
//     const size_t *modes,
//     const real_t *const cgdisp[nfreq], const size_t cdisp_num[nfreq],
//     const char *UC, const size_t cpar_nz, const real_t cpar_zs[cpar_nz], cplx_t (**sense_K)[cpar_nz][GRT_SNSTVTY_MAX]);

// void grt_eigenfn_egy_gdisp_phaseK_util(
//     GRT_MODEL1D *mod1d, const size_t nfreq, const real_t freqs[nfreq], 
//     const size_t ncols, const real_t *const cdisp[nfreq], const size_t *const cdisp_iref[nfreq], const size_t cdisp_num[nfreq], 
//     const bool calc_egfn, const size_t nz, const real_t zsamps[nz], const size_t z_irefs[nz], cplx_t (**eigenfns)[nz][ncols], 
//     const bool calc_egyint, cplx_t (**EgyInt)[GRT_EGYINTS_MAX], 
//     const size_t cpar_nz, const real_t *cpar_zs, const size_t *cpar_z_irefs, cplx_t (**phase_K)[cpar_nz][GRT_SNSTVTY_MAX], 
//     real_t **gdisp, const bool isRayl);

// void grt_group_sensitivity(
//     const size_t nfreq, const real_t freqs[nfreq],
//     const real_t *const cdisp[nfreq], const size_t cdisp_num[nfreq], 
//     const size_t cpar_nz, cplx_t (**phase_K)[cpar_nz][GRT_SNSTVTY_MAX], real_t **gdisp, 
//     cplx_t (**group_K)[cpar_nz][GRT_SNSTVTY_MAX]);