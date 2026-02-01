/**
 * @file   modal_util.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 *     频散相关的辅助函数
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/common/search.h"
#include "grt/common/mynetcdf.h"

#include "grt/modal/modal_util.h"
#include "grt/modal/secular.h"
#include "grt/modal/eigenfn.h"
#include "grt/modal/energy.h"

#include "grt/common/checkerror.h"

void grt_output_cdisp(const char *filepath, const char *full_command, const char *modelpath, EIGENV_INFO *eigmet)
{
    int ncid, f_dimid, n_dimid;
    const int ndims = 2;
    int dimids[ndims];
    int f_varid, n_varid;
    int c_varid, ciref_varid, cnum_varid;

    // 创建 NC 文件
    NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

    // 写入生成频散的命令
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "command", strlen(full_command), full_command));
    // 写入模型路径
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "model", strlen(modelpath), modelpath));

    EIGENV *eigv = eigmet->eigv;

    // 最大零点数
    size_t max_nroots = 0;
    for(size_t iw=0; iw<eigmet->nf; ++iw){
        max_nroots = GRT_MAX(max_nroots, eigv[iw].n);
    }

    // 定义维度
    NC_CHECK(nc_def_dim(ncid, "freq", eigmet->nf, &f_dimid));
    NC_CHECK(nc_def_dim(ncid, "mode", max_nroots, &n_dimid));

    // 定义变量维度数组
    dimids[0] = f_dimid;
    dimids[1] = n_dimid;

    // 定义维度数组
    NC_CHECK(nc_def_var(ncid, "freq", NC_REAL, 1, &f_dimid, &f_varid));
    NC_CHECK(nc_def_var(ncid, "mode", NC_INT,  1, &n_dimid, &n_varid));

    // 定义频散数组
    NC_CHECK(nc_def_var(ncid, "cnum", NC_INT, 1, &f_dimid, &cnum_varid));
    NC_CHECK(nc_def_var(ncid, "c", NC_REAL, ndims, dimids, &c_varid));
    NC_CHECK(nc_def_var(ncid, "ciref", NC_INT, ndims, dimids, &ciref_varid));

    // 填充值
    const real_t fill_value_real_t = 0.0;
    const int fill_value_MYINT = -1;
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, c_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));
    NC_CHECK(NC_FUNC_INT(nc_put_att)  (ncid, ciref_varid, "_FillValue", NC_INT, 1, &fill_value_MYINT));

    // 其它属性
    // Rayleigh or Love
    int isRayl_int = (int)(eigmet->wtype==GRT_DISPERSION_RAYL);
    NC_CHECK(NC_FUNC_INT(nc_put_att) (ncid, NC_GLOBAL, "isRayl", NC_INT, 1, &isRayl_int));
    
    // 结束定义模式
    NC_CHECK(nc_enddef(ncid));

    // 写入数据
    NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, f_varid, eigmet->freqs));
    {
        int *modes = (int *)calloc(max_nroots, sizeof(int));
        for(size_t i=0; i<max_nroots; ++i){
            modes[i] = i;
        } 
        NC_CHECK(NC_FUNC_INT(nc_put_var) (ncid, n_varid, modes));
        GRT_SAFE_FREE_PTR(modes);
    }
    {
        // size_t 数组转为 int 数组保存
        int *tmp = (int *)calloc(eigmet->nf, sizeof(int));
        for(size_t i = 0; i < eigmet->nf; ++i){
            tmp[i] = eigv[i].n;
        }
        NC_CHECK(NC_FUNC_INT(nc_put_var) (ncid, cnum_varid, tmp));
        GRT_SAFE_FREE_PTR(tmp);
    }

    // 分批次写入数据
    size_t startp[ndims];
    size_t countp[ndims];
    for(size_t iw = 0; iw < eigmet->nf; ++iw){
        startp[0] = iw;
        startp[1] = 0;
        countp[0] = 1;
        countp[1] = eigv[iw].n;

        NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, c_varid, startp, countp, eigv[iw].c_roots));
        {
            // size_t 数组转为 int 数组保存
            int *tmp = (int *)calloc(eigv[iw].n, sizeof(int));
            for(size_t i = 0; i < eigv[iw].n; ++i){
                tmp[i] = eigv[iw].c_roots_iref[i];
            }
            NC_CHECK(NC_FUNC_INT(nc_put_vara) (ncid, ciref_varid, startp, countp, tmp));
            GRT_SAFE_FREE_PTR(tmp);
        }
    }
    
    // 关闭文件
    NC_CHECK(nc_close(ncid));

}



// void grt_read_mod1d_cdisp(
//     const char *command, const char *s_filepath, char **pt_modelpath, const bool def_freq_period, 
//     MYINT *pt_nfreq, real_t **pt_freqs, bool *pt_io_period, bool *pt_readin_period,
//     MYINT *pt_nmode, MYINT **pt_modes, 
//     real_t ***pt_cdisp, MYINT ***pt_cdisp_iref, MYINT **pt_cdisp_num,
//     bool *pt_isRayl, GRT_MODEL1D **pt_mod1d, real_t depsrc, real_t deprcv)
// {  
//     int ncid;
//     int f_dimid, n_dimid;
//     int f_varid, n_varid;
//     int c_varid, ciref_varid, cnum_varid;

//     // 打开 NC 文件
//     NC_CHECK(nc_open(s_filepath, NC_NOWRITE, &ncid));

//     // 频散文件中记录的是周期还是频率
//     bool readin_period = false;
//     char *fname = strdup("freq");

//     // 读取一系列全局属性
//     // 读取模型路径
//     {
//         size_t m_len = 0;
//         NC_CHECK(nc_inq_attlen(ncid, NC_GLOBAL, "model", &m_len));
//         *pt_modelpath = (char *)calloc(m_len+1, sizeof(char));
//         NC_CHECK(nc_get_att_text(ncid, NC_GLOBAL, "model", *pt_modelpath));
//         (*pt_modelpath)[m_len] = '\0';

//         if((*pt_mod1d = grt_read_mod1d_from_file(command, *pt_modelpath, depsrc, deprcv, true)) ==NULL){
//             exit(EXIT_FAILURE);
//         }
//     }
//     // 读入的是频率还是周期
//     {
//         int readin_period_int = 0;
//         NC_CHECK(NC_FUNC_MYINT(nc_get_att) (ncid, NC_GLOBAL, "io_period", &readin_period_int));
//         readin_period = (bool)readin_period_int;
//         if(pt_readin_period!=NULL)  *pt_readin_period = readin_period;
//         if(readin_period){
//             GRT_SAFE_ASPRINTF(&fname, "period");
//         }
//     }
//     // Rayleigh or Love
//     {
//         int isRayl_int = 0;
//         NC_CHECK(NC_FUNC_MYINT(nc_get_att) (ncid, NC_GLOBAL, "isRayl", &isRayl_int));
//         *pt_isRayl = (bool)isRayl_int;
//     }
    
//     // 先读取 NC 文件的坐标变量
//     size_t in_nfreq;
//     size_t in_nmode;
//     NC_CHECK(nc_inq_dimid(ncid, fname, &f_dimid));
//     NC_CHECK(nc_inq_dimlen(ncid, f_dimid, &in_nfreq));
//     NC_CHECK(nc_inq_dimid(ncid, "mode", &n_dimid));
//     NC_CHECK(nc_inq_dimlen(ncid, n_dimid, &in_nmode));
//     real_t *in_freqs = (real_t *)calloc(in_nfreq, sizeof(real_t));
//     MYINT *in_modes = (MYINT *)calloc(in_nmode, sizeof(MYINT));
//     NC_CHECK(nc_inq_varid(ncid, fname, &f_varid));
//     NC_CHECK(NC_FUNC_REAL(nc_get_var) (ncid, f_varid, in_freqs));
//     // 如果为周期，则转为频率
//     if(readin_period){
//         for(size_t iw=0; iw < in_nfreq; ++iw){
//             in_freqs[iw] = 1.0 / in_freqs[iw];
//         }
//     }
//     NC_CHECK(nc_inq_varid(ncid, "mode", &n_varid));
//     NC_CHECK(NC_FUNC_MYINT(nc_get_var) (ncid, n_varid, in_modes));

//     // 每个频率下有多少阶
//     MYINT *in_cnum = (MYINT *)calloc(in_nfreq, sizeof(MYINT));
//     NC_CHECK(nc_inq_varid(ncid, "cnum", &cnum_varid));
//     NC_CHECK(NC_FUNC_MYINT(nc_get_var) (ncid, cnum_varid, in_cnum));

//     NC_CHECK(nc_inq_varid(ncid, "c", &c_varid));
//     NC_CHECK(nc_inq_varid(ncid, "ciref", &ciref_varid));


//     // 定义最终要处理的频点
//     // 若函数传入 NULL 频点数组，则使用 NC 文件的全部频点
//     if(*pt_freqs == NULL){
//         *pt_nfreq = in_nfreq;
//         *pt_freqs = (real_t *)calloc(in_nfreq, sizeof(real_t));
//         memcpy(*pt_freqs, in_freqs, sizeof(real_t)*in_nfreq);
//         if(pt_io_period!=NULL) *pt_io_period = readin_period;
//     }
//     // 函数传入了频点数组，但在前两个点定义了频率范围
//     else if(def_freq_period) {
//         real_t freqmin = (*pt_freqs)[0];
//         real_t freqmax = (*pt_freqs)[1];
//         if(pt_io_period!=NULL && *pt_io_period) {
//             freqmin = 1.0/(*pt_freqs)[1];
//             freqmax = 1.0/(*pt_freqs)[0];
//         }
//         GRT_SAFE_FREE_PTR(*pt_freqs);

//         // 动态统计数量和申请内存
//         *pt_freqs = NULL;
//         *pt_nfreq = 0;
//         real_t freq;
//         for(size_t iw=0; iw < in_nfreq; ++iw){
//             freq = in_freqs[iw];
//             if(freq < freqmin || freq > freqmax)  continue;

//             *pt_freqs = (real_t*)realloc(*pt_freqs, sizeof(real_t)*(*pt_nfreq+1));
//             (*pt_freqs)[*pt_nfreq] = freq;
//             (*pt_nfreq)++;
//         }
//     }

//     // 定义最终要处理的阶数
//     // 若函数传入 NULL 阶数数组，则使用 NC 文件的全部阶数
//     if(*pt_modes == NULL){
//         *pt_nmode = in_nmode;
//         *pt_modes = (MYINT *)calloc(in_nmode, sizeof(MYINT));
//         memcpy(*pt_modes, in_modes, sizeof(MYINT)*in_nmode);
//     }

//     // 根据频点数量申请内存
//     *pt_cdisp = (real_t **)calloc(*pt_nfreq, sizeof(real_t*));
//     *pt_cdisp_iref = (MYINT **)calloc(*pt_nfreq, sizeof(MYINT*));
//     *pt_cdisp_num = (MYINT *)calloc(*pt_nfreq, sizeof(MYINT));

//     // 循环每个频率、每个阶数，记录所需结果
//     // 逐个频率读取频散结果
//     for(size_t iw=0; iw < in_nfreq; ++iw){
//         real_t freq = in_freqs[iw];
//         MYINT idx_freq = grt_findClosest_real_t(*pt_freqs, *pt_nfreq, freq);
//         if(fabs((*pt_freqs)[idx_freq] - freq) > 1e-5)  continue;

//         size_t indexp[2];
//         indexp[0] = iw;

//         real_t root_c;
//         MYINT root_ciref;

//         MYINT imode = 0;
//         for(int im=0; im < in_cnum[iw]; ++im){
//             MYINT mode = in_modes[im];
//             MYINT idx_mode = grt_findElement_MYINT(*pt_modes, *pt_nmode, mode);
//             if(idx_mode < 0)  continue;

//             indexp[1] = im;
//             NC_CHECK(NC_FUNC_REAL(nc_get_var1) (ncid, c_varid, indexp, &root_c));
//             NC_CHECK(NC_FUNC_MYINT(nc_get_var1) (ncid, ciref_varid, indexp, &root_ciref));

//             // 为新阶扩容
//             (*pt_cdisp)[idx_freq] = (real_t *)realloc((*pt_cdisp)[idx_freq], sizeof(real_t)*(imode + 1));
//             (*pt_cdisp_iref)[idx_freq] = (MYINT *)realloc((*pt_cdisp_iref)[idx_freq], sizeof(MYINT)*(imode + 1));
//             (*pt_cdisp)[idx_freq][imode] = root_c;
//             (*pt_cdisp_iref)[idx_freq][imode] = root_ciref;

//             // printf("%f, %d, %d, %f\n", freq, mode, root_ciref, root_c);

//             imode++;
//         }

//         (*pt_cdisp_num)[idx_freq] = imode;
//     }

//     GRT_SAFE_FREE_PTR(in_freqs);
//     GRT_SAFE_FREE_PTR(in_modes);
//     GRT_SAFE_FREE_PTR(in_cnum);
//     GRT_SAFE_FREE_PTR(fname);
// }



// /** 根据 grt_eigenfn 的一些参数，计算本征函数、能量积分、群速度、相速度敏感核 */
// void grt_eigenfn_egy_gdisp_phaseK_util(
//     GRT_MODEL1D *mod1d, const MYINT nfreq, const real_t freqs[nfreq], const MYINT ncols, 
//     const real_t *const cdisp[nfreq], const MYINT *const cdisp_iref[nfreq], const MYINT cdisp_num[nfreq], 
//     const bool calc_egfn, const MYINT nz, const real_t zsamps[nz], const MYINT z_irefs[nz], MYCOMPLEX (**eigenfns)[nz][ncols], 
//     const bool calc_egyint, MYCOMPLEX (**EgyInt)[GRT_EGYINTS_MAX], 
//     const MYINT cpar_nz, const real_t *cpar_zs, const MYINT *cpar_z_irefs, MYCOMPLEX (**phase_K)[cpar_nz][GRT_SNSTVTY_MAX], 
//     real_t **gdisp, const bool isRayl)
// {
//     // 抹去mod1d中的Qinv，找本征值暂不考虑衰减
//     for(MYINT i=0; i<mod1d->n; ++i){
//         mod1d->Qa[i] = mod1d->Qb[i] = mod1d->Qainv[i] = mod1d->Qbinv[i] = 0.0;
//     }

//     // 最大零点数
//     MYINT max_nroots = 0;
//     for(MYINT iw=0; iw<nfreq; ++iw){
//         max_nroots = GRT_MAX(max_nroots, cdisp_num[iw]);
//     }
    
//     // 循环每个频率
//     #pragma omp parallel for schedule(guided) default(shared) collapse(2)
//     for(MYINT iw=0; iw<nfreq; ++iw){
//         // 未设置不同频率下的衰减系数，暂不支持
        
//         // 假设这些本征值都是从 iref 层的久期函数算出来的，
//         // 1. 先根据 iref（z_i+）的垂直波函数，计算出每个物理分界面上侧（z_j-）和下侧（z_j+）的垂直波函数
//         // 2. 循环每个深度采样点，使用上侧（z_{j+1}-）和下侧（z_j+）的垂直波函数推导采样点上的垂直波函数以及本征函数
//         for(MYINT ik=0; ik<max_nroots; ++ik){
//             if(ik >= cdisp_num[iw]) continue;

//             real_t omega = PI2*freqs[iw];

//             // 模型每个物理层介质下方 z_j+ 的垂直波函数
//             MYCOMPLEX (*mod_potRaylLove_Down)[ncols] = (MYCOMPLEX (*)[ncols])calloc(mod1d->n, sizeof(MYCOMPLEX)*ncols);
//             // 模型每个物理层介质上方 z_j- 的垂直波函数，申请 n+1 的内存，方便后续不必讨论“半空间中的上行波场”
//             MYCOMPLEX (*mod_potRaylLove_Up)[ncols] = (MYCOMPLEX (*)[ncols])calloc(mod1d->n+1, sizeof(MYCOMPLEX)*ncols);

//             real_t eigenK = omega/cdisp[iw][ik];
//             MYCOMPLEX potRaylLove[ncols];  memset(potRaylLove, 0, sizeof(MYCOMPLEX)*ncols);
//             MYCOMPLEX secRaylLove = 0.0;
//             MYINT stats = GRT_INVERSE_SUCCESS;
//             MYINT iref = cdisp_iref[iw][ik];
//             grt_secular_function_potential(mod1d, omega, eigenK, (iref==0)?0:1, iref, isRayl, &secRaylLove, potRaylLove, &stats);
 
//             // 计算出每个物理分界面上侧（z_{j+1}-）和 下侧（z_j+）的垂直波函数
//             grt_get_mod_potential_Up_Down(mod1d, omega, eigenK, ncols, mod_potRaylLove_Down, mod_potRaylLove_Up, iref, potRaylLove, isRayl, &stats);

//             // 计算不同深度的本征函数
//             if(calc_egfn){
//                 grt_get_eigenfn_depths(
//                     mod1d, omega, eigenK, ncols, 
//                     mod_potRaylLove_Down, mod_potRaylLove_Up, nz, zsamps, z_irefs, eigenfns[iw][ik], isRayl);
//             }

//             // 计算能量积分
//             if(calc_egyint){
//                 grt_energy_integrals(
//                     mod1d, omega, eigenK, ncols, 
//                     mod_potRaylLove_Down, mod_potRaylLove_Up, EgyInt[iw][ik], 
//                     cpar_nz, cpar_zs, cpar_z_irefs, phase_K[iw][ik], isRayl);
//                 // 计算群速度
//                 gdisp[iw][ik] = grt_get_group_velocity(omega, eigenK, EgyInt[iw][ik], isRayl);
//             }

//             GRT_SAFE_FREE_PTR(mod_potRaylLove_Down);
//             GRT_SAFE_FREE_PTR(mod_potRaylLove_Up);

//         }
        
//     }

    
// }


// /** 根据计算好的相速度、群速度和相速度敏感核，计算群速度敏感核 */
// void grt_group_sensitivity(
//     const MYINT nfreq, const real_t freqs[nfreq],
//     const real_t *const cdisp[nfreq], const MYINT cdisp_num[nfreq], 
//     const MYINT cpar_nz, MYCOMPLEX (**phase_K)[cpar_nz][GRT_SNSTVTY_MAX], real_t **gdisp, 
//     MYCOMPLEX (**group_K)[cpar_nz][GRT_SNSTVTY_MAX])
// {
//     // 根据群速度与相速度的关系： U = c^2 / ( c - w*dc/dw )
//     // 推导得到敏感核（其中 X 表示 alpha/beta/rho ）
//     //    dU/dX = U/c * (2 - U/c) * dc/dX + (U/c)^2 * w * d(dc/dX)/dw
//     // 后一项存在对相速度敏感核在频率上做差分，
//     // 因此输入的频率间隔大小会直接影响群速度敏感核精度

//     // 循环每个频率
//     for(MYINT iw=0; iw<nfreq; ++iw){
//         for(MYINT ik=0; ik<cdisp_num[iw]; ++ik){
//             real_t C = cdisp[iw][ik];
//             real_t U = gdisp[iw][ik];
//             real_t UoC = U/C;

//             // 循环每一层
//             for(MYINT iz=0; iz<cpar_nz; ++iz){
//                 // 层内的物理量
//                 for(MYINT ia=0; ia<GRT_SNSTVTY_MAX; ++ia){
//                     // 计算相速度敏感核在频率上的差分
//                     real_t ddCdwdX = 0.0;
//                     if(iw < nfreq-1){
//                         ddCdwdX = phase_K[iw+1][ik][iz][ia] - phase_K[iw][ik][iz][ia];
//                         ddCdwdX /= freqs[iw+1] - freqs[iw];
//                     }
                    
//                     // 合并公式，计算群速度敏感核
//                     group_K[iw][ik][iz][ia] = UoC * (2.0 - UoC) * phase_K[iw][ik][iz][ia] + GRT_SQUARE(UoC) * freqs[iw] * ddCdwdX;
//                 }
//             }
//         }
//     }
// }




// void grt_output_gdisp(
//     const char *filepath, 
//     const MYINT nfreq, const real_t *freqs, const bool io_period, const bool isRayl,
//     const MYINT *modes,
//     const real_t *const cdisp[nfreq], const real_t *const gdisp[nfreq], const MYINT cdisp_num[nfreq])
// {
//     int ncid, f_dimid, n_dimid;
//     const int ndims = 2;
//     int dimids[ndims];
//     int f_varid, n_varid;
//     int c_varid, g_varid, cnum_varid;

//     // 创建 NC 文件
//     NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

//     // 最大零点数
//     MYINT max_nroots = 0;
//     for(MYINT iw=0; iw<nfreq; ++iw){
//         max_nroots = GRT_MAX(max_nroots, cdisp_num[iw]);
//     }

//     // 定义维度
//     const char *fname = (io_period)? "period" : "freq";
//     NC_CHECK(nc_def_dim(ncid, fname, nfreq, &f_dimid));
//     NC_CHECK(nc_def_dim(ncid, "mode", max_nroots, &n_dimid));

//     // 定义变量维度数组
//     dimids[0] = f_dimid;
//     dimids[1] = n_dimid;

//     // 定义维度数组
//     NC_CHECK(nc_def_var(ncid, fname, NC_REAL, 1, &f_dimid, &f_varid));
//     NC_CHECK(nc_def_var(ncid, "mode", GRT_NC_MYINT, 1, &n_dimid, &n_varid));

//     // 定义频散数组
//     NC_CHECK(nc_def_var(ncid, "cnum", GRT_NC_MYINT, 1, &f_dimid, &cnum_varid));
//     NC_CHECK(nc_def_var(ncid, "c", NC_REAL, ndims, dimids, &c_varid));
//     NC_CHECK(nc_def_var(ncid, "g", NC_REAL, ndims, dimids, &g_varid));

//     // 填充值
//     const real_t fill_value_real_t = 0.0;
//     NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, c_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));
//     NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, g_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));

//     // 其它属性
//     // freqs or periods
//     MYINT io_period_int = (MYINT)io_period;
//     NC_CHECK(NC_FUNC_MYINT(nc_put_att) (ncid, NC_GLOBAL, "io_period", GRT_NC_MYINT, 1, &io_period_int));
//     // Rayleigh or Love
//     MYINT isRayl_int = (MYINT)isRayl;
//     NC_CHECK(NC_FUNC_MYINT(nc_put_att) (ncid, NC_GLOBAL, "isRayl", GRT_NC_MYINT, 1, &isRayl_int));
    
//     // 结束定义模式
//     NC_CHECK(nc_enddef(ncid));

//     // 写入数据
//     real_t *freqsOrperiods = (real_t *)calloc(nfreq, sizeof(real_t));
//     for(int iw=0; iw<nfreq; ++iw){
//         freqsOrperiods[iw] = (io_period)? 1.0/freqs[iw] : freqs[iw];
//     }

//     NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, f_varid, freqsOrperiods));
//     NC_CHECK(NC_FUNC_MYINT(nc_put_var) (ncid, n_varid, modes));
//     NC_CHECK(NC_FUNC_MYINT(nc_put_var) (ncid, cnum_varid, cdisp_num));

//     // 分批次写入数据
//     size_t startp[ndims];
//     size_t countp[ndims];
//     for(int iw=0; iw<nfreq; ++iw){
//         startp[0] = iw;
//         startp[1] = 0;
//         countp[0] = 1;
//         countp[1] = cdisp_num[iw];

//         NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, c_varid, startp, countp, cdisp[iw]));
//         NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, g_varid, startp, countp, gdisp[iw]));
//     }
    
//     // 关闭文件
//     NC_CHECK(nc_close(ncid));

//     GRT_SAFE_FREE_PTR(freqsOrperiods);

// }



// /* 输出本征函数结果 */
// void grt_output_eigenfns(
//     const char *filepath, 
//     const MYINT nfreq, const real_t *freqs, const bool io_period, const bool isRayl,
//     const MYINT *modes, const MYINT ncols, 
//     const real_t *const cdisp[nfreq], const MYINT cdisp_num[nfreq], 
//     const MYINT nz, const real_t zsamps[nz],
//     const MYCOMPLEX (*const *eigenfns)[nz][ncols])
// {
//     int ncid, f_dimid, n_dimid, z_dimid, e_dimid;
//     const int ndims = 4;
//     int dimids[ndims];
//     int f_varid, n_varid, z_varid;
//     int c_varid, cnum_varid;
//     int efn_varid;
//     const int nw = 2 * ncols; // 分实部和虚部

//     // 创建 NC 文件
//     NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

//     // 最大零点数
//     MYINT max_nroots = 0;
//     for(MYINT iw=0; iw<nfreq; ++iw){
//         max_nroots = GRT_MAX(max_nroots, cdisp_num[iw]);
//     }

//     // 定义维度
//     const char *fname = (io_period)? "period" : "freq";
//     NC_CHECK(nc_def_dim(ncid, fname, nfreq, &f_dimid));
//     NC_CHECK(nc_def_dim(ncid, "mode", max_nroots, &n_dimid));
//     NC_CHECK(nc_def_dim(ncid, "z", nz, &z_dimid));
//     NC_CHECK(nc_def_dim(ncid, "w", nw, &e_dimid));

//     // 定义变量维度数组
//     dimids[0] = f_dimid;
//     dimids[1] = n_dimid;
//     dimids[2] = z_dimid;
//     dimids[3] = e_dimid;

//     // 定义维度数组
//     NC_CHECK(nc_def_var(ncid, fname, NC_REAL, 1, &f_dimid, &f_varid));
//     NC_CHECK(nc_def_var(ncid, "mode", GRT_NC_MYINT, 1, &n_dimid, &n_varid));
//     NC_CHECK(nc_def_var(ncid, "z", NC_REAL, 1, &z_dimid, &z_varid));

//     // 定义频散数组
//     NC_CHECK(nc_def_var(ncid, "cnum", GRT_NC_MYINT, 1, &f_dimid, &cnum_varid));
//     NC_CHECK(nc_def_var(ncid, "c", NC_REAL, 2, dimids, &c_varid));   // 取前两维

//     // 填充值
//     const real_t fill_value_real_t = 0.0;
//     NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, c_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));

//     // 定义本征函数的实虚部
//     NC_CHECK(nc_def_var(ncid, "eigenfn", NC_REAL, ndims, dimids, &efn_varid));
//     NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, efn_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));
    
//     // 其它属性
//     // freqs or periods
//     MYINT io_period_int = (MYINT)io_period;
//     NC_CHECK(NC_FUNC_MYINT(nc_put_att) (ncid, NC_GLOBAL, "io_period", GRT_NC_MYINT, 1, &io_period_int));
//     // Rayleigh or Love
//     MYINT isRayl_int = (MYINT)isRayl;
//     NC_CHECK(NC_FUNC_MYINT(nc_put_att) (ncid, NC_GLOBAL, "isRayl", GRT_NC_MYINT, 1, &isRayl_int));
//     // 结束定义模式
//     NC_CHECK(nc_enddef(ncid));

//     // 写入数据
//     real_t *freqsOrperiods = (real_t *)calloc(nfreq, sizeof(real_t));
//     for(int iw=0; iw<nfreq; ++iw){
//         freqsOrperiods[iw] = (io_period)? 1.0/freqs[iw] : freqs[iw];
//     }

//     NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, f_varid, freqsOrperiods));
//     NC_CHECK(NC_FUNC_MYINT(nc_put_var) (ncid, n_varid, modes));
//     NC_CHECK(NC_FUNC_MYINT(nc_put_var) (ncid, cnum_varid, cdisp_num));
//     NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, z_varid, zsamps));

//     // 分批次写入数据
//     size_t startp[ndims];
//     size_t countp[ndims];

//     real_t (*realimag_part)[nw] = (real_t (*)[nw])calloc(nz, sizeof(real_t)*nw);

//     for(int iw=0; iw<nfreq; ++iw){
//         startp[0] = iw;
//         startp[1] = 0;
//         startp[2] = 0;
//         startp[3] = 0;
        
//         countp[0] = 1;
//         countp[1] = cdisp_num[iw];
//         countp[2] = nz;
//         countp[3] = nw;

//         NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, c_varid, startp, countp, cdisp[iw]));

//         for(int ik=0; ik < cdisp_num[iw]; ++ik){
//             memset(realimag_part, 0, sizeof(real_t)*nz*nw);

//             const MYCOMPLEX (* const pt_eigenfn)[ncols] = eigenfns[iw][ik];

//             // 将复数数据转为实部和虚部
//             for(int iz=0; iz < nz; ++iz){
//                 for(int j=0; j < ncols; ++j){
//                     realimag_part[iz][2*j]     = creal(pt_eigenfn[iz][j]);
//                     realimag_part[iz][2*j + 1] = cimag(pt_eigenfn[iz][j]);
//                 }
//             }

//             startp[1] = ik;
//             countp[1] = 1;
//             NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, efn_varid, startp, countp, realimag_part[0]));
//         }
//     }
    
//     // 关闭文件
//     NC_CHECK(nc_close(ncid));

//     GRT_SAFE_FREE_PTR(freqsOrperiods);
//     GRT_SAFE_FREE_PTR(realimag_part);

// }


// /* 输出能量积分结果 */
// void grt_output_energy_integrals(
//     const char *filepath, 
//     const MYINT nfreq, const real_t *freqs, const bool io_period, const bool isRayl,
//     const MYINT *modes, 
//     const real_t *const cdisp[nfreq], const MYINT cdisp_num[nfreq], 
//     const MYCOMPLEX (* const *EgyInt)[GRT_EGYINTS_MAX])
// {
//     int ncid, f_dimid, n_dimid, e_dimid;
//     const int ndims = 3;
//     int dimids[ndims];
//     int f_varid, n_varid, eint_varid;
//     int c_varid, cnum_varid;
//     int ne = 2 * (GRT_EGYINTS_MAX + 1);  // + 1 为了验证项

//     // 创建 NC 文件
//     NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

//     // 最大零点数
//     MYINT max_nroots = 0;
//     for(MYINT iw=0; iw<nfreq; ++iw){
//         max_nroots = GRT_MAX(max_nroots, cdisp_num[iw]);
//     }

//     // 定义维度
//     const char *fname = (io_period)? "period" : "freq";
//     NC_CHECK(nc_def_dim(ncid, fname, nfreq, &f_dimid));
//     NC_CHECK(nc_def_dim(ncid, "mode", max_nroots, &n_dimid));
//     NC_CHECK(nc_def_dim(ncid, "e", ne, &e_dimid));

//     // 定义变量维度数组
//     dimids[0] = f_dimid;
//     dimids[1] = n_dimid;
//     dimids[2] = e_dimid;

//     // 定义维度数组
//     NC_CHECK(nc_def_var(ncid, fname, NC_REAL, 1, &f_dimid, &f_varid));
//     NC_CHECK(nc_def_var(ncid, "mode", GRT_NC_MYINT, 1, &n_dimid, &n_varid));

//     // 定义频散数组
//     NC_CHECK(nc_def_var(ncid, "cnum", GRT_NC_MYINT, 1, &f_dimid, &cnum_varid));
//     NC_CHECK(nc_def_var(ncid, "c", NC_REAL, 2, dimids, &c_varid));

//     // 填充值
//     const real_t fill_value_real_t = 0.0;
//     NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, c_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));

//     // 定义能量积分
//     NC_CHECK(nc_def_var(ncid, "egyint", NC_REAL, ndims, dimids, &eint_varid));
//     NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, eint_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));

//     // 其它属性
//     // freqs or periods
//     MYINT io_period_int = (MYINT)io_period;
//     NC_CHECK(NC_FUNC_MYINT(nc_put_att) (ncid, NC_GLOBAL, "io_period", GRT_NC_MYINT, 1, &io_period_int));
//     // Rayleigh or Love
//     MYINT isRayl_int = (MYINT)isRayl;
//     NC_CHECK(NC_FUNC_MYINT(nc_put_att) (ncid, NC_GLOBAL, "isRayl", GRT_NC_MYINT, 1, &isRayl_int));
    
//     // 结束定义模式
//     NC_CHECK(nc_enddef(ncid));

//     // 写入数据
//     real_t *freqsOrperiods = (real_t *)calloc(nfreq, sizeof(real_t));
//     for(int iw=0; iw<nfreq; ++iw){
//         freqsOrperiods[iw] = (io_period)? 1.0/freqs[iw] : freqs[iw];
//     }

//     NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, f_varid, freqsOrperiods));
//     NC_CHECK(NC_FUNC_MYINT(nc_put_var) (ncid, n_varid, modes));
//     NC_CHECK(NC_FUNC_MYINT(nc_put_var) (ncid, cnum_varid, cdisp_num));

//     // 分批次写入数据
//     size_t startp[ndims];
//     size_t countp[ndims];

//     real_t (*realimag_part)[ne] = (real_t (*)[ne])calloc(max_nroots, sizeof(real_t)*ne);

//     for(int iw=0; iw<nfreq; ++iw){
//         startp[0] = iw;
//         startp[1] = 0;
//         startp[2] = 0;
        
//         countp[0] = 1;
//         countp[1] = cdisp_num[iw];
//         countp[2] = ne;

//         NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, c_varid, startp, countp, cdisp[iw]));

//         const MYCOMPLEX (* const pt_egyint)[GRT_EGYINTS_MAX] = EgyInt[iw];

//         memset(realimag_part, 0, sizeof(real_t)*max_nroots*ne);
//         for(int ik=0; ik < cdisp_num[iw]; ++ik){
//             for(int j=0; j < GRT_EGYINTS_MAX; ++j){
//                 realimag_part[ik][2*j] = creal(pt_egyint[ik][j]);
//                 realimag_part[ik][2*j + 1] = cimag(pt_egyint[ik][j]);
//             }

//             // 加上验证项
//             MYCOMPLEX res=-12345.0;
//             res = GRT_SQUARE(PI2*freqs[iw]) * pt_egyint[ik][0]
//                 - GRT_SQUARE(PI2*freqs[iw]/cdisp[iw][ik])*pt_egyint[ik][1]
//                 + (PI2*freqs[iw]/cdisp[iw][ik])*pt_egyint[ik][2]
//                 - pt_egyint[ik][3];

//             realimag_part[ik][ne-2] = creal(res);
//             realimag_part[ik][ne-1] = cimag(res);
//         }

//         NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, eint_varid, startp, countp, realimag_part[0]));
//     }

//     // 关闭文件
//     NC_CHECK(nc_close(ncid));

//     GRT_SAFE_FREE_PTR(freqsOrperiods);
//     GRT_SAFE_FREE_PTR(realimag_part);

// }


// /* 输出敏感核，相/群速度敏感核也应传入对应的相/群速度 */
// void grt_output_sensitivity(
//     const char *filepath, 
//     const MYINT nfreq, const real_t *freqs, const bool io_period, const bool isRayl,
//     const MYINT *modes,
//     const real_t *const cgdisp[nfreq], const MYINT cdisp_num[nfreq],
//     const char *UC, const MYINT cpar_nz, const real_t cpar_zs[cpar_nz], MYCOMPLEX (**sense_K)[cpar_nz][GRT_SNSTVTY_MAX])
// {
//     int ncid, f_dimid, n_dimid, z_dimid, k_dimid;
//     const int ndims = 4;
//     int dimids[ndims];
//     int f_varid, n_varid, z_varid, sens_varid;
//     int cg_varid, cnum_varid;
//     int nk = GRT_SNSTVTY_MAX;  // 敏感核仅取实部

//     // 创建 NC 文件
//     NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

//     // 最大零点数
//     MYINT max_nroots = 0;
//     for(MYINT iw=0; iw<nfreq; ++iw){
//         max_nroots = GRT_MAX(max_nroots, cdisp_num[iw]);
//     }

//     // 定义维度
//     const char *fname = (io_period)? "period" : "freq";
//     NC_CHECK(nc_def_dim(ncid, fname, nfreq, &f_dimid));
//     NC_CHECK(nc_def_dim(ncid, "mode", max_nroots, &n_dimid));
//     NC_CHECK(nc_def_dim(ncid, "z", cpar_nz, &z_dimid));
//     NC_CHECK(nc_def_dim(ncid, "k", nk, &k_dimid));

//     // 定义变量维度数组
//     dimids[0] = f_dimid;
//     dimids[1] = n_dimid;
//     dimids[2] = z_dimid;
//     dimids[3] = k_dimid;

//     // 定义维度数组
//     NC_CHECK(nc_def_var(ncid, fname, NC_REAL, 1, &f_dimid, &f_varid));
//     NC_CHECK(nc_def_var(ncid, "mode", GRT_NC_MYINT, 1, &n_dimid, &n_varid));
//     NC_CHECK(nc_def_var(ncid, "z", NC_REAL, 1, &z_dimid, &z_varid));

//     // 定义频散数组
//     NC_CHECK(nc_def_var(ncid, "cnum", GRT_NC_MYINT, 1, &f_dimid, &cnum_varid));
//     NC_CHECK(nc_def_var(ncid, UC, NC_REAL, 2, dimids, &cg_varid));   // 取前两维

//     // 填充值
//     const real_t fill_value_real_t = 0.0;
//     NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, cg_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));

//     // 敏感核（仅取实部）
//     char *sname = NULL;
//     GRT_SAFE_ASPRINTF(&sname, "%ssens", UC);
//     NC_CHECK(nc_def_var(ncid, sname, NC_REAL, ndims, dimids, &sens_varid));
//     GRT_SAFE_FREE_PTR(sname);
//     NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, sens_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));

//     // 其它属性
//     // freqs or periods
//     MYINT io_period_int = (MYINT)io_period;
//     NC_CHECK(NC_FUNC_MYINT(nc_put_att) (ncid, NC_GLOBAL, "io_period", GRT_NC_MYINT, 1, &io_period_int));
//     // Rayleigh or Love
//     MYINT isRayl_int = (MYINT)isRayl;
//     NC_CHECK(NC_FUNC_MYINT(nc_put_att) (ncid, NC_GLOBAL, "isRayl", GRT_NC_MYINT, 1, &isRayl_int));

//     // 结束定义模式
//     NC_CHECK(nc_enddef(ncid));

//     // 写入数据
//     real_t *freqsOrperiods = (real_t *)calloc(nfreq, sizeof(real_t));
//     for(int iw=0; iw<nfreq; ++iw){
//         freqsOrperiods[iw] = (io_period)? 1.0/freqs[iw] : freqs[iw];
//     }

//     NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, f_varid, freqsOrperiods));
//     NC_CHECK(NC_FUNC_MYINT(nc_put_var) (ncid, n_varid, modes));
//     NC_CHECK(NC_FUNC_MYINT(nc_put_var) (ncid, cnum_varid, cdisp_num));
//     NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, z_varid, cpar_zs));

//     // 分批次写入数据
//     size_t startp[ndims];
//     size_t countp[ndims];

//     real_t (*real_part)[nk] = (real_t (*)[nk])calloc(cpar_nz, sizeof(real_t)*nk);

//     for(int iw=0; iw<nfreq; ++iw){
//         startp[0] = iw;
//         startp[1] = 0;
//         startp[2] = 0;
//         startp[3] = 0;
        
//         countp[0] = 1;
//         countp[1] = cdisp_num[iw];
//         countp[2] = cpar_nz;
//         countp[3] = nk;

//         NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, cg_varid, startp, countp, cgdisp[iw]));

//         for(int ik=0; ik < cdisp_num[iw]; ++ik){
//             memset(real_part, 0, sizeof(real_t)*cpar_nz*nk);

//             const MYCOMPLEX (* const pt_sens)[nk] = sense_K[iw][ik];

//             // 将复数数据转为实部和虚部
//             for(int iz=0; iz < cpar_nz; ++iz){
//                 for(int j=0; j < nk; ++j){
//                     real_part[iz][j] = creal(pt_sens[iz][j]);
//                 }
//             }

//             startp[1] = ik;
//             countp[1] = 1;
//             NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, sens_varid, startp, countp, real_part[0]));
//         }
//     }

//     // 关闭文件
//     NC_CHECK(nc_close(ncid));

//     GRT_SAFE_FREE_PTR(freqsOrperiods);
//     GRT_SAFE_FREE_PTR(real_part);
// }
