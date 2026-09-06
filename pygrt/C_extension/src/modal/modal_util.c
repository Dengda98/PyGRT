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

#include "grt/modal/modal_util.h"
#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/common/mynetcdf.h"
#include "grt/common/checkerror.h"




/** 统计展平后的 freqmode 点数 */
static size_t count_nfreqmode(const EIGENV *eigv, size_t nf)
{
    size_t nfm = 0;
    for(size_t iw = 0; iw < nf; ++iw){
        nfm += eigv[iw].n;
    }
    return nfm;
}

/** 取 mode 维长度：至少为已有 nmode，且不小于各频率实际根数 */
static size_t count_nmode(const EIGENV *eigv, size_t nf, size_t nmode)
{
    for(size_t iw = 0; iw < nf; ++iw){
        nmode = GRT_MAX(nmode, eigv[iw].n);
    }
    return nmode;
}

/** 写入 isRayl 全局属性 */
static void nc_put_isRayl(int ncid, DISPER_TYPE wtype)
{
    int isRayl_int = (int)(wtype == GRT_DISPERSION_RAYL);
    NC_CHECK(NC_FUNC_INT(nc_put_att)(ncid, NC_GLOBAL, "isRayl", NC_INT, 1, &isRayl_int));
}

/** 定义 freq / freqmode / mode 维及 freq、cnum、mode 变量 */
static void nc_def_freq_freqmode(
    int ncid, size_t nf, size_t nfm, size_t nmode,
    int *f_dimid, int *fm_dimid,
    int *f_varid, int *cnum_varid, int *mode_varid)
{
    int n_dimid;
    NC_CHECK(nc_def_dim(ncid, "freq", nf, f_dimid));
    NC_CHECK(nc_def_dim(ncid, "freqmode", nfm, fm_dimid));
    NC_CHECK(nc_def_dim(ncid, "mode", nmode, &n_dimid));
    NC_CHECK(nc_def_var(ncid, "freq", NC_REAL, 1, f_dimid, f_varid));
    NC_CHECK(nc_def_var(ncid, "cnum", NC_INT, 1, f_dimid, cnum_varid));
    NC_CHECK(nc_def_var(ncid, "mode", NC_INT, 1, &n_dimid, mode_varid));
}

/**
 * 写入 freq、cnum、mode(mode) 以及展平后的 c / u / ciref
 * u_varid / ciref_varid 为负表示不写入对应变量
 */
static void nc_put_freq_freqmode(
    int ncid, size_t nf, size_t nfm, size_t nmode,
    const real_t *freqs, const EIGENV *eigv, const size_t *modes,
    int f_varid, int cnum_varid, int mode_varid,
    int c_varid, int u_varid, int ciref_varid)
{
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, f_varid, freqs));

    int *cnum = (int *)calloc(nf, sizeof(int));
    int *mode_vals = (int *)calloc(nmode, sizeof(int));
    real_t *c_flat = (real_t *)calloc(nfm, sizeof(real_t));
    real_t *u_flat = (u_varid >= 0) ? (real_t *)calloc(nfm, sizeof(real_t)) : NULL;
    int *ciref_flat = (ciref_varid >= 0) ? (int *)calloc(nfm, sizeof(int)) : NULL;

    for(size_t i = 0; i < nmode; ++i){
        mode_vals[i] = (modes != NULL) ? (int)modes[i] : (int)i;
    }

    size_t k = 0;
    for(size_t iw = 0; iw < nf; ++iw){
        cnum[iw] = (int)eigv[iw].n;
        for(size_t ic = 0; ic < eigv[iw].n; ++ic){
            c_flat[k] = eigv[iw].c_roots[ic];
            if(u_flat != NULL)     u_flat[k] = eigv[iw].u_roots[ic];
            if(ciref_flat != NULL) ciref_flat[k] = (int)eigv[iw].c_roots_iref[ic];
            k++;
        }
    }

    NC_CHECK(NC_FUNC_INT(nc_put_var)(ncid, cnum_varid, cnum));
    NC_CHECK(NC_FUNC_INT(nc_put_var)(ncid, mode_varid, mode_vals));
    if(c_varid >= 0){
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, c_varid, c_flat));
    }
    if(u_varid >= 0){
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, u_varid, u_flat));
    }
    if(ciref_varid >= 0){
        NC_CHECK(NC_FUNC_INT(nc_put_var)(ncid, ciref_varid, ciref_flat));
    }

    GRT_SAFE_FREE_PTR(cnum);
    GRT_SAFE_FREE_PTR(mode_vals);
    GRT_SAFE_FREE_PTR(c_flat);
    GRT_SAFE_FREE_PTR(u_flat);
    GRT_SAFE_FREE_PTR(ciref_flat);
}

/** 打包模型前四列，供写入 NC */
static real_t (*pack_modarr_nc(size_t nlayer, const real_t (*modarr)[GRT_MODARR_NCOL]))[GRT_MODAL_MODARR_NCOL]
{
    real_t (*out)[GRT_MODAL_MODARR_NCOL] = (real_t (*)[GRT_MODAL_MODARR_NCOL])malloc(
        sizeof(real_t) * GRT_MODAL_MODARR_NCOL * nlayer);
    for(size_t i = 0; i < nlayer; ++i){
        memcpy(out[i], modarr[i], sizeof(real_t) * GRT_MODAL_MODARR_NCOL);
    }
    return out;
}

/** 将 NC 中的四列模型展开为完整矩阵（Qa/Qb 置 0） */
static real_t (*unpack_modarr_nc(size_t nlayer, const real_t (*in)[GRT_MODAL_MODARR_NCOL]))[GRT_MODARR_NCOL]
{
    real_t (*out)[GRT_MODARR_NCOL] = (real_t (*)[GRT_MODARR_NCOL])calloc(
        nlayer, sizeof(real_t) * GRT_MODARR_NCOL);
    for(size_t i = 0; i < nlayer; ++i){
        memcpy(out[i], in[i], sizeof(real_t) * GRT_MODAL_MODARR_NCOL);
    }
    return out;
}


void grt_output_cdisp(
    const char *filepath, const char *full_command, const char *modelname,
    const MODEL1D *mod1d, EIGENV_INFO *eigmet)
{
    if(modelname == NULL || mod1d == NULL || mod1d->nmodarr == 0 || mod1d->modarr == NULL){
        GRTRaiseError("modelname/modarr must be non-empty when writing dispersion.");
    }

    const size_t nlayer = mod1d->nmodarr;
    const real_t (*modarr)[GRT_MODARR_NCOL] = (const real_t (*)[GRT_MODARR_NCOL])mod1d->modarr;

    int ncid, f_dimid, fm_dimid;
    int f_varid, mode_varid;
    int c_varid, ciref_varid, cnum_varid;
    int layer_dimid, param_dimid, model_varid;

    EIGENV *eigv = eigmet->eigv;
    size_t nfm = count_nfreqmode(eigv, eigmet->nf);
    size_t nmode = count_nmode(eigv, eigmet->nf, 0);
    if(nfm == 0){
        GRTRaiseError("No eigenvalues were found, please check.");
    }

    // 创建 NC 文件
    NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

    // 写入生成频散的命令和模型名
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "command", strlen(full_command), full_command));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "modelname", strlen(modelname), modelname));

    // 定义维度和变量
    nc_def_freq_freqmode(ncid, eigmet->nf, nfm, nmode, &f_dimid, &fm_dimid, &f_varid, &cnum_varid, &mode_varid);
    NC_CHECK(nc_def_dim(ncid, "layer", nlayer, &layer_dimid));
    NC_CHECK(nc_def_dim(ncid, "model_param", GRT_MODAL_MODARR_NCOL, &param_dimid));
    {
        int model_dimids[2] = {layer_dimid, param_dimid};
        NC_CHECK(nc_def_var(ncid, "model", NC_REAL, 2, model_dimids, &model_varid));
    }
    NC_CHECK(nc_def_var(ncid, "c", NC_REAL, 1, &fm_dimid, &c_varid));
    NC_CHECK(nc_def_var(ncid, "ciref", NC_INT, 1, &fm_dimid, &ciref_varid));

    nc_put_isRayl(ncid, eigmet->wtype);

    // 结束定义模式
    NC_CHECK(nc_enddef(ncid));

    // 写入数据
    nc_put_freq_freqmode(
        ncid, eigmet->nf, nfm, nmode, eigmet->freqs, eigv, NULL,
        f_varid, cnum_varid, mode_varid, c_varid, -1, ciref_varid);
    {
        real_t (*model_nc)[GRT_MODAL_MODARR_NCOL] = pack_modarr_nc(nlayer, modarr);
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, model_varid, (const real_t *)model_nc));
        GRT_SAFE_FREE_PTR(model_nc);
    }

    // 关闭文件
    NC_CHECK(nc_close(ncid));
}


void grt_group_sensitivity(EIGENFN_INFO *eigfnmet)
{
    // 循环每个频率
    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        for(size_t ic = 0; ic < eigfnmet->eigv[iw].n; ++ic){
            real_t C = eigfnmet->eigv[iw].c_roots[ic];
            real_t U = eigfnmet->eigv[iw].u_roots[ic];
            real_t UoC = U/C;

            EIGENFN *eigfn = &eigfnmet->eigfn[iw][ic];

            // 循环每一层
            for(size_t iz = 0; iz < eigfnmet->cpar_nz; ++iz){
                // 层内的物理量
                for(int ia = 0; ia < GRT_SNSTVTY_MAX; ++ia){
                    // 计算相速度敏感核在频率上的差分
                    real_t ddCdwdX = 0.0;
                    if(iw < eigfnmet->nf-1){
                        EIGENFN *eigfn_w1 = &eigfnmet->eigfn[iw+1][ic];

                        ddCdwdX = eigfn_w1->csens[iz][ia] - eigfn->csens[iz][ia];
                        ddCdwdX /= eigfnmet->freqs[iw+1] - eigfnmet->freqs[iw];
                    }

                    // 合并公式，计算群速度敏感核
                    eigfn->usens[iz][ia] = UoC * (2.0 - UoC) * eigfn->csens[iz][ia] + GRT_SQUARE(UoC) * eigfnmet->freqs[iw] * ddCdwdX;
                }
            }
        }
    }
}


void grt_output_udisp(const char *filepath, EIGENFN_INFO *eigfnmet)
{
    int ncid, f_dimid, fm_dimid;
    int f_varid, mode_varid;
    int c_varid, u_varid, cnum_varid;

    size_t nfm = count_nfreqmode(eigfnmet->eigv, eigfnmet->nf);
    size_t nmode = count_nmode(eigfnmet->eigv, eigfnmet->nf, eigfnmet->nmode);
    if(nfm == 0){
        GRTRaiseError("No eigenvalues were found, please check.");
    }

    // 创建 NC 文件
    NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

    nc_def_freq_freqmode(ncid, eigfnmet->nf, nfm, nmode, &f_dimid, &fm_dimid, &f_varid, &cnum_varid, &mode_varid);
    NC_CHECK(nc_def_var(ncid, "c", NC_REAL, 1, &fm_dimid, &c_varid));
    NC_CHECK(nc_def_var(ncid, "u", NC_REAL, 1, &fm_dimid, &u_varid));

    nc_put_isRayl(ncid, eigfnmet->wtype);

    // 结束定义模式
    NC_CHECK(nc_enddef(ncid));

    nc_put_freq_freqmode(
        ncid, eigfnmet->nf, nfm, nmode, eigfnmet->freqs, eigfnmet->eigv, eigfnmet->modes,
        f_varid, cnum_varid, mode_varid, c_varid, u_varid, -1);

    // 关闭文件
    NC_CHECK(nc_close(ncid));
}


/* 输出本征函数结果 */
void grt_output_eigenfns(const char *filepath, const int ncols, EIGENFN_INFO *eigfnmet)
{
    int ncid, f_dimid, fm_dimid, z_dimid, e_dimid;
    const int ndims = 3;
    int dimids[ndims];
    int f_varid, mode_varid, z_varid;
    int c_varid, cnum_varid;
    int efn_varid;
    const int nw = 2 * ncols; // 分实部和虚部

    size_t nfm = count_nfreqmode(eigfnmet->eigv, eigfnmet->nf);
    size_t nmode = count_nmode(eigfnmet->eigv, eigfnmet->nf, eigfnmet->nmode);
    if(nfm == 0){
        GRTRaiseError("No eigenvalues were found, please check.");
    }

    // 创建 NC 文件
    NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

    nc_def_freq_freqmode(ncid, eigfnmet->nf, nfm, nmode, &f_dimid, &fm_dimid, &f_varid, &cnum_varid, &mode_varid);
    NC_CHECK(nc_def_dim(ncid, "z", eigfnmet->nz, &z_dimid));
    NC_CHECK(nc_def_dim(ncid, "w", nw, &e_dimid));

    dimids[0] = fm_dimid;
    dimids[1] = z_dimid;
    dimids[2] = e_dimid;

    NC_CHECK(nc_def_var(ncid, "z", NC_REAL, 1, &z_dimid, &z_varid));
    NC_CHECK(nc_def_var(ncid, "c", NC_REAL, 1, &fm_dimid, &c_varid));
    NC_CHECK(nc_def_var(ncid, "eigfn", NC_REAL, ndims, dimids, &efn_varid));

    nc_put_isRayl(ncid, eigfnmet->wtype);

    // 结束定义模式
    NC_CHECK(nc_enddef(ncid));

    nc_put_freq_freqmode(
        ncid, eigfnmet->nf, nfm, nmode, eigfnmet->freqs, eigfnmet->eigv, eigfnmet->modes,
        f_varid, cnum_varid, mode_varid, c_varid, -1, -1);
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, z_varid, eigfnmet->zs));

    size_t startp[ndims];
    size_t countp[ndims];
    startp[1] = 0;
    startp[2] = 0;
    countp[0] = 1;
    countp[1] = eigfnmet->nz;
    countp[2] = nw;

    real_t (*realimag_part)[nw] = (real_t (*)[nw])calloc(eigfnmet->nz, sizeof(real_t)*nw);

    size_t k = 0;
    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        for(size_t ic = 0; ic < eigfnmet->eigv[iw].n; ++ic){
            memset(realimag_part, 0, sizeof(real_t)*eigfnmet->nz*nw);

            cplx_t (*fn)[4] = eigfnmet->eigfn[iw][ic].fn;

            // 将复数数据转为实部和虚部
            for(size_t iz = 0; iz < eigfnmet->nz; ++iz){
                for(int j = 0; j < ncols; ++j){
                    realimag_part[iz][2*j]     = creal(fn[iz][j]);
                    realimag_part[iz][2*j + 1] = cimag(fn[iz][j]);
                }
            }

            startp[0] = k;
            NC_CHECK(NC_FUNC_REAL(nc_put_vara)(ncid, efn_varid, startp, countp, realimag_part[0]));
            k++;
        }
    }

    // 关闭文件
    NC_CHECK(nc_close(ncid));

    GRT_SAFE_FREE_PTR(realimag_part);
}


/* 输出能量积分结果 */
void grt_output_energy_integrals(const char *filepath, EIGENFN_INFO *eigfnmet)
{
    int ncid, f_dimid, fm_dimid, e_dimid;
    const int ndims = 2;
    int dimids[ndims];
    int f_varid, mode_varid, eint_varid;
    int c_varid, cnum_varid;
    int ne = 2 * (GRT_EGYINTS_MAX + 1);  // + 1 为了验证项

    size_t nfm = count_nfreqmode(eigfnmet->eigv, eigfnmet->nf);
    size_t nmode = count_nmode(eigfnmet->eigv, eigfnmet->nf, eigfnmet->nmode);
    if(nfm == 0){
        GRTRaiseError("No eigenvalues were found, please check.");
    }

    // 创建 NC 文件
    NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

    nc_def_freq_freqmode(ncid, eigfnmet->nf, nfm, nmode, &f_dimid, &fm_dimid, &f_varid, &cnum_varid, &mode_varid);
    NC_CHECK(nc_def_dim(ncid, "e", ne, &e_dimid));

    dimids[0] = fm_dimid;
    dimids[1] = e_dimid;

    NC_CHECK(nc_def_var(ncid, "c", NC_REAL, 1, &fm_dimid, &c_varid));
    NC_CHECK(nc_def_var(ncid, "egyint", NC_REAL, ndims, dimids, &eint_varid));

    nc_put_isRayl(ncid, eigfnmet->wtype);

    // 结束定义模式
    NC_CHECK(nc_enddef(ncid));

    nc_put_freq_freqmode(
        ncid, eigfnmet->nf, nfm, nmode, eigfnmet->freqs, eigfnmet->eigv, eigfnmet->modes,
        f_varid, cnum_varid, mode_varid, c_varid, -1, -1);

    size_t startp[ndims];
    size_t countp[ndims];
    startp[1] = 0;
    countp[1] = ne;

    real_t (*realimag_part)[ne] = (real_t (*)[ne])calloc(nmode, sizeof(real_t)*ne);

    size_t k = 0;
    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        size_t n = eigfnmet->eigv[iw].n;
        real_t omega = PI2*eigfnmet->freqs[iw];
        real_t *c_roots = eigfnmet->eigv[iw].c_roots;

        memset(realimag_part, 0, sizeof(real_t)*nmode*ne);
        for(size_t ic = 0; ic < n; ++ic){
            cplx_t *egyint = eigfnmet->eigfn[iw][ic].egyint;

            for(int j=0; j < GRT_EGYINTS_MAX; ++j){
                realimag_part[ic][2*j]     = creal(egyint[j]);
                realimag_part[ic][2*j + 1] = cimag(egyint[j]);
            }

            // 加上验证项
            cplx_t res=-12345.0;
            real_t kk = omega/c_roots[ic];
            res = GRT_SQUARE(omega) * egyint[0] - kk*kk*egyint[1] + kk*egyint[2] - egyint[3];
            realimag_part[ic][ne-2] = creal(res);
            realimag_part[ic][ne-1] = cimag(res);
        }

        startp[0] = k;
        countp[0] = n;
        NC_CHECK(NC_FUNC_REAL(nc_put_vara)(ncid, eint_varid, startp, countp, realimag_part[0]));
        k += n;
    }

    // 关闭文件
    NC_CHECK(nc_close(ncid));

    GRT_SAFE_FREE_PTR(realimag_part);
}


void grt_output_sensitivity(const char *filepath, const char *char_uc, EIGENFN_INFO *eigfnmet)
{
    int ncid, f_dimid, fm_dimid, z_dimid, k_dimid;
    const int ndims = 3;
    int dimids[ndims];
    int f_varid, mode_varid, z_varid, sens_varid;
    int cu_varid, cnum_varid;
    int nk = GRT_SNSTVTY_MAX;  // 敏感核仅取实部

    size_t nfm = count_nfreqmode(eigfnmet->eigv, eigfnmet->nf);
    size_t nmode = count_nmode(eigfnmet->eigv, eigfnmet->nf, eigfnmet->nmode);
    if(nfm == 0){
        GRTRaiseError("No eigenvalues were found, please check.");
    }

    // 创建 NC 文件
    NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

    nc_def_freq_freqmode(ncid, eigfnmet->nf, nfm, nmode, &f_dimid, &fm_dimid, &f_varid, &cnum_varid, &mode_varid);
    NC_CHECK(nc_def_dim(ncid, "z", eigfnmet->cpar_nz, &z_dimid));
    NC_CHECK(nc_def_dim(ncid, "k", nk, &k_dimid));

    dimids[0] = fm_dimid;
    dimids[1] = z_dimid;
    dimids[2] = k_dimid;

    NC_CHECK(nc_def_var(ncid, "z", NC_REAL, 1, &z_dimid, &z_varid));
    NC_CHECK(nc_def_var(ncid, char_uc, NC_REAL, 1, &fm_dimid, &cu_varid));

    char *sname = NULL;
    GRT_SAFE_ASPRINTF(&sname, "%ssens", char_uc);
    NC_CHECK(nc_def_var(ncid, sname, NC_REAL, ndims, dimids, &sens_varid));
    GRT_SAFE_FREE_PTR(sname);

    nc_put_isRayl(ncid, eigfnmet->wtype);

    // 结束定义模式
    NC_CHECK(nc_enddef(ncid));

    int u_varid = -1, c_varid = -1;
    if(strcmp(char_uc, "c") == 0){
        c_varid = cu_varid;
    } else if(strcmp(char_uc, "u") == 0){
        u_varid = cu_varid;
    } else {
        GRTRaiseError("Wrong execution.");
    }
    nc_put_freq_freqmode(
        ncid, eigfnmet->nf, nfm, nmode, eigfnmet->freqs, eigfnmet->eigv, eigfnmet->modes,
        f_varid, cnum_varid, mode_varid, c_varid, u_varid, -1);
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, z_varid, eigfnmet->cpar_zs));

    size_t startp[ndims];
    size_t countp[ndims];
    startp[1] = 0;
    startp[2] = 0;
    countp[0] = 1;
    countp[1] = eigfnmet->cpar_nz;
    countp[2] = nk;

    real_t (*real_part)[nk] = (real_t (*)[nk])calloc(eigfnmet->cpar_nz, sizeof(real_t)*nk);

    size_t ifm = 0;
    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        for(size_t ic = 0; ic < eigfnmet->eigv[iw].n; ++ic){
            memset(real_part, 0, sizeof(real_t)*eigfnmet->cpar_nz*nk);

            cplx_t (*cusens)[GRT_SNSTVTY_MAX] = NULL;
            if(strcmp(char_uc, "c") == 0){
                cusens = eigfnmet->eigfn[iw][ic].csens;
            } else {
                cusens = eigfnmet->eigfn[iw][ic].usens;
            }

            // 将复数数据转为实部
            for(size_t iz = 0; iz < eigfnmet->cpar_nz; ++iz){
                for(int j = 0; j < nk; ++j){
                    real_part[iz][j] = creal(cusens[iz][j]);
                }
            }

            startp[0] = ifm;
            NC_CHECK(NC_FUNC_REAL(nc_put_vara)(ncid, sens_varid, startp, countp, real_part[0]));
            ifm++;
        }
    }

    // 关闭文件
    NC_CHECK(nc_close(ncid));

    GRT_SAFE_FREE_PTR(real_part);
}


/** 读取相/群速度频散结果 */
void grt_read_dispersion(
    const char *filepath, EIGENV_INFO *eigmet, char **pt_modelname, MODEL1D **pt_mod1d)
{
    int ncid;
    int f_dimid, fm_dimid;
    int f_varid, mode_varid;
    int c_varid, u_varid, ciref_varid, cnum_varid;

    // 打开 NC 文件
    GRTCheckFileExist(filepath);
    NC_CHECK(nc_open(filepath, NC_NOWRITE, &ncid));

    // 根据是否有变量 u 来判断这个文件记录的是群速度还是相速度
    bool isGroup = false;
    {
        int ret = nc_inq_varid(ncid, "u", &u_varid);
        if(ret == NC_NOERR){
            isGroup = true;
        } else if(ret == NC_ENOTVAR) {
            isGroup = false;
        }
        else {
            NC_CHECK(ret);
        }
    }

    // 在相速度文件中读取模型名和模型数组（仅四列：Thk/Va/Vb/Rho）
    if(! isGroup){
        int layer_dimid, param_dimid, model_varid;
        size_t nlayer = 0, nparam = 0;
        if(nc_inq_dimid(ncid, "layer", &layer_dimid) != NC_NOERR ||
           nc_inq_dimid(ncid, "model_param", &param_dimid) != NC_NOERR){
            NC_CHECK(nc_close(ncid));
            GRTRaiseError(
                "Invalid dispersion nc \"%s\": missing model dimensions "
                "(layer, model_param).", filepath);
        }
        NC_CHECK(nc_inq_dimlen(ncid, layer_dimid, &nlayer));
        NC_CHECK(nc_inq_dimlen(ncid, param_dimid, &nparam));
        if(nlayer == 0 || nparam != GRT_MODAL_MODARR_NCOL){
            NC_CHECK(nc_close(ncid));
            GRTRaiseError(
                "Invalid dispersion nc \"%s\": layer=%zu, model_param=%zu "
                "(expect layer>0, model_param=%d).",
                filepath, nlayer, nparam, GRT_MODAL_MODARR_NCOL);
        }
        NC_CHECK(nc_inq_varid(ncid, "model", &model_varid));
        real_t (*model_nc)[GRT_MODAL_MODARR_NCOL] = (real_t (*)[GRT_MODAL_MODARR_NCOL])malloc(
            sizeof(real_t) * GRT_MODAL_MODARR_NCOL * nlayer);
        NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, model_varid, (real_t *)model_nc));
        real_t (*modarr)[GRT_MODARR_NCOL] = unpack_modarr_nc(nlayer, (const real_t (*)[GRT_MODAL_MODARR_NCOL])model_nc);
        GRT_SAFE_FREE_PTR(model_nc);

        size_t m_len = 0;
        NC_CHECK(nc_inq_attlen(ncid, NC_GLOBAL, "modelname", &m_len));
        char *modelname = (char *)calloc(m_len+1, sizeof(char));
        NC_CHECK(nc_get_att_text(ncid, NC_GLOBAL, "modelname", modelname));
        modelname[m_len] = '\0';

        if(pt_modelname != NULL)  *pt_modelname = modelname;
        else GRT_SAFE_FREE_PTR(modelname);
        if(pt_mod1d != NULL){
            *pt_mod1d = grt_read_mod1d_from_modarr(nlayer, modarr, -1.0, -1.0, true);
        }
        GRT_SAFE_FREE_PTR(modarr);
    } else {
        // 群速度文件不含内嵌模型
        if(pt_modelname != NULL) *pt_modelname = NULL;
        if(pt_mod1d != NULL){
            GRTRaiseError(
                "Group-velocity nc \"%s\" has no embedded model; use a phase-velocity file.",
                filepath);
        }
    }

    // Rayleigh or Love
    {
        int isRayl_int = 0;
        NC_CHECK(NC_FUNC_INT(nc_get_att)(ncid, NC_GLOBAL, "isRayl", &isRayl_int));
        eigmet->wtype = (isRayl_int)? GRT_DISPERSION_RAYL : GRT_DISPERSION_LOVE;
    }

    // 读取 freq / freqmode / mode 维
    int n_dimid;
    size_t nfm = 0, nmode = 0;
    NC_CHECK(nc_inq_dimid(ncid, "freq", &f_dimid));
    NC_CHECK(nc_inq_dimlen(ncid, f_dimid, &eigmet->nf));
    NC_CHECK(nc_inq_dimid(ncid, "freqmode", &fm_dimid));
    NC_CHECK(nc_inq_dimlen(ncid, fm_dimid, &nfm));
    NC_CHECK(nc_inq_dimid(ncid, "mode", &n_dimid));
    NC_CHECK(nc_inq_dimlen(ncid, n_dimid, &nmode));
    if(nmode == 0){
        GRTRaiseError("Invalid dispersion nc \"%s\": empty mode dimension.", filepath);
    }

    eigmet->freqs = (real_t *)calloc(eigmet->nf, sizeof(real_t));
    NC_CHECK(nc_inq_varid(ncid, "freq", &f_varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, f_varid, eigmet->freqs));

    int *cnum = (int *)calloc(eigmet->nf, sizeof(int));
    NC_CHECK(nc_inq_varid(ncid, "cnum", &cnum_varid));
    NC_CHECK(NC_FUNC_INT(nc_get_var)(ncid, cnum_varid, cnum));

    size_t nfm_expect = 0;
    int cnum_max = 0;
    for(size_t iw = 0; iw < eigmet->nf; ++iw){
        if(cnum[iw] < 0){
            GRTRaiseError("Invalid dispersion nc \"%s\": cnum[%zu]=%d.", filepath, iw, cnum[iw]);
        }
        nfm_expect += (size_t)cnum[iw];
        if(cnum[iw] > cnum_max) cnum_max = cnum[iw];
    }
    if(nfm_expect != nfm){
        GRTRaiseError(
            "Invalid dispersion nc \"%s\": sum(cnum)=%zu but freqmode=%zu.",
            filepath, nfm_expect, nfm);
    }
    if((size_t)cnum_max > nmode){
        GRTRaiseError(
            "Invalid dispersion nc \"%s\": max(cnum)=%d exceeds mode dimension %zu.",
            filepath, cnum_max, nmode);
    }

    int *mode_vals = (int *)calloc(nmode, sizeof(int));
    real_t *c_flat = (real_t *)calloc(nfm, sizeof(real_t));
    real_t *u_flat = isGroup ? (real_t *)calloc(nfm, sizeof(real_t)) : NULL;
    int *ciref_flat = (!isGroup) ? (int *)calloc(nfm, sizeof(int)) : NULL;

    NC_CHECK(nc_inq_varid(ncid, "mode", &mode_varid));
    NC_CHECK(NC_FUNC_INT(nc_get_var)(ncid, mode_varid, mode_vals));
    NC_CHECK(nc_inq_varid(ncid, "c", &c_varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, c_varid, c_flat));
    if(isGroup){
        NC_CHECK(nc_inq_varid(ncid, "u", &u_varid));
        NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, u_varid, u_flat));
    } else {
        NC_CHECK(nc_inq_varid(ncid, "ciref", &ciref_varid));
        NC_CHECK(NC_FUNC_INT(nc_get_var)(ncid, ciref_varid, ciref_flat));
    }

    eigmet->eigv = (EIGENV *)calloc(eigmet->nf, sizeof(EIGENV));

    // 按频率还原锯齿状频散；第 iw 频的阶为 mode[0 .. cnum[iw]-1]
    size_t k = 0;
    for(size_t iw = 0; iw < eigmet->nf; ++iw){
        eigmet->eigv[iw].n = (size_t)cnum[iw];
        eigmet->eigv[iw].c_roots = (real_t *)calloc(cnum[iw], sizeof(real_t));
        eigmet->eigv[iw].u_roots = (real_t *)calloc(cnum[iw], sizeof(real_t));
        eigmet->eigv[iw].c_roots_iref = (size_t *)calloc(cnum[iw], sizeof(size_t));

        for(int ic = 0; ic < cnum[iw]; ++ic){
            eigmet->eigv[iw].c_roots[ic] = c_flat[k];
            if(isGroup){
                eigmet->eigv[iw].u_roots[ic] = u_flat[k];
            } else {
                eigmet->eigv[iw].c_roots_iref[ic] = (size_t)ciref_flat[k];
            }
            k++;
        }
    }

    eigmet->nmode = nmode;
    eigmet->modes = (size_t *)calloc(nmode, sizeof(size_t));
    for(size_t i = 0; i < nmode; ++i){
        eigmet->modes[i] = (size_t)mode_vals[i];
    }

    // 关闭文件
    NC_CHECK(nc_close(ncid));

    GRT_SAFE_FREE_PTR(cnum);
    GRT_SAFE_FREE_PTR(mode_vals);
    GRT_SAFE_FREE_PTR(c_flat);
    GRT_SAFE_FREE_PTR(u_flat);
    GRT_SAFE_FREE_PTR(ciref_flat);
}
