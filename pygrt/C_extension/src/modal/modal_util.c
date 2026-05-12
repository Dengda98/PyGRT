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

    if(max_nroots == 0){
        GRTRaiseError("No eigenvalues were found, please check.");
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


void grt_read_cdisp(const char *filepath, EIGENV_INFO *eigmet, char **pt_modelpath)
{  
    int ncid;
    int f_dimid, n_dimid;
    int f_varid;
    int c_varid, ciref_varid, cnum_varid;

    // 打开 NC 文件
    GRTCheckFileExist(filepath);
    NC_CHECK(nc_open(filepath, NC_NOWRITE, &ncid));

    // 读取一系列全局属性
    // 读取模型路径
    {
        size_t m_len = 0;
        NC_CHECK(nc_inq_attlen(ncid, NC_GLOBAL, "model", &m_len));
        char *modelpath = (char *)calloc(m_len+1, sizeof(char));
        NC_CHECK(nc_get_att_text(ncid, NC_GLOBAL, "model", modelpath));
        modelpath[m_len] = '\0';
        if(pt_modelpath != NULL)  *pt_modelpath = modelpath;
    }
    // Rayleigh or Love
    {
        int isRayl_int = 0;
        NC_CHECK(NC_FUNC_INT(nc_get_att) (ncid, NC_GLOBAL, "isRayl", &isRayl_int));
        eigmet->wtype = (isRayl_int)? GRT_DISPERSION_RAYL : GRT_DISPERSION_LOVE;
    }
    
    // 先读取 NC 文件的坐标变量
    NC_CHECK(nc_inq_dimid(ncid, "freq", &f_dimid));
    NC_CHECK(nc_inq_dimlen(ncid, f_dimid, &eigmet->nf));
    NC_CHECK(nc_inq_dimid(ncid, "mode", &n_dimid));
    NC_CHECK(nc_inq_dimlen(ncid, n_dimid, &eigmet->nmode));
    eigmet->freqs = (real_t *)calloc(eigmet->nf, sizeof(real_t));
    NC_CHECK(nc_inq_varid(ncid, "freq", &f_varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var) (ncid, f_varid, eigmet->freqs));

    // 每个频率下有多少阶
    int *cnum = (int *)calloc(eigmet->nf, sizeof(int));
    NC_CHECK(nc_inq_varid(ncid, "cnum", &cnum_varid));
    NC_CHECK(NC_FUNC_INT(nc_get_var) (ncid, cnum_varid, cnum));

    NC_CHECK(nc_inq_varid(ncid, "c", &c_varid));
    NC_CHECK(nc_inq_varid(ncid, "ciref", &ciref_varid));    

    eigmet->eigv = (EIGENV *)calloc(eigmet->nf, sizeof(EIGENV));

    // 读入所有频散
    const int ndims = 2;
    size_t startp[ndims];
    size_t countp[ndims];
    for(size_t iw = 0; iw < eigmet->nf; ++iw){
        eigmet->eigv[iw].n = cnum[iw];
        eigmet->eigv[iw].c_roots = (real_t *)calloc(cnum[iw], sizeof(real_t));
        eigmet->eigv[iw].c_roots_iref = (size_t *)calloc(cnum[iw], sizeof(size_t));

        startp[0] = iw;
        startp[1] = 0;
        countp[0] = 1;
        countp[1] = cnum[iw];

        NC_CHECK(NC_FUNC_REAL(nc_get_vara) (ncid, c_varid, startp, countp, eigmet->eigv[iw].c_roots));
        {
            int *tmp = (int *)calloc(cnum[iw], sizeof(int));
            NC_CHECK(NC_FUNC_INT(nc_get_vara) (ncid, ciref_varid, startp, countp, tmp));
            for(int i = 0; i < cnum[iw]; ++i){
                eigmet->eigv[iw].c_roots_iref[i] = tmp[i];
            }
            GRT_SAFE_FREE_PTR(tmp);
        }
    }

    // 关闭文件
    NC_CHECK(nc_close(ncid));

    GRT_SAFE_FREE_PTR(cnum);
}


/** 根据计算好的相速度、群速度和相速度敏感核，计算群速度敏感核 */
void grt_group_sensitivity(EIGENFN_INFO *eigfnmet)
{
    // 根据群速度与相速度的关系： U = c^2 / ( c - w*dc/dw )
    // 推导得到敏感核（其中 X 表示 alpha/beta/rho ）
    //    dU/dX = U/c * (2 - U/c) * dc/dX + (U/c)^2 * w * d(dc/dX)/dw
    // 后一项存在对相速度敏感核在频率上做差分，
    // 因此输入的频率间隔大小会直接影响群速度敏感核精度

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
    int ncid, f_dimid, n_dimid;
    const int ndims = 2;
    int dimids[ndims];
    int f_varid, n_varid;
    int c_varid, u_varid, cnum_varid;

    // 创建 NC 文件
    NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

    // 定义维度
    NC_CHECK(nc_def_dim(ncid, "freq", eigfnmet->nf, &f_dimid));
    NC_CHECK(nc_def_dim(ncid, "mode", eigfnmet->nmode, &n_dimid));

    // 定义变量维度数组
    dimids[0] = f_dimid;
    dimids[1] = n_dimid;

    // 定义维度数组
    NC_CHECK(nc_def_var(ncid, "freq", NC_REAL, 1, &f_dimid, &f_varid));
    NC_CHECK(nc_def_var(ncid, "mode", NC_INT, 1, &n_dimid, &n_varid));

    // 定义频散数组
    NC_CHECK(nc_def_var(ncid, "cnum", NC_INT, 1, &f_dimid, &cnum_varid));
    NC_CHECK(nc_def_var(ncid, "c", NC_REAL, ndims, dimids, &c_varid));
    NC_CHECK(nc_def_var(ncid, "u", NC_REAL, ndims, dimids, &u_varid));

    // 填充值
    const real_t fill_value_real_t = 0.0;
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, c_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, u_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));

    // 其它属性
    // Rayleigh or Love
    int isRayl_int = (int)(eigfnmet->wtype==GRT_DISPERSION_RAYL);
    NC_CHECK(NC_FUNC_INT(nc_put_att) (ncid, NC_GLOBAL, "isRayl", NC_INT, 1, &isRayl_int));
    
    // 结束定义模式
    NC_CHECK(nc_enddef(ncid));

    // 写入数据
    NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, f_varid, eigfnmet->freqs));
    {
        int *modes = (int *)calloc(eigfnmet->nmode, sizeof(int));
        for(size_t i = 0; i < eigfnmet->nmode; ++i){
            modes[i] = eigfnmet->modes[i];
        } 
        NC_CHECK(NC_FUNC_INT(nc_put_var) (ncid, n_varid, modes));
        GRT_SAFE_FREE_PTR(modes);
    }
    {
        int *tmp = (int *)calloc(eigfnmet->nf, sizeof(int));
        for(size_t i = 0; i < eigfnmet->nf; ++i){
            tmp[i] = eigfnmet->eigv[i].n;
        }
        NC_CHECK(NC_FUNC_INT(nc_put_var) (ncid, cnum_varid, tmp));
        GRT_SAFE_FREE_PTR(tmp);
    }

    // 分批次写入数据
    size_t startp[ndims];
    size_t countp[ndims];
    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        startp[0] = iw;
        startp[1] = 0;
        countp[0] = 1;
        countp[1] = eigfnmet->eigv[iw].n;

        NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, c_varid, startp, countp, eigfnmet->eigv[iw].c_roots));
        NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, u_varid, startp, countp, eigfnmet->eigv[iw].u_roots));
    }
    
    // 关闭文件
    NC_CHECK(nc_close(ncid));

}



/* 输出本征函数结果 */
void grt_output_eigenfns(const char *filepath, const int ncols, EIGENFN_INFO *eigfnmet)
{
    int ncid, f_dimid, n_dimid, z_dimid, e_dimid;
    const int ndims = 4;
    int dimids[ndims];
    int f_varid, n_varid, z_varid;
    int c_varid, cnum_varid;
    int efn_varid;
    const int nw = 2 * ncols; // 分实部和虚部

    // 创建 NC 文件
    NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

    // 定义维度
    NC_CHECK(nc_def_dim(ncid, "freq", eigfnmet->nf, &f_dimid));
    NC_CHECK(nc_def_dim(ncid, "mode", eigfnmet->nmode, &n_dimid));
    NC_CHECK(nc_def_dim(ncid, "z", eigfnmet->nz, &z_dimid));
    NC_CHECK(nc_def_dim(ncid, "w", nw, &e_dimid));

    // 定义变量维度数组
    dimids[0] = f_dimid;
    dimids[1] = n_dimid;
    dimids[2] = z_dimid;
    dimids[3] = e_dimid;

    // 定义维度数组
    NC_CHECK(nc_def_var(ncid, "freq", NC_REAL, 1, &f_dimid, &f_varid));
    NC_CHECK(nc_def_var(ncid, "mode", NC_INT, 1, &n_dimid, &n_varid));
    NC_CHECK(nc_def_var(ncid, "z", NC_REAL, 1, &z_dimid, &z_varid));

    // 定义频散数组
    NC_CHECK(nc_def_var(ncid, "cnum", NC_INT, 1, &f_dimid, &cnum_varid));
    NC_CHECK(nc_def_var(ncid, "c", NC_REAL, 2, dimids, &c_varid));   // 取前两维

    // 填充值
    const real_t fill_value_real_t = 0.0;
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, c_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));

    // 定义本征函数的实虚部
    NC_CHECK(nc_def_var(ncid, "eigfn", NC_REAL, ndims, dimids, &efn_varid));
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, efn_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));
    
    // 其它属性
    // Rayleigh or Love
    int isRayl_int = (int)(eigfnmet->wtype==GRT_DISPERSION_RAYL);
    NC_CHECK(NC_FUNC_INT(nc_put_att) (ncid, NC_GLOBAL, "isRayl", NC_INT, 1, &isRayl_int));
    
    // 结束定义模式
    NC_CHECK(nc_enddef(ncid));

    // 写入数据
    NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, f_varid, eigfnmet->freqs));
    {
        int *modes = (int *)calloc(eigfnmet->nmode, sizeof(int));
        for(size_t i = 0; i < eigfnmet->nmode; ++i){
            modes[i] = eigfnmet->modes[i];
        } 
        NC_CHECK(NC_FUNC_INT(nc_put_var) (ncid, n_varid, modes));
        GRT_SAFE_FREE_PTR(modes);
    }
    {
        int *tmp = (int *)calloc(eigfnmet->nf, sizeof(int));
        for(size_t i = 0; i < eigfnmet->nf; ++i){
            tmp[i] = eigfnmet->eigv[i].n;
        }
        NC_CHECK(NC_FUNC_INT(nc_put_var) (ncid, cnum_varid, tmp));
        GRT_SAFE_FREE_PTR(tmp);
    }
    NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, z_varid, eigfnmet->zs));

    // 分批次写入数据
    size_t startp[ndims];
    size_t countp[ndims];

    real_t (*realimag_part)[nw] = (real_t (*)[nw])calloc(eigfnmet->nz, sizeof(real_t)*nw);

    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        startp[0] = iw;
        startp[1] = 0;
        startp[2] = 0;
        startp[3] = 0;
        
        countp[0] = 1;
        countp[1] = eigfnmet->eigv[iw].n;
        countp[2] = eigfnmet->nz;
        countp[3] = nw;

        NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, c_varid, startp, countp, eigfnmet->eigv[iw].c_roots));

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

            startp[1] = ic;
            countp[1] = 1;
            NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, efn_varid, startp, countp, realimag_part[0]));
        }
    }
    
    // 关闭文件
    NC_CHECK(nc_close(ncid));

    GRT_SAFE_FREE_PTR(realimag_part);
}


/* 输出能量积分结果 */
void grt_output_energy_integrals(const char *filepath, EIGENFN_INFO *eigfnmet)
{
    int ncid, f_dimid, n_dimid, e_dimid;
    const int ndims = 3;
    int dimids[ndims];
    int f_varid, n_varid, eint_varid;
    int c_varid, cnum_varid;
    int ne = 2 * (GRT_EGYINTS_MAX + 1);  // + 1 为了验证项

    // 创建 NC 文件
    NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

    // 定义维度
    NC_CHECK(nc_def_dim(ncid, "freq", eigfnmet->nf, &f_dimid));
    NC_CHECK(nc_def_dim(ncid, "mode", eigfnmet->nmode, &n_dimid));
    NC_CHECK(nc_def_dim(ncid, "e", ne, &e_dimid));

    // 定义变量维度数组
    dimids[0] = f_dimid;
    dimids[1] = n_dimid;
    dimids[2] = e_dimid;

    // 定义维度数组
    NC_CHECK(nc_def_var(ncid, "freq", NC_REAL, 1, &f_dimid, &f_varid));
    NC_CHECK(nc_def_var(ncid, "mode", NC_INT, 1, &n_dimid, &n_varid));

    // 定义频散数组
    NC_CHECK(nc_def_var(ncid, "cnum", NC_INT, 1, &f_dimid, &cnum_varid));
    NC_CHECK(nc_def_var(ncid, "c", NC_REAL, 2, dimids, &c_varid));

    // 填充值
    const real_t fill_value_real_t = 0.0;
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, c_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));

    // 定义能量积分
    NC_CHECK(nc_def_var(ncid, "egyint", NC_REAL, ndims, dimids, &eint_varid));
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, eint_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));

    // 其它属性
    // Rayleigh or Love
    int isRayl_int = (int)(eigfnmet->wtype==GRT_DISPERSION_RAYL);
    NC_CHECK(NC_FUNC_INT(nc_put_att) (ncid, NC_GLOBAL, "isRayl", NC_INT, 1, &isRayl_int));
    
    // 结束定义模式
    NC_CHECK(nc_enddef(ncid));

    // 写入数据
    NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, f_varid, eigfnmet->freqs));
    {
        int *modes = (int *)calloc(eigfnmet->nmode, sizeof(int));
        for(size_t i = 0; i < eigfnmet->nmode; ++i){
            modes[i] = eigfnmet->modes[i];
        } 
        NC_CHECK(NC_FUNC_INT(nc_put_var) (ncid, n_varid, modes));
        GRT_SAFE_FREE_PTR(modes);
    }
    {
        int *tmp = (int *)calloc(eigfnmet->nf, sizeof(int));
        for(size_t i = 0; i < eigfnmet->nf; ++i){
            tmp[i] = eigfnmet->eigv[i].n;
        }
        NC_CHECK(NC_FUNC_INT(nc_put_var) (ncid, cnum_varid, tmp));
        GRT_SAFE_FREE_PTR(tmp);
    }

    // 分批次写入数据
    size_t startp[ndims];
    size_t countp[ndims];

    real_t (*realimag_part)[ne] = (real_t (*)[ne])calloc(eigfnmet->nmode, sizeof(real_t)*ne);

    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        startp[0] = iw;
        startp[1] = 0;
        startp[2] = 0;
        
        countp[0] = 1;
        countp[1] = eigfnmet->eigv[iw].n;
        countp[2] = ne;

        real_t omega = PI2*eigfnmet->freqs[iw];
        real_t *c_roots = eigfnmet->eigv[iw].c_roots;

        NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, c_varid, startp, countp, c_roots));

        memset(realimag_part, 0, sizeof(real_t)*eigfnmet->nmode*ne);
        for(size_t ic = 0; ic < eigfnmet->eigv[iw].n; ++ic){
            cplx_t *egyint = eigfnmet->eigfn[iw][ic].egyint;

            for(int j=0; j < GRT_EGYINTS_MAX; ++j){
                realimag_part[ic][2*j]     = creal(egyint[j]);
                realimag_part[ic][2*j + 1] = cimag(egyint[j]);
            }

            // 加上验证项
            cplx_t res=-12345.0;
            real_t k = omega/c_roots[ic];
            res = GRT_SQUARE(omega) * egyint[0] - k*k*egyint[1] + k*egyint[2] - egyint[3];
            realimag_part[ic][ne-2] = creal(res);
            realimag_part[ic][ne-1] = cimag(res);
        }

        NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, eint_varid, startp, countp, realimag_part[0]));
    }

    // 关闭文件
    NC_CHECK(nc_close(ncid));

    GRT_SAFE_FREE_PTR(realimag_part);

}


void grt_output_sensitivity(const char *filepath, const char *UC, EIGENFN_INFO *eigfnmet)
{
    int ncid, f_dimid, n_dimid, z_dimid, k_dimid;
    const int ndims = 4;
    int dimids[ndims];
    int f_varid, n_varid, z_varid, sens_varid;
    int cg_varid, cnum_varid;
    int nk = GRT_SNSTVTY_MAX;  // 敏感核仅取实部

    // 创建 NC 文件
    NC_CHECK(nc_create(filepath, NC_CLOBBER, &ncid));

    // 定义维度
    NC_CHECK(nc_def_dim(ncid, "freq", eigfnmet->nf, &f_dimid));
    NC_CHECK(nc_def_dim(ncid, "mode", eigfnmet->nmode, &n_dimid));
    NC_CHECK(nc_def_dim(ncid, "z", eigfnmet->cpar_nz, &z_dimid));
    NC_CHECK(nc_def_dim(ncid, "k", nk, &k_dimid));

    // 定义变量维度数组
    dimids[0] = f_dimid;
    dimids[1] = n_dimid;
    dimids[2] = z_dimid;
    dimids[3] = k_dimid;

    // 定义维度数组
    NC_CHECK(nc_def_var(ncid, "freq", NC_REAL, 1, &f_dimid, &f_varid));
    NC_CHECK(nc_def_var(ncid, "mode", NC_INT, 1, &n_dimid, &n_varid));
    NC_CHECK(nc_def_var(ncid, "z", NC_REAL, 1, &z_dimid, &z_varid));

    // 定义频散数组
    NC_CHECK(nc_def_var(ncid, "cnum", NC_INT, 1, &f_dimid, &cnum_varid));
    NC_CHECK(nc_def_var(ncid, UC, NC_REAL, 2, dimids, &cg_varid));   // 取前两维

    // 填充值
    const real_t fill_value_real_t = 0.0;
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, cg_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));

    // 敏感核（仅取实部）
    char *sname = NULL;
    GRT_SAFE_ASPRINTF(&sname, "%ssens", UC);
    NC_CHECK(nc_def_var(ncid, sname, NC_REAL, ndims, dimids, &sens_varid));
    GRT_SAFE_FREE_PTR(sname);
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, sens_varid, "_FillValue", NC_REAL, 1, &fill_value_real_t));

    // 其它属性
    // Rayleigh or Love
    int isRayl_int = (int)(eigfnmet->wtype==GRT_DISPERSION_RAYL);
    NC_CHECK(NC_FUNC_INT(nc_put_att) (ncid, NC_GLOBAL, "isRayl", NC_INT, 1, &isRayl_int));
    
    // 结束定义模式
    NC_CHECK(nc_enddef(ncid));

    // 写入数据
    NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, f_varid, eigfnmet->freqs));
    {
        int *modes = (int *)calloc(eigfnmet->nmode, sizeof(int));
        for(size_t i = 0; i < eigfnmet->nmode; ++i){
            modes[i] = eigfnmet->modes[i];
        } 
        NC_CHECK(NC_FUNC_INT(nc_put_var) (ncid, n_varid, modes));
        GRT_SAFE_FREE_PTR(modes);
    }
    {
        int *tmp = (int *)calloc(eigfnmet->nf, sizeof(int));
        for(size_t i = 0; i < eigfnmet->nf; ++i){
            tmp[i] = eigfnmet->eigv[i].n;
        }
        NC_CHECK(NC_FUNC_INT(nc_put_var) (ncid, cnum_varid, tmp));
        GRT_SAFE_FREE_PTR(tmp);
    }
    NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, z_varid, eigfnmet->cpar_zs));

    // 分批次写入数据
    size_t startp[ndims];
    size_t countp[ndims];

    real_t (*real_part)[nk] = (real_t (*)[nk])calloc(eigfnmet->cpar_nz, sizeof(real_t)*nk);

    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        startp[0] = iw;
        startp[1] = 0;
        startp[2] = 0;
        startp[3] = 0;
        
        countp[0] = 1;
        countp[1] = eigfnmet->eigv[iw].n;
        countp[2] = eigfnmet->cpar_nz;
        countp[3] = nk;

        real_t *cudisp = NULL;
        if(strcmp(UC, "c") == 0){
            cudisp = eigfnmet->eigv[iw].c_roots;
        }
        else if(strcmp(UC, "u") == 0){
            cudisp = eigfnmet->eigv[iw].u_roots;
        }
        else {
            GRTRaiseError("Wrong execution.");
        }

        NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, cg_varid, startp, countp, cudisp));

        for(size_t ic = 0; ic < eigfnmet->eigv[iw].n; ++ic){
            memset(real_part, 0, sizeof(real_t)*eigfnmet->cpar_nz*nk);

            cplx_t (*cusens)[GRT_SNSTVTY_MAX] = NULL;
            
            if(strcmp(UC, "c") == 0){
                cusens = eigfnmet->eigfn[iw][ic].csens;
            }
            else if(strcmp(UC, "u") == 0){
                cusens = eigfnmet->eigfn[iw][ic].usens;
            }
            else {
                GRTRaiseError("Wrong execution.");
            }

            // 将复数数据转为实部和虚部
            for(size_t iz = 0; iz < eigfnmet->cpar_nz; ++iz){
                for(int j = 0; j < nk; ++j){
                    real_part[iz][j] = creal(cusens[iz][j]);
                }
            }

            startp[1] = ic;
            countp[1] = 1;
            NC_CHECK(NC_FUNC_REAL(nc_put_vara) (ncid, sens_varid, startp, countp, real_part[0]));
        }
    }

    // 关闭文件
    NC_CHECK(nc_close(ncid));

    GRT_SAFE_FREE_PTR(real_part);
}
