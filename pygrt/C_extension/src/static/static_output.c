/**
 * @file   static_output.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * 静态位移模块的公共 NetCDF 输出
 *
 */

#include "grt/static/static_output.h"

#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grt/common/mynetcdf.h"

/**
 * 定义位移变量及可选的位移偏导变量
 *
 * @param[in]   ncid         NetCDF 文件 ID
 * @param[in]   ndims        变量维度数
 * @param[in]   dimids       变量维度 ID 数组
 * @param[in]   channels     输出分量编码
 * @param[in]   calc_upar    是否定义位移偏导变量
 * @param[in]   suffix       变量名使用的断层编号后缀，0 表示不添加后缀
 * @param[out]  vars         位移变量 ID 数组
 * @param[out]  dvars        位移偏导变量 ID 数组
 */
static void define_channel_vars(
    int ncid, int ndims, const int *dimids, const char *channels,
    bool calc_upar, size_t suffix,
    int vars[GRT_CHANNEL_NUM], int dvars[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    // 先定义位移变量，再按需定义每个导数方向对应的偏导变量
    for(int c = 0; c < GRT_CHANNEL_NUM; ++c){
        char name[32];
        if(suffix == 0){
            snprintf(name, sizeof(name), "%c", toupper(channels[c]));
        } else {
            snprintf(name, sizeof(name), "%c%zu", toupper(channels[c]), suffix);
        }
        NC_CHECK(nc_def_var(ncid, name, NC_REAL, ndims, dimids, &vars[c]));
        if(calc_upar){
            // 偏导变量按导数方向和位移分量分别组织
            for(int d = 0; d < GRT_CHANNEL_NUM; ++d){
                if(suffix == 0){
                    snprintf(name, sizeof(name), "%c%c", tolower(channels[d]), toupper(channels[c]));
                } else {
                    snprintf(name, sizeof(name), "%c%c%zu",
                        tolower(channels[d]), toupper(channels[c]), suffix);
                }
                NC_CHECK(nc_def_var(ncid, name, NC_REAL, ndims, dimids, &dvars[d][c]));
            }
        }
    }
}

/**
 * 写入位移变量及可选的位移偏导变量
 *
 * @param[in]  ncid         NetCDF 文件 ID
 * @param[in]  npts         当前布局中的接收点数
 * @param[in]  calc_upar    是否写入位移偏导变量
 * @param[in]  syn          位移数组
 * @param[in]  syn_upar     位移偏导数组
 * @param[in]  vars         位移变量 ID 数组
 * @param[in]  dvars        位移偏导变量 ID 数组
 */
static void write_fields(
    int ncid, size_t npts, bool calc_upar,
    const real_t (*syn)[GRT_CHANNEL_NUM],
    const real_t (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    const int vars[GRT_CHANNEL_NUM],
    const int dvars[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    real_t *buffer = (real_t *)calloc(npts, sizeof(real_t));
    // NetCDF 变量按分量写入，临时数组用于提取每个分量的点序列
    for(int c = 0; c < GRT_CHANNEL_NUM; ++c){
        for(size_t i = 0; i < npts; ++i){
            buffer[i] = syn[i][c];
        }
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, vars[c], buffer));
        if(calc_upar){
            // 位移偏导按接收点、导数方向和位移分量的顺序写入
            for(int d = 0; d < GRT_CHANNEL_NUM; ++d){
                for(size_t i = 0; i < npts; ++i){
                    buffer[i] = syn_upar[i][d][c];
                }
                NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, dvars[d][c], buffer));
            }
        }
    }
    GRT_SAFE_FREE_PTR(buffer);
}

/**
 * 定义并写入规则网格布局
 *
 * @param[in]      ncid      NetCDF 文件 ID
 * @param[in]      output    静态位移输出描述
 * @param[out]     vars      位移变量 ID 数组
 * @param[out]     dvars     位移偏导变量 ID 数组
 */
static void write_grid_layout(
    int ncid, const GRT_STATIC_NC_OUTPUT *output,
    int vars[GRT_CHANNEL_NUM], int dvars[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    const GRT_RECV_POINTS *recv = output->recv;
    int dimids[2];
    int north_varid, east_varid;
    real_t deprcv = recv->depths[0];
    real_t rcv_va, rcv_vb, rcv_rho;

    // 规则网格使用统一接收深度，将对应介质参数保存为全局属性
    output->get_medium(
        output->medium_context, deprcv, &rcv_va, &rcv_vb, &rcv_rho);
    NC_CHECK(NC_FUNC_REAL(nc_put_att)(ncid, NC_GLOBAL, "deprcv", NC_REAL, 1, &deprcv));
    NC_CHECK(NC_FUNC_REAL(nc_put_att)(ncid, NC_GLOBAL, "rcv_va", NC_REAL, 1, &rcv_va));
    NC_CHECK(NC_FUNC_REAL(nc_put_att)(ncid, NC_GLOBAL, "rcv_vb", NC_REAL, 1, &rcv_vb));
    NC_CHECK(NC_FUNC_REAL(nc_put_att)(ncid, NC_GLOBAL, "rcv_rho", NC_REAL, 1, &rcv_rho));

    NC_CHECK(nc_def_dim(ncid, "north", recv->nnorth, &dimids[0]));
    NC_CHECK(nc_def_dim(ncid, "east", recv->neast, &dimids[1]));
    NC_CHECK(nc_def_var(ncid, "north", NC_REAL, 1, &dimids[0], &north_varid));
    NC_CHECK(nc_def_var(ncid, "east", NC_REAL, 1, &dimids[1], &east_varid));
    define_channel_vars(
        ncid, 2, dimids, output->channels, output->calc_upar, 0, vars, dvars);
    // 结束定义模式后写入坐标轴和位移数据
    NC_CHECK(nc_enddef(ncid));

    real_t *north_axis = (real_t *)calloc(recv->nnorth, sizeof(real_t));
    real_t *east_axis = (real_t *)calloc(recv->neast, sizeof(real_t));
    // 接收点坐标在内存中按 point 展平，写网格文件时恢复为两个坐标轴
    for(size_t inorth = 0; inorth < recv->nnorth; ++inorth){
        north_axis[inorth] = recv->norths[inorth * recv->neast];
    }
    for(size_t ieast = 0; ieast < recv->neast; ++ieast){
        east_axis[ieast] = recv->easts[ieast];
    }
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, north_varid, north_axis));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, east_varid, east_axis));
    GRT_SAFE_FREE_PTR(north_axis);
    GRT_SAFE_FREE_PTR(east_axis);
    write_fields(ncid, recv->npts, output->calc_upar, output->syn, output->syn_upar, vars, dvars);
}

/**
 * 定义并写入一维接收点布局
 *
 * @param[in]      ncid      NetCDF 文件 ID
 * @param[in]      output    静态位移输出描述
 * @param[out]     vars      位移变量 ID 数组
 * @param[out]     dvars     位移偏导变量 ID 数组
 */
static void write_points_layout(
    int ncid, const GRT_STATIC_NC_OUTPUT *output,
    int vars[GRT_CHANNEL_NUM], int dvars[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    const GRT_RECV_POINTS *recv = output->recv;
    // 有限接收断层已经在 GRT_RECV_POINTS 中展开为一维点数组
    size_t npts = recv->npts;
    const real_t *norths = recv->norths;
    const real_t *easts = recv->easts;
    const real_t *depths = recv->depths;
    int point_dimid;
    int point_dimids[1];
    int north_varid, east_varid, depth_varid;
    int va_varid, vb_varid, rho_varid;
    int strike_varid = -1, dip_varid = -1, rake_varid = -1;
    int offset_varid = -1, stksize_varid = -1, dipsize_varid = -1;

    // 所有接收点变量使用 point 维，有限接收断层再附加 nfault 维
    NC_CHECK(nc_def_dim(ncid, "point", npts, &point_dimid));
    point_dimids[0] = point_dimid;
    NC_CHECK(nc_def_var(ncid, "north", NC_REAL, 1, point_dimids, &north_varid));
    NC_CHECK(nc_def_var(ncid, "east", NC_REAL, 1, point_dimids, &east_varid));
    NC_CHECK(nc_def_var(ncid, "depth", NC_REAL, 1, point_dimids, &depth_varid));
    NC_CHECK(nc_def_var(ncid, "rcv_va", NC_REAL, 1, point_dimids, &va_varid));
    NC_CHECK(nc_def_var(ncid, "rcv_vb", NC_REAL, 1, point_dimids, &vb_varid));
    NC_CHECK(nc_def_var(ncid, "rcv_rho", NC_REAL, 1, point_dimids, &rho_varid));
    if(recv->is_fault){
        int fault_dimid;
        NC_CHECK(nc_def_dim(ncid, "nfault", recv->nfault, &fault_dimid));
        NC_CHECK(nc_def_var(ncid, "strike", NC_REAL, 1, &fault_dimid, &strike_varid));
        NC_CHECK(nc_def_var(ncid, "dip", NC_REAL, 1, &fault_dimid, &dip_varid));
        NC_CHECK(nc_def_var(ncid, "rake", NC_REAL, 1, &fault_dimid, &rake_varid));
        NC_CHECK(nc_def_var(ncid, "offset", NC_INT, 1, &fault_dimid, &offset_varid));
        NC_CHECK(nc_def_var(ncid, "stksize", NC_INT, 1, &fault_dimid, &stksize_varid));
        NC_CHECK(nc_def_var(ncid, "dipsize", NC_INT, 1, &fault_dimid, &dipsize_varid));
    } else if(recv->has_geometry){
        NC_CHECK(nc_def_var(ncid, "strike", NC_REAL, 1, point_dimids, &strike_varid));
        NC_CHECK(nc_def_var(ncid, "dip", NC_REAL, 1, point_dimids, &dip_varid));
        NC_CHECK(nc_def_var(ncid, "rake", NC_REAL, 1, point_dimids, &rake_varid));
    }
    define_channel_vars(
        ncid, 1, point_dimids, output->channels, output->calc_upar, 0, vars, dvars);
    // 完成变量定义后，按点查询接收介质参数并准备写入
    NC_CHECK(nc_enddef(ncid));

    real_t *rcv_va = (real_t *)calloc(npts, sizeof(real_t));
    real_t *rcv_vb = (real_t *)calloc(npts, sizeof(real_t));
    real_t *rcv_rho = (real_t *)calloc(npts, sizeof(real_t));
    for(size_t i = 0; i < npts; ++i){
        output->get_medium(
            output->medium_context, depths[i], &rcv_va[i], &rcv_vb[i], &rcv_rho[i]);
    }

    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, north_varid, norths));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, east_varid, easts));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, depth_varid, depths));
    // 接收介质参数与坐标使用相同的 point 顺序写入
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, va_varid, rcv_va));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, vb_varid, rcv_vb));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, rho_varid, rcv_rho));
    if(recv->is_fault){
        // offset 保存每条有限断层点范围的排他性结束索引，最后一个值等于 point
        int *offsets = (int *)calloc(recv->nfault, sizeof(int));
        int *stksizes = (int *)calloc(recv->nfault, sizeof(int));
        int *dipsizes = (int *)calloc(recv->nfault, sizeof(int));
        for(size_t ifault = 0; ifault < recv->nfault; ++ifault){
            if(recv->offsets[ifault] > (size_t)INT_MAX){
                GRTRaiseError("receiver point offset exceeds the NetCDF integer range.");
            }
            if((recv->stksizes[ifault] > (size_t)INT_MAX) ||
                (recv->dipsizes[ifault] > (size_t)INT_MAX)){
                GRTRaiseError("receiver fault subdivision count exceeds the NetCDF integer range.");
            }
            offsets[ifault] = (int)recv->offsets[ifault];
            stksizes[ifault] = (int)recv->stksizes[ifault];
            dipsizes[ifault] = (int)recv->dipsizes[ifault];
        }
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, strike_varid, recv->fstrikes));
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, dip_varid, recv->fdips));
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, rake_varid, recv->frakes));
        NC_CHECK(nc_put_var_int(ncid, offset_varid, offsets));
        NC_CHECK(nc_put_var_int(ncid, stksize_varid, stksizes));
        NC_CHECK(nc_put_var_int(ncid, dipsize_varid, dipsizes));
        GRT_SAFE_FREE_PTR(offsets);
        GRT_SAFE_FREE_PTR(stksizes);
        GRT_SAFE_FREE_PTR(dipsizes);
    } else if(recv->has_geometry){
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, strike_varid, recv->strikes));
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, dip_varid, recv->dips));
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, rake_varid, recv->rakes));
    }
    GRT_SAFE_FREE_PTR(rcv_va);
    GRT_SAFE_FREE_PTR(rcv_vb);
    GRT_SAFE_FREE_PTR(rcv_rho);
    // 坐标和附加属性写入后，再写入位移及其偏导结果
    write_fields(ncid, npts, output->calc_upar, output->syn, output->syn_upar, vars, dvars);
}

/**
 * 写入公共静态 NetCDF 属性和选定的接收点布局
 *
 * @param[in]  output   静态位移输出描述
 */
void grt_static_save_nc(const GRT_STATIC_NC_OUTPUT *output)
{
    if((output == NULL) || (output->path == NULL)){
        GRTRaiseError("static NetCDF output description is incomplete.");
    }
    if(output->recv == NULL){
        GRTRaiseError("static NetCDF output has no receiver layout.");
    }
    if(output->get_medium == NULL){
        GRTRaiseError("static NetCDF output has no receiver-medium callback.");
    }

    // 有限接收断层和任意接收点使用一维布局，规则网格使用二维布局
    const char *layout;
    if(output->recv->is_grid && !output->recv->is_fault){
        layout = GRT_RECV_LAYOUT_GRID;
    } else {
        layout = GRT_RECV_LAYOUT_POINTS;
    }

    int ncid;
    // 创建文件并写入所有布局共用的全局属性
    NC_CHECK(nc_create(output->path, NC_CLOBBER, &ncid));
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "layout", strlen(layout), layout));
    NC_CHECK(nc_put_att_text(
        ncid, NC_GLOBAL, "computeType",
        strlen(output->compute_type), output->compute_type));
    if(output->coordinate != NULL){
        NC_CHECK(nc_put_att_text(
            ncid, NC_GLOBAL, "coordinate",
            strlen(output->coordinate), output->coordinate));
    }
    if(output->has_elastic_params){
        NC_CHECK(NC_FUNC_REAL(nc_put_att)(
            ncid, NC_GLOBAL, "alpha", NC_REAL, 1, &output->alpha));
        NC_CHECK(NC_FUNC_REAL(nc_put_att)(
            ncid, NC_GLOBAL, "lambda", NC_REAL, 1, &output->lambda));
        NC_CHECK(NC_FUNC_REAL(nc_put_att)(
            ncid, NC_GLOBAL, "mu", NC_REAL, 1, &output->mu));
    }
    if(output->has_depsrc){
        NC_CHECK(NC_FUNC_REAL(nc_put_att)(
            ncid, NC_GLOBAL, "depsrc", NC_REAL, 1, &output->depsrc));
    }
    int calc_upar = (output->calc_upar) ? 1 : 0;
    int rot2ZNE = (output->rot2ZNE) ? 1 : 0;
    NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "calc_upar", NC_INT, 1, &calc_upar));
    NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "rot2ZNE", NC_INT, 1, &rot2ZNE));

    int vars[GRT_CHANNEL_NUM];
    int dvars[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
    // 根据接收点布局定义维度、坐标变量和结果变量并写入数据
    if(output->recv->is_grid && !output->recv->is_fault){
        write_grid_layout(ncid, output, vars, dvars);
    } else {
        write_points_layout(ncid, output, vars, dvars);
    }
    NC_CHECK(nc_close(ncid));
}
