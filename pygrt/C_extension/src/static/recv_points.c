/**
 * @file   recv_points.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * 静态 syn / 后处理用的接收点列表
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "grt/static/recv_points.h"
#include "grt/common/checkerror.h"
#include "grt/common/util.h"
#include "grt/common/mynetcdf.h"

/**
 * 解析一行接收点数据并返回有效数值列数
 *
 * @param[in]   line       待解析的文本行
 * @param[out]  values     解析出的数值数组
 * @param[out]  nvalues    解析出的数值个数
 * @return                 true 表示列数有效，false 表示格式错误
 */
static bool parse_receiver_point_line(
    const char *line, real_t values[6], size_t *nvalues)
{
    const char *cursor = line;
    *nvalues = 0;

    while(true){
        while(isspace((unsigned char)*cursor)) ++cursor;
        if(*cursor == '\0' || *cursor == GRT_COMMENT_HEAD) break;
        if(*nvalues >= 6) return false;

        errno = 0;
        char *end = NULL;
        real_t value = strtod(cursor, &end);
        if(end == cursor || errno == ERANGE || !isfinite(value)) return false;
        values[*nvalues] = value;
        (*nvalues)++;
        cursor = end;
    }

    return *nvalues == 3 || *nvalues == 6;
}

GRT_RECV_POINTS *grt_recv_points_from_grid(
    size_t nnorth, const real_t *norths,
    size_t neast,  const real_t *easts,
    real_t depth)
{
    if(nnorth == 0 || neast == 0 || norths == NULL || easts == NULL){
        GRTRaiseError("empty receiver grid.");
    }
    if(depth < 0.0){
        GRTRaiseError("Negative receiver depth is not supported.");
    }

    GRT_RECV_POINTS *pts = (GRT_RECV_POINTS *)calloc(1, sizeof(GRT_RECV_POINTS));
    pts->is_grid = true;
    pts->nnorth = nnorth;
    pts->neast = neast;
    pts->npts = nnorth * neast;
    pts->norths = (real_t *)calloc(pts->npts, sizeof(real_t));
    pts->easts  = (real_t *)calloc(pts->npts, sizeof(real_t));
    pts->depths = (real_t *)calloc(pts->npts, sizeof(real_t));

    for(size_t inorth = 0; inorth < nnorth; ++inorth){
        for(size_t ieast = 0; ieast < neast; ++ieast){
            size_t ipt = ieast + inorth * neast;
            pts->norths[ipt] = norths[inorth];
            pts->easts[ipt]  = easts[ieast];
            pts->depths[ipt] = depth;
        }
    }
    return pts;
}


GRT_RECV_POINTS *grt_recv_points_from_file(const char *path)
{
    GRTCheckFileExist(path);

    FILE *fp = fopen(path, "r");
    if(fp == NULL){
        GRTRaiseError("Failed to open receiver points file \"%s\".", path);
    }

    // 先统计有效行数并确定文件列数
    size_t npts = 0;
    size_t ncolumns = 0;
    char *line = NULL;
    size_t nlen = 0;
    size_t lineno = 0;
    while(grt_getline(&line, &nlen, fp) != -1){
        lineno++;
        grt_trim_whitespace(line);
        if(grt_is_comment_or_empty(line)) continue;

        real_t values[6];
        size_t nvalues;
        if(!parse_receiver_point_line(line, values, &nvalues)){
            GRT_SAFE_FREE_PTR(line);
            fclose(fp);
            GRTRaiseError(
                "Invalid receiver point at line %zu in \"%s\" "
                "(expect exactly 3 or 6 numeric columns: "
                "north east depth [strike dip rake]).",
                lineno, path);
        }
        if(ncolumns != 0 && nvalues != ncolumns){
            GRT_SAFE_FREE_PTR(line);
            fclose(fp);
            GRTRaiseError(
                "Inconsistent receiver point column count at line %zu in \"%s\" "
                "(all data lines must have either 3 or 6 columns).",
                lineno, path);
        }
        ncolumns = nvalues;
        npts++;
    }
    if(npts == 0){
        GRT_SAFE_FREE_PTR(line);
        fclose(fp);
        GRTRaiseError("No receiver points found in \"%s\".", path);
    }

    GRT_RECV_POINTS *pts = (GRT_RECV_POINTS *)calloc(1, sizeof(GRT_RECV_POINTS));
    pts->is_grid = false;
    pts->nnorth = 0;
    pts->neast = 0;
    pts->npts = npts;
    pts->norths = (real_t *)calloc(npts, sizeof(real_t));
    pts->easts  = (real_t *)calloc(npts, sizeof(real_t));
    pts->depths = (real_t *)calloc(npts, sizeof(real_t));
    pts->has_geometry = ncolumns == 6;
    if(pts->has_geometry){
        pts->strikes = (real_t *)calloc(npts, sizeof(real_t));
        pts->dips    = (real_t *)calloc(npts, sizeof(real_t));
        pts->rakes   = (real_t *)calloc(npts, sizeof(real_t));
    }

    rewind(fp);
    size_t ipt = 0;
    lineno = 0;
    while(grt_getline(&line, &nlen, fp) != -1){
        lineno++;
        grt_trim_whitespace(line);
        if(grt_is_comment_or_empty(line)) continue;

        real_t values[6];
        size_t nvalues;
        if(!parse_receiver_point_line(line, values, &nvalues) || nvalues != ncolumns){
            GRT_SAFE_FREE_PTR(line);
            grt_recv_points_free(pts);
            fclose(fp);
            GRTRaiseError(
                "Invalid receiver point at line %zu in \"%s\" "
                "(expect the same 3 or 6 columns as the other data lines).",
                lineno, path);
        }
        if(values[2] < 0.0){
            GRT_SAFE_FREE_PTR(line);
            grt_recv_points_free(pts);
            fclose(fp);
            GRTRaiseError("Negative receiver depth at line %zu in \"%s\".", lineno, path);
        }
        pts->norths[ipt] = values[0];
        pts->easts[ipt]  = values[1];
        pts->depths[ipt] = values[2];
        if(pts->has_geometry){
            pts->strikes[ipt] = values[3];
            pts->dips[ipt]    = values[4];
            pts->rakes[ipt]   = values[5];
        }
        ipt++;
    }

    GRT_SAFE_FREE_PTR(line);
    fclose(fp);
    return pts;
}


void grt_recv_points_free(GRT_RECV_POINTS *pts)
{
    if(pts == NULL) return;
    // 释放接收点坐标和任意点的逐点几何
    GRT_SAFE_FREE_PTR(pts->norths);
    GRT_SAFE_FREE_PTR(pts->easts);
    GRT_SAFE_FREE_PTR(pts->depths);
    GRT_SAFE_FREE_PTR(pts->strikes);
    GRT_SAFE_FREE_PTR(pts->dips);
    GRT_SAFE_FREE_PTR(pts->rakes);
    // 释放有限接收断层的索引、断层几何和子断层尺度
    GRT_SAFE_FREE_PTR(pts->nsubs);
    GRT_SAFE_FREE_PTR(pts->offsets);
    GRT_SAFE_FREE_PTR(pts->fstrikes);
    GRT_SAFE_FREE_PTR(pts->fdips);
    GRT_SAFE_FREE_PTR(pts->frakes);
    GRT_SAFE_FREE_PTR(pts->stksizes);
    GRT_SAFE_FREE_PTR(pts->dipsizes);
    GRT_SAFE_FREE_PTR(pts);
}


GRT_RECV_POINTS *grt_recv_points_from_faults(
    size_t nfault, const FINITE_FAULT *faults, real_t dL, real_t dW)
{
    if((nfault == 0) || (faults == NULL)){
        GRTRaiseError("empty finite receiver faults.");
    }
    if((dL <= 0.0) != (dW <= 0.0)){
        GRTRaiseError("finite receiver dL and dW must both be positive or both be omitted.");
    }

    GRT_RECV_POINTS *pts = (GRT_RECV_POINTS *)calloc(1, sizeof(*pts));
    pts->is_fault = true;
    pts->nfault = nfault;
    pts->nsubs = (size_t *)calloc(nfault, sizeof(*pts->nsubs));
    pts->offsets = (size_t *)calloc(nfault, sizeof(*pts->offsets));
    pts->fstrikes = (real_t *)calloc(nfault, sizeof(*pts->fstrikes));
    pts->fdips = (real_t *)calloc(nfault, sizeof(*pts->fdips));
    pts->frakes = (real_t *)calloc(nfault, sizeof(*pts->frakes));
    pts->stksizes = (size_t *)calloc(nfault, sizeof(*pts->stksizes));
    pts->dipsizes = (size_t *)calloc(nfault, sizeof(*pts->dipsizes));

    // 遍历每条有限断层，统一调用公共几何函数并追加子断层中心点
    for(size_t ifault = 0; ifault < nfault; ++ifault){
        const FINITE_FAULT *fault = &faults[ifault];
        real_t W, L;
        size_t nW, nL;
        grt_finite_fault_subdiv(fault, dL, dW, &W, &L, &nW, &nL);
        size_t nsubs = nW * nL;
        size_t first_point = pts->npts;
        size_t new_npts = first_point + nsubs;
        pts->nsubs[ifault] = nsubs;
        pts->offsets[ifault] = new_npts;
        pts->fstrikes[ifault] = fault->strike;
        pts->fdips[ifault] = fault->dip;
        pts->frakes[ifault] = fault->rake;
        pts->stksizes[ifault] = nL;
        pts->dipsizes[ifault] = nW;
        pts->norths = (real_t *)realloc(pts->norths, new_npts * sizeof(*pts->norths));
        pts->easts = (real_t *)realloc(pts->easts, new_npts * sizeof(*pts->easts));
        pts->depths = (real_t *)realloc(pts->depths, new_npts * sizeof(*pts->depths));
        pts->npts = new_npts;

        // 将当前有限断层的子断层中心按 point 顺序追加
        size_t ipt = first_point;
        for(size_t iW = 0; iW < nW; ++iW){
            for(size_t iL = 0; iL < nL; ++iL){
                FINITE_SUBFAULT sub;
                grt_finite_fault_subfault(fault, dL, dW, W, L, iW, iL, &sub);
                pts->norths[ipt] = sub.north;
                pts->easts[ipt] = sub.east;
                pts->depths[ipt] = sub.depsrc;
                ++ipt;
            }
        }
    }
    return pts;
}


GRT_RECV_NC_LAYOUT grt_recv_nc_get_layout(int ncid)
{
    size_t len = 0;
    int status = nc_inq_attlen(ncid, NC_GLOBAL, "layout", &len);
    if((status != NC_NOERR) || (len == 0)){
        GRTRaiseError("static receiver layout attribute is missing.");
    }

    // 所有静态输出文件都显式保存 layout 属性，直接按属性确定布局
    char *layout = (char *)calloc(len + 1, 1);
    NC_CHECK(nc_get_att_text(ncid, NC_GLOBAL, "layout", layout));

    GRT_RECV_NC_LAYOUT result;
    if(strcmp(layout, GRT_RECV_LAYOUT_GRID) == 0){
        result = GRT_RECV_NC_LAYOUT_GRID;
    } else if(strcmp(layout, GRT_RECV_LAYOUT_POINTS) == 0){
        result = GRT_RECV_NC_LAYOUT_POINTS;
    } else {
        GRTRaiseError("unsupported static receiver layout \"%s\".", layout);
    }
    GRT_SAFE_FREE_PTR(layout);
    return result;
}


void grt_recv_nc_info_load(int ncid, GRT_RECV_NC_INFO *info)
{
    if(info == NULL){
        GRTRaiseError("receiver NetCDF info is NULL.");
    }
    memset(info, 0, sizeof(*info));
    // 先确定文件布局，再按布局读取并组织接收坐标
    info->layout = grt_recv_nc_get_layout(ncid);

    if(info->layout == GRT_RECV_NC_LAYOUT_GRID){
        size_t nnorth, neast;
        int north_varid, east_varid;
        NC_CHECK(nc_inq_dimid(ncid, "north", &info->dimids[0]));
        NC_CHECK(nc_inq_dimlen(ncid, info->dimids[0], &nnorth));
        NC_CHECK(nc_inq_dimid(ncid, "east", &info->dimids[1]));
        NC_CHECK(nc_inq_dimlen(ncid, info->dimids[1], &neast));
        info->npts = nnorth * neast;
        info->norths = (real_t *)calloc(info->npts, sizeof(real_t));
        info->easts = (real_t *)calloc(info->npts, sizeof(real_t));

        real_t *north_axis = (real_t *)calloc(nnorth, sizeof(real_t));
        real_t *east_axis = (real_t *)calloc(neast, sizeof(real_t));
        // 网格文件按两个坐标轴保存，读取后展开为 point 顺序
        NC_CHECK(nc_inq_varid(ncid, "north", &north_varid));
        NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, north_varid, north_axis));
        NC_CHECK(nc_inq_varid(ncid, "east", &east_varid));
        NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, east_varid, east_axis));
        for(size_t inorth = 0; inorth < nnorth; ++inorth){
            for(size_t ieast = 0; ieast < neast; ++ieast){
                size_t ipt = ieast + inorth * neast;
                info->norths[ipt] = north_axis[inorth];
                info->easts[ipt] = east_axis[ieast];
            }
        }
        GRT_SAFE_FREE_PTR(north_axis);
        GRT_SAFE_FREE_PTR(east_axis);
        return;
    }

    if(info->layout == GRT_RECV_NC_LAYOUT_POINTS){
        int north_varid, east_varid;
        // 一维接收点文件的坐标已经展平，可以直接读取到输出数组
        NC_CHECK(nc_inq_dimid(ncid, "point", &info->dimids[0]));
        NC_CHECK(nc_inq_dimlen(ncid, info->dimids[0], &info->npts));
        info->norths = (real_t *)calloc(info->npts, sizeof(real_t));
        info->easts = (real_t *)calloc(info->npts, sizeof(real_t));
        NC_CHECK(nc_inq_varid(ncid, "north", &north_varid));
        NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, north_varid, info->norths));
        NC_CHECK(nc_inq_varid(ncid, "east", &east_varid));
        NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, east_varid, info->easts));
        return;
    }

}


void grt_recv_nc_info_free(GRT_RECV_NC_INFO *info)
{
    if(info == NULL) return;
    // 坐标数组由该结构体管理，结构体本身由调用方管理
    GRT_SAFE_FREE_PTR(info->norths);
    GRT_SAFE_FREE_PTR(info->easts);
    memset(info, 0, sizeof(*info));
}


bool grt_recv_nc_is_points(int ncid)
{
    return grt_recv_nc_get_layout(ncid) == GRT_RECV_NC_LAYOUT_POINTS;
}
