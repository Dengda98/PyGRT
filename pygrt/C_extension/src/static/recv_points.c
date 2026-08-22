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
#include <math.h>

#include "grt/static/recv_points.h"
#include "grt/common/checkerror.h"
#include "grt/common/util.h"
#include "grt/common/mynetcdf.h"

/** 解析一行接收点数据，并返回有效数值列数 */
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
    GRT_SAFE_FREE_PTR(pts->norths);
    GRT_SAFE_FREE_PTR(pts->easts);
    GRT_SAFE_FREE_PTR(pts->depths);
    GRT_SAFE_FREE_PTR(pts->strikes);
    GRT_SAFE_FREE_PTR(pts->dips);
    GRT_SAFE_FREE_PTR(pts->rakes);
    GRT_SAFE_FREE_PTR(pts);
}


bool grt_recv_nc_is_points(int ncid)
{
    size_t len = 0;
    int status = nc_inq_attlen(ncid, NC_GLOBAL, "layout", &len);
    if(status == NC_NOERR && len > 0){
        char *layout = (char *)calloc(len + 1, 1);
        if(nc_get_att_text(ncid, NC_GLOBAL, "layout", layout) == NC_NOERR){
            bool is_pts = (strcmp(layout, GRT_RECV_LAYOUT_POINTS) == 0);
            GRT_SAFE_FREE_PTR(layout);
            return is_pts;
        }
        GRT_SAFE_FREE_PTR(layout);
    }

    int dimid;
    if(nc_inq_dimid(ncid, "point", &dimid) == NC_NOERR){
        return true;
    }
    return false;
}
