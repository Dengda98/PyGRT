/**
 * @file   rcv_point.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * 动态解和静态解共用的接收点列表
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "grt/common/rcv_point.h"
#include "grt/common/checkerror.h"
#include "grt/common/util.h"


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

    // 逐个读取数值，空白分隔，遇到行尾或注释终止，只接受三列或六列
    while(true){
        while(isspace((unsigned char)*cursor)) {
            ++cursor;
        }
        if(*cursor == '\0' || *cursor == GRT_COMMENT_HEAD) {
            break;
        }
        if(*nvalues >= 6) {
            return false;
        }

        errno = 0;
        char *end = NULL;
        real_t value = strtod(cursor, &end);
        if(end == cursor || errno == ERANGE) {
            return false;
        }
        values[*nvalues] = value;
        (*nvalues)++;
        cursor = end;
    }

    return *nvalues == 3 || *nvalues == 6;
}

RCV_POINT *grt_rcv_points_from_polar(real_t dist, real_t azimuth, real_t depth)
{
    // 极坐标只在输入时转换，后续源台计算使用单点坐标
    RCV_POINT *rcvs = GRT_SAFE_CALLOC(1, sizeof(*rcvs));
    rcvs[0].north = dist * cos(azimuth * DEG1);
    rcvs[0].east = dist * sin(azimuth * DEG1);
    rcvs[0].depth = depth;
    return rcvs;
}

RCV_POINT *grt_rcv_points_from_grid(
    size_t nnorth, const real_t *norths,
    size_t neast,  const real_t *easts,
    real_t depth)
{
    if(nnorth == 0 || neast == 0 || norths == NULL || easts == NULL){
        GRTRaiseError("empty receiver grid.");
    }
    // 网格按 east 方向最快变化的顺序展开为接收点列表
    RCV_POINT *rcvs = GRT_SAFE_CALLOC(nnorth * neast, sizeof(*rcvs));
    for(size_t inorth = 0; inorth < nnorth; ++inorth){
        for(size_t ieast = 0; ieast < neast; ++ieast){
            size_t ipt = ieast + inorth * neast;
            rcvs[ipt].north = norths[inorth];
            rcvs[ipt].east  = easts[ieast];
            rcvs[ipt].depth = depth;
        }
    }
    return rcvs;
}


RCV_POINT *grt_rcv_points_from_file(const char *path, size_t *count)
{
    if(!count) {
        GRTRaiseError("Receiver point count is NULL.");
    }
    GRTCheckFileExist(path);

    FILE *fp = fopen(path, "r");
    if(fp == NULL){
        GRTRaiseError("Failed to open receiver points file \"%s\".", path);
    }

    // 每行只解析一次，统一检查列数后直接追加一个接收点
    RCV_POINT *rcvs = NULL;
    size_t npts = 0;
    size_t capacity = 0;
    size_t ncolumns = 0;
    char *line = NULL;
    size_t nlen = 0;
    size_t lineno = 0;
    while(grt_getline(&line, &nlen, fp) != -1){
        lineno++;
        grt_trim_whitespace(line);
        if(grt_is_comment_or_empty_line(line)) {
            continue;
        }

        real_t values[6];
        size_t nvalues;
        if(!parse_receiver_point_line(line, values, &nvalues)){
            GRTRaiseError(
                "Invalid receiver point at line %zu in \"%s\" "
                "(expect exactly 3 or 6 numeric columns: "
                "north east depth [strike dip rake]).",
                lineno, path);
        }

        // 整个接收文件必须使用相同列数，避免部分点缺失接收机制
        if(ncolumns != 0 && nvalues != ncolumns){
            GRTRaiseError(
                "Inconsistent receiver point column count at line %zu in \"%s\" "
                "(all data lines must have either 3 or 6 columns).",
                lineno, path);
        }
        ncolumns = nvalues;

        // 按需扩大点集容量，避免每读一行都重新分配数组
        if(npts == capacity) {
            capacity = capacity ? 2 * capacity : 64;
            rcvs = GRT_SAFE_REALLOC(rcvs, capacity * sizeof(*rcvs));
        }
        size_t ipt = npts++;
        rcvs[ipt] = (RCV_POINT){0};
        if(values[2] < 0.0){
            GRTRaiseError("Negative receiver depth at line %zu in \"%s\".", lineno, path);
        }
        rcvs[ipt].north = values[0];
        rcvs[ipt].east  = values[1];
        rcvs[ipt].depth = values[2];

        // 六列格式的后三列作为接收机制，三列格式保持机制未定义
        if(ncolumns == 6){
            if(values[3] < 0 || values[3] > 360 || values[4] < 0 || values[4] > 90 || fabs(values[5]) > 180) {
                GRTRaiseError("Invalid receiver mechanism at line %zu in \"%s\" (strike/dip/rake must be in [0,360]/[0,90]/[-180,180]).",
                              lineno, path);
            }
            rcvs[ipt].has_mechanism = true;
            rcvs[ipt].strike = values[3];
            rcvs[ipt].dip    = values[4];
            rcvs[ipt].rake   = values[5];
        }
    }

    if(!npts) {
        GRTRaiseError("No receiver points found in \"%s\".", path);
    }
    *count = npts;
    GRT_SAFE_FREE_PTR(line);
    fclose(fp);
    return rcvs;
}


RCV_POINT *grt_rcv_points_from_faults(size_t nfault, const FINITE_FAULT *faults, size_t *npts)
{
    // 先统计全部剖分点，再按断层及断层内顺序一次性展开
    *npts = 0;
    for(size_t i = 0; i < nfault; ++i) {
        *npts += faults[i].nW * faults[i].nL;
    }
    RCV_POINT *points = GRT_SAFE_MALLOC(*npts * sizeof(*points));
    size_t ipt = 0;
    for(size_t i = 0; i < nfault; ++i) {
        const FINITE_FAULT *fault = &faults[i];
        for(size_t isub = 0; isub < fault->nW * fault->nL; ++isub) {
            points[ipt++] = (RCV_POINT){.north = fault->north[isub], .east = fault->east[isub], .depth = fault->depth[isub],
                                       .fault = fault, .isub = isub};
        }
    }
    return points;
}
