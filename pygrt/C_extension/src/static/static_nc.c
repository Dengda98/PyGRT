/**
 * @file   static_nc.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * 静态解 NetCDF 接收布局读入和结果输出
 *
 */

#include <ctype.h>
#include <stdint.h>
#include <string.h>

#include "grt/static/static_nc.h"
#include "grt/common/model.h"
#include "grt/common/mynetcdf.h"

/** 接收变量的 NetCDF ID，仅在一次写出过程中使用 */
typedef struct {
    int ndims;           ///< 结果变量维度数
    int dimids[2];       ///< 结果变量维度 ID，grid 按 north/east 排列
    int coordinate[3];   ///< north/east/depth 变量 ID，grid 不使用 depth
    int medium[3];       ///< rcv_va/rcv_vb/rcv_rho 变量 ID，用于 points 和 faults
    int geometry[3];     ///< strike/dip/rake 变量 ID，有断层时按 nfault，否则按 point
    int fault_shape[3];  ///< offset/stksize/dipsize 变量 ID，仅用于有限接收断层
} RCV_NC_VARS;

static const char *const layout_names[] = {GRT_RCV_LAYOUT_GRID, GRT_RCV_LAYOUT_POINTS, GRT_RCV_LAYOUT_FAULTS};
static const char *const coordinate_names[] = {"north", "east", "depth"};
static const char *const medium_names[] = {"rcv_va", "rcv_vb", "rcv_rho"};
static const char *const geometry_names[] = {"strike", "dip", "rake"};
static const char *const fault_shape_names[] = {"offset", "stksize", "dipsize"};

/**
 * 定义坐标、介质及可选机制变量，grid 的接收深度和介质保存为属性
 * @param[in]   ncid       NetCDF 文件 ID
 * @param[in]   receivers  接收布局
 * @param[in]   nlayer     原始模型层数
 * @param[in]   modarr     原始模型矩阵
 * @param[out]  vars       接收变量 ID
 */
static void define_receiver_vars(int ncid, const RCV_NC_INFO *receivers, size_t nlayer,
                                 const real_t (*modarr)[GRT_MODARR_NCOL], RCV_NC_VARS *vars)
{
    if(receivers->layout == GRT_RCV_NC_LAYOUT_GRID) {
        real_t depth = receivers->rcvs[0].depth, medium[3];
        grt_modarr_medium_at_depth(nlayer, modarr, depth, &medium[0], &medium[1], &medium[2]);
        NC_CHECK(NC_FUNC_REAL(nc_put_att)(ncid, NC_GLOBAL, "deprcv", NC_REAL, 1, &depth));
        for(int c = 0; c < 3; ++c) {
            NC_CHECK(NC_FUNC_REAL(nc_put_att)(ncid, NC_GLOBAL, medium_names[c], NC_REAL, 1, &medium[c]));
        }
        vars->ndims = 2;
        NC_CHECK(nc_def_dim(ncid, "north", receivers->nnorth, &vars->dimids[0]));
        NC_CHECK(nc_def_dim(ncid, "east", receivers->neast, &vars->dimids[1]));
        for(int c = 0; c < 2; ++c) {
            NC_CHECK(nc_def_var(ncid, coordinate_names[c], NC_REAL, 1, &vars->dimids[c], &vars->coordinate[c]));
        }
        return;
    }

    vars->ndims = 1;
    NC_CHECK(nc_def_dim(ncid, "point", receivers->npts, &vars->dimids[0]));
    for(int c = 0; c < 3; ++c) {
        NC_CHECK(nc_def_var(ncid, coordinate_names[c], NC_REAL, 1, vars->dimids, &vars->coordinate[c]));
    }
    for(int c = 0; c < 3; ++c) {
        NC_CHECK(nc_def_var(ncid, medium_names[c], NC_REAL, 1, vars->dimids, &vars->medium[c]));
    }

    // 普通点机制使用 point 维，有限断层机制使用 nfault 维
    if(receivers->layout == GRT_RCV_NC_LAYOUT_FAULTS || receivers->rcvs[0].has_mechanism) {
        int geometry_dimid = vars->dimids[0];
        if(receivers->layout == GRT_RCV_NC_LAYOUT_FAULTS) {
            NC_CHECK(nc_def_dim(ncid, "nfault", receivers->nfault, &geometry_dimid));
        }
        for(int c = 0; c < 3; ++c) {
            NC_CHECK(nc_def_var(ncid, geometry_names[c], NC_REAL, 1, &geometry_dimid, &vars->geometry[c]));
        }
        if(receivers->layout == GRT_RCV_NC_LAYOUT_FAULTS) {
            for(int c = 0; c < 3; ++c) {
                NC_CHECK(nc_def_var(ncid, fault_shape_names[c], NC_INT, 1, &geometry_dimid, &vars->fault_shape[c]));
            }
        }
    }
}

/**
 * 定义位移及可选位移偏导变量
 * @param[in]   ncid       NetCDF 文件 ID
 * @param[in]   receivers  接收变量 ID 及维度
 * @param[in]   channels   分量编码
 * @param[in]   calc_upar  是否记录位移偏导
 * @param[out]  vars       位移变量 ID
 * @param[out]  dvars      位移偏导变量 ID
 */
static void define_channel_vars(int ncid, const RCV_NC_VARS *receivers, const char *channels, bool calc_upar,
                                int vars[GRT_CHANNEL_NUM], int dvars[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    for(int c = 0; c < GRT_CHANNEL_NUM; ++c) {
        char name[3] = {channels[c], '\0', '\0'};
        NC_CHECK(nc_def_var(ncid, name, NC_REAL, receivers->ndims, receivers->dimids, &vars[c]));
        if(calc_upar) {
            for(int d = 0; d < GRT_CHANNEL_NUM; ++d) {
                name[0] = tolower(channels[d]);
                name[1] = channels[c];
                NC_CHECK(nc_def_var(ncid, name, NC_REAL, receivers->ndims, receivers->dimids, &dvars[d][c]));
            }
        }
    }
}

/**
 * 写入坐标及逐点介质，复用缓冲提取各变量的点序列
 * @param[in]   ncid       NetCDF 文件 ID
 * @param[in]   receivers  接收布局
 * @param[in]   nlayer     原始模型层数
 * @param[in]   modarr     原始模型矩阵
 * @param[in]   vars       接收变量 ID
 * @param[out]  buffer     至少能容纳 npts 个实数的工作缓冲
 */
static void write_receiver_coordinates(int ncid, const RCV_NC_INFO *receivers, size_t nlayer,
                                      const real_t (*modarr)[GRT_MODARR_NCOL], const RCV_NC_VARS *vars, real_t *buffer)
{
    const RCV_POINT *rcvs = receivers->rcvs;
    if(receivers->layout == GRT_RCV_NC_LAYOUT_GRID) {
        for(size_t i = 0; i < receivers->nnorth; ++i) {
            buffer[i] = rcvs[i * receivers->neast].north;
        }
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, vars->coordinate[0], buffer));
        for(size_t i = 0; i < receivers->neast; ++i) {
            buffer[i] = rcvs[i].east;
        }
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, vars->coordinate[1], buffer));
        return;
    }

    for(int c = 0; c < 3; ++c) {
        for(size_t i = 0; i < receivers->npts; ++i) {
            buffer[i] = c == 0 ? rcvs[i].north : c == 1 ? rcvs[i].east : rcvs[i].depth;
        }
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, vars->coordinate[c], buffer));
    }

    // 各点只查询一次介质，再按变量提取；临时介质数组在本阶段结束时释放
    real_t (*medium)[3] = GRT_SAFE_CALLOC(receivers->npts, sizeof(*medium));
    for(size_t i = 0; i < receivers->npts; ++i) {
        grt_modarr_medium_at_depth(nlayer, modarr, rcvs[i].depth, &medium[i][0], &medium[i][1], &medium[i][2]);
    }
    for(int c = 0; c < 3; ++c) {
        for(size_t i = 0; i < receivers->npts; ++i) {
            buffer[i] = medium[i][c];
        }
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, vars->medium[c], buffer));
    }
    GRT_SAFE_FREE_PTR(medium);
}

/**
 * 写入逐点机制或有限断层机制及点范围
 * @param[in]   ncid       NetCDF 文件 ID
 * @param[in]   receivers  接收布局
 * @param[in]   vars       接收变量 ID
 * @param[out]  buffer     至少能容纳 npts 个实数的工作缓冲
 */
static void write_receiver_geometry(int ncid, const RCV_NC_INFO *receivers, const RCV_NC_VARS *vars, real_t *buffer)
{
    if(receivers->layout == GRT_RCV_NC_LAYOUT_GRID ||
       (receivers->layout == GRT_RCV_NC_LAYOUT_POINTS && !receivers->rcvs[0].has_mechanism)) {
        return;
    }
    size_t count = receivers->layout == GRT_RCV_NC_LAYOUT_FAULTS ? receivers->nfault : receivers->npts;
    for(int c = 0; c < 3; ++c) {
        for(size_t i = 0; i < count; ++i) {
            if(receivers->layout == GRT_RCV_NC_LAYOUT_FAULTS) {
                const FINITE_FAULT *fault = &receivers->faults[i];
                buffer[i] = c == 0 ? fault->strike : c == 1 ? fault->dip : fault->rake;
            } else {
                const RCV_POINT *point = &receivers->rcvs[i];
                buffer[i] = c == 0 ? point->strike : c == 1 ? point->dip : point->rake;
            }
        }
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, vars->geometry[c], buffer));
    }
    if(receivers->layout != GRT_RCV_NC_LAYOUT_FAULTS) {
        return;
    }

    // offset 为各断层点范围的排他性结束索引，最后一个值等于 point
    int *shape = GRT_SAFE_CALLOC(count, sizeof(*shape));
    for(int c = 0; c < 3; ++c) {
        size_t offset = 0;
        for(size_t i = 0; i < count; ++i) {
            const FINITE_FAULT *fault = &receivers->faults[i];
            offset += fault->nL * fault->nW;
            shape[i] = c == 0 ? offset : c == 1 ? fault->nL : fault->nW;
        }
        NC_CHECK(nc_put_var_int(ncid, vars->fault_shape[c], shape));
    }
    GRT_SAFE_FREE_PTR(shape);
}

/**
 * 按分量写入位移及可选偏导
 * @param[in]   ncid       NetCDF 文件 ID
 * @param[in]   npts       接收点数
 * @param[in]   calc_upar  是否记录位移偏导
 * @param[in]   syn        位移数组
 * @param[in]   syn_upar   位移偏导数组
 * @param[in]   vars       位移变量 ID
 * @param[in]   dvars      位移偏导变量 ID
 * @param[out]  buffer     至少能容纳 npts 个实数的工作缓冲
 */
static void write_fields(int ncid, size_t npts, bool calc_upar, const real_t (*syn)[GRT_CHANNEL_NUM],
                         const real_t (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM], const int vars[GRT_CHANNEL_NUM],
                         const int dvars[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM], real_t *buffer)
{
    for(int c = 0; c < GRT_CHANNEL_NUM; ++c) {
        for(size_t i = 0; i < npts; ++i) {
            buffer[i] = syn[i][c];
        }
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, vars[c], buffer));
        if(calc_upar) {
            for(int d = 0; d < GRT_CHANNEL_NUM; ++d) {
                for(size_t i = 0; i < npts; ++i) {
                    buffer[i] = syn_upar[i][d][c];
                }
                NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, dvars[d][c], buffer));
            }
        }
    }
}

void grt_static_nc_write(int ncid, const RCV_NC_INFO *receivers, size_t nlayer, const real_t (*modarr)[GRT_MODARR_NCOL],
                         bool rot2ZNE, bool calc_upar, const real_t (*syn)[GRT_CHANNEL_NUM],
                         const real_t (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    if(nlayer == 0 || modarr == NULL || syn == NULL || (calc_upar && syn_upar == NULL)) {
        GRTRaiseError("Static NetCDF model or displacement data are missing.");
    }

    // 公共层只描述接收布局和结果，震源及模块专属属性由调用方写入
    const char *layout = layout_names[receivers->layout];
    const char *channels = rot2ZNE ? GRT_ZNE_CODES : GRT_ZRT_CODES;
    int upar_attr = calc_upar, zne_attr = rot2ZNE;
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "layout", strlen(layout), layout));
    NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "calc_upar", NC_INT, 1, &upar_attr));
    NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "rot2ZNE", NC_INT, 1, &zne_attr));

    // 先完成所有变量定义，再统一结束定义模式并逐阶段写入数据
    RCV_NC_VARS rcv_vars = {0};
    int vars[GRT_CHANNEL_NUM], dvars[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
    define_receiver_vars(ncid, receivers, nlayer, modarr, &rcv_vars);
    define_channel_vars(ncid, &rcv_vars, channels, calc_upar, vars, dvars);
    NC_CHECK(nc_enddef(ncid));

    real_t *buffer = GRT_SAFE_CALLOC(receivers->npts, sizeof(*buffer));
    write_receiver_coordinates(ncid, receivers, nlayer, modarr, &rcv_vars, buffer);
    write_receiver_geometry(ncid, receivers, &rcv_vars, buffer);
    write_fields(ncid, receivers->npts, calc_upar, syn, syn_upar, vars, dvars, buffer);
    GRT_SAFE_FREE_PTR(buffer);
}

GRT_RCV_NC_LAYOUT grt_rcv_nc_get_layout(int ncid)
{
    size_t len = 0;
    int status = nc_inq_attlen(ncid, NC_GLOBAL, "layout", &len);
    if(status != NC_NOERR || len == 0) {
        GRTRaiseError("static receiver layout attribute is missing.");
    }
    char *layout = GRT_SAFE_CALLOC(len + 1, 1);
    NC_CHECK(nc_get_att_text(ncid, NC_GLOBAL, "layout", layout));

    for(size_t i = 0; i < sizeof(layout_names) / sizeof(*layout_names); ++i) {
        if(strcmp(layout, layout_names[i]) == 0) {
            GRT_SAFE_FREE_PTR(layout);
            return (GRT_RCV_NC_LAYOUT)i;
        }
    }
    GRTRaiseError("unsupported static receiver layout \"%s\".", layout);
    return GRT_RCV_NC_LAYOUT_GRID;
}

/**
 * 读取一维坐标变量，核对维度后再写入调用方缓冲
 * @param[in]   ncid    NetCDF 文件 ID
 * @param[in]   name    坐标变量名
 * @param[in]   dimid   预期维度 ID
 * @param[out]  values  能容纳该维度长度的坐标缓冲
 */
static void read_coordinate(int ncid, const char *name, int dimid, real_t *values)
{
    int varid, ndims, actual_dimid;
    NC_CHECK(nc_inq_varid(ncid, name, &varid));
    NC_CHECK(nc_inq_varndims(ncid, varid, &ndims));
    if(ndims != 1) {
        GRTRaiseError("Static receiver coordinate %s must have one dimension.", name);
    }
    NC_CHECK(nc_inq_vardimid(ncid, varid, &actual_dimid));
    if(actual_dimid != dimid) {
        GRTRaiseError("Static receiver coordinate %s has an unexpected dimension.", name);
    }
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, values));
}

void grt_rcv_nc_info_load(int ncid, RCV_NC_INFO *info)
{
    if(info == NULL) {
        GRTRaiseError("receiver NetCDF info is NULL.");
    }
    memset(info, 0, sizeof(*info));
    info->layout = grt_rcv_nc_get_layout(ncid);
    info->ndims = info->layout == GRT_RCV_NC_LAYOUT_GRID ? 2 : 1;

    // 从 layout 统一确定结果维度及点数，所有坐标变量必须使用对应的一维坐标维
    size_t lengths[2] = {1, 1};
    for(int d = 0; d < info->ndims; ++d) {
        const char *name = info->layout == GRT_RCV_NC_LAYOUT_GRID ? coordinate_names[d] : "point";
        NC_CHECK(nc_inq_dimid(ncid, name, &info->dimids[d]));
        NC_CHECK(nc_inq_dimlen(ncid, info->dimids[d], &lengths[d]));
        if(lengths[d] == 0) {
            GRTRaiseError("Static receiver dimension %s is empty.", name);
        }
    }
    if(lengths[0] > SIZE_MAX / lengths[1]) {
        GRTRaiseError("Static receiver grid is too large.");
    }
    info->npts = lengths[0] * lengths[1];
    if(info->layout == GRT_RCV_NC_LAYOUT_GRID) {
        info->nnorth = lengths[0];
        info->neast  = lengths[1];
    }

    // 有限断层由 layout 显式标识，其他布局不允许包含 nfault 维度
    int nfault_dimid;
    int status = nc_inq_dimid(ncid, "nfault", &nfault_dimid);
    if(status != NC_NOERR && status != NC_EBADDIM) {
        NC_CHECK(status);
    }
    if(info->layout == GRT_RCV_NC_LAYOUT_FAULTS) {
        if(status != NC_NOERR) {
            GRTRaiseError("Static receiver faults layout must contain an nfault dimension.");
        }
        info->nfault_dimid = nfault_dimid;
        NC_CHECK(nc_inq_dimlen(ncid, nfault_dimid, &info->nfault));
        if(info->nfault == 0 || info->nfault > info->npts) {
            GRTRaiseError("Static receiver fault count must be nonzero and not exceed the point count.");
        }
    } else if(status == NC_NOERR) {
        GRTRaiseError("Static receiver grid and points layouts must not contain an nfault dimension.");
    }
}

void grt_rcv_nc_info_load_coordinates(int ncid, RCV_NC_INFO *info)
{
    if(info == NULL || info->npts == 0) {
        GRTRaiseError("receiver NetCDF info has no points.");
    }
    GRT_SAFE_FREE_PTR(info->norths);
    GRT_SAFE_FREE_PTR(info->easts);
    info->norths = GRT_SAFE_CALLOC(info->npts, sizeof(*info->norths));
    info->easts = GRT_SAFE_CALLOC(info->npts, sizeof(*info->easts));

    if(info->layout != GRT_RCV_NC_LAYOUT_GRID) {
        read_coordinate(ncid, "north", info->dimids[0], info->norths);
        read_coordinate(ncid, "east", info->dimids[0], info->easts);
        return;
    }

    // grid 的坐标轴逐个读取，再按 east 方向最快变化的点序展开
    real_t *axis = GRT_SAFE_CALLOC(GRT_MAX(info->nnorth, info->neast), sizeof(*axis));
    for(int d = 0; d < 2; ++d) {
        read_coordinate(ncid, coordinate_names[d], info->dimids[d], axis);
        real_t *values = d == 0 ? info->norths : info->easts;
        for(size_t i = 0; i < info->npts; ++i) {
            values[i] = axis[d == 0 ? i / info->neast : i % info->neast];
        }
    }
    GRT_SAFE_FREE_PTR(axis);
}

void grt_rcv_nc_info_free(RCV_NC_INFO *info)
{
    if(info == NULL) {
        return;
    }
    GRT_SAFE_FREE_PTR(info->norths);
    GRT_SAFE_FREE_PTR(info->easts);
    memset(info, 0, sizeof(*info));
}
