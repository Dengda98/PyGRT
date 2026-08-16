/**
 * @file   grt_static_strain.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-08
 * 
 *    根据预先合成的静态位移空间导数，组合成静态应变张量
 * 
 */

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    int dummy;
} GRT_MODULE_CTRL;

/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt static strain] %s\n\n", GRT_VERSION);printf(
"    Combine spatial derivatives of static displacements\n"
"    into strain tensor, and write into the same nc file.\n"
"    Input must be a static syn NetCDF computed with -e.\n"
"    Both grid (-X/-Y) and points (-Q) layouts are supported.\n"
"    Six components are written as strain_ZZ, strain_ZR, ...\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt static strain <ingrid> [-h]\n"
"\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt static syn -Gstgrn.nc -Su1e16 -e -Ostsyn.nc\n"
"    grt static strain stsyn.nc\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    (void)Ctrl;
    int opt;
    while ((opt = getopt(argc, argv, ":h")) != -1) {
        switch (opt) {
            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 暂不支持设置其它参数
}

/** 由静态位移偏导合成应变张量 */
static void compute_strain(
    size_t npts, const real_t *norths, const real_t *easts,
    real_t *const u[GRT_CHANNEL_NUM],
    real_t *const upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    real_t *const res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM], bool rot2ZNE)
{
    const char *chs = rot2ZNE ? GRT_ZNE_CODES : GRT_ZRT_CODES;

    for(size_t ir=0; ir<npts; ++ir){
        real_t dist = hypot(norths[ir], easts[ir]);
        // 联络项（1e-5: km→cm）：r≠0 用 u/r；r=0 改用 ∂_r u，与 syn 中 (1/r)∂_θ 有限部分配套
        real_t ur_over_r = GRT_IS_ZERO(dist) ? upar[1][1][ir] : (u[1][ir] / dist * 1e-5);
        real_t ut_over_r = GRT_IS_ZERO(dist) ? upar[1][2][ir] : (u[2][ir] / dist * 1e-5);

        for(int c=0; c<GRT_CHANNEL_NUM; ++c){
            for(int c2=c; c2<GRT_CHANNEL_NUM; ++c2){
                real_t val = 0.5 * (upar[c2][c][ir] + upar[c][c2][ir]);
                if(chs[c]=='R' && chs[c2]=='T'){
                    val -= 0.5 * ut_over_r;
                }
                else if(chs[c]=='T' && chs[c2]=='T'){
                    val += ur_over_r;
                }
                res[c2][c][ir] = val;
            }
        }
    }
}


/** 子模块主函数 */
int static_strain_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    getopt_from_command(Ctrl, argc, argv);

    // 第二个参数为 nc 文件路径
    char *s_ingrid = strdup(argv[1]);

    // nc 文件相关变量
    int in_ncid;
    int in_syn_varids[GRT_CHANNEL_NUM];
    int in_syn_upar_varids[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
    int out_varids[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];

    // 打开 nc 文件
    GRTCheckFileExist(s_ingrid);
    NC_CHECK(nc_open(s_ingrid, NC_WRITE, &in_ncid));

    // 输出分量格式，即是否需要旋转到ZNE
    bool rot2ZNE = false;

    // 读入数据是否旋转到ZNE
    {
        int rot2ZNE_int;
        NC_CHECK(nc_get_att_int(in_ncid, NC_GLOBAL, "rot2ZNE", &rot2ZNE_int));
        rot2ZNE = !(rot2ZNE_int == 0);
    }

    // 三分量
    const char *chs = (rot2ZNE)? GRT_ZNE_CODES : GRT_ZRT_CODES;

    // 读入的数据是否有位移偏导
    int calc_upar;
    NC_CHECK(nc_get_att_int(in_ncid, NC_GLOBAL, "calc_upar", &calc_upar));
    if(calc_upar == 0){
        GRTRaiseError("Input grid didn't have displacement derivatives.");
    }

    // 识别 grid / points 布局，展开为长度 npts 的坐标
    bool is_points = grt_recv_nc_is_points(in_ncid);
    size_t npts;
    real_t *norths_flat = NULL, *easts_flat = NULL;
    int out_ndims;
    int out_dimids[2];

    if(is_points){
        int point_dimid, north_varid, east_varid;
        NC_CHECK(nc_inq_dimid(in_ncid, "point", &point_dimid));
        NC_CHECK(nc_inq_dimlen(in_ncid, point_dimid, &npts));
        norths_flat = (real_t *)calloc(npts, sizeof(real_t));
        easts_flat  = (real_t *)calloc(npts, sizeof(real_t));
        NC_CHECK(nc_inq_varid(in_ncid, "north", &north_varid));
        NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, north_varid, norths_flat));
        NC_CHECK(nc_inq_varid(in_ncid, "east", &east_varid));
        NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, east_varid, easts_flat));
        out_ndims = 1;
        out_dimids[0] = point_dimid;
    } else {
        int north_dimid, east_dimid, north_varid, east_varid;
        size_t nnorth, neast;
        NC_CHECK(nc_inq_dimid(in_ncid, "north", &north_dimid));
        NC_CHECK(nc_inq_dimlen(in_ncid, north_dimid, &nnorth));
        NC_CHECK(nc_inq_dimid(in_ncid, "east", &east_dimid));
        NC_CHECK(nc_inq_dimlen(in_ncid, east_dimid, &neast));

        real_t *norths = (real_t *)calloc(nnorth, sizeof(real_t));
        real_t *easts  = (real_t *)calloc(neast, sizeof(real_t));
        NC_CHECK(nc_inq_varid(in_ncid, "north", &north_varid));
        NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, north_varid, norths));
        NC_CHECK(nc_inq_varid(in_ncid, "east", &east_varid));
        NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, east_varid, easts));

        npts = nnorth * neast;
        norths_flat = (real_t *)calloc(npts, sizeof(real_t));
        easts_flat  = (real_t *)calloc(npts, sizeof(real_t));
        for(size_t inorth = 0; inorth < nnorth; ++inorth){
            for(size_t ieast = 0; ieast < neast; ++ieast){
                size_t ipt = ieast + inorth * neast;
                norths_flat[ipt] = norths[inorth];
                easts_flat[ipt]  = easts[ieast];
            }
        }
        GRT_SAFE_FREE_PTR(norths);
        GRT_SAFE_FREE_PTR(easts);

        out_ndims = 2;
        out_dimids[0] = north_dimid;
        out_dimids[1] = east_dimid;
    }

    // 读入合成位移偏导 varid
    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        char *s_title = NULL;
        GRT_SAFE_ASPRINTF(&s_title, "%c", toupper(chs[c]));
        NC_CHECK(nc_inq_varid(in_ncid, s_title, &in_syn_varids[c]));

        for(int c2=0; c2<GRT_CHANNEL_NUM; ++c2){
            GRT_SAFE_ASPRINTF(&s_title, "%c%c", tolower(chs[c2]), toupper(chs[c]));
            NC_CHECK(nc_inq_varid(in_ncid, s_title, &in_syn_upar_varids[c2][c]));
        }
        GRT_SAFE_FREE_PTR(s_title);
    }

    // 重新进入定义模式
    NC_CHECK(nc_redef(in_ncid));

    // 定义合成结果 varid
    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        char *s_title = NULL;
        for(int c2=c; c2<GRT_CHANNEL_NUM; ++c2){
            // 这里命名顺序要注意，例如 ZR -> 0.5*(u_{z,r} + u_{r,z})
            GRT_SAFE_ASPRINTF(&s_title, "strain_%c%c", toupper(chs[c]), toupper(chs[c2]));
            NC_CHECK(nc_def_var(in_ncid, s_title, NC_REAL, out_ndims, out_dimids, &out_varids[c2][c]));
        }
        GRT_SAFE_FREE_PTR(s_title);
    }

    // 结束定义模式
    NC_CHECK(nc_enddef(in_ncid));

    // 先读入内存
    real_t *u[GRT_CHANNEL_NUM];
    real_t *upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
    // 计算结果
    real_t *res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        u[c] = (real_t *)calloc(npts, sizeof(real_t));
        NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, in_syn_varids[c], u[c]));
        for(int c2=0; c2<GRT_CHANNEL_NUM; ++c2){
            res[c2][c] = (real_t *)calloc(npts, sizeof(real_t));
            upar[c2][c] = (real_t *)calloc(npts, sizeof(real_t));
            NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, in_syn_upar_varids[c2][c], upar[c2][c]));
        }
    }

    compute_strain(npts, norths_flat, easts_flat, u, upar, res, rot2ZNE);

    // 写入 nc 文件
    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        for(int c2=c; c2<GRT_CHANNEL_NUM; ++c2){
            NC_CHECK(NC_FUNC_REAL(nc_put_var) (in_ncid, out_varids[c2][c], res[c2][c]));
        }
    }

    // 关闭文件
    NC_CHECK(nc_close(in_ncid));

    GRT_SAFE_FREE_PTR(norths_flat);
    GRT_SAFE_FREE_PTR(easts_flat);
    GRT_SAFE_FREE_PTR(s_ingrid);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
