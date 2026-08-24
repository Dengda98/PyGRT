/**
 * @file   grt_static_coulomb.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08-24
 *
 *    Compute Coulomb stress change from sigma_n and tau_s
 *
 */

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    // 输入静态合成文件
    struct {
        bool active;
        char *s_ingrid;
    } G;
    // 等效摩擦系数
    struct {
        bool active;
        real_t friction;
    } F;
} GRT_MODULE_CTRL;


/**
 * 释放参数控制结构体及其动态成员
 *
 * @param[in,out] Ctrl 参数控制结构体
 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    if(Ctrl == NULL) return;
    GRT_SAFE_FREE_PTR(Ctrl->G.s_ingrid);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印模块使用说明 */
static void print_help(void)
{
printf("\n"
"[grt static coulomb] %s\n\n", GRT_VERSION);printf(
"    Compute Coulomb stress change from sigma_n and tau_s.\n"
"    The result is written back to the input nc file as coulomb.\n"
"    The input file must already contain sigma_n and tau_s,\n"
"    which are normally produced by `static sproj`.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt static coulomb -G<ingrid> -F<friction> [-h]\n"
"\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<ingrid>    Input static synthesis file containing sigma_n\n"
"                  and tau_s. The result is written into this file.\n"
"\n"
"    -F<friction>  Nonnegative dimensionless effective friction coefficient.\n"
"\n"
"    -h            Display this help message.\n"
"\n\n"
);
}


/**
 * 从命令行中读取模块选项
 *
 * @param[out] Ctrl 参数控制结构体
 * @param[in]  argc 命令行参数数量
 * @param[in]  argv 命令行参数数组
 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    int opt;

    // 解析输入文件和等效摩擦系数
    while((opt = getopt(argc, argv, ":G:F:h")) != -1){
        switch(opt){
            case 'G':
                Ctrl->G.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->G.s_ingrid);
                Ctrl->G.s_ingrid = strdup(optarg);
                break;

            case 'F': {
                char extra;
                int nscan = sscanf(optarg, "%lf%c", &Ctrl->F.friction, &extra);
                if((nscan != 1) || !isfinite(Ctrl->F.friction) ||
                    (Ctrl->F.friction < 0.0)){
                    GRTBadOptionError(F, "Friction must be a finite nonnegative number.");
                }
                Ctrl->F.active = true;
                break;
            }

            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 输入文件和摩擦系数都是必选项，模块不接受位置参数
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, G);
    GRTCheckOptionActive(Ctrl, F);
    if(optind < argc){
        GRTRaiseError("Unexpected positional argument \"%s\". Use -G<ingrid>.\n", argv[optind]);
    }
}


/**
 * 检查 NetCDF 变量的维度是否与接收布局一致
 *
 * @param[in] ncid            NetCDF 文件 ID
 * @param[in] varid           待检查的变量 ID
 * @param[in] name            变量名称
 * @param[in] expected_ndims  期望的维度数量
 * @param[in] expected_dimids 期望的维度 ID 数组
 */
static void check_var_dimensions(
    int ncid, int varid, const char *name, int expected_ndims,
    const int *expected_dimids)
{
    int ndims;
    int dimids[NC_MAX_DIMS];
    NC_CHECK(nc_inq_varndims(ncid, varid, &ndims));
    if(ndims != expected_ndims){
        GRTRaiseError(
            "Variable \"%s\" has %d dimensions, expected %d.",
            name, ndims, expected_ndims);
    }
    NC_CHECK(nc_inq_vardimid(ncid, varid, dimids));
    for(int i = 0; i < expected_ndims; ++i){
        if(dimids[i] != expected_dimids[i]){
            GRTRaiseError("Variable \"%s\" has an unexpected dimension.", name);
        }
    }
}


/**
 * 获取必须存在的投影应力变量
 *
 * @param[in] ncid       NetCDF 文件 ID
 * @param[in] name       变量名称
 * @param[in] ndims      期望的维度数量
 * @param[in] dimids     期望的维度 ID 数组
 * @return 变量 ID
 */
static int get_projection_var(
    int ncid, const char *name, int ndims, const int *dimids)
{
    int varid;
    int status = nc_inq_varid(ncid, name, &varid);
    if(status == NC_ENOTVAR){
        GRTRaiseError(
            "Input file does not contain variable \"%s\". Run static sproj first.",
            name);
    }
    NC_CHECK(status);
    check_var_dimensions(ncid, varid, name, ndims, dimids);
    return varid;
}


/**
 * 将 Coulomb 应力变化追加或覆盖写入原始 NetCDF 文件
 *
 * @param[in] ncid       输入 NetCDF 文件 ID
 * @param[in] ndims      输出变量维度数量
 * @param[in] dimids     输出变量维度 ID 数组
 * @param[in] coulomb    Coulomb 应力变化
 */
static void write_coulomb_variable(
    int ncid, int ndims, const int *dimids, const real_t *coulomb)
{
    int varid;
    int status = nc_inq_varid(ncid, "coulomb", &varid);
    bool has_coulomb = status == NC_NOERR;
    if((status != NC_NOERR) && (status != NC_ENOTVAR)) NC_CHECK(status);
    if(has_coulomb){
        check_var_dimensions(ncid, varid, "coulomb", ndims, dimids);
        GRTRaiseWarning("Variable \"coulomb\" already exists and will be overwritten.");
    } else {
        // 只定义缺失的结果变量，保留输入文件已有数据
        NC_CHECK(nc_redef(ncid));
        NC_CHECK(nc_def_var(ncid, "coulomb", NC_REAL, ndims, dimids, &varid));
        NC_CHECK(nc_enddef(ncid));
    }

    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, varid, coulomb));
}


/**
 * 执行静态 Coulomb 应力变化计算
 *
 * @param[in] argc 命令行参数数量
 * @param[in] argv 命令行参数数组
 * @return 模块执行状态
 */
int static_coulomb_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);

    GRTCheckFileExist(Ctrl->G.s_ingrid);

    // 以写模式打开输入文件，结果直接写回该文件
    int ncid;
    NC_CHECK(nc_open(Ctrl->G.s_ingrid, NC_WRITE, &ncid));

    GRT_RECV_NC_INFO recv_info;
    grt_recv_nc_info_load(ncid, &recv_info);
    int ndims = (recv_info.layout == GRT_RECV_NC_LAYOUT_POINTS) ? 1 : 2;
    size_t npts = recv_info.npts;

    int sigma_n_varid = get_projection_var(ncid, "sigma_n", ndims, recv_info.dimids);
    int tau_s_varid = get_projection_var(ncid, "tau_s", ndims, recv_info.dimids);
    real_t *sigma_n = (real_t *)calloc(npts, sizeof(real_t));
    real_t *tau_s = (real_t *)calloc(npts, sizeof(real_t));
    real_t *coulomb = (real_t *)calloc(npts, sizeof(real_t));

    // 读取两个投影应力分量，并按点应用 Coulomb 应力公式
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, sigma_n_varid, sigma_n));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, tau_s_varid, tau_s));
    for(size_t i = 0; i < npts; ++i){
        coulomb[i] = tau_s[i] + Ctrl->F.friction * sigma_n[i];
    }

    write_coulomb_variable(ncid, ndims, recv_info.dimids, coulomb);
    NC_CHECK(nc_close(ncid));

    grt_recv_nc_info_free(&recv_info);
    GRT_SAFE_FREE_PTR(sigma_n);
    GRT_SAFE_FREE_PTR(tau_s);
    GRT_SAFE_FREE_PTR(coulomb);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
