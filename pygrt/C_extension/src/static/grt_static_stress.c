/**
 * @file   grt_static_stress.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-08
 * 
 *    根据预先合成的静态位移空间导数，组合成静态应力张量
 * 
 */

#include "grt/common/const.h"
#include "grt/common/util.h"
#include "grt/common/mynetcdf.h"

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
"[grt static stress] %s\n\n", GRT_VERSION);printf(
"    Conbine spatial derivatives of static displacements\n"
"    into stress tensor (unit: dyne/cm^2 = 0.1 Pa), .\n"
"    and write into the same nc file.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt static stress <ingrid>\n"
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


/** 子模块主函数 */
int static_stress_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    getopt_from_command(Ctrl, argc, argv);

    // 第二个参数为 nc 文件路径
    char *s_ingrid = strdup(argv[1]);

    // nc 文件相关变量
    int in_ncid;
    int in_x_dimid, in_y_dimid;
    int in_x_varid, in_y_varid;
    const int ndims = 2;
    int in_dimids[ndims];
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

    // 读入接收点的物性参数，计算 rcv_mu 和 rcv_lambda
    real_t rcv_mu, rcv_lam;
    {
        real_t rcv_va=0.0, rcv_vb=0.0, rcv_rho=0.0;
        NC_CHECK(NC_FUNC_REAL(nc_get_att) (in_ncid, NC_GLOBAL, "rcv_va", &rcv_va));
        NC_CHECK(NC_FUNC_REAL(nc_get_att) (in_ncid, NC_GLOBAL, "rcv_vb", &rcv_vb));
        NC_CHECK(NC_FUNC_REAL(nc_get_att) (in_ncid, NC_GLOBAL, "rcv_rho", &rcv_rho));
        rcv_mu = rcv_vb*rcv_vb*rcv_rho*1e10;
        rcv_lam = rcv_va*rcv_va*rcv_rho*1e10 - 2.0*rcv_mu;
    }

    // 读入坐标变量 dimid, varid
    size_t nx, ny;
    NC_CHECK(nc_inq_dimid(in_ncid, "north", &in_x_dimid));
    NC_CHECK(nc_inq_dimlen(in_ncid, in_x_dimid, &nx));
    NC_CHECK(nc_inq_dimid(in_ncid, "east", &in_y_dimid));
    NC_CHECK(nc_inq_dimlen(in_ncid, in_y_dimid, &ny));
    in_dimids[0] = in_x_dimid;
    in_dimids[1] = in_y_dimid;

    // 读取坐标变量
    real_t *xs = (real_t *)calloc(nx, sizeof(real_t));
    real_t *ys = (real_t *)calloc(ny, sizeof(real_t));
    NC_CHECK(nc_inq_varid(in_ncid, "north", &in_x_varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, in_x_varid, xs));
    NC_CHECK(nc_inq_varid(in_ncid, "east", &in_y_varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, in_y_varid, ys));

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
            // 这里命名顺序要注意，例如 ZR -> mu*(u_{z,r} + u_{r,z})
            GRT_SAFE_ASPRINTF(&s_title, "stress_%c%c", toupper(chs[c]), toupper(chs[c2]));
            NC_CHECK(nc_def_var(in_ncid, s_title, NC_REAL, ndims, in_dimids, &out_varids[c2][c]));
        }
        GRT_SAFE_FREE_PTR(s_title);
    }

    // 结束定义模式
    NC_CHECK(nc_enddef(in_ncid));

    // 总震中距数
    size_t nr = nx * ny;

    // 先读入内存，
    real_t *u[GRT_CHANNEL_NUM];
    real_t *upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
    // 计算结果
    real_t *res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        u[c] = (real_t *)calloc(nr, sizeof(real_t));
        NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, in_syn_varids[c], u[c]));
        for(int c2=0; c2<GRT_CHANNEL_NUM; ++c2){
            res[c2][c] = (real_t *)calloc(nr, sizeof(real_t));
            upar[c2][c] = (real_t *)calloc(nr, sizeof(real_t));
            NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, in_syn_upar_varids[c2][c], upar[c2][c]));
        }
    }
    
    // 每个点逐个处理
    for(size_t ix=0; ix < nx; ++ix){
        real_t x = xs[ix];
        for(size_t iy=0; iy < ny; ++iy){
            real_t y = ys[iy];

            size_t ir = iy + ix*ny;

            // 震中距
            real_t dist = GRT_MAX(sqrt(x*x + y*y), GRT_MIN_DISTANCE);

            // 先计算体积应变u_kk = u_11 + u22 + u33 和 lamda的乘积，ZRT分量需包括协变导数
            real_t lam_ukk = upar[0][0][ir] + upar[1][1][ir] + upar[2][2][ir];
            if(!rot2ZNE)  lam_ukk += u[1][ir] / dist*1e-5;
            lam_ukk *= rcv_lam;

            for(int c=0; c<GRT_CHANNEL_NUM; ++c){
                for(int c2=c; c2<GRT_CHANNEL_NUM; ++c2){
                    real_t val = rcv_mu * (upar[c2][c][ir] + upar[c][c2][ir]);

                    // 对角线分量
                    if(c == c2)   val += lam_ukk;
                    
                    // 特殊情况需加上协变导数，1e-5是因为km->cm
                    if(chs[c]=='R' && chs[c2]=='T'){
                        val -= rcv_mu * u[2][ir] / dist * 1e-5;
                    }
                    else if(chs[c]=='T' && chs[c2]=='T'){
                        val += 2.0 * rcv_mu * u[1][ir] / dist * 1e-5;
                    }

                    res[c2][c][ir] = val;
                }
            }
        }
    }

    // 写入 nc 文件
    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        for(int c2=c; c2<GRT_CHANNEL_NUM; ++c2){
            NC_CHECK(NC_FUNC_REAL(nc_put_var) (in_ncid, out_varids[c2][c], res[c2][c]));
        }
    }
    

    // 关闭文件
    NC_CHECK(nc_close(in_ncid));

    GRT_SAFE_FREE_PTR(s_ingrid);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}