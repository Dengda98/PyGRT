/**
 * @file   grt_static_stress.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-08
 * 
 *    根据预先合成的静态位移空间导数，组合成静态应力张量
 * 
 */

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    int dummy;
} GRT_MODULE_CTRL;

/**
 * 释放模块控制结构体
 *
 * @param[in,out]  Ctrl   模块控制结构体
 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt static stress] %s\n\n", GRT_VERSION);printf(
"    Compute the static stress tensor from a static syn output file\n"
"    (unit: dyne/cm^2 = 0.1 Pa) and write the result back to the\n"
"    same nc file. The input must be generated with -e.\n"
"    Output variables are stress_ZZ, stress_ZR, stress_ZT, stress_RR,\n"
"    stress_RT, and stress_TT for ZRT, or stress_ZZ, stress_ZN,\n"
"    stress_ZE, stress_NN, stress_NE, and stress_EE for ZNE.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt static stress <ingrid> [-h]\n"
"\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt static syn -Gstgrn.nc -Su1e16 -e -Ostsyn.nc\n"
"    grt static stress stsyn.nc\n"
"\n\n\n"
);
}


/**
 * 从命令行中读取选项
 *
 * @param[out]  Ctrl   模块控制结构体
 * @param[in]   argc   命令行参数数量
 * @param[in]   argv   命令行参数数组
 */
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

/**
 * 由静态位移偏导合成应力张量
 *
 * @param[in]   npts       接收点数量
 * @param[in]   norths     接收点 North 坐标 (km)
 * @param[in]   easts      接收点 East 坐标 (km)
 * @param[in]   u          位移数组
 * @param[in]   upar       位移偏导数组
 * @param[out]  res        应力张量数组
 * @param[in]   rot2ZNE    是否使用 ZNE 分量
 * @param[in]   mu         各接收点剪切模量
 * @param[in]   lam        各接收点第一拉梅常数
 */
static void compute_stress(
    size_t npts, const real_t *norths, const real_t *easts,
    real_t *const u[GRT_CHANNEL_NUM],
    real_t *const upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    real_t *const res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    bool rot2ZNE, const real_t *mu, const real_t *lam)
{
    const char *chs = rot2ZNE ? GRT_ZNE_CODES : GRT_ZRT_CODES;

    for(size_t ir=0; ir<npts; ++ir){
        real_t dist = hypot(norths[ir], easts[ir]);
        // 联络项（1e-5: km→cm）：r≠0 用 u/r；r=0 改用 ∂_r u，与 syn 中 (1/r)∂_θ 有限部分配套
        real_t ur_over_r = GRT_IS_ZERO(dist) ? upar[1][1][ir] : (u[1][ir] / dist * 1e-5);
        real_t ut_over_r = GRT_IS_ZERO(dist) ? upar[1][2][ir] : (u[2][ir] / dist * 1e-5);

        real_t mu_i = mu[ir];
        real_t lam_i = lam[ir];
        real_t lam_ukk = upar[0][0][ir] + upar[1][1][ir] + upar[2][2][ir];
        if(!rot2ZNE)  lam_ukk += ur_over_r;
        lam_ukk *= lam_i;

        for(int c=0; c<GRT_CHANNEL_NUM; ++c){
            for(int c2=c; c2<GRT_CHANNEL_NUM; ++c2){
                real_t val = mu_i * (upar[c2][c][ir] + upar[c][c2][ir]);
                if(c == c2)  val += lam_ukk;
                if(chs[c]=='R' && chs[c2]=='T'){
                    val -= mu_i * ut_over_r;
                }
                else if(chs[c]=='T' && chs[c2]=='T'){
                    val += 2.0 * mu_i * ur_over_r;
                }
                res[c2][c][ir] = val;
            }
        }
    }
}


/**
 * 由全局属性读取标量 rcv_va、rcv_vb、rcv_rho 并填充模量数组
 *
 * @param[in]   ncid   NetCDF 文件 ID
 * @param[in]   npts   接收点数量
 * @param[out]  mu     各接收点剪切模量
 * @param[out]  lam    各接收点第一拉梅常数
 */
static void fill_mu_lam_from_global_atts(int ncid, size_t npts, real_t *mu, real_t *lam)
{
    real_t rcv_va=0.0, rcv_vb=0.0, rcv_rho=0.0;
    NC_CHECK(NC_FUNC_REAL(nc_get_att) (ncid, NC_GLOBAL, "rcv_va", &rcv_va));
    NC_CHECK(NC_FUNC_REAL(nc_get_att) (ncid, NC_GLOBAL, "rcv_vb", &rcv_vb));
    NC_CHECK(NC_FUNC_REAL(nc_get_att) (ncid, NC_GLOBAL, "rcv_rho", &rcv_rho));
    real_t rcv_mu  = rcv_vb*rcv_vb*rcv_rho*1e10;
    real_t rcv_lam = rcv_va*rcv_va*rcv_rho*1e10 - 2.0*rcv_mu;
    for(size_t i = 0; i < npts; ++i){
        mu[i]  = rcv_mu;
        lam[i] = rcv_lam;
    }
}


/** 子模块主函数 */
int static_stress_main(int argc, char **argv){
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

    // 识别 grid / points 布局，并将坐标展平供统一计算
    GRT_RECV_NC_INFO recv_info;
    grt_recv_nc_info_load(in_ncid, &recv_info);
    size_t npts = recv_info.npts;
    real_t *norths_flat = recv_info.norths;
    real_t *easts_flat = recv_info.easts;
    int out_ndims = (recv_info.layout == GRT_RECV_NC_LAYOUT_POINTS) ? 1 : 2;
    int out_dimids[2] = {recv_info.dimids[0], recv_info.dimids[1]};

    // 逐点物性参数 mu/lam
    real_t *mu  = (real_t *)calloc(npts, sizeof(real_t));
    real_t *lam = (real_t *)calloc(npts, sizeof(real_t));
    if(recv_info.layout == GRT_RECV_NC_LAYOUT_POINTS){
        int va_varid;
        if(nc_inq_varid(in_ncid, "rcv_va", &va_varid) == NC_NOERR){
            int vb_varid, rho_varid;
            real_t *rcv_va  = (real_t *)calloc(npts, sizeof(real_t));
            real_t *rcv_vb  = (real_t *)calloc(npts, sizeof(real_t));
            real_t *rcv_rho = (real_t *)calloc(npts, sizeof(real_t));
            NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, va_varid, rcv_va));
            NC_CHECK(nc_inq_varid(in_ncid, "rcv_vb", &vb_varid));
            NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, vb_varid, rcv_vb));
            NC_CHECK(nc_inq_varid(in_ncid, "rcv_rho", &rho_varid));
            NC_CHECK(NC_FUNC_REAL(nc_get_var) (in_ncid, rho_varid, rcv_rho));
            for(size_t i = 0; i < npts; ++i){
                mu[i]  = rcv_vb[i]*rcv_vb[i]*rcv_rho[i]*1e10;
                lam[i] = rcv_va[i]*rcv_va[i]*rcv_rho[i]*1e10 - 2.0*mu[i];
            }
            GRT_SAFE_FREE_PTR(rcv_va);
            GRT_SAFE_FREE_PTR(rcv_vb);
            GRT_SAFE_FREE_PTR(rcv_rho);
        } else {
            // 无逐点变量时回退到全局属性
            fill_mu_lam_from_global_atts(in_ncid, npts, mu, lam);
        }
    } else {
        fill_mu_lam_from_global_atts(in_ncid, npts, mu, lam);
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
            // 这里命名顺序要注意，例如 ZR -> mu*(u_{z,r} + u_{r,z})
            GRT_SAFE_ASPRINTF(&s_title, "stress_%c%c", toupper(chs[c]), toupper(chs[c2]));
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

    compute_stress(npts, norths_flat, easts_flat, u, upar, res, rot2ZNE, mu, lam);

    // 写入 nc 文件
    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        for(int c2=c; c2<GRT_CHANNEL_NUM; ++c2){
            NC_CHECK(NC_FUNC_REAL(nc_put_var) (in_ncid, out_varids[c2][c], res[c2][c]));
        }
    }

    // 关闭文件
    NC_CHECK(nc_close(in_ncid));

    grt_recv_nc_info_free(&recv_info);
    GRT_SAFE_FREE_PTR(mu);
    GRT_SAFE_FREE_PTR(lam);
    GRT_SAFE_FREE_PTR(s_ingrid);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
