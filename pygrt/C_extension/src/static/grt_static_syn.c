/**
 * @file   grt_static_syn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 *    根据计算好的静态格林函数，定义震源机制以及方位角等，生成合成的静态三分量位移场
 * 
 */


#include "grt/common/const.h"
#include "grt/common/radiation.h"
#include "grt/common/coord.h"
#include "grt/common/util.h"
#include "grt/common/mynetcdf.h"

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    char *name;
    /** 输入 nc 格式的格林函数 */
    struct {
        bool active;
        char *s_ingrid;
    } G;
    /** 旋转到 Z, N, E */
    struct {
        bool active;
    } N;
    /** 放大系数 */
    struct {
        bool active;
        bool mult_src_mu;
        MYREAL M0;
        MYREAL src_mu;
    } S;  
    /** 剪切源 */
    struct {
        bool active;
    } M;
    /** 单力源 */
    struct {
        bool active;
    } F;
    /** 矩张量源 */
    struct {
        bool active;
    } T;
    /** 输出 nc 文件名 */
    struct {
        bool active;
        char *s_outgrid;
    } O;
    /** 是否计算空间导数 */
    struct {
        bool active;
    } e;

    // 存储不同震源的震源机制相关参数的数组
    MYREAL mchn[GRT_MECHANISM_NUM];

    // 方向因子数组
    MYREAL srcRadi[GRT_SRC_M_NUM][GRT_CHANNEL_NUM];

    // 最终要计算的震源类型
    MYINT computeType;
    char s_computeType[3];

} GRT_MODULE_CTRL;

/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    // G
    GRT_SAFE_FREE_PTR(Ctrl->G.s_ingrid);

    // O
    GRT_SAFE_FREE_PTR(Ctrl->O.s_outgrid);

    GRT_SAFE_FREE_PTR(Ctrl->name);
    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt static syn] %s\n\n", GRT_VERSION);printf(
"    Compute static displacement with the outputs of \n"
"    module `static_greenfn` , output to nc file.\n"
"    Three components are:\n"
"       + Up (Z),\n"
"       + Radial Outward (R),\n"
"       + Transverse Clockwise (T),\n"
"    and the units are cm. You can add -N to rotate ZRT to ZNE.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt static syn -G<ingrid> -S[u]<scale> \n"
"              [-M<strike>/<dip>/<rake>]\n"
"              [-T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>]\n"
"              [-F<fn>/<fe>/<fz>] \n"
"              [-N] [-e]\n"
"              -O<outgrid> \n"
"\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<ingrid>    Filepath to input nc Green's Functions grid.\n"
"\n"
"    -S[u]<scale>  Scale factor to all kinds of source. \n"
"                  + For Explosion, Shear and Moment Tensor,\n"
"                    unit of <scale> is dyne-cm. \n"
"                  + For Single Force, unit of <scale> is dyne.\n"
"                  + Since \"\\mu\" exists in scalar seismic moment\n"
"                    (\\mu*A*D), you can simply set -Su<scale>, <scale>\n"
"                    equals A*D (Area*Slip, [cm^3]), and <scale> will \n"
"                    multiply \\mu automatically in program.\n"
"\n"
"    -O<outgrid>   Filepath to output nc grid.\n"
"\n"
"    For source type, you can only set at most one of\n"
"    '-M', '-T' and '-F'. If none, an Explosion is used.\n"
"\n"
"    -M<strike>/<dip>/<rake>\n"
"                  Three angles to define a fault. \n"
"                  The angles are in degree.\n"
"\n"
"    -T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>\n"
"                  Six elements of Moment Tensor. \n"
"                  x (North), y (East), z (Downward).\n"
"                  Notice they will be scaled by <scale>.\n"
"\n"
"    -F<fn>/<fe>/<fz>\n"
"                  North, East and Vertical(Downward) Forces.\n"
"                  Notice they will be scaled by <scale>.\n"
"\n"
"    -N            Components of results will be Z, N, E.\n"
"\n"
"    -e            Compute the spatial derivatives, ui_z and ui_r,\n"
"                  of displacement u. In filenames, prefix \"r\" means \n"
"                  ui_r and \"z\" means ui_z. \n"
"\n"
"    -h            Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    Say you have computed Static Green's functions with following command:\n"
"        grt static greenfn -Mmilrow -D2/0 -X-5/5/10 -Y-5/5/10 -Ostgrn.nc\n"
"\n"
"    Then you can get static displacement of Explosion\n"
"        grt static syn -Gstgrn.nc -Su1e16 -Ostsyn_ex.nc\n"
"\n"
"    or Shear\n"
"        grt static syn -Gstgrn.nc -Su1e16 -M100/20/80 -Ostsyn_dc.nc\n"
"\n"
"    or Single Force\n"
"        grt static syn -Gstgrn.nc -S1e20 -F0.5/-1.2/3.3 -Ostsyn_sf.nc\n"
"\n"
"    or Moment Tensor\n"
"        grt static syn -Gstgrn.nc -Su1e16 -T2.3/0.2/-4.0/0.3/0.5/1.2 -Ostsyn_mt.nc\n"
"\n\n\n"
"\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    const char *command = Ctrl->name;

    // 先为个别参数设置非0初始值
    Ctrl->computeType = GRT_SYN_COMPUTE_EX;
    sprintf(Ctrl->s_computeType, "%s", "EX");

    int opt;
    while ((opt = getopt(argc, argv, ":G:O:S:M:F:T:Neh")) != -1) {
        switch (opt) {
            // 输入 nc 文件名
            case 'G':
                Ctrl->G.active = true;
                Ctrl->G.s_ingrid = strdup(optarg);
                break;

            // 输出 nc 文件名
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.s_outgrid = strdup(optarg);
                break;

            // 放大系数
            case 'S':
                Ctrl->S.active = true;
                {   
                    // 检查是否存在字符u，若存在表明需要乘上震源处的剪切模量
                    char *upos=NULL;
                    if((upos=strchr(optarg, 'u')) != NULL){
                        Ctrl->S.mult_src_mu = true;
                        *upos = ' ';
                    }
                }
                if(0 == sscanf(optarg, "%lf", &Ctrl->S.M0)){
                    GRTBadOptionError(command, S, "");
                };
                break;

            // 剪切震源
            case 'M':
                Ctrl->M.active = true;
                Ctrl->computeType = GRT_SYN_COMPUTE_DC;
                {
                    double strike, dip, rake;
                    sprintf(Ctrl->s_computeType, "%s", "DC");
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &strike, &dip, &rake)){
                        GRTBadOptionError(command, M, "");
                    };
                    if(strike < 0.0 || strike > 360.0){
                        GRTBadOptionError(command, M, "Strike must be in [0, 360].");
                    }
                    if(dip < 0.0 || dip > 90.0){
                        GRTBadOptionError(command, M, "Dip must be in [0, 90].");
                    }
                    if(rake < -180.0 || rake > 180.0){
                        GRTBadOptionError(command, M, "Rake must be in [-180, 180].");
                    }
                    Ctrl->mchn[0] = strike;
                    Ctrl->mchn[1] = dip;
                    Ctrl->mchn[2] = rake;
                }
                break;

            // 单力源
            case 'F':
                Ctrl->F.active = true;
                Ctrl->computeType = GRT_SYN_COMPUTE_SF;
                {
                    double fn, fe, fz;
                    sprintf(Ctrl->s_computeType, "%s", "SF");
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &fn, &fe, &fz)){
                        GRTBadOptionError(command, F, "");
                    };
                    Ctrl->mchn[0] = fn;
                    Ctrl->mchn[1] = fe;
                    Ctrl->mchn[2] = fz;
                }
                break;

            // 张量震源
            case 'T':
                Ctrl->T.active = true;
                Ctrl->computeType = GRT_SYN_COMPUTE_MT;
                {
                    double Mxx, Mxy, Mxz, Myy, Myz, Mzz;
                    sprintf(Ctrl->s_computeType, "%s", "MT");
                    if(6 != sscanf(optarg, "%lf/%lf/%lf/%lf/%lf/%lf", &Mxx, &Mxy, &Mxz, &Myy, &Myz, &Mzz)){
                        GRTBadOptionError(command, T, "");
                    };
                    Ctrl->mchn[0] = Mxx;
                    Ctrl->mchn[1] = Mxy;
                    Ctrl->mchn[2] = Mxz;
                    Ctrl->mchn[3] = Myy;
                    Ctrl->mchn[4] = Myz;
                    Ctrl->mchn[5] = Mzz;
                }
                break;

            // 是否计算位移空间导数, 影响 calcUTypes 变量
            case 'e':
                Ctrl->e.active = true;
                break;

            // 是否旋转到ZNE, 影响 rot2ZNE 变量
            case 'N':
                Ctrl->N.active = true;
                break;

            GRT_Common_Options_in_Switch(command, (char)(optopt));
        }
    }

    // 检查必选项有没有设置
    GRTCheckOptionSet(command, argc > 1);
    GRTCheckOptionActive(command, Ctrl, G);
    GRTCheckOptionActive(command, Ctrl, O);
    GRTCheckOptionActive(command, Ctrl, S);

    // 只能使用一种震源
    if(Ctrl->M.active + Ctrl->F.active + Ctrl->T.active > 1){
        GRTRaiseError("[%s] Error! Only support at most one of \"-M\", \"-F\" and \"-T\". Use \"-h\" for help.\n", command);
    }
}




/** 子模块主函数 */
int static_syn_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    Ctrl->name = strdup(argv[0]);
    const char *command = Ctrl->name;

    getopt_from_command(Ctrl, argc, argv);

    // 输出分量格式，即是否需要旋转到ZNE
    bool rot2ZNE = Ctrl->N.active;

    // 根据参数设置，选择分量名
    const char *chs = (rot2ZNE)? GRT_ZNE_CODES : GRT_ZRT_CODES;

    // nc 文件相关变量
    int in_ncid;
    int in_x_dimid, in_y_dimid;
    int in_x_varid, in_y_varid;
    const int ndims = 2;
    int in_u_varids[GRT_SRC_M_NUM][GRT_CHANNEL_NUM];
    int in_uiz_varids[GRT_SRC_M_NUM][GRT_CHANNEL_NUM];
    int in_uir_varids[GRT_SRC_M_NUM][GRT_CHANNEL_NUM];
    int out_ncid;
    int out_x_dimid, out_y_dimid;
    int out_x_varid, out_y_varid;
    int out_dimids[ndims];
    int out_syn_varids[GRT_CHANNEL_NUM];
    int out_syn_upar_varids[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
    
    // 打开 nc 文件
    GRT_NC_CHECK(nc_open(Ctrl->G.s_ingrid, NC_NOWRITE, &in_ncid));
    GRT_NC_CHECK(nc_create(Ctrl->O.s_outgrid, NC_CLOBBER, &out_ncid));

    // 读取全局属性，视情况计算 src_mu
    MYREAL src_va=0.0, src_vb=0.0, src_rho=0.0, src_mu=0.0;
    GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_get_att) (in_ncid, NC_GLOBAL, "src_va", &src_va));
    GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_get_att) (in_ncid, NC_GLOBAL, "src_vb", &src_vb));
    GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_get_att) (in_ncid, NC_GLOBAL, "src_rho", &src_rho));
    src_mu = src_vb*src_vb*src_rho*1e10;
    if(Ctrl->S.mult_src_mu) Ctrl->S.M0 *= src_mu;

    // 读入的数据是否有位移偏导
    MYINT calc_upar;
    GRT_NC_CHECK(GRT_NC_FUNC_MYINT(nc_get_att) (in_ncid, NC_GLOBAL, "calc_upar", &calc_upar));
    if(Ctrl->e.active && calc_upar == 0){
        GRTRaiseError("[%s] Input grid didn't have displacement derivatives, you can't set -e.", command);
    }

    // 复制属性
    GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_put_att) (out_ncid, NC_GLOBAL, "src_va", GRT_NC_MYREAL, 1, &src_va));
    GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_put_att) (out_ncid, NC_GLOBAL, "src_vb", GRT_NC_MYREAL, 1, &src_vb));
    GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_put_att) (out_ncid, NC_GLOBAL, "src_rho", GRT_NC_MYREAL, 1, &src_rho));
    GRT_NC_CHECK(GRT_NC_FUNC_MYINT(nc_put_att) (out_ncid, NC_GLOBAL, "calc_upar", GRT_NC_MYINT, 1, &calc_upar));
    {
        MYREAL rcv_va=0.0, rcv_vb=0.0, rcv_rho=0.0;
        GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_get_att) (in_ncid, NC_GLOBAL, "rcv_va", &rcv_va));
        GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_get_att) (in_ncid, NC_GLOBAL, "rcv_vb", &rcv_vb));
        GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_get_att) (in_ncid, NC_GLOBAL, "rcv_rho", &rcv_rho));
        GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_put_att) (out_ncid, NC_GLOBAL, "rcv_va", GRT_NC_MYREAL, 1, &rcv_va));
        GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_put_att) (out_ncid, NC_GLOBAL, "rcv_vb", GRT_NC_MYREAL, 1, &rcv_vb));
        GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_put_att) (out_ncid, NC_GLOBAL, "rcv_rho", GRT_NC_MYREAL, 1, &rcv_rho)); 
    }

    // 是否旋转到ZNE记录到全局属性
    {
        MYINT rot2ZNE_int = rot2ZNE;
        GRT_NC_CHECK(GRT_NC_FUNC_MYINT(nc_put_att) (out_ncid, NC_GLOBAL, "rot2ZNE", GRT_NC_MYINT, 1, &rot2ZNE_int));
    }

    // 震源类型写入全局属性
    GRT_NC_CHECK(nc_put_att_text(out_ncid, NC_GLOBAL, "computeType", strlen(Ctrl->s_computeType), Ctrl->s_computeType));
    
    // 读入坐标变量 dimid, varid
    size_t nx, ny;
    GRT_NC_CHECK(nc_inq_dimid(in_ncid, "north", &in_x_dimid));
    GRT_NC_CHECK(nc_inq_dimlen(in_ncid, in_x_dimid, &nx));
    GRT_NC_CHECK(nc_inq_dimid(in_ncid, "east", &in_y_dimid));
    GRT_NC_CHECK(nc_inq_dimlen(in_ncid, in_y_dimid, &ny));

    // 写入坐标变量 dimid, varid
    GRT_NC_CHECK(nc_def_dim(out_ncid, "north", nx, &out_x_dimid));
    GRT_NC_CHECK(nc_def_dim(out_ncid, "east", ny, &out_y_dimid));
    GRT_NC_CHECK(nc_def_var(out_ncid, "north", GRT_NC_MYREAL, 1, &out_x_dimid, &out_x_varid));
    GRT_NC_CHECK(nc_def_var(out_ncid, "east",  GRT_NC_MYREAL, 1, &out_y_dimid, &out_y_varid));
    out_dimids[0] = out_x_dimid;
    out_dimids[1] = out_y_dimid;

    // 读入格林函数 varid
    for(int i=0; i<GRT_SRC_M_NUM; ++i){
        int modr = GRT_SRC_M_ORDERS[i];
        char *s_title = NULL;
        for(int c=0; c<GRT_CHANNEL_NUM; ++c){
            if(modr==0 && GRT_ZRT_CODES[c]=='T')  continue;

            GRT_SAFE_ASPRINTF(&s_title, "%s%c", GRT_SRC_M_NAME_ABBR[i], GRT_ZRT_CODES[c]);
            GRT_NC_CHECK(nc_inq_varid(in_ncid, s_title, &in_u_varids[i][c]));

            // 位移偏导
            if(Ctrl->e.active){
                GRT_SAFE_ASPRINTF(&s_title, "z%s%c", GRT_SRC_M_NAME_ABBR[i], GRT_ZRT_CODES[c]);
                GRT_NC_CHECK(nc_inq_varid(in_ncid, s_title, &in_uiz_varids[i][c]));
                GRT_SAFE_ASPRINTF(&s_title, "r%s%c", GRT_SRC_M_NAME_ABBR[i], GRT_ZRT_CODES[c]);
                GRT_NC_CHECK(nc_inq_varid(in_ncid, s_title, &in_uir_varids[i][c]));
            }
        }
        GRT_SAFE_FREE_PTR(s_title);
    }

    // 定义合成结果 varid
    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        char *s_title = NULL;
        GRT_SAFE_ASPRINTF(&s_title, "%c", toupper(chs[c]));
        GRT_NC_CHECK(nc_def_var(out_ncid, s_title, GRT_NC_MYREAL, ndims, out_dimids, &out_syn_varids[c]));
        // 位移偏导
        if(Ctrl->e.active){
            for(int c2=0; c2<GRT_CHANNEL_NUM; ++c2){
                GRT_SAFE_ASPRINTF(&s_title, "%c%c", tolower(chs[c2]), toupper(chs[c]));
                GRT_NC_CHECK(nc_def_var(out_ncid, s_title, GRT_NC_MYREAL, ndims, out_dimids, &out_syn_upar_varids[c2][c]));
            }
        }
        GRT_SAFE_FREE_PTR(s_title);
    }

    // 结束定义模式
    GRT_NC_CHECK(nc_enddef(out_ncid));

    // 读取坐标变量
    MYREAL *xs = (MYREAL *)calloc(nx, sizeof(MYREAL));
    MYREAL *ys = (MYREAL *)calloc(ny, sizeof(MYREAL));
    GRT_NC_CHECK(nc_inq_varid(in_ncid, "north", &in_x_varid));
    GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_get_var) (in_ncid, in_x_varid, xs));
    GRT_NC_CHECK(nc_inq_varid(in_ncid, "east", &in_y_varid));
    GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_get_var) (in_ncid, in_y_varid, ys));

    // 写入坐标变量
    GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_put_var) (out_ncid, out_x_varid, xs));
    GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_put_var) (out_ncid, out_y_varid, ys));

    // 总震中距数
    size_t nr = nx * ny;

    // 先将所有格林函数及其偏导读入内存，
    // 否则连续使用 nc_grt_var1 式读入效率太慢
    MYREAL *u[GRT_SRC_M_NUM][GRT_CHANNEL_NUM];
    MYREAL *uiz[GRT_SRC_M_NUM][GRT_CHANNEL_NUM];
    MYREAL *uir[GRT_SRC_M_NUM][GRT_CHANNEL_NUM];
    for(int i=0; i<GRT_SRC_M_NUM; ++i){
        int modr = GRT_SRC_M_ORDERS[i];
        for(int c=0; c<GRT_CHANNEL_NUM; ++c){
            // 先申请全 0 内存
            u[i][c] = (MYREAL *)calloc(nr, sizeof(MYREAL));
            uiz[i][c] = (MYREAL *)calloc(nr, sizeof(MYREAL));
            uir[i][c] = (MYREAL *)calloc(nr, sizeof(MYREAL));

            if(modr==0 && GRT_ZRT_CODES[c]=='T')  continue;

            GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_get_var) (in_ncid, in_u_varids[i][c], u[i][c]));

            if(Ctrl->e.active){
                GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_get_var) (in_ncid, in_uiz_varids[i][c], uiz[i][c]));
                GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_get_var) (in_ncid, in_uir_varids[i][c], uir[i][c]));
            }
        }
    }
    
    // 最终计算的结果
    MYREAL (*syn)[GRT_CHANNEL_NUM] = (MYREAL (*)[GRT_CHANNEL_NUM])calloc(nr, sizeof(MYREAL)*GRT_CHANNEL_NUM);
    MYREAL (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM] = (MYREAL (*)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])calloc(nr, sizeof(MYREAL)*GRT_CHANNEL_NUM*GRT_CHANNEL_NUM);
    
    // 每个点逐个处理
    for(size_t ix=0; ix < nx; ++ix){
        MYREAL x = xs[ix];
        for(size_t iy=0; iy < ny; ++iy){
            MYREAL y = ys[iy];

            size_t ir = iy + ix*ny;

            // 方位角
            MYREAL azrad = atan2(y, x);

            // 震中距
            MYREAL dist = GRT_MAX(sqrt(x*x + y*y), GRT_MIN_DISTANCE);

            // 计算和位移相关量的种类（1-位移，2-ui_z，3-ui_r，4-ui_t）
            int calcUTypes = (Ctrl->e.active)? 4 : 1;
            MYREAL upar_scale = 1.0;

            MYREAL *(*up)[GRT_CHANNEL_NUM];  // 使用对应类型的格林函数
            MYREAL tmpsyn[GRT_CHANNEL_NUM];

            for(int ityp=0; ityp<calcUTypes; ++ityp){

                upar_scale=1.0;

                // 求位移空间导数时，需调整比例系数
                if(ityp > 0){
                    switch (GRT_ZRT_CODES[ityp-1]){
                        // 合成 ui_z, uir
                        case 'Z': case 'R': upar_scale = 1e-5; break;
                        // 合成 ui_t
                        case 'T': upar_scale = 1e-5 / dist; break;
                        default: break;
                    }
                }
                

                if(ityp==1){
                    up = uiz;
                } else if(ityp==2){
                    up = uir;
                } else {
                    up = u;
                }

                memset(tmpsyn, 0, sizeof(MYREAL)*GRT_CHANNEL_NUM);

                // 计算震源辐射因子
                grt_set_source_radiation(Ctrl->srcRadi, Ctrl->computeType, (ityp > 0) && GRT_ZRT_CODES[ityp-1]=='T', Ctrl->S.M0, upar_scale, azrad, Ctrl->mchn);

                // 合成
                for(int i=0; i<GRT_SRC_M_NUM; ++i){
                    int modr = GRT_SRC_M_ORDERS[i];
                    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
                        if(modr==0 && GRT_ZRT_CODES[c]=='T')  continue;

                        tmpsyn[c] += up[i][c][ir] * Ctrl->srcRadi[i][c];
                    }
                }

                // 记录数据
                for(int i=0; i<GRT_CHANNEL_NUM; ++i){
                    if(ityp == 0){
                        syn[ir][i] = tmpsyn[i];
                    } else {
                        syn_upar[ir][ityp-1][i] = tmpsyn[i];
                    }
                }

            } // END loop calcUTypes

            // 是否要转到ZNE
            if(rot2ZNE){
                if(Ctrl->e.active){
                    grt_rot_zrt2zxy_upar(azrad, syn[ir], syn_upar[ir], dist*1e5);
                } else {
                    grt_rot_zxy2zrt_vec(-azrad, syn[ir]);
                }
            }

        }
    }

    // 写入 nc 文件
    MYREAL *tmpdata = (MYREAL *)calloc(nr, sizeof(MYREAL));
    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        for(size_t ir=0; ir < nr; ++ir){
            tmpdata[ir] = syn[ir][c];
        }
        GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_put_var) (out_ncid, out_syn_varids[c], tmpdata));

        // 位移偏导
        if(Ctrl->e.active){
            for(int c2=0; c2<GRT_CHANNEL_NUM; ++c2){
                for(size_t ir=0; ir < nr; ++ir){
                    tmpdata[ir] = syn_upar[ir][c2][c];
                }
                GRT_NC_CHECK(GRT_NC_FUNC_MYREAL(nc_put_var) (out_ncid, out_syn_upar_varids[c2][c], tmpdata));
            }
        }
    }
    GRT_SAFE_FREE_PTR(tmpdata);
    

    // 关闭文件
    GRT_NC_CHECK(nc_close(in_ncid));
    GRT_NC_CHECK(nc_close(out_ncid));

    // 释放内存
    for(int i=0; i<GRT_SRC_M_NUM; ++i){
        for(int c=0; c<GRT_CHANNEL_NUM; ++c){
            GRT_SAFE_FREE_PTR(u[i][c]);
            GRT_SAFE_FREE_PTR(uiz[i][c]);
            GRT_SAFE_FREE_PTR(uir[i][c]);
        }
    }
    GRT_SAFE_FREE_PTR(syn);
    GRT_SAFE_FREE_PTR(syn_upar);


    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}