/**
 * @file   grt_syn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12-2
 * 
 *    根据计算好的格林函数，定义震源机制以及方位角等，生成合成的三分量地震图
 * 
 */

#include "grt.h"

// 防止被替换为虚数单位
#undef I

/** 该子模块的参数控制结构体 */
typedef struct {
    /** 格林函数路径 */
    struct {
        bool active;
        char *s_grnpath;
    } G;
    /** 输出目录 */
    struct {
        bool active;
        char *s_output_dir;
    } O;
    /** 方位角 */
    struct {
        bool active;
        real_t azimuth;
        real_t azrad;
        real_t backazimuth;
    } A;
    /** 旋转到 Z, N, E */
    struct {
        bool active;
    } N;
    /** 放大系数 */
    struct {
        bool active;
        bool mult_src_mu;
        real_t M0;
        real_t src_mu;
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
    /** 积分次数 */
    struct {
        bool active;
        int int_times;
    } I;
    /** 求导次数 */
    struct {
        bool active;
        int dif_times;
    } J;
    /** 时间函数 */
    struct {
        bool active;
        char tftype;
        char *tfparams;
    } D;
    /** 根目录检索时的震源和台站深度 */
    struct {
        bool s_active;
        bool r_active;
        real_t depsrc;
        real_t deprcv;
    } Depth;
    /** 根目录检索时的震中距 */
    struct {
        bool active;
        real_t dist;
    } R;
    /** 静默输出 */
    struct {
        bool active;
    } s;
    /** 是否计算空间导数 */
    struct {
        bool active;
    } e;

    // 存储不同震源的震源机制相关参数的数组
    real_t mchn[GRT_MECHANISM_NUM];

    // 震中距
    real_t dist;

    // 震源层 Vp/Vs
    real_t VpVs_ratio;

    // 方向因子数组
    realChnlGrid srcRadi;

    // 最终要计算的震源类型
    GRT_SYN_TYPE computeType;
    char s_computeType[3];

} GRT_MODULE_CTRL;

static void resolve_syn_grn_path(GRT_MODULE_CTRL *Ctrl);


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    // G
    GRT_SAFE_FREE_PTR(Ctrl->G.s_grnpath);
    // O
    GRT_SAFE_FREE_PTR(Ctrl->O.s_output_dir);
    // D
    GRT_SAFE_FREE_PTR(Ctrl->D.tfparams);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt syn] %s\n\n", GRT_VERSION);printf(
"    A Supplementary Tool of GRT to Compute Three-Component \n"
"    Displacement with the outputs of module `greenfn`.\n"
"    Three components are:\n"
"       + Up (Z),\n"
"       + Radial Outward (R),\n"
"       + Transverse Clockwise (T),\n"
"    and the units are cm. You can add -N to rotate ZRT to ZNE.\n"
"\n"
"    + Default outputs (without -I and -J) are impulse-like displacements.\n"
"    + -D, -I and -J are applied in the time domain.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt syn -G<grn_path> -A<azimuth> -S[u]<scale> -O<outdir> \n"
"            [-Ds<depsrc>] [-Dr<deprcv>] [-R<dist>]\n"
"            [-M<strike>/<dip>[/<rake>]]\n"
"            [-T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>]\n"
"            [-F<fn>/<fe>/<fz>] \n"
"            [-D<tftype>/<tfparams>] [-I<odr>] [-J<odr>]\n" 
"            [-N] [-e] [-s]\n"
"\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<grn_path>  Green's Functions output directory of module `greenfn`.\n"
"                  A single-distance subdirectory may be given directly.\n"
"                  In that mode, -Ds, -Dr and -R cannot be used.\n"
"                  When the output root directory is given, each of -Ds, -Dr\n"
"                  and -R is required when the corresponding library dimension\n"
"                  has multiple values, and optional when it has one value.\n"
"                  Explicit selectors must exactly match the corresponding\n"
"                  values in the library; no depth or distance interpolation\n"
"                  is done.\n"

"    -Ds<depsrc>   Source depth (km) used with a Green's-function root.\n"
"                  Required for multiple source depths and optional for one.\n"
"                  It cannot be used when -G points to a subdirectory.\n"

"    -Dr<deprcv>   Receiver depth (km) used with a Green's-function root.\n"
"                  Required for multiple receiver depths and optional for one.\n"
"                  It cannot be used when -G points to a subdirectory.\n"

"    -R<dist>      Exact epicentral distance (km) used with a root directory.\n"
"                  Required for multiple distances and optional for one.\n"
"                  It cannot be used when -G points to a subdirectory.\n"
"\n"
"    -A<azimuth>   Azimuth in degree, from source to station.\n"
"                  Ignored (forced to 0°) when Green's Functions\n"
"                  have zero epicentral distance.\n"
"\n"
"    -S[u]<scale>  Scale factor to all kinds of source. \n"
"                  + For Explosion, Shear and Moment Tensor,\n"
"                    unit of <scale> is dyne-cm.\n"
"                  + For Single Force, unit of <scale> is dyne.\n"
"                  + Since \"\\mu\" exists in scalar seismic moment\n"
"                    (\\mu*A*D), you can simply set -Su<scale>, <scale>\n"
"                    equals A*D (Area*Slip, [cm^3]), and <scale> will \n"
"                    multiply \\mu automatically in program.\n"
"\n"
"    For source type, you can only set at most one of\n"
"    '-M', '-T' and '-F'. If none, an Explosion is used.\n"
"\n"
"    -M<strike>/<dip>[/<rake>]\n"
"                  Three angles to define a shear fault. \n"
"                  The angles are in degree.\n"
"                  If <rake> not set, then define a tensile fault.\n"
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
"    -O<outdir>    Directory of output for saving. Default is\n"
"                  current directory.\n"
"\n"
"    -D<tftype>/<tfparams>\n"
"                  Convolve a Time Function. All time functions use area\n"
"                  normalization except Ricker wavelet, which has a peak\n"
"                  amplitude of 1.0.\n"
"                  There are several options:\n"
"                  + Parabolic wave (y = a*x^2 + b*x)\n"
"                    set -D%c/<t0>, <t0> (secs) is the duration of wave.\n", GRT_SIG_PARABOLA); printf(
"                    e.g. \n"
"                         -D%c/1.3\n", GRT_SIG_PARABOLA); printf(
"                  + Trapezoidal wave\n"
"                    set -D%c/<t1>/<t2>/<t3>, <t1> is the end time of\n", GRT_SIG_TRAPEZOID); printf(
"                    Rising, <t2> is the end time of Platform, and\n"
"                    <t3> is the end time of Falling.\n"
"                    e.g. \n"
"                         -D%c/0.1/0.2/0.4\n", GRT_SIG_TRAPEZOID); printf(
"                         -D%c/0.4/0.4/0.6 (become a triangle)\n", GRT_SIG_TRAPEZOID); printf(
"                  + Ricker wavelet\n"
"                    set -D%c/<f0>, <f0> (Hz) is the dominant frequency.\n", GRT_SIG_RICKER); printf(
"                    e.g. \n"
"                         -D%c/0.5 \n", GRT_SIG_RICKER); printf(
"                  + Custom wave\n"
"                    set -D%c/<path>, <path> is the filepath to a custom\n", GRT_SIG_CUSTOM); printf(
"                    Time Function ASCII file. The file has just one column\n"
"                    of amplitude and no other columns. Its sequence sum should\n"
"                    be 1.0; the program only issues a warning when it is not.\n"
"                    The file can contain unlimited comment lines with prefix\n"
"                    \"#\".\n"
"                    e.g. \n"
"                         -D%c/tfunc.txt \n", GRT_SIG_CUSTOM); printf(
"                  To match the time interval in Green's Functions, \n"
"                  parameters of Time Function will be slightly modified.\n"
"                  The corresponding Time Function will be saved\n"
"                  as a SAC file under <outdir>.\n"
"\n"
"    -I<odr>       Order of integration. Default not use\n"
"\n"
"    -J<odr>       Order of differentiation. Default not use\n"
"\n"
"    -N            Components of results will be Z, N, E.\n"
"\n"
"    -e            Compute the spatial derivatives, ui_z and ui_r,\n"
"                  of displacement u. In filenames, prefix \"r\" means \n"
"                  ui_r and \"z\" means ui_z. \n"
"\n"
"    -s            Silence all outputs.\n"
"\n"
"    -h            Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    Say you have computed Green's functions with following command:\n"
"        grt greenfn -Mmilrow -N1000/0.01 -D2/0 -Ores -R2,4,6,8,10\n"
"\n"
"    Then you can get synthetic seismograms of Explosion at epicentral\n"
"    distance of 10 km and an azimuth of 30° by running:\n"
"        grt syn -Gres/milrow_2_0_10 -Osyn_ex -A30 -S1e24\n"
"\n"
"    or Shear\n"
"        grt syn -Gres/milrow_2_0_10 -Osyn_dc -A30 -S1e24 -M100/20/80\n"
"\n"
"    or Tension\n"
"        grt syn -Gres/milrow_2_0_10 -Osyn_dc -A30 -S1e24 -M100/20\n"
"\n"
"    or Single Force\n"
"        grt syn -Gres/milrow_2_0_10 -Osyn_sf -A30 -S1e24 -F0.5/-1.2/3.3\n"
"\n"
"    or Moment Tensor\n"
"        grt syn -Gres/milrow_2_0_10 -Osyn_mt -A30 -S1e24 -T2.3/0.2/-4.0/0.3/0.5/1.2\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    // 先为个别参数设置非0初始值
    Ctrl->computeType = GRT_SYN_EX;
    sprintf(Ctrl->s_computeType, "%s", "EX");

    int opt;
    while ((opt = getopt(argc, argv, ":G:A:S:M:F:T:O:D:I:J:R:Nehs")) != -1) {
        switch (opt) {
            // 格林函数路径
            case 'G':
                Ctrl->G.active = true;
                Ctrl->G.s_grnpath = strdup(optarg);
                // 检查是否存在该目录
                GRTCheckDirExist(Ctrl->G.s_grnpath);
                break;

            // 方位角
            case 'A':
                Ctrl->A.active = true;
                if(0 == sscanf(optarg, "%lf", &Ctrl->A.azimuth)){
                    GRTBadOptionError(A, "");
                };
                if(Ctrl->A.azimuth < 0.0 || Ctrl->A.azimuth > 360.0){
                    GRTBadOptionError(A, "Azimuth must be in [0, 360].");
                }
                Ctrl->A.backazimuth = 180.0 + Ctrl->A.azimuth;
                if(Ctrl->A.backazimuth >= 360.0)   Ctrl->A.backazimuth -= 360.0;
                Ctrl->A.azrad = Ctrl->A.azimuth * DEG1;
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
                    GRTBadOptionError(S, "");
                };
                break;
            
            // 剪切震源， 张裂源
            case 'M':
                Ctrl->M.active = true;
                {
                    real_t strike=0.0, dip=0.0, rake=0.0;
                    int nscan = sscanf(optarg, "%lf/%lf/%lf", &strike, &dip, &rake);
                    if(nscan >= 2){
                        Ctrl->computeType = GRT_SYN_TS;
                        sprintf(Ctrl->s_computeType, "%s", "TS");
                        if(strike < 0.0 || strike > 360.0){
                            GRTBadOptionError(M, "Strike must be in [0, 360].");
                        }
                        if(dip < 0.0 || dip > 90.0){
                            GRTBadOptionError(M, "Dip must be in [0, 90].");
                        }
                        if(nscan == 3){
                            Ctrl->computeType = GRT_SYN_DC;
                            sprintf(Ctrl->s_computeType, "%s", "DC");
                            if(rake < -180.0 || rake > 180.0){
                                GRTBadOptionError(M, "Rake must be in [-180, 180].");
                            }
                        }
                    } else {
                        GRTBadOptionError(M, "");
                    };
                    
                    
                    Ctrl->mchn[0] = strike;
                    Ctrl->mchn[1] = dip;
                    Ctrl->mchn[2] = rake;
                }
                break;

            // 单力源
            case 'F':
                Ctrl->F.active = true;
                Ctrl->computeType = GRT_SYN_SF;
                {
                    real_t fn, fe, fz;
                    sprintf(Ctrl->s_computeType, "%s", "SF");
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &fn, &fe, &fz)){
                        GRTBadOptionError(F, "");
                    };
                    Ctrl->mchn[0] = fn;
                    Ctrl->mchn[1] = fe;
                    Ctrl->mchn[2] = fz;
                }
                break;

            // 张量震源
            case 'T':
                Ctrl->T.active = true;
                Ctrl->computeType = GRT_SYN_MT;
                {
                    real_t Mxx, Mxy, Mxz, Myy, Myz, Mzz;
                    sprintf(Ctrl->s_computeType, "%s", "MT");
                    if(6 != sscanf(optarg, "%lf/%lf/%lf/%lf/%lf/%lf", &Mxx, &Mxy, &Mxz, &Myy, &Myz, &Mzz)){
                        GRTBadOptionError(T, "");
                    };
                    Ctrl->mchn[0] = Mxx;
                    Ctrl->mchn[1] = Mxy;
                    Ctrl->mchn[2] = Mxz;
                    Ctrl->mchn[3] = Myy;
                    Ctrl->mchn[4] = Myz;
                    Ctrl->mchn[5] = Mzz;
                }
                break;

            // 输出路径
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.s_output_dir = strdup(optarg);
                break;

            // 卷积时间函数
            case 'D':
                if(optarg[0] == 's'){
                    Ctrl->Depth.s_active = true;
                    if(strchr(optarg + 1, '/') != NULL || 1 != sscanf(optarg + 1, "%lf", &Ctrl->Depth.depsrc)){
                        GRTBadOptionError(Ds, "");
                    }
                    if(Ctrl->Depth.depsrc < 0.0){
                        GRTBadOptionError(Ds, "Negative source depth is not supported.");
                    }
                } else if(optarg[0] == 'r' && optarg[1] != '/'){
                    Ctrl->Depth.r_active = true;
                    if(strchr(optarg + 1, '/') != NULL || 1 != sscanf(optarg + 1, "%lf", &Ctrl->Depth.deprcv)){
                        GRTBadOptionError(Dr, "");
                    }
                    if(Ctrl->Depth.deprcv < 0.0){
                        GRTBadOptionError(Dr, "Negative receiver depth is not supported.");
                    }
                } else {
                    Ctrl->D.active = true;
                    Ctrl->D.tfparams = (char*)malloc(sizeof(char) * (strlen(optarg) + 1));
                    if(optarg[1] != '/' || 1 != sscanf(optarg, "%c", &Ctrl->D.tftype) || 1 != sscanf(optarg + 2, "%s", Ctrl->D.tfparams)){
                        GRTBadOptionError(D, "");
                    }
                    if(! grt_check_tftype_tfparams(Ctrl->D.tftype, Ctrl->D.tfparams)){
                        GRTBadOptionError(D, "");
                    }
                }
                break;

            // 根目录检索时指定震中距
            case 'R':
                Ctrl->R.active = true;
                if(1 != sscanf(optarg, "%lf", &Ctrl->R.dist) || Ctrl->R.dist < 0.0){
                    GRTBadOptionError(R, "Nonnegative epicentral distance is required.");
                }
                break;

            // 对结果做积分
            case 'I':
                Ctrl->I.active = true;
                if(1 != sscanf(optarg, "%d", &Ctrl->I.int_times)){
                    GRTBadOptionError(I, "");
                }
                if(Ctrl->I.int_times <= 0){
                    GRTBadOptionError(I, "Order should be positive.");
                }
                break;

            // 对结果做微分
            case 'J':
                Ctrl->J.active = true;
                if(1 != sscanf(optarg, "%d", &Ctrl->J.dif_times)){
                    GRTBadOptionError(J, "");
                }
                if(Ctrl->J.dif_times <= 0){
                    GRTBadOptionError(J, "Order should be positive.");
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

            // 不打印在终端
            case 's':
                Ctrl->s.active = true;
                break;

            GRT_Common_Options_in_Switch((char)(optopt));
        }

    }

    // 检查必选项有没有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, G);
    GRTCheckOptionActive(Ctrl, A);
    GRTCheckOptionActive(Ctrl, S);
    GRTCheckOptionActive(Ctrl, O);

    // 只能使用一种震源
    if(Ctrl->M.active + Ctrl->F.active + Ctrl->T.active > 1){
        GRTRaiseError("Only support at most one of \"-M\", \"-F\" and \"-T\". Use \"-h\" for help.\n");
    }

    resolve_syn_grn_path(Ctrl);

    // 建立保存目录
    GRTCheckMakeDir(Ctrl->O.s_output_dir);

    // 随机读取一个 sac，确定 dist 和 src_mu
    {
        struct dirent *entry;
        DIR *dp = opendir(Ctrl->G.s_grnpath);
        while ((entry = readdir(dp))) {
            if (strlen(entry->d_name) <= 4)  continue;
            if (strcmp(entry->d_name + strlen(entry->d_name) - 3, "sac") != 0)  continue;

            char *s_filepath = NULL;
            GRT_SAFE_ASPRINTF(&s_filepath, "%s/%s", Ctrl->G.s_grnpath, entry->d_name);
            SACTRACE *sac = grt_read_SACTRACE(s_filepath, true);
            GRT_SAFE_FREE_PTR(s_filepath);

            Ctrl->dist = sac->hd.dist;

            if(!Ctrl->s.active){
                GRTRaiseInfo(
                    "Selected Green's function: depsrc = %.6g km, deprcv = %.6g km, dist = %.6g km.",
                    sac->hd.evdp, -sac->hd.stel * 1e-3, sac->hd.dist);
            }

            float va, vb, rho;  
            va  = sac->hd.user6;
            vb  = sac->hd.user7;
            rho = sac->hd.user8;
            if(va <= 0.0 || vb < 0.0 || rho <= 0.0){
                GRTRaiseError("Bad src_va, src_vb or src_rho in \"%s\" header.\n", entry->d_name);
            }

            Ctrl->VpVs_ratio = (vb == 0.0)? 0.0 : ((real_t)va)/vb;

            if (Ctrl->S.mult_src_mu) {
                if(vb == 0.0){
                    GRTRaiseError("Zero src_vb in \"%s\" header. "
                        "Maybe you try to use -Su<scale> but the source is in the liquid. "
                        "Use -S<scale> instead.\n" , entry->d_name);
                }
                Ctrl->S.src_mu = vb*vb*rho*1e10;
                Ctrl->S.M0 *= Ctrl->S.src_mu;
            }
            
            grt_free_SACTRACE(sac);
            
            break;
        }

        closedir(dp);
    }
}


/** 返回目录中的任意一个 SAC 文件路径 */
static char *find_sac_in_dir(const char *dirpath)
{
    DIR *dir = opendir(dirpath);
    if(dir == NULL) return NULL;

    struct dirent *entry;
    while((entry = readdir(dir)) != NULL){
        size_t len = strlen(entry->d_name);
        if(len <= 4 || strcmp(entry->d_name + len - 4, ".sac") != 0) continue;

        char *filepath = NULL;
        GRT_SAFE_ASPRINTF(&filepath, "%s/%s", dirpath, entry->d_name);
        struct stat st;
        if(stat(filepath, &st) == 0 && S_ISREG(st.st_mode)){
            closedir(dir);
            return filepath;
        }
        GRT_SAFE_FREE_PTR(filepath);
    }

    closedir(dir);
    return NULL;
}


static bool real_value_exists(const real_t *values, size_t nvalues, real_t value, const real_t tolerance)
{
    for(size_t i = 0; i < nvalues; ++i){
        if(fabs(values[i] - value) <= tolerance) return true;
    }
    return false;
}


static void append_unique_real(real_t **values, size_t *nvalues, real_t value, const real_t tolerance)
{
    if(real_value_exists(*values, *nvalues, value, tolerance)) return;

    real_t *new_values = (real_t *)realloc(*values, sizeof(real_t) * (*nvalues + 1));
    if(new_values == NULL){
        GRTRaiseError("Failed to allocate Green's-function directory metadata.");
    }
    new_values[*nvalues] = value;
    *values = new_values;
    *nvalues += 1;
}


/** 根据根目录和库形态选取唯一的动态格林函数子目录 */
static void resolve_syn_grn_path(GRT_MODULE_CTRL *Ctrl)
{
    const real_t match_tolerance = 1e-5;
    char *sample = find_sac_in_dir(Ctrl->G.s_grnpath);
    if(sample != NULL){
        GRT_SAFE_FREE_PTR(sample);
        if(Ctrl->Depth.s_active || Ctrl->Depth.r_active || Ctrl->R.active){
            GRTRaiseError("-Ds/-Dr/-R cannot be used when -G points to a Green's-function subdirectory.");
        }
        return;
    }

    DIR *root = opendir(Ctrl->G.s_grnpath);
    if(root == NULL){
        GRTRaiseError("Unable to open Green's-function root \"%s\".", Ctrl->G.s_grnpath);
    }

    real_t *depsrcs = NULL;
    real_t *deprcvs = NULL;
    real_t *dists = NULL;
    size_t ndepsrcs = 0;
    size_t ndeprcvs = 0;
    size_t ndists = 0;
    struct dirent *entry;
    while((entry = readdir(root)) != NULL){
        if(entry->d_name[0] == '.') continue;

        char *dirpath = NULL;
        GRT_SAFE_ASPRINTF(&dirpath, "%s/%s", Ctrl->G.s_grnpath, entry->d_name);
        struct stat st;
        if(stat(dirpath, &st) != 0 || !S_ISDIR(st.st_mode)){
            GRT_SAFE_FREE_PTR(dirpath);
            continue;
        }

        sample = find_sac_in_dir(dirpath);
        if(sample == NULL){
            GRT_SAFE_FREE_PTR(dirpath);
            continue;
        }

        SACTRACE *sac = grt_read_SACTRACE(sample, true);
        append_unique_real(&depsrcs, &ndepsrcs, sac->hd.evdp, match_tolerance);
        append_unique_real(&deprcvs, &ndeprcvs, -sac->hd.stel * 1e-3, match_tolerance);
        append_unique_real(&dists, &ndists, sac->hd.dist, match_tolerance);
        grt_free_SACTRACE(sac);
        GRT_SAFE_FREE_PTR(sample);
        GRT_SAFE_FREE_PTR(dirpath);
    }
    closedir(root);

    if(ndepsrcs == 0 || ndeprcvs == 0 || ndists == 0){
        GRTRaiseError("Green's-function root contains no valid subdirectories.");
    }

    if(ndepsrcs == 1){
        if(!Ctrl->Depth.s_active){
            Ctrl->Depth.depsrc = depsrcs[0];
        }
    } else if(!Ctrl->Depth.s_active){
        GRTRaiseError("Green's-function root has multiple source depths; -Ds<depsrc> is required.");
    }

    if(ndeprcvs == 1){
        if(!Ctrl->Depth.r_active){
            Ctrl->Depth.deprcv = deprcvs[0];
        }
    } else if(!Ctrl->Depth.r_active){
        GRTRaiseError("Green's-function root has multiple receiver depths; -Dr<deprcv> is required.");
    }

    if(ndists == 1){
        if(!Ctrl->R.active){
            Ctrl->R.dist = dists[0];
        }
    } else if(!Ctrl->R.active){
        GRTRaiseError("Green's-function root has multiple epicentral distances; -R<dist> is required.");
    }

    GRT_SAFE_FREE_PTR(depsrcs);
    GRT_SAFE_FREE_PTR(deprcvs);
    GRT_SAFE_FREE_PTR(dists);

    root = opendir(Ctrl->G.s_grnpath);
    if(root == NULL){
        GRTRaiseError("Unable to open Green's-function root \"%s\".", Ctrl->G.s_grnpath);
    }

    char *match = NULL;
    size_t nmatch = 0;
    while((entry = readdir(root)) != NULL){
        if(entry->d_name[0] == '.') continue;

        char *dirpath = NULL;
        GRT_SAFE_ASPRINTF(&dirpath, "%s/%s", Ctrl->G.s_grnpath, entry->d_name);
        struct stat st;
        if(stat(dirpath, &st) != 0 || !S_ISDIR(st.st_mode)){
            GRT_SAFE_FREE_PTR(dirpath);
            continue;
        }

        sample = find_sac_in_dir(dirpath);
        if(sample == NULL){
            GRT_SAFE_FREE_PTR(dirpath);
            continue;
        }

        SACTRACE *sac = grt_read_SACTRACE(sample, true);
        bool depth_match = fabs(sac->hd.evdp - Ctrl->Depth.depsrc) <= match_tolerance
            && fabs((-sac->hd.stel * 1e-3) - Ctrl->Depth.deprcv) <= match_tolerance;
        bool dist_match = fabs(sac->hd.dist - Ctrl->R.dist) <= match_tolerance;
        grt_free_SACTRACE(sac);
        GRT_SAFE_FREE_PTR(sample);

        if(depth_match && dist_match){
            nmatch++;
            GRT_SAFE_FREE_PTR(match);
            match = dirpath;
        } else {
            GRT_SAFE_FREE_PTR(dirpath);
        }
    }
    closedir(root);

    if(nmatch == 0){
        GRTRaiseError("No Green's-function subdirectory exactly matches depsrc=%.9g, deprcv=%.9g, dist=%.9g.",
            Ctrl->Depth.depsrc, Ctrl->Depth.deprcv, Ctrl->R.dist);
    }
    if(nmatch > 1){
        GRTRaiseError("Multiple Green's-function subdirectories match the requested depths and distance.");
    }

    GRT_SAFE_FREE_PTR(Ctrl->G.s_grnpath);
    Ctrl->G.s_grnpath = match;
}


/** 将某一道合成地震图保存到sac文件 */
static void save_to_sac(GRT_MODULE_CTRL *Ctrl, const char *pfx, const char ch, SACTRACE *sac){
    sac->hd.az = Ctrl->A.azimuth;
    sac->hd.baz = Ctrl->A.backazimuth;
    char *buffer = NULL;
    snprintf(sac->hd.kcmpnm, sizeof(sac->hd.kcmpnm), "%s%c", pfx, ch);
    GRT_SAFE_ASPRINTF(&buffer, "%s/%s%c.sac", Ctrl->O.s_output_dir, pfx, ch);
    grt_write_SACTRACE(buffer, sac);
    GRT_SAFE_FREE_PTR(buffer);
}


/** 判断该震源类型是否参与合成 */
static bool syn_need_src(GRT_SYN_TYPE computeType, int im)
{
    if (computeType == GRT_SYN_EX) {
        return im == GRT_SRC_M_EX_INDEX;
    } else if (computeType == GRT_SYN_SF) {
        return im == GRT_SRC_M_VF_INDEX || im == GRT_SRC_M_HF_INDEX;
    } else if (computeType == GRT_SYN_DC) {
        return im >= GRT_SRC_M_DD_INDEX;
    } else if (computeType == GRT_SYN_TS || computeType == GRT_SYN_MT) {
        return im >= GRT_SRC_M_DD_INDEX || im == GRT_SRC_M_EX_INDEX;
    }
    return false;
}


/** 线性叠加：out += coef * gf，跳过零系数；非零系数时 gf 不可为 NULL */
static void syn_accum_from_gf(
    size_t npts, const pfloatChnlGrid gf,
    const realChnlGrid srcRadi, float *const out[GRT_CHANNEL_NUM])
{
    GRT_LOOP_ChnlGrid(im, c) {
        int modr = GRT_SRC_M_ORDERS[im];
        if(modr == 0 && GRT_ZRT_CODES[c] == 'T') continue;

        const real_t coef = srcRadi[im][c];
        if(coef == 0.0) continue;
        if(gf[im][c] == NULL){
            GRTRaiseError("Missing Green function for %s%c.\n",
                GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
        }

        float *dst = out[c];
        const float *src = gf[im][c];
        for(size_t n = 0; n < npts; ++n){
            dst[n] += (float)(src[n] * coef);
        }
    }
}


/**
 * 由动态格林函数合成三分量地震图（及可选空间偏导）
 *
 * 数组布局：gf[震源][分量][采样点]、syn[分量][采样点]、
 * syn_upar[偏导方向][分量][采样点]。gf_uiz/gf_uir 在 calc_upar=false
 * 时可传 NULL；单个分量指针为 NULL 时跳过该道
 *
 * r=0 时强制 *azrad=0（e_r→N、e_θ→E）并告警；
 * 并用 ∂_r 格林函数合成 (1/r)∂_θ 的有限部分（见函数内注释）
 *
 * @param[in,out]  azrad     方位角（弧度）；r=0 时写回 0
 */
static void syn_from_gf(
    size_t npts, float dist,
    const pfloatChnlGrid gf, const pfloatChnlGrid gf_uiz, const pfloatChnlGrid gf_uir,
    GRT_SYN_TYPE computeType, real_t M0, real_t VpVs_ratio, real_t *azrad,
    const real_t mchn[GRT_MECHANISM_NUM],
    bool rot2ZNE, bool calc_upar,
    float *const syn[GRT_CHANNEL_NUM], float *const syn_upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    // r=0：方位角无定义，约定 e_r→N、e_θ→E，故强制 az=0
    if(GRT_IS_ZERO(dist)){
        GRTRaiseWarning(
            "Zero epicentral distance: azimuth is ignored "
            "(forced to 0°; e_r→N, e_θ→E).");
        *azrad = 0.0;
    }
    const real_t az = *azrad;

    int calcUTypes = calc_upar ? 4 : 1;
    realChnlGrid srcRadi = {0};

    // 清零输出
    for(int c = 0; c < GRT_CHANNEL_NUM; ++c){
        memset(syn[c], 0, npts * sizeof(float));
        if(calc_upar){
            for(int c2 = 0; c2 < GRT_CHANNEL_NUM; ++c2){
                memset(syn_upar[c][c2], 0, npts * sizeof(float));
            }
        }
    }

    for(int ityp = 0; ityp < calcUTypes; ++ityp){
        real_t upar_scale = 1.0;
        // 求位移空间导数时，需调整比例系数（1e-5: km→cm）
        // ZRT 协变导数拆两步：此处合成 (1/r)∂_θ u，后处理再补 ±u/r。
        // r=0 时协变组合仍有限，两项分别换成有限极限（见下方 gf 与后处理）。
        if(ityp > 0){
            switch (GRT_ZRT_CODES[ityp-1]){
                case 'Z': case 'R':
                    upar_scale = 1e-5;
                    break;
                // (1/r)∂_θ：r≠0 时 scale∝1/r；r=0 时改用 ∂_r GF，scale 仅留 km→cm
                case 'T':
                    upar_scale = GRT_IS_ZERO(dist) ? 1e-5 : (1e-5 / dist);
                    break;
                default:
                    break;
            }
        }

        float *const (*up)[GRT_CHANNEL_NUM] = gf;
        if(ityp == 1){
            up = gf_uiz;
        } else if(ityp == 2){
            up = gf_uir;
        } else if(ityp == 3 && GRT_IS_ZERO(dist)){
            // r=0: 用 ∂_r GF（par_θ 辐射因子）合成 (1/r)∂_θ 的有限部分；
            // 后处理中 u/r 联络项改用 ∂_r u，二者合并得有限直角/柱坐标导数。
            // 对 u_z：m≥1 ⇒ u_z(0)=0，lim u_z/r=∂_r u_z；m=0 无 ∂_θ；无 ±u_z/r 联络项。
            up = gf_uir;
        }

        memset(srcRadi, 0, sizeof(srcRadi));
        grt_set_source_radiation(srcRadi, computeType, (ityp == 3), M0, upar_scale, VpVs_ratio, az, mchn);

        float *out_ptrs[GRT_CHANNEL_NUM];
        if(ityp == 0){
            for(int c = 0; c < GRT_CHANNEL_NUM; ++c) out_ptrs[c] = syn[c];
        } else {
            for(int c = 0; c < GRT_CHANNEL_NUM; ++c) out_ptrs[c] = syn_upar[ityp-1][c];
        }
        syn_accum_from_gf(npts, up, srcRadi, out_ptrs);
    }

    // 是否转到 ZNE
    if(rot2ZNE){
        real_t dblsyn[3] = {0};
        real_t dblupar[3][3] = {0};
        for(size_t n = 0; n < npts; ++n){
            for(int i1 = 0; i1 < GRT_CHANNEL_NUM; ++i1){
                dblsyn[i1] = syn[i1][n];
                if(calc_upar){
                    for(int i2 = 0; i2 < GRT_CHANNEL_NUM; ++i2){
                        dblupar[i1][i2] = syn_upar[i1][i2][n];
                    }
                }
            }
            if(calc_upar){
                grt_rot_zrt2zxy_upar(az, dblsyn, dblupar, dist * 1e5);  // 1e5 km→cm
            } else {
                grt_rot_zxy2zrt_vec(-az, dblsyn);
            }
            for(int i1 = 0; i1 < GRT_CHANNEL_NUM; ++i1){
                syn[i1][n] = (float)dblsyn[i1];
                if(calc_upar){
                    for(int i2 = 0; i2 < GRT_CHANNEL_NUM; ++i2){
                        syn_upar[i1][i2][n] = (float)dblupar[i1][i2];
                    }
                }
            }
        }
    }
}


/** 读取一道格林函数 SAC；不存在则返回 NULL（允许 m=0 的 T） */
static SACTRACE *syn_load_one_gf(const char *dirpath, const char *prefix, int im, int c)
{
    int modr = GRT_SRC_M_ORDERS[im];
    if(modr == 0 && GRT_ZRT_CODES[c] == 'T') return NULL;

    char *grnpath = NULL;
    GRT_SAFE_ASPRINTF(&grnpath, "%s/%s%s%c.sac", dirpath, prefix, GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
    if(access(grnpath, F_OK) != 0){
        GRT_SAFE_FREE_PTR(grnpath);
        return NULL;
    }
    SACTRACE *sac = grt_read_SACTRACE(grnpath, false);
    GRT_SAFE_FREE_PTR(grnpath);
    return sac;
}


/** 对一道合成结果做时间函数卷积 / 积分 / 微分 */
static void syn_postprocess_trace(SACTRACE *sac, SACTRACE *tfsac, int int_times, int dif_times)
{
    float dt = sac->hd.delta;
    int nt = sac->hd.npts;

    if(tfsac != NULL){
        // 卷积时间函数前先把虚频率的补偿撤回，这样似乎会更稳定
        float wI = sac->hd.user0;
        float fac = 1.0f;
        float dfac = expf(-wI * dt);
        for(int n = 0; n < nt; ++n){
            sac->data[n] *= fac;
            if(n < tfsac->hd.npts) tfsac->data[n] *= fac;
            fac *= dfac;
        }

        float *convarr = (float *)calloc(nt, sizeof(float));
        grt_oaconvolve(sac->data, nt, tfsac->data, tfsac->hd.npts, convarr, nt, true);
        fac = 1.0f;
        dfac = expf(wI * dt);
        for(int n = 0; n < nt; ++n){
            sac->data[n] = convarr[n] * fac;
            if(n < tfsac->hd.npts) tfsac->data[n] *= fac;
            fac *= dfac;
        }
        GRT_SAFE_FREE_PTR(convarr);
    }

    for(int i = 0; i < int_times; ++i){
        grt_trap_integral(sac->data, nt, dt);
    }
    for(int i = 0; i < dif_times; ++i){
        grt_differential(sac->data, nt, dt);
    }
}


/** 子模块主函数 */
int syn_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    getopt_from_command(Ctrl, argc, argv);

    bool rot2ZNE = Ctrl->N.active;
    bool calc_upar = Ctrl->e.active;
    const char *chs = rot2ZNE ? GRT_ZNE_CODES : GRT_ZRT_CODES;

    // 读入格林函数到内存
    SACTRACE *gf_sac[GRT_SRC_M_NUM][GRT_CHANNEL_NUM] = {{0}};
    SACTRACE *gf_uiz_sac[GRT_SRC_M_NUM][GRT_CHANNEL_NUM] = {{0}};
    SACTRACE *gf_uir_sac[GRT_SRC_M_NUM][GRT_CHANNEL_NUM] = {{0}};
    pfloatChnlGrid gf = {{0}};
    pfloatChnlGrid gf_uiz = {{0}};
    pfloatChnlGrid gf_uir = {{0}};

    SACTRACE *tmpl = NULL;
    GRT_LOOP_ChnlGrid(im, c) {
        if(!syn_need_src(Ctrl->computeType, im)) continue;
        int modr = GRT_SRC_M_ORDERS[im];
        if(modr == 0 && GRT_ZRT_CODES[c] == 'T') continue;

        gf_sac[im][c] = syn_load_one_gf(Ctrl->G.s_grnpath, "", im, c);
        if(gf_sac[im][c] == NULL){
            GRTRaiseError("Failed to read Green function %s%s%c.sac\n",
                "", GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
        }
        gf[im][c] = gf_sac[im][c]->data;
        if(tmpl == NULL) tmpl = gf_sac[im][c];

        if(calc_upar){
            gf_uiz_sac[im][c] = syn_load_one_gf(Ctrl->G.s_grnpath, "z", im, c);
            gf_uir_sac[im][c] = syn_load_one_gf(Ctrl->G.s_grnpath, "r", im, c);
            if(gf_uiz_sac[im][c] == NULL || gf_uir_sac[im][c] == NULL){
                GRTRaiseError("Failed to read Green function derivatives for %s%c\n",
                    GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
            }
            gf_uiz[im][c] = gf_uiz_sac[im][c]->data;
            gf_uir[im][c] = gf_uir_sac[im][c]->data;
        }
    }
    if(tmpl == NULL){
        GRTRaiseError("No Green functions loaded.\n");
    }

    int npts = tmpl->hd.npts;
    float dt = tmpl->hd.delta;

    // 分配合成结果（与格林函数同头段，数据清零）
    SACTRACE *synsac[GRT_CHANNEL_NUM] = {0};
    SACTRACE *synparsac[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM] = {{0}};
    float *syn[GRT_CHANNEL_NUM];
    float *syn_upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM] = {{0}};
    for(int c = 0; c < GRT_CHANNEL_NUM; ++c){
        synsac[c] = grt_copy_SACTRACE(tmpl, true);
        syn[c] = synsac[c]->data;
        if(calc_upar){
            for(int c2 = 0; c2 < GRT_CHANNEL_NUM; ++c2){
                synparsac[c][c2] = grt_copy_SACTRACE(tmpl, true);
                syn_upar[c][c2] = synparsac[c][c2]->data;
            }
        }
    }

    syn_from_gf(
        (size_t)npts, Ctrl->dist,
        gf, calc_upar ? gf_uiz : NULL, calc_upar ? gf_uir : NULL,
        Ctrl->computeType, Ctrl->S.M0, Ctrl->VpVs_ratio, &Ctrl->A.azrad, Ctrl->mchn,
        false, calc_upar, syn, syn_upar);

    // C 可能因 r=0 强制 azrad=0，同步方位角头段
    Ctrl->A.azimuth = Ctrl->A.azrad / DEG1;
    Ctrl->A.backazimuth = Ctrl->A.azimuth + 180.0;
    if(Ctrl->A.backazimuth >= 360.0) Ctrl->A.backazimuth -= 360.0;

    // 时间函数 / 积分 / 微分（在旋转前，与历史行为一致）
    SACTRACE *tfsac = NULL;
    if(Ctrl->D.active){
        int tfnt;
        float *tfarr = grt_get_time_function(&tfnt, dt, Ctrl->D.tftype, Ctrl->D.tfparams);
        if(tfarr == NULL){
            GRTRaiseError("get time function error.\n");
        }
        tfsac = grt_new_SACTRACE(dt, tfnt, 0.0);
        memcpy(tfsac->data, tfarr, sizeof(float) * tfnt);
        GRT_SAFE_FREE_PTR(tfarr);
    }

    for(int c = 0; c < GRT_CHANNEL_NUM; ++c){
        syn_postprocess_trace(synsac[c], tfsac, Ctrl->I.int_times, Ctrl->J.dif_times);
        if(calc_upar){
            for(int c2 = 0; c2 < GRT_CHANNEL_NUM; ++c2){
                syn_postprocess_trace(synparsac[c][c2], tfsac, Ctrl->I.int_times, Ctrl->J.dif_times);
            }
        }
    }

    // 旋转到 ZNE（在时域处理后）
    if(rot2ZNE){
        real_t dblsyn[3] = {0};
        real_t dblupar[3][3] = {0};
        for(int n = 0; n < npts; ++n){
            for(int i1 = 0; i1 < GRT_CHANNEL_NUM; ++i1){
                dblsyn[i1] = synsac[i1]->data[n];
                if(calc_upar){
                    for(int i2 = 0; i2 < GRT_CHANNEL_NUM; ++i2){
                        dblupar[i1][i2] = synparsac[i1][i2]->data[n];
                    }
                }
            }
            if(calc_upar){
                grt_rot_zrt2zxy_upar(Ctrl->A.azrad, dblsyn, dblupar, Ctrl->dist * 1e5);
            } else {
                grt_rot_zxy2zrt_vec(-Ctrl->A.azrad, dblsyn);
            }
            for(int i1 = 0; i1 < GRT_CHANNEL_NUM; ++i1){
                synsac[i1]->data[n] = (float)dblsyn[i1];
                if(calc_upar){
                    for(int i2 = 0; i2 < GRT_CHANNEL_NUM; ++i2){
                        synparsac[i1][i2]->data[n] = (float)dblupar[i1][i2];
                    }
                }
            }
        }
    }

    // 保存到 SAC
    for(int i1 = 0; i1 < GRT_CHANNEL_NUM; ++i1){
        char pfx[20] = "";
        save_to_sac(Ctrl, pfx, chs[i1], synsac[i1]);
        if(calc_upar){
            for(int i2 = 0; i2 < GRT_CHANNEL_NUM; ++i2){
                sprintf(pfx, "%c", tolower(chs[i1]));
                save_to_sac(Ctrl, pfx, chs[i2], synparsac[i1][i2]);
            }
        }
    }

    if(tfsac != NULL){
        char *buffer = NULL;
        GRT_SAFE_ASPRINTF(&buffer, "%s/sig.sac", Ctrl->O.s_output_dir);
        grt_write_SACTRACE(buffer, tfsac);
        GRT_SAFE_FREE_PTR(buffer);
    }

    if(!Ctrl->s.active){
        GRTRaiseInfo("Under \"%s\"", Ctrl->O.s_output_dir);
        GRTRaiseInfo("Synthetic Seismograms of %-13s source done.", srcTypeFullName[Ctrl->computeType]);
        if(tfsac != NULL) GRTRaiseInfo("Time Function saved.");
    }

    if(tfsac != NULL) grt_free_SACTRACE(tfsac);
    for(int i = 0; i < GRT_CHANNEL_NUM; ++i){
        grt_free_SACTRACE(synsac[i]);
        for(int j = 0; j < GRT_CHANNEL_NUM; ++j){
            if(synparsac[i][j] != NULL) grt_free_SACTRACE(synparsac[i][j]);
        }
    }
    GRT_LOOP_ChnlGrid(im, c) {
        if(gf_sac[im][c] != NULL) grt_free_SACTRACE(gf_sac[im][c]);
        if(gf_uiz_sac[im][c] != NULL) grt_free_SACTRACE(gf_uiz_sac[im][c]);
        if(gf_uir_sac[im][c] != NULL) grt_free_SACTRACE(gf_uir_sac[im][c]);
    }

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
