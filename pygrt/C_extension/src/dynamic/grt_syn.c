/**
 * @file   grt_syn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12-2
 * 
 *    根据计算好的格林函数，定义震源机制以及方位角等，生成合成的三分量地震图
 * 
 */

#include "grt/dynamic/signals.h"
#include "grt/common/sacio2.h"
#include "grt/common/const.h"
#include "grt/common/radiation.h"
#include "grt/common/coord.h"

#include "grt.h"

// 防止被替换为虚数单位
#undef I

// 和宏命令对应的震源类型全称
static const char *sourceTypeFullName[] = {"Explosion", "Single Force", "Shear", "Moment Tensor"};

/** 该子模块的参数控制结构体 */
typedef struct {
    char *name;
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
        int tfnt;
        float *tfarr;
    } D;
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

    // 方向因子数组
    realChnlGrid srcRadi;

    // 最终要计算的震源类型
    int computeType;
    char s_computeType[3];

} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl->name);
    // G
    GRT_SAFE_FREE_PTR(Ctrl->G.s_grnpath);
    // O
    GRT_SAFE_FREE_PTR(Ctrl->O.s_output_dir);
    // D
    GRT_SAFE_FREE_PTR(Ctrl->D.tfparams);
    GRT_SAFE_FREE_PTR(Ctrl->D.tfarr);
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
"            [-M<strike>/<dip>/<rake>]\n"
"            [-T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>]\n"
"            [-F<fn>/<fe>/<fz>] \n"
"            [-D<tftype>/<tfparams>] [-I<odr>] [-J<odr>]\n" 
"            [-N] [-e] [-s]\n"
"\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<grn_path>  Green's Functions output directory of module `greenfn`.\n"
"\n"
"    -A<azimuth>   Azimuth in degree, from source to station.\n"
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
"    -O<outdir>    Directory of output for saving. Default is\n"
"                  current directory.\n"
"\n"
"    -D<tftype>/<tfparams>\n"
"                  Convolve a Time Function with a maximum value of 1.0.\n"
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
"                    of the amplitude. File header can write unlimited lines\n"
"                    of comments with prefix \"#\".\n"
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
"    or Single Force\n"
"        grt syn -Gres/milrow_2_0_10 -Osyn_sf -A30 -S1e24 -F0.5/-1.2/3.3\n"
"\n"
"    or Moment Tensor\n"
"        grt syn -Gres/milrow_2_0_10 -Osyn_mt -A30 -S1e24 -T2.3/0.2/-4.0/0.3/0.5/1.2\n"
"\n\n\n"
);
}


/**
 * 检查格林函数文件是否存在
 * 
 * @param    name    格林函数文件名（不含父级目录）
 */
static void check_grn_exist(GRT_MODULE_CTRL *Ctrl, const char *name){
    const char *command = Ctrl->name;
    char *buffer = NULL;
    GRT_SAFE_ASPRINTF(&buffer, "%s/%s", Ctrl->G.s_grnpath, name);
    GRTCheckFileExist(command, buffer);

    // 检查文件的同时将src_mu计算出来
    if(Ctrl->S.src_mu == 0.0 && Ctrl->S.mult_src_mu){
        SACHEAD hd;
        grt_read_SAC_HEAD(command, buffer, &hd);
        real_t va, vb, rho;
        va = hd.user6;
        vb = hd.user7;
        rho = hd.user8;
        if(va <= 0.0 || vb < 0.0 || rho <= 0.0){
            GRTRaiseError("[%s] Error! Bad src_va, src_vb or src_rho in \"%s\" header.\n", command, buffer);
        }
        if(vb == 0.0){
            GRTRaiseError("[%s] Error! Zero src_vb in \"%s\" header. "
                "Maybe you try to use -Su<scale> but the source is in the liquid. "
                "Use -S<scale> instead.\n" , command, buffer);
        }
        Ctrl->S.src_mu = vb*vb*rho*1e10;
    }
    GRT_SAFE_FREE_PTR(buffer);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    const char *command = Ctrl->name;

    // 先为个别参数设置非0初始值
    Ctrl->computeType = GRT_SYN_COMPUTE_EX;
    sprintf(Ctrl->s_computeType, "%s", "EX");

    int opt;
    while ((opt = getopt(argc, argv, ":G:A:S:M:F:T:O:D:I:J:Nehs")) != -1) {
        switch (opt) {
            // 格林函数路径
            case 'G':
                Ctrl->G.active = true;
                Ctrl->G.s_grnpath = strdup(optarg);
                // 检查是否存在该目录
                GRTCheckDirExist(command, Ctrl->G.s_grnpath);
                break;

            // 方位角
            case 'A':
                Ctrl->A.active = true;
                if(0 == sscanf(optarg, "%lf", &Ctrl->A.azimuth)){
                    GRTBadOptionError(command, A, "");
                };
                if(Ctrl->A.azimuth < 0.0 || Ctrl->A.azimuth > 360.0){
                    GRTBadOptionError(command, A, "Azimuth must be in [0, 360].");
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
                    GRTBadOptionError(command, S, "");
                };
                break;
            
            // 剪切震源
            case 'M':
                Ctrl->M.active = true;
                Ctrl->computeType = GRT_SYN_COMPUTE_DC;
                {
                    real_t strike, dip, rake;
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
                    real_t fn, fe, fz;
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
                    real_t Mxx, Mxy, Mxz, Myy, Myz, Mzz;
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

            // 输出路径
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.s_output_dir = strdup(optarg);
                break;

            // 卷积时间函数
            case 'D':
                Ctrl->D.active = true;
                Ctrl->D.tfparams = (char*)malloc(sizeof(char)*strlen(optarg));
                if(optarg[1] != '/' || 1 != sscanf(optarg, "%c", &Ctrl->D.tftype) || 1 != sscanf(optarg+2, "%s", Ctrl->D.tfparams)){
                    GRTBadOptionError(command, D, "");
                }
                // 检查测试
                if(! grt_check_tftype_tfparams(Ctrl->D.tftype, Ctrl->D.tfparams)){
                    GRTBadOptionError(command, D, "");
                }
                break;

            // 对结果做积分
            case 'I':
                Ctrl->I.active = true;
                if(1 != sscanf(optarg, "%d", &Ctrl->I.int_times)){
                    GRTBadOptionError(command, I, "");
                }
                if(Ctrl->I.int_times <= 0){
                    GRTBadOptionError(command, I, "Order should be positive.");
                }
                break;

            // 对结果做微分
            case 'J':
                Ctrl->J.active = true;
                if(1 != sscanf(optarg, "%d", &Ctrl->J.dif_times)){
                    GRTBadOptionError(command, J, "");
                }
                if(Ctrl->J.dif_times <= 0){
                    GRTBadOptionError(command, J, "Order should be positive.");
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

            GRT_Common_Options_in_Switch(command, (char)(optopt));
        }

    }

    // 检查必选项有没有设置
    GRTCheckOptionSet(command, argc > 1);
    GRTCheckOptionActive(command, Ctrl, G);
    GRTCheckOptionActive(command, Ctrl, A);
    GRTCheckOptionActive(command, Ctrl, S);
    GRTCheckOptionActive(command, Ctrl, O);

    // 只能使用一种震源
    if(Ctrl->M.active + Ctrl->F.active + Ctrl->T.active > 1){
        GRTRaiseError("[%s] Error! Only support at most one of \"-M\", \"-F\" and \"-T\". Use \"-h\" for help.\n", command);
    }

    #define EX_GRN_List \
        X(EXR.sac)   \
        X(EXZ.sac)   \

    #define SF_GRN_List \
        X(VFR.sac)   \
        X(VFZ.sac)   \
        X(HFR.sac)   \
        X(HFT.sac)   \
        X(HFZ.sac)   \

    #define DC_GRN_List \
        X(DDR.sac)   \
        X(DDZ.sac)   \
        X(DSR.sac)   \
        X(DST.sac)   \
        X(DSZ.sac)   \
        X(SSR.sac)   \
        X(SST.sac)   \
        X(SSZ.sac)   \

    // 检查对应震源的格林函数文件在不在
    if( (!Ctrl->M.active && !Ctrl->F.active && !Ctrl->T.active) || Ctrl->T.active){
        #define X(name)  check_grn_exist(Ctrl, #name);
            EX_GRN_List
        #undef X
        if(Ctrl->e.active) {
            #define X(name)  check_grn_exist(Ctrl, "z" #name);
                EX_GRN_List
            #undef X
            #define X(name)  check_grn_exist(Ctrl, "r" #name);
                EX_GRN_List
            #undef X
        }
    }
    if(Ctrl->M.active){
        #define X(name)  check_grn_exist(Ctrl, #name);
            DC_GRN_List
        #undef X
        if(Ctrl->e.active) {
            #define X(name)  check_grn_exist(Ctrl, "z" #name);
                DC_GRN_List
            #undef X
            #define X(name)  check_grn_exist(Ctrl, "r" #name);
                DC_GRN_List
            #undef X
        }
    }
    if(Ctrl->F.active){
        #define X(name)  check_grn_exist(Ctrl, #name);
            SF_GRN_List
        #undef X
        if(Ctrl->e.active) {
            #define X(name)  check_grn_exist(Ctrl, "z" #name);
                SF_GRN_List
            #undef X
            #define X(name)  check_grn_exist(Ctrl, "r" #name);
                SF_GRN_List
            #undef X
        }
    }

    // 建立保存目录
    GRTCheckMakeDir(command, Ctrl->O.s_output_dir);

    if(Ctrl->S.mult_src_mu)  Ctrl->S.M0 *= Ctrl->S.src_mu;
}


/**
 * 将某一道合成地震图保存到sac文件
 * 
 * @param      pfx         通道名前缀
 * @param      ch          分量名， Z/R/T
 * @param      arr         数据指针
 * @param      hd          SAC头段变量
 */
static void save_to_sac(GRT_MODULE_CTRL *Ctrl, const char *pfx, const char ch, float *arr, SACHEAD hd){
    hd.az = Ctrl->A.azimuth;
    hd.baz = Ctrl->A.backazimuth;
    char *buffer = NULL;
    snprintf(hd.kcmpnm, sizeof(hd.kcmpnm), "%s%c", pfx, ch);
    GRT_SAFE_ASPRINTF(&buffer, "%s/%s%c.sac", Ctrl->O.s_output_dir, pfx, ch);
    write_sac(buffer, hd, arr);
    GRT_SAFE_FREE_PTR(buffer);
}

/**
 * 将时间函数保存到sac文件
 * 
 * @param      tfarr       时间函数数据指针
 * @param      tfnt        点数
 * @param      dt          采样间隔
 */
static void save_tf_to_sac(GRT_MODULE_CTRL *Ctrl, float *tfarr, int tfnt, float dt){
    char *buffer = NULL;
    SACHEAD hd = new_sac_head(dt, tfnt, 0.0);
    GRT_SAFE_ASPRINTF(&buffer, "%s/sig.sac", Ctrl->O.s_output_dir);
    write_sac(buffer, hd, tfarr);
    GRT_SAFE_FREE_PTR(buffer);
}


/**
 * 将不同ZRT分量的位移以及位移空间导数旋转到ZNE分量
 * 
 * @param    syn       位移
 * @param    syn_upar  位移空间导数
 * @param    nt        时间点数
 * @param    azrad     方位角弧度
 * @param    dist      震中距(km)
 */
static void data_zrt2zne(float *syn[3], float *syn_upar[3][3], int nt, real_t azrad, real_t dist){
    real_t dblsyn[3];
    real_t dblupar[3][3];

    bool doupar = (syn_upar[0][0]!=NULL);

    // 对每一个时间点
    for(int n=0; n<nt; ++n){
        // 复制数据，以调用函数
        for(int i1=0; i1<3; ++i1){
            dblsyn[i1] = syn[i1][n];
            for(int i2=0; i2<3; ++i2){
                if(doupar) dblupar[i1][i2] = syn_upar[i1][i2][n];
            }
        }

        if(doupar) {
            grt_rot_zrt2zxy_upar(azrad, dblsyn, dblupar, dist*1e5);
        } else {
            grt_rot_zxy2zrt_vec(-azrad, dblsyn);
        }
        

        // 将结果写入原数组
        for(int i1=0; i1<3; ++i1){
            syn[i1][n] = dblsyn[i1];
            for(int i2=0; i2<3; ++i2){
                if(doupar)  syn_upar[i1][i2][n] = dblupar[i1][i2];
            }
        }
    }
}



/** 子模块主函数 */
int syn_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    Ctrl->name = strdup(argv[0]);
    const char *command = Ctrl->name;

    getopt_from_command(Ctrl, argc, argv);

    // 输出分量格式，即是否需要旋转到ZNE
    bool rot2ZNE = Ctrl->N.active;

    // 根据参数设置，选择分量名
    const char *chs = (rot2ZNE)? GRT_ZNE_CODES : GRT_ZRT_CODES;

    float **ptarrout=NULL, *arrout=NULL;
    float *arrsyn[3] = {NULL, NULL, NULL};
    float *arrsyn_upar[3][3] = {NULL};
    SACHEAD hdsyn[3], hdsyn_upar[3][3], hd0;
    SACHEAD *pthd=NULL;
    char ch;
    float coef=0.0, fac=0.0, dfac=0.0;
    float wI=0.0; // 虚频率
    int nt=0;
    float dt=0.0;
    float dist=-12345.0; // 震中距

    real_t upar_scale=1.0;

    // 计算和位移相关量的种类（1-位移，2-ui_z，3-ui_r，4-ui_t）
    int calcUTypes = (Ctrl->e.active)? 4 : 1;

    for(int ityp=0; ityp<calcUTypes; ++ityp){
        // 求位移空间导数时，需调整比例系数
        switch (ityp){
            // 合成位移
            case 0:
                upar_scale=1.0;
                break;

            // 合成ui_z
            case 1:
            // 合成ui_r
            case 2:
                upar_scale=1e-5;
                break;

            // 合成ui_t，其中dist会在ityp<3之前从sac文件中读出
            case 3:
                upar_scale=1e-5 / dist;
                break;
                
            default:
                break;
        }
        
        // 重新计算方向因子
        grt_set_source_radiation(Ctrl->srcRadi, Ctrl->computeType, (ityp==3), Ctrl->S.M0, upar_scale, Ctrl->A.azrad, Ctrl->mchn);

        for(int c=0; c<GRT_CHANNEL_NUM; ++c){
            ch = GRT_ZRT_CODES[c];
            
            // 定义SACHEAD指针
            if(ityp==0){
                pthd = &hdsyn[c];
                ptarrout = &arrsyn[c];
            } else {
                pthd = &hdsyn_upar[ityp-1][c];
                ptarrout = &arrsyn_upar[ityp-1][c];
            }
            arrout = *ptarrout;

            for(int k=0; k<GRT_SRC_M_NUM; ++k){
                coef = Ctrl->srcRadi[k][c];
                if(coef == 0.0) continue;

                char *buffer = NULL;
                if(ityp==0 || ityp==3){
                    GRT_SAFE_ASPRINTF(&buffer, "%s/%s%c.sac", Ctrl->G.s_grnpath, GRT_SRC_M_NAME_ABBR[k], ch);
                } else {
                    GRT_SAFE_ASPRINTF(&buffer, "%s/%c%s%c.sac", Ctrl->G.s_grnpath, tolower(GRT_ZRT_CODES[ityp-1]), GRT_SRC_M_NAME_ABBR[k], ch);
                }
                
                float *arr = grt_read_SAC(command, buffer, pthd, NULL);
                hd0 = *pthd; // 备份一份

                nt = pthd->npts;
                dt = pthd->delta;
                dist = pthd->dist;
                // dw = PI2/(nt*dt);

                // 第一次读入元信息，申请数组内存
                if(arrout==NULL){
                    arrout = *ptarrout = (float*)calloc(nt, sizeof(float));
                }    
    
                // 使用虚频率将序列压制，卷积才会稳定
                // 读入虚频率 
                wI = pthd->user0;
                fac = 1.0;
                dfac = expf(-wI*dt);
                for(int n=0; n<nt; ++n){
                    arrout[n] += arr[n]*coef * fac;
                    fac *= dfac;
                }
    
                GRT_SAFE_FREE_PTR(arr);
                GRT_SAFE_FREE_PTR(buffer);
            } // ENDFOR 不同震源
            
            // 再次检查内存，例如爆炸源的T分量，不会进入上述for循环，导致arrout没有分配内存
            if(arrout==NULL){
                arrout = *ptarrout = (float*)calloc(nt, sizeof(float));
                *pthd = hd0;
                continue;  // 直接跳过，认为这一分量全为0
            }
    
            if(Ctrl->D.active && Ctrl->D.tfarr==NULL){
                // 获得时间函数 
                Ctrl->D.tfarr = grt_get_time_function(&Ctrl->D.tfnt, dt, Ctrl->D.tftype, Ctrl->D.tfparams);
                if(Ctrl->D.tfarr==NULL){
                    GRTRaiseError("[%s] get time function error.\n", command);
                }
                fac = 1.0;
                dfac = expf(-wI*dt);
                for(int i=0; i<Ctrl->D.tfnt; ++i){
                    Ctrl->D.tfarr[i] = Ctrl->D.tfarr[i]*fac;
                    fac *= dfac;
                }
            } 

    
            // 时域循环卷积
            if(Ctrl->D.tfarr!=NULL){
                float *convarr = (float*)calloc(nt, sizeof(float));
                grt_oaconvolve(arrout, nt, Ctrl->D.tfarr, Ctrl->D.tfnt, convarr, nt, true);
                for(int i=0; i<nt; ++i){
                    arrout[i] = convarr[i] * dt; // dt是连续卷积的系数
                }
                GRT_SAFE_FREE_PTR(convarr);
            }
    
            // 处理虚频率
            fac = 1.0;
            dfac = expf(wI*dt);
            for(int i=0; i<nt; ++i){
                arrout[i] *= fac;
                fac *= dfac;
            }
    
            // 时域积分或求导
            for(int i=0; i<Ctrl->I.int_times; ++i){
                grt_trap_integral(arrout, nt, dt);
            }
            for(int i=0; i<Ctrl->J.dif_times; ++i){
                grt_differential(arrout, nt, dt);
            }
    
        } // ENDFOR 三分量
    }
    

    // 是否需要旋转
    if(rot2ZNE){
        data_zrt2zne(arrsyn, arrsyn_upar, nt, Ctrl->A.azrad, dist);
    }

    // 保存到SAC文件
    for(int i1=0; i1<GRT_CHANNEL_NUM; ++i1){
        char pfx[20]="";
        save_to_sac(Ctrl, pfx, chs[i1], arrsyn[i1], hdsyn[i1]);
        if(Ctrl->e.active){
            for(int i2=0; i2<GRT_CHANNEL_NUM; ++i2){
                sprintf(pfx, "%c", tolower(chs[i1]));
                save_to_sac(Ctrl, pfx, chs[i2], arrsyn_upar[i1][i2], hdsyn_upar[i1][i2]);
            }
        }
    }

    // 保存时间函数
    if(Ctrl->D.tfnt > 0){
        // 处理虚频率
        // 保存前恢复幅值
        fac = 1.0;
        dfac = expf(wI*dt);
        for(int i=0; i<Ctrl->D.tfnt; ++i){
            Ctrl->D.tfarr[i] *= fac;
            fac *= dfac;
        }
        save_tf_to_sac(Ctrl, Ctrl->D.tfarr, Ctrl->D.tfnt, dt);
    }  
    
    if(! Ctrl->s.active) {
        printf("[%s] Under \"%s\"\n", command, Ctrl->O.s_output_dir);
        printf("[%s] Synthetic Seismograms of %-13s source done.\n", command, sourceTypeFullName[Ctrl->computeType]);
        if(Ctrl->D.tfarr!=NULL) printf("[%s] Time Function saved.\n", command);
    }
    
    
    for(int i=0; i<3; ++i){
        GRT_SAFE_FREE_PTR(arrsyn[i]);
        for(int j=0; j<3; ++j){
            GRT_SAFE_FREE_PTR(arrsyn_upar[i][j]);
        }
    }

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}

