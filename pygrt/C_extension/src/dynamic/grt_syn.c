/**
 * @file   grt_syn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12-2
 * 
 *    根据计算好的格林函数，定义震源机制以及方位角等，生成合成的三分量地震图
 * 
 */

#include <dirent.h>

#include "grt/dynamic/signals.h"
#include "grt/common/sacio2.h"
#include "grt/common/const.h"
#include "grt/common/radiation.h"
#include "grt/common/coord.h"
#include "grt/dynamic/syn.h"

#include "grt.h"

// 防止被替换为虚数单位
#undef I

// 和宏命令对应的震源类型全称
static const char *sourceTypeFullName[] = {"Explosion", "Single Force", "Shear", "Moment Tensor"};

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

    // 方向因子数组
    realChnlGrid srcRadi;

    // 最终要计算的震源类型
    int computeType;
    char s_computeType[3];

} GRT_MODULE_CTRL;


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


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
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
            
            // 剪切震源
            case 'M':
                Ctrl->M.active = true;
                Ctrl->computeType = GRT_SYN_COMPUTE_DC;
                {
                    real_t strike, dip, rake;
                    sprintf(Ctrl->s_computeType, "%s", "DC");
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &strike, &dip, &rake)){
                        GRTBadOptionError(M, "");
                    };
                    if(strike < 0.0 || strike > 360.0){
                        GRTBadOptionError(M, "Strike must be in [0, 360].");
                    }
                    if(dip < 0.0 || dip > 90.0){
                        GRTBadOptionError(M, "Dip must be in [0, 90].");
                    }
                    if(rake < -180.0 || rake > 180.0){
                        GRTBadOptionError(M, "Rake must be in [-180, 180].");
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
                Ctrl->computeType = GRT_SYN_COMPUTE_MT;
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
                Ctrl->D.active = true;
                Ctrl->D.tfparams = (char*)malloc(sizeof(char)*strlen(optarg));
                if(optarg[1] != '/' || 1 != sscanf(optarg, "%c", &Ctrl->D.tftype) || 1 != sscanf(optarg+2, "%s", Ctrl->D.tfparams)){
                    GRTBadOptionError(D, "");
                }
                // 检查测试
                if(! grt_check_tftype_tfparams(Ctrl->D.tftype, Ctrl->D.tfparams)){
                    GRTBadOptionError(D, "");
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

            if (Ctrl->S.mult_src_mu) {
                float va, vb, rho;  
                va  = sac->hd.user6;
                vb  = sac->hd.user7;
                rho = sac->hd.user8;
                if(va <= 0.0 || vb < 0.0 || rho <= 0.0){
                    GRTRaiseError("Bad src_va, src_vb or src_rho in \"%s\" header.\n", entry->d_name);
                }
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


static void data_zrt2zne(SACTRACE *synsac[3], SACTRACE *synparsac[3][3], real_t azrad)
{
    real_t dblsyn[3] = {0};
    real_t dblupar[3][3] = {0};

    bool doupar = (synparsac[0][0]!=NULL);

    float dist = synsac[0]->hd.dist;
    int nt = synsac[0]->hd.npts;

    // 对每一个时间点
    for(int n = 0; n < nt; ++n){
        // 复制数据，以调用函数
        for(int i1=0; i1<3; ++i1){
            dblsyn[i1] = synsac[i1]->data[n];
            for(int i2=0; i2<3; ++i2){
                if(doupar) dblupar[i1][i2] = synparsac[i1][i2]->data[n];
            }
        }

        if(doupar) {
            grt_rot_zrt2zxy_upar(azrad, dblsyn, dblupar, dist*1e5);   // 1e5 km 转为 cm
        } else {
            grt_rot_zxy2zrt_vec(-azrad, dblsyn);
        }

        // 将结果写入原数组
        for(int i1=0; i1<3; ++i1){
            synsac[i1]->data[n] = dblsyn[i1];
            for(int i2=0; i2<3; ++i2){
                if(doupar)  synparsac[i1][i2]->data[n] = dblupar[i1][i2];
            }
        }
    }
}



/** 子模块主函数 */
int syn_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    getopt_from_command(Ctrl, argc, argv);

    // 输出分量格式，即是否需要旋转到ZNE
    bool rot2ZNE = Ctrl->N.active;

    // 根据参数设置，选择分量名
    const char *chs = (rot2ZNE)? GRT_ZNE_CODES : GRT_ZRT_CODES;

    SACTRACE *synsac[GRT_CHANNEL_NUM] = {0};
    SACTRACE *synparsac[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM] = {0};
    SACTRACE **sacs = NULL;
    SACTRACE *tfsac = NULL;

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

            // 合成ui_t
            case 3:
                upar_scale=1e-5 / Ctrl->dist;
                break;
                
            default:
                break;
        }

        if (ityp == 0) {
            sacs = &synsac[0];
        } else {
            sacs = &synparsac[ityp-1][0];
        }

        // 重新计算方向因子
        grt_set_source_radiation(Ctrl->srcRadi, Ctrl->computeType, (ityp==3), Ctrl->S.M0, upar_scale, Ctrl->A.azrad, Ctrl->mchn);

        // 合成地震图
        if (ityp==0 || ityp==3) {
            grt_syn(Ctrl->srcRadi, Ctrl->computeType, Ctrl->G.s_grnpath, "", sacs);
        } else {
            char prefix[] = {tolower(GRT_ZRT_CODES[ityp-1]), '\0'};
            grt_syn(Ctrl->srcRadi, Ctrl->computeType, Ctrl->G.s_grnpath, prefix, sacs);
        }        

        // 首次读取获得时间函数 
        if(Ctrl->D.active && tfsac == NULL){
            int tfnt;
            float *tfarr = grt_get_time_function(&tfnt, sacs[0]->hd.delta, Ctrl->D.tftype, Ctrl->D.tfparams);
            if(tfarr == NULL){
                GRTRaiseError("get time function error.\n");
            }
            tfsac = grt_new_SACTRACE(sacs[0]->hd.delta, tfnt, 0.0);
            memcpy(tfsac->data, tfarr, sizeof(float)*tfnt);
            GRT_SAFE_FREE_PTR(tfarr);
        } 

        for(int c = 0; c < GRT_CHANNEL_NUM; ++c){
            float dt = sacs[0]->hd.delta;
            int nt = sacs[0]->hd.npts;

            // 时域循环卷积
            if(tfsac != NULL){
                float wI;
                float fac, dfac;
                // 虚频率幅值压制
                wI = sacs[0]->hd.user0;
                fac = 1.0;
                dfac = expf(- wI*dt);
                for(int n = 0; n < nt; ++n){
                    sacs[c]->data[n] *= fac;
                    if (n < tfsac->hd.npts)  tfsac->data[n] *= fac;
                    fac *= dfac;
                }

                float *convarr = (float *)calloc(nt, sizeof(float));
                grt_oaconvolve(sacs[c]->data, nt, tfsac->data, tfsac->hd.npts, convarr, nt, true);
                // 虚频率振幅恢复
                fac = 1.0;
                dfac = expf(wI*dt);
                for(int n = 0; n < nt; ++n){
                    sacs[c]->data[n] = convarr[n] * fac * dt; // dt是连续卷积的系数
                    if (n < tfsac->hd.npts)  tfsac->data[n] *= fac;
                    fac *= dfac;
                }
                GRT_SAFE_FREE_PTR(convarr);
            }

            // 时域积分或求导
            for(int i=0; i<Ctrl->I.int_times; ++i){
                grt_trap_integral(sacs[c]->data, nt, dt);
            }
            for(int i=0; i<Ctrl->J.dif_times; ++i){
                grt_differential(sacs[c]->data, nt, dt);
            }
        }
    }

    // 是否需要旋转
    if(rot2ZNE){
        data_zrt2zne(synsac, synparsac, Ctrl->A.azrad);
    }

    // 保存到SAC文件
    for(int i1=0; i1<GRT_CHANNEL_NUM; ++i1){
        char pfx[20]="";
        save_to_sac(Ctrl, pfx, chs[i1], synsac[i1]);
        if(Ctrl->e.active){
            for(int i2=0; i2<GRT_CHANNEL_NUM; ++i2){
                sprintf(pfx, "%c", tolower(chs[i1]));
                save_to_sac(Ctrl, pfx, chs[i2], synparsac[i1][i2]);
            }
        }
    }

    // 保存时间函数
    if(tfsac != NULL){
        char *buffer = NULL;
        GRT_SAFE_ASPRINTF(&buffer, "%s/sig.sac", Ctrl->O.s_output_dir);
        grt_write_SACTRACE(buffer, tfsac);
        GRT_SAFE_FREE_PTR(buffer);
    }
        
    if(! Ctrl->s.active) {
        GRTRaiseInfo("Under \"%s\"", Ctrl->O.s_output_dir);
        GRTRaiseInfo("Synthetic Seismograms of %-13s source done.", sourceTypeFullName[Ctrl->computeType]);
        if(tfsac != NULL) GRTRaiseInfo("Time Function saved.");
    }

    if(tfsac != NULL) grt_free_SACTRACE(tfsac);
    for(int i=0; i<3; ++i){
        grt_free_SACTRACE(synsac[i]);
        for(int j=0; j<3; ++j){
            grt_free_SACTRACE(synparsac[i][j]);
        }
    }

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}

