/**
 * @file   grt_syn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12-2
 * 
 *    根据计算好的格林函数，定义震源机制以及方位角等，生成合成的三分量地震图
 * 
 */


#include <stdio.h>
#include <unistd.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <dirent.h>
#include <ctype.h>

#include "dynamic/signals.h"
#include "common/sacio2.h"
#include "common/const.h"
#include "common/logo.h"
#include "common/colorstr.h"
#include "common/radiation.h"
#include "common/coord.h"



extern char *optarg;
extern int optind;
extern int optopt;

//****************** 在该文件以内的全局变量 ***********************//
// 命令名称
static char *command = NULL;

// 格林函数目录路径
static char *s_grnpath = NULL;
// 输出目录路径
static char *s_output_dir = NULL;
// 保存文件前缀 
static char *s_prefix = NULL;
static const char *s_prefix_default = "out";
// 方位角，以及对应弧度制
static double azimuth = 0.0, azrad = 0.0, backazimuth=0.0;
// 放大系数，对于剪切源、爆炸源、张量震源，M0是标量地震矩；对于单力源，M0是放大系数
static double M0 = 0.0;
// 在放大系数上是否需要乘上震源处的剪切模量
static bool mult_src_mu = false;
// 存储不同震源的震源机制相关参数的数组
static double mchn[MECHANISM_NUM] = {0};
// 最终要计算的震源类型
static int computeType=GRT_SYN_COMPUTE_EX;
static char s_computeType[3] = "EX";
// 和宏命令对应的震源类型全称
static const char *sourceTypeFullName[] = {"Explosion", "Single Force", "Shear", "Moment Tensor"};
// 不打印输出
static bool silenceInput=false;

// 积分次数
static int int_times = 0;
// 求导次数
static int dif_times = 0;

// 是否计算位移空间导数
static bool calc_upar=false;

// 输出分量格式，即是否需要旋转到ZNE
static bool rot2ZNE = false;

// 各选项的标志变量，初始化为0，定义了则为1
static int G_flag=0, O_flag=0, A_flag=0,
           S_flag=0, M_flag=0, F_flag=0,
           T_flag=0, P_flag=0, s_flag=0,
           D_flag=0, I_flag=0, J_flag=0, 
           e_flag=0, N_flag=0;

// 计算和位移相关量的种类（1-位移，2-ui_z，3-ui_r，4-ui_t）
static int calcUTypes=1;

// 震源名称数组，以及方向因子数组
static double srcRadi[SRC_M_NUM][CHANNEL_NUM] = {0};

// 卷积的时间函数类型
static char tftype = GRT_SIG_CUSTOM;
static char *tfparams = NULL;

// 震源处的剪切模量
static double src_mu = 0.0;




/**
 * 打印使用说明
 */
static void print_help(){
print_logo();
printf("\n"
"[grt.syn]\n\n"
"    A Supplementary Tool of GRT to Compute Three-Component \n"
"    Displacement with the outputs of command `grt`.\n"
"    Three components are:\n"
"       + Up (Z),\n"
"       + Radial Outward (R),\n"
"       + Transverse Clockwise (T),\n"
"    and the units are cm. You can add -N to rotate ZRT to ZNE.\n"
"\n"
"    + Default outputs (without -I and -J) are impulse-like displacements.\n"
"    + -D, -I and -J are applied in the frequency domain.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt.syn -G<grn_path> -A<azimuth> -S[u]<scale> -O<outdir> \n"
"            [-M<strike>/<dip>/<rake>]\n"
"            [-T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>]\n"
"            [-F<fn>/<fe>/<fz>] \n"
"            [-D<tftype>/<tfparams>] [-I<odr>] [-J<odr>]\n" 
"            [-P<prefix>] [-N] [-e] [-s]\n"
"\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<grn_path>  Green's Functions output directory of command `grt`.\n"
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
"    -P<prefix>    Prefix for single SAC file. Default is \"%s\".\n", s_prefix_default); printf(
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
"        grt -Mmilrow -N1000/0.01 -D2/0 -Ores -R2,4,6,8,10\n"
"\n"
"    Then you can get synthetic seismograms of Explosion at epicentral\n"
"    distance of 10 km and an azimuth of 30° by running:\n"
"        grt.syn -Gres/milrow_2_0_10 -Osyn_ex -A30 -S1e24\n"
"\n"
"    or Shear\n"
"        grt.syn -Gres/milrow_2_0_10 -Osyn_dc -A30 -S1e24 -M100/20/80\n"
"\n"
"    or Single Force\n"
"        grt.syn -Gres/milrow_2_0_10 -Osyn_sf -A30 -S1e24 -F0.5/-1.2/3.3\n"
"\n"
"    or Moment Tensor\n"
"        grt.syn -Gres/milrow_2_0_10 -Osyn_mt -A30 -S1e24 -T2.3/0.2/-4.0/0.3/0.5/1.2\n"
"\n\n\n"
);
}


/**
 * 检查格林函数文件是否存在
 * 
 * @param    name    格林函数文件名（不含父级目录）
 */
static void check_grn_exist(const char *name){
    char *buffer = (char*)malloc(sizeof(char)*(strlen(s_grnpath)+strlen(name)+100));
    sprintf(buffer, "%s/%s", s_grnpath, name);
    if(access(buffer, F_OK) == -1){
        fprintf(stderr, "[%s] " BOLD_RED "Error! %s not exists.\n" DEFAULT_RESTORE, command, buffer);
        exit(EXIT_FAILURE);
    }
    // 检查文件的同时将src_mu计算出来
    if(src_mu == 0.0 && mult_src_mu){
        SACHEAD hd;
        read_SAC_HEAD(command, buffer, &hd);
        double va, vb, rho;
        va = hd.user6;
        vb = hd.user7;
        rho = hd.user8;
        if(va <= 0.0 || vb < 0.0 || rho <= 0.0){
            fprintf(stderr, "[%s] " BOLD_RED "Error! Bad src_va, src_vb or src_rho in \"%s\" header.\n" DEFAULT_RESTORE, command, buffer);
            exit(EXIT_FAILURE);
        }
        if(vb == 0.0){
            fprintf(stderr, "[%s] " BOLD_RED 
                "Error! Zero src_vb in \"%s\" header. "
                "Maybe you try to use -Su<scale> but the source is in the liquid. "
                "Use -S<scale> instead.\n" 
                DEFAULT_RESTORE, command, buffer);
            exit(EXIT_FAILURE);
        }
        src_mu = vb*vb*rho*1e10;
    }
    free(buffer);
}


/**
 * 从命令行中读取选项，处理后记录到全局变量中
 * 
 * @param     argc      命令行的参数个数
 * @param     argv      多个参数字符串指针
 */
static void getopt_from_command(int argc, char **argv){
    int opt;
    while ((opt = getopt(argc, argv, ":G:A:S:M:F:T:O:P:D:I:J:Nehs")) != -1) {
        switch (opt) {
            // 格林函数路径
            case 'G':
                G_flag = 1;
                s_grnpath = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                strcpy(s_grnpath, optarg);
                // 检查是否存在该目录
                DIR *dir = opendir(s_grnpath);
                if (dir == NULL) {
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Directory \"%s\" set by -G not exists.\n" DEFAULT_RESTORE, command, s_grnpath);
                    exit(EXIT_FAILURE);
                } 
                closedir(dir);
                break;

            // 方位角
            case 'A':
                A_flag = 1;
                if(0 == sscanf(optarg, "%lf", &azimuth)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -A.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                if(azimuth < 0.0 || azimuth > 360.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Azimuth in -A must be in [0, 360].\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                backazimuth = 180.0 + azimuth;
                if(backazimuth >= 360.0) backazimuth -= 360.0;
                azrad = azimuth * DEG1;
                break;

            // 放大系数
            case 'S':
                S_flag = 1;
                {   
                    // 检查是否存在字符u，若存在表明需要乘上震源处的剪切模量
                    char *upos=NULL;
                    if((upos=strchr(optarg, 'u')) != NULL){
                        mult_src_mu = true;
                        *upos = ' ';
                    }
                }
                
                if(0 == sscanf(optarg, "%lf", &M0)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -S.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                break;
            
            // 剪切震源
            case 'M':
                M_flag = 1; 
                computeType = GRT_SYN_COMPUTE_DC;
                double strike, dip, rake;
                sprintf(s_computeType, "%s", "DC");
                if(3 != sscanf(optarg, "%lf/%lf/%lf", &strike, &dip, &rake)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -M.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                if(strike < 0.0 || strike > 360.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Strike in -M must be in [0, 360].\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                if(dip < 0.0 || dip > 90.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Dip in -M must be in [0, 90].\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                if(rake < -180.0 || rake > 180.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Rake in -M must be in [-180, 180].\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                mchn[0] = strike;
                mchn[1] = dip;
                mchn[2] = rake;
                break;

            // 单力源
            case 'F':
                F_flag = 1;
                computeType = GRT_SYN_COMPUTE_SF;
                double fn, fe, fz;
                sprintf(s_computeType, "%s", "SF");
                if(3 != sscanf(optarg, "%lf/%lf/%lf", &fn, &fe, &fz)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -F.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                mchn[0] = fn;
                mchn[1] = fe;
                mchn[2] = fz;
                break;

            // 张量震源
            case 'T':
                T_flag = 1;
                computeType = GRT_SYN_COMPUTE_MT;
                double Mxx, Mxy, Mxz, Myy, Myz, Mzz;
                sprintf(s_computeType, "%s", "MT");
                if(6 != sscanf(optarg, "%lf/%lf/%lf/%lf/%lf/%lf", &Mxx, &Mxy, &Mxz, &Myy, &Myz, &Mzz)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -T.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                mchn[0] = Mxx;
                mchn[1] = Mxy;
                mchn[2] = Mxz;
                mchn[3] = Myy;
                mchn[4] = Myz;
                mchn[5] = Mzz;
                break;

            // 输出路径
            case 'O':
                O_flag = 1;
                s_output_dir = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                strcpy(s_output_dir, optarg);
                break;

            // 保存文件前缀 
            case 'P':
                P_flag = 1; 
                s_prefix = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                strcpy(s_prefix, optarg);
                break;

            // 卷积时间函数
            case 'D':
                D_flag = 1;
                tfparams = (char*)malloc(sizeof(char)*strlen(optarg));
                if(optarg[1] != '/' || 1 != sscanf(optarg, "%c", &tftype) || 1 != sscanf(optarg+2, "%s", tfparams)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -D.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                // 检查测试
                if(! check_tftype_tfparams(tftype, tfparams)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -D.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                break;

            // 对结果做积分
            case 'I':
                I_flag = 1;
                if(1 != sscanf(optarg, "%d", &int_times)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -I.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                if(int_times <= 0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Order in -I should be positive.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                break;

            // 对结果做微分
            case 'J':
                J_flag = 1;
                if(1 != sscanf(optarg, "%d", &dif_times)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -J.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                if(dif_times <= 0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Order in -J should be positive.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                break;

            // 是否计算位移空间导数
            case 'e':
                e_flag = 1;
                calc_upar = true;
                calcUTypes = 4;
                break;

            // 是否旋转到ZNE
            case 'N':
                N_flag = 1;
                rot2ZNE = true;
                break;

            // 不打印在终端
            case 's':
                s_flag = 1;
                silenceInput = true;
                break;

            // 帮助
            case 'h':
                print_help();
                exit(EXIT_SUCCESS);
                break;

            // 参数缺失
            case ':':
                fprintf(stderr, "[%s] " BOLD_RED "Error! Option '-%c' requires an argument. Use '-h' for help.\n" DEFAULT_RESTORE, command, optopt);
                exit(EXIT_FAILURE);
                break;

            // 非法选项
            case '?':
            default:
                fprintf(stderr, "[%s] " BOLD_RED "Error! Option '-%c' is invalid. Use '-h' for help.\n" DEFAULT_RESTORE, command, optopt);
                exit(EXIT_FAILURE);
                break;
        }

    }

    // 检查必选项有没有设置
    if(argc == 1){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set options. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(G_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -G. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(A_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -A. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(S_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -S. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(O_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -O. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }

    // 只能使用一种震源
    if(M_flag + F_flag + T_flag > 1){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Only support at most one of '-M', '-F' and '-T'. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }

    // 检查对应震源的格林函数文件在不在
    if( (M_flag==0&&F_flag==0&&T_flag==0) || T_flag == 1){
        check_grn_exist("EXR.sac");
        check_grn_exist("EXZ.sac");
        if(calc_upar) {
            check_grn_exist("zEXR.sac");  check_grn_exist("rEXR.sac");
            check_grn_exist("zEXZ.sac");  check_grn_exist("rEXZ.sac");
        }
    }
    if(M_flag == 1){
        check_grn_exist("DDR.sac");
        check_grn_exist("DDZ.sac");
        check_grn_exist("DSR.sac");
        check_grn_exist("DST.sac");
        check_grn_exist("DSZ.sac");
        check_grn_exist("SSR.sac");
        check_grn_exist("SST.sac");
        check_grn_exist("SSZ.sac");
        if(calc_upar){
            check_grn_exist("zDDR.sac");  check_grn_exist("rDDR.sac");
            check_grn_exist("zDDZ.sac");  check_grn_exist("rDDZ.sac");  
            check_grn_exist("zDSR.sac");  check_grn_exist("rDSR.sac");
            check_grn_exist("zDST.sac");  check_grn_exist("rDST.sac");
            check_grn_exist("zDSZ.sac");  check_grn_exist("rDSZ.sac");  
            check_grn_exist("zSSR.sac");  check_grn_exist("rSSR.sac");  
            check_grn_exist("zSST.sac");  check_grn_exist("rSST.sac");  
            check_grn_exist("zSSZ.sac");  check_grn_exist("rSSZ.sac");  
        }
    }
    if(F_flag == 1){
        check_grn_exist("VFR.sac");
        check_grn_exist("VFZ.sac");
        check_grn_exist("HFR.sac");
        check_grn_exist("HFT.sac");
        check_grn_exist("HFZ.sac");
        if(calc_upar){
            check_grn_exist("zVFR.sac");  check_grn_exist("rVFR.sac");  
            check_grn_exist("zVFZ.sac");  check_grn_exist("rVFZ.sac");  
            check_grn_exist("zHFR.sac");  check_grn_exist("rHFR.sac");  
            check_grn_exist("zHFT.sac");  check_grn_exist("rHFT.sac");  
            check_grn_exist("zHFZ.sac");  check_grn_exist("rHFZ.sac");  
        }
    }
    

    // 建立保存目录
    if(mkdir(s_output_dir, 0777) != 0){
        if(errno != EEXIST){
            fprintf(stderr, "[%s] " BOLD_RED "Error! Unable to create folder %s. Error code: %d\n" DEFAULT_RESTORE, command, s_output_dir, errno);
            exit(EXIT_FAILURE);
        }
    }

    if(P_flag == 0){
        s_prefix = (char*)malloc(sizeof(char)*100);
        strcpy(s_prefix, s_prefix_default);
    }

    if(mult_src_mu)  M0 *= src_mu;
}


/**
 * 将某一道合成地震图保存到sac文件
 * 
 * @param      buffer      输出文件夹字符串(重复使用)
 * @param      pfx         通道名前缀
 * @param      ch          分量名， Z/R/T
 * @param      arr         数据指针
 * @param      hd          SAC头段变量
 */
static void save_to_sac(char *buffer, const char *pfx, const char ch, float *arr, SACHEAD hd){
    hd.az = azimuth;
    hd.baz = backazimuth;
    snprintf(hd.kcmpnm, sizeof(hd.kcmpnm), "%s%s%c", pfx, s_computeType, ch);
    sprintf(buffer, "%s/%s%s%c.sac", s_output_dir, pfx, s_prefix, ch);
    write_sac(buffer, hd, arr);
}

/**
 * 将时间函数保存到sac文件
 * 
 * @param      buffer      输出文件夹
 * @param      tfarr       时间函数数据指针
 * @param      tfnt        点数
 * @param      dt          采样间隔
 */
static void save_tf_to_sac(char *buffer, float *tfarr, int tfnt, float dt){
    SACHEAD hd = new_sac_head(dt, tfnt, 0.0);
    sprintf(buffer, "%s/sig.sac", s_output_dir);
    write_sac(buffer, hd, tfarr);
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
static void data_zrt2zne(float *syn[3], float *syn_upar[3][3], int nt, double azrad, double dist){
    double dblsyn[3];
    double dblupar[3][3];

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
            rot_zrt2zxy_upar(azrad, dblsyn, dblupar, dist*1e5);
        } else {
            rot_zxy2zrt_vec(-azrad, dblsyn);
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



//====================================================================================
//====================================================================================
//====================================================================================
int main(int argc, char **argv){
    command = argv[0];
    getopt_from_command(argc, argv);

    // 根据参数设置，选择分量名
    const char *chs = (rot2ZNE)? ZNEchs : ZRTchs;

    char *buffer = (char*)malloc(sizeof(char)*(strlen(s_grnpath)+strlen(s_output_dir)+strlen(s_prefix)+100));
    float **ptarrout=NULL, *arrout=NULL;
    float *arrsyn[3] = {NULL, NULL, NULL};
    float *arrsyn_upar[3][3] = {NULL};
    SACHEAD hdsyn[3], hdsyn_upar[3][3], hd0;
    SACHEAD *pthd=NULL;
    float *tfarr = NULL;
    int tfnt = 0;
    char ch;
    float coef=0.0, fac=0.0, dfac=0.0;
    float wI=0.0; // 虚频率
    int nt=0;
    float dt=0.0;
    float dist=-12345.0; // 震中距

    double upar_scale=1.0;

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
        set_source_radiation(srcRadi, computeType, (ityp==3), M0, upar_scale, azrad, mchn);

        for(int c=0; c<CHANNEL_NUM; ++c){
            ch = ZRTchs[c];
            
            // 定义SACHEAD指针
            if(ityp==0){
                pthd = &hdsyn[c];
                ptarrout = &arrsyn[c];
            } else {
                pthd = &hdsyn_upar[ityp-1][c];
                ptarrout = &arrsyn_upar[ityp-1][c];
            }
            arrout = *ptarrout;

            for(int k=0; k<SRC_M_NUM; ++k){
                coef = srcRadi[k][c];
                if(coef == 0.0) continue;

                if(ityp==0 || ityp==3){
                    sprintf(buffer, "%s/%s%c.sac", s_grnpath, SRC_M_NAME_ABBR[k], ch);
                } else {
                    sprintf(buffer, "%s/%c%s%c.sac", s_grnpath, tolower(ZRTchs[ityp-1]), SRC_M_NAME_ABBR[k], ch);
                }
                
                float *arr = read_SAC(command, buffer, pthd, NULL);
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
    
                free(arr);
            } // ENDFOR 不同震源
            
            // 再次检查内存，例如爆炸源的T分量，不会进入上述for循环，导致arrout没有分配内存
            if(arrout==NULL){
                arrout = *ptarrout = (float*)calloc(nt, sizeof(float));
                *pthd = hd0;
                continue;  // 直接跳过，认为这一分量全为0
            }
    
            if(D_flag == 1 && tfarr==NULL){
                // 获得时间函数 
                tfarr = get_time_function(&tfnt, dt, tftype, tfparams);
                if(tfarr==NULL){
                    fprintf(stderr, "[%s] " BOLD_RED "get time function error.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                fac = 1.0;
                dfac = expf(-wI*dt);
                for(int i=0; i<tfnt; ++i){
                    tfarr[i] = tfarr[i]*fac;
                    fac *= dfac;
                }
            } 

    
            // 时域循环卷积
            if(tfarr!=NULL){
                float *convarr = (float*)calloc(nt, sizeof(float));
                oaconvolve(arrout, nt, tfarr, tfnt, convarr, nt, true);
                for(int i=0; i<nt; ++i){
                    arrout[i] = convarr[i] * dt; // dt是连续卷积的系数
                }
                free(convarr);
            }
    
            // 处理虚频率
            fac = 1.0;
            dfac = expf(wI*dt);
            for(int i=0; i<nt; ++i){
                arrout[i] *= fac;
                fac *= dfac;
            }
    
            // 时域积分或求导
            for(int i=0; i<int_times; ++i){
                trap_integral(arrout, nt, dt);
            }
            for(int i=0; i<dif_times; ++i){
                differential(arrout, nt, dt);
            }
    
        } // ENDFOR 三分量
    }
    

    // 是否需要旋转
    if(rot2ZNE){
        data_zrt2zne(arrsyn, arrsyn_upar, nt, azrad, dist);
    }

    // 保存到SAC文件
    for(int i1=0; i1<CHANNEL_NUM; ++i1){
        char pfx[20]="";
        save_to_sac(buffer, pfx, chs[i1], arrsyn[i1], hdsyn[i1]);
        if(calc_upar){
            for(int i2=0; i2<CHANNEL_NUM; ++i2){
                sprintf(pfx, "%c", tolower(chs[i1]));
                save_to_sac(buffer, pfx, chs[i2], arrsyn_upar[i1][i2], hdsyn_upar[i1][i2]);
            }
        }
    }


    // 保存时间函数
    if(tfnt > 0){
        // 处理虚频率
        // 保存前恢复幅值
        fac = 1.0;
        dfac = expf(wI*dt);
        for(int i=0; i<tfnt; ++i){
            tfarr[i] *= fac;
            fac *= dfac;
        }
        save_tf_to_sac(buffer, tfarr, tfnt, dt);
    }  

    free(buffer);


    if(!silenceInput) {
        printf("[%s] Under \"%s\"\n", command, s_output_dir);
        printf("[%s] Synthetic Seismograms of %-13s source done.\n", command, sourceTypeFullName[computeType]);
        if(tfarr!=NULL) printf("[%s] Time Function saved.\n", command);
    }

    free(s_output_dir);
    free(s_prefix);
    free(s_grnpath);
    if(tfparams!=NULL) free(tfparams);
    if(tfarr!=NULL)  free(tfarr);

}

