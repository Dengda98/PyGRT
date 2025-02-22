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

#include "sacio.h"
#include "const.h"
#include "logo.h"
#include "colorstr.h"
#include "signals.h"

#define DEG1 0.017453292519943295
#define COMPUTE_EXP 0
#define COMPUTE_SF 1
#define COMPUTE_DC 2
#define COMPUTE_MT 3

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
static double azimuth = 0.0, azrad = 0.0;
static double caz = 0.0, saz = 0.0;
static double caz2 = 0.0, saz2 = 0.0;
// 放大系数，对于位错源、爆炸源、张量震源，M0是标量地震矩；对于单力源，M0是放大系数
static double M0 = 0.0;
// 位错震源机制
static double strike=0.0, dip=0.0, rake=0.0;
// 单力源
static double fn=0.0, fe=0.0, fz=0.0;
// 张量震源 
static double Mxx=0.0, Mxy=0.0, Mxz=0.0, Myy=0.0, Myz=0.0, Mzz=0.0;
// 最终要计算的震源类型
static int computeType=COMPUTE_EXP;
// 和宏命令对应的震源类型全称
static const char *sourceTypeFullName[] = {"Explosion", "Single Force", "Double Couple", "Moment Tensor"};
// 不打印输出
static bool silenceInput=false;

// 积分次数
static int int_times = 0;
// 求导次数
static int dif_times = 0;

// 各选项的标志变量，初始化为0，定义了则为1
static int G_flag=0, O_flag=0, A_flag=0,
           S_flag=0, M_flag=0, F_flag=0,
           T_flag=0, P_flag=0, s_flag=0,
           D_flag=0, I_flag=0, J_flag=0;

// 三分量代号
const char chs[3] = {'Z', 'R', 'T'};

// 震源名称数组，以及方向因子数组
static const int srcnum = 6;
static const char *srcName[] = {"EX", "VF", "HF", "DD", "DS", "SS"};
static double srcCoef[3][6] = { // 三分量和chs数组对应
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
};

// 卷积的时间函数类型
static char tftype = GRT_SIG_CUSTOM;
static char *tfparams = NULL;




/**
 * 打印使用说明
 */
static void print_help(){
print_logo();
printf("\n"
"[grt.syn]\n\n"
"    A Supplementary Tool of GRT to Compute Three-Component \n"
"    Displacement with the Green's Functions from command `grt`.\n"
"    Three components are:\n"
"       + Up (Z),\n"
"       + Radial Outward (R),\n"
"       + Transverse Clockwise (T),\n"
"    and the units are cm.\n"
"\n"
"    + Default outputs (without -I and -J) are impulse-like displacements.\n"
"    + -D, -I and -J are applied in the frequency domain.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt.syn -G<grn_path> -A<azimuth> -S<scale> \n"
"            [-M<strike>/<dip>/<rake>]\n"
"            [-T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>]\n"
"            [-F<fn>/<fe>/<fz>] [-O<outdir>] \n"
"            [-D<tftype>/<tfparams>] [-I<odr>] [-J<odr>]\n" 
"            [-P<prefix>] [-s]\n"
"\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<grn_path>  Green's Functions output directory of command `grt`.\n"
"\n"
"    -A<azimuth>   Azimuth in degree, from source to station.\n"
"\n"
"    -S<scale>     Scale factor to all kinds of source. \n"
"                  + For Explosion, Double Couple and Moment Tensor,\n"
"                    unit of <scale> is dyne-cm.\n"
"                  + For Single Force, unit of <scale> is dyne.\n"
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
"    -s            Silence all outputs.\n"
"\n"
"    -h            Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    Say you have computed computed Green's functions with following \n"
"    command:\n"
"        grt -Mmilrow -N1000/0.01 -D2/0 -Ores -R2,4,6,8,10\n"
"\n"
"    Then you can get synthetic seismograms of Explosion at epicentral\n"
"    distance of 10 km and an azimuth of 30° by running:\n"
"        grt.syn -Gres/milrow_2_0_10 -Osyn_ex -A30 -S1e24\n"
"\n"
"    or Double Couple\n"
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
    while ((opt = getopt(argc, argv, ":G:A:S:M:F:T:O:P:D:I:J:hs")) != -1) {
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
                azrad = azimuth * DEG1;
                saz = sin(azrad);
                caz = cos(azrad);
                saz2 = 2.0*saz*caz;
                caz2 = 2.0*caz*caz - 1.0;
                break;

            // 放大系数
            case 'S':
                S_flag = 1;
                if(0 == sscanf(optarg, "%lf", &M0)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -S.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                break;
            
            // 位错震源
            case 'M':
                M_flag = 1; 
                computeType = COMPUTE_DC;
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
                break;

            // 单力源
            case 'F':
                F_flag = 1;
                computeType = COMPUTE_SF;
                if(3 != sscanf(optarg, "%lf/%lf/%lf", &fn, &fe, &fz)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -F.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                break;

            // 张量震源
            case 'T':
                T_flag = 1;
                computeType = COMPUTE_MT;
                if(6 != sscanf(optarg, "%lf/%lf/%lf/%lf/%lf/%lf", &Mxx, &Mxy, &Mxz, &Myy, &Myz, &Mzz)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -T.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
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

    // 只能使用一种震源
    if(M_flag + F_flag + T_flag > 1){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Only support at most one of '-M', '-F' and '-T'. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }

    // 检查对应震源的格林函数文件在不在
    if( (M_flag==0&&F_flag==0&&T_flag==0) || T_flag == 1){
        check_grn_exist("EXR.sac");
        check_grn_exist("EXZ.sac");
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
    }
    if(F_flag == 1){
        check_grn_exist("VFR.sac");
        check_grn_exist("VFZ.sac");
        check_grn_exist("HFR.sac");
        check_grn_exist("HFT.sac");
        check_grn_exist("HFZ.sac");
    }


    if(O_flag == 1){
        // 建立保存目录
        if(mkdir(s_output_dir, 0777) != 0){
            if(errno != EEXIST){
                fprintf(stderr, "[%s] " BOLD_RED "Error! Unable to create folder %s. Error code: %d\n" DEFAULT_RESTORE, command, s_output_dir, errno);
                exit(EXIT_FAILURE);
            }
        }
    } else {
        // 使用当前目录
        s_output_dir = (char*)malloc(sizeof(char)*100);
        strcpy(s_output_dir, ".");
    }

    if(P_flag == 0){
        s_prefix = (char*)malloc(sizeof(char)*100);
        strcpy(s_prefix, s_prefix_default);
    }
}


/**
 * 将某一道合成地震图保存到sac文件
 * 
 * @param      buffer      输出文件夹字符串(重复使用)
 * @param      ch          分量名， Z/R/T
 * @param      arr         数据指针
 * @param      hd          SAC头段变量
 */
static void save_to_sac(char *buffer, const char ch, float *arr, SACHEAD hd){
    hd.az = azimuth;
    sprintf(hd.kcmpnm, "HH%c", ch);
    sprintf(buffer, "%s/%s%c.sac", s_output_dir, s_prefix, ch);
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
 * 设置每个震源的方向因子
 */
static void set_source_coef(){
    double mult;
    if(computeType == COMPUTE_SF){
        mult = 1e-15*M0;
    } else {
        mult = 1e-20*M0;
    }

    if(computeType == COMPUTE_EXP){
        srcCoef[0][0] = srcCoef[1][0] = mult; // Z/R
        srcCoef[2][0] = 0.0; // T
    }  
    else if(computeType == COMPUTE_SF){
        // 公式(4.6.20)
        srcCoef[0][1] = srcCoef[1][1] = fz*mult; // VF, Z/R
        srcCoef[0][2] = srcCoef[1][2] = (fn*caz + fe*saz)*mult; // HF, Z/R
        srcCoef[2][1] = 0.0; // VF, T
        srcCoef[2][2] = (- fn*saz + fe*caz)*mult; // HF, T
    }
    else if(computeType == COMPUTE_DC){
        // 公式(4.8.35)
        double stkrad = strike*DEG1;
        double diprad = dip*DEG1;
        double rakrad = rake*DEG1;
        double therad = azrad - stkrad;
        double srak, crak, sdip, cdip, sdip2, cdip2, sthe, cthe, sthe2, cthe2;
        srak = sin(rakrad);     crak = cos(rakrad);
        sdip = sin(diprad);     cdip = cos(diprad);
        sdip2 = 2.0*sdip*cdip;  cdip2 = 2.0*cdip*cdip - 1.0;
        sthe = sin(therad);     cthe = cos(therad);
        sthe2 = 2.0*sthe*cthe;  cthe2 = 2.0*cthe*cthe - 1.0;

        double A0, A1, A2, A4, A5;
        A0 = mult * (0.5*sdip2*srak);
        A1 = mult * (cdip*crak*cthe - cdip2*srak*sthe);
        A2 = mult * (0.5*sdip2*srak*cthe2 + sdip*crak*sthe2);
        A4 = mult * (- cdip2*srak*cthe - cdip*crak*sthe);
        A5 = mult * (sdip*crak*cthe2 - 0.5*sdip2*srak*sthe2);

        srcCoef[0][3] = srcCoef[1][3] = A0; // DD, Z/R
        srcCoef[0][4] = srcCoef[1][4] = A1; // DS, Z/R
        srcCoef[0][5] = srcCoef[1][5] = A2; // SS, Z/R
        srcCoef[2][3] = 0.0; // DD, T
        srcCoef[2][4] = A4;  // DS, T
        srcCoef[2][5] = A5;  // DS, T
    }
    else if(computeType == COMPUTE_MT){
        // 公式(4.9.7)但修改了各向同性的量
        double M11, M12, M13, M22, M23, M33;
        M11 = Mxx;   M12 = Mxy;   M13 = Mxz;
                     M22 = Myy;   M23 = Myz;
                                  M33 = Mzz;
        double Mexp = (M11 + M22 + M33)/3.0;
        M11 -= Mexp;
        M22 -= Mexp;
        M33 -= Mexp;

        double A0, A1, A2, A4, A5;
        A0 = mult * ((2.0*M33 - M11 - M22)/6.0 );
        A1 = mult * (- (M13*caz + M23*saz));
        A2 = mult * (0.5*(M11 - M22)*caz2+ M12*saz2);
        A4 = mult * (M13*saz - M23*caz);
        A5 = mult * (-0.5*(M11 - M22)*saz2 + M12*caz2);

        srcCoef[0][0] = srcCoef[1][0] = mult*Mexp; // EX, Z/R
        srcCoef[0][3] = srcCoef[1][3] = A0; // DD, Z/R
        srcCoef[0][4] = srcCoef[1][4] = A1; // DS, Z/R
        srcCoef[0][5] = srcCoef[1][5] = A2; // SS, Z/R
        srcCoef[2][0] = 0.0; // EX, T
        srcCoef[2][3] = 0.0; // DD, T
        srcCoef[2][4] = A4;  // DS, T
        srcCoef[2][5] = A5;  // DS, T
    }
}



int main(int argc, char **argv){
    command = argv[0];
    getopt_from_command(argc, argv);

    set_source_coef();

    char *buffer = (char*)malloc(sizeof(char)*(strlen(s_grnpath)+strlen(s_output_dir)+strlen(s_prefix)+100));
    float *arr, *arrout=NULL, *convarr=NULL;
    SACHEAD hd;
    float *tfarr = NULL;
    int tfnt = 0;
    char ch;
    float coef=0.0, fac=0.0, dfac=0.0;
    float wI=0.0; // 虚频率
    int nt=0;
    float dt=0.0;


    for(int c=0; c<3; ++c){
        ch = chs[c];
        for(int k=0; k<srcnum; ++k){
            coef = srcCoef[c][k];
            if(coef == 0.0) continue;

            sprintf(buffer, "%s/%s%c.sac", s_grnpath, srcName[k], ch);
            arr = read_sac(buffer, &hd);
            nt = hd.npts;
            dt = hd.delta;
            // dw = PI2/(nt*dt);
            if(arrout==NULL){
                arrout = (float*)calloc(nt, sizeof(float));
            }    

            // 使用虚频率将序列压制，卷积才会稳定
            // 读入虚频率 
            wI = hd.user0;
            fac = 1.0;
            dfac = expf(-wI*dt);
            for(int n=0; n<nt; ++n){
                arrout[n] += arr[n]*coef * fac;
                fac *= dfac;
            }

            free(arr);
        } // ENDFOR 不同震源

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
            convarr = (float*)calloc(nt, sizeof(float));
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

        // 保存成SAC文件
        save_to_sac(buffer, ch, arrout, hd);

        // 置零
        for(int n=0; n<nt; ++n){
            arrout[n] = 0.0f;
        }
    } // ENDFOR 三分量


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
        save_tf_to_sac(buffer, tfarr, tfnt, hd.delta);
    }  

    free(arrout);
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

