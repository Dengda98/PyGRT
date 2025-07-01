/**
 * @file   grt_main.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-11-28
 * 
 *    定义main函数，形成命令行式的用法（不使用python的entry_points，会牺牲性能）
 *    计算不同震源的格林函数
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
#include <omp.h>

#include "dynamic/grt.h"
#include "dynamic/signals.h"
#include "travt/travt.h"
#include "common/const.h"
#include "common/model.h"
#include "common/search.h"
#include "common/sacio.h"
#include "common/logo.h"
#include "common/colorstr.h"


#ifdef GRT_USE_FLOAT 
#define _FFTW_COMPLEX   fftwf_complex
#define _FFTW_PLAN      fftwf_plan
#define _FFTW_EXECUTE   fftwf_execute
#define _FFTW_MALLOC    fftwf_malloc
#define _FFTW_FREE      fftwf_free
#define _FFTW_DESTROY_PLAN   fftwf_destroy_plan
#define _FFTW_PLAN_DFT_C2R_1D   fftwf_plan_dft_c2r_1d
#else 
#define _FFTW_COMPLEX   fftw_complex
#define _FFTW_PLAN      fftw_plan
#define _FFTW_EXECUTE   fftw_execute
#define _FFTW_MALLOC    fftw_malloc
#define _FFTW_FREE      fftw_free
#define _FFTW_DESTROY_PLAN   fftw_destroy_plan
#define _FFTW_PLAN_DFT_C2R_1D   fftw_plan_dft_c2r_1d
#endif


extern char *optarg;
extern int optind;
extern int optopt;

//****************** 在该文件以内的全局变量 ***********************//
// 命令名称
static char *command = NULL;
// 线程数
static int nthreads=-1;
// 模型路径，模型PYMODEL1D指针，全局最大最小速度
static char *s_modelpath = NULL;
static char *s_modelname = NULL;
static PYMODEL1D *pymod;
static double vmax, vmin;
// 震源和场点深度
static double depsrc, deprcv;
static char *s_depsrc = NULL, *s_deprcv = NULL;
// 点数、采样间隔、时窗长度、频率点数、频率间隔
static int nt; 
static double dt;
static double winT;
static int nf1, nf2, nf;
static double df, dw;
static MYREAL *freqs = NULL;
// 输出目录
static char *s_output_dir = NULL;
// 计算频率范围
static double freq1=-1.0, freq2=-1.0;
// 虚频率系数和虚频率
static double zeta=0.8, wI=0.0;
// 波数积分间隔, Filon积分间隔，自适应Filon积分采样精度，Filon积分起始点
static double Length=0.0, filonLength=0.0, safilonTol=0.0, filonCut=0.0;
// 波数积分相关变量
static double keps=-1.0, ampk=1.15, k0=5.0;
// 参考最小速度，小于0表示使用峰谷平均法;
static double vmin_ref=0.0;
static const double min_vmin_ref=0.1;
// 自动使用峰谷平均法的最小厚度差
static const double hs_ptam = MIN_DEPTH_GAP_SRC_RCV;
// 时间延迟量，延迟参考速度。总延迟=T0 + dist/V0;
static double delayT=0.0, delayT0=0.0, delayV0=0.0;
static double tmax; // 时窗最大截止时刻
// 震中距数组以及保存对应初至波走时的数组
static char **s_rs = NULL;
static MYREAL *rs = NULL;
static int nr=0;
// 是否silence整个输出
static bool silenceInput=false;
// 输出单个频率下的波数积分过程文件 -S
static char *s_statsdir = NULL; // 保存目录，和SAC文件目录同级
static char **s_statsidxs = NULL;
static MYINT nstatsidxs=0;
static MYINT *statsidxs = NULL;
// 计算哪些格林函数，确定震源类型, 默认计算全部
static bool doEX=true, doVF=true, doHF=true, doDC=true;

// 是否计算位移空间导数
static bool calc_upar=false;

// 各选项的标志变量，初始化为0，定义了则为1
static int M_flag=0, D_flag=0, N_flag=0, 
            O_flag=0, H_flag=0,
            L_flag=0, V_flag=0, E_flag=0, 
            K_flag=0, s_flag=0, 
            S_flag=0, R_flag=0, P_flag=0,
            G_flag=0, e_flag=0;


/**
 * 打印使用说明
 */
static void print_help(){
print_logo();
printf("\n"
"+ To get more precise results when source and receiver are \n"
"  at a close or same depth, Peak-Trough Average Method(PTAM)\n"
"  (Zhang et al., 2003) will be applied automatically.\n"
"\n"
"+ To use large dk to increase computing speed at a large\n"
"  epicentral distance, Filon's Integration Method(FIM) with \n"
"  2-point linear interpolation(Ji and Yao, 1995) and \n"
"  Self Adaptive FIM (SAFIM) (Chen and Zhang, 2001) can be applied.\n" 
"\n\n"
"The units of output Green's Functions for different sources are: \n"
"    + Explosion:     1e-20 cm/(dyne-cm)\n"
"    + Single Force:  1e-15 cm/(dyne)\n"
"    + Shear:         1e-20 cm/(dyne-cm)\n"
"    + Moment Tensor: 1e-20 cm/(dyne-cm)\n" 
"\n\n"
"The components of Green's Functions are :\n"
"     +------+-----------------------------------------------+\n"
"     | Name |       Description (Source, Component)         |\n"
"     +------+-----------------------------------------------+\n"
"     | EXZ  | Explosion, Vertical Upward                    |\n"
"     | EXR  | Explosion, Radial Outward                     |\n"
"     | VFZ  | Vertical Downward Force, Vertical Upward      |\n"
"     | VFR  | Vertical Downward Force, Radial Outward       |\n"
"     | HFZ  | Horizontal Force, Vertical Upward             |\n"
"     | HFR  | Horizontal Force, Radial Outward              |\n"
"     | HFT  | Horizontal Force, Transverse Clockwise        |\n"
"     | DDZ  | 45° dip slip, Vertical Upward                 |\n"
"     | DDR  | 45° dip slip, Radial Outward                  |\n"
"     | DSZ  | 90° dip slip, Vertical Upward                 |\n"
"     | DSR  | 90° dip slip, Radial Outward                  |\n"
"     | DST  | 90° dip slip, Transverse Clockwise            |\n"
"     | SSZ  | Vertical strike slip, Vertical Upward         |\n"
"     | SSR  | Vertical strike slip, Radial Outward          |\n"
"     | SST  | Vertical strike slip, Transverse Clockwise    |\n"
"     +------+-----------------------------------------------+\n"
"\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt -M<model> -D<depsrc>/<deprcv> -N<nt>/<dt>[/<zeta>] \n"
"        -R<r1>,<r2>[,...]    [-O<outdir>]    [-H<f1>/<f2>] \n"
"        [-L<length>]    [-V<vmin_ref>]     [-E<t0>[/<v0>]] \n" 
"        [-K<k0>[/<ampk>/<keps>]]            [-P<nthreads>]\n"
"        [-G<b1>[/<b2>/<b3>/<b4>]] [-S<i1>,<i2>[,...]] [-e]\n"
"        [-s]\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -M<model>    Filepath to 1D horizontally layered halfspace \n"
"                 model. The model file has 6 columns: \n"
"\n"
"         +-------+----------+----------+-------------+----+----+\n"
"         | H(km) | Vp(km/s) | Vs(km/s) | Rho(g/cm^3) | Qp | Qa |\n"
"         +-------+----------+----------+-------------+----+----+\n"
"\n"
"                 and the number of layers are unlimited.\n"
"\n"
"    -D<depsrc>/<deprcv>\n"
"                 <depsrc>: source depth (km).\n"
"                 <deprcv>: receiver depth (km).\n"
"\n"
"    -N<nt>/<dt>[/<zeta>]\n"
"                 <nt>:   number of points. (NOT requires 2^n).\n"
"                 <dt>:   time interval (secs). \n"
"                 <zeta>: define the coefficient of imaginary \n"
"                         frequency wI=zeta*PI/T, where T=nt*dt.\n"
"                         Default zeta=%.1f.\n", zeta); printf(
"\n"
"    -R<r1>,<r2>[,...]\n"
"                 Multiple epicentral distance (km), \n"
"                 seperated by comma.\n"
"\n"
"    -O<outdir>   Directorypath of output for saving. Default is\n"
"                 current directory.\n"
"\n"
"    -H<f1>/<f2>  Apply bandpass filer with rectangle window, \n"
"                 default no filter.\n"
"                 <f1>: lower frequency (Hz), %.1f means low pass.\n", freq1); printf(
"                 <f2>: upper frequency (Hz), %.1f means high pass.\n", freq2); printf(
"\n"
"    -L[a]<length>[/<Flength>/<Fcut>]\n"
"                 Define the wavenumber integration interval\n"
"                 dk=(2*PI)/(<length>*rmax). rmax is the maximum \n"
"                 epicentral distance. \n"
"                 There are 4 cases:\n"
"                 + (default) not set or set %.1f.\n", Length); printf(
"                   <length> will be determined automatically\n"
"                   in program with the criterion (Bouchon, 1980).\n"
"                 + manually set one POSITIVE value, e.g. -L20\n"
"                 + manually set three POSITIVE values, \n"
"                   e.g. -L20/5/10, means split the integration \n"
"                   into two parts, [0, k*] and [k*, kmax], \n"
"                   in which k*=<Fcut>/rmax, and use DWM with\n"
"                   <length> and FIM with <Flength>, respectively.\n"
"                 + manually set three POSITIVE values, with -La,\n"
"                   in this case, <Flength> will be <Ftol> for Self-\n"
"                   Adaptive FIM.\n"
"\n"
"    -V<vmin_ref> \n"
"                 Minimum velocity (km/s) for reference. This\n"
"                 is designed to define the upper bound \n"
"                 of wavenumber integration, see the\n"
"                 description of -K for the specific formula.\n"
"                 There are 3 cases:\n"
"                 + (default) not set or set %.1f.\n", vmin_ref); printf(
"                   <vmin_ref> will be the minimum velocity\n"
"                   of model, but limited to %.1f. and if the \n", min_vmin_ref); printf(
"                   depth gap between source and receiver is \n"
"                   thinner than %.1f km, PTAM will be appled\n", hs_ptam); printf(
"                   automatically.\n"
"                 + manually set POSITIVE value. \n"
"                 + manually set NEGATIVE value, \n"
"                   and PTAM will be appled.\n"
"\n"
"    -E<t0>[/<v0>]\n"
"                 Introduce the time delay in results. The total \n"
"                 delay = <t0> + dist/<v0>, dist is the\n"
"                 straight-line distance between source and \n"
"                 receiver.\n"
"                 <t0>: reference delay (s), default t0=%.1f\n", delayT0); printf(
"                 <v0>: reference velocity (km/s), \n"
"                       default %.1f not use.\n", delayV0); printf(
"\n"
"    -K<k0>[/<ampk>/<keps>]\n"
"                 Several parameters designed to define the\n"
"                 behavior in wavenumber integration. The upper\n"
"                 bound is \n"
"                 sqrt( (<k0>*mult)^2 + (<ampk>*w/<vmin_ref>)^2 ),\n"
"                 default mult=1.0.\n"
"                 <k0>:   designed to give residual k at\n"
"                         0 frequency, default is %.1f, and \n", k0); printf(
"                         multiply PI/hs in program, \n"
"                         where hs = max(fabs(depsrc-deprcv), %.1f).\n", MIN_DEPTH_GAP_SRC_RCV); printf(
"                 <ampk>: amplification factor, default is %.2f.\n", ampk); printf(
"                 <keps>: a threshold for break wavenumber \n"
"                         integration in advance. See \n"
"                         (Yao and Harkrider, 1983) for details.\n"
"                         Default %.1f not use.\n", keps); printf(
"\n"
"    -P<n>        Number of threads. Default use all cores.\n"
"\n"
"    -G<b1>[/<b2>/<b3>/<b4>]\n"
"                 Designed to choose which kind of source's Green's \n"
"                 functions will be computed, default is all (%d/%d/%d/%d). \n", (int)doEX, (int)doVF, (int)doHF, (int)doDC); printf(
"                 Four bool type (0 or 1) options are\n"
"                 <b1>: Explosion (EX)\n"
"                 <b2>: Vertical Force (VF)\n"
"                 <b3>: Horizontal Force (HF)\n"
"                 <b4>: Shear (DC)\n"
"\n"
"    -S<i1>,<i2>[,...]\n"
"                 Frequency (index) of statsfile in wavenumber\n"
"                 integration to be output, require 0 <= i <= nf-1,\n"
"                 where nf=nt/2+1. These option is designed to check\n"
"                 the trend of kernel with wavenumber.\n"
"                 -1 means all frequency index.\n"
"\n"
"    -e           Compute the spatial derivatives, ui_z and ui_r,\n"
"                 of displacement u. In filenames, prefix \"r\" means \n"
"                 ui_r and \"z\" means ui_z. The units of derivatives\n"
"                 for different sources are: \n"
"                 + Explosion:     1e-25 /(dyne-cm)\n"
"                 + Single Force:  1e-20 /(dyne)\n"
"                 + Shear:         1e-25 /(dyne-cm)\n"
"                 + Moment Tensor: 1e-25 /(dyne-cm)\n" 
"\n"
"    -s           Silence all outputs.\n"
"\n"
"    -h           Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt -Mmilrow -N1000/0.01 -D2/0 -Ores -R2,4,6,8,10\n"
"\n\n\n"
);

}


/**
 * 从路径字符串中找到用/或\\分隔的最后一项
 * 
 * @param    path     路径字符串指针
 * 
 * @return   指向最后一项字符串的指针
 */
static char* get_basename(char* path) {
    // 找到最后一个 '/'
    char* last_slash = strrchr(path, '/'); 
    
#ifdef _WIN32
    char* last_backslash = strrchr(path, '\\');
    if (last_backslash && (!last_slash || last_backslash > last_slash)) {
        last_slash = last_backslash;
    }
#endif
    if (last_slash) {
        // 返回最后一个 '/' 之后的部分
        return last_slash + 1; 
    }
    // 如果没有 '/'，整个路径就是最后一项
    return path; 
}

/**
 * 从命令行中读取选项，处理后记录到全局变量中
 * 
 * @param     argc      命令行的参数个数
 * @param     argv      多个参数字符串指针
 */
static void getopt_from_command(int argc, char **argv){
    int opt;
    while ((opt = getopt(argc, argv, ":M:D:N:O:H:L:V:E:K:shR:S:P:G:e")) != -1) {
        switch (opt) {
            // 模型路径，其中每行分别为 
            //      厚度(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  Qp   Qs
            // 互相用空格隔开即可
            case 'M':
                M_flag = 1;
                s_modelpath = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                strcpy(s_modelpath, optarg);
                if(access(s_modelpath, F_OK) == -1){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! File \"%s\" set by -M not exists.\n" DEFAULT_RESTORE, command, s_modelpath);
                    exit(EXIT_FAILURE);
                }
            
                s_modelname = get_basename(s_modelpath);
                break;

            // 震源和场点深度， -Ddepsrc/deprcv
            case 'D':
                D_flag = 1;
                s_depsrc = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                s_deprcv = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                if(2 != sscanf(optarg, "%[^/]/%s", s_depsrc, s_deprcv)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -D.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                if(1 != sscanf(s_depsrc, "%lf", &depsrc)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -D.\n"  DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                if(1 != sscanf(s_deprcv, "%lf", &deprcv)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -D.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                if(depsrc < 0.0 || deprcv < 0.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Negative value in -D is not supported.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                break;

            // 点数,采样间隔,虚频率 -Nnt/dt/[zeta]
            case 'N':
                N_flag = 1;
                if(2 > sscanf(optarg, "%d/%lf/%lf", &nt, &dt, &zeta)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -N.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                if(nt <= 0 || dt <= 0.0 || zeta <= 0.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Nonpositive value in -N is not supported.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                break;

            // 输出路径 -Ooutput_dir
            case 'O':
                O_flag = 1;
                s_output_dir = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                strcpy(s_output_dir, optarg);
                break;

            // 频带 -H f1/f2
            case 'H':
                H_flag = 1;
                if(2 != sscanf(optarg, "%lf/%lf", &freq1, &freq2)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -H.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                if(freq1>0.0 && freq2>0.0 && freq1 > freq2){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! In -H, positive freq1 should be less than positive freq2.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                break;

            // 波数积分间隔 -L[a]<length>[/<Flength>/<Fcut>]
            case 'L':
                L_flag = 1;
                {
                    // 检查首字母是否为a，表明使用自适应Filon积分
                    int pos=0;
                    bool useSAFIM = false;
                    if(optarg[0] == 'a'){
                        pos++;
                        useSAFIM = true;
                    }
                    double filona = 0.0;
                    int n = sscanf(optarg+pos, "%lf/%lf/%lf", &Length, &filona, &filonCut);
                    if(n != 1 && n != 3){
                        fprintf(stderr, "[%s] " BOLD_RED "Error in -L.\n" DEFAULT_RESTORE, command);
                        exit(EXIT_FAILURE);
                    };
                    if(n == 1 && Length <= 0){
                        fprintf(stderr, "[%s] " BOLD_RED "Error! In -L, length should be positive.\n" DEFAULT_RESTORE, command);
                        exit(EXIT_FAILURE);
                    }
                    if(n == 3 && (filona <= 0 || filonCut < 0)){
                        fprintf(stderr, "[%s] " BOLD_RED "Error! In -L, Flength/Ftol should be positive, Fcut should be nonnegative.\n" DEFAULT_RESTORE, command);
                        exit(EXIT_FAILURE);
                    }
                    if(n == 3){
                        if(useSAFIM){
                            safilonTol = filona;
                        } else {
                            filonLength = filona;
                        }
                    }
                }
                
                break;

            // 参考最小速度 -Vvmin_ref
            case 'V':
                V_flag = 1;
                if(0 == sscanf(optarg, "%lf", &vmin_ref)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -V.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                break;

            // 时间延迟 -ET0/V0
            case 'E':
                E_flag = 1;
                if(0 == sscanf(optarg, "%lf/%lf", &delayT0, &delayV0)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -E.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                if(delayV0 < 0.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Can't set negative v0(%f) in -E.\n" DEFAULT_RESTORE, command, delayV0);
                    exit(EXIT_FAILURE);
                }
                break;

            // 波数积分相关变量 -Kk0/ampk/keps
            case 'K':
                K_flag = 1;
                {
                    if(0 == sscanf(optarg, "%lf/%lf/%lf", &k0, &ampk, &keps)){
                        fprintf(stderr, "[%s] " BOLD_RED "Error in -K.\n" DEFAULT_RESTORE, command);
                        exit(EXIT_FAILURE);
                    };
                }
                
                if(k0 < 0.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Can't set negative k0(%f) in -K.\n" DEFAULT_RESTORE, command, k0);
                    exit(EXIT_FAILURE);
                }
                if(ampk < 0.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Can't set negative ampk(%f) in -K.\n" DEFAULT_RESTORE, command, ampk);
                    exit(EXIT_FAILURE);
                }
                break;

            // 不打印在终端
            case 's':
                s_flag = 1;
                silenceInput = true;
                break;

            // 震中距数组，-Rr1,r2,r3,r4 ...
            case 'R':
                R_flag = 1;
                {
                    char *token;
                    char *str_copy = strdup(optarg);  // 创建字符串副本，以免修改原始字符串
                    token = strtok(str_copy, ",");

                    while(token != NULL){
                        s_rs = (char**)realloc(s_rs, sizeof(char*)*(nr+1));
                        s_rs[nr] = NULL;
                        s_rs[nr] = (char*)realloc(s_rs[nr], sizeof(char)*(strlen(token)+1));
                        rs = (MYREAL*)realloc(rs, sizeof(MYREAL)*(nr+1));
                        strcpy(s_rs[nr], token);
                        rs[nr] = atof(token);
                        if(rs[nr] == 0.0){
                            fprintf(stderr, "[%s] " BOLD_RED "Warning! Add 1e-5 to Zero epicentral distance in -R.\n" DEFAULT_RESTORE, command);
                            rs[nr] += 1e-5;
                        }
                        if(rs[nr] < 0.0){
                            fprintf(stderr, "[%s] " BOLD_RED "Error! Can't set negative epicentral distance(%f) in -R.\n" DEFAULT_RESTORE, command, rs[nr]);
                            exit(EXIT_FAILURE);
                        }


                        token = strtok(NULL, ",");
                        nr++;
                    }
                    free(str_copy);
                }
                break;

            // 多线程数 -Pnthreads
            case 'P':
                P_flag = 1;
                if(1 != sscanf(optarg, "%d", &nthreads)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -P.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                if(nthreads <= 0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Nonpositive value in -P is not supported.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                set_num_threads(nthreads);
                break;

            // 选择要计算的格林函数 -G1/1/1/1
            case 'G': 
                G_flag = 1;
                doEX = doVF = doHF = doDC = false;
                {
                    int i1, i2, i3, i4;
                    i1 = i2 = i3 = i4 = 0;
                    if(0 == sscanf(optarg, "%d/%d/%d/%d", &i1, &i2, &i3, &i4)){
                        fprintf(stderr, "[%s] " BOLD_RED "Error in -G.\n" DEFAULT_RESTORE, command);
                        exit(EXIT_FAILURE);
                    };
                    doEX = (i1!=0);
                    doVF  = (i2!=0);
                    doHF  = (i3!=0);
                    doDC  = (i4!=0);
                }
                
                // 至少要有一个真
                if(!(doEX || doVF || doHF || doDC)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! At least set one true value in -G.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }

                break;

            // 输出波数积分中间文件， -Sidx1,idx2,idx3,...
            case 'S':
                S_flag = 1;
                {
                    char *token;
                    char *str_copy = strdup(optarg);  // 创建字符串副本，以免修改原始字符串
                    token = strtok(str_copy, ",");

                    while(token != NULL){
                        s_statsidxs = (char**)realloc(s_statsidxs, sizeof(char*)*(nstatsidxs+1));
                        s_statsidxs[nstatsidxs] = NULL;
                        s_statsidxs[nstatsidxs] = (char*)realloc(s_statsidxs[nstatsidxs], sizeof(char)*(strlen(token)+1));
                        strcpy(s_statsidxs[nstatsidxs], token);
                        statsidxs = (MYINT*)realloc(statsidxs, sizeof(MYINT)*(nstatsidxs+1));
                        statsidxs[nstatsidxs] = atoi(token);

                        token = strtok(NULL, ",");
                        nstatsidxs++;
                    }
                    free(str_copy);
                }
                break;

            // 是否计算位移空间导数
            case 'e':
                e_flag = 1;
                calc_upar = true;
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
    } // END get options

    // 检查必须设置的参数是否有设置
    if(argc == 1){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set options. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(M_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -M. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(D_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -D. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(N_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -N. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(R_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -R. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(O_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -O. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }


    // 建立保存目录
    if(mkdir(s_output_dir, 0777) != 0){
        if(errno != EEXIST){
            fprintf(stderr, "[%s] " BOLD_RED "Error! Unable to create folder %s. Error code: %d\n" DEFAULT_RESTORE, command, s_output_dir, errno);
            exit(EXIT_FAILURE);
        }
    }
    

    // 在目录中保留命令
    char *dummy = (char*)malloc(sizeof(char)*(strlen(s_output_dir)+100));
    sprintf(dummy, "%s/command", s_output_dir);
    FILE *fp = fopen(dummy, "a");
    for(int i=0; i<argc; ++i){
        fprintf(fp, "%s ", argv[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
    free(dummy);
}

/**
 * 将长数组以换行的方式打印
 * 
 * @param   nlen1     第一列字符宽度
 * @param   nlen2     第二列字符宽度
 * @param   rowname   行名
 * @param   s_arr     字符串数组
 * @param   arrsize   数组长度
 */
static void print_long_array_in_tabel(
    const int nlen1, const int nlen2, const char *rowname,
    const char **s_arr, const int arrsize)
{
    char tmp[nlen2-2];
    printf("| %-*s | ", nlen1-3, rowname);
    {
        tmp[0] = '\0';
        int len = 0;
        for(int m=0; m<arrsize; ++m){
            len = strlen(tmp);
            char s2[strlen(s_arr[m])+2]; // 加上'\0'和逗号的长度
            snprintf(s2, sizeof(s2), "%s,", s_arr[m]);

            snprintf(tmp+len, sizeof(tmp)-len, "%s", s2);
            while(sizeof(tmp) - strlen(tmp) == 1){ // 允许换行
                printf("%-*s |\n", nlen2-3, tmp);
                snprintf(tmp, sizeof(tmp), "%s", s2+(strlen(tmp)-len));
                if(strcmp(tmp, ",") == 0 && m==arrsize-1)  break;
                printf("| %-*s | ", nlen1-3, "");
                len -= strlen(tmp);
            }
            if(strcmp(tmp, ",") == 0 && m==arrsize-1)  break;
        }

        if(strlen(tmp) > 0 && strcmp(tmp, ",") != 0){
            printf("%-*s |\n", nlen2-3, tmp);
        }
    }
}


/**
 * 将某一道做ifft，做时间域处理，保存到sac文件
 * 
 * @param     delay     时间延迟
 * @param     mult      幅值放大系数
 * @param     grncplx   复数形式的格林函数频谱
 * @param     fftw_grn  将频谱写到FFTW_COMPLEX类型中
 * @param     out       ifft后的时域数据
 * @param     float_arr 将时域数据写到float类型的数组中
 * @param     plan      FFTW_PLAN
 * @param     hd        SAC头段变量结构体
 * @param     outpath   sac文件保存路径
 */
static void ifft_one_trace(
    MYREAL delay, MYREAL mult,
    MYCOMPLEX *grncplx, _FFTW_COMPLEX *fftw_grn, MYREAL *out, float *float_arr,
    _FFTW_PLAN plan, SACHEAD *hd, const char *outpath)
{
    // 赋值复数，包括时移
    MYCOMPLEX cfac, ccoef;
    cfac = exp(I*dw*delay);
    ccoef = mult;
    for(int i=0; i<nf; ++i){
        fftw_grn[i] = grncplx[i] * ccoef;
        ccoef *= cfac;
    }

    // 发起fft任务 
    _FFTW_EXECUTE(plan);

    // 归一化，并处理虚频
    double fac, coef;
    coef = df * exp(delay*wI);
    fac = exp(wI*dt);
    for(int i=0; i<nt; ++i){
        out[i] *= coef;
        coef *= fac;
    }

    // 以sac文件保存到本地
    for(int i=0; i<nt; ++i){
        float_arr[i] = out[i];
    }

    write_sac(outpath, *hd, float_arr);
    // FILE *fp = fopen(outpath, "wb");
    // fwrite(out, sizeof(float), nt, fp);
    // fclose(fp);
}


/**
 * 以表格形式打印使用的参数
 */
static void print_parameters(){
    // 模拟打两列表格，第一列参数名，第二列参数值
    print_pymod(pymod);
    const int nlen1=20, nlen2=45; // 两列字符宽度
    // 制作每行分割线
    char splitline[nlen1+nlen2+2];
    splitline[0] = '+';
    for(int i=1; i<nlen1; splitline[i]='-', ++i);
    splitline[nlen1] = '+';
    for(int i=nlen1+1; i<nlen1+nlen2; splitline[i]='-', ++i);
    splitline[nlen1+nlen2] = '+';
    splitline[nlen1+nlen2+1] = '\0';
    // 用于输出值的特殊处理
    char tmp[nlen2-2];

    printf("%s\n", splitline);
    printf("| %-*s | %-*s |\n", nlen1-3, "Parameter", nlen2-3, "Value");
    printf("%s\n", splitline);
    printf("| %-*s | %-*s |\n", nlen1-3, "model_path", nlen2-3, s_modelpath);
    printf("| %-*s | %-*f |\n", nlen1-3, "vmin", nlen2-3, vmin);
    printf("| %-*s | %-*f |\n", nlen1-3, "vmax", nlen2-3, vmax);
    // 特殊处理vmin_ref的输出
    snprintf(tmp, sizeof(tmp), "%f", fabs(vmin_ref));
    if(vmin_ref < 0.0) strncat(tmp, ", using PTAM.", sizeof(tmp)-strlen(tmp)-1);
    printf("| %-*s | %-*s |\n", nlen1-3, "vmin_ref", nlen2-3, tmp);
    // 特殊处理length的输出
    snprintf(tmp, sizeof(tmp), "%f", Length);
    if(filonLength > 0.0){  
        snprintf(tmp, sizeof(tmp), "%f,%f,%f", Length, filonLength, filonCut);
        strncat(tmp, ", using FIM.", sizeof(tmp)-strlen(tmp)-1);
    } else if(safilonTol > 0.0){
        snprintf(tmp, sizeof(tmp), "%f,%f,%f", Length, safilonTol, filonCut);
        strncat(tmp, ", using SAFIM.", sizeof(tmp)-strlen(tmp)-1);
    }
    printf("| %-*s | %-*s |\n", nlen1-3, "length", nlen2-3, tmp);
    // 
    printf("| %-*s | %-*d |\n", nlen1-3, "nt", nlen2-3, nt);
    printf("| %-*s | %-*f |\n", nlen1-3, "dt", nlen2-3, dt);
    printf("| %-*s | %-*f |\n", nlen1-3, "winT", nlen2-3, winT);
    printf("| %-*s | %-*f |\n", nlen1-3, "zeta", nlen2-3, zeta);
    printf("| %-*s | %-*f |\n", nlen1-3, "delayT0", nlen2-3, delayT0);
    printf("| %-*s | %-*f |\n", nlen1-3, "delayV0", nlen2-3, delayV0);
    printf("| %-*s | %-*f |\n", nlen1-3, "tmax", nlen2-3, tmax);
    printf("| %-*s | %-*f |\n", nlen1-3, "k0", nlen2-3, k0);
    printf("| %-*s | %-*f |\n", nlen1-3, "ampk", nlen2-3, ampk);
    printf("| %-*s | %-*f |\n", nlen1-3, "keps", nlen2-3, keps);
    printf("| %-*s | %-*f |\n", nlen1-3, "maxfreq(Hz)", nlen2-3, freqs[nf-1]);
    printf("| %-*s | %-*f |\n", nlen1-3, "f1(Hz)", nlen2-3, freqs[nf1]);
    printf("| %-*s | %-*f |\n", nlen1-3, "f2(Hz)", nlen2-3, freqs[nf2]);
    // 特殊处理震中距的输出
    print_long_array_in_tabel(nlen1, nlen2, "distances(km)", (const char **)s_rs, nr);
    
    // 特殊处理statsfile index的输出
    if(nstatsidxs > 0){
        print_long_array_in_tabel(nlen1, nlen2, "statsfile_index", (const char **)s_statsidxs, nstatsidxs);
    }
        
    printf("| %-*s | ", nlen1-3, "sources");
    tmp[0] = '\0';
    if(doEX) snprintf(tmp+strlen(tmp), sizeof(tmp)-strlen(tmp), "EX,");
    if(doVF)  snprintf(tmp+strlen(tmp), sizeof(tmp)-strlen(tmp), "VF,");
    if(doHF)  snprintf(tmp+strlen(tmp), sizeof(tmp)-strlen(tmp), "HF,");
    if(doDC)  snprintf(tmp+strlen(tmp), sizeof(tmp)-strlen(tmp), "DC,");
    printf("%-*s |\n", nlen2-3, tmp);
    
    // 特殊处理输出路径
    printf("| %-*s | ", nlen1-3, "output_path");
    {
        int len=0;
        tmp[0] = '\0';
        snprintf(tmp+strlen(tmp), sizeof(tmp)-strlen(tmp), "%s", s_output_dir);
        len = strlen(tmp);
        while(sizeof(tmp) - strlen(tmp) == 1){ // 允许换行
            printf("%-*s |\n", nlen2-3, tmp);
            snprintf(tmp, sizeof(tmp), "%s", s_output_dir+(strlen(tmp)-len));
            printf("| %-*s | ", nlen1-3, "");
            len -= strlen(tmp);
        }
        if(strlen(tmp) > 0){
            printf("%-*s |\n", nlen2-3, tmp);
        }
    }
    printf("%s\n", splitline);
    printf("\n");
    
}


/**
 * 以表格形式打印输出的文件夹和走时
 * 
 * @param     s_output_subdir     输出文件夹
 * @param     s_R                 字符串形式的震中距
 * @param     Tp                  P波初至走时
 * @param     Ts                  S波初至走时
 */
static void print_outdir_travt(const char *s_output_subdir, const char *s_R, double Tp, double Ts){
    // 进入该函数的次数 
    static int numin = 0;

    static const int ncols = 4;
    static const int nlens[] = {28, 17, 13, 13};
    int Nlen=0;
    for(int ic=0; ic<ncols; ++ic){
        Nlen += nlens[ic]; 
    }
    // 定义分割线
    char splitline[Nlen+2];
    {
        int n=0;
        for(int ic=0; ic<ncols; ++ic){
            splitline[n] = '+';
            for(int i=1; i<nlens[ic]; ++i){
                splitline[n + i] = '-';
            }
            n += nlens[ic];
        }
        splitline[Nlen] = '+';
        splitline[Nlen+1] = '\0';
    }
    
    // 第一次执行该函数，打印题头
    if(numin == 0){
        printf("\n");
        printf("%s\n", splitline);
        printf("| %-*s ", nlens[0]-3, "Output Directory");
        printf("| %-*s ", nlens[1]-3, "Distance(km)");
        printf("| %-*s ", nlens[2]-3, "Tp(secs)");
        printf("| %-*s ", nlens[2]-3, "Ts(secs)");
        printf("|\n");
        printf("%s\n", splitline);
    }


    // 打印目录、震中距、走时P、走时S
    // 目录和震中距的字符串可能分行，且行数不一
    // 这里不断分行打印，直到打印完
    char dirtmp[nlens[0]-2];
    char Rtmp[nlens[1]-2];
    int lendir=0, lenR=0;
    int LENdir=strlen(s_output_subdir);
    int LENR=strlen(s_R);
    dirtmp[0] = '\0';
    Rtmp[0] = '\0';
    int iline=0;
    int idir=0, iR=0;

    snprintf(dirtmp, sizeof(dirtmp), "%s", s_output_subdir);
    snprintf(Rtmp, sizeof(Rtmp), "%s", s_R);
    lendir = strlen(dirtmp);
    lenR = strlen(Rtmp);
    do{
        printf("| %-*s ", nlens[0]-3, dirtmp);
        idir += strlen(dirtmp);
        if(sizeof(dirtmp) - strlen(dirtmp) == 1 && lendir < LENdir){
            snprintf(dirtmp, sizeof(dirtmp), "%s", s_output_subdir+lendir);
            lendir += strlen(dirtmp);
        } else {
            dirtmp[0] = '\0';
        }
        
        printf("| %-*s ", nlens[1]-3, Rtmp);
        iR += strlen(Rtmp);
        if(sizeof(Rtmp) - strlen(Rtmp) == 1 && lenR < LENR){
            snprintf(Rtmp, sizeof(Rtmp), "%s", s_R+lenR);
            lenR += strlen(Rtmp);
        } else {
            Rtmp[0] = '\0';
        }

        if(iline == 0){
            printf("| %-*.3f ", nlens[2]-3, Tp);
            printf("| %-*.3f ", nlens[3]-3, Ts);
        } else {
            printf("| %-*s ", nlens[2]-3, " ");
            printf("| %-*s ", nlens[3]-3, " ");
        }
        printf("|\n");

        iline++;

        // printf("iR=%d, idir=%d, LENDIR=%d\n", iR, idir, LENdir);

    }while(iR < LENR || idir < LENdir);

    if(numin == nr-1){
        printf("%s\n", splitline);
        printf("\n");
    }
    

    numin++;
}


/**
 * 将一条数据反变换回时间域再进行处理，保存到SAC文件
 * 
 * @param     srcname       震源类型
 * @param     ch            三分量类型（Z,R,T）
 * @param     hd            SAC头段变量结构体指针
 * @param     s_outpath     用于接收保存路径字符串
 * @param     s_output_subdir    保存路径所在文件夹
 * @param     s_prefix           sac文件名以及通道名名称前缀
 * @param     sgn                数据待乘符号(-1/1)
 * @param     grncplx   复数形式的格林函数频谱
 * @param     fftw_grn  将频谱写到FFTW_COMPLEX类型中
 * @param     out       ifft后的时域数据
 * @param     float_arr 将时域数据写到float类型的数组中
 * @param     plan      FFTW_PLAN
 * 
 */
static void write_one_to_sac(
    const char *srcname, const char ch, 
    SACHEAD *hd, char *s_outpath, const char *s_output_subdir, const char *s_prefix,
    const int sgn, MYCOMPLEX *grncplx, fftw_complex *fftw_grn, MYREAL *out, float *float_arr, fftw_plan plan)
{
    char kcmpnm[9];
    snprintf(kcmpnm, sizeof(kcmpnm), "%s%s%c", s_prefix, srcname, ch);
    strcpy(hd->kcmpnm, kcmpnm);
    sprintf(s_outpath, "%s/%s.sac", s_output_subdir, kcmpnm);
    ifft_one_trace(delayT, sgn, grncplx, fftw_grn, out, float_arr, plan, hd, s_outpath);
}



//====================================================================================
//====================================================================================
//====================================================================================
int main(int argc, char **argv) {
    command = argv[0];

    // 传入参数 
    getopt_from_command(argc, argv);

    // 读入模型文件
    if((pymod = read_pymod_from_file(command, s_modelpath, depsrc, deprcv, true)) ==NULL){
        exit(EXIT_FAILURE);
    }

    // 当震源位于液体层中时，仅允许计算爆炸源对应的格林函数
    // 程序结束前会输出对应警告
    if(pymod->Vb[pymod->isrc]==0.0){
        doHF = doVF = doDC = false;
    }

    // 最大最小速度
    get_pymod_vmin_vmax(pymod, &vmin, &vmax);

    // 参考最小速度
    if(vmin_ref == 0.0){
        vmin_ref = vmin;
        if(vmin_ref < min_vmin_ref) vmin_ref = min_vmin_ref;
    } 

    // 如果没有主动设置vmin_ref，则判断是否要自动使用PTAM
    if(V_flag == 0 && fabs(deprcv - depsrc) <= hs_ptam) {
        vmin_ref = - fabs(vmin_ref);
    }

    // 时窗长度 
    winT = nt*dt;

    // 最大震中距
    MYREAL rmax=rs[findMinMax_MYREAL(rs, nr, true)];   

    // 时窗最大截止时刻
    tmax = delayT0 + winT;
    if(delayV0 > 0.0) tmax += rmax/delayV0;

    // 自动选择积分间隔，默认使用传统离散波数积分
    // 自动选择会给出很保守的值（较大的Length）
    if(Length == 0.0){
        Length = 15.0; 
        double jus = (vmax*tmax)*(vmax*tmax) - (deprcv-depsrc)*(deprcv-depsrc);
        if(jus >= 0.0){
            Length = 1.0 + sqrt(jus)/rmax + 0.5; // +0.5为保守值
            if(Length < 15.0) Length = 15.0;
        }
    }

    // 虚频率
    wI = zeta*PI/winT;

    // 定义要计算的频率、时窗等
    nf = nt/2 + 1;
    df = 1.0/winT;
    dw = 2.0*PI*df;
    freqs = (MYREAL*)malloc(nf*sizeof(MYREAL));
    for(int i=0; i<nf; ++i){
        freqs[i] = i*df;
    }

    nf1 = 0; nf2 = nf-1;
    if(freq1 > 0.0){
        nf1 = ceil(freq1/df);
        if(nf1 >= nf-1)    nf1 = nf-1;
    }
    if(freq2 > 0.0){
        nf2 = floor(freq2/df);
        if(nf2 >= nf-1)    nf2 = nf-1;
    }
    if(nf2 < nf1) nf2 = nf1;

    // 波数积分中间文件输出目录
    if(nstatsidxs > 0){
        s_statsdir = (char*)malloc(sizeof(char)*(strlen(s_modelpath)+strlen(s_output_dir)+strlen(s_depsrc)+strlen(s_deprcv)+100));
        sprintf(s_statsdir, "%s_grtstats", s_output_dir);
        // 建立保存目录
        if(mkdir(s_statsdir, 0777) != 0){
            if(errno != EEXIST){
                fprintf(stderr, "[%s] " BOLD_RED "Error! Unable to create folder %s. Error code: %d\n" DEFAULT_RESTORE, command, s_statsdir, errno);
                exit(EXIT_FAILURE);
            }
        }
        sprintf(s_statsdir, "%s/%s_%s_%s", s_statsdir, s_modelname, s_depsrc, s_deprcv);
        if(mkdir(s_statsdir, 0777) != 0){
            if(errno != EEXIST){
                fprintf(stderr, "[%s] " BOLD_RED "Error! Unable to create folder %s. Error code: %d\n" DEFAULT_RESTORE, command, s_statsdir, errno);
                exit(EXIT_FAILURE);
            }
        }
    }
    

    // 建立格林函数的complex数组
    MYCOMPLEX *(*grn)[SRC_M_NUM][CHANNEL_NUM] = (MYCOMPLEX*(*)[SRC_M_NUM][CHANNEL_NUM]) calloc(nr, sizeof(*grn));
    MYCOMPLEX *(*grn_uiz)[SRC_M_NUM][CHANNEL_NUM] = (calc_upar)? (MYCOMPLEX*(*)[SRC_M_NUM][CHANNEL_NUM]) calloc(nr, sizeof(*grn_uiz)) : NULL;
    MYCOMPLEX *(*grn_uir)[SRC_M_NUM][CHANNEL_NUM] = (calc_upar)? (MYCOMPLEX*(*)[SRC_M_NUM][CHANNEL_NUM]) calloc(nr, sizeof(*grn_uir)) : NULL;

    for(int ir=0; ir<nr; ++ir){
        for(int i=0; i<SRC_M_NUM; ++i){
            for(int c=0; c<CHANNEL_NUM; ++c){
                grn[ir][i][c] = (MYCOMPLEX*)calloc(nf, sizeof(MYCOMPLEX));
                if(grn_uiz)  grn_uiz[ir][i][c] = (MYCOMPLEX*)calloc(nf, sizeof(MYCOMPLEX));
                if(grn_uir)  grn_uir[ir][i][c] = (MYCOMPLEX*)calloc(nf, sizeof(MYCOMPLEX));
            }
        }
    }


    // 在计算前打印所有参数
    if(! silenceInput){
        print_parameters();
    }
    

    //==============================================================================
    // 计算格林函数
    integ_grn_spec(
        pymod, nf1, nf2, freqs, nr, rs, wI,
        vmin_ref, keps, ampk, k0, Length, filonLength, safilonTol, filonCut, !silenceInput,
        grn, calc_upar, grn_uiz, grn_uir,
        s_statsdir, nstatsidxs, statsidxs
    );
    //==============================================================================
    

    // 使用fftw3做反傅里叶变换
    // 分配fftw_complex内存
    _FFTW_COMPLEX *fftw_grn = (_FFTW_COMPLEX*)_FFTW_MALLOC(sizeof(_FFTW_COMPLEX)*nf);
    MYREAL *out = (MYREAL*)malloc(sizeof(MYREAL)*nt);
    float *float_arr = (float*)malloc(sizeof(float)*nt);

    // fftw计划
    _FFTW_PLAN plan = _FFTW_PLAN_DFT_C2R_1D(nt, fftw_grn, out, FFTW_ESTIMATE);
    
    // 建立SAC头文件，包含必要的头变量
    SACHEAD hd = new_sac_head(dt, nt, delayT0);
    // 发震时刻作为参考时刻
    hd.o = 0.0; 
    hd.iztype = IO; 
    // 记录震源和台站深度
    hd.evdp = depsrc; // km
    hd.stel = (-1.0)*deprcv*1e3; // m
    // 写入虚频率
    hd.user0 = wI;
    // 写入接受点的Vp,Vs,rho
    hd.user1 = pymod->Va[pymod->ircv];
    hd.user2 = pymod->Vb[pymod->ircv];
    hd.user3 = pymod->Rho[pymod->ircv];
    hd.user4 = RONE/pymod->Qa[pymod->ircv];
    hd.user5 = RONE/pymod->Qb[pymod->ircv];
    // 写入震源点的Vp,Vs,rho
    hd.user6 = pymod->Va[pymod->isrc];
    hd.user7 = pymod->Vb[pymod->isrc];
    hd.user8 = pymod->Rho[pymod->isrc];

    
    // 做反傅里叶变换，保存SAC文件
    for(int ir=0; ir<nr; ++ir){
        hd.dist = rs[ir];

        // 文件保存子目录
        char *s_output_subdir = (char*)malloc(sizeof(char)*(
            strlen(s_output_dir)+strlen(s_modelpath)+
            strlen(s_depsrc)+strlen(s_deprcv)+strlen(s_rs[ir])+100));
        
        sprintf(s_output_subdir, "%s/%s_%s_%s_%s", s_output_dir, s_modelname, s_depsrc, s_deprcv, s_rs[ir]);
        if(mkdir(s_output_subdir, 0777) != 0){
            if(errno != EEXIST){
                fprintf(stderr, "[%s] " BOLD_RED "Error! Unable to create folder %s. Error code: %d\n" DEFAULT_RESTORE, command, s_output_subdir, errno);
                exit(EXIT_FAILURE);
            }
        }
        // 时间延迟 
        delayT = delayT0;
        if(delayV0 > 0.0) delayT += sqrt(rs[ir]*rs[ir] + (deprcv-depsrc)*(deprcv-depsrc))/delayV0;
        // 修改SAC头段时间变量
        hd.b = delayT;

        // 计算理论走时
        hd.t0 = compute_travt1d(pymod->Thk, pymod->Va, pymod->n, pymod->isrc, pymod->ircv, rs[ir]);
        strcpy(hd.kt0, "P");
        hd.t1 = compute_travt1d(pymod->Thk, pymod->Vb, pymod->n, pymod->isrc, pymod->ircv, rs[ir]);
        strcpy(hd.kt1, "S");

        for(int im=0; im<SRC_M_NUM; ++im){
            if(!doEX && im==0)  continue;
            if(!doVF  && im==1)  continue;
            if(!doHF  && im==2)  continue;
            if(!doDC  && im>=3)  continue;

            int modr = SRC_M_ORDERS[im];
            int sgn=1;  // 用于反转Z分量
            for(int c=0; c<CHANNEL_NUM; ++c){
                if(modr==0 && ZRTchs[c]=='T')  continue;  // 跳过输出0阶的T分量

                // 文件保存总路径
                // char *s_outpath = (char*)malloc(sizeof(char)*(strlen(s_output_dir)+100));
                char *s_outpath = (char*)malloc(sizeof(char)*(strlen(s_output_subdir)+100));
                // char *s_suffix = (char*)malloc(sizeof(char)*(strlen(s_depsrc)+strlen(s_deprcv)+strlen(s_rs[ir])+100));
                // sprintf(s_suffix, "%s_%s_%s", s_depsrc, s_deprcv, s_rs[ir]);
                char s_prefix[] = "";

                // Z分量反转
                sgn = (ZRTchs[c]=='Z') ? -1 : 1;

                write_one_to_sac(SRC_M_NAME_ABBR[im], ZRTchs[c], &hd, s_outpath, s_output_subdir, s_prefix, sgn, grn[ir][im][c], fftw_grn, out, float_arr, plan);
                if(calc_upar){
                    write_one_to_sac(SRC_M_NAME_ABBR[im], ZRTchs[c], &hd, s_outpath, s_output_subdir, "z", sgn*(-1), grn_uiz[ir][im][c], fftw_grn, out, float_arr, plan);
                    write_one_to_sac(SRC_M_NAME_ABBR[im], ZRTchs[c], &hd, s_outpath, s_output_subdir, "r", sgn, grn_uir[ir][im][c], fftw_grn, out, float_arr, plan);
                }

                free(s_outpath);
            }
        }


        if(!silenceInput){
            print_outdir_travt(s_output_subdir, s_rs[ir], hd.t0, hd.t1);
        }

        free(s_output_subdir);
    }

    // 输出警告：当震源位于液体层中时，仅允许计算爆炸源对应的格林函数
    if(pymod->Vb[pymod->isrc]==0.0){
        fprintf(stderr, "[%s] " BOLD_YELLOW 
            "The source is located in the liquid layer, "
            "therefore only the Green's Funtions for the Explosion source will be computed.\n" 
            DEFAULT_RESTORE, command);
    }
    

    // 释放内存
    free(s_modelpath);
    free(s_depsrc);
    free(s_deprcv);
    free(s_output_dir);

    for(int ir=0; ir<nr; ++ir){
        free(s_rs[ir]);
        for(int i=0; i<SRC_M_NUM; ++i){
            for(int c=0; c<CHANNEL_NUM; ++c){
                free(grn[ir][i][c]);
                if(grn_uiz)  free(grn_uiz[ir][i][c]);
                if(grn_uir)  free(grn_uir[ir][i][c]);
            }
        }
    }
    free(grn);
    if(grn_uiz)  free(grn_uiz);
    if(grn_uir)  free(grn_uir);

    free(s_rs);
    free(rs);
    free(s_statsdir);
    free(statsidxs);
    for(int i=0; i<nstatsidxs; ++i){
        free(s_statsidxs[i]);
    }
    free(s_statsidxs);
    free(freqs);
    _FFTW_FREE(fftw_grn);
    free(out);
    free(float_arr);
    _FFTW_DESTROY_PLAN(plan);

    free_pymod(pymod);


    return 0;
}


