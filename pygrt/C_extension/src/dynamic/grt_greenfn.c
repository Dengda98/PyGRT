/**
 * @file   grt_greenfn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-11-28
 * 
 *    定义main函数，形成命令行式的用法（不使用python的entry_points，会牺牲性能）
 *    计算不同震源的格林函数
 * 
 */

#include <complex.h>
#include <fftw3.h>
#include <omp.h>

#include "dynamic/grn.h"
#include "dynamic/signals.h"
#include "travt/travt.h"
#include "common/const.h"
#include "common/model.h"
#include "common/search.h"
#include "common/sacio.h"
#include "common/util.h"

#include "grt.h"


// 一些变量的非零默认值
#define GRT_GREENFN_N_ZETA        0.8
#define GRT_GREENFN_H_FREQ1      -1.0
#define GRT_GREENFN_H_FREQ2      -1.0
#define GRT_GREENFN_V_VMIN_REF    0.1
#define GRT_GREENFN_K_K0          5.0
#define GRT_GREENFN_K_AMPK       1.15
#define GRT_GREENFN_G_EX       true
#define GRT_GREENFN_G_VF       true
#define GRT_GREENFN_G_HF       true
#define GRT_GREENFN_G_DC       true



/** 该子模块的参数控制结构体 */
typedef struct {
    char *name;
    /** 输入模型 */
    struct {
        bool active;
        char *s_modelpath;        ///< 模型路径
        const char *s_modelname;  ///< 模型名称
        PYMODEL1D *pymod;         ///< 模型PYMODEL1D结构体指针
    } M;
    /** 震源和接收器深度 */
    struct {
        bool active;
        MYREAL depsrc;
        MYREAL deprcv;
        char *s_depsrc;
        char *s_deprcv;
    } D;
    /** 波形时窗 */
    struct {
        bool active;
        MYINT nt;
        MYINT nf;
        MYREAL dt;
        MYREAL df;
        MYREAL winT;  ///< 时窗长度 
        MYREAL zeta;  ///< 虚频率系数， w <- w - zeta*PI/r* 1j
        MYREAL wI;    ///< 虚频率  zeta*PI/r
        MYREAL *freqs;
    } N;
    /** 输出目录 */
    struct {
        bool active;
        char *s_output_dir;
    } O;
    /** 频段 */
    struct {
        bool active;
        MYREAL freq1;
        MYREAL freq2;
        MYINT nf1;
        MYINT nf2;
    } H;
    /** 波数积分间隔 */
    struct {
        bool active;
        MYREAL Length;
        MYREAL filonLength;
        MYREAL safilonTol;
        MYREAL filonCut;
    } L;
    /** 波数积分上限 */
    struct {
        bool active;
        MYREAL keps;
        MYREAL ampk;
        MYREAL k0;
    } K;
    /** 参考速度 */
    struct {
        bool active;
        MYREAL vmin_ref;
    } V;
    /** 时间延迟 */
    struct {
        bool active;
        MYREAL delayT0;
        MYREAL delayV0;
    } E;
    /** 波数积分过程的核函数文件 */
    struct {
        bool active;
        char *s_raw;
        char **s_statsidxs;
        MYINT *statsidxs;
        MYINT nstatsidxs;
        char *s_statsdir;  ///< 保存目录，和SAC文件目录同级
    } S;
    /** 震中距 */
    struct {
        bool active;
        char *s_raw;
        char **s_rs;
        MYREAL *rs;
        MYINT nr;
    } R;
    /** 多线程 */
    struct {
        bool active;
        MYINT nthreads; ///< 线程数
    } P;
    /** 输出哪些震源的格林函数 */
    struct {
        bool active;
        bool doEX;
        bool doVF;
        bool doHF;
        bool doDC;
    } G;
    /** 是否计算空间导数 */
    struct {
        bool active;
    } e;
    /** 静默输出 */
    struct {
        bool active;
    } s;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    free(Ctrl->name);

    // M
    free(Ctrl->M.s_modelpath);
    free_pymod(Ctrl->M.pymod);
    
    // D
    free(Ctrl->D.s_depsrc);
    free(Ctrl->D.s_deprcv);

    // N
    free(Ctrl->N.freqs);

    // O
    free(Ctrl->O.s_output_dir);

    // R
    free(Ctrl->R.s_raw);
    for(int ir=0; ir<Ctrl->R.nr; ++ir){
        free(Ctrl->R.s_rs[ir]);
    }
    free(Ctrl->R.s_rs);
    free(Ctrl->R.rs);

    // S
    if(Ctrl->S.active){
        free(Ctrl->S.s_raw);
        for(int i=0; i<Ctrl->S.nstatsidxs; ++i){
            free(Ctrl->S.s_statsidxs[i]);
        }
        free(Ctrl->S.s_statsidxs);
        free(Ctrl->S.statsidxs);
        free(Ctrl->S.s_statsdir);
    }

    free(Ctrl);
}


/** 打印结构体中的参数 */
static void print_Ctrl(const GRT_MODULE_CTRL *Ctrl){
    print_pymod(Ctrl->M.pymod);

    const char format[]      = "   \%-20s  \%s\n";
    const char format_real[] = "   \%-20s  \%.3f\n";
    const char format_int[]  = "   \%-20s  \%d\n";
    char line[100];
    printf("------------------------------------------------\n");
    printf(format, "PARAMETER", "VALUE");
    printf(format, "model_path", Ctrl->M.s_modelpath);
    if(Ctrl->V.vmin_ref < 0.0){
        snprintf(line, sizeof(line), "%.3f, Using PTAM", Ctrl->V.vmin_ref);
    } else {
        snprintf(line, sizeof(line), "%.3f", Ctrl->V.vmin_ref);
    }
    printf(format, "vmin_ref", line);
    if(Ctrl->L.filonLength > 0.0){  
        snprintf(line, sizeof(line), "%.3f,%.3f,%.3f, using FIM", Ctrl->L.Length, Ctrl->L.filonLength, Ctrl->L.filonCut);
    } else if(Ctrl->L.safilonTol > 0.0){
        snprintf(line, sizeof(line), "%.3f,%.3e,%.3f, using SAFIM.", Ctrl->L.Length, Ctrl->L.safilonTol, Ctrl->L.filonCut);
    } else {
        snprintf(line, sizeof(line), "%.3f", Ctrl->L.Length);
    }
    printf(format, "Length", line);
    printf(format_int, "nt", Ctrl->N.nt);
    printf(format_real, "dt", Ctrl->N.dt);
    printf(format_real, "winT", Ctrl->N.winT);
    printf(format_real, "zeta", Ctrl->N.zeta);
    printf(format_real, "delayT0", Ctrl->E.delayT0);
    printf(format_real, "delayV0", Ctrl->E.delayV0);
    printf(format_real, "tmax", Ctrl->E.delayT0 + Ctrl->N.winT);
    printf(format_real, "k0", Ctrl->K.k0);
    printf(format_real, "ampk", Ctrl->K.ampk);
    printf(format_real, "keps", Ctrl->K.keps);
    printf(format_real, "maxfreq(Hz)", Ctrl->N.freqs[Ctrl->N.nf-1]);
    printf(format_real, "f1(Hz)", Ctrl->N.freqs[Ctrl->H.nf1]);
    printf(format_real, "f2(Hz)", Ctrl->N.freqs[Ctrl->H.nf2]);
    printf(format, "distances(km)", Ctrl->R.s_raw);
    if(Ctrl->S.nstatsidxs > 0){
        printf(format, "statsfile_index", Ctrl->S.s_raw);
    }
    line[0] = '\0';
    if(Ctrl->G.doEX) snprintf(line+strlen(line), sizeof(line)-strlen(line), "EX,");
    if(Ctrl->G.doVF)  snprintf(line+strlen(line), sizeof(line)-strlen(line), "VF,");
    if(Ctrl->G.doHF)  snprintf(line+strlen(line), sizeof(line)-strlen(line), "HF,");
    if(Ctrl->G.doDC)  snprintf(line+strlen(line), sizeof(line)-strlen(line), "DC,");
    printf(format, "sources", line);
    
    printf("------------------------------------------------\n");

    printf("\n\n");
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt greenfn] %s\n\n", GRT_VERSION);printf(
"    Compute the Green's Functions in Horizontally Layered\n"
"    Halfspace Model.\n"
"\n\n"
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
"    grt greenfn -M<model> -D<depsrc>/<deprcv> -N<nt>/<dt>[/<zeta>] \n"
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
"                         Default zeta=%.1f.\n", GRT_GREENFN_N_ZETA); printf(
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
"                 <f1>: lower frequency (Hz), %.1f means low pass.\n", GRT_GREENFN_H_FREQ1); printf(
"                 <f2>: upper frequency (Hz), %.1f means high pass.\n", GRT_GREENFN_H_FREQ2); printf(
"\n"
"    -L[a]<length>[/<Flength>/<Fcut>]\n"
"                 Define the wavenumber integration interval\n"
"                 dk=(2*PI)/(<length>*rmax). rmax is the maximum \n"
"                 epicentral distance. \n"
"                 There are 4 cases:\n"
"                 + (default) not set or set 0.0.\n"); printf(
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
"                 + (default) not set or set 0.0.\n"); printf(
"                   <vmin_ref> will be the minimum velocity\n"
"                   of model, but limited to %.1f. and if the \n", GRT_GREENFN_V_VMIN_REF); printf(
"                   depth gap between source and receiver is \n"
"                   thinner than %.1f km, PTAM will be appled\n", MIN_DEPTH_GAP_SRC_RCV); printf(
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
"                 <t0>: reference delay (s), default t0=0.0\n"); printf(
"                 <v0>: reference velocity (km/s), \n"
"                       default 0.0 not use.\n"); printf(
"\n"
"    -K<k0>[/<ampk>/<keps>]\n"
"                 Several parameters designed to define the\n"
"                 behavior in wavenumber integration. The upper\n"
"                 bound is \n"
"                 sqrt( (<k0>*mult)^2 + (<ampk>*w/<vmin_ref>)^2 ),\n"
"                 default mult=1.0.\n"
"                 <k0>:   designed to give residual k at\n"
"                         0 frequency, default is %.1f, and \n", GRT_GREENFN_K_K0); printf(
"                         multiply PI/hs in program, \n"
"                         where hs = max(fabs(depsrc-deprcv), %.1f).\n", MIN_DEPTH_GAP_SRC_RCV); printf(
"                 <ampk>: amplification factor, default is %.2f.\n", GRT_GREENFN_K_AMPK); printf(
"                 <keps>: a threshold for break wavenumber \n"
"                         integration in advance. See \n"
"                         (Yao and Harkrider, 1983) for details.\n"
"                         Default 0.0 not use.\n"); printf(
"\n"
"    -P<n>        Number of threads. Default use all cores.\n"
"\n"
"    -G<b1>[/<b2>/<b3>/<b4>]\n"
"                 Designed to choose which kind of source's Green's \n"
"                 functions will be computed, default is all (%d/%d/%d/%d). \n", 
(int)GRT_GREENFN_G_EX, (int)GRT_GREENFN_G_VF, (int)GRT_GREENFN_G_HF, (int)GRT_GREENFN_G_DC); printf(
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
"    grt greenfn -Mmilrow -N1000/0.01 -D2/0 -Ores -R2,4,6,8,10\n"
"\n\n\n"
);

}



/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    char* command = Ctrl->name;

    // 先为个别参数设置非0初始值
    Ctrl->N.zeta = GRT_GREENFN_N_ZETA;
    Ctrl->H.freq1 = GRT_GREENFN_H_FREQ1;
    Ctrl->H.freq2 = GRT_GREENFN_H_FREQ2;
    Ctrl->V.vmin_ref = GRT_GREENFN_V_VMIN_REF;
    Ctrl->K.k0 = GRT_GREENFN_K_K0;
    Ctrl->K.ampk = GRT_GREENFN_K_AMPK;
    Ctrl->G.doEX = GRT_GREENFN_G_EX;
    Ctrl->G.doVF = GRT_GREENFN_G_VF;
    Ctrl->G.doHF = GRT_GREENFN_G_HF;
    Ctrl->G.doDC = GRT_GREENFN_G_DC;

    int opt;
    while ((opt = getopt(argc, argv, ":M:D:N:O:H:L:V:E:K:R:S:P:G:esh")) != -1) {
        switch (opt) {
            // 模型路径，其中每行分别为 
            //      厚度(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  Qp   Qs
            // 互相用空格隔开即可
            case 'M':
                Ctrl->M.active = true;
                Ctrl->M.s_modelpath = strdup(optarg);
                Ctrl->M.s_modelname = get_basename(Ctrl->M.s_modelpath);
                break;

            // 震源和场点深度， -Ddepsrc/deprcv
            case 'D':
                Ctrl->D.active = true;
                Ctrl->D.s_depsrc = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                Ctrl->D.s_deprcv = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                if(2 != sscanf(optarg, "%[^/]/%s", Ctrl->D.s_depsrc, Ctrl->D.s_deprcv)){
                    GRTBadOptionError(command, D, "");
                };
                if(1 != sscanf(Ctrl->D.s_depsrc, "%lf", &Ctrl->D.depsrc)){
                    GRTBadOptionError(command, D, "");
                }
                if(1 != sscanf(Ctrl->D.s_deprcv, "%lf", &Ctrl->D.deprcv)){
                    GRTBadOptionError(command, D, "");
                }
                if(Ctrl->D.depsrc < 0.0 || Ctrl->D.deprcv < 0.0){
                    GRTBadOptionError(command, D, "Negative value in -D is not supported.");
                }
                break;

            // 点数,采样间隔,虚频率 -Nnt/dt/[zeta]
            case 'N':
                Ctrl->N.active = true;
                if(2 > sscanf(optarg, "%d/%lf/%lf", &Ctrl->N.nt, &Ctrl->N.dt, &Ctrl->N.zeta)){
                    GRTBadOptionError(command, N, "");
                };
                if(Ctrl->N.nt <= 0 || Ctrl->N.dt <= 0.0 || Ctrl->N.zeta <= 0.0){
                    GRTBadOptionError(command, N, "Nonpositive value in -N is not supported.");
                }
                break;

            // 输出路径 -Ooutput_dir
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.s_output_dir = strdup(optarg);
                break;

            // 频带 -H f1/f2
            case 'H':
                Ctrl->H.active = true;
                if(2 != sscanf(optarg, "%lf/%lf", &Ctrl->H.freq1, &Ctrl->H.freq2)){
                    GRTBadOptionError(command, H, "");
                };
                if(Ctrl->H.freq1>0.0 && Ctrl->H.freq2>0.0 && Ctrl->H.freq1 > Ctrl->H.freq2){
                    GRTBadOptionError(command, H, "Positive freq1 should be less than positive freq2.");
                }
                break;

            // 波数积分间隔 -L[a]<length>[/<Flength>/<Fcut>]
            case 'L':
                Ctrl->L.active = true;
                {
                    // 检查首字母是否为a，表明使用自适应Filon积分
                    int pos=0;
                    bool useSAFIM = false;
                    if(optarg[0] == 'a'){
                        pos++;
                        useSAFIM = true;
                    }
                    double filona = 0.0;
                    int n = sscanf(optarg+pos, "%lf/%lf/%lf", &Ctrl->L.Length, &filona, &Ctrl->L.filonCut);
                    if(n != 1 && n != 3){
                        GRTBadOptionError(command, L, "");
                    };
                    if(n == 1 && Ctrl->L.Length <= 0){
                        GRTBadOptionError(command, L, "Length should be positive.");
                    }
                    if(n == 3 && (filona <= 0 || Ctrl->L.filonCut < 0)){
                        GRTBadOptionError(command, L, "Flength/Ftol should be positive, Fcut should be nonnegative.");
                    }
                    if(n == 3){
                        useSAFIM ? (Ctrl->L.safilonTol = filona) : (Ctrl->L.filonLength = filona);
                    }
                }
                break;

            // 参考最小速度 -Vvmin_ref
            case 'V':
                Ctrl->V.active = true;
                if(0 == sscanf(optarg, "%lf", &Ctrl->V.vmin_ref)){
                    GRTBadOptionError(command, V, "");
                };
                break;

            // 时间延迟 -ET0/V0
            case 'E':
                Ctrl->E.active = true;
                if(0 == sscanf(optarg, "%lf/%lf", &Ctrl->E.delayT0, &Ctrl->E.delayV0)){
                    GRTBadOptionError(command, E, "");
                };
                if(Ctrl->E.delayV0 < 0.0){
                    GRTBadOptionError(command, E, "Can't set negative v0(%f) in -E.", Ctrl->E.delayV0);
                }
                break;

            // 波数积分相关变量 -Kk0/ampk/keps
            case 'K':
                Ctrl->K.active = true;
                if(0 == sscanf(optarg, "%lf/%lf/%lf", &Ctrl->K.k0, &Ctrl->K.ampk, &Ctrl->K.keps)){
                    GRTBadOptionError(command, K, "");
                };
                if(Ctrl->K.k0 < 0.0){
                    GRTBadOptionError(command, K, "Can't set negative k0(%f).", Ctrl->K.k0);
                }
                if(Ctrl->K.ampk < 0.0){
                    GRTBadOptionError(command, K, "Can't set negative ampk(%f).", Ctrl->K.ampk);
                }
                break;

            // 不打印在终端
            case 's':
                Ctrl->s.active = true;
                break;

            // 震中距数组，-Rr1,r2,r3,r4 ...
            case 'R':
                Ctrl->R.active = true;
                Ctrl->R.s_raw = strdup(optarg);
                Ctrl->R.s_rs = string_split(optarg, ",", &Ctrl->R.nr);
                // 转为浮点数
                Ctrl->R.rs = (MYREAL*)realloc(Ctrl->R.rs, sizeof(MYREAL)*(Ctrl->R.nr));
                for(MYINT i=0; i<Ctrl->R.nr; ++i){
                    Ctrl->R.rs[i] = atof(Ctrl->R.s_rs[i]);
                    if(Ctrl->R.rs[i] < 0.0){
                        GRTBadOptionError(command, R, "Can't set negative epicentral distance(%f).", Ctrl->R.rs[i]);
                    }
                }
                break;

            // 多线程数 -Pnthreads
            case 'P':
                Ctrl->P.active = true;
                if(1 != sscanf(optarg, "%d", &Ctrl->P.nthreads)){
                    GRTBadOptionError(command, P, "");
                };
                if(Ctrl->P.nthreads <= 0){
                    GRTBadOptionError(command, P, "Nonpositive value is not supported.");
                }
                set_num_threads(Ctrl->P.nthreads);
                break;

            // 选择要计算的格林函数 -G1/1/1/1
            case 'G': 
                Ctrl->G.active = true;
                Ctrl->G.doEX = Ctrl->G.doVF = Ctrl->G.doHF = Ctrl->G.doDC = false;
                {
                    int i1, i2, i3, i4;
                    i1 = i2 = i3 = i4 = 0;
                    if(0 == sscanf(optarg, "%d/%d/%d/%d", &i1, &i2, &i3, &i4)){
                        fprintf(stderr, "[%s] " BOLD_RED "Error in -G.\n" DEFAULT_RESTORE, command);
                        exit(EXIT_FAILURE);
                    };
                    Ctrl->G.doEX = (i1!=0);
                    Ctrl->G.doVF  = (i2!=0);
                    Ctrl->G.doHF  = (i3!=0);
                    Ctrl->G.doDC  = (i4!=0);
                }
                // 至少要有一个真
                if(!(Ctrl->G.doEX || Ctrl->G.doVF || Ctrl->G.doHF || Ctrl->G.doDC)){
                    GRTBadOptionError(command, G, "At least set one true value.");
                }
                break;

            // 输出波数积分中间文件， -Sidx1,idx2,idx3,...
            case 'S':
                Ctrl->S.active = true;
                Ctrl->S.s_raw = strdup(optarg);
                Ctrl->S.s_statsidxs = string_split(optarg, ",", &Ctrl->S.nstatsidxs);
                // 转为浮点数
                Ctrl->S.statsidxs = (MYINT*)realloc(Ctrl->S.statsidxs, sizeof(MYINT)*(Ctrl->S.nstatsidxs));
                for(MYINT i=0; i<Ctrl->S.nstatsidxs; ++i){
                    Ctrl->S.statsidxs[i] = atof(Ctrl->S.s_statsidxs[i]);
                }
                break;

            // 是否计算位移空间导数
            case 'e':
                Ctrl->e.active = true;
                break;
            
            GRT_Common_Options_in_Switch(command, (char)(optopt));
        }
    } // END get options

    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(command, argc > 1);
    GRTCheckOptionActive(command, Ctrl, M);
    GRTCheckOptionActive(command, Ctrl, D);
    GRTCheckOptionActive(command, Ctrl, N);
    GRTCheckOptionActive(command, Ctrl, R);
    GRTCheckOptionActive(command, Ctrl, O);

    // 建立保存目录
    GRTCheckMakeDir(command, Ctrl->O.s_output_dir);

    // 在目录中保留命令
    char *dummy = (char*)malloc(sizeof(char)*(strlen(Ctrl->O.s_output_dir)+100));
    sprintf(dummy, "%s/command", Ctrl->O.s_output_dir);
    FILE *fp = GRTCheckOpenFile(command, dummy, "a");
    fprintf(fp, GRT_MAIN_COMMAND " ");  // 主程序名
    for(int i=0; i<argc; ++i){
        fprintf(fp, "%s ", argv[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
    free(dummy);
}


/**
 * 将某一道做ifft，做时间域处理，保存到sac文件
 * 
 * @param     delay     时间延迟
 * @param     mult      幅值放大系数
 * @param     nt        点数
 * @param     dt        采样间隔
 * @param     nf        频率点数
 * @param     df        频率间隔
 * @param     wI        虚频率
 * @param     grncplx   复数形式的格林函数频谱
 * @param     fftw_grn  将频谱写到FFTW_COMPLEX类型中
 * @param     out       ifft后的时域数据
 * @param     float_arr 将时域数据写到float类型的数组中
 * @param     plan      FFTW_PLAN
 * @param     hd        SAC头段变量结构体
 * @param     outpath   sac文件保存路径
 */
static void ifft_one_trace(
    MYREAL delay, MYREAL mult, MYINT nt, MYREAL dt, MYINT nf, MYREAL df, MYREAL wI,
    MYCOMPLEX *grncplx, _FFTW_COMPLEX *fftw_grn, MYREAL *out, float *float_arr,
    _FFTW_PLAN plan, SACHEAD *hd, const char *outpath)
{
    // 赋值复数，包括时移
    MYCOMPLEX cfac, ccoef;
    cfac = exp(I*PI2*df*delay);
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
}





/**
 * 将一条数据反变换回时间域再进行处理，保存到SAC文件
 * 
 * @param     Ctrl          参数控制
 * @param     srcname       震源类型
 * @param     delayT        延迟时间
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
    const GRT_MODULE_CTRL *Ctrl, const char *srcname, const char ch, MYREAL delayT,
    SACHEAD *hd, char *s_outpath, const char *s_output_subdir, const char *s_prefix,
    const int sgn, MYCOMPLEX *grncplx, fftw_complex *fftw_grn, MYREAL *out, float *float_arr, fftw_plan plan)
{
    char kcmpnm[9];
    snprintf(kcmpnm, sizeof(kcmpnm), "%s%s%c", s_prefix, srcname, ch);
    strcpy(hd->kcmpnm, kcmpnm);
    sprintf(s_outpath, "%s/%s.sac", s_output_subdir, kcmpnm);
    ifft_one_trace(
        delayT, sgn, 
        Ctrl->N.nt, Ctrl->N.dt, Ctrl->N.nf, Ctrl->N.df, Ctrl->N.wI,
        grncplx, fftw_grn, out, float_arr, plan, hd, s_outpath);
}


/** 子模块主函数 */
int greenfn_main(int argc, char **argv) {
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    Ctrl->name = strdup(argv[0]);
    const char *command = Ctrl->name;

    // 传入参数 
    getopt_from_command(Ctrl, argc, argv);

    // 读入模型文件
    if((Ctrl->M.pymod = read_pymod_from_file(command, Ctrl->M.s_modelpath, Ctrl->D.depsrc, Ctrl->D.deprcv, true)) == NULL){
        exit(EXIT_FAILURE);
    }
    PYMODEL1D *pymod = Ctrl->M.pymod;

    // 当震源位于液体层中时，仅允许计算爆炸源对应的格林函数
    // 程序结束前会输出对应警告
    if(pymod->Vb[pymod->isrc]==0.0){
        Ctrl->G.doHF = Ctrl->G.doVF = Ctrl->G.doDC = false;
    }

    // 最大最小速度
    MYREAL vmin, vmax;
    get_pymod_vmin_vmax(pymod, &vmin, &vmax);

    // 参考最小速度
    if(!Ctrl->V.active){
        Ctrl->V.vmin_ref = GRT_MAX(vmin, GRT_GREENFN_V_VMIN_REF);
    } 

    // 如果没有主动设置vmin_ref，则判断是否要自动使用PTAM
    if( !Ctrl->V.active && fabs(Ctrl->D.deprcv - Ctrl->D.depsrc) <= MIN_DEPTH_GAP_SRC_RCV) {
        Ctrl->V.vmin_ref = - fabs(Ctrl->V.vmin_ref);
    }

    // 时窗长度 
    Ctrl->N.winT = Ctrl->N.nt*Ctrl->N.dt;

    // 最大震中距
    MYREAL rmax = Ctrl->R.rs[findMinMax_MYREAL(Ctrl->R.rs, Ctrl->R.nr, true)];   

    // 时窗最大截止时刻
    MYREAL tmax = Ctrl->E.delayT0 + Ctrl->N.winT;
    if(Ctrl->E.delayV0 > 0.0)   tmax += rmax/Ctrl->E.delayV0;

    // 自动选择积分间隔，默认使用传统离散波数积分
    // 自动选择会给出很保守的值（较大的Length）
    if(Ctrl->L.Length == 0.0){
        Ctrl->L.Length = 15.0; 
        double jus = GRT_SQUARE(vmax*tmax) - GRT_SQUARE(Ctrl->D.deprcv - Ctrl->D.depsrc);
        if(jus >= 0.0){
            Ctrl->L.Length = GRT_MAX(1.0 + sqrt(jus)/rmax + 0.5, Ctrl->L.Length); // +0.5为保守值
        }
    }

    // 虚频率
    Ctrl->N.wI = Ctrl->N.zeta*PI/Ctrl->N.winT;

    // 定义要计算的频率、时窗等
    Ctrl->N.nf = Ctrl->N.nt/2 + 1;
    Ctrl->N.df = 1.0/Ctrl->N.winT;
    Ctrl->N.freqs = (MYREAL*)malloc(Ctrl->N.nf*sizeof(MYREAL));
    for(int i=0; i<Ctrl->N.nf; ++i){
        Ctrl->N.freqs[i] = i*Ctrl->N.df;
    }

    // 自定义频段
    Ctrl->H.nf1 = 0; Ctrl->H.nf2 = Ctrl->N.nf-1;
    if(Ctrl->H.freq1 > 0.0){
        Ctrl->H.nf1 = GRT_MIN(ceil(Ctrl->H.freq1/Ctrl->N.df), Ctrl->N.nf-1);
    }
    if(Ctrl->H.freq2 > 0.0){
        Ctrl->H.nf2 = GRT_MIN(floor(Ctrl->H.freq2/Ctrl->N.df), Ctrl->N.nf-1);
    }
    Ctrl->H.nf2 = GRT_MAX(Ctrl->H.nf1, Ctrl->H.nf2);

    // 波数积分中间文件输出目录
    if(Ctrl->S.nstatsidxs > 0){
        Ctrl->S.s_statsdir = (char*)malloc(sizeof(char)*(strlen(Ctrl->M.s_modelpath)+strlen(Ctrl->O.s_output_dir)+strlen(Ctrl->D.s_depsrc)+strlen(Ctrl->D.s_deprcv)+100));
        sprintf(Ctrl->S.s_statsdir, "%s_grtstats", Ctrl->O.s_output_dir);
        
        // 建立保存目录
        GRTCheckMakeDir(command, Ctrl->S.s_statsdir);
        sprintf(Ctrl->S.s_statsdir, "%s/%s_%s_%s", Ctrl->S.s_statsdir, Ctrl->M.s_modelname, Ctrl->D.s_depsrc, Ctrl->D.s_deprcv);
        GRTCheckMakeDir(command, Ctrl->S.s_statsdir);
    }

    // 建立格林函数的complex数组
    MYCOMPLEX *(*grn)[SRC_M_NUM][CHANNEL_NUM] = (MYCOMPLEX*(*)[SRC_M_NUM][CHANNEL_NUM]) calloc(Ctrl->R.nr, sizeof(*grn));
    MYCOMPLEX *(*grn_uiz)[SRC_M_NUM][CHANNEL_NUM] = (Ctrl->e.active)? (MYCOMPLEX*(*)[SRC_M_NUM][CHANNEL_NUM]) calloc(Ctrl->R.nr, sizeof(*grn_uiz)) : NULL;
    MYCOMPLEX *(*grn_uir)[SRC_M_NUM][CHANNEL_NUM] = (Ctrl->e.active)? (MYCOMPLEX*(*)[SRC_M_NUM][CHANNEL_NUM]) calloc(Ctrl->R.nr, sizeof(*grn_uir)) : NULL;

    for(int ir=0; ir<Ctrl->R.nr; ++ir){
        for(int i=0; i<SRC_M_NUM; ++i){
            for(int c=0; c<CHANNEL_NUM; ++c){
                grn[ir][i][c] = (MYCOMPLEX*)calloc(Ctrl->N.nf, sizeof(MYCOMPLEX));
                if(grn_uiz)  grn_uiz[ir][i][c] = (MYCOMPLEX*)calloc(Ctrl->N.nf, sizeof(MYCOMPLEX));
                if(grn_uir)  grn_uir[ir][i][c] = (MYCOMPLEX*)calloc(Ctrl->N.nf, sizeof(MYCOMPLEX));
            }
        }
    }


    // 在计算前打印所有参数
    if(! Ctrl->s.active){
        print_Ctrl(Ctrl);
    }
    

    //==============================================================================
    // 计算格林函数
    integ_grn_spec(
        pymod, Ctrl->H.nf1, Ctrl->H.nf2, Ctrl->N.freqs, Ctrl->R.nr, Ctrl->R.rs, Ctrl->N.wI,
        Ctrl->V.vmin_ref, Ctrl->K.keps, Ctrl->K.ampk, Ctrl->K.k0, Ctrl->L.Length, Ctrl->L.filonLength, Ctrl->L.safilonTol, Ctrl->L.filonCut, !Ctrl->s.active,
        grn, Ctrl->e.active, grn_uiz, grn_uir,
        Ctrl->S.s_statsdir, Ctrl->S.nstatsidxs, Ctrl->S.statsidxs
    );
    //==============================================================================
    

    // 使用fftw3做反傅里叶变换
    // 分配fftw_complex内存
    _FFTW_COMPLEX *fftw_grn = (_FFTW_COMPLEX*)_FFTW_MALLOC(sizeof(_FFTW_COMPLEX)*Ctrl->N.nf);
    MYREAL *out = (MYREAL*)malloc(sizeof(MYREAL)*Ctrl->N.nt);
    float *float_arr = (float*)malloc(sizeof(float)*Ctrl->N.nt);

    // fftw计划
    _FFTW_PLAN plan = _FFTW_PLAN_DFT_C2R_1D(Ctrl->N.nt, fftw_grn, out, FFTW_ESTIMATE);
    
    // 建立SAC头文件，包含必要的头变量
    SACHEAD hd = new_sac_head(Ctrl->N.dt, Ctrl->N.nt, Ctrl->E.delayT0);
    // 发震时刻作为参考时刻
    hd.o = 0.0; 
    hd.iztype = IO; 
    // 记录震源和台站深度
    hd.evdp = Ctrl->D.depsrc; // km
    hd.stel = (-1.0)*Ctrl->D.deprcv*1e3; // m
    // 写入虚频率
    hd.user0 = Ctrl->N.wI;
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

    
    // 下面计算的同时也打印走时
    if( ! Ctrl->s.active){
        printf("\n\n");
        printf("------------------------------------------------\n");
        printf(" Distance(km)     Tp(secs)         Ts(secs)     \n");
    }
    
    // 做反傅里叶变换，保存SAC文件
    for(int ir=0; ir<Ctrl->R.nr; ++ir){
        hd.dist = Ctrl->R.rs[ir];

        // 文件保存子目录
        char *s_output_subdir = (char*)malloc(sizeof(char)*(
            strlen(Ctrl->O.s_output_dir)+strlen(Ctrl->M.s_modelpath)+
            strlen(Ctrl->D.s_depsrc)+strlen(Ctrl->D.s_deprcv)+strlen(Ctrl->R.s_rs[ir])+100));
        
        sprintf(s_output_subdir, "%s/%s_%s_%s_%s", Ctrl->O.s_output_dir, Ctrl->M.s_modelname, Ctrl->D.s_depsrc, Ctrl->D.s_deprcv, Ctrl->R.s_rs[ir]);
        GRTCheckMakeDir(command, s_output_subdir);

        // 时间延迟 
        MYREAL delayT = Ctrl->E.delayT0;
        if(Ctrl->E.delayV0 > 0.0)   delayT += sqrt( GRT_SQUARE(Ctrl->R.rs[ir]) + GRT_SQUARE(Ctrl->D.deprcv - Ctrl->D.depsrc) ) / Ctrl->E.delayV0;
        // 修改SAC头段时间变量
        hd.b = delayT;

        // 计算理论走时
        hd.t0 = compute_travt1d(pymod->Thk, pymod->Va, pymod->n, pymod->isrc, pymod->ircv, Ctrl->R.rs[ir]);
        strcpy(hd.kt0, "P");
        hd.t1 = compute_travt1d(pymod->Thk, pymod->Vb, pymod->n, pymod->isrc, pymod->ircv, Ctrl->R.rs[ir]);
        strcpy(hd.kt1, "S");

        for(int im=0; im<SRC_M_NUM; ++im){
            if(!Ctrl->G.doEX  && im==0)  continue;
            if(!Ctrl->G.doVF  && im==1)  continue;
            if(!Ctrl->G.doHF  && im==2)  continue;
            if(!Ctrl->G.doDC  && im>=3)  continue;

            int modr = SRC_M_ORDERS[im];
            int sgn=1;  // 用于反转Z分量
            for(int c=0; c<CHANNEL_NUM; ++c){
                if(modr==0 && ZRTchs[c]=='T')  continue;  // 跳过输出0阶的T分量

                // 文件保存总路径
                char *s_outpath = (char*)malloc(sizeof(char)*(strlen(s_output_subdir)+100));
                char s_prefix[] = "";

                // Z分量反转
                sgn = (ZRTchs[c]=='Z') ? -1 : 1;

                write_one_to_sac(Ctrl, SRC_M_NAME_ABBR[im], ZRTchs[c], delayT, &hd, s_outpath, s_output_subdir, s_prefix, sgn, grn[ir][im][c], fftw_grn, out, float_arr, plan);
                if(Ctrl->e.active){
                    write_one_to_sac(Ctrl, SRC_M_NAME_ABBR[im], ZRTchs[c], delayT, &hd, s_outpath, s_output_subdir, "z", sgn*(-1), grn_uiz[ir][im][c], fftw_grn, out, float_arr, plan);
                    write_one_to_sac(Ctrl, SRC_M_NAME_ABBR[im], ZRTchs[c], delayT, &hd, s_outpath, s_output_subdir, "r", sgn, grn_uir[ir][im][c], fftw_grn, out, float_arr, plan);
                }

                free(s_outpath);
            }
        }


        if( ! Ctrl->s.active){
            printf(" %-15s  %-15.3f  %-15.3f\n", Ctrl->R.s_rs[ir], hd.t0, hd.t1);
        }

        free(s_output_subdir);
    } // End distances loop

    if( ! Ctrl->s.active){
        printf("------------------------------------------------\n");
        printf("\n");
    }

    // 输出警告：当震源位于液体层中时，仅允许计算爆炸源对应的格林函数
    if(pymod->Vb[pymod->isrc]==0.0){
        fprintf(stderr, "[%s] " BOLD_YELLOW 
            "The source is located in the liquid layer, "
            "therefore only the Green's Funtions for the Explosion source will be computed.\n" 
            DEFAULT_RESTORE, command);
    }
    

    // 释放内存
    for(int ir=0; ir<Ctrl->R.nr; ++ir){
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

    _FFTW_FREE(fftw_grn);
    free(out);
    free(float_arr);
    _FFTW_DESTROY_PLAN(plan);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}

