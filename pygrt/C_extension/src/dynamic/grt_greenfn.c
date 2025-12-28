/**
 * @file   grt_greenfn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-11-28
 * 
 *    定义main函数，形成命令行式的用法（不使用python的entry_points，会牺牲性能）
 *    计算不同震源的格林函数
 * 
 */


#include "grt/dynamic/grn.h"
#include "grt/dynamic/signals.h"
#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/common/search.h"
#include "grt/common/sacio2.h"
#include "grt/common/util.h"
#include "grt/common/myfftw.h"
#include "grt/common/travt.h"
#include "grt/integral/integ_method.h"
#include "grt.h"


// 一些变量的非零默认值
#define GRT_GREENFN_N_ZETA        0.8
#define GRT_GREENFN_N_UPSAMPLE    1
#define GRT_GREENFN_H_FREQ1      -1.0
#define GRT_GREENFN_H_FREQ2      -1.0
#define GRT_GREENFN_K_VMIN        0.1
#define GRT_GREENFN_K_K0          5.0
#define GRT_GREENFN_K_AMPK       1.15
#define GRT_GREENFN_G_EX       true
#define GRT_GREENFN_G_VF       true
#define GRT_GREENFN_G_HF       true
#define GRT_GREENFN_G_DC       true



/** 该子模块的参数控制结构体 */
typedef struct {
    /** 输入模型 */
    struct {
        bool active;
        char *s_modelpath;        ///< 模型路径
        const char *s_modelname;  ///< 模型名称
        GRT_MODEL1D *mod1d;         ///< 模型结构体指针
    } M;
    /** 震源和接收器深度 */
    struct {
        bool active;
        real_t depsrc;
        real_t deprcv;
        char *s_depsrc;
        char *s_deprcv;
    } D;
    /** 波形时窗 */
    struct {
        bool active;
        size_t nt;
        size_t nf;
        real_t dt;
        real_t df;
        real_t winT;  ///< 时窗长度 
        real_t zeta;  ///< 虚频率系数， w <- w - zeta*PI/r* 1j
        real_t wI;    ///< 虚频率  zeta*PI/r
        real_t *freqs;
        size_t upsample_n;  ///< 升采样倍数
        bool keepAllFreq;
    } N;
    /** 输出目录 */
    struct {
        bool active;
        char *s_output_dir;
    } O;
    /** 频段 */
    struct {
        bool active;
        real_t freq1;
        real_t freq2;
        size_t nf1;
        size_t nf2;
    } H;
    /** 波数积分间隔以及方法 */
    struct {
        bool active;
        real_t Length;
        real_t kcut;
        struct {
            bool active;
            real_t Length;
        } FIM;
        struct {
            bool active;
            real_t tol;
        } SAFIM;
    } L;
    /** 波数积分收敛方法 */
    struct {
        bool active;
        bool applyDCM;
        bool applyPTAM;
        bool applyNoConverg;
    } C;
    /** 波数积分上限 */
    struct {
        bool active;
        real_t keps;
        real_t ampk;
        real_t k0;
        real_t vmin;
    } K;
    /** 时间延迟 */
    struct {
        bool active;
        real_t delayT0;
        real_t delayV0;
    } E;
    /** 波数积分过程的核函数文件 */
    struct {
        bool active;
        char *s_raw;
        char **s_statsidxs;
        size_t *statsidxs;
        size_t nstatsidxs;
        char *s_statsdir;  ///< 保存目录，和SAC文件目录同级
    } S;
    /** 震中距 */
    struct {
        bool active;
        char *s_raw;
        char **s_rs;
        real_t *rs;
        size_t nr;
    } R;
    /** 多线程 */
    struct {
        bool active;
        int nthreads; ///< 线程数
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
    // M
    GRT_SAFE_FREE_PTR(Ctrl->M.s_modelpath);
    grt_free_mod1d(Ctrl->M.mod1d);
    
    // D
    GRT_SAFE_FREE_PTR(Ctrl->D.s_depsrc);
    GRT_SAFE_FREE_PTR(Ctrl->D.s_deprcv);

    // N
    GRT_SAFE_FREE_PTR(Ctrl->N.freqs);

    // O
    GRT_SAFE_FREE_PTR(Ctrl->O.s_output_dir);

    // R
    GRT_SAFE_FREE_PTR_ARRAY(Ctrl->R.s_rs, Ctrl->R.nr);
    GRT_SAFE_FREE_PTR(Ctrl->R.s_raw);
    GRT_SAFE_FREE_PTR(Ctrl->R.rs);

    // S
    if(Ctrl->S.active){
        GRT_SAFE_FREE_PTR(Ctrl->S.s_raw);
        GRT_SAFE_FREE_PTR_ARRAY(Ctrl->S.s_statsidxs, Ctrl->S.nstatsidxs);
        GRT_SAFE_FREE_PTR(Ctrl->S.statsidxs);
        GRT_SAFE_FREE_PTR(Ctrl->S.s_statsdir);
    }

    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印结构体中的参数 */
static void print_Ctrl(const GRT_MODULE_CTRL *Ctrl){
    grt_print_mod1d(Ctrl->M.mod1d);

    const char format[]      = "   \%-20s  \%s\n";
    const char format_real[] = "   \%-20s  \%.3f\n";
    const char format_size[]  = "   \%-20s  \%zu\n";
    char line[100];
    printf("------------------------------------------------\n");
    printf(format, "PARAMETER", "VALUE");
    printf(format, "model_path", Ctrl->M.s_modelpath);
    printf(format_real, "k0", Ctrl->K.k0);
    printf(format_real, "ampk", Ctrl->K.ampk);
    printf(format_real, "keps", Ctrl->K.keps);
    printf(format_real, "vmin", Ctrl->K.vmin);
    printf(format_real, "Length", Ctrl->L.Length);
    printf(format_real, "kcut",   Ctrl->L.kcut);
    printf(format_real, "filonLength",   Ctrl->L.FIM.Length);
    printf(format_real, "safilonTol",   Ctrl->L.SAFIM.tol);
    printf(format, "applyDCM", (Ctrl->C.applyDCM)? "true" : "false");
    printf(format, "applyPTAM", (Ctrl->C.applyPTAM)? "true" : "false");

    printf(format_size, "nt", Ctrl->N.nt);
    printf(format_real, "dt", Ctrl->N.dt);
    printf(format_real, "winT", Ctrl->N.winT);
    printf(format_real, "zeta", Ctrl->N.zeta);
    printf(format_real, "delayT0", Ctrl->E.delayT0);
    printf(format_real, "delayV0", Ctrl->E.delayV0);
    printf(format_real, "tmax", Ctrl->E.delayT0 + Ctrl->N.winT);
    
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
"  at a close or same depth, Direct Convergence Method (DCM)\n"
"  (Zhu et al., 2026) will be applied automatically.\n"
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
"    grt greenfn -M<model> -D<depsrc>/<deprcv> \n"
"        -N<nt>/<dt>[+w<zeta>][+n<fac>][+a] \n"
"        -R<r1>,<r2>[,...]     -O<outdir>     [-H<f1>/<f2>] \n"
"        [-L<length>] [-C[d|p|n]] [-E<t0>[/<v0>]] \n" 
"        [-K[+k<k0>][+s<ampk>][+e<keps>][+v<vmin>]]\n"
"        [-P<nthreads>] [-Ge|v|h|s] \n"
"        [-S[<i1>,<i2>,...]] [-e] [-s]\n"
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
"    -N<nt>/<dt>[+w<zeta>][+n<fac>][+a] \n"
"                 <nt>:   number of points. (NOT requires 2^n).\n"
"                 <dt>:   time interval (secs). \n"
"                 +w<zeta>: define the coefficient of imaginary \n"
"                           frequency wI=zeta*PI/T, where T=nt*dt.\n"
"                           Default zeta=%.1f.\n", GRT_GREENFN_N_ZETA); printf(
"                 +n<fac>:  upsampling factor (integer)\n"
"                           i.e.  nt <-- nt * <fac>\n"
"                                 dt <-- dt / <fac>\n"
"                           and calculated frequencies stay unchanged.\n"
"                 +a:       All frequencies are calculated regardless of\n"
"                           how low the frequency is.\n"
"\n"
"    -R<r1>,<r2>[,...]\n"
"                 Multiple epicentral distance (km), \n"
"                 seperated by comma.\n"
"\n"
"    -O<outdir>   Directorypath of output for saving.\n"
"\n"
"    -H<f1>/<f2>  Apply bandpass filer with rectangle window, \n"
"                 default no filter.\n"
"                 <f1>: lower frequency (Hz), %.1f means low pass.\n", GRT_GREENFN_H_FREQ1); printf(
"                 <f2>: upper frequency (Hz), %.1f means high pass.\n", GRT_GREENFN_H_FREQ2); printf(
"\n"
"    -L[<length>][+l<Flength>][+a<Ftol>][+o<offset>]\n"
"                 Define the wavenumber integration interval\n"
"                 dk=(2*PI)/(<length>*rmax) and methods. \n"
"                 rmax is the maximum epicentral distance. \n"
"                 For DWM:\n"
"                 + (default) not set or set 0.0. \n"
"                   <length> will be determined automatically\n"
"                   in program with the criterion (Bouchon, 1980).\n"
"                 + manually set one POSITIVE <length>, e.g. -L20\n"
"                 For FIM or SAFIM:\n"
"                 + +l<Flength> defines the dk of the FIM.\n"
"                 + +a<Ftol> defines the tolerance of the SAFIM.\n"
"                   you can't set both.\n"
"                 + +o<offset> split the integration into two parts,\n"
"                   [0, k*] and [k*, kmax], in which k*=<offset>/rmax,\n"
"                   the former uses DWM and the latter uses FIM/SAFIM.\n"
"\n"
"    -Cd|p|n      Set convergence method,\n"
"                 + d: Direct Convergence Method (DCM).\n"
"                 + p: Peak-Trough Averaging Method (PTAM).\n"
"                 + n: None.\n"
"                 Default use +cd when fabs(depsrc-deprcv) <= %.1f.\n", GRT_MIN_DEPTH_GAP_SRC_RCV); printf(
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
"    -K[+k<k0>][+s<ampk>][+e<keps>][+v<vmin>]\n"
"                 Several parameters designed to define the\n"
"                 behavior in wavenumber integration. The upper\n"
"                 bound is \n"
"                 sqrt( <k0>^2 + (<ampk>*w/<vmin_ref>)^2 ),\n"
"                 <k0>:   designed to give residual k at\n"
"                         0 frequency, default is %.1f, and \n", GRT_GREENFN_K_K0); printf(
"                         multiply PI/hs in program, \n"
"                         where hs = max(fabs(depsrc-deprcv), %.1f).\n", GRT_MIN_DEPTH_GAP_SRC_RCV); printf(
"                 <ampk>: amplification factor, default is %.2f.\n", GRT_GREENFN_K_AMPK); printf(
"                 <keps>: a threshold for break wavenumber \n"
"                         integration in advance. See \n"
"                         (Yao and Harkrider, 1983) for details.\n"
"                         Default 0.0 not use.\n"
"                 <vmin>: Minimum velocity (km/s) for reference. This\n"
"                         is designed to define the upper bound \n"
"                         of wavenumber integration.\n"
"                         There are 2 cases:\n"
"                         + (default) not set or set 0.0.\n"); printf(
"                           <vmin> will be the minimum velocity\n"
"                           of model, but limited to %.1f.\n", GRT_GREENFN_K_VMIN); printf(
"                         + manually set POSITIVE value. \n"
"\n"
"    -P<n>        Number of threads. Default use all cores.\n"
"\n"
"    -Ge|v|h|s\n"
"                 Designed to choose which kind of source's Green's \n"
"                 functions will be computed, default is all (-Gevhs). \n"); printf(
"                 Four bool type (0 or 1) options are\n"
"                 <b1>: Explosion (EX)\n"
"                 <b2>: Vertical Force (VF)\n"
"                 <b3>: Horizontal Force (HF)\n"
"                 <b4>: Shear (DC)\n"
"\n"
"    -S[<i1>,<i2>,...]\n"
"                 Frequency (index) of statsfile in wavenumber\n"
"                 integration to be output, require 0 <= i <= nf-1,\n"
"                 where nf=nt/2+1. These option is designed to check\n"
"                 the trend of kernel with wavenumber.\n"
"                 Empty -S means all frequency index.\n"
"\n"
"    -e           Compute the spatial derivatives, ui_z and ui_r,\n"
"                 of displacement u. In filenames, prefix \"r\" means \n"
"                 ui_r and \"z\" means ui_z. The units of derivatives\n"
"                 for different sources are: \n"
"                 + Explosion:     1e-25 /(dyne-cm)\n"
"                 + Single Force:  1e-20 /(dyne)\n"
"                 + Shear:         1e-25 /(dyne-cm)\n"
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
    // 先为个别参数设置非0初始值
    Ctrl->N.zeta = GRT_GREENFN_N_ZETA;
    Ctrl->N.upsample_n = GRT_GREENFN_N_UPSAMPLE;
    Ctrl->H.freq1 = GRT_GREENFN_H_FREQ1;
    Ctrl->H.freq2 = GRT_GREENFN_H_FREQ2;
    Ctrl->K.k0 = GRT_GREENFN_K_K0;
    Ctrl->K.ampk = GRT_GREENFN_K_AMPK;
    Ctrl->G.doEX = GRT_GREENFN_G_EX;
    Ctrl->G.doVF = GRT_GREENFN_G_VF;
    Ctrl->G.doHF = GRT_GREENFN_G_HF;
    Ctrl->G.doDC = GRT_GREENFN_G_DC;

    int opt;
    while ((opt = getopt(argc, argv, ":M:D:N:O:H:L:C:E:K:R:S::P:G:esh")) != -1) {
        switch (opt) {
            // 模型路径，其中每行分别为 
            //      厚度(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  Qp   Qs
            // 互相用空格隔开即可
            case 'M':
                Ctrl->M.active = true;
                Ctrl->M.s_modelpath = strdup(optarg);
                Ctrl->M.s_modelname = grt_get_basename(Ctrl->M.s_modelpath);
                break;

            // 震源和场点深度， -Ddepsrc/deprcv
            case 'D':
                Ctrl->D.active = true;
                Ctrl->D.s_depsrc = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                Ctrl->D.s_deprcv = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                if(2 != sscanf(optarg, "%[^/]/%s", Ctrl->D.s_depsrc, Ctrl->D.s_deprcv)){
                    GRTBadOptionError(D, "");
                };
                if(1 != sscanf(Ctrl->D.s_depsrc, "%lf", &Ctrl->D.depsrc)){
                    GRTBadOptionError(D, "");
                }
                if(1 != sscanf(Ctrl->D.s_deprcv, "%lf", &Ctrl->D.deprcv)){
                    GRTBadOptionError(D, "");
                }
                if(Ctrl->D.depsrc < 0.0 || Ctrl->D.deprcv < 0.0){
                    GRTBadOptionError(D, "Negative value in -D is not supported.");
                }
                break;

            // 点数,采样间隔,虚频率 -Nnt/dt[+w<zeta>][+n<scale>][+a]
            case 'N':
                Ctrl->N.active = true;
                {
                    char *string = strdup(optarg);
                    char *token = strtok(string, "+");
                    if(2 != sscanf(token, "%zu/%lf", &Ctrl->N.nt, &Ctrl->N.dt)){
                        GRTBadOptionError(N, "");
                    };
                    if(Ctrl->N.nt <= 0 || Ctrl->N.dt <= 0.0){
                        GRTBadOptionError(N, "Nonpositive value in -N is not supported.");
                    }

                    // 处理 + 号指令
                    token = strtok(NULL, "+");
                    while(token != NULL){
                        switch (token[0]){
                            case 'n':
                                if(1 != sscanf(token+1, "%zu", &Ctrl->N.upsample_n)){
                                    GRTBadOptionError(N, "");
                                }
                                if(Ctrl->N.upsample_n <= 0){
                                    GRTBadOptionError(N, "+%s need positive integer, but get (%zu).", token, Ctrl->N.upsample_n);
                                }
                                break;

                            case 'w':
                                if(1 != sscanf(token+1, "%lf", &Ctrl->N.zeta)){
                                    GRTBadOptionError(N, "");
                                }
                                if(Ctrl->N.zeta <= 0.0){
                                    GRTBadOptionError(N, "+%s need positive float, but get (%lf).", token, Ctrl->N.zeta);
                                }
                                break;

                            case 'a':
                                Ctrl->N.keepAllFreq = true;
                                break;
                            
                            default:
                                GRTBadOptionError(N, "+%s is not supported.", token);
                                break;
                        }

                        token = strtok(NULL, "+");
                    }

                    GRT_SAFE_FREE_PTR(string);
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
                    GRTBadOptionError(H, "");
                };
                if(Ctrl->H.freq1>0.0 && Ctrl->H.freq2>0.0 && Ctrl->H.freq1 > Ctrl->H.freq2){
                    GRTBadOptionError(H, "Positive freq1 should be less than positive freq2.");
                }
                break;

            // 波数积分间隔 -L[<length>][+l<Flength>][+a<Ftol>][+o<offset>]
            case 'L':
                Ctrl->L.active = true;
                {
                    char *string = strdup(optarg);
                    char *token = strtok(string, "+");
                    // 如果首先不是加号，则先读取DWM的length
                    if(optarg[0] != '+'){
                        if(1 != sscanf(optarg, "%lf", &Ctrl->L.Length)){
                            GRTBadOptionError(L, "");
                        }
                        token = strtok(NULL, "+");
                    }

                    while(token != NULL){
                        switch(token[0]) {
                            case 'l':
                                Ctrl->L.FIM.active = true;
                                if(1 != sscanf(token+1, "%lf", &Ctrl->L.FIM.Length)){
                                    GRTBadOptionError(L+l, "");
                                }
                                break;
                            
                            case 'a':
                                Ctrl->L.SAFIM.active = true;
                                if(1 != sscanf(token+1, "%lf", &Ctrl->L.SAFIM.tol)){
                                    GRTBadOptionError(L+a, "");
                                }
                                break;
                            
                            case 'o':
                                if(1 != sscanf(token+1, "%lf", &Ctrl->L.kcut)){
                                    GRTBadOptionError(L+o, "");
                                }
                                break;

                            default:
                                GRTBadOptionError(L, "-L+%s is not supported.", token);
                                break;
                        }

                        token = strtok(NULL, "+");
                    }

                    if(Ctrl->L.FIM.active && Ctrl->L.SAFIM.active){
                        GRTBadOptionError(L, "You can't set -L+a and -L+l both.");
                    }
                    
                    GRT_SAFE_FREE_PTR(string);
                }
                break;

            // 波数积分收敛方法  -Cd|p|n
            case 'C':
                Ctrl->C.active = true;
                if(strlen(optarg) == 0){
                    GRTBadOptionError(C, "");
                }
                switch (optarg[0]){
                    case 'p':
                        Ctrl->C.applyPTAM = true;
                        break;
                    case 'd':
                        Ctrl->C.applyDCM = true;
                        break;
                    case 'n':
                        Ctrl->C.applyNoConverg = true;
                        break;
                    default:
                        GRTBadOptionError(C, "-C+%s is not supported.", optarg);
                        break;
                }
                if(Ctrl->C.applyPTAM && Ctrl->C.applyDCM){
                    GRTBadOptionError(C, "You can't set -Cd and -Cp both.");
                }
                break;

            // 时间延迟 -ET0/V0
            case 'E':
                Ctrl->E.active = true;
                if(0 == sscanf(optarg, "%lf/%lf", &Ctrl->E.delayT0, &Ctrl->E.delayV0)){
                    GRTBadOptionError(E, "");
                };
                if(Ctrl->E.delayV0 < 0.0){
                    GRTBadOptionError(E, "Can't set negative v0(%f) in -E.", Ctrl->E.delayV0);
                }
                break;

            // 波数积分相关变量 -K[+k<k0>][+s<ampk>][+e<keps>][+v<vmin>]
            case 'K':
                Ctrl->K.active = true;
                {
                char *line = strdup(optarg);
                char *token = strtok(line, "+");
                while(token != NULL){
                    switch(token[0]) {
                        case 'k':
                            if(1 != sscanf(token+1, "%lf", &Ctrl->K.k0)){
                                GRTBadOptionError(K+k, "");
                            }
                            if(Ctrl->K.k0 < 0.0){
                                GRTBadOptionError(K, "Can't set negative k0(%f).", Ctrl->K.k0);
                            }
                            break;

                        case 's':
                            if(1 != sscanf(token+1, "%lf", &Ctrl->K.ampk)){
                                GRTBadOptionError(K+s, "");
                            }
                            if(Ctrl->K.ampk < 0.0){
                                GRTBadOptionError(K, "Can't set negative ampk(%f).", Ctrl->K.ampk);
                            }
                            break;

                        case 'e':
                            if(1 != sscanf(token+1, "%lf", &Ctrl->K.keps)){
                                GRTBadOptionError(K+e, "");
                            }
                            break;

                        case 'v':
                            if(1 != sscanf(token+1, "%lf", &Ctrl->K.vmin)){
                                GRTBadOptionError(K+v, "");
                            }
                            if(Ctrl->K.vmin <= 0.0){
                                GRTBadOptionError(K, "Illegal vmin (%f).", Ctrl->K.vmin);
                            }
                            break;

                        default:
                            GRTBadOptionError(K, "-K+%s is not supported.", token);
                            break;
                    }

                    token = strtok(NULL, "+");
                }

                GRT_SAFE_FREE_PTR(line);
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
                // 如果输入仅由数字、小数点和间隔符组成，则直接读取
                if(grt_string_composed_of(optarg, GRT_NUM_STR ".,")){
                    Ctrl->R.s_rs = grt_string_split(optarg, ",", &Ctrl->R.nr);
                } 
                // 否则从文件读取
                else {
                    FILE *fp = GRTCheckOpenFile(optarg, "r");
                    Ctrl->R.s_rs = grt_string_from_file(fp, &Ctrl->R.nr);
                    fclose(fp);
                }
                // 转为浮点数
                Ctrl->R.rs = (real_t*)realloc(Ctrl->R.rs, sizeof(real_t)*(Ctrl->R.nr));
                for(size_t i=0; i<Ctrl->R.nr; ++i){
                    Ctrl->R.rs[i] = atof(Ctrl->R.s_rs[i]);
                    if(Ctrl->R.rs[i] < 0.0){
                        GRTBadOptionError(R, "Can't set negative epicentral distance(%f).", Ctrl->R.rs[i]);
                    }
                }
                break;

            // 多线程数 -Pnthreads
            case 'P':
                Ctrl->P.active = true;
                if(1 != sscanf(optarg, "%d", &Ctrl->P.nthreads)){
                    GRTBadOptionError(P, "");
                };
                if(Ctrl->P.nthreads <= 0){
                    GRTBadOptionError(P, "Nonpositive value is not supported.");
                }
                grt_set_num_threads(Ctrl->P.nthreads);
                break;

            // 选择要计算的格林函数 -Ge|v|h|s
            case 'G': 
                Ctrl->G.active = true;
                // 先全部置否
                Ctrl->G.doEX = Ctrl->G.doVF = Ctrl->G.doHF = Ctrl->G.doDC = false;
                for(size_t i=0; i < strlen(optarg); ++i){
                    switch (optarg[i]) {
                        case 'e':
                            Ctrl->G.doEX = true;
                            break;
                        case 'v':
                            Ctrl->G.doVF = true;
                            break;
                        case 'h':
                            Ctrl->G.doHF = true;
                            break;
                        case 's':
                            Ctrl->G.doDC = true;
                            break;
                        default:
                            GRTBadOptionError(G, "unknown type %c.", optarg[i]);
                            break;
                    }
                }
                // 至少要有一个真
                if(!(Ctrl->G.doEX || Ctrl->G.doVF || Ctrl->G.doHF || Ctrl->G.doDC)){
                    GRTBadOptionError(G, "At least set one true value.");
                }
                break;

            // 输出波数积分中间文件， -Sidx1,idx2,idx3,...
            case 'S':
                Ctrl->S.active = true;
                // 如果非空，则读取对应索引，要求不能为负数；否则在 switch 外手动创建所有索引值
                if(optarg != NULL){
                    Ctrl->S.s_raw = strdup(optarg);
                    Ctrl->S.s_statsidxs = grt_string_split(optarg, ",", &Ctrl->S.nstatsidxs);
                    // 转为浮点数
                    Ctrl->S.statsidxs = (size_t*)realloc(Ctrl->S.statsidxs, sizeof(size_t)*(Ctrl->S.nstatsidxs));
                    for(size_t i=0; i<Ctrl->S.nstatsidxs; ++i){
                        int tmp = atoi(Ctrl->S.s_statsidxs[i]);
                        if(tmp < 0){
                            GRTBadOptionError(S, "index (%d) can't negative.", tmp);
                        }
                        Ctrl->S.statsidxs[i] = (size_t)tmp;
                    }
                }
                break;

            // 是否计算位移空间导数
            case 'e':
                Ctrl->e.active = true;
                break;
            
            GRT_Common_Options_in_Switch((char)(optopt));
        }
    } // END get options

    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, M);
    GRTCheckOptionActive(Ctrl, D);
    GRTCheckOptionActive(Ctrl, N);
    GRTCheckOptionActive(Ctrl, R);
    GRTCheckOptionActive(Ctrl, O);

    // 建立保存目录
    GRTCheckMakeDir(Ctrl->O.s_output_dir);

    // 在目录中保留命令
    char *dummy = NULL;
    GRT_SAFE_ASPRINTF(&dummy, "%s/command", Ctrl->O.s_output_dir);
    FILE *fp = GRTCheckOpenFile(dummy, "a");
    fprintf(fp, GRT_MAIN_COMMAND " ");  // 主程序名
    for(int i=0; i<argc; ++i){
        fprintf(fp, "%s ", argv[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
    GRT_SAFE_FREE_PTR(dummy);

}



/** 将一条数据反变换回时间域再进行处理，保存到SAC文件 */
static void write_one_to_sac(
    const char *srcname, const char ch, GRT_FFTW_HOLDER *fh, const real_t wI, 
    SACTRACE *sac, const char *s_output_subdir, const char *s_prefix,
    const int sgn, const cplx_t *grncplx)
{
    snprintf(sac->hd.kcmpnm, sizeof(sac->hd.kcmpnm), "%s%s%c", s_prefix, srcname, ch);
    
    char *s_outpath = NULL;
    GRT_SAFE_ASPRINTF(&s_outpath, "%s/%s.sac", s_output_subdir, sac->hd.kcmpnm);

    // 执行fft任务会修改数组，需重新置零
    grt_reset_fftw_holder_zero(fh);
    
    // 赋值复数，包括时移
    cplx_t cfac, ccoef;
    cfac = exp(I*PI2*fh->df*sac->hd.b);
    ccoef = sgn;
    // 只赋值有效长度，其余均为0
    for(size_t i = 0; i < fh->nf_valid; ++i){
        fh->W_f[i] = grncplx[i] * ccoef;
        ccoef *= cfac;
    }


    if(! fh->naive_inv){
        // 发起fft任务 
        fftw_execute(fh->plan);
    } else {
        grt_naive_inverse_transform_double(fh);
    }

    // 归一化，并处理虚频
    // 并转为 SAC 需要的单精度类型
    real_t fac, coef;
    coef = fh->df * exp(sac->hd.b * wI);
    fac = exp(wI*fh->dt);
    for(size_t i = 0; i < fh->nt; ++i){
        sac->data[i] = fh->w_t[i] * coef;
        coef *= fac;
    }

    // 以sac文件保存到本地
    grt_write_SACTRACE(s_outpath, sac);

    GRT_SAFE_FREE_PTR(s_outpath);
}



/** 子模块主函数 */
int greenfn_main(int argc, char **argv) {
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    // 传入参数 
    getopt_from_command(Ctrl, argc, argv);

    // 读入模型文件
    if((Ctrl->M.mod1d = grt_read_mod1d_from_file(Ctrl->M.s_modelpath, Ctrl->D.depsrc, Ctrl->D.deprcv, true)) == NULL){
        exit(EXIT_FAILURE);
    }
    GRT_MODEL1D *mod1d = Ctrl->M.mod1d;

    // 当震源位于液体层中时，仅允许计算爆炸源对应的格林函数
    // 程序结束前会输出对应警告
    if(mod1d->Vb[mod1d->isrc]==0.0){
        Ctrl->G.doHF = Ctrl->G.doVF = Ctrl->G.doDC = false;
    }

    // 最大最小速度
    real_t vmin, vmax;
    grt_get_mod1d_vmin_vmax(mod1d, &vmin, &vmax);

    // 参考最小速度
    if(Ctrl->K.vmin == 0.0){
        Ctrl->K.vmin = GRT_MAX(vmin, GRT_GREENFN_K_VMIN);
    } 

    // 判断是否要自动使用收敛方法
    if( ! Ctrl->C.active && fabs(Ctrl->D.deprcv - Ctrl->D.depsrc) <= GRT_MIN_DEPTH_GAP_SRC_RCV) {
        Ctrl->C.applyDCM = true;
    }

    // 时窗长度 
    Ctrl->N.winT = Ctrl->N.nt*Ctrl->N.dt;

    // 最大震中距
    real_t rmax = Ctrl->R.rs[grt_findMax_real_t(Ctrl->R.rs, Ctrl->R.nr)];   

    // 时窗最大截止时刻
    real_t tmax = Ctrl->E.delayT0 + Ctrl->N.winT;
    if(Ctrl->E.delayV0 > 0.0)   tmax += rmax/Ctrl->E.delayV0;

    // 自动选择积分间隔，默认使用传统离散波数积分
    // 自动选择会给出很保守的值（较大的Length）
    if(Ctrl->L.Length == 0.0){
        Ctrl->L.Length = 15.0; 
        real_t jus = GRT_SQUARE(vmax*tmax) - GRT_SQUARE(Ctrl->D.deprcv - Ctrl->D.depsrc);
        if(jus >= 0.0){
            Ctrl->L.Length = GRT_MAX(1.0 + sqrt(jus)/rmax + 0.5, Ctrl->L.Length); // +0.5为保守值
        }
    }

    // 虚频率
    Ctrl->N.wI = Ctrl->N.zeta*PI/Ctrl->N.winT;

    // 定义要计算的频率、时窗等
    Ctrl->N.nf = Ctrl->N.nt/2 + 1;
    Ctrl->N.df = 1.0/Ctrl->N.winT;
    Ctrl->N.freqs = (real_t*)malloc(Ctrl->N.nf*sizeof(real_t));
    for(size_t i=0; i<Ctrl->N.nf; ++i){
        Ctrl->N.freqs[i] = i*Ctrl->N.df;
    }

    // 如果只传入了 -S, 未指定索引，则默认所有频率索引
    if(Ctrl->S.active && Ctrl->S.statsidxs == NULL){
        // 另外两个字符相关的指针仍指向 NULL
        Ctrl->S.nstatsidxs = Ctrl->N.nf;
        Ctrl->S.statsidxs = (size_t*)realloc(Ctrl->S.statsidxs, sizeof(size_t)*(Ctrl->S.nstatsidxs));
        for(size_t i=0; i < Ctrl->S.nstatsidxs; ++i){
            Ctrl->S.statsidxs[i] = i;
        }
        Ctrl->S.s_raw = strdup("(all)");
    }

    // 自定义频段
    Ctrl->H.nf1 = 0; 
    Ctrl->H.nf2 = Ctrl->N.nf-1;
    if(Ctrl->H.freq1 > 0.0){
        Ctrl->H.nf1 = GRT_MIN(ceil(Ctrl->H.freq1/Ctrl->N.df), Ctrl->N.nf-1);
    }
    if(Ctrl->H.freq2 > 0.0){
        Ctrl->H.nf2 = GRT_MIN(floor(Ctrl->H.freq2/Ctrl->N.df), Ctrl->N.nf-1);
    }
    Ctrl->H.nf2 = GRT_MAX(Ctrl->H.nf1, Ctrl->H.nf2);

    // 波数积分中间文件输出目录
    if(Ctrl->S.active){
        Ctrl->S.s_statsdir = NULL;
        GRT_SAFE_ASPRINTF(&Ctrl->S.s_statsdir, "%s_grtstats", Ctrl->O.s_output_dir);
        
        // 建立保存目录
        GRTCheckMakeDir(Ctrl->S.s_statsdir);
        GRT_SAFE_ASPRINTF(&Ctrl->S.s_statsdir, "%s/%s_%s_%s", Ctrl->S.s_statsdir, Ctrl->M.s_modelname, Ctrl->D.s_depsrc, Ctrl->D.s_deprcv);
        GRTCheckMakeDir(Ctrl->S.s_statsdir);
    }

    // 建立格林函数的complex数组
    pt_cplxChnlGrid *grn = (pt_cplxChnlGrid *) calloc(Ctrl->R.nr, sizeof(*grn));
    pt_cplxChnlGrid *grn_uiz = (Ctrl->e.active)? (pt_cplxChnlGrid *) calloc(Ctrl->R.nr, sizeof(*grn_uiz)) : NULL;
    pt_cplxChnlGrid *grn_uir = (Ctrl->e.active)? (pt_cplxChnlGrid *) calloc(Ctrl->R.nr, sizeof(*grn_uir)) : NULL;

    for(size_t ir=0; ir<Ctrl->R.nr; ++ir){
        GRT_LOOP_ChnlGrid(im, c){
            grn[ir][im][c] = (cplx_t*)calloc(Ctrl->N.nf, sizeof(cplx_t));
            if(grn_uiz)  grn_uiz[ir][im][c] = (cplx_t*)calloc(Ctrl->N.nf, sizeof(cplx_t));
            if(grn_uir)  grn_uir[ir][im][c] = (cplx_t*)calloc(Ctrl->N.nf, sizeof(cplx_t));
        }
    }

        
    // 波数积分方法
    K_INTEG_METHOD KMET = {0};
    {   
        real_t hs = GRT_MAX(fabs(mod1d->depsrc - mod1d->deprcv), GRT_MIN_DEPTH_GAP_SRC_RCV);
        KMET.k0 = Ctrl->K.k0 * PI / hs;
        KMET.ampk = Ctrl->K.ampk;
        KMET.keps = (Ctrl->C.applyPTAM || Ctrl->C.applyDCM)? 0.0 : Ctrl->K.keps;  // 如果使用了显式收敛方法，则不使用keps进行收敛判断
        KMET.vmin = Ctrl->K.vmin;
        
        KMET.kcut = Ctrl->L.kcut / rmax;

        KMET.dk = PI2 / (Ctrl->L.Length * rmax);

        KMET.applyFIM = Ctrl->L.FIM.active;
        KMET.filondk = (Ctrl->L.FIM.active) ? PI2 / (Ctrl->L.FIM.Length * rmax) : 0.0;

        KMET.applySAFIM = Ctrl->L.SAFIM.active;
        KMET.sa_tol = Ctrl->L.SAFIM.tol;
        
        KMET.applyDCM = Ctrl->C.applyDCM;
        KMET.applyPTAM = Ctrl->C.applyPTAM;
    }

    // 在计算前打印所有参数
    if(! Ctrl->s.active){
        print_Ctrl(Ctrl);
    }
    

    //==============================================================================
    // 计算格林函数
    grt_integ_grn_spec(
        mod1d, Ctrl->H.nf1, Ctrl->H.nf2, Ctrl->N.freqs, Ctrl->R.nr, Ctrl->R.rs, Ctrl->N.wI, Ctrl->N.keepAllFreq,
        &KMET, !Ctrl->s.active, Ctrl->e.active, 
        grn, grn_uiz, grn_uir,
        Ctrl->S.s_statsdir, Ctrl->S.nstatsidxs, Ctrl->S.statsidxs
    );
    //==============================================================================

    // 使用fftw3做反傅里叶变换，并保存到 SAC 
    // 其中考虑了升采样倍数
    GRT_FFTW_HOLDER *fh = grt_create_fftw_holder_C2R_1D(
        Ctrl->N.nt*Ctrl->N.upsample_n, Ctrl->N.dt/Ctrl->N.upsample_n, Ctrl->N.nf, Ctrl->N.df);

    real_t (* travtPS)[2] = (real_t (*)[2])calloc(Ctrl->R.nr, sizeof(real_t)*2);

    char *s_output_dirprefx = NULL;
    GRT_SAFE_ASPRINTF(&s_output_dirprefx, "%s/%s_%s_%s", Ctrl->O.s_output_dir, Ctrl->M.s_modelname, Ctrl->D.s_depsrc, Ctrl->D.s_deprcv);
    
    // 建立SAC头文件，包含必要的头变量
    SACTRACE *sac = grt_new_SACTRACE(fh->dt, fh->nt, Ctrl->E.delayT0);
    // 发震时刻作为参考时刻
    sac->hd.o = 0.0; 
    sac->hd.iztype = IO; 
    // 记录震源和台站深度
    sac->hd.evdp = Ctrl->D.depsrc; // km
    sac->hd.stel = (-1.0)*Ctrl->D.deprcv*1e3; // m
    // 写入虚频率
    sac->hd.user0 = Ctrl->N.wI;
    // 写入接受点的Vp,Vs,rho
    sac->hd.user1 = mod1d->Va[mod1d->ircv];
    sac->hd.user2 = mod1d->Vb[mod1d->ircv];
    sac->hd.user3 = mod1d->Rho[mod1d->ircv];
    sac->hd.user4 = mod1d->Qainv[mod1d->ircv];
    sac->hd.user5 = mod1d->Qbinv[mod1d->ircv];
    // 写入震源点的Vp,Vs,rho
    sac->hd.user6 = mod1d->Va[mod1d->isrc];
    sac->hd.user7 = mod1d->Vb[mod1d->isrc];
    sac->hd.user8 = mod1d->Rho[mod1d->isrc];
    
    // 做反傅里叶变换，保存SAC文件
    for(size_t ir = 0; ir < Ctrl->R.nr; ++ir){
        real_t dist = Ctrl->R.rs[ir];
        sac->hd.dist = dist;

        // 文件保存子目录
        char *s_output_subdir = NULL;
        
        GRT_SAFE_ASPRINTF(&s_output_subdir, "%s_%s", s_output_dirprefx, Ctrl->R.s_rs[ir]);
        GRTCheckMakeDir(s_output_subdir);

        // 时间延迟 
        real_t delayT = Ctrl->E.delayT0;
        if(Ctrl->E.delayV0 > 0.0)   delayT += sqrt( GRT_SQUARE(dist) + GRT_SQUARE(Ctrl->D.deprcv - Ctrl->D.depsrc) ) / Ctrl->E.delayV0;
        // 修改SAC头段时间变量
        sac->hd.b = delayT;

        // 计算理论走时
        sac->hd.t0 = grt_compute_travt1d(mod1d->Thk, mod1d->Va, mod1d->n, mod1d->isrc, mod1d->ircv, dist);
        strcpy(sac->hd.kt0, "P");
        sac->hd.t1 = grt_compute_travt1d(mod1d->Thk, mod1d->Vb, mod1d->n, mod1d->isrc, mod1d->ircv, dist);
        strcpy(sac->hd.kt1, "S");

        GRT_LOOP_ChnlGrid(im, c){
            if(! Ctrl->G.doEX  && im==GRT_SRC_M_EX_INDEX)  continue;
            if(! Ctrl->G.doVF  && im==GRT_SRC_M_VF_INDEX)  continue;
            if(! Ctrl->G.doHF  && im==GRT_SRC_M_HF_INDEX)  continue;
            if(! Ctrl->G.doDC  && im>=GRT_SRC_M_DD_INDEX)  continue;

            int modr = GRT_SRC_M_ORDERS[im];
            int sgn=1;  // 用于反转Z分量

            if(modr==0 && GRT_ZRT_CODES[c]=='T')  continue;  // 跳过输出0阶的T分量

            // Z分量反转
            sgn = (GRT_ZRT_CODES[c]=='Z') ? -1 : 1;

            char ch = GRT_ZRT_CODES[c];

            write_one_to_sac(GRT_SRC_M_NAME_ABBR[im], ch, fh, Ctrl->N.wI, sac, s_output_subdir, "", sgn, grn[ir][im][c]);

            if(Ctrl->e.active){
                write_one_to_sac(GRT_SRC_M_NAME_ABBR[im], ch, fh, Ctrl->N.wI, sac, s_output_subdir, "z", sgn*(-1), grn_uiz[ir][im][c]);
                write_one_to_sac(GRT_SRC_M_NAME_ABBR[im], ch, fh, Ctrl->N.wI, sac, s_output_subdir, "r", sgn,      grn_uir[ir][im][c]);
            }
        }

        GRT_SAFE_FREE_PTR(s_output_subdir);
        
        // 记录走时
        if(travtPS != NULL){
            travtPS[ir][0] = sac->hd.t0;
            travtPS[ir][1] = sac->hd.t1;
        }
    }

    // 输出警告：当震源位于液体层中时，仅允许计算爆炸源对应的格林函数
    if(mod1d->Vb[mod1d->isrc]==0.0){
        GRTRaiseWarning(
            "The source is located in the liquid layer, "
            "therefore only the Green's Funtions for the Explosion source will be computed.\n");
    }

    grt_free_SACTRACE(sac);

    
    // 打印走时
    if( ! Ctrl->s.active){
        printf("\n\n");
        printf("------------------------------------------------\n");
        printf(" Distance(km)     Tp(secs)         Ts(secs)     \n");
        for(size_t ir=0; ir<Ctrl->R.nr; ++ir){
            printf(" %-15s  %-15.3f  %-15.3f\n", Ctrl->R.s_rs[ir], travtPS[ir][0], travtPS[ir][1]);
        }
        printf("------------------------------------------------\n");
        printf("\n");
    }

    // 释放内存
    for(size_t ir=0; ir<Ctrl->R.nr; ++ir){
        GRT_LOOP_ChnlGrid(im, c){
            GRT_SAFE_FREE_PTR(grn[ir][im][c]);
            if(grn_uiz) GRT_SAFE_FREE_PTR(grn_uiz[ir][im][c]);
            if(grn_uir) GRT_SAFE_FREE_PTR(grn_uir[ir][im][c]);
        }
    }
    GRT_SAFE_FREE_PTR(grn);
    GRT_SAFE_FREE_PTR(grn_uiz);
    GRT_SAFE_FREE_PTR(grn_uir);
    GRT_SAFE_FREE_PTR(travtPS);

    grt_destroy_fftw_holder(fh);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}

