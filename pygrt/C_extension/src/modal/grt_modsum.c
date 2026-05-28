/**
 * @file   grt_modsum.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 *     根据计算的相速度频散(本征值)，使用模态叠加（modal summation）计算面波格林函数
 * 
 */

#include "grt.h"

// 一些变量的非零默认值
#define GRT_GREENFN_N_UPSAMPLE    1
#define GRT_GREENFN_G_EX       true
#define GRT_GREENFN_G_VF       true
#define GRT_GREENFN_G_HF       true
#define GRT_GREENFN_G_DC       true

/** 该子模块的参数控制结构体 */
typedef struct {
    /* 输入频散文件 */
    struct {
        bool active;
        char *s_filepath;
    } C;
    /* 阶数范围 */
    struct {
        bool active;
        size_t *modes;
        size_t nmode;
    } N;
    /* 频率范围 */
    struct {
        bool active;
        size_t nf;
        real_t *freqs;
    } F;
    /* 输入震源和台站深度，计算面波格林函数 */
    struct {
        bool active;
        real_t depsrc;
        real_t deprcv;
        char *s_depsrc;
        char *s_deprcv;
    } D;
    /* 面波格林函数输出路径 */
    struct {
        bool active;
        char *s_output_dir;
    } O;
    /** 时间延迟 */
    struct {
        bool active;
        real_t delayT0;
        real_t delayV0;
        bool refFirstP;  ///< 是否参考初至P
    } E;
    /* 震中距数组 */
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
    /* 升采样倍数 */
    struct {
        bool active;
        int upsample_n;
    } W;

} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    // C
    GRT_SAFE_FREE_PTR(Ctrl->C.s_filepath);
    // N
    GRT_SAFE_FREE_PTR(Ctrl->N.modes);

    // D
    GRT_SAFE_FREE_PTR(Ctrl->D.s_depsrc);
    GRT_SAFE_FREE_PTR(Ctrl->D.s_deprcv);
    // O
    GRT_SAFE_FREE_PTR(Ctrl->O.s_output_dir);
    // R
    GRT_SAFE_FREE_PTR_ARRAY(Ctrl->R.s_rs, Ctrl->R.nr);
    GRT_SAFE_FREE_PTR(Ctrl->R.s_raw);
    GRT_SAFE_FREE_PTR(Ctrl->R.rs);

    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt modsum] %s\n\n", GRT_VERSION);printf(
"    Compute the surface-wave Green's Functions using the Modal Summation.\n"
"\n"
"    The modelpath was stored in the dispersion result from \n"
"    module `eigenv`, therefore this module requires the modelpath\n"
"    to exist and will automatically read the model.\n"
"\n"
"    To apply IFFT, `eigenv` must use the form of -F<f1>/<f2>/<df>\n"
"    to define the equidistant frequency points.\n"
"\n\n"
"Usage:\n"
"---------------------------------------------------------------------\n"
"    grt modsum -C<path> -D<depsrc>/<deprcv> -O<outdir> \n"
"         -R<r1>,<r2>[,...]|<r1>/<r2>/<dr>|<file> [-F<f1>/<f2>]\n"
"         [-N[<n1>][/<n2>][/<dn>]] [-E[p]<t0>[/<v0>]] [-Ge|v|h|s]\n"
"         [-P<nthreads>] [-W<fac>] [-e] [-h]\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -C<path>    Input dispersion result file from module `eigenv` \n"
"                (.nc format)\n"
"\n"
"    -D<depsrc>/<deprcv>\n"
"                 <depsrc>: source depth (km).\n"
"                 <deprcv>: receiver depth (km).\n"
"\n"
"    -O<outdir>   Directorypath of output for saving.\n"
"\n"
"    -R<r1>,<r2>[,...]|<r1>/<r2>/<dr>|<file>\n"
"                 Multiple epicentral distances (km), support three ways:\n"
"                 + <r1>,<r2>[,...]: seperated by comma.\n"
"                 + <r1>/<r2>/<dr>:  equal distance <dr> within [r1,r2].\n"
"                 + <file>: each line contains a distance value.\n"
"\n"
"    -F<f1>/<f2>\n"
"                Select the frequency range from the input file.\n"
"                <f1>: start frequency (Hz)\n"
"                <f2>: end frequency (Hz)\n"
"\n"
"    -N[<n1>][/<n2>][/<dn>]\n"
"                Select the order range from the input file, \n"
"                0 means fundamental mode.\n"
"                <n1>: start order\n"
"                <n2>: end order\n"
"                <dn>: order interval\n"
"                + If -N is not given, equal to -N0.\n"
"                + If set an empty -N, means all orders will be considered.\n"
"                + If set -N<n1>/<n2>, means set min/max range.\n"
"\n"
"    -E[p]<t0>[/<v0>]\n"
"                 Introduce the time delay in results. The total \n"
"                 delay = <t0> + dist/<v0>, dist is the\n"
"                 straight-line distance between source and \n"
"                 receiver.\n"
"                 <t0>: reference delay (s), default t0=0.0\n"); printf(
"                 <v0>: reference velocity (km/s), \n"
"                       default 0.0 not use. \n"
"                 -Ep<t0>: \n"
"                       the delay (the begining time) will be <t0> + first P, \n"
"                       e.g., -Ep-10.\n"
"\n"
"    -Ge|v|h|s\n"
"                 Designed to choose which kind of source's Green's \n"
"                 functions will be computed, default is all (-Gevhs). \n"
"                 <e>: Explosion (EX)\n"
"                 <v>: Vertical Force (VF)\n"
"                 <h>: Horizontal Force (HF)\n"
"                 <s>: Shear (DC)\n"
"\n"
"    -P<n>        Number of threads. Default use all cores.\n"
"\n"
"    -W<fac>      upsampling factor (integer)\n"
"                 i.e.  nt <-- nt * <fac>\n"
"                       dt <-- dt / <fac>\n"
"\n"
"    -e           Compute the spatial derivatives, ui_z and ui_r,\n"
"                 of displacement u. In filenames, prefix \"r\" means \n"
"                 ui_r and \"z\" means ui_z. The units of derivatives\n"
"                 for different sources are: \n"
"                 + Explosion:     1e-25 /(dyne-cm)\n"
"                 + Single Force:  1e-20 /(dyne)\n"
"                 + Shear:         1e-25 /(dyne-cm)\n"
"\n"
"    -h           Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt modsum -Cphase_R.nc -D2/0 -R5/100/5 -OGRN_NM -W5 -e\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    // 先为个别参数设置非0初始值
    Ctrl->W.upsample_n = GRT_GREENFN_N_UPSAMPLE;

    Ctrl->G.doEX = GRT_GREENFN_G_EX;
    Ctrl->G.doVF = GRT_GREENFN_G_VF;
    Ctrl->G.doHF = GRT_GREENFN_G_HF;
    Ctrl->G.doDC = GRT_GREENFN_G_DC;

    int opt;
    while ((opt = getopt(argc, argv, ":C:N::F:D:O:E:R:G:W:P:esh" )) != -1) {
        switch (opt) {
            // 频散文件
            case 'C':
                Ctrl->C.active = true;
                Ctrl->C.s_filepath = strdup(optarg);
                break;

            // 阶数范围 -N[n1][/n2][/dn]
            case 'N':
                Ctrl->N.active = true;
                // 没有参数，直接跳过
                if(optarg==NULL)  break;
                {
                    ssize_t a1, a2, dif;
                    a1 = a2 = 0;
                    dif = 1;
                    int nscan = sscanf(optarg, "%zd/%zd/%zd", &a1, &a2, &dif);
                    if( !(nscan == 1 || nscan == 2 || nscan == 3)){
                        GRTBadOptionError(N, "");
                    }
                    if(nscan == 1 && a1 < 0){
                        GRTBadOptionError(N, "Can't set a single negative value.");
                    }
                    if(nscan == 1){
                        a2 = a1;
                        dif = 1;
                    }
                    else if(nscan >=2){
                        if(dif <= 0){
                            GRTBadOptionError(N, "Can't set negative dn(%zd).", dif);
                        }
                        if(a1 < 0 || a2 < 0){
                            GRTBadOptionError(N, "Can't set nonpositive n1(%zd), n2(%zd).", a1, a2);
                        }
                        if(a1 > a2){
                            GRTBadOptionError(N, "n1(%zd) > n2(%zd).", a1, a2);
                        }
                    }
                    Ctrl->N.nmode = floor((a2-a1)/dif) + 1;
                    Ctrl->N.modes = (size_t *)calloc(Ctrl->N.nmode, sizeof(size_t));
                    for(size_t i=0; i < Ctrl->N.nmode; ++i){
                        Ctrl->N.modes[i] = a1 + dif*i;
                    }
                }
                break;

            // 频率值 -Ff1/f2  指定频率范围
            case 'F':
                Ctrl->F.active = true;
                Ctrl->F.freqs = (real_t *)calloc(2, sizeof(real_t));
                if(2 != sscanf(optarg, "%lf/%lf", Ctrl->F.freqs, Ctrl->F.freqs+1)){
                    GRTBadOptionError(F, "");
                };
                if(Ctrl->F.freqs[0] < 0.0){
                    GRTBadOptionError(F, "Can't set negative freq1(%f).", Ctrl->F.freqs[0]);
                }
                if(Ctrl->F.freqs[1] < 0.0){
                    GRTBadOptionError(F, "Can't set negative freq2(%f).", Ctrl->F.freqs[1]);
                }
                if(Ctrl->F.freqs[0] > Ctrl->F.freqs[1]){
                    GRTBadOptionError(F, "Positive freq1 should be less than positive freq2.");
                }
                break;

            // 输入震源和台站深度、震中距和保存目录，计算面波格林函数
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
            
            // 输出路径 -Ooutput_dir
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.s_output_dir = strdup(optarg);
                break;

            // 时间延迟 -E[p]<t0>[/<v0>]
            case 'E':
                Ctrl->E.active = true;
                if(optarg[0] == 'p'){
                    if(1 != sscanf(optarg+1, "%lf", &Ctrl->E.delayT0)){
                        GRTBadOptionError(E, "");
                    };
                    if(Ctrl->E.delayT0 >= 0.0){
                        GRTBadOptionError(E, "Can't set positive t0(%f) in -Ep.", Ctrl->E.delayT0);
                    }
                    Ctrl->E.refFirstP = true;
                } else {
                    if(0 == sscanf(optarg+1, "%lf/%lf", &Ctrl->E.delayT0, &Ctrl->E.delayV0)){
                        GRTBadOptionError(E, "");
                    };
                    if(Ctrl->E.delayV0 < 0.0){
                        GRTBadOptionError(E, "Can't set negative v0(%f) in -E.", Ctrl->E.delayV0);
                    }
                }
                break;

            // 震中距数组，-R<r1>,<r2>[,...]|<r1>/<r2>/<dr>|<file>
            case 'R':
                Ctrl->R.active = true;
                Ctrl->R.s_raw = strdup(optarg);
                {
                    real_t a1, a2, delta;
                    
                    // 如果输入仅由数字、小数点和间隔符组成，则直接读取
                    if(grt_string_composed_of(optarg, GRT_NUM_STR "eE+-" ".,")){
                        Ctrl->R.s_rs = grt_string_split(optarg, ",", &Ctrl->R.nr);
                    }
                    // 尝试按照 <r1>/<r2>/<dr> 读取
                    else if(3 == sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &delta)){
                        if(delta <= 0){
                            GRTBadOptionError(R, "Can't set nonpositive dr(%f)", delta);
                        }
                        if(a1 > a2){
                            GRTBadOptionError(R, "r1(%f) > r2(%f).", a1, a2);
                        }

                        // 根据最大的小数位数来设置 s_rs
                        int place = 0;
                        real_t test[3] = {0};
                        test[0] = a1;
                        test[1] = a2;
                        test[2] = delta;
                        do {
                            bool agree = true;
                            for(size_t i = 0; i < 3; ++i)  {
                                real_t new_frac = test[i] - (int)test[i];
                                test[i] *= 10.0;
                                // 如果乘以10后变成整数（或非常接近整数）
                                if (fabs(new_frac) > 1e-8) {
                                    agree = false;
                                }
                            }
                            if(agree) break;
                            place++;
                        } while (place < 8);  // 最多小数点后 8 位

                        Ctrl->R.nr = floor((a2-a1)/delta) + 1;
                        Ctrl->R.s_rs = (char **)calloc(Ctrl->R.nr, sizeof(char*) * Ctrl->R.nr);
                        for(size_t ir = 0; ir < Ctrl->R.nr; ++ir){
                            GRT_SAFE_ASPRINTF(&Ctrl->R.s_rs[ir], "%.*f", place, a1 + delta*ir);
                        }
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
                }
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

            // 升采样倍数
            case 'W':
                Ctrl->W.active = true;
                if(1 != sscanf(optarg, "%d", &Ctrl->W.upsample_n)){
                    GRTBadOptionError(W, "");
                };
                if(Ctrl->W.upsample_n <= 0){
                    GRTBadOptionError(W, "need positive integer, but get (%d).", Ctrl->W.upsample_n);
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

            // 是否计算位移空间导数
            case 'e':
                Ctrl->e.active = true;
                break;

            GRT_Common_Options_in_Switch((char)(optopt));

        }
    }


    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, C);
    GRTCheckOptionActive(Ctrl, D);
    GRTCheckOptionActive(Ctrl, O);
    GRTCheckOptionActive(Ctrl, R);

    // -N 的默认项，仅处理基阶
    if( ! Ctrl->N.active ){
        Ctrl->N.nmode = 1;
        Ctrl->N.modes = (size_t *)calloc(Ctrl->N.nmode, sizeof(size_t));
        Ctrl->N.modes[0] = 0;
    }
    
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


/* 子模块主函数 */
int modsum_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    // 传入参数 
    getopt_from_command(Ctrl, argc, argv);

    // 读取频散
    char *modelpath = NULL;
    EIGENV_INFO *eigmet = (EIGENV_INFO *)calloc(1, sizeof(EIGENV_INFO));
    grt_read_dispersion(Ctrl->C.s_filepath, eigmet, &modelpath);

    // 读取模型（不插入震源和台站的虚拟层）
    MODEL1D *mod1d = NULL;
    if((mod1d = grt_read_mod1d_from_file(modelpath, -1.0, -1.0, true)) ==NULL){
        exit(EXIT_FAILURE);
    }

    // 当震源位于液体层中时，仅允许计算爆炸源对应的格林函数
    // 程序结束前会输出对应警告
    if(mod1d->isLiquid[mod1d->isrc]){
        Ctrl->G.doHF = Ctrl->G.doVF = Ctrl->G.doDC = false;
    }

    // 不再插入虚拟层，直接记录震源和台站所在的真实层位
    mod1d->depsrc = Ctrl->D.depsrc;
    mod1d->deprcv = Ctrl->D.deprcv;
    mod1d->ircvup = (Ctrl->D.depsrc >= Ctrl->D.deprcv);
    for(mod1d->isrc=0; mod1d->isrc < mod1d->n-1; ++mod1d->isrc){
        if(mod1d->Dep[mod1d->isrc+1] >= Ctrl->D.depsrc)  break; 
    }
    for(mod1d->ircv=0; mod1d->ircv < mod1d->n-1; ++mod1d->ircv){
        if(mod1d->Dep[mod1d->ircv+1] >= Ctrl->D.deprcv)  break; 
    }

    // 根据命令行参数确定出所需的部分频散信息
    EIGENFN_INFO *eigfnmet = (EIGENFN_INFO *)calloc(1, sizeof(EIGENFN_INFO));
    grt_filter_eigenfn_info(
        Ctrl->F.nf, Ctrl->F.freqs, true,
        Ctrl->N.nmode, Ctrl->N.modes, eigmet, eigfnmet);
    // 使用模型层计算能量积分
    eigfnmet->cpar_nz = mod1d->n;
    
    // 如果要计算面波格林函数，对生成频散文件的参数有限制，在 eigenv 必须以 -F0/fmax/df 的形式使用
    // 不能只传入一个频率值
    if(eigfnmet->nf <= 1){
        GRTRaiseError("Modal summation need at least two sampled freqs, but only found %zu.\n", eigfnmet->nf);
    }
    // 原频散文件中，起始频率必须等于频率间隔，且不能以周期的形式给定
    if(fabs(2.0*eigmet->freqs[0] - eigmet->freqs[1]) > 1e-8){
        GRTRaiseError("Modal summation need dispersion file generated like `grt eigenv -F0/fmax/df`.\n");
    }

    // 初始化 eigfn
    eigfnmet->eigfn = (EIGENFN **)calloc(eigfnmet->nf, sizeof(EIGENFN *));
    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        size_t cnum = eigfnmet->eigv[iw].n;
        eigfnmet->eigfn[iw] = (EIGENFN *)calloc(cnum, sizeof(EIGENFN));
    }

    // FFT需要的参数，包括零频和区域外的（从原频散文件确定，因此本模块中 -F 的作用仅用于筛选频段）
    real_t fft_df = eigmet->freqs[1] - eigmet->freqs[0];
    size_t fft_nf = eigmet->nf + 1; // + 1 是为零频
    size_t fft_nt = 2*(fft_nf - 1);
    real_t fft_dt = 0.5 / eigmet->freqs[eigmet->nf-1];

    // 格林函数频谱
    GRNSPEC *grn = &(GRNSPEC){0};
    {
        grn->nf = fft_nf;
        grn->freqs = (real_t *)calloc(fft_nf, sizeof(real_t));   // 单独增加 0 频
        memcpy(grn->freqs+1, eigmet->freqs, sizeof(real_t)*eigmet->nf);
        grn->nf1 = 0;
        grn->nf2 = fft_nf-1;
        grn->nr = Ctrl->R.nr;
        grn->rs = Ctrl->R.rs;
        grn->wI = 0.0;
        grn->keepAllFreq = true;
        grn->calc_upar = Ctrl->e.active;

        grt_grnspec_allocate_u(grn);
        grn->statsstr = NULL;
        grn->nstatsidxs = 0;
        grn->statsidxs = NULL;
    }

    grt_modsum_grn_spec(mod1d, eigmet->wtype, eigfnmet, grn);


    // 使用fftw3做反傅里叶变换，并保存到 SAC 
    // 其中考虑了升采样倍数
    GRT_FFTW_HOLDER *fh = grt_create_fftw_holder_C2R_1D(fft_nt*Ctrl->W.upsample_n, fft_dt/Ctrl->W.upsample_n, fft_nf, fft_df);
    
    // 建立SAC头文件，包含必要的头变量
    SACTRACE *sac = grt_new_SACTRACE(fh->dt, fh->nt, 0.0);
    // 发震时刻作为参考时刻
    sac->hd.o = 0.0; 
    sac->hd.iztype = IO; 
    // 记录震源和台站深度
    sac->hd.evdp = Ctrl->D.depsrc; // km
    sac->hd.stel = (-1.0)*Ctrl->D.deprcv*1e3; // m
    // 写入虚频率
    sac->hd.user0 = 0.0;
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
    
    // 为每个震中距设置对应的变量
    real_t (*travtPS)[2] = (real_t (*)[2])calloc(grn->nr, sizeof(real_t)*2);
    real_t *begintimes = (real_t *)calloc(grn->nr, sizeof(real_t));
    char **outputdirs = (char **)calloc(grn->nr, sizeof(char *));
    for(size_t ir = 0; ir < Ctrl->R.nr; ++ir){
        real_t dist = Ctrl->R.rs[ir];

        outputdirs[ir] = NULL;
        GRT_SAFE_ASPRINTF(&outputdirs[ir], "%s/%s_%s_%s_%s", 
            Ctrl->O.s_output_dir, grt_get_basename(modelpath), Ctrl->D.s_depsrc, Ctrl->D.s_deprcv, Ctrl->R.s_rs[ir]);

        // 计算理论走时
        travtPS[ir][0] = grt_compute_travt1d(mod1d->Thk, mod1d->Va, mod1d->n, mod1d->isrc, mod1d->ircv, dist);
        travtPS[ir][1] = grt_compute_travt1d(mod1d->Thk, mod1d->Vb, mod1d->n, mod1d->isrc, mod1d->ircv, dist);

        // 时间延迟
        {
            real_t delayT = 0.0;
            if (! Ctrl->E.refFirstP){
                delayT = Ctrl->E.delayT0;
                if(Ctrl->E.delayV0 > 0.0)   delayT += sqrt( GRT_SQUARE(dist) + GRT_SQUARE(Ctrl->D.deprcv - Ctrl->D.depsrc) ) / Ctrl->E.delayV0;
            } else {
                delayT = Ctrl->E.delayT0 + sac->hd.t0;
            }
            begintimes[ir] = delayT;
        }
    }

    const char *validChnls = (eigmet->wtype == GRT_DISPERSION_RAYL)? "ZR" : "T";

    // 反傅里叶变换后保存到SAC文件
    grt_grnspec_write_sac(
        grn, travtPS, begintimes, outputdirs, fh, sac, 
        validChnls, true, Ctrl->G.doEX, Ctrl->G.doVF, Ctrl->G.doHF, Ctrl->G.doDC);

    // 输出警告：当震源位于液体层中时，仅允许计算爆炸源对应的格林函数
    if(mod1d->isLiquid[mod1d->isrc]){
        GRTRaiseWarning(
            "The source is located in the liquid layer, "
            "therefore only the Green's Funtions for the Explosion source will be computed.\n");
    }

    // 释放内存
    GRT_SAFE_FREE_PTR(grn->freqs);
    grt_grnspec_free_u(grn);
    grt_free_SACTRACE(sac);
    grt_destroy_fftw_holder(fh);
    
    GRT_SAFE_FREE_PTR(modelpath);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
