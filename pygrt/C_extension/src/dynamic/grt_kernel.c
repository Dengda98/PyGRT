/**
 * @file   grt_kernel.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-12
 * 
 *    计算不同频率、不同相速度的核函数
 * 
 */

#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/integral/kernel.h"
#include "grt/integral/iostats.h"
#include "grt/common/util.h"
#include "grt/common/progressbar.h"

#include "grt.h"

// 一些变量的非零默认值
#define GRT_GREENFN_F_ZETA        20.0

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
    /* 相速度搜索范围 */
    struct {
        bool active;
        real_t cmin;
        real_t cmax;
        real_t dc;
        real_t *c_phases;
        size_t nc;
        bool isset;
    } C;
    /* 频率范围 */
    struct {
        bool active;
        real_t *freqs;
        size_t nf;
        real_t zeta;  ///< 虚频率系数， w <- w - zeta*PI/r* 1j
        real_t wI;    ///< 虚频率  zeta*PI/r
    } F;
    /** 输出目录 */
    struct {
        bool active;
        char *s_output_dir;
    } O;
    /** 是否计算空间导数 */
    struct {
        bool active;
    } e;
    /** 多线程 */
    struct {
        bool active;
        int nthreads; ///< 线程数
    } P;

} GRT_MODULE_CTRL;

/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    // M
    GRT_SAFE_FREE_PTR(Ctrl->M.s_modelpath);
    grt_free_mod1d(Ctrl->M.mod1d);

    // D
    GRT_SAFE_FREE_PTR(Ctrl->D.s_depsrc);
    GRT_SAFE_FREE_PTR(Ctrl->D.s_deprcv);

    // O
    GRT_SAFE_FREE_PTR(Ctrl->O.s_output_dir);

    // F
    GRT_SAFE_FREE_PTR(Ctrl->F.freqs);

    // C
    GRT_SAFE_FREE_PTR(Ctrl->C.c_phases);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt kernel] %s\n\n", GRT_VERSION);printf(
"    Compute kernel functions at different frequencies and\n"
"    phase velocities.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt kernel -M<model> -D<depsrc>/<deprcv> -C[<cmin>/<cmax>/]<dc>\n"
"               -F<f1>/<f2>/<df>[+w<zeta>] -O<outdir> [-P<nthreads>]\n"
"               [-e] [-h]\n"
"\n"
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
"    -C[<cmin>/<cmax>/]<dc>\n"
"                 <cmin>: minimum phase velocity.\n"
"                 <cmax>: maximum phase velocity.\n"
"                 <dc>:   interval. \n"
"                 If only provide <dc>, then\n"
"                     <cmin> = 0.8*MIN{Vel in model && Vel!=0}\n"
"                     <cmax> = MAX{Vel in model}.\n"
"                 All units are km/s.\n"
"\n"
"    -F<f1>/<f2>/<df>[+w<zeta>]\n"
"                 <f1>: minimum frequency (auto skip 0).\n"
"                 <f2>: maximum frequency.\n"
"                 <df>: interval.\n"
"                 +w<zeta>: define the coefficient of imaginary \n"
"                           frequency wI=zeta*PI*<df>.\n"
"                           Default zeta=%.1f.\n", GRT_GREENFN_F_ZETA); printf(
"\n"
"    -O<outdir>   Directorypath of output for saving.\n"
"\n"
"    -P<n>        Number of threads. Default use all cores.\n"
"\n"
"    -e           Compute the z-derivatives of kernel functions.\n"
"\n"
"    -h           Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt kernel -Mmod -D0.05/0 -F0/50/0.01 -C0.01 -OKERNEL\n"
"\n\n\n"
);
}

/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    // 先为个别参数设置非0初始值
    Ctrl->F.zeta = GRT_GREENFN_F_ZETA;

    int opt;
    while ((opt = getopt(argc, argv, ":M:D:F:C:O:P:eh")) != -1) 
    {
        switch(opt)
        {
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

            // -F<f1>/<f2>/<df>[+w<zeta>]
            case 'F':
                Ctrl->F.active = true;
                {
                    char *string = strdup(optarg);
                    char *token = strtok(string, "+");
                    real_t a1, a2, df;
                    a1 = a2 = df = 0;
                    if(3 != sscanf(token, "%lf/%lf/%lf", &a1, &a2, &df)){
                        GRTBadOptionError(F, "");
                    };
                    if(df <= 0){
                        GRTBadOptionError(F, "Can't set nonpositive df(%lf).", df);
                    }
                    if(a1 < 0.0 || a2 <= 0.0){
                        GRTBadOptionError(F, "Can't set nonpositive f1(%lf), f2(%lf).", a1, a2);
                    }
                    if(a1 > a2){
                        GRTBadOptionError(F, "f1(%lf) > f2(%lf).", a1, a2);
                    }

                    // 跳过零频
                    a1 = GRT_MAX(a1, df);
                    a2 = GRT_MAX(a1, a2);
                    
                    // 处理 + 号指令
                    token = strtok(NULL, "+");
                    if(token != NULL){
                        switch (token[0]){
                            case 'w':
                                if(1 != sscanf(token+1, "%lf", &Ctrl->F.zeta)){
                                    GRTBadOptionError(F, "");
                                }
                                if(Ctrl->F.zeta <= 0.0){
                                    GRTBadOptionError(F, "+%s need positive float, but get (%lf).", token, Ctrl->F.zeta);
                                }
                                break;
                            default:
                                GRTBadOptionError(F, "+%s is not supported.", token);
                                break;
                        }
                    }

                    Ctrl->F.nf = floor((a2-a1)/df) + 1;
                    // 至少要有两个频率点，才能在程序中用差分计算 df
                    if(Ctrl->F.nf < 2){
                        GRTBadOptionError(F, "Too few frequency points, only %zu.", Ctrl->F.nf);
                    }
                    Ctrl->F.freqs = (real_t *)calloc(Ctrl->F.nf, sizeof(real_t));
                    for(size_t i=0; i<Ctrl->F.nf; ++i){
                        Ctrl->F.freqs[i] = a1 + df*i;
                    }

                    GRT_SAFE_FREE_PTR(string);
                }
                break;

            // 输出路径 -Ooutput_dir
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.s_output_dir = strdup(optarg);
                break;

            // -C[<cmin>/<cmax>/]<dc>
            case 'C':
                Ctrl->C.active = true;
                {
                    real_t a1, a2, df;
                    a1 = a2 = df = 0;
                    int nscan = sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &df);
                    if( !(nscan == 1 || nscan == 3)){
                        GRTBadOptionError(C, "");
                    };
                    if(nscan == 1 && a1 <= 0.0){
                        GRTBadOptionError(C, "Can't set a single nonpositive value.");
                    }
                    if(nscan == 1){
                        Ctrl->C.dc = a1;
                        Ctrl->C.isset = false;
                    }
                    else{
                        Ctrl->C.isset = true;
                        Ctrl->C.cmin = a1;
                        Ctrl->C.cmax = a2;
                        Ctrl->C.dc = df;
                        if(Ctrl->C.cmin <= 0.0 || Ctrl->C.cmax <= 0.0){
                            GRTBadOptionError(C, "Can't set nonpositive cmin(%lf), cmax(%lf).", Ctrl->C.cmin, Ctrl->C.cmax);
                        }
                        if(Ctrl->C.cmin >= Ctrl->C.cmax){
                            GRTBadOptionError(C, "cmin(%lf) >= cmax(%lf).", Ctrl->C.cmin, Ctrl->C.cmax);
                        }

                        Ctrl->C.nc = floor((a2-a1)/df) + 1;
                        Ctrl->C.c_phases = (real_t *)calloc(Ctrl->C.nc, sizeof(real_t));
                        for(size_t i=0; i<Ctrl->C.nc; ++i){
                            Ctrl->C.c_phases[i] = a1 + df*i;
                        }
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
            
            // 是否计算位移空间导数
            case 'e':
                Ctrl->e.active = true;
                break;

            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, M);
    GRTCheckOptionActive(Ctrl, D);
    GRTCheckOptionActive(Ctrl, F);
    GRTCheckOptionActive(Ctrl, C);
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



/** 子模块主函数 */
int kernel_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    // 传入参数 
    getopt_from_command(Ctrl, argc, argv);

    // 读入模型文件
    if((Ctrl->M.mod1d = grt_read_mod1d_from_file(Ctrl->M.s_modelpath, Ctrl->D.depsrc, Ctrl->D.deprcv, true)) == NULL){
        exit(EXIT_FAILURE);
    }
    GRT_MODEL1D *mod1d = Ctrl->M.mod1d;

    // 最大最小速度
    real_t vmin, vmax;
    grt_get_mod1d_vmin_vmax(mod1d, &vmin, &vmax);
    if(! Ctrl->C.isset)
    {
        Ctrl->C.cmin = 0.8 * vmin;
        Ctrl->C.cmax = vmax;
        Ctrl->C.nc = floor((Ctrl->C.cmax - Ctrl->C.cmin)/Ctrl->C.dc) + 1;
        Ctrl->C.c_phases = (real_t *)calloc(Ctrl->C.nc, sizeof(real_t));
        for(size_t i=0; i<Ctrl->C.nc; ++i){
            Ctrl->C.c_phases[i] = Ctrl->C.cmin + Ctrl->C.dc*i;
        }
    }

    // 虚频率
    Ctrl->F.wI = Ctrl->F.zeta * PI * (Ctrl->F.freqs[1] - Ctrl->F.freqs[0]);

    // 保存目录
    char *s_output_dir = NULL;
    GRT_SAFE_ASPRINTF(&s_output_dir, "%s/%s_%s_%s", Ctrl->O.s_output_dir, Ctrl->M.s_modelname, Ctrl->D.s_depsrc, Ctrl->D.s_deprcv);
    GRTCheckMakeDir(s_output_dir);

    const real_t Rho = mod1d->Rho[mod1d->isrc]; // 震源区密度
    const real_t fac = 1.0/(4.0*PI*Rho);

    // 频率循环
    #pragma omp parallel for schedule(guided) default(shared) 
    for(size_t iw = 0; iw < Ctrl->F.nf; ++iw)
    {
        real_t freq = Ctrl->F.freqs[iw];
        real_t w = PI2*freq;
        cplx_t omega = w - I*Ctrl->F.wI;

        GRT_MODEL1D *local_mod1d = NULL;
    #ifdef _OPENMP 
        // 定义局部模型对象
        local_mod1d = grt_copy_mod1d(mod1d);
    #else 
        local_mod1d = mod1d;
    #endif

        // 将 omega 计入模型结构体
        local_mod1d->omega = omega;

        grt_attenuate_mod1d(local_mod1d, omega);

        // 为当前频率创建波数积分记录文件
        FILE *fstats = NULL;
        char *fname = NULL;
        GRT_SAFE_ASPRINTF(&fname, "%s/C_%04zu_%.5e", s_output_dir, iw, freq);
        fstats = fopen(fname, "wb");
        GRT_SAFE_FREE_PTR(fname);

        // 不同震源不同阶数的核函数 F(k, w) 
        cplxChnlGrid QWV = {0};
        cplxChnlGrid QWV_uiz = {0};

        // 相速度循环
        for(size_t ic = 0; ic < Ctrl->C.nc; ++ic)
        {
            real_t cc = Ctrl->C.c_phases[ic];
            real_t k = w/cc;
            grt_kernel(local_mod1d, k, QWV, Ctrl->e.active, QWV_uiz);

            // 系数
            GRT_LOOP_ChnlGrid(im, c){
                cplx_t tmp = - fac/(omega*omega);
                QWV[im][c] *= tmp;
                if(Ctrl->e.active){
                    QWV_uiz[im][c] *= tmp;
                }
            }

            // 记录积分核函数，注意这里传入的是相速度
            grt_write_stats(fstats, cc, (Ctrl->e.active)? QWV_uiz : QWV);
        }

        fclose(fstats);

    #ifdef _OPENMP
        grt_free_mod1d(local_mod1d);
    #endif

    }

    printf("Number of frequencies: %zu\n", Ctrl->F.nf);
    printf("Number of phase velocities: %zu\n", Ctrl->C.nc);

    GRT_SAFE_FREE_PTR(s_output_dir);
    free_Ctrl(Ctrl);

    return EXIT_SUCCESS;
}