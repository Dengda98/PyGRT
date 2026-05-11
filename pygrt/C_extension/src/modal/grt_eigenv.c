/**
 * @file   grt_eigenv.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-05
 * 
 *     求解Rayleigh波和Love波的相速度频散
 * 
 */

#include "grt.h"

#define GRT_EIGENV_SATOL  1e-3
#define GRT_EIGENV_CGAP   1e-7
#define GRT_EIGENV_RTOL   1e-4
#define GRT_EIGENV_VGAP   1e-7


/** 该子模块的参数控制结构体 */
typedef struct {
    /** 输入模型 */
    struct {
        bool active;
        char *s_modelpath;        ///< 模型路径
        const char *s_modelname;  ///< 模型名称
        MODEL1D *mod1d;         ///< 模型结构体指针
    } M;
    /* 相速度频散结果输出路径 */
    struct {
        bool active;
        char *s_phasepath;
    } C;
    /* 频散类型 */
    struct {
        bool active;
        DISPER_TYPE wtype;
    } S;
    /* 最大阶数 */
    struct {
        bool active;
        size_t nmode;
    } N;
    /* 频率范围 */
    struct {
        bool active;
        real_t *freqs;
        size_t nf;
    } F;
    /** 输出久期函数 */
    struct {
        bool active;
        bool manual_crange;
        real_t cmin;
        real_t cmax;
    } X;
    /* 零点阈值 */
    struct {
        bool active;
        real_t satol;
        real_t uniform_dc;
        real_t cgap;
        real_t rtol;
        real_t vgap;
    } T;    
    /** 多线程 */
    struct {
        bool active;
        int nthreads; ///< 线程数
    } P;
    /* 静默输出 */
    struct {
        bool active;
    } s;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    // M
    GRT_SAFE_FREE_PTR(Ctrl->M.s_modelpath);
    grt_free_mod1d(Ctrl->M.mod1d);
    // C
    GRT_SAFE_FREE_PTR(Ctrl->C.s_phasepath);
    // F
    GRT_SAFE_FREE_PTR(Ctrl->F.freqs);

    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt eigenv] %s\n\n", GRT_VERSION);printf(
"    Compute eigenvalues (dispersion) of surface waves.\n"
"\n\n"
"Usage:\n"
"---------------------------------------------------------------------\n"
"    grt eigenv -M<model> -C<path> -SR|L -F<f1>[/<f2>/<df>][+p]\n"
"         [-N[<max_order>]] [-X<freq>[+c<cmin>/<cmax>]] [-T+<params>] \n"
"         [-P<nthreads>] [-s] [-h]\n"
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
"         The 1st column H is the layer thickness. \n"
"         For compatibility, if the \"thickness\" of the first layer \n"
"         set 0.0, then the first column will be the layer top depth, \n"
"         which must be in an ascending order.\n"
"         The number of layers are unlimited.\n"
"\n"
"    -C<path>    Dispersion result file (.nc format)\n"
"\n"
"    -SR|L       Surface wave type:\n"
"                + R: Rayleigh wave\n"
"                + L: Love wave\n"
"\n"
"    -F<f1>[/<f2>/<df>][+p]\n"
"                Set the frequency points.\n"
"                <f1>: start frequency (Hz)\n"
"                <f2>: end frequency (Hz)\n"
"                <df>: frequency interval (Hz)\n"
"                + If set -F<f1>, means set only one point. \n"
"                + If add +p at the end, <f1>, <f2>\n"
"                  and <df> will be period (sec).\n"
"\n"
"    -N[<max_order>]\n"
"                Set the max order of roots at each frequency,\n"
"                0 means fundamental mode.\n"
"                + If -N is not given, equal to -N0.\n"
"                + If set an empty -N, means find all roots\n"
"                  as much as possible.\n"
"\n"
"    -X<freq>[+c<cmin>/<cmax>]\n"
"                Debug option. You can set -X to output the\n"
"                secular function at frequency <freq> within \n"
"                phase velocity [<cmin>, <cmax>], and\n"
"                at this time, -F will be omitted.\n"
"\n"
"    -T+<params> set some control parameters in root searching,\n"
"                e.g. -T+t3+c7\n"
"                + t<tol>: tolerance (10^{-<tol>}) for adaptive searching,\n"
"                          default is %.1f\n", - log10(GRT_EIGENV_SATOL)); printf(
"                + c<cgap>: min gap (10^{-<cgap>} km/s) between different\n"
"                           roots at the same frequency, default is %.1f\n", - log10(GRT_EIGENV_CGAP)); printf(
"                + r<thrd>: amplitude threshold (10^{-<thrd>}) to identify\n"
"                           a root of secular function, default is %.1f\n", - log10(GRT_EIGENV_RTOL)); printf(
"                + v<vgap>: min gap (10^{-<vgap>} km/s) between root\n"
"                           and layer velocity, default is %.1f\n", - log10(GRT_EIGENV_VGAP)); printf(
"                + u<dc>: searching step (10^{-<dc>} km/s) for\n"
"                         fixed-interval searching, this is only for test.\n"); printf(
"\n"
"    -P<n>        Number of threads. Default use all cores.\n"
"\n"
"    -s           Silence all outputs.\n"
"\n"
"    -h           Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt eigenv -Mmod -F0/100/0.5 -SR -N -Cphase_R.nc\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    // 先为个别参数设置非0初始值
    Ctrl->T.satol = GRT_EIGENV_SATOL;

    // 以下浮点数不易过大或过小
    Ctrl->T.cgap = GRT_EIGENV_CGAP;
    Ctrl->T.rtol = GRT_EIGENV_RTOL;
    Ctrl->T.vgap = GRT_EIGENV_VGAP;

    int opt;
    while ((opt = getopt(argc, argv, ":M:C:S:N::F:X:T:P:sh")) != -1) {
        switch(opt) {
            // 模型路径，其中每行分别为 
            //      厚度(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  Qp   Qs
            // 互相用空格隔开即可
            case 'M':
                Ctrl->M.active = true;
                Ctrl->M.s_modelpath = strdup(optarg);
                Ctrl->M.s_modelname = grt_get_basename(Ctrl->M.s_modelpath);
                break;

            // -C<path>
            case 'C':
                Ctrl->C.active = true;
                Ctrl->C.s_phasepath = strdup(optarg);
                break;
            
            // -S[R|L]
            case 'S':
                Ctrl->S.active = true;
                {
                    char x;
                    if(1 != sscanf(optarg, "%c", &x)){
                        GRTBadOptionError(S, "");
                    }
                    if(x == 'R'){
                        Ctrl->S.wtype = GRT_DISPERSION_RAYL;
                    }
                    else if(x == 'L'){
                        Ctrl->S.wtype = GRT_DISPERSION_LOVE;
                    }
                    else {
                        GRTBadOptionError(S, "only support -SR and -SL.");
                    }
                }
                break;

            // -N<max_mode>
            case 'N':
                Ctrl->N.active = true;
                // 没有参数，直接跳过
                if(optarg==NULL)  break;
                {
                    ssize_t tmp = 0;
                    if(1 != sscanf(optarg, "%zd", &tmp)){
                        GRTBadOptionError(N, "");
                    }
                    if(tmp < 0){
                        GRTBadOptionError(N, "Can't set negative nmode(%zd).", tmp);
                    }
                    Ctrl->N.nmode = tmp;
                }
                
                Ctrl->N.nmode++;
                break;

            // -Ff1[/f2/df][+p]
            case 'F':
                Ctrl->F.active = true;
                {
                    char *string = strdup(optarg);
                    char *token = strtok(string, "+");
                    real_t a1, a2, df;
                    a1 = a2 = df = 0;
                    int nscan = sscanf(token, "%lf/%lf/%lf", &a1, &a2, &df);
                    if( !(nscan == 1 || nscan == 3)){
                        GRTBadOptionError(F, "");
                    };
                    if(nscan == 1 && a1 <= 0.0){
                        GRTBadOptionError(F, "Can't set a single nonpositive value.");
                    }
                    if(nscan == 1){
                        a2 = a1;
                        df = 1;
                    }
                    else if(nscan == 3){
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
                    }
                    
                    bool isperiod = false;
                    // 处理 + 号指令
                    token = strtok(NULL, "+");
                    if(token != NULL){
                        switch (token[0]){
                            case 'p':
                                isperiod = true;
                                break;
                            default:
                                GRTBadOptionError(F, "+%s is not supported.", token);
                                break;
                        }
                    }

                    Ctrl->F.nf = floor((a2-a1)/df) + 1;
                    Ctrl->F.freqs = (real_t*)calloc(Ctrl->F.nf, sizeof(real_t));
                    for(size_t i = 0; i < Ctrl->F.nf; ++i){
                        if(isperiod){
                            Ctrl->F.freqs[Ctrl->F.nf-1 - i] = 1.0/(a1 + df*i);
                        } else {
                            Ctrl->F.freqs[i] = a1 + df*i;
                        }
                    }

                    GRT_SAFE_FREE_PTR(string);
                }
                break;
            
            // -X<freq>[+c<cmin>/<cmax>]
            case 'X':
                Ctrl->X.active = true;
                {
                    char *string = strdup(optarg);
                    char *token = strtok(string, "+");
                    real_t fx = -1.0;
                    if(1 != sscanf(token, "%lf", &fx)){
                        GRTBadOptionError(X, "");
                    }
                    if(fx <= 0){
                        GRTBadOptionError(X, "Can't set nonpositive freq(%lf).", fx);
                    }

                    // 处理 + 号指令
                    token = strtok(NULL, "+");
                    if(token != NULL){
                        switch (token[0]){
                            case 'c':
                                Ctrl->X.manual_crange = true;
                                if(2 != sscanf(token+1, "%lf/%lf", &Ctrl->X.cmin, &Ctrl->X.cmax)){
                                    GRTBadOptionError(X, "");
                                }
                                if(Ctrl->X.cmin < 0.0 || Ctrl->X.cmax < 0.0){
                                    GRTBadOptionError(X, "Can't set negative velocity.");
                                }
                                if(Ctrl->X.cmin >= Ctrl->X.cmax){
                                    GRTBadOptionError(X, "requires cmax < cmax.");
                                }
                                break;
                            default:
                                GRTBadOptionError(X, "+%s is not supported.", token);
                                break;
                        }
                    }

                    // 如果freqs已经在-F中申请了内存，则先释放
                    if(Ctrl->F.freqs != NULL)  free(Ctrl->F.freqs);
                    // 仅申请单一频率，用于输出久期函数值
                    Ctrl->F.nf = 1;
                    Ctrl->F.freqs = (real_t*)calloc(Ctrl->F.nf, sizeof(real_t));
                    Ctrl->F.freqs[0] = fx;

                    GRT_SAFE_FREE_PTR(string);
                }
                break;

            // -T+t<tol>+c<cgap>+r<thrd>+v<vgap>+u<dc>
            case 'T':
                Ctrl->T.active = true;
                {
                    char *line = strdup(optarg);
                    char *token = strtok(line, "+");
                    real_t p=0;
                    while(token != NULL){
                        switch(token[0]){
                            case 't':
                                if(1 != sscanf(token+1, "%lf", &p)){
                                    GRTBadOptionError(T+t, "");
                                }
                                Ctrl->T.satol = pow(10, -p);
                                break;
                            case 'c':
                                if(1 != sscanf(token+1, "%lf", &p)){
                                    GRTBadOptionError(T+k, "");
                                }
                                Ctrl->T.cgap = pow(10, -p);
                                break;
                            case 'r':
                                if(1 != sscanf(token+1, "%lf", &p)){
                                    GRTBadOptionError(T+r, "");
                                }
                                Ctrl->T.rtol = pow(10, -p);
                                break;
                            case 'v':
                                if(1 != sscanf(token+1, "%lf", &p)){
                                    GRTBadOptionError(T+v, "");
                                }
                                Ctrl->T.vgap = pow(10, -p);
                                break;
                            case 'u':
                                if(1 != sscanf(token+1, "%lf", &p)){
                                    GRTBadOptionError(T+v, "");
                                }
                                Ctrl->T.uniform_dc = pow(10, -p);
                                break;
                            default:
                                GRTBadOptionError(T, "-T+%s is not supported.", token);
                                break;
                        }
                        token = strtok(NULL, "+");
                    } // END -T loop

                    GRT_SAFE_FREE_PTR(line);
                }
                break;

            // 多线程数 -Pnthreads
            case 'P':
                Ctrl->P.active = true;
                if(1 != sscanf(optarg, "%d", &Ctrl->P.nthreads)){
                    GRTBadOptionError(P, "");
                };
                if(Ctrl->P.nthreads <= 0){
                    GRTBadOptionError(P, "Nonpositive value in -P is not supported.");
                }
                grt_set_num_threads(Ctrl->P.nthreads);
                break;

            // 静默输出
            case 's':
                Ctrl->s.active = true;
                break;

            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, M);
    GRTCheckOptionActive(Ctrl, S);
    if( ! Ctrl->X.active){
        GRTCheckOptionActive(Ctrl, C);
    }
    if(Ctrl->F.active + Ctrl->X.active != 1){
        GRTRaiseError("Need set and only set one of -F and -X. Use '-h' for help.\n");
    }

    // 如果设置了 -X 输出久期函数，就实行静默输出信息，防止和 -X 输出撞车
    if(Ctrl->X.active){
        Ctrl->s.active = false;
    }

    // 默认只找基阶
    // 如果注释这句，则根据找零点函数中的限制，
    // 会寻找当前频段内所有根
    if(! Ctrl->N.active){
        Ctrl->N.nmode = 1;
    }

}


/** 子模块主函数 */
int eigenv_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    // 传入参数 
    getopt_from_command(Ctrl, argc, argv);

    // 读入模型文件（暂设置为不支持液体层）
    if((Ctrl->M.mod1d = grt_read_mod1d_from_file(Ctrl->M.s_modelpath, -1.0, -1.0, true)) ==NULL){
        exit(EXIT_FAILURE);
    }
    MODEL1D *mod1d = Ctrl->M.mod1d;

    // 目前边界条件暂有限制
    if( ! (mod1d->topbound==GRT_BOUND_FREE && mod1d->botbound==GRT_BOUND_HALFSPACE)){
        GRTRaiseError("Currently only support top-free and bottom-halfspace in Surface wave problem.");
    }

    // 将信息转入结构体
    EIGENV_INFO *eigmet = (EIGENV_INFO *)calloc(1, sizeof(EIGENV_INFO));
    eigmet->nf = Ctrl->F.nf;
    eigmet->freqs = (real_t *)calloc(eigmet->nf, sizeof(real_t));
    memcpy(eigmet->freqs, Ctrl->F.freqs, sizeof(real_t)*eigmet->nf);
    eigmet->nmode = Ctrl->N.nmode;
    eigmet->wtype = Ctrl->S.wtype;
    eigmet->print_sec = Ctrl->X.active;
    eigmet->satol = Ctrl->T.satol;
    eigmet->uniform_dc = Ctrl->T.uniform_dc;
    eigmet->cgap = Ctrl->T.cgap;
    eigmet->rtol = Ctrl->T.rtol;
    eigmet->vgap = Ctrl->T.vgap;
    // 存储频散值的结构体
    eigmet->eigv = (EIGENV *)calloc(eigmet->nf, sizeof(EIGENV));

    // 相速度搜索范围
    if(Ctrl->X.manual_crange)
    {
        eigmet->manual_crange = true;
        eigmet->cmin = Ctrl->X.cmin;
        eigmet->cmax = Ctrl->X.cmax;
    }
    else
    {
        // 取 S 波速度范围为搜索范围
        real_t vbmin = 9.9e30;  // 最小 S 波速度
        for(size_t i = 0; i < mod1d->n; ++i){
            real_t va = mod1d->Va[i];
            real_t vb = mod1d->Vb[i];
            vbmin = GRT_MIN(vbmin, (vb > 0.0)? vb : va);
        }
        real_t vbn = mod1d->Vb[mod1d->n-1];  // 半空间速度
        if(vbn == 0.0) vbn = mod1d->Va[mod1d->n-1];

        // 如果存在上层固体、下层液体的情况，设置 vbmin 为极小非零值
        for(size_t i = 0; i < mod1d->n-1; ++i){
            if( ! mod1d->isLiquid[i] && mod1d->isLiquid[i+1]){
                vbmin *= 1e-3;
                break;
            }
        }

        if(vbmin >= vbn){
            GRTRaiseError("minimum velocity (vbmin=%f) >= halfspace velocity (vbn=%f).", vbmin, vbn);
        }

        eigmet->cmin = vbmin;
        eigmet->cmax = vbn;
    }
    
    // 寻找久期函数零点
    grt_get_secular_roots(mod1d, eigmet, !Ctrl->s.active);

    if(! Ctrl->s.active) printf("# Number of evaluation: %zu\n", eigmet->neval);

    if(Ctrl->X.active) goto FINISH;

    // 生成总命令
    char *full_command = NULL;
    GRT_SAFE_ASPRINTF(&full_command, GRT_MAIN_COMMAND);
    for(int i = 0; i < argc; ++i){
        char *tmp1=NULL;
        GRT_SAFE_ASPRINTF(&tmp1, "%s %s", full_command, argv[i]);
        GRT_SAFE_FREE_PTR(full_command);
        full_command = tmp1;
    }

    // 输出频散结果
    grt_output_cdisp(Ctrl->C.s_phasepath, full_command, Ctrl->M.s_modelpath, eigmet);
    GRT_SAFE_FREE_PTR(full_command);
    
    grt_free_eigenv_info(eigmet);

FINISH:
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}




