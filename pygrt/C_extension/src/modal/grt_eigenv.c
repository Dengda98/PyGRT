/**
 * @file   grt_eigenv.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-05
 * 
 *     求解Rayleigh波和Love波的相速度频散
 * 
 */
#include "grt.h"

#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/common/search.h"
#include "grt/common/util.h"

#include "grt/modal/modal_def.h"
#include "grt/modal/modal_util.h"
#include "grt/modal/root.h"


/** 该子模块的参数控制结构体 */
typedef struct {
    /** 输入模型 */
    struct {
        bool active;
        char *s_modelpath;        ///< 模型路径
        const char *s_modelname;  ///< 模型名称
        GRT_MODEL1D *mod1d;         ///< 模型结构体指针
    } M;
    /* 相速度搜索范围 */
    struct {
        bool active;
        real_t cmin;
        real_t cmax;
    } V;
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
    } X;
    /* 零点阈值 */
    struct {
        bool active;
        real_t satol;
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
"[grt eigenval]\n\n"
"    Compute eigenvalue of surface waves (dispersion).\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt eigenval <syn_dir>/<name>\n"
"    [WAITING TO FINISH]\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    // 先为个别参数设置非0初始值
    Ctrl->T.satol = 1e-3;

    // 以下浮点数不易过大或过小
    Ctrl->T.cgap = 1e-8;
    Ctrl->T.rtol = 1e-6;
    Ctrl->T.vgap = 1e-7;

    Ctrl->V.cmin = -1.0;
    Ctrl->V.cmax = -1.0;

    int opt;
    while ((opt = getopt(argc, argv, ":M:V:C:S:N::F:X:T:P:sh")) != -1) {
        switch(opt) {
            // 模型路径，其中每行分别为 
            //      厚度(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  Qp   Qs
            // 互相用空格隔开即可
            case 'M':
                Ctrl->M.active = true;
                Ctrl->M.s_modelpath = strdup(optarg);
                Ctrl->M.s_modelname = grt_get_basename(Ctrl->M.s_modelpath);
                break;

            // -V<cmin>/<cmax>
            case 'V':
                Ctrl->V.active = true;
                if(2 != sscanf(optarg, "%lf/%lf", &Ctrl->V.cmin, &Ctrl->V.cmax)){
                    GRTBadOptionError(V, "");
                }
                if(Ctrl->V.cmin <= 0.0 || Ctrl->V.cmax <= 0.0){
                    GRTBadOptionError(V, "Can't set nonpositive cmin(%lf), cmax(%lf).", Ctrl->V.cmin, Ctrl->V.cmax);
                }
                if(Ctrl->V.cmin >= Ctrl->V.cmax){
                    GRTBadOptionError(V, "cmin(%lf) >= cmax(%lf).", Ctrl->V.cmin, Ctrl->V.cmax);
                }
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
            
            // -X<freq>
            case 'X':
                Ctrl->X.active = true;
                {
                    real_t fx = -1.0;
                    if(1 != sscanf(optarg, "%lf", &fx)){
                        GRTBadOptionError(X, "");
                    }
                    if(fx <= 0){
                        GRTBadOptionError(X, "Can't set nonpositive freq(%lf).", fx);
                    }

                    // 如果freqs已经在-F中申请了内存，则先释放
                    if(Ctrl->F.freqs != NULL)  free(Ctrl->F.freqs);
                    // 仅申请单一频率，用于输出久期函数值
                    Ctrl->F.nf = 1;
                    Ctrl->F.freqs = (real_t*)calloc(Ctrl->F.nf, sizeof(real_t));
                    Ctrl->F.freqs[0] = fx;
                }
                break;

            // -T+t<tol>+c<cgap>+r<tol2>+v<vgap>
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
    GRT_MODEL1D *mod1d = Ctrl->M.mod1d;

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
    eigmet->cmin = Ctrl->V.cmin;
    eigmet->cmax = Ctrl->V.cmax;
    eigmet->satol = Ctrl->T.satol;
    eigmet->cgap = Ctrl->T.cgap;
    eigmet->rtol = Ctrl->T.rtol;
    eigmet->vgap = Ctrl->T.vgap;
    // 存储频散值的结构体
    eigmet->eigv = (EIGENV *)calloc(eigmet->nf, sizeof(EIGENV));

    // 自动取 S 波速度范围为搜索范围
    if(eigmet->cmin < 0.0 || eigmet->cmax < 0.0){
        real_t vbmin = mod1d->Vb[0];  // 最小 S 波速度
        for(size_t i = 0; i < mod1d->n; ++i){
            vbmin = GRT_MIN(vbmin, mod1d->Vb[i]);
        }
        real_t vbn = mod1d->Vb[mod1d->n-1];  // 半空间的 S 波速度
        
        if(vbmin >= vbn){
            GRTRaiseError("minimum Vs(%f) >= halfspace Vs(%f), you should set -C manually.", vbmin, vbn);
        }
        if(vbmin == 0.0){
            GRTRaiseError("minimum Vs == 0.0 , you should set -C manually.");
        }

        eigmet->cmin = vbmin;
        eigmet->cmax = vbn;
        if(eigmet->wtype == GRT_DISPERSION_RAYL)  eigmet->cmin *= 0.9;
    }
    
    // 寻找久期函数零点
    grt_get_secular_roots(mod1d, eigmet, !Ctrl->s.active);

    if(Ctrl->X.active) goto FINISH;

    // 生成总命令
    char *full_command;
    {
        char *tmp=NULL;
        GRT_SAFE_ASPRINTF(&tmp, GRT_MAIN_COMMAND);
        for(int i=0; i<argc; ++i){
            char *tmp1=NULL;
            GRT_SAFE_ASPRINTF(&tmp1, "%s %s", tmp, argv[i]);
            GRT_SAFE_FREE_PTR(tmp);
            tmp = tmp1;
        }
        
        // 对 -N 特殊处理，因为 -N 可能没有传入，需要根据程序自动设置
        // 这里在最后再加上一个 -N ，功能上会覆盖之前可能设置的 -N
        // 相当于之前设置的最大阶数有冗余

        // 最大零点数
        size_t max_nroots = 0;
        for(size_t iw=0; iw < Ctrl->F.nf; ++iw){
            max_nroots = GRT_MAX(max_nroots, eigmet->eigv[iw].n);
        }
        char *tmp1=NULL;
        GRT_SAFE_ASPRINTF(&tmp1, "%s -N%zu", tmp, max_nroots-1);
        GRT_SAFE_FREE_PTR(tmp);
        tmp = tmp1;

        full_command = tmp;
    }

    // 输出频散结果
    grt_output_cdisp(Ctrl->C.s_phasepath, full_command, Ctrl->M.s_modelpath, eigmet);
    GRT_SAFE_FREE_PTR(full_command);
    
    grt_free_eigenv_info(eigmet);

FINISH:
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}




