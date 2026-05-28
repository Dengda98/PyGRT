/**
 * @file   grt_disp2asc.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-02
 * 
 *     将 eigenv 模块输出的相速度频散结果或 eigenfn 模块输出的群速度频散结果
 *     转为文本格式以供阅读
 * 
 */

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    /* 相速度频散结果输出路径 */
    struct {
        bool active;
        char *s_phasepath;
    } C;
    /* 群速度频散结果输出路径 */
    struct {
        bool active;
        char *s_grouppath;
    } U;
    /* 阶数范围 */
    struct {
        bool active;
        size_t *modes;
        size_t nmode;
    } N;
    /* 频率范围 */
    struct {
        bool active;
        real_t *freqs;
        size_t nf;
        bool def_range;  ///< 定义了搜索范围而不是离散的一些点
        bool io_period;  ///< 以周期的形式输入输出
    } F;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    // C
    GRT_SAFE_FREE_PTR(Ctrl->C.s_phasepath);
    // U
    GRT_SAFE_FREE_PTR(Ctrl->U.s_grouppath);
    // N
    GRT_SAFE_FREE_PTR(Ctrl->N.modes);
    // F
    GRT_SAFE_FREE_PTR(Ctrl->F.freqs);
    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt disp2asc] %s\n\n", GRT_VERSION);printf(
"    Convert a nc-format dispersion result from `eigenv` module \n"
"    or `eigenfn` module into an ASCII file, write to \n"
"    standard output (ignore some vars).\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt disp2asc -C<path> -U<path> [-F<f1>[/<f2>][/<df>][+p]] \n"
"         [-N[<n1>][/<n2>][/<dn>]] [-h]\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -C<path>    Input phase-velocity dispersion result file \n"
"                from module `eigenv` (.nc format)\n"
"\n"
"    -U<path>    Input group-velocity dispersion result file \n"
"                from module `eigenfn` (.nc format)\n"
"\n"
"    -F<f1>[/<f2>][/<df>][+p]\n"
"                Select the frequency range from the input file.\n"
"                <f1>: start frequency (Hz)\n"
"                <f2>: end frequency (Hz)\n"
"                <df>: frequency interval (Hz)\n"
"                + If set -F<f1>, means set only one point. \n"
"                + If set -F<f1>/<f2>, means set min/max range. \n"
"                + If add +p at the end, <f1>, <f2>\n"
"                  and <df> will be period (sec).\n"
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
"    -h           Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt disp2asc -Cphase_R.nc -N > phase_R.txt\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    int opt;
    while ((opt = getopt(argc, argv, ":C:U:F:N::h")) != -1) {
        switch (opt) {
            // -C<path>
            case 'C':
                Ctrl->C.active = true;
                Ctrl->C.s_phasepath = strdup(optarg);
                break;

            // -U<path>
            case 'U':
                Ctrl->U.active = true;
                Ctrl->U.s_grouppath = strdup(optarg);
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

            // 频率值 -Ff1[/f2][/df][+p]  可使用
            case 'F':
                Ctrl->F.active = true;
                {
                    char *string = strdup(optarg);
                    char *token = strtok(string, "+");
                    real_t a1, a2, df;
                    a1 = a2 = 0;
                    df = 1;
                    int nscan = sscanf(token, "%lf/%lf/%lf", &a1, &a2, &df);
                    if( !(nscan == 1 || nscan == 2 || nscan == 3)){
                        GRTBadOptionError(F, "");
                    };
                    if(nscan == 1 && a1 <= 0.0){
                        GRTBadOptionError(F, "Can't set a single nonpositive value.");
                    }
                    if(nscan == 1){
                        a2 = a1;
                        df = 1;
                    }
                    else if(nscan >= 2){
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
                        if(nscan == 3){
                            a1 = GRT_MAX(a1, df);
                            a2 = GRT_MAX(a1, a2);
                        }
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

                    if(nscan == 3 || nscan == 1){
                        Ctrl->F.nf = floor((a2-a1)/df) + 1;
                        Ctrl->F.freqs = (real_t*)calloc(Ctrl->F.nf, sizeof(real_t));
                        for(size_t i=0; i<Ctrl->F.nf; ++i){
                            if(isperiod){
                                Ctrl->F.freqs[Ctrl->F.nf-1 - i] = 1.0/(a1 + df*i);
                            } else {
                                Ctrl->F.freqs[i] = a1 + df*i;
                            }
                        }
                    }
                    else if(nscan == 2){
                        // 只约束最大值和最小值，内存申请由读取频散文件的过程中进行
                        Ctrl->F.nf = 2;
                        Ctrl->F.freqs = (real_t*)calloc(Ctrl->F.nf, sizeof(real_t));
                        Ctrl->F.freqs[0] = a1;
                        Ctrl->F.freqs[1] = a2;
                        Ctrl->F.def_range = true;
                    }
                    
                    Ctrl->F.io_period = isperiod;

                    GRT_SAFE_FREE_PTR(string);
                }
                break;

            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 检查必选项有没有设置
    GRTCheckOptionSet(argc > 1);

    // -C 和 -U 必须选择其中一个
    if(Ctrl->C.active + Ctrl->U.active != 1){
        GRTRaiseError("Need set and only set one of -C and -U. Use '-h' for help.\n");
    }

    // -N 的默认项，仅处理基阶
    if( ! Ctrl->N.active ){
        Ctrl->N.nmode = 1;
        Ctrl->N.modes = (size_t *)calloc(Ctrl->N.nmode, sizeof(size_t));
        Ctrl->N.modes[0] = 0;
    }
}


/** 子模块主函数 */
int disp2asc_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    // 传入参数 
    getopt_from_command(Ctrl, argc, argv);

    // 读取频散
    EIGENV_INFO *eigmet = (EIGENV_INFO *)calloc(1, sizeof(EIGENV_INFO));

    if(Ctrl->C.active){
        char *modelpath = NULL;
        grt_read_dispersion(Ctrl->C.s_phasepath, eigmet, &modelpath);
        GRT_SAFE_FREE_PTR(modelpath);
    } else if(Ctrl->U.active) {
        grt_read_dispersion(Ctrl->U.s_grouppath, eigmet, NULL);
    } else {
        GRTRaiseError("Wrong Execution.");
    }
    

    // 根据命令行参数确定出所需的部分频散信息
    EIGENFN_INFO *eigfnmet = (EIGENFN_INFO *)calloc(1, sizeof(EIGENFN_INFO));
    grt_filter_eigenfn_info(
        Ctrl->F.nf, Ctrl->F.freqs, Ctrl->F.def_range,
        Ctrl->N.nmode, Ctrl->N.modes, eigmet, eigfnmet);
    
    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        real_t *c_roots = eigfnmet->eigv[iw].c_roots;
        real_t *u_roots = eigfnmet->eigv[iw].u_roots;
        size_t *c_roots_iref = eigfnmet->eigv[iw].c_roots_iref;
        size_t cnum = eigfnmet->eigv[iw].n;

        // 打印相速度
        if(Ctrl->C.active){
            for(size_t ik = 0; ik < cnum; ++ik){
                printf("%15.5e %10zu %20.8e %10zu\n", eigfnmet->freqs[iw], eigfnmet->modes[ik], c_roots[ik], c_roots_iref[ik]);
            }
        // 打印群速度
        } else if(Ctrl->U.active){
            for(size_t ik = 0; ik < cnum; ++ik){
                printf("%15.5e %10zu %20.8e %20.8e\n", eigfnmet->freqs[iw], eigfnmet->modes[ik], c_roots[ik], u_roots[ik]);
            }
        }
        
    }

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}