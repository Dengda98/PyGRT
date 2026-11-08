/**
 * @file   grt_eigenfn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-07
 * 
 *     根据计算的相速度频散(本征值)，求解本征函数（包括位移和应力）
 * 
 */

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    /* 输入频散文件 */
    struct {
        bool active;
        char *s_filepath;
        bool readin_period;
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
        real_t *freqs;
        size_t nf;
        bool def_range;  ///< 定义了搜索范围而不是离散的一些点
        bool io_period;  ///< 以周期的形式输入输出
    } F;
    /* 群速度 */
    struct {
        bool active;
        char *s_udisp_filepath;
    } U;
    /* 敏感核 */
    struct {
        bool active;
        char *s_cpar_filepath;    ///< 相速度对模型参数的偏导数（敏感核）
        char *s_upar_filepath;    ///< 群速度对模型参数的偏导数（敏感核）
        char *s_egyint_filepath;  ///< 输出能量积分
        real_t cpar_dz;           ///< 计算敏感核时对模型深度进行采样
        real_t freq_shift;        ///< 群速度核的正频率偏移量，Hz，0 表示使用相邻频点差分
    } K;
    /* 本征函数 */
    struct {
        bool active;
        char *s_egn_filepath;
        real_t *zs;
        size_t nz;
    } W;

    bool calc_egyint;   ///< 是否需要计算能量积分

} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    // C
    GRT_SAFE_FREE_PTR(Ctrl->C.s_filepath);
    // N
    GRT_SAFE_FREE_PTR(Ctrl->N.modes);
    // F
    GRT_SAFE_FREE_PTR(Ctrl->F.freqs);
    // U
    GRT_SAFE_FREE_PTR(Ctrl->U.s_udisp_filepath);
    // K
    GRT_SAFE_FREE_PTR(Ctrl->K.s_cpar_filepath);
    GRT_SAFE_FREE_PTR(Ctrl->K.s_upar_filepath);
    GRT_SAFE_FREE_PTR(Ctrl->K.s_egyint_filepath);
    // W
    GRT_SAFE_FREE_PTR(Ctrl->W.s_egn_filepath);
    GRT_SAFE_FREE_PTR(Ctrl->W.zs);

    // GRT_SAFE_FREE_PTR_ARRAY(Ctrl->cdisp, Ctrl->F.nf);
    // GRT_SAFE_FREE_PTR_ARRAY(Ctrl->cdisp_iref, Ctrl->F.nf);
    // GRT_SAFE_FREE_PTR(Ctrl->cdisp_num);

    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt eigenfn] %s\n\n", GRT_VERSION);printf(
"    Compute eigenfunctions, group velocity, energy integrals and \n"
"    dispersion sensitivity based on eigenvalues of surface waves.\n"
"\n"
"    The 1-D model is stored in the dispersion result from \n"
"    module `eigenv`, therefore this module will automatically\n"
"    read the model from the nc file.\n"
"\n\n"
"Usage:\n"
"---------------------------------------------------------------------\n"
"    grt eigenfn -C<path> [-F<f1>[/<f2>][/<df>][+p]]\n"
"         [-N[<n1>][/<n2>][/<dn>]] [-W<path>[+z<z1>[/<z2>/<dz>]]]\n"
"         [-U<path>] [-K[+c<path>][+u<path>][+x<path>][+z<dz>][+s<shift>]] [-h]\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -C<path>    Input dispersion result file from module `eigenv` \n"
"                (.nc format)\n"
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
"    -W<path>[+z<z1>[/<z2>/<dz>]]\n"
"                Set the depth points and output eigenfunctions\n"
"                result file (.nc format). \n"
"                <z1>: start depth (km)\n"
"                <z2>: end depth (km)\n"
"                <dz>: depth interval (km)\n"
"                If set +z<z1>, means set only one point. \n"
"                If not set +z, means set the layer depths. \n"
"\n"
"    -U<path>    Output group velocity (.nc format)\n"
"\n"
"    -K[+c<path>][+u<path>][+x<path>][+z<dz>][+s<shift>]\n"
"                Output some quantities related to sensitivity,\n"
"                e.g. -K+ccsens.nc+uusens.nc\n"
"                + c<path>: phase-velocity dispersion sensitivity (.nc format).\n"
"                + u<path>: group-velocity dispersion sensitivity (.nc format).\n"
"                + z<dz>:   divide depth with interval dz (km),\n"
"                           equivalent to performing stepwise \n"
"                           sampling on the model before computing \n"
"                           the sensitivity.\n"
"                + x<path>: energy integrals (.nc format)\n"
"                + s<shift>: positive frequency shift (Hz) for group kernels,\n"
"                            requires +u. Recompute at shifted frequencies\n"
"                            for finite differences. If omitted, use existing\n"
"                            frequencies for finite differences.\n"
"\n"
"    -h           Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt eigenfn -Cphase_R.nc -F20 -N0/10 -Wegn_R.nc+z0/0.06/1e-3\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    int opt;
    while ((opt = getopt(argc, argv, ":C:N::F:U:K:W:sh" )) != -1) {
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
                    Ctrl->N.nmode = (a2 - a1) / dif + 1;
                    Ctrl->N.modes = GRT_SAFE_CALLOC(Ctrl->N.nmode, sizeof(size_t));
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
                        Ctrl->F.nf = floor((a2 - a1) / df + 1e-8) + 1;
                        Ctrl->F.freqs = GRT_SAFE_CALLOC(Ctrl->F.nf, sizeof(real_t));
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
                        Ctrl->F.freqs = GRT_SAFE_CALLOC(Ctrl->F.nf, sizeof(real_t));
                        Ctrl->F.freqs[0] = a1;
                        Ctrl->F.freqs[1] = a2;
                        Ctrl->F.def_range = true;
                    }
                    
                    Ctrl->F.io_period = isperiod;

                    GRT_SAFE_FREE_PTR(string);
                }
                break;
            
            // 输出群速度 -U<path>
            case 'U':
                Ctrl->U.active = true;
                Ctrl->U.s_udisp_filepath = strdup(optarg);
                break;

            // 输出能量积分和敏感核 -K+c<path>+u<path>+z<dz>+x<path>+s<shift>
            case 'K':
                Ctrl->K.active = true;
                {
                    char *string = strdup(optarg);
                    char *token = strtok(string, "+");
                    while(token != NULL){
                        switch (token[0]){
                            case 'c':
                                Ctrl->K.s_cpar_filepath = strdup(token+1);
                                if(strlen(token)==1){
                                    GRTBadOptionError(K, "+%s wrong usage.", token);
                                }
                                break;
                            case 'u':
                                Ctrl->K.s_upar_filepath = strdup(token+1);
                                if(strlen(token)==1){
                                    GRTBadOptionError(K, "+%s wrong usage.", token);
                                }
                                break;                     
                            case 'z':
                                if(1 != sscanf(token+1, "%lf", &Ctrl->K.cpar_dz)){
                                    GRTBadOptionError(K, "+%s wrong usage.", token);
                                }
                                if(Ctrl->K.cpar_dz <= 0.0){
                                    GRTBadOptionError(K, "+%s can't set nonpositive dz(%lf).", token, Ctrl->K.cpar_dz);
                                }
                                break;
                            case 'x':
                                Ctrl->K.s_egyint_filepath = strdup(token+1);
                                if(strlen(token)==1){
                                    GRTBadOptionError(K, "+%s wrong usage.", token);
                                }
                                break;
                            case 's': {
                                char extra;
                                if(sscanf(token+1, "%lf%c", &Ctrl->K.freq_shift, &extra) != 1 || !(Ctrl->K.freq_shift > 0.0)){
                                    GRTBadOptionError(K, "+s requires a positive frequency shift (Hz).");
                                }
                                break;
                            }
                            default:
                                GRTBadOptionError(K, "+%s is not supported.", token);
                                break;
                        }
                        token = strtok(NULL, "+");
                    }
                    GRT_SAFE_FREE_PTR(string);
                }
                break;

            // 输出深度采样的本征函数  -Wpath[+zz1[/z2/dz]]
            case 'W':
                Ctrl->W.active = true;
                {
                    char *string = strdup(optarg);
                    char *token = strtok(string, "+");
                    if(strlen(token) == 0){
                        GRTBadOptionError(W, "");
                    }
                    Ctrl->W.s_egn_filepath = strdup(token);
                    
                    token = strtok(NULL, "+");
                    if(token != NULL){
                        switch (token[0]){
                            case 'z': {
                                real_t a1, a2, dif;
                                a1 = a2 = dif = 0;
                                int nscan = sscanf(token+1, "%lf/%lf/%lf", &a1, &a2, &dif);
                                if( !(nscan == 1 || nscan == 3)){
                                    GRTBadOptionError(Z, "");
                                };
                                if(nscan == 1 && a1 <= 0.0){
                                    GRTBadOptionError(Z, "Can't set a single nonpositive value.");
                                }
                                if(nscan == 1){
                                    a2 = a1;
                                    dif = 1;
                                }
                                else if(nscan == 3){
                                    if(dif <= 0){
                                        GRTBadOptionError(Z, "Can't set nonpositive dz(%lf).", dif);
                                    }
                                    if(a1 < 0.0 || a2 <= 0.0){
                                        GRTBadOptionError(Z, "Can't set nonpositive z1(%lf), z2(%lf)", a1, a2);
                                    }
                                    if(a1 > a2){
                                        GRTBadOptionError(Z, "z1(%lf) > z2(%lf).", a1, a2);
                                    }
                                }

                                Ctrl->W.nz = floor((a2 - a1) / dif + 1e-8) + 1;
                                Ctrl->W.zs = GRT_SAFE_CALLOC(Ctrl->W.nz, sizeof(real_t));
                                for(size_t i=0; i<Ctrl->W.nz; ++i){
                                    Ctrl->W.zs[i] = a1 + dif*i;
                                }
                                break;
                            }
                            default:
                                GRTBadOptionError(W, "+%s is not supported.", token);
                                break;
                        }
                    }
                    GRT_SAFE_FREE_PTR(string);
                }
                break;

            GRT_Common_Options_in_Switch((char)(optopt));

        }
    }


    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, C);
    if(Ctrl->K.freq_shift > 0.0 && Ctrl->K.s_upar_filepath == NULL){
        GRTBadOptionError(K, "+s requires +u<path>.");
    }

    // 是否需要计算能量积分
    Ctrl->calc_egyint = Ctrl->U.active || Ctrl->K.active;

    // -N 的默认项，仅处理基阶
    if( ! Ctrl->N.active ){
        Ctrl->N.nmode = 1;
        Ctrl->N.modes = GRT_SAFE_CALLOC(Ctrl->N.nmode, sizeof(size_t));
        Ctrl->N.modes[0] = 0;
    }

}



/* 子模块主函数 */
int eigenfn_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));

    // 传入参数 
    getopt_from_command(Ctrl, argc, argv);

    // 读取频散及其中保存的模型
    MODEL1D *mod1d = NULL;
    EIGENV_INFO *eigmet = GRT_SAFE_CALLOC(1, sizeof(EIGENV_INFO));
    grt_read_dispersion(Ctrl->C.s_filepath, eigmet, &mod1d);

    // 根据命令行参数确定出所需的部分频散信息
    EIGENFN_INFO *eigfnmet = GRT_SAFE_CALLOC(1, sizeof(EIGENFN_INFO));
    grt_filter_eigenfn_info(
        Ctrl->F.nf, Ctrl->F.freqs, Ctrl->F.def_range,
        Ctrl->N.nmode, Ctrl->N.modes, eigmet, eigfnmet);
    
    // 命令行只映射输出和采样选项，计算与内存管理统一交给公共接口
    EIGENFN_OPTIONS options = {
        .eigenfunctions     = Ctrl->W.active,
        .energy             = Ctrl->calc_egyint,
        .phase_kernel       = Ctrl->K.s_cpar_filepath != NULL,
        .group_kernel       = Ctrl->K.s_upar_filepath != NULL,
        .group_kernel_shift = Ctrl->K.freq_shift,
        .kernel_dz          = Ctrl->K.cpar_dz,
        .nz                 = Ctrl->W.nz,
        .zs                 = Ctrl->W.zs,
    };
    grt_compute_eigenfn(mod1d, eigfnmet, &options);

    // 打印本征函数结果
    if(Ctrl->W.active){
        int ncols = eigfnmet->wtype == GRT_DISPERSION_RAYL ? GRT_RAYL_DIM : GRT_LOVE_DIM;
        grt_output_eigenfns(Ctrl->W.s_egn_filepath, ncols, eigfnmet);
    }

    // 输出群速度
    if(Ctrl->U.active){
        grt_output_udisp(Ctrl->U.s_udisp_filepath, eigfnmet);
    }

    // 输出能量积分
    if(Ctrl->K.s_egyint_filepath!=NULL){
        grt_output_energy_integrals(Ctrl->K.s_egyint_filepath, eigfnmet);
    }

    // 输出相速度敏感核
    if(Ctrl->K.s_cpar_filepath!=NULL){
        grt_output_sensitivity(Ctrl->K.s_cpar_filepath, "c", eigfnmet);
    }
    

    // 计算并输出群速度敏感核
    if(Ctrl->K.s_upar_filepath!=NULL){
        grt_output_sensitivity(Ctrl->K.s_upar_filepath, "u", eigfnmet);
    }


    grt_free_eigenfn_info(eigfnmet);
    grt_free_eigenv_info(eigmet);
    grt_free_mod1d(mod1d);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
