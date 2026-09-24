/**
 * @file   grt_lamb1.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    求解第一类 Lamb 问题的主函数
 */

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    /** 模型参数 */
    struct {
        bool active;
        real_t nu;    ///<  泊松比
    } P;

    /** 无量纲时间序列 tbar=t/(r/beta) */
    struct {
        bool active;
        real_t *ts;
        int nt;
    } T;

    /** 方位角 */
    struct {
        bool active;
        real_t azimuth;  ///<  方位角，单位为度
    } A;

    /** 无量纲运动源速度 */
    struct {
        bool active;
        real_t cbar;  ///<  c/beta
    } C;

    /** 选择输出的震相 */
    struct {
        bool active;
        char *phase_list;
    } Q;

} GRT_MODULE_CTRL;



/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl->T.ts);
    GRT_SAFE_FREE_PTR(Ctrl->Q.phase_list);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/**
 * 打印使用说明
 */
static void print_help(){
printf("\n"
"[grt lamb1] %s\n\n", GRT_VERSION);printf(
"    Compute the exact generalized closed-form solution for the first-kind Lamb problem\n"
"    (both the source and receiver are on the surface).\n"
"\n"
"    Without -C, output contains dimensionless time and the 9 step-convolved fixed-source Green functions Gij.\n"
"    With -C, output contains dimensionless time and the 3 vertical-force displacements u1, u2, u3.\n"
"    To recover the physical quantities, you can:\n"
"       + G_{ij} <- G_{ij} / (pi^2*mu*r) for fixed sources\n"
"       + u_i <- u_i / (pi^2*mu*r) for moving sources\n"
"    where mu is the shear modulus and r is the source-receiver distance.\n"
"\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt lamb1 -P<nu> -T<t1>/<t2>/<dt> -A<azimuth> [-C<cbar>]\n"
"               [-Q<phases>]\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -P<nu>         Poisson ratio of the halfspace, (0, 0.5).\n"
"\n"
"    -T<t1>/<t2>/<dt>\n"
"                   Dimensionless time tbar = t/(r/beta) = beta*t/r.\n"
"                   Here t is physical time, r is the direct source-receiver distance,\n"
"                   and beta is the S-wave speed.\n"
"                   <t1>: start time.\n"
"                   <t2>: end time.\n"
"                   <dt>: time interval.\n"
"\n"
"    -A<azimuth>    Azimuth in degree, from source to station.\n"
"\n"
"    -C<cbar>      Enable the moving vertical point-force mode, where cbar = c/beta.\n"
"                   cbar must be positive and smaller than vR/beta; the station must be\n"
"                   off the x1 axis. Output contains only u1, u2, u3. Without -C, output\n"
"                   contains the 9 fixed-source Green functions.\n"
"\n"
"    -Q<phases>     Keep only selected phase terms in fixed-source mode. Use a\n"
"                   comma-separated list of P, S and R. It is ignored with -C.\n"
"                   If no valid phase remains, a warning is issued and the output is all zeros.\n"
"\n"
"    -h             Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt lamb1 -P0.25 -T0/2/1e-3 -A30\n"
"    grt lamb1 -P0.25 -T0/2/1e-3 -A30 -C0.1\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    int opt;

    while ((opt = getopt(argc, argv, ":P:T:A:C:Q:h")) != -1) {
        switch (opt) {
            // 模型参数， -P<nu>
            case 'P':
                Ctrl->P.active = true;
                {
                    char extra;
                    if(1 != sscanf(optarg, "%lf%c", &Ctrl->P.nu, &extra)){
                        GRTBadOptionError(P, "expected nu.");
                    }
                }
                if(Ctrl->P.nu <= 0.0 || Ctrl->P.nu >= 0.5){
                    GRTBadOptionError(P, "poisson ratio (%lf) is out of bound.", Ctrl->P.nu);
                }
                break;

            // 无量纲时间序列 tbar, -Tt1/t2/dt
            case 'T':
                Ctrl->T.active = true;
                {
                    real_t a1, a2, delta;
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &delta)){
                        GRTBadOptionError(T, "");
                    };
                    if(a1 < 0.0 || a2 < 0.0){
                        GRTBadOptionError(T, "t1 < 0.0 or t2 < 0.0.");
                    }
                    if(delta <= 0.0){
                        GRTBadOptionError(T, "dt <= 0.0.");
                    }
                    if(a1 > a2){
                        GRTBadOptionError(T, "t1(%f) > t2(%f).", a1, a2);
                    }

                    Ctrl->T.nt = floor((a2-a1)/delta) + 1;
                    Ctrl->T.ts = GRT_SAFE_CALLOC(Ctrl->T.nt, sizeof(real_t));
                    for(int i=0; i<Ctrl->T.nt; ++i){
                        Ctrl->T.ts[i] = a1 + delta*i;
                    }
                }
                break;

            // 方位角，  -Aazimuth
            case 'A':
                Ctrl->A.active = true;
                if(1 != sscanf(optarg, "%lf", &Ctrl->A.azimuth)){
                    GRTBadOptionError(A, "");
                }
                if(Ctrl->A.azimuth < 0.0 || Ctrl->A.azimuth > 360){
                    GRTBadOptionError(A, "azimuth should be in [0, 360].");
                }
                break;

            // 无量纲运动源速度，-Ccbar
            case 'C':
                Ctrl->C.active = true;
                {
                    char extra;
                    if(1 != sscanf(optarg, "%lf%c", &Ctrl->C.cbar, &extra)){
                        GRTBadOptionError(C, "expected a positive cbar value.");
                    }
                }
                if(!isfinite(Ctrl->C.cbar) || Ctrl->C.cbar <= 0.0){
                    GRTBadOptionError(C, "cbar should be finite and positive.");
                }
                break;

            case 'Q':
                GRT_SAFE_FREE_PTR(Ctrl->Q.phase_list);
                Ctrl->Q.phase_list = strdup(optarg);
                Ctrl->Q.active = true;
                break;
            
            GRT_Common_Options_in_Switch((char)(optopt)); 
        }
    }

    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, P);
    GRTCheckOptionActive(Ctrl, T);
    GRTCheckOptionActive(Ctrl, A);
}



/** 模块主函数 */
int lamb1_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));

    // 传入参数 
    getopt_from_command(Ctrl, argc, argv);

    // 运动源模式只输出竖向力源对应的三个位移分量
    if(Ctrl->C.active){
        real_t (*u)[3][3] = GRT_SAFE_CALLOC(Ctrl->T.nt, sizeof(*u));
        grt_solve_lamb1(Ctrl->P.nu, Ctrl->T.ts, Ctrl->T.nt, Ctrl->A.azimuth,
            Ctrl->C.cbar, Ctrl->Q.active ? Ctrl->Q.phase_list : NULL, u);

        printf("#%13s%14s%14s%14s\n", "tbar", "u1", "u2", "u3");
        for(int i=0; i<Ctrl->T.nt; ++i){
            printf("%14.6e", Ctrl->T.ts[i]);
            for(int j=0; j<3; ++j){
                printf("%14.6e", u[i][j][2]);
            }
            printf("\n");
        }
        GRT_SAFE_FREE_PTR(u);
    }
    else{
        grt_solve_lamb1(Ctrl->P.nu, Ctrl->T.ts, Ctrl->T.nt, Ctrl->A.azimuth,
            0.0, Ctrl->Q.active ? Ctrl->Q.phase_list : NULL, NULL);
    }

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
