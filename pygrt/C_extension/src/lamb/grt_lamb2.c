/**
 * @file   grt_lamb2.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 *    求解第二类 Lamb 问题的主函数
 */

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    /** 模型参数 */
    struct {
        bool active;
        real_t nu; ///< 泊松比
    } P;

    /** 无量纲时间序列 tbar=t/(r/beta) */
    struct {
        bool active;
        real_t *ts;
        int nt;
    } T;

    /** 水平震中距 */
    struct {
        bool active;
        real_t distance;
    } R;

    /** 源点或接收点深度：-Ds 与 -Dr 互斥，未设置的一侧为 0 */
    struct {
        bool s_active;       ///< -Ds
        bool r_active;       ///< -Dr
        real_t source_depth;
        real_t receiver_depth;
    } D;

    /** 空间导数输出路径 */
    struct {
        bool active;
        char *source_path;
        char *receiver_path;
    } S;

    /** 方位角 */
    struct {
        bool active;
        real_t azimuth;
    } A;

} GRT_MODULE_CTRL;


static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    GRT_SAFE_FREE_PTR(Ctrl->T.ts);
    GRT_SAFE_FREE_PTR(Ctrl->S.source_path);
    GRT_SAFE_FREE_PTR(Ctrl->S.receiver_path);
    GRT_SAFE_FREE_PTR(Ctrl);
}


static void print_help(void)
{
printf("\n"
"[grt lamb2] %s\n\n", GRT_VERSION);printf(
"    Compute the exact generalized closed-form solution for the second-kind Lamb problem\n"
"    (exactly one of the source and receiver is on the free surface).\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt lamb2 -P<nu> -T<t1>/<t2>/<dt> -R<dist>\n"
"              (-Ds<depsrc> | -Dr<deprcv>) -A<azimuth>\n"
"              [-S+s<source-path>+r<receiver-path>]\n"
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
"    -R<dist>       Horizontal epicentral distance from source to receiver, positive.\n"
"                   Values with R/r <= 1e-3 trigger a numerical-stability warning.\n"
"\n"
"    -Ds<depsrc>    Source depth, strictly positive. The receiver is then on the surface.\n"
"                   Mutually exclusive with -Dr. Values with depsrc/r < 1e-3 trigger a\n"
"                   numerical-stability warning.\n"
"\n"
"    -Dr<deprcv>    Receiver depth, strictly positive. The source is then on the surface.\n"
"                   Mutually exclusive with -Ds. Values with deprcv/r < 1e-3 trigger a\n"
"                   numerical-stability warning. This case is obtained from the buried-source\n"
"                   solution by reciprocity.\n"
"\n"
"    -S+s<path>+r<path>\n"
"                   Save source-coordinate derivatives to +s<path> and receiver-coordinate\n"
"                   derivatives to +r<path>. Either suboption may be omitted.\n"
"\n"
"    -A<azimuth>    Azimuth in degree, from source to receiver, [0, 360].\n"
"\n"
"    -h             Display this help message.\n"
"\n\n"
"Output:\n"
"----------------------------------------------------------------\n"
"    Standard output contains dimensionless time and the 9 dimensionless step-force\n"
"    displacement Green functions Gij. When -S is used, each specified derivative file\n"
"    contains dimensionless time and the corresponding 27 dimensionless derivatives.\n"
"    Here r=sqrt(R^2+h^2) is the straight source-receiver distance, \n"
"    and h is the underground depth. Divide Gij by\n"
"    pi^2*mu*r and derivatives by pi^2*mu*r^2 to recover physical quantities.\n"
"\n\n"
"Example:\n"
"----------------------------------------------------------------\n"
"    grt lamb2 -P0.25 -T0/3/1e-3 -R10 -Ds5 -A30\n"
"    grt lamb2 -P0.25 -T0/3/1e-3 -R10 -Dr5 -A30\n"
"    grt lamb2 -P0.25 -T0/3/1e-3 -R10 -Ds5 -A30 -S+slamb2_source.txt+rlamb2_receiver.txt\n"
"\n\n\n");
}


static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    int opt;
    while ((opt = getopt(argc, argv, ":P:T:R:D:S:A:h")) != -1) {
        switch (opt) {
            case 'P':
                Ctrl->P.active = true;
                {
                    char extra;
                    if (1 != sscanf(optarg, "%lf%c", &Ctrl->P.nu, &extra)) {
                        GRTBadOptionError(P, "expected nu.");
                    }
                }
                if (Ctrl->P.nu <= 0.0 || Ctrl->P.nu >= 0.5) {
                    GRTBadOptionError(P, "poisson ratio (%lf) is out of bound.", Ctrl->P.nu);
                }
                break;

            case 'T':
                Ctrl->T.active = true;
                {
                    real_t t1, t2, dt;
                    if (3 != sscanf(optarg, "%lf/%lf/%lf", &t1, &t2, &dt)) {
                        GRTBadOptionError(T, "");
                    }
                    if (t1 < 0.0 || t2 < 0.0) {
                        GRTBadOptionError(T, "t1 < 0.0 or t2 < 0.0.");
                    }
                    if (dt <= 0.0) {
                        GRTBadOptionError(T, "dt <= 0.0.");
                    }
                    if (t1 > t2) {
                        GRTBadOptionError(T, "t1(%f) > t2(%f).", t1, t2);
                    }
                    Ctrl->T.nt = (int)floor((t2 - t1) / dt) + 1;
                    Ctrl->T.ts = (real_t *)calloc((size_t)Ctrl->T.nt, sizeof(*Ctrl->T.ts));
                    for (int i = 0; i < Ctrl->T.nt; ++i) {
                        Ctrl->T.ts[i] = t1 + dt * i;
                    }
                }
                break;

            case 'R':
                Ctrl->R.active = true;
                {
                    char extra;
                    if (1 != sscanf(optarg, "%lf%c", &Ctrl->R.distance, &extra) || Ctrl->R.distance <= 0.0) {
                        GRTBadOptionError(R, "horizontal distance should be positive.");
                    }
                }
                break;

            case 'D':
                if (optarg[0] == 's') {
                    Ctrl->D.s_active = true;
                    {
                        char extra;
                        if (1 != sscanf(optarg + 1, "%lf%c", &Ctrl->D.source_depth, &extra) || Ctrl->D.source_depth <= 0.0) {
                            GRTBadOptionError(Ds, "source depth should be strictly positive.");
                        }
                    }
                } else if (optarg[0] == 'r') {
                    Ctrl->D.r_active = true;
                    {
                        char extra;
                        if (1 != sscanf(optarg + 1, "%lf%c", &Ctrl->D.receiver_depth, &extra) || Ctrl->D.receiver_depth <= 0.0) {
                            GRTBadOptionError(Dr, "receiver depth should be strictly positive.");
                        }
                    }
                } else {
                    GRTBadOptionError(D, "use -Ds<depsrc> or -Dr<deprcv>.");
                }
                break;

            case 'S':
                Ctrl->S.active = true;
                grt_lamb_parse_derivative_paths(optarg, &Ctrl->S.source_path, &Ctrl->S.receiver_path);
                if (Ctrl->S.source_path != NULL && Ctrl->S.receiver_path != NULL &&
                    strcmp(Ctrl->S.source_path, Ctrl->S.receiver_path) == 0) {
                    GRTBadOptionError(S, "source and receiver derivative paths must be different.");
                }
                break;

            case 'A':
                Ctrl->A.active = true;
                {
                    char extra;
                    if (1 != sscanf(optarg, "%lf%c", &Ctrl->A.azimuth, &extra)) {
                        GRTBadOptionError(A, "");
                    }
                    if (Ctrl->A.azimuth < 0.0 || Ctrl->A.azimuth > 360.0) {
                        GRTBadOptionError(A, "azimuth should be in [0, 360].");
                    }
                }
                break;

            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, P);
    GRTCheckOptionActive(Ctrl, T);
    GRTCheckOptionActive(Ctrl, R);
    if (!Ctrl->D.s_active && !Ctrl->D.r_active) {
        GRTRaiseError("Need set one of options \"-Ds\" and \"-Dr\". Use \"-h\" for help.\n");
    }
    if (Ctrl->D.s_active && Ctrl->D.r_active) {
        GRTRaiseError("Options -Ds and -Dr are mutually exclusive in lamb2.\n");
    }
    GRTCheckOptionActive(Ctrl, A);
}


static void run_lamb2_with_derivative_outputs(const GRT_MODULE_CTRL *Ctrl)
{
    const size_t nt = (size_t)Ctrl->T.nt;
    real_t (*G)[3][3] = calloc(nt, sizeof(*G));
    real_t (*dG_source)[3][3][3] = calloc(nt, sizeof(*dG_source));
    real_t (*dG_receiver)[3][3][3] = calloc(nt, sizeof(*dG_receiver));
    if (G == NULL || dG_source == NULL || dG_receiver == NULL) {
        GRT_SAFE_FREE_PTR(G);
        GRT_SAFE_FREE_PTR(dG_source);
        GRT_SAFE_FREE_PTR(dG_receiver);
        GRTRaiseError("Cannot allocate lamb2 output arrays.\n");
    }

    FILE *source_file = NULL;
    FILE *receiver_file = NULL;
    if (Ctrl->S.source_path != NULL) {
        source_file = GRTCheckOpenFile(Ctrl->S.source_path, "w");
    }
    if (Ctrl->S.receiver_path != NULL) {
        receiver_file = GRTCheckOpenFile(Ctrl->S.receiver_path, "w");
    }

    grt_solve_lamb2(Ctrl->P.nu, Ctrl->T.ts, Ctrl->T.nt, Ctrl->R.distance,
        Ctrl->D.source_depth, Ctrl->D.receiver_depth, Ctrl->A.azimuth,
        G, dG_source, dG_receiver);
    grt_lamb_print_green_series(stdout, Ctrl->T.ts, Ctrl->T.nt, G);
    if (source_file != NULL) {
        grt_lamb_print_derivative_series(source_file, Ctrl->T.ts, Ctrl->T.nt, dG_source, true);
        fclose(source_file);
    }
    if (receiver_file != NULL) {
        grt_lamb_print_derivative_series(receiver_file, Ctrl->T.ts, Ctrl->T.nt, dG_receiver, false);
        fclose(receiver_file);
    }

    GRT_SAFE_FREE_PTR(G);
    GRT_SAFE_FREE_PTR(dG_source);
    GRT_SAFE_FREE_PTR(dG_receiver);
}


int lamb2_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);
    if (Ctrl->S.active) {
        run_lamb2_with_derivative_outputs(Ctrl);
    } else {
        grt_solve_lamb2(Ctrl->P.nu, Ctrl->T.ts, Ctrl->T.nt, Ctrl->R.distance,
            Ctrl->D.source_depth, Ctrl->D.receiver_depth, Ctrl->A.azimuth, NULL, NULL, NULL);
    }
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
