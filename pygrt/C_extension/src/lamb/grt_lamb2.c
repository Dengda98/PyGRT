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
        real_t depsrc;
        real_t deprcv;
    } D;

    /** 空间导数输出路径 */
    struct {
        bool active;
        char *source_path;
        char *receiver_path;
        char *mixed_path;
    } S;

    /** 方位角 */
    struct {
        bool active;
        real_t azimuth;
    } A;

    /** 选择输出的震相 */
    struct {
        bool active;
        char *phase_list;
    } Q;

} GRT_MODULE_CTRL;


static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    GRT_SAFE_FREE_PTR(Ctrl->T.ts);
    GRT_SAFE_FREE_PTR(Ctrl->S.source_path);
    GRT_SAFE_FREE_PTR(Ctrl->S.receiver_path);
    GRT_SAFE_FREE_PTR(Ctrl->S.mixed_path);
    GRT_SAFE_FREE_PTR(Ctrl->Q.phase_list);
    GRT_SAFE_FREE_PTR(Ctrl);
}


static void print_help(void)
{
printf("\n"
"[grt lamb2] %s\n\n", GRT_VERSION);printf(
"    Compute the exact generalized closed-form solution for the second-kind Lamb problem\n"
"    (exactly one of the source and receiver is on the free surface).\n"
"\n"
"    All outputs are dimensionless and convolved with the step function.\n"
"    Standard output contains dimensionless time and the 9 displacement Green functions Gij.\n"
"    Each requested derivative file contains dimensionless time and the corresponding\n"
"    27 first derivatives or 81 mixed second derivatives.\n"
"    To recover the physical quantities, you can:\n"
"       + G_{ij} <- G_{ij} / (pi^2*mu*r)\n"
"       + G_{ij,k'} <- G_{ij,k'} / (pi^2*mu*r^2)\n"
"       + G_{ij,k}  <- G_{ij,k} / (pi^2*mu*r^2)\n"
"       + G_{ij,k,k'}  <- G_{ij,k,k'} / (pi^2*mu*r^3)\n"
"    where mu is the shear modulus and r is the source-receiver distance.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt lamb2 -P<nu> -T<t1>/<t2>/<dt> -R<dist>\n"
"              (-Ds<depsrc> | -Dr<deprcv>) -A<azimuth>\n"
"              [-Q<phases>]\n"
"              [-S[+s<source-path>][+r<receiver-path>][+m<mixed-path>]]\n"
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
"                   Mutually exclusive with -Dr. Values with depsrc/r <= 1e-3 trigger a\n"
"                   numerical-stability warning.\n"
"\n"
"    -Dr<deprcv>    Receiver depth, strictly positive. The source is then on the surface.\n"
"                   Mutually exclusive with -Ds. Values with deprcv/r <= 1e-3 trigger a\n"
"                   numerical-stability warning. This case is obtained from the buried-source\n"
"                   solution by reciprocity.\n"
"\n"
"    -S[+s<source-path>][+r<receiver-path>][+m<mixed-path>]\n"
"                   +s<source-path>: save source-coordinate derivatives G_{ij,k'}\n"
"                   +r<receiver-path>: save receiver-coordinate derivatives G_{ij,k}\n"
"                   +m<mixed-path>: save mixed second derivatives G_{ij,k,k'}\n"
"                   At least one suboption is required; multiple suboptions may be combined.\n"
"\n"
"    -A<azimuth>    Azimuth in degree, from source to receiver, [0, 360].\n"
"\n"
"    -Q<phases>     Keep only selected phase terms. Use a comma-separated list\n"
"                   of P, S, SP and PS. Only SP is\n"
"                   available for an underground source; only PS is available\n"
"                   for an underground receiver.\n"
"                   If no valid phase remains, a warning is issued and the output is all zeros.\n"
"\n"
"    -h             Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt lamb2 -P0.25 -T0/3/1e-3 -R10 -Ds5 -A30\n"
"    grt lamb2 -P0.25 -T0/3/1e-3 -R10 -Dr5 -A30\n"
"    grt lamb2 -P0.25 -T0/3/1e-3 -R10 -Ds5 -A30 -S+slamb2_source.txt+rlamb2_receiver.txt+mlamb2_mixed.txt\n"
"\n\n\n");
}


static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    int opt;
    while ((opt = getopt(argc, argv, ":P:T:R:D:S:A:Q:h")) != -1) {
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
                    Ctrl->T.ts = GRT_SAFE_CALLOC((size_t)Ctrl->T.nt, sizeof(*Ctrl->T.ts));
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
                        if (1 != sscanf(optarg + 1, "%lf%c", &Ctrl->D.depsrc, &extra) || Ctrl->D.depsrc <= 0.0) {
                            GRTBadOptionError(Ds, "source depth should be strictly positive.");
                        }
                    }
                } else if (optarg[0] == 'r') {
                    Ctrl->D.r_active = true;
                    {
                        char extra;
                        if (1 != sscanf(optarg + 1, "%lf%c", &Ctrl->D.deprcv, &extra) || Ctrl->D.deprcv <= 0.0) {
                            GRTBadOptionError(Dr, "receiver depth should be strictly positive.");
                        }
                    }
                } else {
                    GRTBadOptionError(D, "use -Ds<depsrc> or -Dr<deprcv>.");
                }
                break;

            case 'S':
                Ctrl->S.active = true;
                grt_lamb_parse_derivative_paths_with_mixed(
                    optarg, &Ctrl->S.source_path, &Ctrl->S.receiver_path, &Ctrl->S.mixed_path);
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

            case 'Q':
                GRT_SAFE_FREE_PTR(Ctrl->Q.phase_list);
                Ctrl->Q.phase_list = strdup(optarg);
                Ctrl->Q.active = true;
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
    real_t (*G)[3][3] = GRT_SAFE_CALLOC(nt, sizeof(*G));
    real_t (*dG_source)[3][3][3] = GRT_SAFE_CALLOC(nt, sizeof(*dG_source));
    real_t (*dG_receiver)[3][3][3] = GRT_SAFE_CALLOC(nt, sizeof(*dG_receiver));
    real_t (*dG_mixed)[3][3][3][3] = Ctrl->S.mixed_path != NULL ? GRT_SAFE_CALLOC(nt, sizeof(*dG_mixed)) : NULL;

    FILE *source_file = NULL;
    FILE *receiver_file = NULL;
    if (Ctrl->S.source_path != NULL) {
        source_file = GRTCheckOpenFile(Ctrl->S.source_path, "w");
    }
    if (Ctrl->S.receiver_path != NULL) {
        receiver_file = GRTCheckOpenFile(Ctrl->S.receiver_path, "w");
    }
    FILE *mixed_file = NULL;
    if (Ctrl->S.mixed_path != NULL) {
        mixed_file = GRTCheckOpenFile(Ctrl->S.mixed_path, "w");
    }

    grt_solve_lamb2(Ctrl->P.nu, Ctrl->T.ts, Ctrl->T.nt, Ctrl->R.distance,
        Ctrl->D.depsrc, Ctrl->D.deprcv, Ctrl->A.azimuth,
        Ctrl->Q.active ? Ctrl->Q.phase_list : NULL,
        G, dG_source, dG_receiver, dG_mixed);
    grt_lamb_print_green_series(stdout, Ctrl->T.ts, Ctrl->T.nt, G);
    if (source_file != NULL) {
        grt_lamb_print_derivative_series(source_file, Ctrl->T.ts, Ctrl->T.nt, dG_source, true);
        fclose(source_file);
    }
    if (receiver_file != NULL) {
        grt_lamb_print_derivative_series(receiver_file, Ctrl->T.ts, Ctrl->T.nt, dG_receiver, false);
        fclose(receiver_file);
    }
    if (mixed_file != NULL) {
        grt_lamb_print_mixed_derivative_series(mixed_file, Ctrl->T.ts, Ctrl->T.nt, dG_mixed);
        fclose(mixed_file);
    }

    GRT_SAFE_FREE_PTR(G);
    GRT_SAFE_FREE_PTR(dG_source);
    GRT_SAFE_FREE_PTR(dG_receiver);
    GRT_SAFE_FREE_PTR(dG_mixed);
}


int lamb2_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);
    if (Ctrl->S.active) {
        run_lamb2_with_derivative_outputs(Ctrl);
    } else {
        grt_solve_lamb2(Ctrl->P.nu, Ctrl->T.ts, Ctrl->T.nt, Ctrl->R.distance,
            Ctrl->D.depsrc, Ctrl->D.deprcv, Ctrl->A.azimuth,
            Ctrl->Q.active ? Ctrl->Q.phase_list : NULL, NULL, NULL, NULL, NULL);
    }
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
