/**
 * @file   grt_lamb3.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-09
 *
 *    求解第三类 Lamb 问题的主函数
 */

#include "grt.h"


typedef struct {
    struct {
        bool active;
        real_t nu;
    } P;

    /** 无量纲时间序列 tbar=t/(r/beta) */
    struct {
        bool active;
        real_t *ts;
        int nt;
    } T;

    struct {
        bool active;
        real_t distance;
    } R;

    struct {
        bool active;
        real_t source_depth;
        real_t receiver_depth;
    } D;

    struct {
        bool active;
        char *source_path;
        char *receiver_path;
    } S;

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
"[grt lamb3] %s\n\n", GRT_VERSION);printf(
"    Compute the generalized closed-form solution for the third-kind Lamb problem\n"
"    (both the source and receiver are underground).\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt lamb3 -P<nu> -T<t1>/<t2>/<dt> -R<dist>\n"
"               -D<depsrc>/<deprcv> -A<azimuth>\n"
"               [-S+s<source-path>+r<receiver-path>]\n"
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
"    -R<dist>       Horizontal source-receiver distance, positive.\n"
"                   Values with R/r' <= 1e-2 trigger a numerical-stability warning,\n"
"                   where r'=sqrt(R^2+(depsrc+deprcv)^2).\n"
"\n"
"    -D<depsrc>/<deprcv>\n"
"                   Source and receiver depths, both strictly positive.\n"
"                   Values with depsrc/r or deprcv/r < 1e-3 trigger a numerical-stability warning.\n"
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
"    Here r is the straight source-receiver distance sqrt(R^2+(depsrc-deprcv)^2).\n"
"    Divide Gij by pi^2*mu*r and derivatives by pi^2*mu*r^2 to recover physical quantities.\n"
"\n\n"
"Example:\n"
"----------------------------------------------------------------\n"
"    grt lamb3 -P0.25 -T0/2/1e-3 -R10 -D2/1 -A30\n"
"    grt lamb3 -P0.25 -T0/2/1e-3 -R10 -D2/1 -A30 -S+slamb3_source.txt+rlamb3_receiver.txt\n"
"\n\n\n"
);
}


static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    int opt;
    while((opt = getopt(argc, argv, ":P:T:R:D:S:A:h")) != -1){
        switch(opt){
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

            case 'T':
                Ctrl->T.active = true;
                {
                    real_t t1, t2, dt;
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &t1, &t2, &dt)){
                        GRTBadOptionError(T, "");
                    }
                    if(t1 < 0.0 || t2 < 0.0){
                        GRTBadOptionError(T, "t1 < 0.0 or t2 < 0.0.");
                    }
                    if(dt <= 0.0){
                        GRTBadOptionError(T, "dt <= 0.0.");
                    }
                    if(t1 > t2){
                        GRTBadOptionError(T, "t1(%f) > t2(%f).", t1, t2);
                    }
                    Ctrl->T.nt = (int)floor((t2 - t1) / dt) + 1;
                    Ctrl->T.ts = calloc((size_t)Ctrl->T.nt, sizeof(*Ctrl->T.ts));
                    for(int i = 0; i < Ctrl->T.nt; ++i){
                        Ctrl->T.ts[i] = t1 + dt * i;
                    }
                }
                break;

            case 'R':
                Ctrl->R.active = true;
                {
                    char extra;
                    if(1 != sscanf(optarg, "%lf%c", &Ctrl->R.distance, &extra) || Ctrl->R.distance <= 0.0){
                        GRTBadOptionError(R, "horizontal distance should be positive.");
                    }
                }
                break;

            case 'S':
                Ctrl->S.active = true;
                grt_lamb_parse_derivative_paths(optarg, &Ctrl->S.source_path, &Ctrl->S.receiver_path);
                if(Ctrl->S.source_path != NULL && Ctrl->S.receiver_path != NULL &&
                    strcmp(Ctrl->S.source_path, Ctrl->S.receiver_path) == 0){
                    GRTBadOptionError(S, "source and receiver derivative paths must be different.");
                }
                break;

            case 'D':
                Ctrl->D.active = true;
                {
                    char extra;
                    if(2 != sscanf(optarg, "%lf/%lf%c", &Ctrl->D.source_depth, &Ctrl->D.receiver_depth, &extra)){
                        GRTBadOptionError(D, "expected source-depth/receiver-depth.");
                    }
                }
                if(Ctrl->D.source_depth <= 0.0 || Ctrl->D.receiver_depth <= 0.0){
                    GRTBadOptionError(D, "source and receiver depths should be strictly positive.");
                }
                break;

            case 'A':
                Ctrl->A.active = true;
                {
                    char extra;
                    if(1 != sscanf(optarg, "%lf%c", &Ctrl->A.azimuth, &extra)){
                        GRTBadOptionError(A, "");
                    }
                    if(Ctrl->A.azimuth < 0.0 || Ctrl->A.azimuth > 360.0){
                        GRTBadOptionError(A, "azimuth is out of bound.");
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
    GRTCheckOptionActive(Ctrl, D);
    GRTCheckOptionActive(Ctrl, A);
}


static void run_lamb3_with_derivative_outputs(const GRT_MODULE_CTRL *Ctrl)
{
    const size_t nt = (size_t)Ctrl->T.nt;
    real_t (*G)[3][3] = calloc(nt, sizeof(*G));
    real_t (*dG_source)[3][3][3] = calloc(nt, sizeof(*dG_source));
    real_t (*dG_receiver)[3][3][3] = calloc(nt, sizeof(*dG_receiver));
    if (G == NULL || dG_source == NULL || dG_receiver == NULL) {
        GRT_SAFE_FREE_PTR(G);
        GRT_SAFE_FREE_PTR(dG_source);
        GRT_SAFE_FREE_PTR(dG_receiver);
        GRTRaiseError("Cannot allocate lamb3 output arrays.\n");
    }

    FILE *source_file = NULL;
    FILE *receiver_file = NULL;
    if (Ctrl->S.source_path != NULL) {
        source_file = GRTCheckOpenFile(Ctrl->S.source_path, "w");
    }
    if (Ctrl->S.receiver_path != NULL) {
        receiver_file = GRTCheckOpenFile(Ctrl->S.receiver_path, "w");
    }

    grt_solve_lamb3(Ctrl->P.nu, Ctrl->T.ts, Ctrl->T.nt, Ctrl->R.distance, Ctrl->D.source_depth,
        Ctrl->D.receiver_depth, Ctrl->A.azimuth, G, dG_source, dG_receiver);
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


int lamb3_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);
    if (Ctrl->S.active) {
        run_lamb3_with_derivative_outputs(Ctrl);
    } else {
        grt_solve_lamb3(Ctrl->P.nu, Ctrl->T.ts, Ctrl->T.nt, Ctrl->R.distance, Ctrl->D.source_depth,
            Ctrl->D.receiver_depth, Ctrl->A.azimuth, NULL, NULL, NULL);
    }
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
