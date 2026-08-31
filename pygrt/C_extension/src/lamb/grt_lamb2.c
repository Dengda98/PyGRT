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
        real_t nu;  ///< 泊松比
    } P;

    /** 归一化时间序列 */
    struct {
        bool active;
        real_t *ts;
        int nt;
    } T;

    /** 震源射线角 */
    struct {
        bool active;
        real_t theta;  ///< 震源射线角，单位度
    } D;

    /** 方位角 */
    struct {
        bool active;
        real_t azimuth;  ///< 方位角，单位度
    } A;

} GRT_MODULE_CTRL;


static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    GRT_SAFE_FREE_PTR(Ctrl->T.ts);
    GRT_SAFE_FREE_PTR(Ctrl);
}

static void print_help(void)
{
printf("\n"
"[grt lamb2] %s\n\n", GRT_VERSION);printf(
"    Compute the exact generalized closed-form solution for the second-kind Lamb problem\n"
"    (the receiver is on the surface and the source is underground).\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt lamb2 -P<nu> -T<t1>/<t2>/<dt> -D<theta> -A<azimuth>\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -P<nu>         Poisson ratio of the halfspace, (0, 0.5).\n"
"\n"
"    -T<t1>/<t2>/<dt>\n"
"                   Dimensionless time.\n"
"                   <t1>: start time.\n"
"                   <t2>: end time.\n"
"                   <dt>: time interval.\n"
"\n"
"    -D<theta>      Source ray angle in degree, (0, 90).\n"
"                   theta is measured from the upward vertical.\n"
"\n"
"    -A<azimuth>    Azimuth in degree, from source to station.\n"
"\n"
"    -h             Display this help message.\n"
"\n\n"
"Output:\n"
"----------------------------------------------------------------\n"
"    The output is dimensionless step-force displacement Green functions.\n"
"    The first 9 columns are Gij. They are followed by 27 source-coordinate\n"
"    derivatives Gij,k' and 27 receiver-coordinate derivatives Gij,k, all in\n"
"    row-major order. The physical Green function is G^H = Gbar^H/(pi^2*mu*r);\n"
"    both physical spatial derivatives are obtained by dividing their normalized\n"
"    results by pi^2*mu*r^2, where mu is the shear modulus.\n"
"\n\n"
"Example:\n"
"----------------------------------------------------------------\n"
"    grt lamb2 -P0.25 -T0/3/1e-3 -D60 -A30\n"
"\n\n\n"
);
}

static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    int opt;
    while((opt = getopt(argc, argv, ":P:T:D:A:h")) != -1){
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
                    Ctrl->T.nt = (int)floor((t2 - t1)/dt) + 1;
                    Ctrl->T.ts = (real_t*)calloc((size_t)Ctrl->T.nt, sizeof(real_t));
                    for(int i=0; i<Ctrl->T.nt; ++i){
                        Ctrl->T.ts[i] = t1 + dt*i;
                    }
                }
                break;

            case 'D':
                Ctrl->D.active = true;
                if(1 != sscanf(optarg, "%lf", &Ctrl->D.theta)){
                    GRTBadOptionError(D, "");
                }
                if(Ctrl->D.theta <= 0.0 || Ctrl->D.theta >= 90.0){
                    GRTBadOptionError(D, "theta should be in (0, 90).\n");
                }
                break;

            case 'A':
                Ctrl->A.active = true;
                if(1 != sscanf(optarg, "%lf", &Ctrl->A.azimuth)){
                    GRTBadOptionError(A, "");
                }
                if(Ctrl->A.azimuth < 0.0 || Ctrl->A.azimuth > 360.0){
                    GRTBadOptionError(A, "azimuth should be in [0, 360].");
                }
                break;

            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, P);
    GRTCheckOptionActive(Ctrl, T);
    GRTCheckOptionActive(Ctrl, D);
    GRTCheckOptionActive(Ctrl, A);
}

int lamb2_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);
    grt_solve_lamb2(Ctrl->P.nu, Ctrl->T.ts, Ctrl->T.nt, Ctrl->D.theta, Ctrl->A.azimuth, NULL, NULL, NULL);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
