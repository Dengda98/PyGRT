/**
 * @file   grt_lamb.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-09
 *
 *    根据 Lamb 闭合解计算半空间中的物理 Green 函数
 */

#include "grt.h"

/* 防止复数虚数单位宏替换积分次数成员名 */
#undef I


/** 该子模块的参数控制结构体 */
typedef struct {
    /** 半空间参数 */
    struct {
        bool active;
        real_t vp;
        real_t vs;
        real_t rho;
        real_t nu;
    } H;

    /** 输出目录 */
    struct {
        bool active;
        char *s_output_dir;
    } O;

    /** 方位角 */
    struct {
        bool active;
        real_t azimuth;
        real_t azrad;
        real_t backazimuth;
    } A;

    /** 源强缩放 */
    struct {
        bool active;
        bool mult_src_mu;
        real_t M0;
        real_t src_mu;
    } S;

    /** 剪切源 */
    struct {
        bool active;
    } M;

    /** 单力源 */
    struct {
        bool active;
    } F;

    /** 矩张量源 */
    struct {
        bool active;
    } T;

    /** 时间积分次数 */
    struct {
        bool active;
        int int_times;
    } I;

    /** 时间微分次数 */
    struct {
        bool active;
        int dif_times;
    } J;

    /** 时间函数 */
    struct {
        bool active;
        char tftype;
        char *tfparams;
    } D;

    /** 时间延迟 */
    struct {
        bool active;       ///< 是否设置时间延迟
        real_t delayT0;    ///< 参考时间偏移，s
        real_t delayV0;    ///< 参考速度，km/s
        bool refFirstP;    ///< 是否参考初至 P
    } E;

    /** 物理时间序列 */
    struct {
        bool active;
        int nt;
        real_t dt;
        real_t *tbar;
    } N;

    /** 源点和接收点深度 */
    struct {
        bool s_active;
        bool r_active;
        real_t depsrc;
        real_t deprcv;
    } Depth;

    /** 水平震中距 */
    struct {
        bool active;
        real_t dist;
    } R;

    /** 旋转到 Z、N、E */
    struct {
        bool active;
    } n;

    /** 是否静默输出 */
    struct {
        bool active;
    } s;

    /** 是否计算空间导数 */
    struct {
        bool active;
    } e;

    /** 震源机制相关参数 */
    real_t mchn[GRT_MECHANISM_NUM];

    /** 最终要计算的震源类型 */
    GRT_SYN_TYPE computeType;
    char s_computeType[3];
} GRT_MODULE_CTRL;


/** Lamb 求解器输出的各类数组 */
typedef struct {
    real_t (*G)[3][3];                    ///< 位移 Green 函数
    real_t (*dG_source)[3][3][3];         ///< 源点一阶导数
    real_t (*dG_receiver)[3][3][3];       ///< 接收点一阶导数
    real_t (*dG_mixed)[3][3][3][3];       ///< 混合二阶导数
} LAMB_RESULT;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    GRT_SAFE_FREE_PTR(Ctrl->N.tbar);
    GRT_SAFE_FREE_PTR(Ctrl->D.tfparams);
    GRT_SAFE_FREE_PTR(Ctrl->O.s_output_dir);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印 Lamb 模块的命令行使用说明 */
static void print_help(void)
{
printf("\n"
"[grt lamb] %s\n\n", GRT_VERSION);printf(
"    Compute physical synthetic seismograms from Lamb's closed-form solutions.\n"
"    The selected source components are combined into physical displacement records.\n"
"    Output filenames follow module `syn`: Z, R, T or Z, N, E.  Prefixes z, r and t\n"
"    (or z, n and e) contain spatial derivatives when -e is used.\n"
"\n"
"    The Lamb solvers use the global ZNE coordinate system with Z positive downward.\n"
"    This module converts it to Z-up/R/T by default; -n keeps the horizontal N/E\n"
"    components in the output.\n"
"\n"
"    The analytic Lamb solutions are normalized by pi^2*mu*r.  The module applies\n"
"    the corresponding physical normalization and source scale before writing SAC.\n"
"    The absolute halfspace parameters are specified by -H<vp>/<vs>/<rho>; the\n"
"    Poisson ratio required by the Lamb solvers is derived from vp/vs.\n"
"    The -D, -I and -J options are applied in the time domain.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt lamb -H<vp>/<vs>/<rho> -N<nt>/<dt> -R<dist>\n"
"              -Ds<depsrc> -Dr<deprcv>\n"
"              -A<azimuth> -S[u]<scale> -O<outdir>\n"
"              [-M<strike>/<dip>[/<rake>]]\n"
"              [-T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>] [-F<fn>/<fe>/<fz>]\n"
"              [-D<tftype>/<tfparams>] [-E[p]<t0>[/<v0>]] [-I<odr>] [-J<odr>]\n"
"              [-n] [-e] [-s]\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -H<vp>/<vs>/<rho>\n"
"                  Homogeneous halfspace parameters. vp and vs are P- and S-wave\n"
"                  speeds in km/s, and rho is density in g/cm^3. The Poisson ratio\n"
"                  is calculated from vp/vs, which must be greater than sqrt(2).\n"
"\n"
"    -N<nt>/<dt>   Number of samples and physical time interval in seconds.\n"
"                  The time series starts at the origin time.\n"
"\n"
"    -R<dist>       Horizontal source-receiver distance in km, positive.\n"
"\n"
"    -Ds<depsrc>    Required source depth in km, nonnegative.\n"
"\n"
"    -Dr<deprcv>    Required receiver depth in km, nonnegative.\n"
"\n"
"    -A<azimuth>    Azimuth from source to receiver, in degree, [0, 360].\n"
"\n"
"    -S[u]<scale>   Source scale.  Moment sources use dyne-cm and force sources\n"
"                  use dyne.  `-Su` multiplies the scale by source shear modulus.\n"
"\n"
"    -M<strike>/<dip>[/<rake>]\n"
"                  Shear source, or tensile source when rake is omitted.\n"
"\n"
"    -T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>\n"
"                  Moment tensor source in North-East-Down coordinates.\n"
"\n"
"    -F<fn>/<fe>/<fz>\n"
"                  Single force in North-East-Down coordinates.\n"
"\n"
"                  When both the source and receiver are on the free surface,\n"
"                  only -F is supported; -e is ignored in this case.\n"
"\n"
"    -O<outdir>     Output directory.\n"
"\n"
"    -D<tftype>/<tfparams>\n"
"                  Convolve a time function. All time functions use area\n"
"                  normalization except Ricker wavelet, which has a peak\n"
"                  amplitude of 1.0.\n"
"                  There are several options:\n"
"                  + Parabolic wave (y = a*x^2 + b*x)\n"
"                    set -D%c/<t0>, <t0> (secs) is the duration of wave.\n", GRT_SIG_PARABOLA); printf(
"                    e.g.\n"
"                         -D%c/1.3\n", GRT_SIG_PARABOLA); printf(
"                  + Trapezoidal wave\n"
"                    set -D%c/<t1>/<t2>/<t3>, <t1> is the end time of\n", GRT_SIG_TRAPEZOID); printf(
"                    Rising, <t2> is the end time of Platform, and\n"
"                    <t3> is the end time of Falling.\n"
"                    e.g.\n"
"                         -D%c/0.1/0.2/0.4\n", GRT_SIG_TRAPEZOID); printf(
"                         -D%c/0.4/0.4/0.6 (become a triangle)\n", GRT_SIG_TRAPEZOID); printf(
"                  + Ricker wavelet\n"
"                    set -D%c/<f0>, <f0> (Hz) is the dominant frequency.\n", GRT_SIG_RICKER); printf(
"                    e.g.\n"
"                         -D%c/0.5\n", GRT_SIG_RICKER); printf(
"                  + Custom wave\n"
"                    set -D%c/<path>, <path> is the filepath to a custom\n", GRT_SIG_CUSTOM); printf(
"                    Time Function ASCII file. The file has just one column\n"
"                    of amplitude and no other columns. Its sequence sum should\n"
"                    be 1/dt, where dt is the sampling interval; the program\n"
"                    only issues a warning when it is not.\n"
"                    The file can contain unlimited comment lines with prefix\n"
"                    \"#\".\n"
"                    e.g.\n"
"                         -D%c/tfunc.txt\n", GRT_SIG_CUSTOM); printf(
"                  To match the physical time interval, parameters of the time\n"
"                  function may be slightly modified. The corresponding time\n"
"                  function is saved as a SAC file under <outdir>.\n"
"\n"
"    -E[p]<t0>[/<v0>]\n"
"                  Introduce a time shift in the output SAC records. The time\n"
"                  series starts at <t0> + distance/<v0>, where distance is the\n"
"                  straight-line source-receiver distance. <v0> is a reference\n"
"                  velocity in km/s; when omitted or zero, no distance correction\n"
"                  is applied.\n"
"                  -Ep<t0> starts the series at <t0> plus the direct P arrival\n"
"                  in the homogeneous halfspace. For example, use -Ep-10.\n"
"                  Without -E, the first time sample is the origin time.\n"
"\n"
"    -I<odr>        Apply odr time integrations after physical normalization.\n"
"\n"
"    -J<odr>        Apply odr time differentiations after physical normalization.\n"
"\n"
"    -n             Write receiver components as Z, N and E; default is Z, R and T.\n"
"\n"
"    -e             Also write spatial derivatives with direction prefixes matching\n"
"                  the output coordinates.\n"
"\n"
"    -s             Do not print completion information.\n"
"\n"
"    -h             Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt lamb -H8.0/4.62/3.3 -N6000/0.001 -R10 -A30 -S1e24 -Ds5 -Dr0 -M100/20/80 -e -n -Ores\n"
"\n\n\n"
);
}


/**
 * 解析 Lamb 模块的命令行选项
 *
 * @param[in,out]  Ctrl      Lamb 模块参数控制结构体
 * @param[in]      argc      命令行参数个数
 * @param[in]      argv      命令行参数数组
 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    Ctrl->computeType = GRT_SYN_EX;
    snprintf(Ctrl->s_computeType, sizeof(Ctrl->s_computeType), "%s", "EX");

    int opt;
    while ((opt = getopt(argc, argv, ":H:N:R:A:S:M:F:T:O:D:E:I:J:nesh")) != -1) {
        switch (opt) {
            /* 半空间参数 */
            case 'H': {
                char extra;
                real_t vp, vs, rho;
                real_t vp_vs_squared;
                if (sscanf(optarg, "%lf/%lf/%lf%c", &vp, &vs, &rho, &extra) != 3) {
                    GRTBadOptionError(H, "expected vp/vs/rho.");
                }
                if (!isfinite(vp) || !isfinite(vs) || !isfinite(rho) || vp <= 0.0 || vs <= 0.0 || rho <= 0.0) {
                    GRTBadOptionError(H, "vp, vs and rho should be positive.");
                }
                vp_vs_squared = GRT_SQUARE(vp / vs);
                Ctrl->H.nu = (vp_vs_squared - 2.0) / (2.0 * (vp_vs_squared - 1.0));
                if (!isfinite(Ctrl->H.nu) || Ctrl->H.nu <= 0.0 || Ctrl->H.nu >= 0.5) {
                    GRTBadOptionError(H, "vp/vs gives an invalid Poisson ratio.");
                }
                Ctrl->H.active = true;
                Ctrl->H.vp = vp;
                Ctrl->H.vs = vs;
                Ctrl->H.rho = rho;
                break;
            }

            /* 物理时间序列 */
            case 'N': {
                char extra;
                int nt;
                real_t dt;
                if (sscanf(optarg, "%d/%lf%c", &nt, &dt, &extra) != 2 || nt <= 0 || dt <= 0.0 || !isfinite(dt)) {
                    GRTBadOptionError(N, "expected positive nt/dt.");
                }
                Ctrl->N.active = true;
                Ctrl->N.nt = nt;
                Ctrl->N.dt = dt;
                break;
            }

            /* 水平震中距 */
            case 'R': {
                char extra;
                if (sscanf(optarg, "%lf%c", &Ctrl->R.dist, &extra) != 1 || Ctrl->R.dist <= 0.0 ||
                    !isfinite(Ctrl->R.dist)) {
                    GRTBadOptionError(R, "horizontal distance should be positive.");
                }
                Ctrl->R.active = true;
                break;
            }

            /* 方位角 */
            case 'A': {
                char extra;
                if (sscanf(optarg, "%lf%c", &Ctrl->A.azimuth, &extra) != 1 ||
                    !isfinite(Ctrl->A.azimuth) || Ctrl->A.azimuth < 0.0 || Ctrl->A.azimuth > 360.0) {
                    GRTBadOptionError(A, "azimuth should be in [0, 360].");
                }
                Ctrl->A.backazimuth = Ctrl->A.azimuth + 180.0;
                if (Ctrl->A.backazimuth >= 360.0) {
                    Ctrl->A.backazimuth -= 360.0;
                }
                Ctrl->A.azrad = Ctrl->A.azimuth * DEG1;
                Ctrl->A.active = true;
                break;
            }

            /* 源强缩放 */
            case 'S': {
                char extra;
                const char *scale = optarg;
                if (scale[0] == 'u') {
                    Ctrl->S.mult_src_mu = true;
                    ++scale;
                }
                if (sscanf(scale, "%lf%c", &Ctrl->S.M0, &extra) != 1 || !isfinite(Ctrl->S.M0)) {
                    GRTBadOptionError(S, "expected a numeric source scale.");
                }
                Ctrl->S.active = true;
                break;
            }

            /* 走向、倾角和滑动角 */
            case 'M': {
                if (Ctrl->M.active || Ctrl->F.active || Ctrl->T.active) {
                    GRTBadOptionError(M, "only one of -M, -T and -F may be used.");
                }
                char extra;
                real_t strike, dip, rake;
                int count = sscanf(optarg, "%lf/%lf/%lf%c", &strike, &dip, &rake, &extra);
                if (count != 2 && count != 3) {
                    GRTBadOptionError(M, "expected strike/dip[/rake].");
                }
                if (!isfinite(strike) || !isfinite(dip) || (count == 3 && !isfinite(rake)) ||
                    strike < 0.0 || strike > 360.0 || dip < 0.0 || dip > 90.0 ||
                    (count == 3 && (rake < -180.0 || rake > 180.0))) {
                    GRTBadOptionError(M, "strike, dip or rake is out of bound.");
                }
                Ctrl->mchn[0] = strike;
                Ctrl->mchn[1] = dip;
                Ctrl->mchn[2] = count == 3 ? rake : 0.0;
                Ctrl->computeType = count == 3 ? GRT_SYN_DC : GRT_SYN_TS;
                snprintf(Ctrl->s_computeType, sizeof(Ctrl->s_computeType), "%s",
                    count == 3 ? "DC" : "TS");
                Ctrl->M.active = true;
                break;
            }

            /* 单力源 */
            case 'F': {
                if (Ctrl->M.active || Ctrl->F.active || Ctrl->T.active) {
                    GRTBadOptionError(F, "only one of -M, -T and -F may be used.");
                }
                char extra;
                int count = sscanf(optarg, "%lf/%lf/%lf%c",
                    &Ctrl->mchn[0], &Ctrl->mchn[1], &Ctrl->mchn[2], &extra);
                if (count != 3 || !isfinite(Ctrl->mchn[0]) || !isfinite(Ctrl->mchn[1]) ||
                    !isfinite(Ctrl->mchn[2])) {
                    GRTBadOptionError(F, "expected fn/fe/fz.");
                }
                Ctrl->computeType = GRT_SYN_SF;
                snprintf(Ctrl->s_computeType, sizeof(Ctrl->s_computeType), "%s", "SF");
                Ctrl->F.active = true;
                break;
            }

            /* 矩张量源 */
            case 'T': {
                if (Ctrl->M.active || Ctrl->F.active || Ctrl->T.active) {
                    GRTBadOptionError(T, "only one of -M, -T and -F may be used.");
                }
                char extra;
                int count = sscanf(optarg, "%lf/%lf/%lf/%lf/%lf/%lf%c",
                    &Ctrl->mchn[0], &Ctrl->mchn[1], &Ctrl->mchn[2],
                    &Ctrl->mchn[3], &Ctrl->mchn[4], &Ctrl->mchn[5], &extra);
                if (count != 6 || !isfinite(Ctrl->mchn[0]) || !isfinite(Ctrl->mchn[1]) ||
                    !isfinite(Ctrl->mchn[2]) || !isfinite(Ctrl->mchn[3]) ||
                    !isfinite(Ctrl->mchn[4]) || !isfinite(Ctrl->mchn[5])) {
                    GRTBadOptionError(T, "expected six moment-tensor components.");
                }
                Ctrl->computeType = GRT_SYN_MT;
                snprintf(Ctrl->s_computeType, sizeof(Ctrl->s_computeType), "%s", "MT");
                Ctrl->T.active = true;
                break;
            }

            /* 输出目录 */
            case 'O':
                GRT_SAFE_FREE_PTR(Ctrl->O.s_output_dir);
                Ctrl->O.s_output_dir = strdup(optarg);
                Ctrl->O.active = true;
                break;

            /* 深度或时间函数 */
            case 'D':
                if (optarg[0] == 's' && optarg[1] != '/') {
                    char extra;
                    if (Ctrl->Depth.s_active) {
                        GRTBadOptionError(Ds, "the option is duplicated.");
                    }
                    if (sscanf(optarg + 1, "%lf%c", &Ctrl->Depth.depsrc, &extra) != 1 ||
                        !isfinite(Ctrl->Depth.depsrc) || Ctrl->Depth.depsrc < 0.0) {
                        GRTBadOptionError(Ds, "source depth should be nonnegative.");
                    }
                    Ctrl->Depth.s_active = true;
                } else if (optarg[0] == 'r' && optarg[1] != '/') {
                    char extra;
                    if (Ctrl->Depth.r_active) {
                        GRTBadOptionError(Dr, "the option is duplicated.");
                    }
                    if (sscanf(optarg + 1, "%lf%c", &Ctrl->Depth.deprcv, &extra) != 1 ||
                        !isfinite(Ctrl->Depth.deprcv) || Ctrl->Depth.deprcv < 0.0) {
                        GRTBadOptionError(Dr, "receiver depth should be nonnegative.");
                    }
                    Ctrl->Depth.r_active = true;
                } else {
                    if (optarg[0] == '\0' || optarg[1] != '/' || optarg[2] == '\0') {
                        GRTBadOptionError(D, "expected tftype/tfparams.");
                    }
                    GRT_SAFE_FREE_PTR(Ctrl->D.tfparams);
                    Ctrl->D.tftype = optarg[0];
                    Ctrl->D.tfparams = strdup(optarg + 2);
                    if (!grt_check_tftype_tfparams(Ctrl->D.tftype, Ctrl->D.tfparams)) {
                        GRTBadOptionError(D, "invalid time function.");
                    }
                    Ctrl->D.active = true;
                }
                break;

            /* 时间延迟 */
            case 'E': {
                Ctrl->E.active = true;
                if (optarg[0] == 'p') {
                    char extra;
                    real_t delayT0;
                    if (sscanf(optarg + 1, "%lf%c", &delayT0, &extra) != 1 || !isfinite(delayT0)) {
                        GRTBadOptionError(E, "expected a finite t0 after -Ep.");
                    }
                    if (delayT0 >= 0.0) {
                        GRTBadOptionError(E, "Can't set positive t0(%f) in -Ep.", delayT0);
                    }
                    Ctrl->E.delayT0 = delayT0;
                    Ctrl->E.delayV0 = 0.0;
                    Ctrl->E.refFirstP = true;
                } else {
                    char extra;
                    real_t delayT0;
                    real_t delayV0 = 0.0;
                    int count = sscanf(optarg, "%lf/%lf%c", &delayT0, &delayV0, &extra);
                    if ((count != 1 && count != 2) ||
                        (count == 1 && sscanf(optarg, "%lf%c", &delayT0, &extra) != 1) ||
                        !isfinite(delayT0) || !isfinite(delayV0)) {
                        GRTBadOptionError(E, "expected t0[/v0].");
                    }
                    if (delayV0 < 0.0) {
                        GRTBadOptionError(E, "Can't set negative v0(%f) in -E.", delayV0);
                    }
                    Ctrl->E.delayT0 = delayT0;
                    Ctrl->E.delayV0 = delayV0;
                    Ctrl->E.refFirstP = false;
                }
                break;
            }

            /* 时间积分次数 */
            case 'I': {
                char extra;
                if (sscanf(optarg, "%d%c", &Ctrl->I.int_times, &extra) != 1 || Ctrl->I.int_times <= 0) {
                    GRTBadOptionError(I, "order should be positive.");
                }
                Ctrl->I.active = true;
                break;
            }

            /* 时间微分次数 */
            case 'J': {
                char extra;
                if (sscanf(optarg, "%d%c", &Ctrl->J.dif_times, &extra) != 1 || Ctrl->J.dif_times <= 0) {
                    GRTBadOptionError(J, "order should be positive.");
                }
                Ctrl->J.active = true;
                break;
            }

            /* 是否旋转到 ZNE */
            case 'n':
                Ctrl->n.active = true;
                break;

            /* 是否计算空间导数 */
            case 'e':
                Ctrl->e.active = true;
                break;

            /* 是否静默输出 */
            case 's':
                Ctrl->s.active = true;
                break;

            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, H);
    if (!Ctrl->N.active) {
        GRTRaiseError("Need set options \"-N\". Use \"-h\" for help.\n");
    }
    GRTCheckOptionActive(Ctrl, R);
    if (!Ctrl->Depth.s_active) {
        GRTRaiseError("Need set options \"-Ds\". Use \"-h\" for help.\n");
    }
    if (!Ctrl->Depth.r_active) {
        GRTRaiseError("Need set options \"-Dr\". Use \"-h\" for help.\n");
    }
    GRTCheckOptionActive(Ctrl, A);
    GRTCheckOptionActive(Ctrl, S);
    GRTCheckOptionActive(Ctrl, O);
}


/**
 * 判断源点和接收点是否同时位于自由表面
 *
 * @param[in]      source_depth    源点深度，km
 * @param[in]      receiver_depth  接收点深度，km
 *
 * @return 同时位于 z=0 时返回 true，否则返回 false
 */
static bool is_surface_source_receiver(const real_t source_depth, const real_t receiver_depth)
{
    /* 只有源点和接收点同时位于 z=0 时才使用第一类 Lamb 解 */
    return source_depth == 0.0 && receiver_depth == 0.0;
}


/**
 * 判断源项是否为单力源
 *
 * @param[in]      source_index  震源在源项数组中的索引
 *
 * @return 源项为垂直力或水平力时返回 true，否则返回 false
 */
static bool source_is_force(const int source_index)
{
    /* 垂直力和水平力使用位移 Green 函数，其他源使用其空间导数 */
    return GRT_SRC_M_INDEX_IS_FORCE(source_index);
}


/**
 * 判断当前源项是否参与合成
 *
 * @param[in]      computeType   当前计算的震源类型
 * @param[in]      source_index  震源在源项数组中的索引
 * @param[in]      surface       源点和接收点是否均位于自由表面
 *
 * @return 需要输出当前源项时返回 true，否则返回 false
 */
static bool lamb_need_src(
    const GRT_SYN_TYPE computeType, const int source_index, const bool surface)
{
    /* 根据几何类型和用户选择的震源类型筛选需要写出的源 */
    if (surface) {
        return computeType == GRT_SYN_SF && source_is_force(source_index);
    }
    if (computeType == GRT_SYN_EX) {
        return source_index == GRT_SRC_M_EX_INDEX;
    }
    if (computeType == GRT_SYN_SF) {
        return source_is_force(source_index);
    }
    if (computeType == GRT_SYN_DC) {
        return source_index >= GRT_SRC_M_DD_INDEX;
    }
    return source_index == GRT_SRC_M_EX_INDEX || source_index >= GRT_SRC_M_DD_INDEX;
}


/** Lamb 模块使用的坐标变换矩阵 */
typedef struct {
    real_t local_from_global[3][3];
    real_t component_from_global[3][3];
} LAMB_COORDINATES;


/**
 * 创建 Lamb 求解器与输出所需的坐标变换矩阵
 *
 * @param[in]      azrad        源点到接收点的方位角，弧度
 * @param[in]      rot2ZNE      是否输出 Z、N、E 坐标
 * @param[out]     coordinates  坐标变换矩阵
 *
 */
static void make_lamb_coordinates(
    const real_t azrad, const bool rot2ZNE, LAMB_COORDINATES *coordinates)
{
    const real_t sin_azimuth = sin(azrad);
    const real_t cos_azimuth = cos(azrad);
    coordinates->local_from_global[0][0] = cos_azimuth;
    coordinates->local_from_global[0][1] = sin_azimuth;
    coordinates->local_from_global[0][2] = 0.0;
    coordinates->local_from_global[1][0] = -sin_azimuth;
    coordinates->local_from_global[1][1] = cos_azimuth;
    coordinates->local_from_global[1][2] = 0.0;
    coordinates->local_from_global[2][0] = 0.0;
    coordinates->local_from_global[2][1] = 0.0;
    coordinates->local_from_global[2][2] = 1.0;

    if (rot2ZNE) {
        memcpy(coordinates->component_from_global, (real_t[3][3]){
            {0.0, 0.0, -1.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
        }, sizeof(coordinates->component_from_global));
    } else {
        memcpy(coordinates->component_from_global, (real_t[3][3]){
            {0.0, 0.0, -1.0},
            {cos_azimuth, sin_azimuth, 0.0},
            {-sin_azimuth, cos_azimuth, 0.0},
        }, sizeof(coordinates->component_from_global));
    }
}


/** 计算单力源与位移 Green 函数的收缩结果
 *
 * @param[in]      G                       位移 Green 函数
 * @param[in]      component_from_global  输出分量坐标变换矩阵
 * @param[in]      output_component       输出分量索引
 * @param[in]      source                 全局坐标下的源力向量
 *
 * @return 收缩后的 Green 函数值
 */
static real_t evaluate_force(
    const real_t G[3][3], const real_t component_from_global[3][3],
    const int output_component, const real_t source[3])
{
    real_t value = 0.0;
    for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
        value += component_from_global[output_component][i] * G[i][j] * source[j];
    }}
    return value;
}


/** 计算矩源与源点一阶导数 Green 函数的收缩结果
 *
 * @param[in]      dG                      源点一阶导数 Green 函数
 * @param[in]      component_from_global  输出分量坐标变换矩阵
 * @param[in]      output_component       输出分量索引
 * @param[in]      source                  全局坐标下的源矩张量
 *
 * @return 收缩后的 Green 函数值
 */
static real_t evaluate_moment(
    const real_t dG[3][3][3], const real_t component_from_global[3][3],
    const int output_component, const real_t source[3][3])
{
    real_t value = 0.0;
    for (int k = 0; k < 3; ++k) {
    for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
        value += component_from_global[output_component][i] * dG[k][i][j] * source[j][k];
    }}}
    return value;
}


/** 计算单力源接收点空间导数与源力的收缩结果
 *
 * @param[in]      dG                      接收点一阶导数 Green 函数
 * @param[in]      component_from_global  坐标变换矩阵
 * @param[in]      derivative_direction   导数方向索引
 * @param[in]      output_component       输出分量索引
 * @param[in]      source                  全局坐标下的源力向量
 *
 * @return 收缩后的 Green 函数值
 */
static real_t evaluate_force_derivative(
    const real_t dG[3][3][3], const real_t component_from_global[3][3],
    const int derivative_direction, const int output_component, const real_t source[3])
{
    real_t value = 0.0;
    for (int k = 0; k < 3; ++k) {
    for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
        value += component_from_global[derivative_direction][k] *
            component_from_global[output_component][i] * dG[k][i][j] * source[j];
    }}}
    return value;
}


/** 计算矩源混合空间导数与源矩的收缩结果
 *
 * @param[in]      dG                      混合二阶导数 Green 函数
 * @param[in]      component_from_global  坐标变换矩阵
 * @param[in]      derivative_direction   接收点导数方向索引
 * @param[in]      output_component       输出分量索引
 * @param[in]      source                  全局坐标下的源矩张量
 *
 * @return 收缩后的 Green 函数值
 */
static real_t evaluate_moment_derivative(
    const real_t dG[3][3][3][3], const real_t component_from_global[3][3],
    const int derivative_direction, const int output_component, const real_t source[3][3])
{
    real_t value = 0.0;
    for (int k = 0; k < 3; ++k) {
    for (int kp = 0; kp < 3; ++kp) {
    for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
        value += component_from_global[derivative_direction][k] *
            component_from_global[output_component][i] * dG[k][kp][i][j] * source[j][kp];
    }}}}
    return value;
}


/**
 * 计算当前距离阶数下各源项的辐射系数
 *
 * @param[in]      computeType   当前计算的震源类型
 * @param[in]      M0            源强缩放因子
 * @param[in]      nu            半空间泊松比
 * @param[in]      azrad         源点到接收点的方位角，弧度
 * @param[in]      mchn          震源机制参数数组
 * @param[in]      par_theta     是否计算方位角导数
 * @param[in]      coef          当前计算项的距离缩放因子
 * @param[out]     srcRadi       各源项和各输出分量的辐射系数
 */
static void make_source_radiation(
    const GRT_SYN_TYPE computeType, const real_t M0,
    const real_t nu, const real_t azrad, const real_t mchn[GRT_MECHANISM_NUM],
    const bool par_theta, const real_t coef, realChnlGrid srcRadi)
{
    /* 先清零所有源项，再根据源类型生成当前距离阶数的辐射系数 */
    memset(srcRadi, 0, sizeof(realChnlGrid));
    /* grt_set_source_radiation 使用 vp/vs 作为水平分量和垂直分量的换算比 */
    grt_set_source_radiation(
        srcRadi, computeType, par_theta, M0, coef,
        sqrt(2.0 * (1.0 - nu) / (1.0 - 2.0 * nu)), azrad, mchn);
}


/**
 * 计算 Lamb 位移 Green 函数及按需计算其空间导数
 *
 * @param[in]      nu                    半空间泊松比
 * @param[in]      tbar                  无量纲时间序列
 * @param[in]      nt                    时间序列长度
 * @param[in]      horizontal_distance   源点与接收点的水平距离，km
 * @param[in]      source_depth          源点深度，km
 * @param[in]      receiver_depth        接收点深度，km
 * @param[in]      azimuth_degree        源点到接收点的方位角，度
 * @param[in]      computeType           当前计算的震源类型
 * @param[in]      calculate_derivatives 是否计算空间导数
 * @param[out]     result                Lamb Green 函数及其导数的结果结构体
 */
static void make_lamb_result(
    const real_t nu, const real_t *tbar, const int nt, const real_t horizontal_distance,
    const real_t source_depth, const real_t receiver_depth, const real_t azimuth_degree,
    const GRT_SYN_TYPE computeType, const bool calculate_derivatives,
    LAMB_RESULT *result)
{
    /* 只有需要空间导数的源才分配相应的导数数组 */
    const bool surface = is_surface_source_receiver(source_depth, receiver_depth);
    const bool moment = computeType != GRT_SYN_SF;
    const bool force = computeType == GRT_SYN_SF;
    const bool need_upar = calculate_derivatives && !surface;
    const size_t nt_size = (size_t)nt;
    result->G = GRT_SAFE_CALLOC(nt_size, sizeof(*result->G));
    result->dG_source = moment ? GRT_SAFE_CALLOC(nt_size, sizeof(*result->dG_source)) : NULL;
    result->dG_receiver = need_upar && force ? GRT_SAFE_CALLOC(nt_size, sizeof(*result->dG_receiver)) : NULL;
    result->dG_mixed = need_upar && moment ? GRT_SAFE_CALLOC(nt_size, sizeof(*result->dG_mixed)) : NULL;

    /* 地表、单侧地下和双侧地下分别对应三类 Lamb 求解器 */
    if (surface) {
        grt_solve_lamb1(nu, tbar, nt, azimuth_degree, result->G);
    } else if (source_depth > 0.0 && receiver_depth > 0.0) {
        grt_solve_lamb3(nu, tbar, nt, horizontal_distance,
            source_depth, receiver_depth, azimuth_degree,
            result->G, result->dG_source, result->dG_receiver, result->dG_mixed);
    } else {
        grt_solve_lamb2(nu, tbar, nt, horizontal_distance,
            source_depth, receiver_depth, azimuth_degree,
            result->G, result->dG_source, result->dG_receiver, result->dG_mixed);
    }
}


/**
 * 创建并初始化一个 Lamb 模块的 SAC 记录
 *
 * @param[in]      nt                   记录采样点数
 * @param[in]      dt                   记录采样间隔，s
 * @param[in]      horizontal_distance  源点与接收点的水平距离，km
 * @param[in]      source_depth         源点深度，km
 * @param[in]      receiver_depth       接收点深度，km
 * @param[in]      azimuth              源点到接收点的方位角，度
 * @param[in]      begin_time           SAC 记录起始时刻，s
 * @param[in]      vp                   半空间 P 波速度，km/s
 * @param[in]      vs                   半空间 S 波速度，km/s
 * @param[in]      rho                  半空间密度，g/cm^3
 *
 * @return 新创建的 SAC 记录
 */
static SACTRACE *new_lamb_trace(
    const int nt, const real_t dt, const real_t horizontal_distance,
    const real_t source_depth, const real_t receiver_depth,
    const real_t azimuth, const real_t begin_time,
    const real_t vp, const real_t vs, const real_t rho)
{
    /* 头段字段与 greenfn 输出的 SAC 原型保持一致 */
    SACTRACE *sac = grt_new_SACTRACE((float)dt, nt, (float)begin_time);
    sac->hd.o = 0.0;
    sac->hd.iztype = IO;
    sac->hd.dist = horizontal_distance;
    sac->hd.evdp = source_depth;
    sac->hd.stel = -receiver_depth * 1e3;
    /* 直接时域解不使用 greenfn 的虚频率补偿，也不包含衰减 */
    sac->hd.user0 = 0.0;
    sac->hd.user1 = vp;
    sac->hd.user2 = vs;
    sac->hd.user3 = rho;
    sac->hd.user4 = 0.0;
    sac->hd.user5 = 0.0;
    sac->hd.user6 = vp;
    sac->hd.user7 = vs;
    sac->hd.user8 = rho;

    sac->hd.az = azimuth;
    sac->hd.baz = azimuth + 180.0;
    if (sac->hd.baz >= 360.0) {
        sac->hd.baz -= 360.0;
    }
    return sac;
}


enum {
    LAMB_PHASE_P,                         ///< t0/kt0：直达 P 波
    LAMB_PHASE_S,                         ///< t1/kt1：直达 S 波
    LAMB_PHASE_R,                         ///< t2/kt2：Rayleigh 波参考到时
    LAMB_PHASE_sP,                        ///< t3/kt3：滑行 sP 波
    LAMB_PHASE_PP,                        ///< t4/kt4：反射 PP 波
    LAMB_PHASE_SS,                        ///< t5/kt5：反射 SS 波
    LAMB_PHASE_PS,                        ///< t6/kt6：PS 转换波
    LAMB_PHASE_SP,                        ///< t7/kt7：SP 转换波
    LAMB_PHASE_sPs,                       ///< t8/kt8：滑行 sPs 波
    LAMB_PHASE_COUNT,                     ///< 已定义的 Lamb 震相数量
    LAMB_SAC_PICK_COUNT = 10,             ///< SAC 用户震相槽位总数
};


/** 按 SAC 的 8 字节定长格式写入震相名称 */
static void copy_lamb_phase_name(char target[9], const char *name)
{
    const size_t length = GRT_MIN(strlen(name), (size_t)SAC_HEADER_STRING_LENGTH_FILE);
    memset(target, ' ', SAC_HEADER_STRING_LENGTH_FILE);
    memcpy(target, name, length);
    target[SAC_HEADER_STRING_LENGTH_FILE] = '\0';
}


/** 清空 Lamb SAC 记录中的震相到时及名称 */
static void clear_lamb_arrivals(SACTRACE *sac)
{
    float *times[LAMB_SAC_PICK_COUNT] = {
        &sac->hd.t0, &sac->hd.t1, &sac->hd.t2, &sac->hd.t3,
        &sac->hd.t4, &sac->hd.t5, &sac->hd.t6, &sac->hd.t7,
        &sac->hd.t8, &sac->hd.t9,
    };
    char *names[LAMB_SAC_PICK_COUNT] = {
        sac->hd.kt0, sac->hd.kt1, sac->hd.kt2, sac->hd.kt3,
        sac->hd.kt4, sac->hd.kt5, sac->hd.kt6, sac->hd.kt7,
        sac->hd.kt8, sac->hd.kt9,
    };
    for (int i = 0; i < LAMB_SAC_PICK_COUNT; ++i) {
        *times[i] = SAC_FLOAT_UNDEF;
        copy_lamb_phase_name(names[i], SAC_CHAR8_UNDEF);
    }
}


/** 将一个无量纲震相到时写入 SAC 头段 */
static void set_lamb_arrival(
    SACTRACE *sac, const int phase, const real_t tbar, const real_t time_scale, const char *name)
{
    if (phase < 0 || phase >= LAMB_PHASE_COUNT || tbar < 0.0 || !isfinite(tbar)) {
        return;
    }
    float *times[LAMB_SAC_PICK_COUNT] = {
        &sac->hd.t0, &sac->hd.t1, &sac->hd.t2, &sac->hd.t3,
        &sac->hd.t4, &sac->hd.t5, &sac->hd.t6, &sac->hd.t7,
        &sac->hd.t8, &sac->hd.t9,
    };
    char *names[LAMB_SAC_PICK_COUNT] = {
        sac->hd.kt0, sac->hd.kt1, sac->hd.kt2, sac->hd.kt3,
        sac->hd.kt4, sac->hd.kt5, sac->hd.kt6, sac->hd.kt7,
        sac->hd.kt8, sac->hd.kt9,
    };
    *times[phase] = (float)(tbar * time_scale);
    copy_lamb_phase_name(names[phase], name);
}


/** 根据源点和接收点位置写入 Lamb SAC 记录的震相到时 */
static void set_lamb_arrivals(
    SACTRACE *sac, const real_t nu, const real_t horizontal_distance,
    const real_t source_depth, const real_t receiver_depth,
    const real_t direct_distance, const real_t vs)
{
    clear_lamb_arrivals(sac);
    const real_t time_scale = direct_distance / vs;
    real_t tP;
    real_t tR;
    set_lamb_arrival(sac, LAMB_PHASE_S, 1.0, time_scale, "S");
    grt_compute_lamb1_travt(nu, &tP, &tR);

    if (source_depth > 0.0 && receiver_depth > 0.0) {
        const real_t reflected_distance = hypot(horizontal_distance, source_depth + receiver_depth);
        tR *= reflected_distance / direct_distance;
    }
    set_lamb_arrival(sac, LAMB_PHASE_R, tR, time_scale, "R");

    if (is_surface_source_receiver(source_depth, receiver_depth)) {
        set_lamb_arrival(sac, LAMB_PHASE_P, tP, time_scale, "P");
        return;
    }

    if (source_depth > 0.0 && receiver_depth > 0.0) {
        real_t tPP;
        real_t tSS;
        real_t tPS;
        real_t tSP;
        real_t t_sPs;
        grt_compute_lamb3_travt(
            nu, horizontal_distance, source_depth, receiver_depth,
            &tP, &tPP, &tSS, &tPS, &tSP, &t_sPs);
        set_lamb_arrival(sac, LAMB_PHASE_P, tP, time_scale, "P");
        set_lamb_arrival(sac, LAMB_PHASE_PP, tPP, time_scale, "PP");
        set_lamb_arrival(sac, LAMB_PHASE_SS, tSS, time_scale, "SS");
        set_lamb_arrival(sac, LAMB_PHASE_PS, tPS, time_scale, "PS");
        set_lamb_arrival(sac, LAMB_PHASE_SP, tSP, time_scale, "SP");
        set_lamb_arrival(sac, LAMB_PHASE_sPs, t_sPs, time_scale, "sPs");
        return;
    }

    real_t t_sP;
    grt_compute_lamb2_travt(
        nu, horizontal_distance, source_depth, receiver_depth, &tP, &t_sP);
    set_lamb_arrival(sac, LAMB_PHASE_P, tP, time_scale, "P");
    set_lamb_arrival(sac, LAMB_PHASE_sP, t_sP, time_scale, "sP");
}


/**
 * 去除 Lamb 闭合解自带的一次时间积分
 *
 * @param[in,out]  sac     待恢复的 SAC 记录
 */
static void lamb_normalize_trace(SACTRACE *sac)
{
    /* 闭合解和混合导数都带有一次时间积分，先恢复为脉冲型物理解 */
    if (sac->hd.npts <= 1) {
        sac->data[0] = 0.0f;
        return;
    }
    grt_differential(sac->data, sac->hd.npts, sac->hd.delta);
}


/**
 * 对 SAC 记录执行时间函数卷积及时间积分或微分
 *
 * @param[in,out]  sac                待处理的 SAC 记录
 * @param[in]      integration_order  时间积分次数
 * @param[in]      differential_order 时间微分次数
 * @param[in]      time_function      可选的时间函数记录，NULL 表示不进行卷积
 */
static void lamb_postprocess_trace(
    SACTRACE *sac, const int integration_order, const int differential_order,
    const SACTRACE *time_function)
{
    /* 单点序列无法进行稳定的差分，直接保留零值 */
    if (sac->hd.npts <= 1) {
        sac->data[0] = 0.0f;
        return;
    }

    /* 时间函数卷积先于用户要求的积分和微分 */
    if (time_function != NULL) {
        float *convolution = GRT_SAFE_CALLOC(sac->hd.npts, sizeof(*convolution));
        grt_oaconvolve(sac->data, sac->hd.npts, time_function->data, time_function->hd.npts,
            convolution, sac->hd.npts, false);
        /* 时间函数样本表示物理时间函数，连续卷积的离散积分因子为 dt */
        for (int n = 0; n < sac->hd.npts; ++n) {
            sac->data[n] = convolution[n] * sac->hd.delta;
        }
        GRT_SAFE_FREE_PTR(convolution);
    }

    for (int i = 0; i < integration_order; ++i) {
        grt_trap_integral(sac->data, sac->hd.npts, sac->hd.delta);
    }
    for (int i = 0; i < differential_order; ++i) {
        grt_differential(sac->data, sac->hd.npts, sac->hd.delta);
    }
}


/**
 * 按模块命名规则写出一个 Lamb SAC 记录
 *
 * @param[in]      output_path       输出目录
 * @param[in]      derivative_prefix 空间导数方向前缀
 * @param[in]      channel           接收分量名称
 * @param[in,out]  sac               待写出的 SAC 记录
 */
static void save_to_sac(
    const char *output_path, const char *derivative_prefix,
    const char channel, SACTRACE *sac)
{
    /* 文件名由导数方向和接收分量组成，与 syn 模块一致 */
    char component[32];
    char *path = NULL;
    snprintf(component, sizeof(component), "%s%c", derivative_prefix, channel);
    snprintf(sac->hd.kcmpnm, sizeof(sac->hd.kcmpnm), "%.8s", component);
    GRT_SAFE_ASPRINTF(&path, "%s/%s.sac", output_path, component);
    grt_write_SACTRACE(path, sac);
    GRT_SAFE_FREE_PTR(path);
}


/**
 * 将一个 SAC 记录叠加到目标记录
 *
 * @param[in,out]  target  目标 SAC 记录
 * @param[in]      source  待叠加的 SAC 记录
 */
static void accumulate_trace(SACTRACE *target, const SACTRACE *source)
{
    /* 将一个源分量的记录叠加到最终合成记录 */
    for (int n = 0; n < target->hd.npts; ++n) {
        target->data[n] += source->data[n];
    }
}


/** 创建一组 SAC 记录及其空间导数记录
 *
 * @param[in]      prototype   SAC 记录原型
 * @param[in]      calc_upar   是否创建空间导数记录
 * @param[out]     base        三个基本接收分量的 SAC 记录
 * @param[out]     derivative  三个导数方向和三个接收分量的 SAC 记录
 */
static void allocate_lamb_traces(
    SACTRACE *prototype, const bool calc_upar, SACTRACE *base[3], SACTRACE *derivative[3][3])
{
    for (int component = 0; component < 3; ++component) {
        base[component] = grt_copy_SACTRACE(prototype, true);
        if (calc_upar) {
            for (int direction = 0; direction < 3; ++direction) {
                derivative[direction][component] = grt_copy_SACTRACE(prototype, true);
            }
        }
    }
}


/** 释放一组 SAC 记录及其空间导数记录
 *
 * @param[in]      calc_upar   是否释放空间导数记录
 * @param[in,out]  base        三个基本接收分量的 SAC 记录
 * @param[in,out]  derivative  三个导数方向和三个接收分量的 SAC 记录
 */
static void free_lamb_traces(
    const bool calc_upar, SACTRACE *base[3], SACTRACE *derivative[3][3])
{
    for (int component = 0; component < 3; ++component) {
        grt_free_SACTRACE(base[component]);
        if (calc_upar) {
            for (int direction = 0; direction < 3; ++direction) {
                grt_free_SACTRACE(derivative[direction][component]);
            }
        }
    }
}


/** 对一组 SAC 记录执行 Lamb 归一化、卷积、积分和微分
 *
 * @param[in]      calc_upar      是否处理空间导数记录
 * @param[in]      int_times      时间积分次数
 * @param[in]      dif_times      时间微分次数
 * @param[in]      time_function  可选的时间函数记录
 * @param[in,out]  base           三个基本接收分量的 SAC 记录
 * @param[in,out]  derivative     三个导数方向和三个接收分量的 SAC 记录
 */
static void postprocess_lamb_traces(
    const bool calc_upar, const int int_times, const int dif_times,
    const SACTRACE *time_function, SACTRACE *base[3], SACTRACE *derivative[3][3])
{
    for (int component = 0; component < 3; ++component) {
        lamb_normalize_trace(base[component]);
        if (calc_upar) {
            for (int direction = 0; direction < 3; ++direction) {
                lamb_normalize_trace(derivative[direction][component]);
            }
        }
    }

    for (int component = 0; component < 3; ++component) {
        lamb_postprocess_trace(base[component], int_times, dif_times, time_function);
        if (calc_upar) {
            for (int direction = 0; direction < 3; ++direction) {
                lamb_postprocess_trace(derivative[direction][component], int_times, dif_times, time_function);
            }
        }
    }
}


/** 将一组源项记录累加到合成记录
 *
 * @param[in]      calc_upar          是否累加空间导数记录
 * @param[in]      source_base        当前源项的基本接收分量记录
 * @param[in]      source_derivative 当前源项的空间导数记录
 * @param[in,out]  base               合成的基本接收分量记录
 * @param[in,out]  derivative         合成的空间导数记录
 */
static void accumulate_lamb_traces(
    const bool calc_upar, SACTRACE *source_base[3], SACTRACE *source_derivative[3][3],
    SACTRACE *base[3], SACTRACE *derivative[3][3])
{
    for (int component = 0; component < 3; ++component) {
        accumulate_trace(base[component], source_base[component]);
        if (calc_upar) {
            for (int direction = 0; direction < 3; ++direction) {
                accumulate_trace(derivative[direction][component], source_derivative[direction][component]);
            }
        }
    }
}


/**
 * 将单力源辐射系数转换为源力向量
 *
 * @param[in]      srcRadi        各源项和各分量的辐射系数
 * @param[in]      source_index   震源在源项数组中的索引
 * @param[out]     source         R、T、Z_down 坐标下的源力向量
 */
static void make_source_force(const realChnlGrid srcRadi, const int source_index, real_t source[3])
{
    /* 源辐射在内部柱坐标中按 R、T、Z_down 排列 */
    memset(source, 0, sizeof(real_t) * 3);
    if (source_index == GRT_SRC_M_VF_INDEX) {
        source[2] = srcRadi[source_index][0];
    } else if (source_index == GRT_SRC_M_HF_INDEX) {
        source[0] = srcRadi[source_index][0];
        source[1] = srcRadi[source_index][2];
    }
}


/**
 * 将矩源辐射系数转换为源矩张量
 *
 * @param[in]      srcRadi        各源项和各分量的辐射系数
 * @param[in]      source_index   震源在源项数组中的索引
 * @param[out]     source         R、T、Z_down 坐标下的源矩张量
 */
static void make_source_moment(const realChnlGrid srcRadi, const int source_index, real_t source[3][3])
{
    /* 将 EX、DD、DS 和 SS 的辐射系数还原为 R、T、Z_down 下的矩张量 */
    memset(source, 0, sizeof(real_t) * 3 * 3);
    if (source_index == GRT_SRC_M_EX_INDEX) {
        source[0][0] = source[1][1] = source[2][2] = srcRadi[source_index][0];
    } else if (source_index == GRT_SRC_M_DD_INDEX) {
        source[0][0] = source[1][1] = -srcRadi[source_index][0];
        source[2][2] = 2.0 * srcRadi[source_index][0];
    } else if (source_index == GRT_SRC_M_DS_INDEX) {
        source[0][2] = source[2][0] = -srcRadi[source_index][0];
        source[1][2] = source[2][1] = -srcRadi[source_index][2];
    } else if (source_index == GRT_SRC_M_SS_INDEX) {
        source[0][0] = srcRadi[source_index][0];
        source[1][1] = -srcRadi[source_index][0];
        source[0][1] = source[1][0] = srcRadi[source_index][2];
    }
}


/** 将局部坐标下的源力转换为全局坐标
 *
 * @param[in]      local_from_global  全局坐标到局部坐标的变换矩阵
 * @param[in]      source_local       局部坐标下的源力
 * @param[out]     source_global      全局坐标下的源力
 */
static void transform_source_force(
    const real_t local_from_global[3][3], const real_t source_local[3], real_t source_global[3])
{
    memset(source_global, 0, sizeof(real_t) * 3);
    for (int i = 0; i < 3; ++i) {
    for (int a = 0; a < 3; ++a) {
        source_global[i] += local_from_global[a][i] * source_local[a];
    }}
}


/** 将局部坐标下的源矩转换为全局坐标
 *
 * @param[in]      local_from_global  全局坐标到局部坐标的变换矩阵
 * @param[in]      source_local       局部坐标下的源矩
 * @param[out]     source_global      全局坐标下的源矩
 */
static void transform_source_moment(
    const real_t local_from_global[3][3], const real_t source_local[3][3], real_t source_global[3][3])
{
    memset(source_global, 0, sizeof(real_t) * 3 * 3);
    for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
    for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
        source_global[i][j] += local_from_global[a][i] * local_from_global[b][j] * source_local[a][b];
    }}}}
}


/** 将一个源项转换为全局坐标下的源力或源矩
 *
 * @param[in]      srcRadi            各源项和各分量的辐射系数
 * @param[in]      source_index       震源在源项数组中的索引
 * @param[in]      force              当前源项是否为单力源
 * @param[in]      local_from_global  全局坐标到局部坐标的变换矩阵
 * @param[out]     source_force      全局坐标下的源力
 * @param[out]     source_moment     全局坐标下的源矩
 */
static void make_source_terms(
    const realChnlGrid srcRadi, const int source_index, const bool force,
    const real_t local_from_global[3][3], real_t source_force[3], real_t source_moment[3][3])
{
    memset(source_force, 0, sizeof(real_t) * 3);
    memset(source_moment, 0, sizeof(real_t) * 3 * 3);
    if (force) {
        real_t source_local[3];
        make_source_force(srcRadi, source_index, source_local);
        transform_source_force(local_from_global, source_local, source_force);
    } else {
        real_t source_local[3][3];
        make_source_moment(srcRadi, source_index, source_local);
        transform_source_moment(local_from_global, source_local, source_moment);
    }
}


/** 将一个源项合成为请求坐标下的记录及其空间导数
 *
 * @param[in]      rot2ZNE                  是否输出 Z、N、E 坐标
 * @param[in]      horizontal_distance      源点与接收点的水平距离，km
 * @param[in]      computeType              当前计算的震源类型
 * @param[in]      M0                       源强缩放因子
 * @param[in]      nu                       半空间泊松比
 * @param[in]      azrad                    源点到接收点的方位角，弧度
 * @param[in]      mchn                     震源机制参数数组
 * @param[in]      calc_upar                是否计算空间导数
 * @param[in]      nt                       时间序列长度
 * @param[in]      result                   Lamb Green 函数及其导数
 * @param[in]      source_index             震源在源项数组中的索引
 * @param[in]      force_factor             单力源位移的物理归一化因子
 * @param[in]      moment_factor            矩源位移的物理归一化因子
 * @param[in]      force_derivative_factor  单力源空间导数的物理归一化因子
 * @param[in]      moment_derivative_factor 矩源空间导数的物理归一化因子
 * @param[out]     base                     三个基本接收分量的 SAC 记录
 * @param[out]     derivative               三个导数方向和三个接收分量的 SAC 记录
 */
static void fill_source_traces(
    const bool rot2ZNE, const real_t horizontal_distance,
    const GRT_SYN_TYPE computeType, const real_t M0, const real_t nu,
    const real_t azrad, const real_t mchn[GRT_MECHANISM_NUM],
    const bool calc_upar, const int nt, const LAMB_RESULT *result,
    const int source_index, const real_t force_factor, const real_t moment_factor,
    const real_t force_derivative_factor, const real_t moment_derivative_factor,
    SACTRACE *base[3], SACTRACE *derivative[3][3])
{
    /* 统一在全局 N、E、Z_down 坐标中收缩 Green 函数 */
    const bool force = source_is_force(source_index);
    const real_t base_factor = force ? force_factor : moment_factor;
    const real_t derivative_factor = force ? force_derivative_factor : moment_derivative_factor;
    const bool skip_transverse = !rot2ZNE && GRT_SRC_M_ORDERS[source_index] == 0;
    LAMB_COORDINATES coordinates;
    realChnlGrid baseRadiation;
    make_lamb_coordinates(azrad, rot2ZNE, &coordinates);
    make_source_radiation(computeType, M0, nu, azrad, mchn, false, 1.0, baseRadiation);

    real_t source_force_global[3];
    real_t source_moment_global[3][3];
    make_source_terms(baseRadiation, source_index, force, coordinates.local_from_global,
        source_force_global, source_moment_global);

    for (int output_component = 0; output_component < 3; ++output_component) {
        if (skip_transverse && output_component == 2) {
            continue;
        }
        for (int n = 0; n < nt; ++n) {
            const real_t value = force
                ? evaluate_force(result->G[n], coordinates.component_from_global, output_component, source_force_global)
                : evaluate_moment(result->dG_source[n], coordinates.component_from_global, output_component, source_moment_global);
            base[output_component]->data[n] = (float)(base_factor * value);
        }
    }

    if (!calc_upar) {
        return;
    }

    realChnlGrid receiverRadiation;
    make_source_radiation(computeType, M0, nu, azrad, mchn, false, 1e-5, receiverRadiation);
    make_source_terms(receiverRadiation, source_index, force, coordinates.local_from_global,
        source_force_global, source_moment_global);

    realChnlGrid thetaRadiation;
    real_t theta_force_global[3] = {0.0};
    real_t theta_moment_global[3][3] = {{0.0}};
    if (!rot2ZNE) {
        make_source_radiation(computeType, M0, nu, azrad, mchn, true,
            1e-5 / horizontal_distance, thetaRadiation);
        make_source_terms(thetaRadiation, source_index, force, coordinates.local_from_global,
            theta_force_global, theta_moment_global);
    }

    for (int direction = 0; direction < 3; ++direction) {
        for (int output_component = 0; output_component < 3; ++output_component) {
            if (skip_transverse && output_component == 2) {
                continue;
            }
            for (int n = 0; n < nt; ++n) {
                real_t value;
                real_t scale;
                if (!rot2ZNE && direction == 2) {
                    value = force
                        ? evaluate_force(result->G[n], coordinates.component_from_global, output_component, theta_force_global)
                        : evaluate_moment(result->dG_source[n], coordinates.component_from_global, output_component, theta_moment_global);
                    scale = base_factor;
                } else {
                    value = force
                        ? evaluate_force_derivative(result->dG_receiver[n], coordinates.component_from_global,
                            direction, output_component, source_force_global)
                        : evaluate_moment_derivative(result->dG_mixed[n], coordinates.component_from_global,
                            direction, output_component, source_moment_global);
                    scale = derivative_factor;
                }
                derivative[direction][output_component]->data[n] = (float)(scale * value);
            }
        }
    }
}


/** 模块主函数 */
int lamb_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);

    const real_t depsrc = Ctrl->Depth.depsrc;
    const real_t deprcv = Ctrl->Depth.deprcv;
    const real_t horizontal_distance = Ctrl->R.dist;
    const bool surface = is_surface_source_receiver(depsrc, deprcv);
    const bool rot2ZNE = Ctrl->n.active;
    const bool calc_upar = Ctrl->e.active && !surface;
    const char *chs = rot2ZNE ? GRT_ZNE_CODES : GRT_ZRT_CODES;
    if (surface) {
        /* 第一类 Lamb 解只支持由 -F 指定的地表单力源 */
        if (Ctrl->computeType != GRT_SYN_SF) {
            GRTRaiseError(
                "When both source and receiver are on the free surface, only the single force source specified by -F is supported.\n");
        }
        if (Ctrl->e.active) {
            GRTRaiseWarning("Both source and receiver are on the free surface; -e is ignored.");
        }
    }

    GRTCheckMakeDir(Ctrl->O.s_output_dir);
    const real_t distance = hypot(horizontal_distance, depsrc - deprcv);
    real_t begin_time = Ctrl->E.delayT0;
    if (Ctrl->E.refFirstP) {
        begin_time += distance / Ctrl->H.vp;
    } else if (Ctrl->E.delayV0 > 0.0) {
        begin_time += distance / Ctrl->E.delayV0;
    }
    /* 按记录的实际物理时刻构造无量纲时间 tbar=vs*t/r */
    Ctrl->N.tbar = GRT_SAFE_CALLOC((size_t)Ctrl->N.nt, sizeof(*Ctrl->N.tbar));
    for (int n = 0; n < Ctrl->N.nt; ++n) {
        Ctrl->N.tbar[n] = (begin_time + n * Ctrl->N.dt) * Ctrl->H.vs / distance;
    }

    LAMB_RESULT result = {0};
    make_lamb_result(
        Ctrl->H.nu, Ctrl->N.tbar, Ctrl->N.nt, horizontal_distance, depsrc, deprcv,
        Ctrl->A.azimuth, Ctrl->computeType, calc_upar, &result);

    /* 将无量纲闭合解恢复为物理量，导数阶数每增加一阶再除以一个 r */
    const real_t mu = Ctrl->H.vs * Ctrl->H.vs * Ctrl->H.rho;
    const real_t factor = 1.0 / (PI * PI * mu);
    const real_t force_factor = factor / distance;
    const real_t moment_factor = factor / (distance * distance);
    const real_t force_derivative_factor = factor / (distance * distance);
    const real_t moment_derivative_factor = factor / (distance * distance * distance);
    real_t source_scale = Ctrl->S.M0;
    SACTRACE *time_function = NULL;

    if (Ctrl->S.mult_src_mu) {
        /* -Su 的输入单位是 potency，换算为带剪切模量的源强 */
        Ctrl->S.src_mu = mu * 1e10;
        source_scale *= Ctrl->S.src_mu;
    }
    if (Ctrl->D.active) {
        /* 时间函数与输出序列使用相同的物理采样间隔 */
        int time_function_nt;
        float *values = grt_get_time_function(&time_function_nt, (float)Ctrl->N.dt,
            Ctrl->D.tftype, Ctrl->D.tfparams);
        if (values == NULL) {
            GRTRaiseError("get time function error.\n");
        }
        time_function = grt_new_SACTRACE((float)Ctrl->N.dt, time_function_nt, 0.0f);
        memcpy(time_function->data, values, sizeof(*values) * time_function_nt);
        GRT_SAFE_FREE_PTR(values);
    }

    SACTRACE *prototype = new_lamb_trace(
        Ctrl->N.nt, Ctrl->N.dt, horizontal_distance, depsrc, deprcv,
        Ctrl->A.azimuth, begin_time, Ctrl->H.vp, Ctrl->H.vs, Ctrl->H.rho);
    set_lamb_arrivals(
        prototype, Ctrl->H.nu, horizontal_distance, depsrc, deprcv,
        distance, Ctrl->H.vs);
    SACTRACE *base[3] = {0};
    SACTRACE *derivative[3][3] = {{0}};
    allocate_lamb_traces(prototype, calc_upar, base, derivative);

    for (int source_index = 0; source_index < GRT_SRC_M_NUM; ++source_index) {
        if (!lamb_need_src(Ctrl->computeType, source_index, surface)) {
            continue;
        }
        SACTRACE *source_base[3] = {0};
        SACTRACE *source_derivative[3][3] = {{0}};
        allocate_lamb_traces(prototype, calc_upar, source_base, source_derivative);

        fill_source_traces(
            rot2ZNE, horizontal_distance, Ctrl->computeType, source_scale, Ctrl->H.nu,
            Ctrl->A.azrad, Ctrl->mchn, calc_upar, Ctrl->N.nt, &result, source_index,
            force_factor, moment_factor, force_derivative_factor, moment_derivative_factor,
            source_base, source_derivative);
        /* 先去除 Lamb 解自带的时间积分，再执行用户指定的时间操作 */
        postprocess_lamb_traces(
            calc_upar, Ctrl->I.int_times, Ctrl->J.dif_times, time_function,
            source_base, source_derivative);
        accumulate_lamb_traces(
            calc_upar, source_base, source_derivative, base, derivative);
        free_lamb_traces(calc_upar, source_base, source_derivative);
    }

    /* 填充阶段已经按照请求的坐标系生成结果 */
    for (int component = 0; component < 3; ++component) {
        save_to_sac(Ctrl->O.s_output_dir, "", chs[component], base[component]);
        if (calc_upar) {
            for (int direction = 0; direction < 3; ++direction) {
                char prefix[2] = {(char)tolower(chs[direction]), '\0'};
                save_to_sac(Ctrl->O.s_output_dir, prefix, chs[component], derivative[direction][component]);
            }
        }
    }
    free_lamb_traces(calc_upar, base, derivative);

    if (time_function != NULL) {
        char *path = NULL;
        GRT_SAFE_ASPRINTF(&path, "%s/sig.sac", Ctrl->O.s_output_dir);
        grt_write_SACTRACE(path, time_function);
        GRT_SAFE_FREE_PTR(path);
        grt_free_SACTRACE(time_function);
    }

    if (!Ctrl->s.active) {
        GRTRaiseInfo("Under \"%s\".", Ctrl->O.s_output_dir);
        GRTRaiseInfo("Synthetic Seismograms of %-13s source done.", srcTypeFullName[Ctrl->computeType]);
        if (Ctrl->D.active) {
            GRTRaiseInfo("Time Function saved.");
        }
    }

    GRT_SAFE_FREE_PTR(result.G);
    GRT_SAFE_FREE_PTR(result.dG_source);
    GRT_SAFE_FREE_PTR(result.dG_receiver);
    GRT_SAFE_FREE_PTR(result.dG_mixed);
    grt_free_SACTRACE(prototype);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
