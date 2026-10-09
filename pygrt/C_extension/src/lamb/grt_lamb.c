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
    } A;

    /** 源强缩放 */
    struct {
        bool active;
        bool mult_src_mu;  ///< 是否将单点源强度乘以震源处的剪切模量
        real_t scale;      ///< 源强
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
        char *option;  ///< 完整 -D 时间函数选项
    } D;

    /** 时间延迟 */
    struct {
        bool active;     ///< 是否设置时间延迟
        real_t delayT0;  ///< 参考时间偏移，s
        real_t delayV0;  ///< 参考速度，km/s
        bool refFirstP;  ///< 是否参考初至 P
    } E;

    /** 物理时间序列 */
    struct {
        bool active;
        int nt;
        real_t dt;
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

    /** 参数解析时确定的震相选择 */
    struct {
        bool active;
        unsigned int phase_mask;
    } L;

    /** 接收点文件 */
    struct {
        bool active;
        char *path;  ///< 接收点文件
    } Q;

    /** 有限震源断层 */
    struct {
        bool active;
        char *option;         ///< 有限震源选项，确定采样间隔后读入
        FINITE_FAULT *faults;  ///< 震源断层数组
        size_t nfault;         ///< 震源断层数量
        real_t dL;             ///< 走向剖分间隔，km
        real_t dW;             ///< 倾向剖分间隔，km
    } C;

    /** 有限接收断层 */
    struct {
        bool active;
        FINITE_FAULT *faults;  ///< 接收断层数组
        size_t nfault;         ///< 接收断层数量
        real_t dL;             ///< 走向剖分间隔，km
        real_t dW;             ///< 倾向剖分间隔，km
    } U;

    /** 子源并行线程数 */
    struct {
        bool active;
        int nthreads;  ///< 子源并行线程数
    } P;

    /** 单点源临时参数 */
    GRT_SYN_TYPE source_type;             ///< 单点源类型
    real_t mechanism[GRT_MECHANISM_NUM];  ///< 单点源机制参数

} GRT_MODULE_CTRL;


/** 每个接收点的时间窗及多源初至 */
typedef struct {
    real_t begin;       ///< 采样对齐后的共同起点，s
    real_t travtPS[2];  ///< 含破裂延迟的最早 P、S 初至，s，无有效源时为 INFINITY
} LAMB_RECEIVER_TIMING;

/** Lamb 求解器输出的各类数组 */
typedef struct {
    real_t (*G)[3][3];               ///< 位移 Green 函数
    real_t (*dG_source)[3][3][3];    ///< 源点一阶导数
    real_t (*dG_receiver)[3][3][3];  ///< 接收点一阶导数
    real_t (*dG_mixed)[3][3][3][3];  ///< 混合二阶导数
} LAMB_RESULT;

/** Lamb 模块使用的坐标变换矩阵 */
typedef struct {
    real_t local_from_global[3][3];
    real_t component_from_global[3][3];
} LAMB_COORDINATES;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    GRT_SAFE_FREE_PTR(Ctrl->D.option);
    GRT_SAFE_FREE_PTR(Ctrl->O.s_output_dir);
    GRT_SAFE_FREE_PTR(Ctrl->Q.path);
    GRT_SAFE_FREE_PTR(Ctrl->C.option);
    grt_finite_fault_free(Ctrl->C.nfault, Ctrl->C.faults);
    grt_finite_fault_free(Ctrl->U.nfault, Ctrl->U.faults);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印 Lamb 模块的命令行使用说明 */
static void print_help(void)
{
printf("\n"
"[grt lamb] %s\n\n", GRT_VERSION);printf(
"    Compute three-component dynamic displacement from Lamb closed-form solutions.\n"
"\n"
"    Define five groups of information:\n"
"      1. Medium and sampling: -H gives vp/vs/rho for the homogeneous halfspace;\n"
"         -N gives the minimum output sample count and time interval. No GF library is needed.\n"
"      2. Source location: a point source is at the horizontal origin; -Ds\n"
"         gives its depth. For finite sources, -C supplies all locations.\n"
"      3. Receiver locations: choose polar coordinates (-R, -A, -Dr), a point\n"
"         file (-Q), or receiver faults (-U). These three modes are exclusive.\n"
"      4. Source mechanism and strength: a point source requires -S; add at most\n"
"         one of -M (shear/tensile), -F (force), or -T (moment tensor).\n"
"         Without -M/-F/-T, use an explosion. -C supplies both location and\n"
"         mechanism/strength, so omit -Ds/-S/-M/-F/-T when using source faults.\n"
"      5. Output: -O sets the SAC directory; -n selects ZNE and -e adds spatial\n"
"         derivatives. Finite sources always output ZNE.\n"
"\n"
"    All coordinates share one horizontal origin; with -C, -R/-A locate the\n"
"    receiver relative to that origin. Rectangular sources require explicit\n"
"    subdivision sizes through +i on -C. Receiver subdivision with -U is optional.\n"
"    Optional time controls: -D sets a time function, -E sets the output start,\n"
"    -I/-J integrate or differentiate in time; -L selects phases and -P threads.\n"
"    The source and receiver depths select the appropriate Lamb solution.\n"
"    When both are on the free surface, only -F is supported and -e is ignored.\n"
"    Displacements are in cm, with Z upward, R radial outward and T clockwise\n"
"    from R by default; output filenames follow `syn`.\n"
"\n"
"\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    # Point source\n"
"    grt lamb -H<vp>/<vs>/<rho> -N<nt>/<dt> -Ds<depsrc> -S[u]<scale> -O<outdir>\n"
"             [-M<strike>/<dip>[/<rake>] | -T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>\n"
"              | -F<fn>/<fe>/<fz>] <receiver options> [common options]\n"
"\n"
"    # Finite source\n"
"    grt lamb -H<vp>/<vs>/<rho> -N<nt>/<dt> -C<fault>[+i<dL>/<dW>] -O<outdir>\n"
"             <receiver options> [common options]\n"
"\n"
"    Receiver options (choose one mode):\n"
"      -R<dist> -A<azimuth> -Dr<deprcv>\n"
"      -Q<points>\n"
"      -U<fault>+i<dL>/<dW>\n"
"\n"
"    Common options:\n"
"      [-D<tftype>[/<tfparams>][+d<delay>]] [-E[p]<t0>[/<v0>]]\n"
"      [-I<odr>] [-J<odr>] [-L<P,S,R,PP,SS,PS,SP,sPs>] [-P<nthreads>]\n"
"      [-n] [-e] [-s] [-h]\n"
"\n"
"\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -H<vp>/<vs>/<rho>\n"
"                  Homogeneous halfspace parameters. vp and vs are P- and S-wave\n"
"                  speeds in km/s, and rho is density in g/cm^3. The Poisson ratio\n"
"                  is calculated from vp/vs, which must be greater than sqrt(2).\n"
"\n"
"    -N<nt>/<dt>   Minimum sample count and physical time interval in seconds.\n"
"                  Long source processes (including delay) extend nt for all receivers.\n"
"                  Time functions use linear convolution.\n"
"                  The time series starts at the origin time.\n"
"\n"
"    -R<dist>       Horizontal source-receiver distance in km, positive.\n"
"                  With -C, measure from the common horizontal origin.\n"
"                  -R/-A/-Dr, -Q and -U are mutually exclusive receiver modes.\n"
"\n"
"    -Ds<depsrc>    Required point-source depth in km, nonnegative; omit with -C.\n"
"\n"
"    -Dr<deprcv>    Required polar-receiver depth in km, nonnegative; omit with -Q/-U.\n"
"\n"
"    -A<azimuth>    Azimuth from source to receiver, in degree, [0, 360].\n"
"                  With finite sources, this is the polar receiver azimuth\n"
"                  from the common horizontal origin. Omit with -Q/-U.\n"
"\n"
"    -S[u]<scale>   Source scale.  Moment sources use dyne-cm and force sources\n"
"                  use dyne.  `-Su` multiplies the scale by source shear modulus.\n"
"                  Use -C instead of -S/-M/-F/-T for finite sources.\n"
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
"                  A polar receiver writes here directly; -Q/-U always use\n"
"                  index_north_east_depth subdirectories. Finite-source sig.sac\n"
"                  saves the complete scalar moment rate (dyne-cm/s) only\n"
"                  when a time function is explicitly specified.\n"
"                  Receiver coordinates/angles use SAC unused1..5; receiver index\n"
"                  and fault grouping use unused10..14 (project C field names).\n"
"\n"
"    -D<tftype>[/<tfparams>][+d<delay>]\n"
"                  Convolve a time function. Source time functions use area\n"
"                  normalization.\n"
"                  There are several options:\n"
"                  + Impulse\n"
"                    set -D%c.\n", GRT_SIG_IMPULSE); printf(
"                  + Parabolic wave (y = a*x^2 + b*x)\n"
"                    set -D%c/<t0>, <t0> (secs) is the duration of wave.\n", GRT_SIG_PARABOLA); printf(
"                    e.g.\n"
"                         -D%c/1.3\n", GRT_SIG_PARABOLA); printf(
"                  + Trapezoidal wave\n"
"                    set -D%c/<t1>/<t2>/<t3>, rise/plateau/fall durations in seconds.\n", GRT_SIG_TRAPEZOID); printf(
"                    Durations must be nonnegative, with a positive total.\n"
"                    t2=0 gives a triangle; t1=t3=0 gives a rectangle.\n"
"                    e.g.\n"
"                         -D%c/0.1/0.1/0.2\n", GRT_SIG_TRAPEZOID); printf(
"                         -D%c/0.4/0/0.2 (become a triangle)\n", GRT_SIG_TRAPEZOID); printf(
"                         -D%c/0/0.5/0 (become a rectangle)\n", GRT_SIG_TRAPEZOID); printf(
"                  + AsymmetricCosine\n"
"                    set -D%c/<t1>/<t2>, positive rise/fall durations in seconds.\n", GRT_SIG_ASYMMETRIC_COSINE); printf(
"                    The peak is at t1 and the end is at t1+t2.\n"
"                    e.g. -D%c/4.5/1.5+d43\n", GRT_SIG_ASYMMETRIC_COSINE); printf(
"                  + Custom wave\n"
"                    set -D%c/<path>, <path> is the filepath to a custom\n", GRT_SIG_CUSTOM); printf(
"                    Time Function ASCII file. The file has just one column\n"
"                    of amplitude and no other columns. Its sequence sum should\n"
"                    be 1/dt, where dt is the sampling interval; the program\n"
"                    normalizes it with a warning when it is not.\n"
"                    The file can contain unlimited comment lines with prefix\n"
"                    \"#\".\n"
"                    e.g.\n"
"                         -D%c/tfunc.txt\n", GRT_SIG_CUSTOM); printf(
"                  Also accepts a signed Ricker convolution wavelet:\n"
"                  -D%c/<f0>, <f0> is the dominant frequency in Hz.\n", GRT_SIG_RICKER); printf(
"                  Its analytic peak amplitude is 1, without area normalization;\n"
"                  it is not a unit-slip source process.\n"
"                  To match the physical time interval, parameters of the time\n"
"                  function may be slightly modified. The corresponding time\n"
"                  function is saved as a SAC file under <outdir>.\n"
"\n"
"                  Append +d<delay> for rupture delay in seconds.\n"
"                  Append a complete -D option at the end of each fault row\n"
"                  to specify its rupture process.\n"
"\n"
"    -E[p]<t0>[/<v0>]\n"
"                  Introduce a time shift in the output SAC records. The time\n"
"                  series starts at <t0> + distance/<v0>, where distance is the\n"
"                  straight-line source-receiver distance. <v0> is a reference\n"
"                  velocity in km/s; when omitted or zero, no distance correction\n"
"                  is applied.\n"
"                  -Ep<t0> starts each receiver at <t0> plus its earliest P\n"
"                  arrival, including the sampled rupture delay. Use -Ep-10.\n"
"                  Without -E, the first time sample is the origin time.\n"
"                  Faults are subdivided before computing arrivals and starts.\n"
"                  Starts are rounded to the output sampling interval.\n"
"                  -E sets the solve end time; the internal start is fixed.\n"
"                  Convolution precedes cropping to the requested window.\n"
"                  Multiple sources save only the earliest P/S arrival picks.\n"
"\n"
"    -I<odr>        Apply odr time integrations after physical normalization.\n"
"\n"
"    -J<odr>        Apply odr time differentiations after physical normalization.\n"
"\n"
"    -n             Write receiver components as Z, N and E; default is Z, R and T.\n"
"                  Finite sources always output ZNE.\n"
"\n"
"    -e             Also write spatial derivatives with direction prefixes matching\n"
"                  the output coordinates.\n"
"\n"
"    -L<phases>     Keep only selected phase terms. Use a comma-separated list of\n"
"                  P, S, R, PP, SS, PS, SP and sPs.\n"
"                  The applicable phase names depend on whether the selected\n"
"                  Lamb problem is of the first, second or third kind.\n"
"                  If no valid phase remains, a warning is issued and the output is all zeros.\n"
"\n"
"    -s             Do not print completion information.\n"
"\n"
"    -C<fault>[+i<dL>/<dW>]\n"
"                  Coulomb source faults, including point-source Kode records.\n"
"                  The file supplies location, mechanism and signed slip/potency;\n"
"                  do not set -Ds/-S/-M/-F/-T. Finite sources always output ZNE.\n"
"                  Kode 100/200/300: rectangular shear/tensile sources.\n"
"                  Kode 400: point double couple; Kode 500: point tensile/inflation.\n"
"                  Two header lines precede the 11 numeric columns. An exact\n"
"                  \"rake\" token in the seventh header column selects rake/net\n"
"                  slip for Kode 100; the filename suffix does not select format.\n"
"                  +i gives along-strike/dip subdivision sizes (km).\n"
"                  Point-source Kode records remain single points at fault centers.\n"
"                  Rectangular sources require explicit +i subdivision sizes.\n"
"\n"
"    -Q<points>    ASCII rows: north east depth (km) [strike dip rake (degrees)].\n"
"                  Lines starting with # are comments. Optional angles are saved\n"
"                  for later stress projection and do not affect synthesis.\n"
"\n"
"    -U<fault>[+i<dL>/<dW>]\n"
"                  Coulomb receiver faults. Receivers are subfault centers;\n"
"                  no receiver-area averaging is performed.\n"
"                  Without +i, use each fault center. Slip magnitude is ignored.\n"
"                  Only Kode=100 is supported;\n"
"                  receiver angles and grouping are saved.\n"
"\n"
"    -P<nthreads>  OpenMP source-point threads. Receivers are processed serially.\n"
"\n"
"    -h             Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt lamb -H8.0/4.62/3.3 -N6000/0.001 -R10 -A30 -S1e24 -Ds5 -Dr0 -M100/20/80 -e -n -Ores\n"
"\n\n\n"
"    Finite sources with explicit receivers:\n"
"        grt lamb -H6/3.464/2.7 -N1000/0.01 -Cfaults.inp+i1/1 -Qreceivers.txt -Osyn_ff\n"

);
}

/**
 * 解析 Lamb 模块的命令行选项
 *
 * @param[in,out]  Ctrl  Lamb 模块参数控制结构体
 * @param[in]      argc  命令行参数个数
 * @param[in]      argv  命令行参数数组
 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    const char *phase_list = NULL;

    GRTCheckOptionSet(argc > 1);
    Ctrl->source_type = GRT_SYN_EX;

    int opt;
    char extra;
    while ((opt = getopt(argc, argv, ":H:N:R:A:S:M:F:T:O:D:E:I:J:C:U:Q:L:P:nesh")) != -1) {
        switch (opt) {
            /* 半空间参数 */
            case 'H': {
                char extra;
                real_t vp, vs, rho;
                real_t vp_vs_squared;
                if (sscanf(optarg, "%lf/%lf/%lf%c", &vp, &vs, &rho, &extra) != 3) {
                    GRTBadOptionError(H, "expected vp/vs/rho.");
                }
                if (vp <= 0.0 || vs <= 0.0 || rho <= 0.0) {
                    GRTBadOptionError(H, "vp, vs and rho should be positive.");
                }
                vp_vs_squared = GRT_SQUARE(vp / vs);
                Ctrl->H.nu = (vp_vs_squared - 2.0) / (2.0 * (vp_vs_squared - 1.0));
                if (Ctrl->H.nu <= 0.0 || Ctrl->H.nu >= 0.5) {
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
                if (sscanf(optarg, "%d/%lf%c", &nt, &dt, &extra) != 2 || nt <= 0 || dt <= 0.0) {
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
                if (sscanf(optarg, "%lf%c", &Ctrl->R.dist, &extra) != 1 || Ctrl->R.dist <= 0.0) {
                    GRTBadOptionError(R, "horizontal distance should be positive.");
                }
                Ctrl->R.active = true;
                break;
            }

            /* 方位角 */
            case 'A': {
                char extra;
                if (sscanf(optarg, "%lf%c", &Ctrl->A.azimuth, &extra) != 1 || Ctrl->A.azimuth < 0.0 || Ctrl->A.azimuth > 360.0) {
                    GRTBadOptionError(A, "azimuth should be in [0, 360].");
                }
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
                if (sscanf(scale, "%lf%c", &Ctrl->S.scale, &extra) != 1) {
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
                if ((count != 2 && count != 3) || (count == 2 && sscanf(optarg, "%lf/%lf%c", &strike, &dip, &extra) != 2)) {
                    GRTBadOptionError(M, "expected strike/dip[/rake].");
                }
                if (strike < 0.0 || strike > 360.0 || dip < 0.0 || dip > 90.0 || (count == 3 && (rake < -180.0 || rake > 180.0))) {
                    GRTBadOptionError(M, "strike, dip or rake is out of bound.");
                }
                Ctrl->mechanism[0] = strike;
                Ctrl->mechanism[1] = dip;
                Ctrl->mechanism[2] = count == 3 ? rake : 0.0;
                Ctrl->source_type = count == 3 ? GRT_SYN_DC : GRT_SYN_TS;
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
                    &Ctrl->mechanism[0], &Ctrl->mechanism[1], &Ctrl->mechanism[2], &extra);
                if (count != 3) {
                    GRTBadOptionError(F, "expected fn/fe/fz.");
                }
                Ctrl->source_type = GRT_SYN_SF;
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
                    &Ctrl->mechanism[0], &Ctrl->mechanism[1], &Ctrl->mechanism[2],
                    &Ctrl->mechanism[3], &Ctrl->mechanism[4], &Ctrl->mechanism[5], &extra);
                if (count != 6) {
                    GRTBadOptionError(T, "expected six moment-tensor components.");
                }
                Ctrl->source_type = GRT_SYN_MT;
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
                    if (sscanf(optarg + 1, "%lf%c", &Ctrl->Depth.depsrc, &extra) != 1 || Ctrl->Depth.depsrc < 0.0) {
                        GRTBadOptionError(Ds, "source depth should be nonnegative.");
                    }
                    Ctrl->Depth.s_active = true;
                } else if (optarg[0] == 'r' && optarg[1] != '/') {
                    char extra;
                    if (Ctrl->Depth.r_active) {
                        GRTBadOptionError(Dr, "the option is duplicated.");
                    }
                    if (sscanf(optarg + 1, "%lf%c", &Ctrl->Depth.deprcv, &extra) != 1 || Ctrl->Depth.deprcv < 0.0) {
                        GRTBadOptionError(Dr, "receiver depth should be nonnegative.");
                    }
                    Ctrl->Depth.r_active = true;
                } else {
                    GRT_SAFE_FREE_PTR(Ctrl->D.option);
                    Ctrl->D.active = true;
                    GRT_SAFE_ASPRINTF(&Ctrl->D.option, "-D%s", optarg);
                }
                break;

            /* 时间延迟 */
            case 'E': {
                Ctrl->E.active = true;
                if (optarg[0] == 'p') {
                    char extra;
                    real_t delayT0;
                    if (sscanf(optarg + 1, "%lf%c", &delayT0, &extra) != 1) {
                        GRTBadOptionError(E, "expected t0 after -Ep.");
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
                    if ((count != 1 && count != 2) || (count == 1 && sscanf(optarg, "%lf%c", &delayT0, &extra) != 1)) {
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

            /* 选择输出的震相 */
            case 'L':
                Ctrl->L.active = true;
                phase_list = optarg;
                break;

            /* 是否静默输出 */
            case 's':
                Ctrl->s.active = true;
                break;

            // 任意接收点文件，包含坐标和可选接收机制
            case 'Q':
                Ctrl->Q.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->Q.path);
                Ctrl->Q.path = strdup(optarg);
                break;

            // 有限震源文件，时间函数在确定采样间隔后统一读取
            case 'C':
                Ctrl->C.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->C.option);
                Ctrl->C.option = strdup(optarg);
                break;

            // 有限接收断层文件，解析几何和可选剖分尺寸
            case 'U':
                Ctrl->U.active = true;
                grt_finite_fault_free(Ctrl->U.nfault, Ctrl->U.faults);
                Ctrl->U.faults = grt_finite_fault_from_option(optarg, &Ctrl->U.nfault, &Ctrl->U.dL, &Ctrl->U.dW, false, 0, 1, NULL);
                break;

            // 设置子源合成的线程数
            case 'P':
                Ctrl->P.active = true;
                if(sscanf(optarg, "%d%c", &Ctrl->P.nthreads, &extra) != 1 || Ctrl->P.nthreads <= 0) {
                    GRTBadOptionError(P, "Expected a positive thread count.");
                }
                grt_set_num_threads(Ctrl->P.nthreads);
                break;

            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 输出目录、半空间介质和时间采样参数均为必选项
    GRTCheckOptionActive(Ctrl, O);
    GRTCheckOptionActive(Ctrl, H);
    GRTCheckOptionActive(Ctrl, N);

    // 有限震源文件不能与点源强度、机制或源深度选项同时使用
    if(Ctrl->C.active && (Ctrl->S.active || (Ctrl->M.active + Ctrl->F.active + Ctrl->T.active) || Ctrl->Depth.s_active)) {
        GRTRaiseError("Finite sources and point-source parameters are mutually exclusive.");
    }

    // 点源必须显式指定震源强度
    if (!Ctrl->C.active && !Ctrl->S.active) {
        GRTRaiseError("Point sources require -S.");
    }

    // 点源机制最多只能选用 -M、-F、-T 中的一种
    if((Ctrl->M.active + Ctrl->F.active + Ctrl->T.active) > 1) {
        GRTRaiseError("Only one point-source mechanism may be selected.");
    }

    // 逐点文件、有限接收断层和极坐标接收点只能选择一种
    int explicit_rcv = Ctrl->Q.active + Ctrl->U.active;

    // -Q 和 -U 互斥，且不能与极坐标选项 -R/-A 混用
    if(explicit_rcv > 1 || (explicit_rcv && (Ctrl->R.active || Ctrl->A.active))) {
        GRTRaiseError("Receiver geometry options are mutually exclusive.");
    }

    // 接收文件已经包含深度，不能再用 -Dr 指定
    if((Ctrl->Q.active || Ctrl->U.active) && Ctrl->Depth.r_active) {
        GRTRaiseError("-Q/-U provide receiver depths; do not set -Dr.");
    }

    // 使用极坐标接收点时必须给出方位角
    if(!explicit_rcv && !Ctrl->A.active) {
        GRTRaiseError("A single polar receiver requires -A.");
    }

    // 点源必须给出源深度，有限震源的深度来自断层文件
    if (!Ctrl->C.active && !Ctrl->Depth.s_active) {
        GRTRaiseError("Lamb point source requires -Ds.");
    }

    // 极坐标接收点必须给出接收深度
    if (!Ctrl->Q.active && !Ctrl->U.active && !Ctrl->Depth.r_active) {
        GRTRaiseError("Lamb polar receiver requires -Dr.");
    }

    // Lamb 极坐标接收点必须显式给出震中距
    if(!explicit_rcv && !Ctrl->R.active) {
        GRTRaiseError("Lamb polar receiver requires -R.");
    }

    // 有限震源的各子源方位不同，叠加时统一使用 ZNE 坐标系
    if(Ctrl->C.active) {
        Ctrl->n.active = true;
    }
    Ctrl->L.phase_mask = grt_lamb_parse_phase_list(phase_list, GRT_LAMB_ALL_PHASES);
}


/**
 * 判断源点和接收点是否同时位于自由表面
 *
 * @param[in]  source_depth    源点深度，km
 * @param[in]  receiver_depth  接收点深度，km
 *
 * @return 同时位于 z=0 时返回 true，否则返回 false
 */
static bool is_surface_source_receiver(const real_t source_depth, const real_t receiver_depth)
{
    /* 只有源点和接收点同时位于 z=0 时才使用第一类 Lamb 解 */
    return source_depth == 0.0 && receiver_depth == 0.0;
}


/**
 * 创建 Lamb 求解器与输出所需的坐标变换矩阵
 *
 * @param[in]   azrad        源点到接收点的方位角，弧度
 * @param[in]   rot2ZNE      是否输出 Z、N、E 坐标
 * @param[out]  coordinates  坐标变换矩阵
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
 * @param[in]  G                      位移 Green 函数
 * @param[in]  component_from_global  输出分量坐标变换矩阵
 * @param[in]  output_component       输出分量索引
 * @param[in]  source                 全局坐标下的源力向量
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
 * @param[in]  dG                     源点一阶导数 Green 函数
 * @param[in]  component_from_global  输出分量坐标变换矩阵
 * @param[in]  output_component       输出分量索引
 * @param[in]  source                 全局坐标下的源矩张量
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
 * @param[in]  dG                     接收点一阶导数 Green 函数
 * @param[in]  component_from_global  坐标变换矩阵
 * @param[in]  derivative_direction   导数方向索引
 * @param[in]  output_component       输出分量索引
 * @param[in]  source                 全局坐标下的源力向量
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
 * @param[in]  dG                     混合二阶导数 Green 函数
 * @param[in]  component_from_global  坐标变换矩阵
 * @param[in]  derivative_direction   接收点导数方向索引
 * @param[in]  output_component       输出分量索引
 * @param[in]  source                 全局坐标下的源矩张量
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
 * @param[in]   source_type  当前计算的震源类型
 * @param[in]   scale        源强缩放因子
 * @param[in]   nu           半空间泊松比
 * @param[in]   azrad        源点到接收点的方位角，弧度
 * @param[in]   mchn         震源机制参数数组
 * @param[in]   par_theta    是否计算方位角导数
 * @param[in]   coef         当前计算项的距离缩放因子
 * @param[out]  srcRadi      各源项和各输出分量的辐射系数
 */
static void make_source_radiation(
    const GRT_SYN_TYPE source_type, const real_t scale,
    const real_t nu, const real_t azrad, const real_t mchn[GRT_MECHANISM_NUM],
    const bool par_theta, const real_t coef, realChnlGrid srcRadi)
{
    /* 先清零所有源项，再根据源类型生成当前距离阶数的辐射系数 */
    memset(srcRadi, 0, sizeof(realChnlGrid));
    /* grt_set_source_radiation 使用 vp/vs 作为水平分量和垂直分量的换算比 */
    grt_set_source_radiation(
        srcRadi, source_type, par_theta, scale, coef,
        sqrt(2.0 * (1.0 - nu) / (1.0 - 2.0 * nu)), azrad, mchn);
}


/**
 * 计算 Lamb 位移 Green 函数及按需计算其空间导数
 *
 * @param[in]   nu                     半空间泊松比
 * @param[in]   tbar                   无量纲时间序列
 * @param[in]   nt                     时间序列长度
 * @param[in]   horizontal_distance    源点与接收点的水平距离，km
 * @param[in]   source_depth           源点深度，km
 * @param[in]   receiver_depth         接收点深度，km
 * @param[in]   azimuth_degree         源点到接收点的方位角，度
 * @param[in]   source_type            当前计算的震源类型
 * @param[in]   calculate_derivatives  是否计算空间导数
 * @param[in]   phase_mask             震相筛选掩码
 * @param[out]  result                 Lamb Green 函数及其导数的结果结构体
 */
static void make_lamb_result(
    const real_t nu, const real_t *tbar, const int nt, const real_t horizontal_distance,
    const real_t source_depth, const real_t receiver_depth, const real_t azimuth_degree,
    const GRT_SYN_TYPE source_type, const bool calculate_derivatives, const unsigned int phase_mask,
    LAMB_RESULT *result)
{
    /* 只有需要空间导数的源才分配相应的导数数组 */
    const bool surface = is_surface_source_receiver(source_depth, receiver_depth);
    const bool moment = source_type != GRT_SYN_SF;
    const bool force = source_type == GRT_SYN_SF;
    const bool need_upar = calculate_derivatives && !surface;
    const size_t nt_size = (size_t)nt;
    result->G = GRT_SAFE_CALLOC(nt_size, sizeof(*result->G));
    result->dG_source = moment ? GRT_SAFE_CALLOC(nt_size, sizeof(*result->dG_source)) : NULL;
    result->dG_receiver = need_upar && force ? GRT_SAFE_CALLOC(nt_size, sizeof(*result->dG_receiver)) : NULL;
    result->dG_mixed = need_upar && moment ? GRT_SAFE_CALLOC(nt_size, sizeof(*result->dG_mixed)) : NULL;

    /* 地表、单侧地下和双侧地下分别对应三类 Lamb 求解器 */
    if (surface) {
        grt_solve_lamb1(nu, tbar, nt, azimuth_degree, 0.0, phase_mask, result->G);
    } else if (source_depth > 0.0 && receiver_depth > 0.0) {
        grt_solve_lamb3(nu, tbar, nt, horizontal_distance,
            source_depth, receiver_depth, azimuth_degree,
            phase_mask,
            result->G, result->dG_source, result->dG_receiver, result->dG_mixed);
    } else {
        grt_solve_lamb2(nu, tbar, nt, horizontal_distance,
            source_depth, receiver_depth, azimuth_degree,
            phase_mask,
            result->G, result->dG_source, result->dG_receiver, result->dG_mixed);
    }
}


enum {
    LAMB_PHASE_P,                         ///< t0/kt0：直达 P 波
    LAMB_PHASE_S,                         ///< t1/kt1：直达 S 波
    LAMB_PHASE_R,                         ///< t2/kt2：Rayleigh 波参考到时
    LAMB_PHASE_PP,                        ///< t3/kt3：反射 PP 波
    LAMB_PHASE_SS,                        ///< t4/kt4：反射 SS 波
    LAMB_PHASE_PS,                        ///< t5/kt5：PS 转换波
    LAMB_PHASE_SP,                        ///< t6/kt6：SP 转换波
    LAMB_PHASE_sPs,                       ///< t7/kt7：滑行 sPs 波
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
    SACTRACE *sac, const int phase, const real_t tbar, const real_t time_scale, const real_t delay, const char *name)
{
    if (phase < 0 || phase >= LAMB_PHASE_COUNT || tbar < 0.0) {
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
    *times[phase] = (float)(tbar * time_scale + delay);
    copy_lamb_phase_name(names[phase], name);
}


/** 根据源点和接收点位置写入 Lamb SAC 记录的震相到时 */
static void set_lamb_arrivals(
    SACTRACE *sac, const real_t nu, const real_t horizontal_distance,
    const real_t source_depth, const real_t receiver_depth,
    const real_t direct_distance, const real_t vs, const real_t delay)
{
    clear_lamb_arrivals(sac);
    const real_t time_scale = direct_distance / vs;
    real_t tP;
    real_t tR;
    set_lamb_arrival(sac, LAMB_PHASE_S, 1.0, time_scale, delay, "S");
    grt_compute_lamb1_travt(nu, &tP, &tR);

    if (source_depth > 0.0 && receiver_depth > 0.0) {
        const real_t reflected_distance = hypot(horizontal_distance, source_depth + receiver_depth);
        tR *= reflected_distance / direct_distance;
    }
    set_lamb_arrival(sac, LAMB_PHASE_R, tR, time_scale, delay, "R");

    if (is_surface_source_receiver(source_depth, receiver_depth)) {
        set_lamb_arrival(sac, LAMB_PHASE_P, tP, time_scale, delay, "P");
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
        set_lamb_arrival(sac, LAMB_PHASE_P, tP, time_scale, delay, "P");
        set_lamb_arrival(sac, LAMB_PHASE_PP, tPP, time_scale, delay, "PP");
        set_lamb_arrival(sac, LAMB_PHASE_SS, tSS, time_scale, delay, "SS");
        set_lamb_arrival(sac, LAMB_PHASE_PS, tPS, time_scale, delay, "PS");
        set_lamb_arrival(sac, LAMB_PHASE_SP, tSP, time_scale, delay, "SP");
        set_lamb_arrival(sac, LAMB_PHASE_sPs, t_sPs, time_scale, delay, "sPs");
        return;
    }

    real_t t_sliding;
    grt_compute_lamb2_travt(
        nu, horizontal_distance, source_depth, receiver_depth, &tP, &t_sliding);
    set_lamb_arrival(sac, LAMB_PHASE_P, tP, time_scale, delay, "P");
    /* 地下源的滑行项沿互易路径为 SP，互易回地表源问题后为 PS */
    if (source_depth > 0.0) {
        set_lamb_arrival(sac, LAMB_PHASE_SP, t_sliding, time_scale, delay, "SP");
    } else {
        set_lamb_arrival(sac, LAMB_PHASE_PS, t_sliding, time_scale, delay, "PS");
    }
}


/**
 * 将单力源辐射系数转换为源力向量
 *
 * @param[in]   srcRadi       各源项和各分量的辐射系数
 * @param[in]   source_index  震源在源项数组中的索引
 * @param[out]  source        R、T、Z_down 坐标下的源力向量
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
 * @param[in]   srcRadi       各源项和各分量的辐射系数
 * @param[in]   source_index  震源在源项数组中的索引
 * @param[out]  source        R、T、Z_down 坐标下的源矩张量
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
 * @param[in]   local_from_global  全局坐标到局部坐标的变换矩阵
 * @param[in]   source_local       局部坐标下的源力
 * @param[out]  source_global      全局坐标下的源力
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
 * @param[in]   local_from_global  全局坐标到局部坐标的变换矩阵
 * @param[in]   source_local       局部坐标下的源矩
 * @param[out]  source_global      全局坐标下的源矩
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
 * @param[in]   srcRadi            各源项和各分量的辐射系数
 * @param[in]   source_index       震源在源项数组中的索引
 * @param[in]   force              当前源项是否为单力源
 * @param[in]   local_from_global  全局坐标到局部坐标的变换矩阵
 * @param[out]  source_force       全局坐标下的源力，单力源时写入
 * @param[out]  source_moment      全局坐标下的源矩，矩源时写入
 */
static void make_source_terms(
    const realChnlGrid srcRadi, const int source_index, const bool force,
    const real_t local_from_global[3][3], real_t source_force[3], real_t source_moment[3][3])
{
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


/** 将一个震源项的基本源型累加到共同响应
 *
 * @param[in]      rot2ZNE                   是否输出 Z、N、E 坐标
 * @param[in]      horizontal_distance       源点与接收点的水平距离，km
 * @param[in]      term                      震源项
 * @param[in]      scale                     当前子源的源强
 * @param[in]      nu                        半空间泊松比
 * @param[in]      azrad                     源点到接收点的方位角，弧度
 * @param[in]      calc_upar                 是否计算空间导数
 * @param[in]      nt                        时间序列长度
 * @param[in]      result                    Lamb Green 函数及其导数
 * @param[in]      factor                    介质物理归一化因子
 * @param[in]      distance                  源台直线距离，km
 * @param[in,out]  traces                    输出累加数组，先位移、再按导数方向排列
 */
static void accumulate_lamb_term(
    const bool rot2ZNE, const real_t horizontal_distance, const FINITE_SOURCE_TERM *term,
    const real_t scale, const real_t nu, const real_t azrad, const bool calc_upar,
    const int nt, const LAMB_RESULT *result, const real_t factor, const real_t distance, real_t **traces)
{
    LAMB_COORDINATES coordinates;
    realChnlGrid baseRadiation, receiverRadiation, thetaRadiation;
    make_lamb_coordinates(azrad, rot2ZNE, &coordinates);
    make_source_radiation(term->type, scale, nu, azrad, term->mechanism, false, 1.0, baseRadiation);
    if(calc_upar) {
        make_source_radiation(term->type, scale, nu, azrad, term->mechanism, false, 1e-5, receiverRadiation);
        if(!rot2ZNE) {
            make_source_radiation(term->type, scale, nu, azrad, term->mechanism, true, 1e-5 / horizontal_distance, thetaRadiation);
        }
    }

    for(int source_index = 0; source_index < GRT_SRC_M_NUM; ++source_index) {
        if(!grt_source_has_component(term->type, source_index)) continue;
        const bool force = GRT_SRC_M_INDEX_IS_FORCE(source_index);
        const real_t base_factor = force ? factor / distance : factor / (distance * distance);
        const real_t derivative_factor = base_factor / distance;
        const bool skip_transverse = !rot2ZNE && GRT_SRC_M_ORDERS[source_index] == 0;
        real_t source_force_global[3];
        real_t source_moment_global[3][3];
        make_source_terms(baseRadiation, source_index, force, coordinates.local_from_global,
            source_force_global, source_moment_global);

        for (int output_component = 0; output_component < 3; ++output_component) {
            if (skip_transverse && output_component == 2) continue;
            for (int n = 0; n < nt; ++n) {
                const real_t value = force
                    ? evaluate_force(result->G[n], coordinates.component_from_global, output_component, source_force_global)
                    : evaluate_moment(result->dG_source[n], coordinates.component_from_global, output_component, source_moment_global);
                traces[output_component][n] += base_factor * value;
            }
        }

        if (!calc_upar) {
            continue;
        }

        make_source_terms(receiverRadiation, source_index, force, coordinates.local_from_global,
            source_force_global, source_moment_global);

        real_t theta_force_global[3] = {0.0};
        real_t theta_moment_global[3][3] = {{0.0}};
        if (!rot2ZNE) {
            make_source_terms(thetaRadiation, source_index, force, coordinates.local_from_global,
                theta_force_global, theta_moment_global);
        }

        for (int direction = 0; direction < 3; ++direction) {
            for (int output_component = 0; output_component < 3; ++output_component) {
                if (skip_transverse && output_component == 2) continue;
                for (int n = 0; n < nt; ++n) {
                    real_t value;
                    real_t component_factor;
                    if (!rot2ZNE && direction == 2) {
                        value = force
                            ? evaluate_force(result->G[n], coordinates.component_from_global, output_component, theta_force_global)
                            : evaluate_moment(result->dG_source[n], coordinates.component_from_global, output_component, theta_moment_global);
                        component_factor = base_factor;
                    } else {
                        value = force
                            ? evaluate_force_derivative(result->dG_receiver[n], coordinates.component_from_global,
                                direction, output_component, source_force_global)
                            : evaluate_moment_derivative(result->dG_mixed[n], coordinates.component_from_global,
                                direction, output_component, source_moment_global);
                        component_factor = derivative_factor;
                    }
                    traces[3 + 3 * direction + output_component][n] += component_factor * value;
                }
            }
        }
    }
}

/**
 * 计算一对源台的物理 Lamb 脉冲响应
 *
 * @param[in]      nt          样本数
 * @param[in]      dt          采样间隔，s
 * @param[in]      begin       实际采样起点，s
 * @param[in]      medium      vp、vs、rho
 * @param[in]      nu          半空间泊松比
 * @param[in]      zs          源深度，km
 * @param[in]      zr          接收深度，km
 * @param[in]      dist        震中距，km
 * @param[in]      az          方位角，弧度
 * @param[in]      fault       震源断层
 * @param[in]      area        子断层面积，点源传 1
 * @param[in]      zne         是否输出 ZNE
 * @param[in]      upar        是否计算空间导数
 * @param[in]      phase_mask  震相筛选掩码
 * @param[in,out]  traces      输出累加数组，先位移、再按导数方向排列
 */
static void lamb_synthesis_samples(int nt, real_t dt, real_t begin, const real_t medium[3], real_t nu,
                                   real_t zs, real_t zr, real_t dist, real_t az, const FINITE_FAULT *fault,
                                   real_t area, bool zne, bool upar, unsigned int phase_mask, real_t **traces)
{
    real_t distance = hypot(dist, zs - zr);
    real_t *tbar = GRT_SAFE_MALLOC(nt * sizeof(*tbar));
    long long first_sample = llround(begin / dt);
    for(int i = 0; i < nt; ++i) {
        tbar[i] = (first_sample + i) * dt * medium[1] / distance;
    }

    // 同一断层的源项共用 Green 函数，有限断层的各源项均为矩源
    LAMB_RESULT result = {0};
    make_lamb_result(nu, tbar, nt, dist, zs, zr, az / DEG1, fault->terms[0].type, upar, phase_mask, &result);
    real_t factor = 1.0 / (PI * PI * medium[1] * medium[1] * medium[2]);
    real_t src_mu = medium[1] * medium[1] * medium[2] * 1e10;
    for(int t = 0; t < fault->nterms; ++t) {
        const FINITE_SOURCE_TERM *term = &fault->terms[t];
        real_t scale = term->scale * area * (term->with_mu ? src_mu : 1);
        accumulate_lamb_term(zne, dist, term, scale, nu, az, upar, nt, &result,
                             factor, distance, traces);
    }

    // 闭合解和混合导数都带有一次时间积分，各源项累加后统一恢复为物理解
    int nc = upar ? 12 : 3;
    for(int c = 0; c < nc; ++c) {
        grt_differential(traces[c], nt, dt);
    }
    GRT_SAFE_FREE_PTR(result.G);
    GRT_SAFE_FREE_PTR(result.dG_source);
    GRT_SAFE_FREE_PTR(result.dG_receiver);
    GRT_SAFE_FREE_PTR(result.dG_mixed);
    GRT_SAFE_FREE_PTR(tbar);
}

/**
 * 检查震源剖分尺寸并展开为统一源点集合
 *
 * @param[in,out]  Ctrl    命令行参数，保存断层
 * @param[in]      modarr  均匀半空间模型矩阵
 * @param[out]     nsrc    展开后的震源点数
 * @return         源点集合
 */
static SRC_POINT *build_lamb_sources(GRT_MODULE_CTRL *Ctrl, const real_t (*modarr)[GRT_MODARR_NCOL], size_t *nsrc)
{
    // Lamb 没有格林函数库提供默认间隔，矩形震源必须显式指定剖分尺寸
    for(size_t i = 0; i < Ctrl->C.nfault; ++i) {
        if(KODE_IS_FINITE(Ctrl->C.faults[i].kode) && Ctrl->C.dL <= 0) {
            GRTRaiseError("Lamb finite rectangular sources require +idL/dW.");
        }
    }

    // 在均匀半空间中逐条剖分震源，点源保持单个子源
    for(size_t i = 0; i < Ctrl->C.nfault; ++i) {
        FINITE_FAULT *fault = &Ctrl->C.faults[i];
        grt_finite_fault_subdiv(fault, KODE_IS_POINT(fault->kode) ? 0 : Ctrl->C.dL, KODE_IS_POINT(fault->kode) ? 0 : Ctrl->C.dW, 1, modarr);
    }
    return grt_src_points_from_faults(Ctrl->C.nfault, Ctrl->C.faults, nsrc);
}

/**
 * 将单点、逐点文件或有限接收断层统一为接收点集合
 *
 * @param[in,out]  Ctrl  命令行参数
 * @param[out]     nrcv  展开后的接收点数
 * @return         接收点集合
 */
static RCV_POINT *build_lamb_receivers(GRT_MODULE_CTRL *Ctrl, size_t *nrcv)
{
    // 逐点接收文件直接提供坐标和可选机制
    if(Ctrl->Q.active) {
        return grt_rcv_points_from_file(Ctrl->Q.path, nrcv);
    }

    // 接收断层在读入选项时已统一完成几何剖分
    if(Ctrl->U.active) {
        return grt_rcv_points_from_faults(Ctrl->U.nfault, Ctrl->U.faults, nrcv);
    }

    // 极坐标接收点转换为统一的北向、东向和深度坐标
    *nrcv = 1;
    return grt_rcv_points_from_polar(Ctrl->R.dist, Ctrl->A.azimuth, Ctrl->Depth.deprcv);
}

/**
 * 检查全部源台几何的 Lamb 求解支持，先计算各接收点初至，再确定共同起点
 *
 * @param[in,out]  Ctrl  命令行参数，地表单力模式可能关闭空间导数
 * @param[in]      nsrc  震源点数
 * @param[in]      srcs  源点集合
 * @param[in]      nrcv  接收点数
 * @param[in]      rcvs  接收点集合
 * @return         各接收点的时间窗及初至数组，调用方负责释放
 */
static LAMB_RECEIVER_TIMING *prepare_lamb_geometry_and_timing(GRT_MODULE_CTRL *Ctrl, size_t nsrc, const SRC_POINT *srcs, size_t nrcv, const RCV_POINT *rcvs)
{
    // 地表点力且所有接收点均在地表时使用第一类 Lamb 解，关闭不支持的空间导数
    if(!Ctrl->C.active && srcs[0].fault->terms[0].type == GRT_SYN_SF && srcs[0].depth == 0) {
        bool all_surface = true;
        for(size_t i = 0; i < nrcv; ++i) {
            all_surface &= rcvs[i].depth == 0;
        }
        if(all_surface && Ctrl->e.active) {
            GRTRaiseWarning("Surface-force Lamb derivatives are not supported; -e is ignored.");
            Ctrl->e.active = false;
        }
    }

    // 震源读入和剖分已完成，按实际几何分别计算 P、S 最早到时
    LAMB_RECEIVER_TIMING *timings = GRT_SAFE_CALLOC(nrcv, sizeof(*timings));
    for(size_t ir = 0; ir < nrcv; ++ir) {
        LAMB_RECEIVER_TIMING *timing = &timings[ir];
        timing->travtPS[0] = timing->travtPS[1] = INFINITY;
        real_t min_distance = INFINITY;
        for(size_t is = 0; is < nsrc; ++is) {
            real_t dist = hypot(rcvs[ir].north - srcs[is].north, rcvs[ir].east - srcs[is].east);
            bool surface = srcs[is].depth == 0 && rcvs[ir].depth == 0;

            // 所有源台组合必须有正的水平距离，以避开 Lamb 解的轴上奇点
            if(dist <= 0) {
                GRTRaiseError("Lamb horizontal distance must be positive at source %zu receiver %zu.", is, ir);
            }

            // 源台均在地表时，只支持通过 -F 指定的单力源
            if(surface && (Ctrl->C.active || srcs[is].fault->terms[0].type != GRT_SYN_SF)) {
                GRTRaiseError("When both source and receiver are on the free surface, only the single force source specified by -F is supported.");
            }

            // 含地表源台组合时，不允许计算尚未支持的空间导数
            if(surface && Ctrl->e.active) {
                GRTRaiseError("Lamb derivatives are unavailable for a surface source/receiver pair.");
            }

            // 只有非零源项参与最早初至和参考距离的计算，初至包含破裂延迟
            if(srcs[is].fault->nterms) {
                real_t distance = hypot(dist, srcs[is].depth - rcvs[ir].depth);
                real_t delay = srcs[is].fault->stf_delay;
                timing->travtPS[0] = fmin(timing->travtPS[0], distance / Ctrl->H.vp + delay);
                timing->travtPS[1] = fmin(timing->travtPS[1], distance / Ctrl->H.vs + delay);
                min_distance = fmin(min_distance, distance);
            }
        }

        // -Ep 参考含破裂延迟的最早 P，普通 -E 的参考距离仍取非零源的最短直线距离
        real_t begin = Ctrl->E.delayT0;
        if(Ctrl->E.refFirstP && timing->travtPS[0] != INFINITY) {
            begin += timing->travtPS[0];
        } else if(Ctrl->E.delayV0 > 0 && min_distance != INFINITY) {
            begin += min_distance / Ctrl->E.delayV0;
        }
        timing->begin = grt_sample_aligned_time(begin, Ctrl->N.dt);
    }
    return timings;
}

/**
 * 在接收点共同时间窗内合成并卷积一个子源
 *
 * @param[in]      Ctrl      命令行参数
 * @param[in]      source    当前源点
 * @param[in]      receiver  当前接收点
 * @param[in]      begin     当前接收点共同起点，s
 * @param[in,out]  data      当前线程的接收波形，按分量连续排列
 */
static void lamb_one_pair(const GRT_MODULE_CTRL *Ctrl, const SRC_POINT *source, const RCV_POINT *receiver,
                          real_t begin, real_t *data)
{
    int nt = Ctrl->N.nt, nc = Ctrl->e.active ? 12 : 3;
    real_t dt = Ctrl->N.dt;

    // 按实际源台坐标计算距离和方位角，各子源独立求解后统一叠加
    real_t north = receiver->north - source->north, east = receiver->east - source->east;
    real_t dist = hypot(north, east), az = atan2(east, north);
    if(az < 0) {
        az += 2 * PI;
    }

    // 矩形子源按面积换算强度，点源保持原强度
    real_t area = KODE_IS_FINITE(source->fault->kode) ? source->fault->width[source->isub] * source->fault->length[source->isub] : 1;

    const real_t *stf = source->fault->stfd;
    int stf_npts = source->fault->stf_npts;

    // 输出窗口只决定求解末时刻，内部起点始终固定，保留所有震相及卷积历史
    // 最多两次求导恢复混合 Green 函数，再一次求导恢复物理解，两端各需三个额外采样点
    long long first_output = llround((begin - source->fault->stf_delay) / dt);
    long long last_output = first_output + nt - 1;
    if(last_output < -3) {
        return;
    }
    int margin = 3;
    int tail = margin + Ctrl->J.dif_times;
    int work_nt = GRT_MAX(margin + tail + 1, last_output + margin + tail + 1);
    real_t b = -margin * dt;
    real_t *response = GRT_SAFE_CALLOC((size_t)nc * work_nt, sizeof(*response));
    real_t *convolution = GRT_SAFE_MALLOC(work_nt * sizeof(*convolution));

    // 一个子源的全部源项先累加到同一响应，再统一卷积和执行时间算子
    real_t medium[3] = {Ctrl->H.vp, Ctrl->H.vs, Ctrl->H.rho};
    real_t *traces[12];
    for(int c = 0; c < nc; ++c) {
        traces[c] = response + (size_t)c * work_nt;
    }
    lamb_synthesis_samples(work_nt, dt, b, medium, Ctrl->H.nu, source->depth, receiver->depth, dist, az, source->fault,
                           area, Ctrl->n.active, Ctrl->e.active, Ctrl->L.phase_mask, traces);
    for(int c = 0; c < nc; ++c) {
        grt_oaconvolve(traces[c], work_nt, stf, stf_npts, convolution, work_nt, false);

        // 时间算子在完整响应上执行，最后才裁剪到输出窗口
        for(int j = 0; j < Ctrl->I.int_times; ++j) {
            grt_trap_integral(convolution, work_nt, dt);
        }
        for(int j = 0; j < Ctrl->J.dif_times; ++j) {
            grt_differential(convolution, work_nt, dt);
        }
        for(int n = 0; n < nt; ++n) {
            long long sample = first_output + n + margin;
            if(sample >= 0 && sample < work_nt - tail) {
                data[(size_t)c * nt + n] += convolution[sample] * dt;
            }
        }
    }
    GRT_SAFE_FREE_PTR(response);
    GRT_SAFE_FREE_PTR(convolution);
}

/**
 * 合成多点震源到单个接收点，并立即保存该点的结果
 *
 * @param[in]  Ctrl     命令行参数
 * @param[in]  output   动态合成输出设置
 * @param[in]  nsrc     震源点数
 * @param[in]  srcs     源点集合
 * @param[in]  ir       当前接收点索引
 * @param[in]  timing   当前接收点的时间窗及最早初至
 * @param[in]  threads  源点线程数
 */
static void lamb_one_receiver(const GRT_MODULE_CTRL *Ctrl, const DY_SYN_OUTPUT *output, size_t nsrc, const SRC_POINT *srcs,
                              size_t ir, const LAMB_RECEIVER_TIMING *timing, int threads)
{
    const RCV_POINT *receiver = &output->rcvs[ir];
    int nt = Ctrl->N.nt, nc = output->calc_upar ? 12 : 3;

    // 每个线程独占一段完整接收波形，避免子源叠加时竞争写入
    size_t samples = (size_t)nc * nt;
    real_t *data = GRT_SAFE_CALLOC((size_t)threads * samples, sizeof(*data));

    // 各线程独立求解和累加子源，共用只读几何及震源时间函数
    #pragma omp parallel num_threads(threads) if(nsrc > 1)
    {
        int tid = grt_get_thread_index();
        real_t *local = data + (size_t)tid * samples;

        #pragma omp for schedule(guided)
        for(size_t is = 0; is < nsrc; ++is) {
            // 无有效源项的子源不参与计算和叠加
            if(!srcs[is].fault->nterms) {
                continue;
            }
            lamb_one_pair(Ctrl, &srcs[is], receiver, timing->begin, local);
        }
    }

    // 按固定线程顺序归约，第一段缓冲保存最终结果
    for(int t = 1; t < threads; ++t) {
        for(size_t i = 0; i < samples; ++i) {
            data[i] += data[(size_t)t * samples + i];
        }
    }

    // 多源只记录 P、S 最早初至，单源保留其后续震相并计入相同破裂延迟
    SACTRACE *trace = grt_new_SACTRACE(Ctrl->N.dt, nt, timing->begin);
    GRT_SACHEAD_SET_IMAG_FREQ(&trace->hd, 0);
    real_t modarr[1][GRT_MODARR_NCOL] = {{0, Ctrl->H.vp, Ctrl->H.vs, Ctrl->H.rho, 0, 0}};
    grt_syn_output_set_receiver_header(output, ir, &trace->hd, Ctrl->C.active ? SAC_FLOAT_UNDEF : Ctrl->Depth.depsrc, 1, modarr);
    clear_lamb_arrivals(trace);

    // 单个有效子源保留完整震相信息，多源仅记录 P、S 最早初至
    if(nsrc == 1 && srcs[0].fault->nterms) {
        real_t dist = hypot(receiver->north - srcs[0].north, receiver->east - srcs[0].east);
        real_t distance = hypot(dist, srcs[0].depth - receiver->depth);
        set_lamb_arrivals(trace, Ctrl->H.nu, dist, srcs[0].depth, receiver->depth, distance, Ctrl->H.vs, srcs[0].fault->stf_delay);
    }

    // P 初至有效时写入头段，否则保留 SAC 未定义值
    if(timing->travtPS[0] != INFINITY) {
        set_lamb_arrival(trace, LAMB_PHASE_P, timing->travtPS[0], 1, 0, "P");
    }

    // S 初至有效时写入头段，否则保留 SAC 未定义值
    if(timing->travtPS[1] != INFINITY) {
        set_lamb_arrival(trace, LAMB_PHASE_S, timing->travtPS[1], 1, 0, "S");
    }
    grt_syn_output_save_receiver(output, ir, trace, data, Ctrl->N.dt, 0, 0);
    grt_free_SACTRACE(trace);
    GRT_SAFE_FREE_PTR(data);
}

/** 模块主函数 */
int lamb_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);

    // 读入时直接选择全局或行内时间函数，点源与有限源使用相同接口
    if(Ctrl->C.active) {
        Ctrl->C.faults = grt_finite_fault_from_option(Ctrl->C.option, &Ctrl->C.nfault, &Ctrl->C.dL, &Ctrl->C.dW, true, Ctrl->N.dt, 1, Ctrl->D.option);
    } else {
        Ctrl->C.nfault = 1;
        Ctrl->C.faults = grt_finite_fault_from_point(Ctrl->Depth.depsrc, Ctrl->source_type, Ctrl->S.scale,
                                                    Ctrl->S.mult_src_mu, Ctrl->mechanism, Ctrl->N.dt, 1, Ctrl->D.option);
    }

    // 统一展开震源和接收点，预检全部源台几何
    size_t nsrc, nrcv;
    real_t modarr[1][GRT_MODARR_NCOL] = {{0, Ctrl->H.vp, Ctrl->H.vs, Ctrl->H.rho, 0, 0}};
    SRC_POINT *srcs = build_lamb_sources(Ctrl, modarr, &nsrc);

    // 长震源可扩展指定的输出点数，所有接收点共用同一长度，卷积仍使用线性卷积
    Ctrl->N.nt = grt_syn_output_npts(Ctrl->N.nt, Ctrl->N.dt, Ctrl->C.nfault, Ctrl->C.faults);
    RCV_POINT *rcvs = build_lamb_receivers(Ctrl, &nrcv);
    LAMB_RECEIVER_TIMING *timings = prepare_lamb_geometry_and_timing(Ctrl, nsrc, srcs, nrcv, rcvs);
    DY_SYN_OUTPUT output = {
        .root = Ctrl->O.s_output_dir, .rcv_subdirs = Ctrl->Q.active || Ctrl->U.active, .rcvs = rcvs, .rcv_faults = Ctrl->U.faults,
        .channels = Ctrl->n.active ? "ZNE" : "ZRT", .calc_upar = Ctrl->e.active,
    };

    // 逐个接收点合成并保存，仅内部的源点循环并行
    int threads = grt_get_num_threads(nsrc);
    for(size_t ir = 0; ir < nrcv; ++ir) {
        lamb_one_receiver(Ctrl, &output, nsrc, srcs, ir, &timings[ir], threads);
    }

    // 总震源时间函数只保存一次，避免按接收点重复累加矩率
    if(Ctrl->C.faults[0].stf_explicit) {
        grt_syn_output_save_signal(output.root, Ctrl->N.dt, 1, Ctrl->C.nfault, Ctrl->C.faults, Ctrl->C.active);
    }

    // 非静默模式下报告展开后的源点数、接收点数和线程数
    if (!Ctrl->s.active) {
        GRTRaiseInfo("Synthesized %zu source point(s), %zu receiver(s), %d source thread(s).", nsrc, nrcv, threads);
    }

    // 所有接收结果保存完成后，释放时间窗和源台点集
    GRT_SAFE_FREE_PTR(timings);
    GRT_SAFE_FREE_PTR(srcs);
    GRT_SAFE_FREE_PTR(rcvs);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
