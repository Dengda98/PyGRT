/**
 * @file   grt_syn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12-2
 * 
 *    根据计算好的格林函数，定义震源机制以及方位角等，生成合成的三分量地震图
 * 
 */

#include "grt.h"

// 防止被替换为虚数单位
#undef I

/** 该子模块的参数控制结构体 */
typedef struct {
    /** 格林函数路径 */
    struct {
        bool active;
        char *s_grnpath;
    } G;
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
    /** 旋转到 Z, N, E */
    struct {
        bool active;
    } N;
    /** 放大系数 */
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
    /** 积分次数 */
    struct {
        bool active;
        int int_times;
    } I;
    /** 求导次数 */
    struct {
        bool active;
        int dif_times;
    } J;
    /** 时间函数 */
    struct {
        bool active;
        char *option;  ///< 完整 -D 时间函数选项
    } D;
    /** 根目录检索时的震源和台站深度 */
    struct {
        bool s_active;
        bool r_active;
        real_t depsrc;
        real_t deprcv;
    } Depth;
    /** 根目录检索时的震中距 */
    struct {
        bool active;
        real_t dist;
    } R;
    /** 静默输出 */
    struct {
        bool active;
    } s;
    /** 是否计算空间导数 */
    struct {
        bool active;
    } e;

    /** 接收点文件 */
    struct {
        bool active;
        char *path;  ///< 接收点文件
    } Q;

    /** 有限震源断层 */
    struct {
        bool active;
        bool has_i;           ///< 是否显式指定 +i
        char *option;         ///< 有限震源选项，确定采样间隔后读入
        FINITE_FAULT *faults;  ///< 震源断层数组
        size_t nfault;         ///< 震源断层数量
        real_t dL;             ///< 走向剖分间隔，km
        real_t dW;             ///< 倾向剖分间隔，km
    } C;

    /** 有限接收断层 */
    struct {
        bool active;
        bool has_i;           ///< 是否显式指定 +i
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

    /** 格林函数库查询方式 */
    struct {
        bool active;
        bool interpolate;  ///< 是否对角点合成结果插值
    } i;

    /** 单点源临时参数 */
    GRT_SYN_TYPE source_type;             ///< 单点源类型
    real_t mechanism[GRT_MECHANISM_NUM];  ///< 单点源机制参数

} GRT_MODULE_CTRL;

/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    // G
    GRT_SAFE_FREE_PTR(Ctrl->G.s_grnpath);

    // O
    GRT_SAFE_FREE_PTR(Ctrl->O.s_output_dir);

    // D
    GRT_SAFE_FREE_PTR(Ctrl->D.option);

    // 有限震源和接收点
    GRT_SAFE_FREE_PTR(Ctrl->C.option);
    grt_finite_fault_free(Ctrl->C.nfault, Ctrl->C.faults);
    grt_finite_fault_free(Ctrl->U.nfault, Ctrl->U.faults);
    GRT_SAFE_FREE_PTR(Ctrl->Q.path);
    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt syn] %s\n\n", GRT_VERSION);printf(
"    Compute three-component dynamic displacement from `greenfn` outputs.\n"
"\n"
"    Define five groups of information:\n"
"      1. Green functions: -G selects the library root or one node directory.\n"
"      2. Source location: a point source is at the horizontal origin; -Ds\n"
"         selects its depth. For finite sources, -C supplies all locations.\n"
"      3. Receiver locations: choose polar coordinates (-A, -R, -Dr), a point\n"
"         file (-Q), or receiver faults (-U). These three modes are exclusive.\n"
"      4. Source mechanism and strength: a point source requires -S; add at most\n"
"         one of -M (shear/tensile), -F (force), or -T (moment tensor).\n"
"         Without -M/-F/-T, use an explosion. -C supplies both location and\n"
"         mechanism/strength, so omit -Ds/-S/-M/-F/-T when using source faults.\n"
"      5. Output: -O sets the SAC directory; -N selects ZNE and -e adds spatial\n"
"         derivatives. Finite sources always output ZNE.\n"
"\n"
"    With a library root, -Ds/-Dr/-R may be inferred from singleton dimensions.\n"
"    With one node directory, omit -Ds/-Dr/-R; supply -A and the source mechanism.\n"
"    A node directory is used directly, and -i is ignored in that mode.\n"
"    Finite sources and -Q/-U receivers require a library root.\n"
"    -C requires at least two GF nodes, including point-source Kode records.\n"
"    All coordinates share one horizontal origin; with -C, -R/-A locate the\n"
"    receiver relative to that origin. +i on -C/-U controls fault subdivision.\n"
"    Time processing is optional: -D sets a time function, -I/-J integrate or\n"
"    differentiate in time. -i controls GF interpolation and -P sets threads.\n"
"    Default outputs are impulse-like displacements in cm, with Z upward,\n"
"    R radial outward and T clockwise from R.\n"
"\n"
"\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    # Point source\n"
"    grt syn -G<grn_path> -S[u]<scale> -O<outdir> [-Ds<depsrc>]\n"
"            [-M<strike>/<dip>[/<rake>] | -T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>\n"
"             | -F<fn>/<fe>/<fz>] <receiver options> [common options]\n"
"\n"
"    # Finite source\n"
"    grt syn -G<grn_root> -C<fault>[+i<dL>/<dW>] -O<outdir>\n"
"            <receiver options> [common options]\n"
"\n"
"    Receiver options (choose one mode):\n"
"      -A<azimuth> [-R<dist>] [-Dr<deprcv>]\n"
"      -Q<points>\n"
"      -U<fault>[+i<dL>/<dW>]\n"
"\n"
"    Common options:\n"
"      [-D<tftype>[/<tfparams>][+d<delay>]] [-I<odr>] [-J<odr>]\n"
"      [-i<0|1>] [-P<nthreads>] [-N] [-e] [-s] [-h]\n"
"\n"
"\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<grn_path>  Green's Functions output directory of module `greenfn`.\n"
"                  Output length is the maximum of GF and source-process lengths.\n"
"                  All receivers share that length; convolution is circular in time.\n"
"                  A single-distance subdirectory may be given directly.\n"
"                  In that mode, -Ds, -Dr and -R cannot be used.\n"
"                  The node is used directly; -i has no effect in that mode.\n"
"                  Finite sources and explicit receivers require a library root.\n"
"                  The root must contain every combination of source depth,\n"
"                  receiver depth and distance; each axis may be nonuniform.\n"
"                  With a root, -Ds selects point-source depth; -Dr/-R select\n"
"                  polar-receiver depth/distance. Each may be omitted when its\n"
"                  library dimension has one value. -C/-Q/-U supply their own\n"
"                  coordinates. Root queries default to linear interpolation\n"
"                  for all source/receiver types; -i0 selects nearest nodes.\n"
"\n"
"    -Ds<depsrc>   Source depth (km) used with a Green's-function root.\n"
"                  Required for multiple source depths and optional for one.\n"
"                  It cannot be used when -G points to a subdirectory.\n"
"                  It cannot be combined with finite-source option -C.\n"
"\n"

"    -Dr<deprcv>   Receiver depth (km) used with a Green's-function root.\n"
"                  Required for multiple receiver depths and optional for one.\n"
"                  It cannot be used when -G points to a subdirectory.\n"
"                  Omit with -Q/-U, which supply receiver depths.\n"
"\n"

"    -R<dist>      Polar receiver distance (km) used with a root directory.\n"
"                  Required for multiple distances and optional for one.\n"
"                  It cannot be used when -G points to a subdirectory.\n"
"                  With -C, this is the receiver distance from the horizontal\n"
"                  origin; source-receiver distances are computed from coordinates.\n"
"                  -R/-A/-Dr, -Q and -U are mutually exclusive receiver modes.\n"
"\n"
"    -A<azimuth>   Azimuth in degree, from source to station.\n"
"                  Ignored (forced to 0°) when Green's Functions\n"
"                  have zero epicentral distance.\n"
"                  With finite sources, this is the polar receiver azimuth\n"
"                  from the common horizontal origin. Omit with -Q/-U.\n"
"\n"
"    -S[u]<scale>  Scale factor to all kinds of source. \n"
"                  + For Explosion, Shear and Moment Tensor,\n"
"                    unit of <scale> is dyne-cm.\n"
"                  + For Single Force, unit of <scale> is dyne.\n"
"                  + Since \"\\mu\" exists in scalar seismic moment\n"
"                    (\\mu*A*D), you can simply set -Su<scale>, <scale>\n"
"                    equals A*D (Area*Slip, [cm^3]), and <scale> will \n"
"                    multiply \\mu automatically in program.\n"
"                  Use -C instead of -S/-M/-F/-T for finite sources.\n"
"\n"
"    For source type, you can only set at most one of\n"
"    '-M', '-T' and '-F'. If none, an Explosion is used.\n"
"\n"
"    -M<strike>/<dip>[/<rake>]\n"
"                  Three angles to define a shear fault. \n"
"                  The angles are in degree.\n"
"                  If <rake> not set, then define a tensile fault.\n"
"\n"
"    -T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>\n"
"                  Six elements of Moment Tensor. \n"
"                  x (North), y (East), z (Downward).\n"
"                  Notice they will be scaled by <scale>.\n"
"\n"
"    -F<fn>/<fe>/<fz>\n"
"                  North, East and Vertical(Downward) Forces.\n"
"                  Notice they will be scaled by <scale>.\n"
"\n"
"    -O<outdir>    Directory of output for saving. Default is\n"
"                  current directory.\n"
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
"                    e.g. \n"
"                         -D%c/1.3\n", GRT_SIG_PARABOLA); printf(
"                  + Trapezoidal wave\n"
"                    set -D%c/<t1>/<t2>/<t3>, rise/plateau/fall durations in seconds.\n", GRT_SIG_TRAPEZOID); printf(
"                    Durations must be nonnegative, with a positive total.\n"
"                    t2=0 gives a triangle; t1=t3=0 gives a rectangle.\n"
"                    e.g. \n"
"                         -D%c/0.1/0.1/0.2\n", GRT_SIG_TRAPEZOID); printf(
"                         -D%c/0.4/0/0.2 (become a triangle)\n", GRT_SIG_TRAPEZOID); printf(
"                         -D%c/0/0.5/0 (become a rectangle)\n", GRT_SIG_TRAPEZOID); printf(
"                  + AsymmetricCosine\n"
"                    set -D%c/<t1>/<t2>, positive rise/fall durations in seconds.\n", GRT_SIG_ASYMMETRIC_COSINE); printf(
"                    The peak is at t1 and the end is at t1+t2.\n"
"                    e.g. -D%c/4.5/1.5+d43\n", GRT_SIG_ASYMMETRIC_COSINE); printf(
"                  + Custom wave\n"
"                    set -D%c/<path>, <path> is the filepath to a custom\n", GRT_SIG_CUSTOM); printf(
"                    Time Function ASCII file with one amplitude column sampled\n"
"                    at dt, or two columns: time (s) and amplitude. Two-column\n"
"                    data are linearly interpolated to dt. Custom time functions\n"
"                    are automatically area-normalized.\n"
"                    e.g. \n"
"                         -D%c/tfunc.txt \n", GRT_SIG_CUSTOM); printf(
"                  Also accepts a signed Ricker convolution wavelet:\n"
"                  -D%c/<f0>, <f0> is the dominant frequency in Hz.\n", GRT_SIG_RICKER); printf(
"                  Its analytic peak amplitude is 1, without area normalization;\n"
"                  it is not a unit-slip source process.\n"
"                  To match the time interval in Green's Functions, \n"
"                  parameters of Time Function will be slightly modified.\n"
"                  The corresponding Time Function will be saved\n"
"                  as a SAC file under <outdir>.\n"
"\n"
"                  Append +d<delay> for rupture delay in seconds.\n"
"                  Append a complete -D option at the end of each fault row\n"
"                  to specify its rupture process.\n"
"\n"
"    -I<odr>       Order of integration. Default not use\n"
"\n"
"    -J<odr>       Order of differentiation. Default not use\n"
"\n"
"    -N            Components of results will be Z, N, E.\n"
"                  Finite sources always output ZNE; -N is forced with a warning if omitted.\n"
"\n"
"    -e            Compute the spatial derivatives, ui_z and ui_r,\n"
"                  of displacement u. In filenames, prefix \"r\" means \n"
"                  ui_r and \"z\" means ui_z. \n"
"                  The azimuthal derivative is also synthesized. Direction prefixes\n"
"                  are z/r/t in ZRT or z/n/e in ZNE.\n"
"\n"
"    -s            Silence all outputs.\n"
"\n"
"    -C<fault>[+i<dL>/<dW>]\n"
"                  Coulomb source faults, including point-source Kode records.\n"
"                  Requires a library root with at least two nodes; +i does not\n"
"                  allow finite sources with a single-node library.\n"
"                  The file supplies location, mechanism and signed slip/potency;\n"
"                  do not set -Ds/-S/-M/-F/-T. Finite sources always output ZNE.\n"
"                  Kode 100/200/300: rectangular shear/tensile sources.\n"
"                  Kode 400: point double couple; Kode 500: point tensile/inflation.\n"
"                  Two header lines precede the 11 numeric columns. An exact\n"
"                  \"rake\" token in the seventh header column selects rake/net\n"
"                  slip for Kode 100; the filename suffix does not select format.\n"
"                  +i gives along-strike/dip subdivision sizes (km).\n"
"                  Sizes must both be positive or both be zero. +i0/0 uses\n"
"                  one center point source weighted by the whole fault area.\n"
"                  Point-source Kode records remain single points at fault centers.\n"
"                  Without +i, rectangular sizes default to the smallest positive\n"
"                  distance/source-depth/receiver-depth interval in the GF library.\n"
"\n"
"    -Q<points>    ASCII rows: north east depth (km) [strike dip rake (degrees)].\n"
"                  Lines starting with # are comments. Optional angles are saved\n"
"                  for later stress projection and do not affect synthesis.\n"
"\n"
"    -U<fault>[+i<dL>/<dW>]\n"
"                  Coulomb receiver faults. Receivers are subfault centers;\n"
"                  no receiver-area averaging is performed.\n"
"                  Without +i or with +i0/0, use each fault center.\n"
"                  Slip magnitude is ignored.\n"
"                  Only Kode=100 is supported;\n"
"                  receiver angles and grouping are saved.\n"
"\n"
"    -i<0|1>       Library-root query method, independent of source/receiver count.\n"
"                  Default 1: linearly combine synthesized corner results.\n"
"                  0: use nearest nodes, including for one source/one receiver.\n"
"                  Ignored when -G points directly to a node directory.\n"
"                  All geometries must lie inside the library. Combined results\n"
"                  require consistent nt/dt and imaginary frequency.\n"
"\n"
"    -P<nthreads>  OpenMP source-point threads. Receivers are processed serially.\n"
"\n"
"    -h            Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    Say you have computed Green's functions with following command:\n"
"        grt greenfn -Mmilrow -N1000/0.01 -D2/0 -Ores -R2,4,6,8,10\n"
"\n"
"    Then you can get synthetic seismograms of Explosion at epicentral\n"
"    distance of 10 km and an azimuth of 30° by running:\n"
"        grt syn -Gres/milrow_2_0_10 -Osyn_ex -A30 -S1e24\n"
"\n"
"    or Shear\n"
"        grt syn -Gres/milrow_2_0_10 -Osyn_dc -A30 -S1e24 -M100/20/80\n"
"\n"
"    or Tension\n"
"        grt syn -Gres/milrow_2_0_10 -Osyn_dc -A30 -S1e24 -M100/20\n"
"\n"
"    or Single Force\n"
"        grt syn -Gres/milrow_2_0_10 -Osyn_sf -A30 -S1e24 -F0.5/-1.2/3.3\n"
"\n"
"    or Moment Tensor\n"
"        grt syn -Gres/milrow_2_0_10 -Osyn_mt -A30 -S1e24 -T2.3/0.2/-4.0/0.3/0.5/1.2\n"
"\n\n\n"
"    Finite sources with explicit receivers:\n"
"        grt syn -Gres -Cfaults.inp+i1/1 -Qreceivers.txt -Osyn_ff\n"

);
}

/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    // 先为个别参数设置非0初始值
    Ctrl->source_type = GRT_SYN_EX;
    Ctrl->i.interpolate = true;

    int opt, mechanisms = 0;
    char extra;
    while ((opt = getopt(argc, argv, ":G:R:A:S:M:F:T:O:D:I:J:C:U:Q:i:P:Nesh")) != -1) {
        switch (opt) {
            // 格林函数路径
            case 'G':
                Ctrl->G.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->G.s_grnpath);
                Ctrl->G.s_grnpath = strdup(optarg);

                // 检查是否存在该目录
                GRTCheckDirExist(Ctrl->G.s_grnpath);
                break;

            // 方位角
            case 'A':
                Ctrl->A.active = true;
                if(1 != sscanf(optarg, "%lf%c", &Ctrl->A.azimuth, &extra)){
                    GRTBadOptionError(A, "");
                };
                if(Ctrl->A.azimuth < 0.0 || Ctrl->A.azimuth > 360.0){
                    GRTBadOptionError(A, "Azimuth must be in [0, 360].");
                }
                break;

            // 放大系数
            case 'S':
                Ctrl->S.active = true;

                // 检查是否存在字符u，若存在表明需要乘上震源处的剪切模量
                Ctrl->S.mult_src_mu = optarg[0] == 'u';
                if(1 != sscanf(optarg + Ctrl->S.mult_src_mu, "%lf%c", &Ctrl->S.scale, &extra)){
                    GRTBadOptionError(S, "");
                };
                break;

            // 剪切震源， 张裂源
            case 'M':
                Ctrl->M.active = true;
                {
                    real_t strike=0.0, dip=0.0, rake=0.0;
                    int nscan = sscanf(optarg, "%lf/%lf/%lf%c", &strike, &dip, &rake, &extra);
                    if((nscan == 2 || nscan == 3) && (nscan != 2 || sscanf(optarg, "%lf/%lf%c", &strike, &dip, &extra) == 2)){
                        Ctrl->source_type = GRT_SYN_TS;
                        if(strike < 0.0 || strike > 360.0){
                            GRTBadOptionError(M, "Strike must be in [0, 360].");
                        }
                        if(dip < 0.0 || dip > 90.0){
                            GRTBadOptionError(M, "Dip must be in [0, 90].");
                        }
                        if(nscan == 3){
                            Ctrl->source_type = GRT_SYN_DC;
                            if(rake < -180.0 || rake > 180.0){
                                GRTBadOptionError(M, "Rake must be in [-180, 180].");
                            }
                        }
                    } else {
                        GRTBadOptionError(M, "");
                    };


                    Ctrl->mechanism[0] = strike;
                    Ctrl->mechanism[1] = dip;
                    Ctrl->mechanism[2] = rake;
                }
                ++mechanisms;
                break;

            // 单力源
            case 'F':
                Ctrl->F.active = true;
                Ctrl->source_type = GRT_SYN_SF;
                {
                    real_t fn, fe, fz;
                    if(3 != sscanf(optarg, "%lf/%lf/%lf%c", &fn, &fe, &fz, &extra)){
                        GRTBadOptionError(F, "");
                    };
                    Ctrl->mechanism[0] = fn;
                    Ctrl->mechanism[1] = fe;
                    Ctrl->mechanism[2] = fz;
                }
                ++mechanisms;
                break;

            // 张量震源
            case 'T':
                Ctrl->T.active = true;
                Ctrl->source_type = GRT_SYN_MT;
                {
                    real_t Mxx, Mxy, Mxz, Myy, Myz, Mzz;
                    if(6 != sscanf(optarg, "%lf/%lf/%lf/%lf/%lf/%lf%c", &Mxx, &Mxy, &Mxz, &Myy, &Myz, &Mzz, &extra)){
                        GRTBadOptionError(T, "");
                    };
                    Ctrl->mechanism[0] = Mxx;
                    Ctrl->mechanism[1] = Mxy;
                    Ctrl->mechanism[2] = Mxz;
                    Ctrl->mechanism[3] = Myy;
                    Ctrl->mechanism[4] = Myz;
                    Ctrl->mechanism[5] = Mzz;
                }
                ++mechanisms;
                break;

            // 输出路径
            case 'O':
                Ctrl->O.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->O.s_output_dir);
                Ctrl->O.s_output_dir = strdup(optarg);
                break;

            // 点源深度、接收深度或震源时间函数，按 -D 的前缀区分
            case 'D':
                if(optarg[0] == 's'){
                    Ctrl->Depth.s_active = true;
                    if(strchr(optarg + 1, '/') != NULL || 1 != sscanf(optarg + 1, "%lf%c", &Ctrl->Depth.depsrc, &extra)){
                        GRTBadOptionError(Ds, "");
                    }
                    if(Ctrl->Depth.depsrc < 0.0){
                        GRTBadOptionError(Ds, "Negative source depth is not supported.");
                    }
                } else if(optarg[0] == 'r' && optarg[1] != '/'){
                    Ctrl->Depth.r_active = true;
                    if(strchr(optarg + 1, '/') != NULL || 1 != sscanf(optarg + 1, "%lf%c", &Ctrl->Depth.deprcv, &extra)){
                        GRTBadOptionError(Dr, "");
                    }
                    if(Ctrl->Depth.deprcv < 0.0){
                        GRTBadOptionError(Dr, "Negative receiver depth is not supported.");
                    }
                } else {
                    GRT_SAFE_FREE_PTR(Ctrl->D.option);
                    Ctrl->D.active = true;
                    GRT_SAFE_ASPRINTF(&Ctrl->D.option, "-D%s", optarg);
                }
                break;

            // 根目录检索时指定震中距
            case 'R':
                Ctrl->R.active = true;
                if(1 != sscanf(optarg, "%lf%c", &Ctrl->R.dist, &extra) || Ctrl->R.dist < 0.0){
                    GRTBadOptionError(R, "Nonnegative epicentral distance is required.");
                }
                break;

            // 对结果做积分
            case 'I':
                Ctrl->I.active = true;
                if(1 != sscanf(optarg, "%d%c", &Ctrl->I.int_times, &extra)){
                    GRTBadOptionError(I, "");
                }
                if(Ctrl->I.int_times <= 0){
                    GRTBadOptionError(I, "Order should be positive.");
                }
                break;

            // 对结果做微分
            case 'J':
                Ctrl->J.active = true;
                if(1 != sscanf(optarg, "%d%c", &Ctrl->J.dif_times, &extra)){
                    GRTBadOptionError(J, "");
                }
                if(Ctrl->J.dif_times <= 0){
                    GRTBadOptionError(J, "Order should be positive.");
                }
                break;

            // 是否计算位移空间导数, 影响 calcUTypes 变量
            case 'e':
                Ctrl->e.active = true;
                break;

            // 是否旋转到ZNE, 影响 rot2ZNE 变量
            case 'N':
                Ctrl->N.active = true;
                break;

            // 不打印在终端
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
                Ctrl->U.faults = grt_finite_fault_from_option(optarg, &Ctrl->U.nfault, &Ctrl->U.has_i, &Ctrl->U.dL, &Ctrl->U.dW,
                                                              false, 0, 1, NULL);
                break;

            // 是否对格林函数库进行三维线性插值
            case 'i':
                Ctrl->i.active = true;
                if(strcmp(optarg, "0") && strcmp(optarg, "1")) {
                    GRTRaiseError("-i expects 0 or 1.");
                }
                Ctrl->i.interpolate = atoi(optarg);
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

    // 检查必选项有没有设置
    GRTCheckOptionActive(Ctrl, G);
    GRTCheckOptionActive(Ctrl, O);

    // 有限震源文件不能与点源强度、机制或源深度选项同时使用
    if(Ctrl->C.active && (Ctrl->S.active || mechanisms || Ctrl->Depth.s_active)) {
        GRTRaiseError("Finite sources and point-source parameters are mutually exclusive.");
    }

    // 点源必须显式指定震源强度
    if(!Ctrl->C.active && !Ctrl->S.active) {
        GRTRaiseError("Point sources require -S.");
    }

    // 点源机制最多只能选用 -M、-F、-T 中的一种
    if(mechanisms > 1) {
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

    // 有限震源的各子源方位不同，叠加时统一使用 ZNE 坐标系
    if(Ctrl->C.active && !Ctrl->N.active) {
        GRTRaiseWarning("Finite sources (-C) require ZNE components; -N is enabled automatically.");
        Ctrl->N.active = true;
    }
}

/**
 * 合成一个源项，在角点处完成联络项和坐标转换后累加
 *
 * 数组布局：gf[基本源型][分量][采样点]、syn[分量][采样点]、
 * syn_upar[偏导方向][分量][采样点]
 *
 * @param[in]      npts         时间样本数
 * @param[in]      dist         库角点震中距，km
 * @param[in]      target_dist  实际源台震中距，km
 * @param[in]      gf           位移格林函数
 * @param[in]      gf_uiz       深度偏导格林函数，未计算导数时可为 NULL
 * @param[in]      gf_uir       距离偏导格林函数，未计算导数时可为 NULL
 * @param[in]      source_type  源类型
 * @param[in]      scale        源强
 * @param[in]      VpVs_ratio   源点速度比
 * @param[in]      azrad        源台方位角，弧度
 * @param[in]      mchn         源机制参数
 * @param[in]      rot2ZNE      是否输出 ZNE 分量
 * @param[in]      calc_upar    是否计算空间导数
 * @param[in,out]  syn          位移累加数组
 * @param[in,out]  syn_upar     导数累加数组
 */
static void syn_accumulate_term(
    size_t npts, real_t dist, real_t target_dist,
    const prealChnlGrid gf, const prealChnlGrid gf_uiz, const prealChnlGrid gf_uir,
    GRT_SYN_TYPE source_type, real_t scale, real_t VpVs_ratio, real_t azrad,
    const real_t mchn[GRT_MECHANISM_NUM], bool rot2ZNE, bool calc_upar,
    real_t *const syn[GRT_CHANNEL_NUM], real_t *const syn_upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    int ntypes = calc_upar ? 4 : 1;
    realChnlGrid radiation[4] = {0};
    real_t *const (*grns[4])[GRT_CHANNEL_NUM] = {gf, gf_uiz, gf_uir, GRT_IS_ZERO(dist) ? gf_uir : gf};
    for(int ityp = 0; ityp < ntypes; ++ityp) {
        real_t derivative_scale = ityp == 0 ? 1 : 1e-5;
        if(ityp == 3 && !GRT_IS_ZERO(dist)) {
            derivative_scale /= dist;
        }
        // r=0 时用距离偏导 GF 合成角向导数的有限部分，联络项的有限极限由坐标转换统一处理
        grt_set_source_radiation(radiation[ityp], source_type, ityp == 3, scale, derivative_scale, VpVs_ratio, azrad, mchn);
    }

    for(size_t n = 0; n < npts; ++n) {
        real_t sample[4][GRT_CHANNEL_NUM] = {{0}};
        for(int ityp = 0; ityp < ntypes; ++ityp) {
            GRT_LOOP_ChnlGrid(im, c) {
                real_t coef = radiation[ityp][im][c];
                if(coef != 0) {
                    sample[ityp][c] += grns[ityp][im][c][n] * coef;
                }
            }
        }

        if(calc_upar) {
            // 先用角点距离得到直角坐标导数，再转换为实际位置的 ZRT 表示
            // ZRT 输出仍保存裸角向偏导，后处理在实际位置补联络项
            grt_rot_zrt2zxy_upar(azrad, sample[0], sample + 1, dist * 1e5);
            if(!rot2ZNE) {
                grt_rot_zxy2zrt_upar(azrad, sample[0], sample + 1, target_dist * 1e5);
            }
        } else if(rot2ZNE) {
            grt_rot_zxy2zrt_vec(-azrad, sample[0]);
        }

        for(int c = 0; c < GRT_CHANNEL_NUM; ++c) {
            syn[c][n] += sample[0][c];
            if(calc_upar) {
                for(int d = 0; d < GRT_CHANNEL_NUM; ++d) {
                    syn_upar[d][c][n] += sample[d + 1][c];
                }
            }
        }
    }
}

/**
 * 根据输入路径确定查询方式，并补全可省略的几何选项
 * @param[in,out]  Ctrl  命令行参数
 * @param[in]      lib   动态格林函数库
 * @return         单节点直接匹配，根目录使用最近邻或线性插值
 */
static GRT_SAMPLE_MODE select_library_geometry(GRT_MODULE_CTRL *Ctrl, const DYGRNLIB *lib)
{
    // 按实际节点数量拒绝有限源，不区分单节点目录和仅含一个节点的根目录
    if(Ctrl->C.active && lib->nnode == 1) {
        GRTRaiseError("Finite sources (-C) require a Green-function library with at least two nodes.");
    }

    // 单节点目录直接确定源台几何，查询选项不影响匹配
    if(lib->direct_node_input) {
        if(Ctrl->Q.active || Ctrl->U.active) {
            GRTRaiseError("Receiver files (-Q/-U) require a Green-library root.");
        }
        if(Ctrl->Depth.s_active || Ctrl->Depth.r_active || Ctrl->R.active) {
            GRTRaiseError("Do not set depth/distance selectors with a GF subdirectory.");
        }
        Ctrl->Depth.depsrc = lib->nodes[0].zs;
        Ctrl->Depth.deprcv = lib->nodes[0].zr;
        Ctrl->R.dist = lib->nodes[0].r;
        return GRT_SAMPLE_EXACT;
    }

    // 点源未指定深度时只允许库中有一个源深度，并使用该默认值
    if(!Ctrl->C.active) {
        if(!Ctrl->Depth.s_active && lib->ndepsrc != 1) {
            GRTRaiseError("Multiple source depths require -Ds.");
        }
        if(!Ctrl->Depth.s_active) {
            Ctrl->Depth.depsrc = lib->depsrcs[0];
        }
    }

    // 极坐标接收点未指定深度时，只允许库中有一个接收深度
    if(!Ctrl->Q.active && !Ctrl->U.active) {
        if(!Ctrl->Depth.r_active && lib->ndeprcv != 1) {
            GRTRaiseError("Multiple receiver depths require -Dr.");
        }
        if(!Ctrl->Depth.r_active) {
            Ctrl->Depth.deprcv = lib->deprcvs[0];
        }
    }

    // 未指定接收文件和震中距时，只允许库中有一个距离采样点
    if(!Ctrl->Q.active && !Ctrl->U.active && !Ctrl->R.active) {
        if(lib->nr != 1) {
            GRTRaiseError("Multiple distances require -R or explicit receiver geometry.");
        }
        Ctrl->R.dist = lib->rs[0];
    }
    return Ctrl->i.interpolate ? GRT_SAMPLE_LINEAR : GRT_SAMPLE_NEAREST;
}

/**
 * 确定震源剖分尺寸并展开为统一源点集合
 *
 * @param[in]      Ctrl  命令行参数
 * @param[out]     nsrc  展开后的震源点数
 * @return         源点集合
 */
static SRC_POINT *build_syn_sources(GRT_MODULE_CTRL *Ctrl, const DYGRNLIB *lib, size_t *nsrc)
{
    // 只有矩形断层需要剖分尺寸，未指定时使用库内三个采样轴的最小正间隔
    bool rectangular = false;
    for(size_t i = 0; i < Ctrl->C.nfault; ++i) {
        rectangular |= KODE_IS_FINITE(Ctrl->C.faults[i].kode);
    }
    if(rectangular && !Ctrl->C.has_i) {
        Ctrl->C.dL = Ctrl->C.dW = grt_dygrnlib_default_subfault_size(lib);
    }

    // 逐条断层剖分，保存矩率时计算子源标量矩，点源保持单个子源
    for(size_t i = 0; i < Ctrl->C.nfault; ++i) {
        FINITE_FAULT *fault = &Ctrl->C.faults[i];
        grt_finite_fault_subdiv(fault, KODE_IS_POINT(fault->kode) ? 0 : Ctrl->C.dL, KODE_IS_POINT(fault->kode) ? 0 : Ctrl->C.dW,
                                lib->nlayer, lib->modarr);
    }
    return grt_src_points_from_faults(Ctrl->C.nfault, Ctrl->C.faults, nsrc);
}

/**
 * 将单点、逐点文件或有限接收断层统一为接收点集合
 *
 * @param[in,out]  Ctrl  命令行参数，补全默认剖分尺寸
 * @param[in]      lib   格林函数库
 * @param[out]     nrcv  展开后的接收点数
 * @return         接收点集合
 */
static RCV_POINT *build_syn_receivers(GRT_MODULE_CTRL *Ctrl, size_t *nrcv)
{
    // 逐点接收文件直接提供坐标和可选机制，无需剖分
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
 * 预检源台几何并确定各接收点的共同起点，角点查询结果仅在当前源台对内使用
 *
 * @param[in]   lib     格林函数库
 * @param[in]   nsrc    震源点数量
 * @param[in]   srcs    源点集合
 * @param[in]   nrcv    接收点数量
 * @param[in]   rcvs    接收点集合
 * @param[in]   mode    几何查询方式
 * @param[out]  begins  各接收点的共同起点，长度为 nrcv
 */
static void prepare_syn_timing(const DYGRNLIB *lib, size_t nsrc, const SRC_POINT *srcs,
                               size_t nrcv, const RCV_POINT *rcvs, GRT_SAMPLE_MODE mode, real_t *begins)
{
    for(size_t ir = 0; ir < nrcv; ++ir) {
        begins[ir] = INFINITY;
        for(size_t is = 0; is < nsrc; ++is) {
            const SRC_POINT *source = &srcs[is];
            const RCV_POINT *receiver = &rcvs[ir];
            real_t dist = hypot(receiver->north - source->north, receiver->east - source->east);
            size_t indices[8];
            real_t weights[8];
            int ncorners = grt_dygrnlib_find_corners(lib, source->depth, receiver->depth, dist, mode, indices, weights);
            if(ncorners <= 0) {
                GRTRaiseError("Out-of-range or unmatched GF geometry at source %zu receiver %zu: depths %.9g/%.9g, distance %.9g.",
                              is, ir, source->depth, receiver->depth, dist);
            }

            // 零强度子源仍检查几何，但不参与共同起点的选择
            for(int k = 0; k < ncorners && source->fault->nterms; ++k) {
                const DYGRNLIB_NODE *node = &lib->nodes[indices[k]];
                for(int t = 0; t < source->fault->nterms; ++t) {
                    const FINITE_SOURCE_TERM *term = &source->fault->terms[t];
                    grt_check_source_medium(term->type, term->with_mu, node->vs);
                }
                begins[ir] = GRT_MIN(begins[ir], node->b);
            }
        }
        if(begins[ir] == INFINITY) {
            begins[ir] = 0;
        }
    }
}

/**
 * 合成一个源台对，按角点累加到当前线程的接收波形
 *
 * @param[in]      Ctrl         命令行参数
 * @param[in]      lib          动态格林函数库
 * @param[in]      source       震源点
 * @param[in]      receiver     接收点
 * @param[in]      mode         几何查询方式
 * @param[in]      begin        当前接收点共同起点，s
 * @param[in]      nt           所有接收点共用的输出点数
 * @param[out]     travtPS      P、S 初至，含破裂延迟，缺失为 INFINITY
 * @param[in,out]  data         当前线程的接收波形，按分量连续排列
 * @param[in,out]  response     当前角点响应缓冲
 * @param[in,out]  convolution  当前线程复用的卷积缓冲，长度为 nt
 */
static void syn_one_pair(const GRT_MODULE_CTRL *Ctrl, const DYGRNLIB *lib, const SRC_POINT *source,
                         const RCV_POINT *receiver, GRT_SAMPLE_MODE mode,
                         real_t begin, int nt, real_t travtPS[2], real_t *data, real_t *response, real_t *convolution)
{
    real_t north = receiver->north - source->north;
    real_t east = receiver->east - source->east;
    real_t dist = hypot(north, east);
    real_t azrad = GRT_IS_ZERO(dist) ? 0 : atan2(east, north);
    size_t indices[8];
    real_t weights[8];
    int ncorners = grt_dygrnlib_find_corners(lib, source->depth, receiver->depth, dist, mode, indices, weights);
    const FINITE_FAULT *fault = source->fault;
    int grn_nt = lib->nt, nc = lib->calc_upar ? 12 : 3;
    real_t dt = lib->dt;
    real_t area = KODE_IS_FINITE(fault->kode) ? fault->width[source->isub] * fault->length[source->isub] : 1;

    real_t growth = 1 / lib->stf_decay;

    real_t *traces[12], *syn[3], *upar[3][3] = {{NULL}};
    for(int c = 0; c < nc; ++c) {
        traces[c] = response + (size_t)c * grn_nt;
    }
    for(int c = 0; c < 3; ++c) {
        syn[c] = traces[c];
        if(lib->calc_upar) {
            for(int d = 0; d < 3; ++d) {
                upar[d][c] = traces[3 + 3 * d + c];
            }
        }
    }

    // 实际源台几何决定初至，不随最近邻或插值方式改变
    if(lib->modarr) {
        MODEL1D *mod1d = grt_read_mod1d_from_modarr(lib->nlayer, lib->modarr, "syn", source->depth, receiver->depth, true);
        travtPS[0] = grt_compute_travt1d(mod1d->Thk, mod1d->Va, mod1d->n, mod1d->isrc, mod1d->ircv, dist);
        travtPS[1] = grt_compute_travt1d(mod1d->Thk, mod1d->Vb, mod1d->n, mod1d->isrc, mod1d->ircv, dist);
        grt_free_mod1d(mod1d);
    } else {
        travtPS[0] = lib->refhead.t0;
        travtPS[1] = lib->refhead.t1;
    }
    for(int phase = 0; phase < 2; ++phase) {
        travtPS[phase] = travtPS[phase] != INFINITY && travtPS[phase] != SAC_FLOAT_UNDEF ? travtPS[phase] + fault->stf_delay : INFINITY;
    }

    for(int k = 0; k < ncorners; ++k) {
        const DYGRNLIB_NODE *node = &lib->nodes[indices[k]];
        DYGRNLIB_WAVEFORMS *waveforms = grt_dygrnlib_read_node(lib, indices[k], fault);
        long long shift = llround((node->b - begin + fault->stf_delay) / dt);
        real_t src_mu = node->vs * node->vs * node->rho * 1e10;
        real_t vpvs = node->vs > 0 ? node->vp / node->vs : 0;

        // 各源项在角点处完成坐标转换，再统一卷积和对齐
        memset(response, 0, (size_t)nc * grn_nt * sizeof(*response));
        for(int t = 0; t < fault->nterms; ++t) {
            const FINITE_SOURCE_TERM *term = &fault->terms[t];
            real_t scale = term->scale * area * (term->with_mu ? src_mu : 1);
            syn_accumulate_term(grn_nt, node->r, dist, waveforms->gf[0], waveforms->gf[1], waveforms->gf[2],
                                term->type, scale, vpvs, azrad, term->mechanism, Ctrl->N.active, lib->calc_upar, syn, upar);
        }
        grt_dygrnlib_free_waveforms(waveforms);

        for(int c = 0; c < nc; ++c) {
            // 卷积前去除虚频率补偿，较短输入补零，以统一输出点数为周期做时域循环卷积
            real_t damping = 1;
            for(int n = 0; n < grn_nt; ++n) {
                traces[c][n] *= damping;
                damping *= lib->stf_decay;
            }
            grt_oaconvolve(traces[c], grn_nt, fault->stfd, fault->stf_npts, convolution, nt, true);
            // 恢复虚频率补偿并乘连续卷积的 dt 因子，破裂延迟仍按输出窗口平移和裁剪
            damping = 1;
            for(int n = 0; n < nt; ++n) {
                long long out = n + shift;
                if(out >= 0 && out < nt) {
                    data[(size_t)c * nt + out] += weights[k] * convolution[n] * damping * dt;
                }
                damping *= growth;
            }
        }
    }
}

/**
 * 合成多点震源到单个接收点，并立即保存该点的结果
 *
 * @param[in]  Ctrl     命令行参数
 * @param[in]  lib      格林函数库
 * @param[in]  output   动态合成输出设置
 * @param[in]  nsrc     震源点数
 * @param[in]  srcs     震源点集合
 * @param[in]  mode     几何查询方式
 * @param[in]  ir       当前接收点索引
 * @param[in]  begin    当前接收点的共同起点，s
 * @param[in]  nt       所有接收点共用的输出点数
 * @param[in]  threads  源点线程数
 */
static void syn_one_receiver(const GRT_MODULE_CTRL *Ctrl, const DYGRNLIB *lib, const DY_SYN_OUTPUT *output,
                             size_t nsrc, const SRC_POINT *srcs, GRT_SAMPLE_MODE mode, size_t ir, real_t begin, int nt, int threads)
{
    int nc = output->calc_upar ? 12 : 3;

    // 每个线程独占一段完整接收波形，避免子源叠加时竞争写入
    size_t samples = (size_t)nc * nt;
    real_t *data = GRT_SAFE_CALLOC((size_t)threads * samples, sizeof(*data));
    real_t firstP = INFINITY, firstS = INFINITY;

    // 各线程独立求解和累加子源，共用只读几何及震源时间函数
    #pragma omp parallel num_threads(threads) if(nsrc > 1) reduction(min:firstP,firstS)
    {
        // 缺失震相或线程未处理有效源时，保持无穷初值以区别于合法到时
        firstP = firstS = INFINITY;
        int tid = grt_get_thread_index();
        real_t *local = data + (size_t)tid * samples;
        real_t *response = GRT_SAFE_MALLOC((size_t)nc * lib->nt * sizeof(*response));
        real_t *convolution = GRT_SAFE_MALLOC(nt * sizeof(*convolution));

        #pragma omp for schedule(guided)
        for(size_t is = 0; is < nsrc; ++is) {
            // 跳过无有效源项的子源，避免更新波形和最早初至
            if(!srcs[is].fault->nterms) {
                continue;
            }
            real_t travtPS[2];
            syn_one_pair(Ctrl, lib, &srcs[is], &output->rcvs[ir], mode, begin, nt, travtPS, local, response, convolution);
            firstP = fmin(firstP, travtPS[0]);
            firstS = fmin(firstS, travtPS[1]);
        }
        GRT_SAFE_FREE_PTR(response);
        GRT_SAFE_FREE_PTR(convolution);
    }

    // 按固定线程顺序归约，第一段缓冲保存最终结果
    for(int t = 1; t < threads; ++t) {
        for(size_t i = 0; i < samples; ++i) {
            data[i] += data[(size_t)t * samples + i];
        }
    }

    // 爆炸点源和纯垂直单力点源具有轴对称性，ZRT 输出中的理论零分量不保留坐标转换的舍入残差
    bool axisymmetric = !Ctrl->C.active && (Ctrl->source_type == GRT_SYN_EX ||
                         (Ctrl->source_type == GRT_SYN_SF && Ctrl->mechanism[0] == 0 && Ctrl->mechanism[1] == 0));
    if(axisymmetric && !Ctrl->N.active) {
        const int zero_channels[] = {2, 5, 8, 9, 10, 11};  // T、zT、rT、tZ、tR、tT
        int nzero = output->calc_upar ? 6 : 1;
        for(int i = 0; i < nzero; ++i) {
            memset(data + (size_t)zero_channels[i] * nt, 0, nt * sizeof(*data));
        }
    }

    // 使用共同时间窗创建输出头段，再填写接收位置、介质和最早初至
    SACTRACE *trace = grt_new_SACTRACE(lib->dt, nt, begin);
    GRT_SACHEAD_SET_IMAG_FREQ(&trace->hd, 0);
    if(!lib->modarr) {
        // 单节点目录没有原始模型，介质直接沿用库的参考头段
        GRT_SACHEAD_SET_RCV_VP(&trace->hd,     GRT_SACHEAD_GET_RCV_VP(&lib->refhead));
        GRT_SACHEAD_SET_RCV_VS(&trace->hd,     GRT_SACHEAD_GET_RCV_VS(&lib->refhead));
        GRT_SACHEAD_SET_RCV_RHO(&trace->hd,    GRT_SACHEAD_GET_RCV_RHO(&lib->refhead));
        GRT_SACHEAD_SET_RCV_QP_INV(&trace->hd, GRT_SACHEAD_GET_RCV_QP_INV(&lib->refhead));
        GRT_SACHEAD_SET_RCV_QS_INV(&trace->hd, GRT_SACHEAD_GET_RCV_QS_INV(&lib->refhead));
        GRT_SACHEAD_SET_SRC_VP(&trace->hd,     GRT_SACHEAD_GET_SRC_VP(&lib->refhead));
        GRT_SACHEAD_SET_SRC_VS(&trace->hd,     GRT_SACHEAD_GET_SRC_VS(&lib->refhead));
        GRT_SACHEAD_SET_SRC_RHO(&trace->hd,    GRT_SACHEAD_GET_SRC_RHO(&lib->refhead));
    }
    grt_syn_output_set_receiver_header(output, ir, &trace->hd, Ctrl->C.active ? SAC_FLOAT_UNDEF : Ctrl->Depth.depsrc, lib->nlayer, lib->modarr);

    // 只记录实际存在的初至，缺失震相保留 SAC 未定义值
    if(firstP != INFINITY) {
        trace->hd.t0 = firstP;
        strcpy(trace->hd.kt0, "P");
    }

    // S 初至独立检查，避免缺失震相被记录为有效到时
    if(firstS != INFINITY) {
        trace->hd.t1 = firstS;
        strcpy(trace->hd.kt1, "S");
    }
    grt_syn_output_save_receiver(output, ir, trace, data, lib->dt, Ctrl->I.int_times, Ctrl->J.dif_times);
    grt_free_SACTRACE(trace);
    GRT_SAFE_FREE_PTR(data);
}



/** 子模块主函数 */
int syn_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));

    getopt_from_command(Ctrl, argc, argv);

    // 读取格林函数库，统一展开震源和接收点
    DYGRNLIB *lib = grt_dygrnlib_load(Ctrl->G.s_grnpath, Ctrl->e.active);
    GRT_SAMPLE_MODE mode = select_library_geometry(Ctrl, lib);

    size_t nsrc, nrcv;
    RCV_POINT *rcvs = build_syn_receivers(Ctrl, &nrcv);

    // 读入时直接选择全局或行内时间函数，点源与有限源使用相同接口
    if(Ctrl->C.active) {
        Ctrl->C.faults = grt_finite_fault_from_option(Ctrl->C.option, &Ctrl->C.nfault, &Ctrl->C.has_i, &Ctrl->C.dL, &Ctrl->C.dW,
                                                      true, lib->dt, lib->stf_decay, Ctrl->D.option);
    } else {
        Ctrl->C.nfault = 1;
        Ctrl->C.faults = grt_finite_fault_from_point(Ctrl->Depth.depsrc, Ctrl->source_type, Ctrl->S.scale,
                                                    Ctrl->S.mult_src_mu, Ctrl->mechanism, lib->dt, lib->stf_decay, Ctrl->D.option);
    }
    SRC_POINT *srcs = build_syn_sources(Ctrl, lib, &nsrc);

    // 输出点数覆盖格林函数和完整震源过程，只计算一次，不随接收点几何改变
    int nt = grt_syn_output_npts(lib->nt, lib->dt, Ctrl->C.nfault, Ctrl->C.faults);

    // 单点源台重合时无法定义方位角，使用固定的零度方向
    if(!Ctrl->C.active && nrcv == 1 && GRT_IS_ZERO(hypot(rcvs[0].north, rcvs[0].east))) {
        GRTRaiseWarning("Zero epicentral distance: azimuth is forced to 0 degrees.");
    }

    // 有限源统一输出 ZNE，位移与空间导数共用接收信息和保存设置
    DY_SYN_OUTPUT output = {
        .root = Ctrl->O.s_output_dir, .rcv_subdirs = Ctrl->Q.active || Ctrl->U.active, .rcvs = rcvs, .rcv_faults = Ctrl->U.faults,
        .channels = Ctrl->N.active ? "ZNE" : "ZRT", .calc_upar = Ctrl->e.active,
    };

    // 预检源台几何并确定各接收点的共同时间窗
    real_t *begins = GRT_SAFE_MALLOC(nrcv * sizeof(*begins));
    prepare_syn_timing(lib, nsrc, srcs, nrcv, rcvs, mode, begins);
    int threads = grt_get_num_threads(nsrc);

    // 逐个接收点合成并保存，仅内部的源点循环并行
    for(size_t ir = 0; ir < nrcv; ++ir) {
        syn_one_receiver(Ctrl, lib, &output, nsrc, srcs, mode, ir, begins[ir], nt, threads);
    }

    // 总震源时间函数只保存一次，避免按接收点重复累加矩率
    if(Ctrl->C.faults[0].stf_explicit) {
        grt_syn_output_save_signal(output.root, lib->dt, lib->stf_decay, Ctrl->C.nfault, Ctrl->C.faults, Ctrl->C.active);
    }

    // 非静默模式下报告展开后的源点数、接收点数和线程数
    if(!Ctrl->s.active) {
        GRTRaiseInfo("Synthesized %zu source point(s), %zu receiver(s), %d source thread(s).", nsrc, nrcv, threads);
    }

    GRT_SAFE_FREE_PTR(begins);
    GRT_SAFE_FREE_PTR(srcs);
    grt_dygrnlib_free(lib);
    GRT_SAFE_FREE_PTR(rcvs);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
