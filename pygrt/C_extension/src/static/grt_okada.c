/**
 * @file   grt_okada.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * 使用 Okada 均匀半空间解析解计算静态位移场
 *
 */

#include "grt.h"

#undef I

/** 该子模块的参数控制结构体 */
typedef struct {
    /** 均匀半空间介质参数 */
    struct {
        bool active;
        real_t vp;
        real_t vs;
        real_t rho;
    } I;
    /** 点源放大系数 */
    struct {
        bool active;
        bool mult_src_mu;
        real_t value;
    } S;
    /** 点源走向、倾角和滑动角 */
    struct {
        bool active;
        real_t strike;
        real_t dip;
        real_t rake;
        bool has_rake;
    } M;
    /** Coulomb 格式有限断层 */
    struct {
        bool active;
        real_t dL;
        real_t dW;
        FINITE_FAULT *faults;
        size_t nfault;
    } C;
    /** North 方向接收点坐标 */
    struct {
        bool active;
        size_t n;
        real_t *values;
    } X;
    /** East 方向接收点坐标 */
    struct {
        bool active;
        size_t n;
        real_t *values;
    } Y;
    /** 任意接收点文件 */
    struct {
        bool active;
        char *path;
    } Q;
    /** Coulomb 格式有限接收断层 */
    struct {
        bool active;
        real_t dL;
        real_t dW;
        size_t nfault;
        FINITE_FAULT *faults;
    } R;
    /** 输出 nc 文件 */
    struct {
        bool active;
        char *path;
    } O;
    /** 源深度、接收深度、坐标旋转、导数和静默输出选项 */
    struct {
        bool active;
    } Dsrc, Drcv, N, e, s;
    /** 点源深度，单位为 km */
    real_t depsrc;
    /** 规则网格接收点深度，单位为 km */
    real_t deprcv;
} GRT_MODULE_CTRL;

/** 均匀半空间的介质参数 */
typedef struct {
    real_t vp;
    real_t vs;
    real_t rho;
    real_t alpha;
    real_t lambda;
    real_t mu;
} OKADA_MEDIUM_PARAMS;

// 与 Coulomb 中用于避开点源奇异点的偏移量保持一致，单位为 km
static const real_t OKADA_SINGULAR_OFFSET = 1.0e-4;

/** 释放命令行参数结构体及其动态分配的成员
 *
 * @param[in,out] Ctrl  命令行参数结构体
 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    if(Ctrl == NULL) return;
    GRT_SAFE_FREE_PTR(Ctrl->Q.path);
    GRT_SAFE_FREE_PTR(Ctrl->O.path);
    GRT_SAFE_FREE_PTR(Ctrl->X.values);
    GRT_SAFE_FREE_PTR(Ctrl->Y.values);
    grt_finite_fault_free(Ctrl->C.faults);
    grt_finite_fault_free(Ctrl->R.faults);
    free(Ctrl);
}

/** 打印使用说明 */
static void print_help(void){
printf("\n"
"[grt okada] %s\n\n", GRT_VERSION);printf(
"    Compute static displacement and displacement derivatives with the Okada analytical solution.\n"
"    The input medium is a homogeneous elastic half-space. The output displacement is in cm\n"
"    and the displacement derivatives are dimensionless.\n"
"    The output format follows module `static_syn`.\n"
"\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    # Point source on a regular grid\n"
"    grt okada -I<vp>/<vs>/<rho> -S[u]<scale> -Ds<depsrc> -Dr<deprcv> -O<outgrid>\n"
"              [-M<strike>/<dip>[/<rake>]]\n"
"              -X<x1>/<x2>/<dx> -Y<y1>/<y2>/<dy> [-N] [-e] [-s]\n"
"\n"
"    # Point source at arbitrary receiver points\n"
"    grt okada -I<vp>/<vs>/<rho> -S[u]<scale> -Ds<depsrc> -Q<file> -O<outgrid>\n"
"              [-M<strike>/<dip>[/<rake>]] [-N] [-e] [-s]\n"
"\n"
"    # Finite faults in Coulomb format\n"
"    grt okada -I<vp>/<vs>/<rho> -C<path> -O<outgrid>\n"
"              [-X<x1>/<x2>/<dx> -Y<y1>/<y2>/<dy> | -Q<file> | -R<fault>[+i<dL>/<dW>]] [-N] [-e] [-s]\n"
"\n"
"    -I specifies <vp>/<vs>/<rho> for a homogeneous elastic half-space.\n"
"    -C evaluates each Coulomb fault row as one exact rectangular Okada fault.\n"
"    Point-source -Ds is required. Grid receivers require -Dr; -Q/-R provide receiver depths.\n"
"\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -I<vp>/<vs>/<rho>\n"
"                  Homogeneous half-space parameters. vp and vs are in km/s; rho is in g/cm^3.\n"
"\n"
"    -S[u]<scale> Point-source scale factor. Without u, <scale> is in dyne-cm.\n"
"                  With u, <scale> is potency or area-slip in cm^3 and is multiplied by source moduli\n"
"                  where required by the point-source type.\n"
"\n"
"    -Ds<depsrc>  Point-source depth in km.\n"
"\n"
"    -Dr<deprcv>  Regular-grid receiver depth in km. Forbidden with -Q/-R.\n"
"\n"
"    -M<strike>/<dip>[/<rake>]\n"
"                  Point-source geometry in degrees. Without rake, use a tensile source; with rake,\n"
"                  use a double-couple source. If -M is omitted, use an explosive source.\n"
"\n"
"    -C<path>      Coulomb finite-fault input file. Each row is evaluated as one rectangular fault.\n"
"                  Kode 100/200/300 select rectangular shear/tensile components;\n"
"                  Kode 400/500 select point sources at the fault-plane center.\n"
"                  Positive Coulomb right-lateral slip is converted to negative Okada strike-slip\n"
"\n"
"    -X<x1>/<x2>/<dx>\n"
"                  Set equidistant receiver points in the north direction, in km.\n"
"                  Mutually exclusive with -Q/-R.\n"
"\n"
"    -Y<y1>/<y2>/<dy>\n"
"                  Set equidistant receiver points in the east direction, in km.\n"
"                  Mutually exclusive with -Q/-R.\n"
"\n"
"    -Q<file>      Arbitrary receiver points from an ASCII file. Each line contains north east depth\n"
"                  in km, optionally followed by strike dip rake in degrees; lines beginning with\n"
"                  # are comments. If present, the three angles are saved as point variables but\n"
"                  are not used in the current calculation. Mutually exclusive with -X/-Y/-Dr/-R.\n"
"\n"
"    -R<fault>[+i<dL>/<dW>]\n"
"                  Coulomb-format finite receiver faults. Without +i, each\n"
"                  fault contributes one point at its rectangular center.\n"
"                  With +i, each fault is subdivided along strike/dip and\n"
"                  the receiver points are the subfault centers. The output\n"
"                  uses one point dimension for all receivers and adds\n"
"                  nfault-dimensional strike/dip/rake/offset/stksize/dipsize variables.\n"
"                  Mutually exclusive with -Q, -X/-Y and -Dr.\n"
"\n"
"    -N            Output components are Z, N, E. Without -N, output is Z, R, T.\n"
"\n"
"    -e            Also output displacement derivatives as nc variables with prefixes z/r/t or z/n/e.\n"
"                  The derivative output can be used by static strain / stress / rotation modules.\n"
"\n"
"    -s            Silence all informational output.\n"
"\n"
"    -h            Display this help message.\n"
"\n"
"\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    Explosion point source on a north/east grid:\n"
"        grt okada -I6/3.464/2.7 -Su1e12 -Ds50 -Dr0 -X-5/5/0.5 -Y-5/5/0.5 -Ookada_ex.nc\n"
"\n"
"    Double-couple and tensile point sources:\n"
"        grt okada -I6/3.464/2.7 -Su1e16 -Ds10 -Dr0 -M100/20/80 -N -Ookada_dc.nc\n"
"        grt okada -I6/3.464/2.7 -Su1e16 -Ds10 -Dr0 -M100/20 -N -Ookada_ts.nc\n"
"\n"
"    Arbitrary receiver points:\n"
"        grt okada -I6/3.464/2.7 -Su1e16 -Ds10 -Qrcv.txt -N -Ookada_q.nc\n"
"\n"
"    Coulomb finite faults and derivatives:\n"
"        grt okada -I6/3.464/2.7 -Cfaults.inp -Dr0 -X-5/5/0.5 -Y-5/5/0.5 -e -Ookada_ff.nc\n"
"\n\n\n"
"\n"
);
}

/** 解析规则接收点坐标轴，并生成等间隔坐标数组
 *
 * @param[in]   text    x1/x2/dx 格式的坐标范围
 * @param[in]   option  当前命令行选项字符
 * @param[out]  n       坐标数组长度
 * @param[out]  values  动态分配的坐标数组
 */
static void parse_axis(const char *text, char option, size_t *n, real_t **values)
{
    real_t a1, a2, delta;
    if((sscanf(text, "%lf/%lf/%lf", &a1, &a2, &delta) != 3) || (delta <= 0.0) || (a1 > a2)){
        GRTRaiseError("Error in \"-%c\". expected x1/x2/dx with x1 <= x2 and dx > 0. Use \"-h\" for help.\n", option);
    }
    *n = (size_t)floor((a2 - a1) / delta) + 1;
    *values = (real_t *)calloc(*n, sizeof(real_t));
    if(*values == NULL) GRTRaiseError("failed to allocate receiver axis.");
    for(size_t i = 0; i < *n; ++i) (*values)[i] = a1 + i * delta;
}

/** 读取并检查 Okada 模块的命令行参数
 *
 * @param[out]  Ctrl    保存解析结果的参数结构体
 * @param[in]   argc    命令行参数数量
 * @param[in]   argv    命令行参数数组
 */
static void parse_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    int opt;
    while((opt = getopt(argc, argv, ":I:O:S:M:C:X:Y:D:Q:R:Nesh")) != -1){
        switch(opt){
            // 读取均匀半空间介质参数
            case 'I': {
                char extra;
                if(sscanf(optarg, "%lf/%lf/%lf%c", &Ctrl->I.vp, &Ctrl->I.vs, &Ctrl->I.rho, &extra) != 3){
                    GRTBadOptionError(I, "expected vp/vs/rho.");
                }
                Ctrl->I.active = true;
                break;
            }
            // 设置输出 nc 文件
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.path = strdup(optarg);
                break;
            // 设置点源放大系数
            case 'S': {
                Ctrl->S.active = true;
                char *value = optarg;
                if(value[0] == 'u'){
                    Ctrl->S.mult_src_mu = true;
                    value++;
                }
                if(*value == '\0' || sscanf(value, "%lf", &Ctrl->S.value) != 1){
                    GRTBadOptionError(S, "");
                }
                break;
            }
            // 设置点源走向、倾角和滑动角
            case 'M': {
                int nscan = sscanf(optarg, "%lf/%lf/%lf", &Ctrl->M.strike, &Ctrl->M.dip, &Ctrl->M.rake);
                if(nscan != 2 && nscan != 3) GRTBadOptionError(M, "expected strike/dip[/rake].");
                Ctrl->M.active = true;
                Ctrl->M.has_rake = nscan == 3;
                if(Ctrl->M.strike < 0.0 || Ctrl->M.strike > 360.0){
                    GRTBadOptionError(M, "strike must be in [0, 360].");
                }
                if(Ctrl->M.dip < 0.0 || Ctrl->M.dip > 90.0){
                    GRTBadOptionError(M, "dip must be in [0, 90].");
                }
                if(Ctrl->M.has_rake && (Ctrl->M.rake < -180.0 || Ctrl->M.rake > 180.0)){
                    GRTBadOptionError(M, "rake must be in [-180, 180].");
                }
                break;
            }
            // 读取 Coulomb 格式有限断层
            case 'C':
                Ctrl->C.active = true;
                grt_finite_fault_free(Ctrl->C.faults);
                Ctrl->C.faults = grt_finite_fault_from_option(
                    optarg, &Ctrl->C.nfault, &Ctrl->C.dL, &Ctrl->C.dW);
                if(Ctrl->C.dL > 0.0){
                    GRTBadOptionError(C, "subdivision suffix is not used by the direct rectangular Okada solution.");
                }
                break;
            // 设置 North 方向规则坐标轴
            case 'X':
                Ctrl->X.active = true;
                parse_axis(optarg, 'X', &Ctrl->X.n, &Ctrl->X.values);
                break;
            // 设置 East 方向规则坐标轴
            case 'Y':
                Ctrl->Y.active = true;
                parse_axis(optarg, 'Y', &Ctrl->Y.n, &Ctrl->Y.values);
                break;
            // 读取任意接收点文件
            case 'Q':
                Ctrl->Q.active = true;
                Ctrl->Q.path = strdup(optarg);
                break;
            // 读取 Coulomb 格式有限接收断层
            case 'R': {
                Ctrl->R.active = true;
                grt_finite_fault_free(Ctrl->R.faults);
                Ctrl->R.faults = grt_finite_fault_from_option(
                    optarg, &Ctrl->R.nfault, &Ctrl->R.dL, &Ctrl->R.dW);
                break;
            }
            // 设置点源深度或规则网格接收点深度
            case 'D':
                if(optarg[0] == 's'){
                    Ctrl->Dsrc.active = true;
                    if(sscanf(optarg + 1, "%lf", &Ctrl->depsrc) != 1 || Ctrl->depsrc < 0.0){
                        GRTBadOptionError(Ds, "source depth must be nonnegative.");
                    }
                } else if(optarg[0] == 'r'){
                    Ctrl->Drcv.active = true;
                    if(sscanf(optarg + 1, "%lf", &Ctrl->deprcv) != 1 || Ctrl->deprcv < 0.0){
                        GRTBadOptionError(Dr, "receiver depth must be nonnegative.");
                    }
                } else {
                    GRTBadOptionError(D, "use -Ds<depth> or -Dr<depth>.");
                }
                break;
            // 输出 ZNE 分量
            case 'N': Ctrl->N.active = true; break;
            // 输出位移偏导
            case 'e': Ctrl->e.active = true; break;
            // 静默输出信息
            case 's': Ctrl->s.active = true; break;
            GRT_Common_Options_in_Switch((char)optopt);
        }
    }

    GRTCheckOptionSet(argc > 1);
    if(!Ctrl->I.active) GRTRaiseError("Okada requires -I<vp>/<vs>/<rho>.\n");
    if(!Ctrl->O.active) GRTRaiseError("Okada requires -O<out>.\n");
    if(Ctrl->X.active ^ Ctrl->Y.active){
        GRTRaiseError("-X and -Y must be specified together.\n");
    }
    if((Ctrl->Q.active || Ctrl->R.active) && (Ctrl->X.active || Ctrl->Y.active || Ctrl->Drcv.active)){
        GRTRaiseError("-Q and -R are mutually exclusive with -X/-Y/-Dr.\n");
    }
    if(Ctrl->Q.active && Ctrl->R.active){
        GRTRaiseError("-Q and -R are mutually exclusive.\n");
    }

    bool point = Ctrl->S.active || Ctrl->M.active;
    if(point == Ctrl->C.active){
        GRTRaiseError("Specify either a point source (-S/-M) or a finite fault (-C).\n");
    }
    if(point){
        if(!Ctrl->S.active) GRTRaiseError("Point source requires -S<scale>.\n");
        if(!Ctrl->Dsrc.active) GRTRaiseError("Point source requires -Ds<depth>.\n");
    } else {
        if(Ctrl->Dsrc.active) GRTRaiseError("-Ds is not used for finite faults.\n");
    }
    if((!Ctrl->Q.active) && (!Ctrl->R.active) && ((!Ctrl->X.active) || (!Ctrl->Y.active) || (!Ctrl->Drcv.active))){
        GRTRaiseError("Grid receivers require -X, -Y and -Dr, or use -Q/-R.\n");
    }
}

/**
 * 根据规则网格或接收点文件构造接收点数组
 *
 * @param[in]  Ctrl   Okada 命令行控制结构体
 * @return            新分配的接收点结构体
 */
static GRT_RECV_POINTS *build_receivers(const GRT_MODULE_CTRL *Ctrl)
{
    if(Ctrl->Q.active) return grt_recv_points_from_file(Ctrl->Q.path);
    if(Ctrl->R.active){
        return grt_recv_points_from_faults(
            Ctrl->R.nfault, Ctrl->R.faults, Ctrl->R.dL, Ctrl->R.dW);
    }
    return grt_recv_points_from_grid(Ctrl->X.n, Ctrl->X.values, Ctrl->Y.n, Ctrl->Y.values, Ctrl->deprcv);
}

/** 将 Okada 局部坐标中的位移和偏导转换到 PyGRT ZNE 坐标
 *
 * @param[in]   strike  断层走向，单位为度
 * @param[in]   local_u Okada 局部位移
 * @param[in]   local_d Okada 局部位移偏导
 * @param[out]  world_u PyGRT ZNE 位移
 * @param[out]  world_d PyGRT ZNE 位移偏导
 */
static void local_to_zne(real_t strike, const real_t local_u[3], const real_t local_d[3][3],
    real_t world_u[3], real_t world_d[3][3])
{
    // Q 的行表示世界坐标分量，列表示 Okada 局部分量
    real_t cs = cos(strike * DEG1), ss = sin(strike * DEG1);
    const real_t q[3][3] = {
        {0.0, 0.0, 1.0},
        {cs, ss, 0.0},
        {ss, -cs, 0.0}
    };

    for(int c = 0; c < 3; ++c){
        world_u[c] = 0.0;
        for(int b = 0; b < 3; ++b) world_u[c] += q[c][b] * local_u[b];
    }
    for(int d = 0; d < 3; ++d){
        for(int c = 0; c < 3; ++c){
            world_d[d][c] = 0.0;
            for(int a = 0; a < 3; ++a){
                for(int b = 0; b < 3; ++b) world_d[d][c] += q[d][a] * q[c][b] * local_d[a][b];
            }
        }
    }
}

/** 计算一个点源在全部接收点上的位移和位移偏导
 *
 * @param[in]       Ctrl      命令行参数结构体
 * @param[in]       medium    均匀半空间介质参数
 * @param[in]       npts      接收点数量
 * @param[in]       norths    接收点 North 坐标
 * @param[in]       easts     接收点 East 坐标
 * @param[in]       depths    接收点深度
 * @param[in,out]   syn       累加后的位移
 * @param[in,out]   syn_d     累加后的位移偏导
 */
static void add_point_source(const GRT_MODULE_CTRL *Ctrl, const OKADA_MEDIUM_PARAMS *medium,
    size_t npts, const real_t *norths, const real_t *easts, const real_t *depths,
    real_t (*syn)[3], real_t (*syn_d)[3][3])
{
    // 将 PyGRT 点源参数转换为 Okada 的四类 potency
    real_t strike = (Ctrl->M.active) ? Ctrl->M.strike : 0.0;
    real_t dip = (Ctrl->M.active) ? Ctrl->M.dip : 0.0;
    real_t pot1 = 0.0, pot2 = 0.0, pot3 = 0.0, pot4 = 0.0;
    real_t scalar = (Ctrl->S.mult_src_mu) ? Ctrl->S.value : Ctrl->S.value / medium->mu;

    if(!Ctrl->M.active){
        pot4 = scalar;
    } else if(Ctrl->M.has_rake){
        pot1 = scalar * cos(Ctrl->M.rake * DEG1);
        pot2 = scalar * sin(Ctrl->M.rake * DEG1);
    } else {
        pot3 = (Ctrl->S.mult_src_mu) ? Ctrl->S.value * medium->mu / medium->lambda : Ctrl->S.value / medium->lambda;
    }

    // Okada 点源输出的位移和偏导分别需要乘以 1e-10 和 1e-15
    for(size_t i = 0; i < npts; ++i){
        real_t cs = cos(strike * DEG1), ss = sin(strike * DEG1);
        real_t x = norths[i] * cs + easts[i] * ss;
        real_t y = norths[i] * ss - easts[i] * cs;
        real_t z = -depths[i];
        real_t u[3], up[3][3], uw[3], dw[3][3];
        int iret = grt_okada_dc3d0(medium->alpha, x, y, z, Ctrl->depsrc, dip,
            pot1, pot2, pot3, pot4, u, up);
        if(iret != 0){
            const char *reason = (iret == 1) ? "singular point" :
                ((iret == 2) ? "receiver is above the free surface" : "unknown error");
            x += OKADA_SINGULAR_OFFSET;
            GRTRaiseWarning(
                "Okada point-source evaluation reached %s (return code %d) at receiver %zu/%zu "
                "(north=%.6g km, east=%.6g km, depth=%.6g km); retry after shifting local X "
                "by %.6g km.",
                reason, iret, i + 1, npts, norths[i], easts[i], depths[i], OKADA_SINGULAR_OFFSET);
            iret = grt_okada_dc3d0(medium->alpha, x, y, z, Ctrl->depsrc, dip,
                pot1, pot2, pot3, pot4, u, up);
            if(iret != 0){
                reason = (iret == 1) ? "singular point" :
                    ((iret == 2) ? "receiver is above the free surface" : "unknown error");
                GRTRaiseWarning(
                    "Okada point-source retry still reached %s (return code %d) at receiver %zu/%zu; "
                    "the contribution is set to zero.",
                    reason, iret, i + 1, npts);
            }
        }
        for(int c = 0; c < 3; ++c){
            u[c] *= 1e-10;
            for(int d = 0; d < 3; ++d) up[d][c] *= 1e-15;
        }
        local_to_zne(strike, u, up, uw, dw);
        if(Ctrl->N.active){
            for(int c = 0; c < 3; ++c) syn[i][c] += uw[c];
            if(Ctrl->e.active){
                for(int d = 0; d < 3; ++d) for(int c = 0; c < 3; ++c) syn_d[i][d][c] += dw[d][c];
            }
        } else {
            // N、E 分量对应公共坐标变换中的 X、Y 分量
            // 与 static_syn 使用相同的零震中距约定，避免网格浮点误差造成任意方位角
            real_t dist = hypot(norths[i], easts[i]);
            real_t theta = (GRT_IS_ZERO(dist)) ? 0.0 : atan2(easts[i], norths[i]);
            real_t radius = dist * 1e5;
            grt_rot_zxy2zrt_upar(theta, uw, dw, radius);
            for(int c = 0; c < 3; ++c) syn[i][c] += uw[c];
            if(Ctrl->e.active){
                for(int d = 0; d < 3; ++d) for(int c = 0; c < 3; ++c) syn_d[i][d][c] += dw[d][c];
            }
        }
    }
}

/** 计算 Coulomb 有限断层在全部接收点上的位移和位移偏导
 *
 * @param[in]       Ctrl      命令行参数结构体
 * @param[in]       medium    均匀半空间介质参数
 * @param[in]       npts      接收点数量
 * @param[in]       norths    接收点 North 坐标
 * @param[in]       easts     接收点 East 坐标
 * @param[in]       depths    接收点深度
 * @param[in,out]   syn       累加后的位移
 * @param[in,out]   syn_d     累加后的位移偏导
 */
static void add_finite_faults(const GRT_MODULE_CTRL *Ctrl, const OKADA_MEDIUM_PARAMS *medium,
    size_t npts, const real_t *norths, const real_t *easts, const real_t *depths,
    real_t (*syn)[3], real_t (*syn_d)[3][3])
{
    for(size_t nf = 0; nf < Ctrl->C.nfault; ++nf){
        const FINITE_FAULT *fault = &Ctrl->C.faults[nf];
        real_t strike = fault->strike;
        real_t dip = fault->dip;
        real_t width, length;
        size_t nW, nL;
        FINITE_SUBFAULT center;
        grt_finite_fault_subdiv(fault, 0.0, 0.0, &width, &length, &nW, &nL);
        (void)nW;
        (void)nL;
        grt_finite_fault_subfault(fault, 0.0, 0.0, width, length, 0, 0, &center);
        real_t coss = cos(strike * DEG1);
        real_t sins = sin(strike * DEG1);
        real_t north_mid = 0.5 * (fault->north_begin + fault->north_end);
        real_t east_mid = 0.5 * (fault->east_begin + fault->east_end);
        real_t depsrc = center.depsrc;

        // Coulomb 的 Kode=400/500 点源位于矩形断层面中心，而不是水平投影中心
        real_t point_north = center.north;
        real_t point_east = center.east;

        real_t disl1 = 0.0;
        real_t disl2 = 0.0;
        real_t disl3 = 0.0;
        real_t pot1 = 0.0;
        real_t pot2 = 0.0;
        real_t pot3 = 0.0;
        real_t pot4 = 0.0;
        if(fault->kode == KODE_RTLAT_REVERSE){
            disl1 = -100.0 * fault->right_lateral;
            disl2 = 100.0 * fault->reverse;
        } else if(fault->kode == KODE_RTLAT_TENSILE){
            disl1 = -100.0 * fault->right_lateral;
            disl3 = 100.0 * fault->tensile;
        } else if(fault->kode == KODE_TENSILE_REVERSE){
            disl2 = 100.0 * fault->reverse;
            disl3 = 100.0 * fault->tensile;
        } else if(fault->kode == KODE_POINT_DC){
            pot1 = -1e6 * fault->right_lateral;
            pot2 = 1e6 * fault->reverse;
        } else if(fault->kode == KODE_POINT_TENSILE_INFLATE){
            pot3 = 1e6 * fault->tensile;
            pot4 = 1e6 * fault->inflate;
        } else {
            GRTRaiseError("unsupported Coulomb Kode=%u.", fault->kode);
        }

        // 同一断层的各接收点只写入自身结果，断层之间保持串行累加
        #pragma omp parallel for schedule(guided) default(shared) if(npts > 1)
        for(size_t i = 0; i < npts; ++i){
            real_t u[3], up[3][3], uw[3], dw[3][3];
            int iret;
            if(KODE_IS_FINITE(fault->kode)){
                // 用顶边水平投影作为 DC3D 参考点，沿上倾方向的范围为 [-width, 0]
                real_t dn = norths[i] - north_mid;
                real_t de = easts[i] - east_mid;
                real_t x = dn * coss + de * sins;
                real_t y = dn * sins - de * coss;
                iret = grt_okada_dc3d(medium->alpha, x, y, -depths[i], fault->top, dip,
                    -0.5 * length, 0.5 * length, -width, 0.0, disl1, disl2, disl3, u, up);
            } else {
                real_t dn = norths[i] - point_north;
                real_t de = easts[i] - point_east;
                real_t x = dn * coss + de * sins;
                real_t y = dn * sins - de * coss;
                iret = grt_okada_dc3d0(medium->alpha, x, y, -depths[i], depsrc, dip,
                    pot1, pot2, pot3, pot4, u, up);
            }
            if(iret != 0){
                const char *reason = (iret == 1) ? "singular point" :
                    ((iret == 2) ? "receiver is above the free surface" : "unknown error");
                GRTRaiseError(
                    "Okada finite-fault evaluation failed: %s (return code %d) at Coulomb fault row %zu/%zu "
                    "(Kode=%u) and receiver %zu/%zu (north=%.6g km, east=%.6g km, depth=%.6g km).",
                    reason, iret, nf + 1, Ctrl->C.nfault, fault->kode,
                    i + 1, npts, norths[i], easts[i], depths[i]);
            }
            if(KODE_IS_FINITE(fault->kode)){
                for(int c = 0; c < 3; ++c){
                    for(int d = 0; d < 3; ++d) up[d][c] *= 1e-5;
                }
            } else {
                // 点源 potency 为 cm^3，水平坐标为 km
                for(int c = 0; c < 3; ++c){
                    u[c] *= 1e-10;
                    for(int d = 0; d < 3; ++d) up[d][c] *= 1e-15;
                }
            }
            local_to_zne(strike, u, up, uw, dw);
            for(int c = 0; c < 3; ++c) syn[i][c] += uw[c];
            if(Ctrl->e.active){
                for(int d = 0; d < 3; ++d) for(int c = 0; c < 3; ++c) syn_d[i][d][c] += dw[d][c];
            }
        }

        if(!Ctrl->s.active){
            GRTRaiseInfo("finite fault[%zu/%zu]", nf + 1, Ctrl->C.nfault);
        }
    }
}

/**
 * 查询均匀介质中的接收介质参数
 *
 * @param[in]   context   均匀介质参数
 * @param[in]   depth     接收点深度 (km)，此处不参与查询
 * @param[out]  va        P 波速度 (km/s)
 * @param[out]  vb        S 波速度 (km/s)
 * @param[out]  rho       密度 (g/cm^3)
 */
static void okada_get_medium(
    void *context, real_t depth, real_t *va, real_t *vb, real_t *rho)
{
    const OKADA_MEDIUM_PARAMS *medium = (const OKADA_MEDIUM_PARAMS *)context;
    (void)depth;
    *va = medium->vp;
    *vb = medium->vs;
    *rho = medium->rho;
}


/**
 * 使用公共静态 NetCDF 输出函数写出 Okada 结果
 *
 * @param[in]  Ctrl         Okada 命令行控制结构体
 * @param[in]  medium       均匀半空间介质参数
 * @param[in]  recv         规则网格、任意点或有限接收断层点列表
 * @param[in]  syn          位移数组
 * @param[in]  syn_d        位移偏导数组
 */
static void save_nc(
    const GRT_MODULE_CTRL *Ctrl, const OKADA_MEDIUM_PARAMS *medium,
    const GRT_RECV_POINTS *recv,
    const real_t (*syn)[3], const real_t (*syn_d)[3][3])
{
    const char *channels;
    if(Ctrl->N.active){
        channels = GRT_ZNE_CODES;
    } else {
        channels = GRT_ZRT_CODES;
    }

    const char *compute_type;
    if(Ctrl->C.active){
        compute_type = "FF";
    } else if((Ctrl->M.active) && (Ctrl->M.has_rake)){
        compute_type = "DC";
    } else if(Ctrl->M.active){
        compute_type = "TS";
    } else {
        compute_type = "EX";
    }

    GRT_STATIC_NC_OUTPUT output = {
        .path = Ctrl->O.path,
        .channels = channels,
        .compute_type = compute_type,
        .coordinate = "Okada X=strike,Y=up-dip horizontal,Z=up",
        .calc_upar = Ctrl->e.active,
        .rot2ZNE = Ctrl->N.active,
        .has_depsrc = (!Ctrl->C.active),
        .depsrc = Ctrl->depsrc,
        .has_elastic_params = true,
        .alpha = medium->alpha,
        .lambda = medium->lambda,
        .mu = medium->mu,
        .recv = recv,
        .get_medium = okada_get_medium,
        .medium_context = (void *)medium,
        .syn = (const real_t (*)[GRT_CHANNEL_NUM])syn,
        .syn_upar = (const real_t (*)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])syn_d,
    };
    grt_static_save_nc(&output);
}
/** Okada 子模块主函数 */
int okada_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = (GRT_MODULE_CTRL *)calloc(1, sizeof(*Ctrl));
    parse_command(Ctrl, argc, argv);
    OKADA_MEDIUM_PARAMS medium = {
        .vp = Ctrl->I.vp,
        .vs = Ctrl->I.vs,
        .rho = Ctrl->I.rho,
    };
    // 速度单位为 km/s，密度单位为 g/cm^3，模量转换为 dyne/cm^2
    medium.alpha = 1.0 - (medium.vs / medium.vp) * (medium.vs / medium.vp);
    medium.mu = medium.rho * medium.vs * medium.vs * 1e10;
    medium.lambda = medium.rho * (medium.vp * medium.vp - 2.0 * medium.vs * medium.vs) * 1e10;
    GRT_RECV_POINTS *recv = build_receivers(Ctrl);
    size_t npts = recv->npts;
    const real_t *norths = recv->norths;
    const real_t *easts = recv->easts;
    const real_t *depths = recv->depths;
    real_t (*syn)[3] = (real_t (*)[3])calloc(npts, sizeof(*syn));
    real_t (*syn_d)[3][3] = (real_t (*)[3][3])calloc(npts, sizeof(*syn_d));
    if((syn == NULL) || (syn_d == NULL)) GRTRaiseError("failed to allocate Okada output.");

    // 按源类型计算位移场
    if(Ctrl->C.active){
        add_finite_faults(Ctrl, &medium, npts, norths, easts, depths, syn, syn_d);
        // 有限断层先在全局 ZNE 中累加，再按 -N 决定最终保存的坐标系
        if(!Ctrl->N.active){
            for(size_t i = 0; i < npts; ++i){
                real_t dist = hypot(norths[i], easts[i]);
                real_t theta = (GRT_IS_ZERO(dist)) ? 0.0 : atan2(easts[i], norths[i]);
                if(Ctrl->e.active){
                    grt_rot_zxy2zrt_upar(theta, syn[i], syn_d[i], dist * 1e5);
                } else {
                    grt_rot_zxy2zrt_vec(theta, syn[i]);
                }
            }
        }
    } else {
        add_point_source(Ctrl, &medium, npts, norths, easts, depths, syn, syn_d);
    }
    save_nc(Ctrl, &medium, recv, syn, syn_d);

    if(!Ctrl->s.active) GRTRaiseInfo("Okada static displacements saved in \"%s\".", Ctrl->O.path);
    GRT_SAFE_FREE_PTR(syn);
    GRT_SAFE_FREE_PTR(syn_d);
    grt_recv_points_free(recv);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
