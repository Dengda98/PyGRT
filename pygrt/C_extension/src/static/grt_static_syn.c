/**
 * @file   grt_static_syn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 *    根据计算好的静态格林函数，定义震源机制以及方位角等，生成合成的静态三分量位移场
 * 
 */

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    /** 输入 nc 格式的格林函数 */
    struct {
        bool active;
        char *s_ingrid;
    } G;
    /** 旋转到 Z, N, E */
    struct {
        bool active;
    } N;
    /** 放大系数 */
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
    /** 对应 Coulomb 程序版本的有限断层 */
    struct {
        bool active;
        real_t dL;   // 沿走向方向的剖分间隔，<=0 表示用默认
        real_t dW;   // 沿倾向方向的剖分间隔，<=0 表示用默认
        size_t nfault;
        FINITE_FAULT *faults;
    } C;
    /** -X: north 坐标 */
    struct {
        bool active;
        size_t nnorth;
        real_t *norths;
    } X;
    /** -Y: east 坐标 */
    struct {
        bool active;
        size_t neast;
        real_t *easts;
    } Y;
    /** -Q: 任意接收点文件 */
    struct {
        bool active;
        char *s_path;
    } Q;
    /** -R: Coulomb 格式有限接收断层 */
    struct {
        bool active;
        real_t dL;
        real_t dW;
        size_t nfault;
        FINITE_FAULT *faults;
    } R;
    /** 输出 nc 文件名 */
    struct {
        bool active;
        char *s_outgrid;
    } O;
    /** 静默输出 */
    struct {
        bool active;
    } s;
    /** 是否计算空间导数 */
    struct {
        bool active;
    } e;
    /** 点源震源深度 / 接收深度：-Ds / -Dr（单值） */
    struct {
        bool s_active;
        bool r_active;
        real_t depsrc;
        real_t deprcv;
    } D;

    // 是否使用 -X/-Y 指定的新接收点网格
    bool isnewNEgrid;

    // 存储不同震源的震源机制相关参数的数组
    real_t mchn[GRT_MECHANISM_NUM];

    // 最终要计算的震源类型
    GRT_SYN_TYPE computeType;
    char s_computeType[3];

    bool isPointSource;
    bool isFiniteFault;

} GRT_MODULE_CTRL;

/**
 * 释放 static_syn 命令行参数结构体及其动态成员
 *
 * @param[in,out]  Ctrl   命令行参数结构体
 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    // G
    GRT_SAFE_FREE_PTR(Ctrl->G.s_ingrid);

    // C
    grt_finite_fault_free(Ctrl->C.faults);
    Ctrl->C.faults = NULL;

    // X
    GRT_SAFE_FREE_PTR(Ctrl->X.norths);

    // Y
    GRT_SAFE_FREE_PTR(Ctrl->Y.easts);

    // Q
    GRT_SAFE_FREE_PTR(Ctrl->Q.s_path);

    // R
    grt_finite_fault_free(Ctrl->R.faults);
    Ctrl->R.faults = NULL;

    // O
    GRT_SAFE_FREE_PTR(Ctrl->O.s_outgrid);

    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt static syn] %s\n\n", GRT_VERSION);printf(
"    Compute static displacement with the outputs of \n"
"    module `static_greenfn` , output to nc file.\n"
"    Three components are:\n"
"       + Up (Z),\n"
"       + Radial Outward (R),\n"
"       + Transverse Clockwise (T),\n"
"    and the units are cm. You can add -N to rotate ZRT to ZNE.\n"
"\n"
"    Receivers: by default reuse the library north/east grid (from\n"
"    static greenfn -X/-Y or -R). Optionally redefine with -X/-Y (uniform\n"
"    depth via -Dr when needed), -Q<file> for arbitrary points, or\n"
"    -R<fault>[+i<dL>/<dW>] for finite receiver faults. -Q/-R are\n"
"    mutually exclusive with -X/-Y and -Dr. If the library was built with -R, the default grid is\n"
"    a 1-D line (north=0, east=R); set -X/-Y or -Q to get a 2-D field.\n"
"    For each receiver, the module first synthesizes results at the\n"
"    surrounding epicentral-distance samples in the library, then combines\n"
"    those results with weights from the requested epicentral distance.\n"
"    When source or receiver depth-based weighting is needed, it performs the\n"
"    same synthesis for each surrounding depth combination and combines those\n"
"    synthesized results with depth-based weights. It combines synthesized\n"
"    results, rather than directly interpolating Green's-function arrays.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    # Point source\n"
"    grt static syn -G<ingrid.nc> -S[u]<scale> -O<outgrid> \n"
"              [-Ds<depsrc>] [-Dr<deprcv>]\n"
"              [-M<strike>/<dip>[/<rake>]]\n"
"              [-T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>]\n"
"              [-F<fn>/<fe>/<fz>] \n"
"              [-X<x1>/<x2>/<dx>] [-Y<y1>/<y2>/<dy>] | [-Q<file>] | [-R<fault>[+i<dL>/<dW>]]\n"
"              [-N] [-e] [-s]\n"
"\n"
"    # Finite faults (Coulomb format)\n"
"    grt static syn -G<ingrid.nc> -C<path>[+i<dL>/<dW>] -O<outgrid>\n"
"              [-Dr<deprcv>] [-X<x1>/<x2>/<dx>] [-Y<y1>/<y2>/<dy>] | [-Q<file>] | [-R<fault>[+i<dL>/<dW>]]\n"
"              [-N] [-e] [-s]\n"
"\n"
"    -G always points to a single 4D STGRNLIB nc file.\n"
"    Depth options (without -Q/-R) depend on the library shape:\n"
"      ndepsrc=1, ndeprcv=1: -Ds/-Dr optional; finite faults forbidden\n"
"      ndepsrc=1, ndeprcv>1: -Dr required; -Ds optional;\n"
"                            finite faults forbidden\n"
"      ndepsrc>1, ndeprcv=1: -Ds required for point source;\n"
"                            -Dr optional; finite faults allowed\n"
"      ndepsrc>1, ndeprcv>1: -Ds required for point source;\n"
"                            -Dr required; finite faults allowed\n"
"    When an optional depth is set, it must be within the library range.\n"
"    With -Q/-R, receiver depths come from the input; do not set -Dr.\n"
"\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<ingrid>    Filepath to a single STGRNLIB nc Green's function\n"
"                  library (dims depsrc×deprcv×north×east).\n"
"\n"
"    -Ds<depsrc>   Point-source source depth (km). Required when the\n"
"                  library has multiple source depths; optional when\n"
"                  it has one source depth. Forbidden for finite faults.\n"
"\n"
"    -Dr<deprcv>   Receiver depth (km) for grid receivers. Required\n"
"                  when the library has multiple receiver depths and -Q/-R\n"
"                  is not used; optional for one depth, but forbidden with -Q/-R.\n"
"\n"
"    -S[u]<scale>  Scale factor to all kinds of point source. \n"
"                  + For Explosion, Shear and Moment Tensor,\n"
"                    unit of <scale> is dyne-cm. \n"
"                  + For Single Force, unit of <scale> is dyne.\n"
"                  + Since \"\\mu\" exists in scalar seismic moment\n"
"                    (\\mu*A*D), you can simply set -Su<scale>, <scale>\n"
"                    equals A*D (Area*Slip, [cm^3]), and <scale> will \n"
"                    multiply \\mu automatically in program.\n"
"\n"
"    -O<outgrid>   Filepath to output nc grid.\n"
"\n"
"    For point source, you can only set at most one of\n"
"    '-M', '-T' and '-F'. If none, an Explosion is used.\n"
"    For finite faults, use '-C' instead (mutually exclusive with\n"
"    point-source options '-S'/'-M'/'-F'/'-T').\n"
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
"    -C<path>[+i<dL>/<dW>]\n"
"                  Finite faults in Coulomb input format.\n"
"                  <path>: Coulomb fault table with 11 numeric columns.\n"
"                  The exact token \"rake\" in the seventh header column\n"
"                  selects the rake/net-slip interpretation; the filename\n"
"                  suffix is not used to select the format.\n"
"                  Kode 100/200/300 are rectangular shear/tensile sources;\n"
"                  Kode 400 is a point double couple; Kode 500 is a point\n"
"                  tensile/inflation source with potency in the two fields.\n"
"                  The rake/net-slip interpretation is restricted to Kode 100.\n"
"                  Optional +i<dL>/<dW>: along-strike / along-dip\n"
"                  subfault size (km). If omitted, both default to\n"
"                  the smallest positive interval among epicentral-distance,\n"
"                  source-depth and receiver-depth sampling in the library.\n"
"                  Each fault: dip in (0, 90], bot > top (km).\n"
"                  Receiver locations default to the library grid;\n"
"                  optional -X/-Y, -Q or -R to redefine.\n"
"                  Requires a library with ndepsrc>1.\n"
"\n"
"    -X<x1>/<x2>/<dx>\n"
"                 Set the equidistant points in the north direction.\n"
"                 <x1>: start coordinate (km).\n"
"                 <x2>: end coordinate (km).\n"
"                 <dx>: sampling interval (km).\n"
"                 Mutually exclusive with -Q/-R.\n"
"\n"
"    -Y<y1>/<y2>/<dy>\n"
"                 Set the equidistant points in the east direction.\n"
"                 <y1>: start coordinate (km).\n"
"                 <y2>: end coordinate (km).\n"
"                 <dy>: sampling interval (km).\n"
"                 Mutually exclusive with -Q/-R.\n"
"\n"
"    -Q<file>      Arbitrary receiver points from an ASCII file.\n"
"                  Each line contains north east depth (km), optionally\n"
"                  followed by strike dip rake (degrees); lines starting\n"
"                  with # are comments. If the three angles are present,\n"
"                  they are saved as point variables but are not used in\n"
"                  synthesis. Mutually exclusive with -X/-Y/-Dr/-R\n"
"                  (depths come from the file).\n"
"\n"
"    -R<fault>[+i<dL>/<dW>]\n"
"                  Coulomb-format finite receiver faults. Without +i, dL=dW\n"
"                  defaults to the smallest positive interval among the\n"
"                  epicentral-distance, source-depth and receiver-depth\n"
"                  sampling in the Green's function library.\n"
"                  With +i, each fault is subdivided along strike/dip and\n"
"                  the receiver points are the subfault centers. The output\n"
"                  uses one point dimension for all receivers and adds\n"
"                  nfault-dimensional strike/dip/rake/offset/stksize/dipsize variables.\n"
"                  -R is mutually exclusive with -Q, -X/-Y and -Dr.\n"
"\n"
"    -N            Components of results will be Z, N, E.\n"
"\n"
"    -e            Also synthesize spatial derivatives of displacement.\n"
"                  Written as nc variables with prefixes z/r/t (ZRT)\n"
"                  or z/n/e (ZNE), e.g. zZ, rR. Required later for\n"
"                  static strain / stress / rotation.\n"
"\n"
"    -s            Silence all outputs.\n"
"\n"
"    -h            Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    2-D north/east grid (syn may omit -X/-Y and reuse the grid):\n"
"        grt static greenfn -Mmilrow -D2/0 -X-5/5/1 -Y-5/5/1 -Ostgrn.nc\n"
"        grt static syn -Gstgrn.nc -Su1e16 -Ostsyn_ex.nc\n"
"\n"
"    Epicentral distances (recommended GF library). Syn should set\n"
"    receivers with -X/-Y or -Q; otherwise output is a line along east:\n"
"        grt static greenfn -Mmilrow -D2/0 -R0/7/0.1 -Ostgrn.nc\n"
"        grt static syn -Gstgrn.nc -Su1e16 -X-5/5/0.5 -Y-5/5/0.5 -Ostsyn_ex.nc\n"
"\n"
"    Other point sources (same -G file):\n"
"        grt static syn -Gstgrn.nc -Su1e16 -M100/20/80 -Ostsyn_dc.nc\n"
"        grt static syn -Gstgrn.nc -Su1e16 -M100/20 -Ostsyn_ts.nc\n"
"        grt static syn -Gstgrn.nc -S1e20 -F0.5/-1.2/3.3 -Ostsyn_sf.nc\n"
"        grt static syn -Gstgrn.nc -Su1e16 -T2.3/0.2/-4.0/0.3/0.5/1.2 -Ostsyn_mt.nc\n"
"\n"
"    Arbitrary receiver points (each line: north east depth in km):\n"
"        grt static syn -Gstgrn.nc -Su1e16 -Qrcv.txt -Ostsyn_q.nc\n"
"\n"
"    Multi-depth library and interpolated source depth:\n"
"        grt static greenfn -Mmilrow -Ds1,2,3 -Dr0 -R0/7/0.1 -Ostgrn.nc\n"
"        grt static syn -Gstgrn.nc -Su1e16 -Ds1.5 -X-5/5/0.5 -Y-5/5/0.5 -Ostsyn.nc\n"
"\n"
"    Finite faults (Coulomb format; library must have ndepsrc>1):\n"
"        grt static syn -Gstgrn.nc -Cfaults.inp+i1/1 -X-5/5/0.5 -Y-5/5/0.5 -Ostsyn_ff.nc\n"
"\n"
"    Spatial derivatives for later strain/stress/rotation:\n"
"        grt static syn -Gstgrn.nc -Su1e16 -e -N -X-5/5/0.5 -Y-5/5/0.5 -Ostsyn.nc\n"
"        grt static strain stsyn.nc\n"
"\n\n\n"
"\n"
);
}


/**
 * 从命令行中读取选项并记录到控制结构体
 *
 * @param[out]  Ctrl   保存解析结果的控制结构体
 * @param[in]   argc   命令行参数数量
 * @param[in]   argv   命令行参数数组
 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    // 先为个别参数设置非0初始值
    Ctrl->computeType = GRT_SYN_EX;
    sprintf(Ctrl->s_computeType, "%s", "EX");

    int opt;
    while ((opt = getopt(argc, argv, ":G:O:S:M:F:T:C:X:Y:D:Q:R:Nesh")) != -1) {
        switch (opt) {
            // 输入 nc 文件名
            case 'G':
                Ctrl->G.active = true;
                Ctrl->G.s_ingrid = strdup(optarg);
                break;

            // -Ds<depsrc> 或 -Dr<deprcv>（单值）
            case 'D':
                if(optarg[0] == 's'){
                    Ctrl->D.s_active = true;
                    if(1 != sscanf(optarg + 1, "%lf", &Ctrl->D.depsrc)){
                        GRTBadOptionError(Ds, "");
                    }
                    if(Ctrl->D.depsrc < 0.0){
                        GRTBadOptionError(Ds, "Negative source depth is not supported.");
                    }
                } else if(optarg[0] == 'r'){
                    Ctrl->D.r_active = true;
                    if(1 != sscanf(optarg + 1, "%lf", &Ctrl->D.deprcv)){
                        GRTBadOptionError(Dr, "");
                    }
                    if(Ctrl->D.deprcv < 0.0){
                        GRTBadOptionError(Dr, "Negative receiver depth is not supported.");
                    }
                } else {
                    GRTBadOptionError(D, "use -Ds<depsrc> or -Dr<deprcv>.");
                }
                break;

            // 输出 nc 文件名
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.s_outgrid = strdup(optarg);
                break;

            // 放大系数
            case 'S':
                Ctrl->S.active = true;
                {   
                    // 检查是否存在字符u，若存在表明需要乘上震源处的剪切模量
                    char *upos=NULL;
                    if((upos=strchr(optarg, 'u')) != NULL){
                        Ctrl->S.mult_src_mu = true;
                        *upos = ' ';
                    }
                }
                if(0 == sscanf(optarg, "%lf", &Ctrl->S.M0)){
                    GRTBadOptionError(S, "");
                };
                break;

            // 剪切震源， 张裂源
            case 'M':
                Ctrl->M.active = true;
                {
                    real_t strike=0.0, dip=0.0, rake=0.0;
                    int nscan = sscanf(optarg, "%lf/%lf/%lf", &strike, &dip, &rake);
                    if(nscan >= 2){
                        Ctrl->computeType = GRT_SYN_TS;
                        sprintf(Ctrl->s_computeType, "%s", "TS");
                        if(strike < 0.0 || strike > 360.0){
                            GRTBadOptionError(M, "Strike must be in [0, 360].");
                        }
                        if(dip < 0.0 || dip > 90.0){
                            GRTBadOptionError(M, "Dip must be in [0, 90].");
                        }
                        if(nscan == 3){
                            Ctrl->computeType = GRT_SYN_DC;
                            sprintf(Ctrl->s_computeType, "%s", "DC");
                            if(rake < -180.0 || rake > 180.0){
                                GRTBadOptionError(M, "Rake must be in [-180, 180].");
                            }
                        }
                    } else {
                        GRTBadOptionError(M, "");
                    };
                    
                    
                    Ctrl->mchn[0] = strike;
                    Ctrl->mchn[1] = dip;
                    Ctrl->mchn[2] = rake;
                }
                break;

            // 单力源
            case 'F':
                Ctrl->F.active = true;
                Ctrl->computeType = GRT_SYN_SF;
                {
                    real_t fn, fe, fz;
                    sprintf(Ctrl->s_computeType, "%s", "SF");
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &fn, &fe, &fz)){
                        GRTBadOptionError(F, "");
                    };
                    Ctrl->mchn[0] = fn;
                    Ctrl->mchn[1] = fe;
                    Ctrl->mchn[2] = fz;
                }
                break;

            // 张量震源
            case 'T':
                Ctrl->T.active = true;
                Ctrl->computeType = GRT_SYN_MT;
                {
                    real_t Mxx, Mxy, Mxz, Myy, Myz, Mzz;
                    sprintf(Ctrl->s_computeType, "%s", "MT");
                    if(6 != sscanf(optarg, "%lf/%lf/%lf/%lf/%lf/%lf", &Mxx, &Mxy, &Mxz, &Myy, &Myz, &Mzz)){
                        GRTBadOptionError(T, "");
                    };
                    Ctrl->mchn[0] = Mxx;
                    Ctrl->mchn[1] = Mxy;
                    Ctrl->mchn[2] = Mxz;
                    Ctrl->mchn[3] = Myy;
                    Ctrl->mchn[4] = Myz;
                    Ctrl->mchn[5] = Mzz;
                }
                break;

            // 从文件中读取有限断层（Coulomb 程序所用格式）-C<path>[+i<dL>/<dW>]
            case 'C':
                Ctrl->C.active = true;
                Ctrl->computeType = GRT_SYN_DC;
                sprintf(Ctrl->s_computeType, "%s", "FF");
                grt_finite_fault_free(Ctrl->C.faults);
                Ctrl->C.faults = grt_finite_fault_from_option(
                    optarg, &Ctrl->C.nfault, &Ctrl->C.dL, &Ctrl->C.dW);
                break;

            // X坐标数组，-Xx1/x2/dx
            case 'X':
                Ctrl->X.active = true;
                {
                    real_t a1, a2, delta;
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &delta)){
                        GRTBadOptionError(X, "");
                    };
                    if(delta <= 0){
                        GRTBadOptionError(X, "Can't set nonpositive dx(%f)", delta);
                    }
                    if(a1 > a2){
                        GRTBadOptionError(X, "x1(%f) > x2(%f).", a1, a2);
                    }

                    Ctrl->X.nnorth = floor((a2-a1)/delta) + 1;
                    Ctrl->X.norths = (real_t*)calloc(Ctrl->X.nnorth, sizeof(real_t));
                    for(size_t i=0; i<Ctrl->X.nnorth; ++i){
                        Ctrl->X.norths[i] = a1 + delta*i;
                    }
                }
                break;

            // Y坐标数组，-Yy1/y2/dy
            case 'Y':
                Ctrl->Y.active = true;
                {
                    real_t a1, a2, delta;
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &delta)){
                        GRTBadOptionError(Y, "");
                    };
                    if(delta <= 0){
                        GRTBadOptionError(Y, "Can't set nonpositive dy(%f)", delta);
                    }
                    if(a1 > a2){
                        GRTBadOptionError(Y, "y1(%f) > y2(%f).", a1, a2);
                    }

                    Ctrl->Y.neast = floor((a2-a1)/delta) + 1;
                    Ctrl->Y.easts = (real_t*)calloc(Ctrl->Y.neast, sizeof(real_t));
                    for(size_t i=0; i<Ctrl->Y.neast; ++i){
                        Ctrl->Y.easts[i] = a1 + delta*i;
                    }
                }
                break;

            // 任意接收点文件，-Q<path>
            case 'Q':
                Ctrl->Q.active = true;
                Ctrl->Q.s_path = strdup(optarg);
                break;

            // 有限接收断层文件，-R<path>[+i<dL>/<dW>]
            case 'R':
                Ctrl->R.active = true;
                grt_finite_fault_free(Ctrl->R.faults);
                Ctrl->R.faults = grt_finite_fault_from_option(
                    optarg, &Ctrl->R.nfault, &Ctrl->R.dL, &Ctrl->R.dW);
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

            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 要么通过 -S, -M/-F/-T 来指定点源，要么通过 -C 来指定有限断层
    bool isPointSource = Ctrl->S.active || Ctrl->M.active || Ctrl->F.active || Ctrl->T.active;
    bool isFiniteFault = Ctrl->C.active;
    if(isPointSource == isFiniteFault){
        GRTRaiseError("You may set either a point source or finite faults — not both, and not neither. Use \"-h\" for help.\n");
    }

    // 检查必选项有没有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, G);
    GRTCheckOptionActive(Ctrl, O);
    if(isPointSource) GRTCheckOptionActive(Ctrl, S);

    // 点源只能使用一种震源
    if(isPointSource && (Ctrl->M.active + Ctrl->F.active + Ctrl->T.active > 1)){
        GRTRaiseError("Only support at most one of \"-M\", \"-F\" and \"-T\". Use \"-h\" for help.\n");
    }

    // -Q/-R 与 -X/-Y 互斥
    if((Ctrl->Q.active || Ctrl->R.active) && (Ctrl->X.active || Ctrl->Y.active)){
        GRTRaiseError("\"-Q\" and \"-R\" are mutually exclusive with \"-X\"/\"-Y\". Use \"-h\" for help.\n");
    }
    if(Ctrl->Q.active && Ctrl->R.active){
        GRTRaiseError("\"-Q\" and \"-R\" are mutually exclusive. Use \"-h\" for help.\n");
    }

    // 指定新接收点网格时必须同时指定 -X 和 -Y
    if(Ctrl->X.active ^ Ctrl->Y.active){
        GRTRaiseError("If you want to set a new north/east grid, you need set \"-X\" and \"-Y\" both.\n");
    }
    Ctrl->isnewNEgrid = Ctrl->X.active;

    // -Q/-R 时深度来自文件，禁止 -Dr
    if((Ctrl->Q.active || Ctrl->R.active) && Ctrl->D.r_active){
        GRTRaiseError("Do not set -Dr with -Q/-R; receiver depths come from the input file.\n");
    }

    Ctrl->isPointSource = isPointSource;
    Ctrl->isFiniteFault = isFiniteFault;
}


/**
 * 在单个震中距采样点上由静态格林函数合成三分量及可选空间偏导
 *
 * @param[in]   azrad       接收点方位角 (rad)
 * @param[in]   ir_pick     震中距采样点索引
 * @param[in]   dist0       震中距采样值 (km)
 * @param[in]   u           位移格林函数
 * @param[in]   uiz         位移对深度的偏导格林函数
 * @param[in]   uir         位移对震中距的偏导格林函数
 * @param[in]   computeType 震源类型
 * @param[in]   M0          标量矩或矩势
 * @param[in]   VpVs_ratio  P 波与 S 波速度比
 * @param[in]   mchn        震源机制参数
 * @param[in]   rot2ZNE     是否输出 ZNE 分量
 * @param[in]   calc_upar   是否计算位移偏导
 * @param[out]  syn         位移结果
 * @param[out]  syn_upar    位移偏导结果
 */
static void static_syn_from_gf_one(
    real_t azrad, size_t ir_pick, real_t dist0,
    const realChnlGrid *u, const realChnlGrid *uiz, const realChnlGrid *uir,
    GRT_SYN_TYPE computeType, real_t M0, real_t VpVs_ratio, const real_t mchn[GRT_MECHANISM_NUM],
    bool rot2ZNE, bool calc_upar,
    real_t syn[GRT_CHANNEL_NUM], real_t syn_upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    int calcUTypes = calc_upar ? 4 : 1;
    realChnlGrid srcRadi = {0};
    real_t tmpsyn[GRT_CHANNEL_NUM];

    for(int ityp = 0; ityp < calcUTypes; ++ityp){
        real_t upar_scale = 1.0;
        // 求位移空间导数时，需调整比例系数（1e-5: km→cm）
        // ZRT 协变导数拆两步：此处合成 (1/r)∂_θ u，后处理再补 ±u/r。
        // r=0 时协变组合仍有限，两项分别换成有限极限（见下方 up 与后处理）。
        if(ityp > 0){
            switch (GRT_ZRT_CODES[ityp-1]){
                case 'Z': case 'R':
                    upar_scale = 1e-5;
                    break;
                // (1/r)∂_θ：r≠0 时 scale∝1/r；r=0 时改用 ∂_r GF，scale 仅留 km→cm
                case 'T':
                    upar_scale = GRT_IS_ZERO(dist0) ? 1e-5 : (1e-5 / dist0);
                    break;
                default:
                    break;
            }
        }

        const realChnlGrid *up = u;
        if(ityp == 1){
            up = uiz;
        } else if(ityp == 2){
            up = uir;
        } else if(ityp == 3 && GRT_IS_ZERO(dist0)){
            // r=0: 用 ∂_r GF（par_θ 辐射因子）合成 (1/r)∂_θ 的有限部分；
            // 后处理中 u/r 联络项改用 ∂_r u，二者合并得有限直角/柱坐标导数。
            // 对 u_z：m≥1 ⇒ u_z(0)=0，lim u_z/r=∂_r u_z；m=0 无 ∂_θ；无 ±u_z/r 联络项。
            up = uir;
        }

        memset(tmpsyn, 0, sizeof(tmpsyn));
        grt_set_source_radiation(srcRadi, computeType, (ityp == 3), M0, upar_scale, VpVs_ratio, azrad, mchn);

        GRT_LOOP_ChnlGrid(im, c){
            int modr = GRT_SRC_M_ORDERS[im];
            if(modr == 0 && GRT_ZRT_CODES[c] == 'T') continue;
            tmpsyn[c] += up[ir_pick][im][c] * srcRadi[im][c];
        }

        for(int i = 0; i < GRT_CHANNEL_NUM; ++i){
            if(ityp == 0){
                syn[i] = tmpsyn[i];
            } else {
                syn_upar[ityp-1][i] = tmpsyn[i];
            }
        }
    }

    if(rot2ZNE){
        if(calc_upar){
            // dist0: km→cm；r=0 时 coord 内 u/r 联络项改用 ∂_r u
            grt_rot_zrt2zxy_upar(azrad, syn, syn_upar, dist0 * 1e5);
        } else {
            grt_rot_zxy2zrt_vec(-azrad, syn);
        }
    }
}


/**
 * 由静态格林函数合成三分量位移场及可选空间偏导
 *
 * GF 侧使用已准备好的升序震中距元数据
 * （sort_rs0 / sort_rs0_idx / isUniform / dr）；查询点为平坦 north/east 列表
 * r=0 时强制方位角为 0（e_r→N、e_θ→E）
 *
 * 数组布局：u[采样点][震源][分量]、syn[接收点][分量]、
 * syn_upar[接收点][偏导方向][分量]。uiz/uir 在 calc_upar=false 时可传 NULL
 *
 * @param[in]      nr0          震中距采样点数
 * @param[in]      sort_rs0     升序震中距采样值
 * @param[in]      sort_rs0_idx 升序采样值对应的原始索引
 * @param[in]      isUniform    震中距采样是否等间隔
 * @param[in]      dr           等间隔震中距步长 (km)
 * @param[in]      npts         接收点数量
 * @param[in]      norths       接收点 North 坐标 (km)
 * @param[in]      easts        接收点 East 坐标 (km)
 * @param[in]      u            位移格林函数
 * @param[in]      uiz          位移对深度的偏导格林函数
 * @param[in]      uir          位移对震中距的偏导格林函数
 * @param[in]      computeType  震源类型
 * @param[in]      M0            标量矩或矩势
 * @param[in]      VpVs_ratio    P 波与 S 波速度比
 * @param[in]      mchn          震源机制参数
 * @param[in]      rot2ZNE       是否输出 ZNE 分量
 * @param[in]      calc_upar     是否计算位移偏导
 * @param[in,out]  syn           位移累加结果
 * @param[in,out]  syn_upar      位移偏导累加结果
 */
static void static_syn_from_gf(
    size_t nr0, const real_t *sort_rs0, const size_t *sort_rs0_idx,
    bool isUniform, real_t dr,
    size_t npts, const real_t *norths, const real_t *easts,
    const realChnlGrid *u, const realChnlGrid *uiz, const realChnlGrid *uir,
    GRT_SYN_TYPE computeType, real_t M0, real_t VpVs_ratio, const real_t mchn[GRT_MECHANISM_NUM],
    bool rot2ZNE, bool calc_upar,
    real_t (*syn)[GRT_CHANNEL_NUM], real_t (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    // 每个接收点逐个处理
    for(size_t ir = 0; ir < npts; ++ir){
        real_t north = norths[ir];
        real_t east = easts[ir];

        real_t dist = hypot(north, east);

        // 方位角；r=0 时 atan2(0,0) 无定义，约定 e_r→N、e_θ→E ⇒ az=0
        real_t azrad = GRT_IS_ZERO(dist) ? 0.0 : atan2(east, north);

        // syn/syn_upar 不在此处清零，以便多次调用（有限断层子源）时对 syn 做累加
        // 调用方需保证首次调用前 syn/syn_upar 已清零（如 calloc）

        // 检查是否越界（允许查询点为精确的零震中距）
        bool r_OutofBound = (dist < sort_rs0[0] - 1e-8 || dist > sort_rs0[nr0-1] + 1e-8);
        if(r_OutofBound){
            GRTRaiseWarning("(north, east)=(%.3e, %.3e) is out of distance bounds, skip.", north, east);
            continue;
        }

        size_t sort_ir_pick = 0, sort_ir_pick1 = 0;
        if(isUniform){
            sort_ir_pick = (size_t)((dist - sort_rs0[0]) / dr);
        } else {
            for(sort_ir_pick = 0; sort_ir_pick < nr0-1; ++sort_ir_pick)   if(sort_rs0[sort_ir_pick+1] > dist)  break;
        }
        sort_ir_pick1 = GRT_MIN(sort_ir_pick + 1, nr0-1);

        // 重复震中距时避免除零（-X/-Y 建库时对称点可能 r 相同）
        real_t r0 = sort_rs0[sort_ir_pick];
        real_t r1 = sort_rs0[sort_ir_pick1];
        real_t drs = (sort_ir_pick == sort_ir_pick1 || fabs(r1 - r0) < 1e-15)
            ? 0.0 : (dist - r0) / (r1 - r0);

        real_t syn2[GRT_CHANNEL_NUM] = {0.0}, syn2_upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM] = {{0.0}};

        size_t iir[2] = {sort_ir_pick, sort_ir_pick1};
        real_t facr[2] = {1.0 - drs, drs};
        for(int j = 0; j < 2; ++j){
            if(j==1 && (sort_ir_pick == sort_ir_pick1 || fabs(r1 - r0) < 1e-15))  continue;

            size_t ir_pick = sort_rs0_idx[iir[j]];
            real_t dist0 = sort_rs0[iir[j]];
            static_syn_from_gf_one(azrad, ir_pick, dist0, u, uiz, uir, computeType, M0, VpVs_ratio, mchn, rot2ZNE, calc_upar, syn2, syn2_upar);

            for(int c = 0; c < GRT_CHANNEL_NUM; ++c){
                syn[ir][c] += facr[j] * syn2[c];
                for(int c2 = 0; c2 < GRT_CHANNEL_NUM; ++c2){
                    syn_upar[ir][c][c2] += facr[j] * syn2_upar[c][c2];
                }
            }
        }
    }
}


/**
 * depsrc×deprcv 双线性角点循环：对给定水平坐标子集做震中距合成，按 fac 累加到 syn
 *
 * 角点权 fac = w_src * w_rcv；结果写入 syn[ipt0 + i]（i = 0..nloc-1）
 *
 * @param[in]      lib              静态格林函数库
 * @param[in]      is0              depsrc 括号左端下标
 * @param[in]      is1              depsrc 括号右端下标
 * @param[in]      ws               depsrc 插值权（落在 is1 侧）
 * @param[in]      na               depsrc 角点个数（1 或 2）
 * @param[in]      ir0              deprcv 括号左端下标
 * @param[in]      ir1              deprcv 括号右端下标
 * @param[in]      wr               deprcv 插值权（落在 ir1 侧）
 * @param[in]      ipt0             写回 syn / syn_upar 的起始下标
 * @param[in]      nloc             本次合成的接收点数
 * @param[in]      loc_n            水平 north 坐标 (km)，长度 nloc
 * @param[in]      loc_e            水平 east 坐标 (km)，长度 nloc
 * @param[in]      computeType      震源类型
 * @param[in]      M0               标量矩或 potency（见 scale_by_src_mu）
 * @param[in]      scale_by_src_mu  为真时用角点 μ 将 M0 转为矩
 * @param[in]      mchn             震源机制参数数组
 * @param[in]      rot2ZNE          是否输出 ZNE
 * @param[in]      calc_upar        是否合成位移空间偏导
 * @param[out]     tmp              单角点位移缓冲，容量 >= nloc
 * @param[out]     tmp_upar         单角点偏导缓冲，容量 >= nloc
 * @param[in,out]  syn              位移累加输出，下标从 ipt0 起
 * @param[in,out]  syn_upar         偏导累加输出，下标从 ipt0 起；仅在 calc_upar 为真时写入
 */
static void static_syn_ps_depth_corners(
    const STGRNLIB *lib,
    size_t is0, size_t is1, real_t ws, int na,
    size_t ir0, size_t ir1, real_t wr,
    size_t ipt0, size_t nloc, const real_t *loc_n, const real_t *loc_e,
    GRT_SYN_TYPE computeType, real_t M0, bool scale_by_src_mu,
    const real_t mchn[GRT_MECHANISM_NUM],
    bool rot2ZNE, bool calc_upar,
    real_t (*tmp)[GRT_CHANNEL_NUM],
    real_t (*tmp_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    real_t (*syn)[GRT_CHANNEL_NUM],
    real_t (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    size_t is_idx[2] = {is0, is1};
    size_t ir_idx[2] = {ir0, ir1};
    int nb = (ir0 == ir1) ? 1 : 2;

    for(int a = 0; a < na; ++a){
        for(int b = 0; b < nb; ++b){
            real_t fac = ((a == 0) ? (1.0 - ws) : ws) * ((b == 0) ? (1.0 - wr) : wr);
            if(fabs(fac) < 1e-15) continue;

            memset(tmp, 0, nloc * sizeof(real_t) * GRT_CHANNEL_NUM);
            if(calc_upar){
                memset(tmp_upar, 0, nloc * sizeof(real_t) * GRT_CHANNEL_NUM * GRT_CHANNEL_NUM);
            }

            size_t is = is_idx[a];
            size_t ir = ir_idx[b];

            // 角点介质取 depsrcs[is]
            real_t va = lib->src_va[is];
            real_t vb = lib->src_vb[is];
            real_t rho = lib->src_rho[is];
            real_t VpVs_ratio = (vb == 0.0) ? 0.0 : (va / vb);
            real_t M0_use = M0;
            if(scale_by_src_mu){
                M0_use = M0 * (vb * vb * rho * 1e10); // dyne/cm^2 * potency
            }

            static_syn_from_gf(
                lib->nr, lib->sort_rs, lib->sort_rs_idx, lib->isUniform, lib->dr,
                nloc, loc_n, loc_e,
                lib->u[is][ir],
                calc_upar ? lib->uiz[is][ir] : NULL,
                calc_upar ? lib->uir[is][ir] : NULL,
                computeType, M0_use, VpVs_ratio, mchn,
                rot2ZNE, calc_upar,
                tmp, tmp_upar
            );

            for(size_t i = 0; i < nloc; ++i){
                size_t ipt = ipt0 + i;
                for(int c = 0; c < GRT_CHANNEL_NUM; ++c){
                    syn[ipt][c] += fac * tmp[i][c];
                    if(calc_upar){
                        for(int c2 = 0; c2 < GRT_CHANNEL_NUM; ++c2){
                            syn_upar[ipt][c][c2] += fac * tmp_upar[i][c][c2];
                        }
                    }
                }
            }
        }
    }
}


/**
 * 基于 STGRNLIB 对各 depsrc×deprcv 邻点做震中距合成，再对结果做二维组合
 *
 * shared_depth 为真（-X/-Y 网格或延用库水平网格）：全部接收点共面，用 depths[0] 求一次 deprcv 括号后批量合成
 * shared_depth 为假（-Q 任意点）：逐点求 deprcv 括号并合成，不做深度归组
 * 各角点的 Vp/Vs（及可选 μ）取自对应 depsrcs 采样；
 * scale_by_src_mu 为真时，M0 为 potency，角点矩为 M0 * μ[is]
 *
 * 输出 syn / syn_upar 按原 npts 下标累加，调用方需事先清零（如 calloc）
 *
 * @param[in]      lib            静态格林函数库
 * @param[in]      depsrc        点源深度 (km)
 * @param[in]      npts          接收点数量
 * @param[in]      norths        接收点 North 坐标 (km)
 * @param[in]      easts         接收点 East 坐标 (km)
 * @param[in]      depths        接收点深度 (km)
 * @param[in]      shared_depth  是否所有接收点共面
 * @param[in]      computeType   震源类型
 * @param[in]      M0            标量矩或矩势
 * @param[in]      scale_by_src_mu 是否使用震源处剪切模量缩放矩势
 * @param[in]      mchn          震源机制参数
 * @param[in]      rot2ZNE       是否输出 ZNE 分量
 * @param[in]      calc_upar     是否计算位移偏导
 * @param[in,out]  syn           位移累加结果
 * @param[in,out]  syn_upar      位移偏导累加结果
 */
static void static_syn_from_gf_PS(
    const STGRNLIB *lib, real_t depsrc,
    size_t npts, const real_t *norths, const real_t *easts, const real_t *depths,
    bool shared_depth,
    GRT_SYN_TYPE computeType, real_t M0, bool scale_by_src_mu,
    const real_t mchn[GRT_MECHANISM_NUM],
    bool rot2ZNE, bool calc_upar,
    real_t (*syn)[GRT_CHANNEL_NUM],
    real_t (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    if(lib == NULL || lib->ndepsrc == 0 || lib->ndeprcv == 0){
        GRTRaiseError("empty STGRNLIB.");
    }
    if(npts == 0 || norths == NULL || easts == NULL || depths == NULL){
        GRTRaiseError("empty receiver points.");
    }
    if(calc_upar && !lib->calc_upar){
        GRTRaiseError("STGRNLIB has no displacement derivatives, cannot set calc_upar.");
    }

    // 震源深度括号：depsrc 落在 depsrcs[is0], depsrcs[is1] 之间，权为 ws
    size_t is0, is1;
    real_t ws;
    if(!grt_locateLinearInterp(lib->depsrcs, lib->ndepsrc, depsrc, &is0, &is1, &ws)){
        GRTRaiseError(
            "Source depth %.6g km is out of Green's function depsrc range [%.6g, %.6g].",
            depsrc, lib->depsrcs[0], lib->depsrcs[lib->ndepsrc - 1]);
    }
    int na = (is0 == is1) ? 1 : 2;  // 落在采样点上时只取一侧

    // 单角点合成缓冲：共面时容量 npts，逐点时容量 1
    size_t nbuf = shared_depth ? npts : 1;
    real_t (*tmp)[GRT_CHANNEL_NUM] =
        (real_t (*)[GRT_CHANNEL_NUM])calloc(nbuf, sizeof(real_t) * GRT_CHANNEL_NUM);
    real_t (*tmp_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM] =
        (real_t (*)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])calloc(
            nbuf, sizeof(real_t) * GRT_CHANNEL_NUM * GRT_CHANNEL_NUM);

    if(shared_depth){
        // 网格：统一深度 depths[0]，一次括号 + 批量震中距合成
        real_t zr = depths[0];
        size_t ir0, ir1;
        real_t wr;
        if(!grt_locateLinearInterp(lib->deprcvs, lib->ndeprcv, zr, &ir0, &ir1, &wr)){
            GRTRaiseError(
                "Receiver depth %.6g km is out of Green's function deprcv range [%.6g, %.6g].",
                zr, lib->deprcvs[0], lib->deprcvs[lib->ndeprcv - 1]);
        }
        static_syn_ps_depth_corners(
            lib, is0, is1, ws, na, ir0, ir1, wr,
            0, npts, norths, easts,
            computeType, M0, scale_by_src_mu, mchn,
            rot2ZNE, calc_upar, tmp, tmp_upar, syn, syn_upar);
    } else {
        // 任意点：逐点括号与合成
        for(size_t ipt = 0; ipt < npts; ++ipt){
            real_t zr = depths[ipt];
            size_t ir0, ir1;
            real_t wr;
            if(!grt_locateLinearInterp(lib->deprcvs, lib->ndeprcv, zr, &ir0, &ir1, &wr)){
                GRTRaiseError(
                    "Receiver depth %.6g km is out of Green's function deprcv range [%.6g, %.6g].",
                    zr, lib->deprcvs[0], lib->deprcvs[lib->ndeprcv - 1]);
            }
            static_syn_ps_depth_corners(
                lib, is0, is1, ws, na, ir0, ir1, wr,
                ipt, 1, &norths[ipt], &easts[ipt],
                computeType, M0, scale_by_src_mu, mchn,
                rot2ZNE, calc_upar, tmp, tmp_upar, syn, syn_upar);
        }
    }

    GRT_SAFE_FREE_PTR(tmp);
    GRT_SAFE_FREE_PTR(tmp_upar);
}


/**
 * 将单个有限断层按 Kode 拆分为 DC、TS、EX 源并累加各子源结果
 *
 * @param[in]      lib           静态格林函数库
 * @param[in]      fault         当前有限断层
 * @param[in]      dL            沿走向子断层尺寸 (km)
 * @param[in]      dW            沿倾向子断层尺寸 (km)
 * @param[in]      W             当前断层沿倾向总长 (km)
 * @param[in]      L             当前断层沿走向总长 (km)
 * @param[in]      nW            倾向子断层数
 * @param[in]      nL            走向子断层数
 * @param[in]      npts          接收点数量
 * @param[in]      norths        接收点 North 坐标 (km)
 * @param[in]      easts         接收点 East 坐标 (km)
 * @param[in]      depths        接收点深度 (km)
 * @param[in]      shared_depth  是否所有接收点共面
 * @param[in]      calc_upar     是否计算位移偏导
 * @param[in,out]  syn           位移累加结果
 * @param[in,out]  syn_upar      位移偏导累加结果
 */
static void static_syn_one_finite_fault(
    const STGRNLIB *lib, const FINITE_FAULT *fault,
    real_t dL, real_t dW, real_t W, real_t L, size_t nW, size_t nL,
    size_t npts, const real_t *norths, const real_t *easts, const real_t *depths,
    bool shared_depth, bool calc_upar,
    real_t (*syn)[GRT_CHANNEL_NUM],
    real_t (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    bool point_source = KODE_IS_POINT(fault->kode);
    size_t loop_nW = point_source ? 1 : nW;
    size_t loop_nL = point_source ? 1 : nL;

    // 按子源并行：各线程累加到私有缓冲，最后归约到 syn，避免对同一接收点写竞争
    size_t nsub = loop_nW * loop_nL;
    #pragma omp parallel default(shared) if(nsub > 1)
    {
        real_t *rcv_norths = (real_t *)calloc(npts, sizeof(real_t));
        real_t *rcv_easts = (real_t *)calloc(npts, sizeof(real_t));
        real_t (*local_syn)[GRT_CHANNEL_NUM] =
            (real_t (*)[GRT_CHANNEL_NUM])calloc(npts, sizeof(real_t) * GRT_CHANNEL_NUM);
        real_t (*local_syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM] =
            (real_t (*)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])calloc(
                npts, sizeof(real_t) * GRT_CHANNEL_NUM * GRT_CHANNEL_NUM);

        #pragma omp for collapse(2) schedule(guided)
        for(size_t iW = 0; iW < loop_nW; ++iW){
            for(size_t iL = 0; iL < loop_nL; ++iL){
                FINITE_SUBFAULT sub;
                if(point_source){
                    // 点源不随 +i 剖分，仍用有限断层面中心确定其位置和深度
                    grt_finite_fault_subfault(fault, L, W, W, L, 0, 0, &sub);
                } else {
                    grt_finite_fault_subfault(fault, dL, dW, W, L, iW, iL, &sub);
                }

                for(size_t ipt = 0; ipt < npts; ++ipt){
                    rcv_norths[ipt] = norths[ipt] - sub.north;
                    rcv_easts[ipt] = easts[ipt] - sub.east;
                }

                real_t mchn[GRT_MECHANISM_NUM] = {0};
                mchn[0] = fault->strike;
                mchn[1] = fault->dip;
                real_t area_scale = sub.width * sub.length * 1e12;

                #define CALL_STATIC_SYN_FROM_GF_PS(source_type, m0) \
                    static_syn_from_gf_PS( \
                        lib, sub.depsrc, \
                        npts, rcv_norths, rcv_easts, depths, \
                        shared_depth, \
                        source_type, m0, true, mchn, \
                        true, calc_upar, \
                        local_syn, local_syn_upar \
                    )

                if(fault->kode == KODE_RTLAT_REVERSE){
                    if(fault->slip != 0.0){
                        mchn[2] = fault->rake;
                        CALL_STATIC_SYN_FROM_GF_PS(GRT_SYN_DC, sub.potency);
                    }
                } else if(fault->kode == KODE_RTLAT_TENSILE){
                    if(fault->right_lateral != 0.0){
                        mchn[2] = 180.0;
                        CALL_STATIC_SYN_FROM_GF_PS(GRT_SYN_DC, fault->right_lateral * area_scale);
                    }
                    if(fault->tensile != 0.0){
                        CALL_STATIC_SYN_FROM_GF_PS(GRT_SYN_TS, fault->tensile * area_scale);
                    }
                } else if(fault->kode == KODE_TENSILE_REVERSE){
                    if(fault->tensile != 0.0){
                        CALL_STATIC_SYN_FROM_GF_PS(GRT_SYN_TS, fault->tensile * area_scale);
                    }
                    if(fault->reverse != 0.0){
                        mchn[2] = 90.0;
                        CALL_STATIC_SYN_FROM_GF_PS(GRT_SYN_DC, fault->reverse * area_scale);
                    }
                } else if(fault->kode == KODE_POINT_DC){
                    real_t potency = hypot(fault->right_lateral, fault->reverse) * 1e6;
                    if(potency != 0.0){
                        mchn[2] = fault->rake;
                        CALL_STATIC_SYN_FROM_GF_PS(GRT_SYN_DC, potency);
                    }
                } else if(fault->kode == KODE_POINT_TENSILE_INFLATE){
                    if(fault->tensile != 0.0){
                        CALL_STATIC_SYN_FROM_GF_PS(GRT_SYN_TS, fault->tensile * 1e6);
                    }
                    if(fault->inflate != 0.0){
                        CALL_STATIC_SYN_FROM_GF_PS(GRT_SYN_EX, fault->inflate * 1e6);
                    }
                } else {
                    GRTRaiseError("unsupported Coulomb Kode=%u.", fault->kode);
                }

                #undef CALL_STATIC_SYN_FROM_GF_PS
            }
        }

        #pragma omp critical(static_syn_ff_reduce)
        {
            for(size_t ipt = 0; ipt < npts; ++ipt){
                for(int c = 0; c < GRT_CHANNEL_NUM; ++c){
                    syn[ipt][c] += local_syn[ipt][c];
                    if(calc_upar){
                        for(int c2 = 0; c2 < GRT_CHANNEL_NUM; ++c2){
                            syn_upar[ipt][c][c2] += local_syn_upar[ipt][c][c2];
                        }
                    }
                }
            }
        }

        GRT_SAFE_FREE_PTR(rcv_norths);
        GRT_SAFE_FREE_PTR(rcv_easts);
        GRT_SAFE_FREE_PTR(local_syn);
        GRT_SAFE_FREE_PTR(local_syn_upar);
    }
}


/**
 * 基于 Coulomb 有限断层和 STGRNLIB 合成静态位移
 * 调用方应先统一确定 dL/dW
 * shared_depth 含义同 static_syn_from_gf_PS
 *
 * @param[in]      lib           静态格林函数库
 * @param[in]      nfault       有限断层数量
 * @param[in]      faults       有限断层数组
 * @param[in]      dL            沿走向子断层尺寸 (km)
 * @param[in]      dW            沿倾向子断层尺寸 (km)
 * @param[in]      npts          接收点数量
 * @param[in]      norths        接收点 North 坐标 (km)
 * @param[in]      easts         接收点 East 坐标 (km)
 * @param[in]      depths        接收点深度 (km)
 * @param[in]      shared_depth  是否所有接收点共面
 * @param[in]      calc_upar     是否计算位移偏导
 * @param[in,out]  syn           位移累加结果
 * @param[in,out]  syn_upar      位移偏导累加结果
 */
static void static_syn_from_gf_FF(
    const STGRNLIB *lib,
    size_t nfault, const FINITE_FAULT *faults,
    real_t dL, real_t dW,
    size_t npts, const real_t *norths, const real_t *easts, const real_t *depths,
    bool shared_depth, bool calc_upar,
    real_t (*syn)[GRT_CHANNEL_NUM],
    real_t (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    if(lib == NULL || lib->ndepsrc == 0 || lib->ndeprcv == 0){
        GRTRaiseError("empty STGRNLIB.");
    }
    if(lib->ndepsrc <= 1){
        GRTRaiseError("Finite faults require a Green's function library with ndepsrc > 1.");
    }
    if(npts == 0 || norths == NULL || easts == NULL || depths == NULL){
        GRTRaiseError("empty receiver points.");
    }
    if((dL <= 0.0) || (dW <= 0.0)){
        GRTRaiseError("finite fault subdivision dL and dW must be positive.");
    }
    if(calc_upar && !lib->calc_upar){
        GRTRaiseError("STGRNLIB has no displacement derivatives, cannot set calc_upar.");
    }

    // 预先检查所有接收深度均在库范围内
    {
        size_t i0, i1;
        real_t w;
        for(size_t ipt = 0; ipt < npts; ++ipt){
            if(!grt_locateLinearInterp(lib->deprcvs, lib->ndeprcv, depths[ipt], &i0, &i1, &w)){
                GRTRaiseError(
                    "Receiver depth %.6g km is out of Green's function deprcv range [%.6g, %.6g].",
                    depths[ipt], lib->deprcvs[0], lib->deprcvs[lib->ndeprcv - 1]);
            }
        }
    }

    for(size_t ifault = 0; ifault < nfault; ++ifault){
        const FINITE_FAULT *f = &faults[ifault];

        real_t W, L;
        size_t nW, nL;
        grt_finite_fault_subdiv(f, dL, dW, &W, &L, &nW, &nL);
        size_t nsub = KODE_IS_POINT(f->kode) ? 1 : nW * nL;
        GRTRaiseInfo("finite fault[%zu/%zu]: nsubfaults = %zu", ifault + 1, nfault, nsub);

        static_syn_one_finite_fault(
            lib, f, dL, dW, W, L, nW, nL,
            npts, norths, easts, depths,
            shared_depth, calc_upar, syn, syn_upar
        );
    }
}


/**
 * 按格林函数库形态校验 -Ds、-Dr 和 -C 选项
 *
 * @param[in]  Ctrl   static_syn 命令行控制结构体
 * @param[in]  lib    静态格林函数库
 */
static void check_syn_depth_options(const GRT_MODULE_CTRL *Ctrl, const STGRNLIB *lib)
{
    bool multi_src = (lib->ndepsrc > 1);
    bool multi_rcv = (lib->ndeprcv > 1);

    if(Ctrl->isFiniteFault && !multi_src){
        GRTRaiseError("Finite faults require a Green's function library with ndepsrc > 1.");
    }
    if(Ctrl->isFiniteFault && Ctrl->D.s_active){
        GRTRaiseError("Do not set -Ds for finite faults; source depths come from the fault geometry.");
    }

    if(Ctrl->Q.active || Ctrl->R.active){
        // -Q/-R：深度来自文件，禁止 -Dr（getopt 已拦一道，此处再保险）
        if(Ctrl->D.r_active){
            GRTRaiseError("Do not set -Dr with -Q/-R; receiver depths come from the input file.");
        }
    } else {
        // 网格接收：多台站深度库必须 -Dr，单台站深度库可省略或显式设置
        if(multi_rcv && !Ctrl->D.r_active){
            GRTRaiseError("Library has multiple receiver depths; -Dr<deprcv> is required.");
        }
    }

    // 点源：多震源深度必须 -Ds，单震源深度可省略或显式设置
    if(Ctrl->isPointSource){
        if(multi_src && !Ctrl->D.s_active){
            GRTRaiseError("Library has multiple source depths; -Ds<depsrc> is required for point source.");
        }
    }
}


/**
 * 统一确定有限震源断层或有限接收断层的剖分尺寸
 *
 * @param[in]      lib      静态格林函数库
 * @param[in,out]  dL       沿走向剖分尺寸 (km)
 * @param[in,out]  dW       沿倾向剖分尺寸 (km)
 * @param[in]      label    输出信息中的断层类型
 */
static void resolve_finite_fault_subdiv(
    const STGRNLIB *lib, real_t *dL, real_t *dW, const char *label)
{
    if((lib == NULL) || (dL == NULL) || (dW == NULL)){
        GRTRaiseError("finite fault subdivision arguments are incomplete.");
    }
    if((*dL <= 0.0) != (*dW <= 0.0)){
        GRTRaiseError("finite fault dL and dW must both be positive or both be omitted.");
    }
    if(*dL <= 0.0){
        // 未指定 +i 时统一使用格林函数库三个采样方向的最小正间隔
        *dL = *dW = grt_stgrnlib_default_subfault_size(lib);
        GRTRaiseInfo("%s: use default dL = dW = %.6g km", label, *dL);
    }
}


/**
 * 构建接收点列表
 *
 * -Q：任意点（各点自有深度，is_grid=false）
 * -R：有限接收断层子断层中心（is_fault=true）
 * 否则：-X/-Y 或延用库水平网格，统一深度（-Dr 或库 deprcvs[0]），is_grid=true
 *
 * @param[in]  Ctrl   static_syn 命令行控制结构体
 * @param[in]  lib    静态格林函数库
 * @return            新分配的接收点结构体
 */
static GRT_RECV_POINTS *build_syn_recv(const GRT_MODULE_CTRL *Ctrl, const STGRNLIB *lib)
{
    if(Ctrl->Q.active){
        return grt_recv_points_from_file(Ctrl->Q.s_path);
    }
    if(Ctrl->R.active){
        return grt_recv_points_from_faults(
            Ctrl->R.nfault, Ctrl->R.faults, Ctrl->R.dL, Ctrl->R.dW);
    }

    size_t nnorth = Ctrl->isnewNEgrid ? Ctrl->X.nnorth : lib->nnorth;
    size_t neast  = Ctrl->isnewNEgrid ? Ctrl->Y.neast  : lib->neast;
    const real_t *norths = Ctrl->isnewNEgrid ? Ctrl->X.norths : lib->norths;
    const real_t *easts  = Ctrl->isnewNEgrid ? Ctrl->Y.easts  : lib->easts;
    real_t deprcv = Ctrl->D.r_active ? Ctrl->D.deprcv : lib->deprcvs[0];
    return grt_recv_points_from_grid(nnorth, norths, neast, easts, deprcv);
}


/**
 * 查询分层介质中的接收介质参数
 *
 * @param[in]   context   静态格林函数库
 * @param[in]   depth     接收点深度 (km)
 * @param[out]  va        P 波速度 (km/s)
 * @param[out]  vb        S 波速度 (km/s)
 * @param[out]  rho       密度 (g/cm^3)
 */
static void static_syn_get_medium(
    void *context, real_t depth, real_t *va, real_t *vb, real_t *rho)
{
    const STGRNLIB *lib = (const STGRNLIB *)context;
    grt_modarr_medium_at_depth(lib->nlayer, lib->modarr, depth, va, vb, rho);
}


/**
 * 使用公共静态 NetCDF 输出函数写出 static_syn 结果
 *
 * @param[in]  path         输出 NetCDF 文件路径
 * @param[in]  Ctrl         static_syn 命令行控制结构体
 * @param[in]  lib          静态格林函数库
 * @param[in]  recv         规则网格、任意点或有限接收断层点列表
 * @param[in]  depsrc       点源深度 (km)
 * @param[in]  syn          位移数组
 * @param[in]  syn_upar     位移偏导数组
 */
static void save_syn_nc(
    const char *path,
    const GRT_MODULE_CTRL *Ctrl,
    const STGRNLIB *lib,
    const GRT_RECV_POINTS *recv,
    real_t depsrc,
    const real_t (*syn)[GRT_CHANNEL_NUM],
    const real_t (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])
{
    const char *channels;
    if(Ctrl->N.active){
        channels = GRT_ZNE_CODES;
    } else {
        channels = GRT_ZRT_CODES;
    }

    GRT_STATIC_NC_OUTPUT output = {
        .path = path,
        .channels = channels,
        .compute_type = Ctrl->s_computeType,
        .calc_upar = Ctrl->e.active,
        .rot2ZNE = Ctrl->N.active,
        .has_depsrc = Ctrl->isPointSource,
        .depsrc = depsrc,
        .has_elastic_params = false,
        .recv = recv,
        .get_medium = static_syn_get_medium,
        .medium_context = (void *)lib,
        .syn = syn,
        .syn_upar = syn_upar,
    };
    grt_static_save_nc(&output);
}
/** 子模块主函数 */
int static_syn_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);

    STGRNLIB *lib = grt_stgrnlib_load_nc(Ctrl->G.s_ingrid);
    check_syn_depth_options(Ctrl, lib);
    if(Ctrl->e.active && !lib->calc_upar){
        GRTRaiseError("Input grid didn't have displacement derivatives, you can't set -e.");
    }

    // -C 和 -R 共用同一套格林函数库剖分间隔规则
    if(Ctrl->C.active){
        resolve_finite_fault_subdiv(lib, &Ctrl->C.dL, &Ctrl->C.dW, "finite fault");
    }
    if(Ctrl->R.active){
        resolve_finite_fault_subdiv(lib, &Ctrl->R.dL, &Ctrl->R.dW, "finite receiver fault");
    }

    // 接收点：网格、-Q 逐点或 -R 有限断层子断层中心
    GRT_RECV_POINTS *recv = build_syn_recv(Ctrl, lib);
    size_t npts = recv->npts;
    const real_t *norths = recv->norths;
    const real_t *easts = recv->easts;
    const real_t *depths = recv->depths;
    bool shared_depth = recv->is_grid;

    real_t (*syn)[GRT_CHANNEL_NUM] = (real_t (*)[GRT_CHANNEL_NUM])calloc(npts, sizeof(real_t) * GRT_CHANNEL_NUM);
    real_t (*syn_upar)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM] =
        (real_t (*)[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM])calloc(npts, sizeof(real_t) * GRT_CHANNEL_NUM * GRT_CHANNEL_NUM);

    // depsrc 仅点源需要；有限断层由各子断层几何提供
    real_t depsrc = 0.0;
    if(Ctrl->isPointSource){
        depsrc = Ctrl->D.s_active ? Ctrl->D.depsrc : lib->depsrcs[0];
        static_syn_from_gf_PS(
            lib, depsrc,
            npts, norths, easts, depths,
            shared_depth,
            Ctrl->computeType, Ctrl->S.M0, Ctrl->S.mult_src_mu, Ctrl->mchn,
            Ctrl->N.active, Ctrl->e.active,
            syn, syn_upar);
    } else {
        static_syn_from_gf_FF(
            lib,
            Ctrl->C.nfault, Ctrl->C.faults,
            Ctrl->C.dL, Ctrl->C.dW,
            npts, norths, easts, depths,
            shared_depth, Ctrl->e.active,
            syn, syn_upar);

        // 有限断层先在全局 ZNE 中累加，再按 -N 决定最终保存的坐标系
        if(!Ctrl->N.active){
            for(size_t i = 0; i < npts; ++i){
                real_t dist = hypot(norths[i], easts[i]);
                real_t theta = (GRT_IS_ZERO(dist)) ? 0.0 : atan2(easts[i], norths[i]);
                if(Ctrl->e.active){
                    grt_rot_zxy2zrt_upar(theta, syn[i], syn_upar[i], dist * 1e5);
                } else {
                    grt_rot_zxy2zrt_vec(theta, syn[i]);
                }
            }
        }
    }

    // 写出 nc（接收介质只在此处按布局查询，不参与合成）
    save_syn_nc(Ctrl->O.s_outgrid, Ctrl, lib, recv, depsrc, syn, syn_upar);

    if(!Ctrl->s.active){
        if(Ctrl->isFiniteFault){
            GRTRaiseInfo("Synthetic static displacements of Coulomb finite faults saved in \"%s\".", Ctrl->O.s_outgrid);
        } else {
            GRTRaiseInfo(
                "Synthetic static displacements of %s source saved in \"%s\".",
                srcTypeFullName[Ctrl->computeType], Ctrl->O.s_outgrid);
        }
    }

    GRT_SAFE_FREE_PTR(syn);
    GRT_SAFE_FREE_PTR(syn_upar);
    grt_recv_points_free(recv);
    grt_stgrnlib_free(lib);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
