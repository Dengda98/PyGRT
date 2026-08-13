/**
 * @file   grt_static_greenfn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 *    计算静态位移
 * 
 */

#include "grt.h"

// 一些变量的非零默认值
#define GRT_GREENFN_K_K0         50.0
#define GRT_GREENFN_L_LENGTH     15.0


/** 该子模块的参数控制结构体 */
typedef struct {
    /** 输入模型 */
    struct {
        bool active;
        char *s_modelpath;        ///< 模型路径
    } M;
    /** 震源和接收器深度：-Dsrc/rcv 或 -Ds/-Dr */
    struct {
        bool active;       ///< 旧式 -D<depsrc>/<deprcv>
        bool s_active;     ///< -Ds
        bool r_active;     ///< -Dr
        size_t ndepsrc;
        real_t *depsrcs;
        size_t ndeprcv;
        real_t *deprcvs;
    } D;
    /** 顶层和底层的边界条件 */
    struct {
        bool active;
        GRT_BOUND_TYPE topbound;
        GRT_BOUND_TYPE botbound;
    } B;
    /** 波数积分间隔以及方法 */
    struct {
        bool active;
        real_t Length;
        real_t kcut;
        struct {
            bool active;
            real_t Length;
        } FIM;
        struct {
            bool active;
            real_t tol;
        } SAFIM;
    } L;
    /** 波数积分收敛方法 */
    struct {
        bool active;
        K_INTEG_CONVERG_METHOD convmet;
    } C;
    /** 波数积分上限 */
    struct {
        bool active;
        real_t keps;
        real_t k0;
        bool use_kmax_ref;
    } K;
    /** 波数积分过程的核函数文件 */
    struct {
        bool active;
        char *s_statsdir;  ///< 保存目录，和当前目录同级
    } S;
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
    /** 输出 nc 文件名 */
    struct {
        bool active;
        char *s_outgrid;
    } O;
    /** 是否计算空间导数 */
    struct {
        bool active;
    } e;

    size_t nr;
    real_t *rs;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    // M
    GRT_SAFE_FREE_PTR(Ctrl->M.s_modelpath);

    // D
    GRT_SAFE_FREE_PTR(Ctrl->D.depsrcs);
    GRT_SAFE_FREE_PTR(Ctrl->D.deprcvs);

    // X
    GRT_SAFE_FREE_PTR(Ctrl->X.norths);

    // Y
    GRT_SAFE_FREE_PTR(Ctrl->Y.easts);

    // O
    GRT_SAFE_FREE_PTR(Ctrl->O.s_outgrid);

    GRT_SAFE_FREE_PTR(Ctrl->rs);

    // S
    GRT_SAFE_FREE_PTR(Ctrl->S.s_statsdir);

    GRT_SAFE_FREE_PTR(Ctrl);
}


/**
 * 打印使用说明
 */
static void print_help(){
printf("\n"
"[grt static greenfn] %s\n\n", GRT_VERSION);printf(
"    Compute static Green's Functions, output to nc file. \n"
"    The units and components are consistent with the dynamics, \n"
"    check \"grt greenfn -h\" for details.\n"
"\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt static greenfn -M<model> -O<outgrid>\n"
"           (-D<depsrc>/<deprcv> | -Ds<source> -Dr<receiver>) \n"
"           [-X<x1>/<x2>/<dx>] [-Y<y1>/<y2>/<dy>] \n"
"           [-R<r1>,<r2>[,...]|<r1>/<r2>/<dr>|<file>]\n" 
"           [-L<length>]   [-C[d|p|n]]  [-Bf|F|r|R|h|H] \n"
"           [-K[+k<k0>][+f][+e<keps>]] [-S]  [-e]\n"
"\n"
"    Output is always one nc file with dims\n"
"    [depsrc][deprcv][north][east] (even if each depth size is 1).\n"
"\n"
"    There're two ways to define the \"epicentral distances\":\n"
"    1. set both -X and -Y (north/east). The kernel still depends\n"
"       only on r=hypot(north,east); results are stored with north/east\n"
"       dims for compatibility.\n"
"    2. set -R for a 1D distance list. Stored as nnorth=1, norths=[0],\n"
"       easts=R (recommended for building a GF library).\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -M<model>    Filepath to 1D horizontally layered halfspace \n"
"                 model. The model file has 6 columns: \n"
"\n"
"         +-------+----------+----------+-------------+----+----+\n"
"         | H(km) | Vp(km/s) | Vs(km/s) | Rho(g/cm^3) | Qp | Qa |\n"
"         +-------+----------+----------+-------------+----+----+\n"
"\n"
"         The 1st column H is the layer thickness. \n"
"         For compatibility, if the \"thickness\" of the first layer \n"
"         set 0.0, then the first column will be the layer top depth, \n"
"         which must be in an ascending order.\n"
"         The number of layers are unlimited.\n"
"\n"
"    -D<depsrc>/<deprcv>\n"
"                 Single source/receiver depth (km). Compatible form.\n"
"                 Mutually exclusive with -Ds/-Dr.\n"
"\n"
"    -Ds<source>  Source depth list (km), same syntax as -R:\n"
"                 + z1,z2[,...]\n"
"                 + z1/z2/dz\n"
"                 + <file>\n"
"                 Must be paired with -Dr.\n"
"\n"
"    -Dr<receiver>\n"
"                 Receiver depth list (km), same syntax as -Ds.\n"
"                 Must be paired with -Ds.\n"
"\n"
"    -X<x1>/<x2>/<dx>\n"
"                 Set the equidistant points in the north direction.\n"
"                 <x1>: start coordinate (km).\n"
"                 <x2>: end coordinate (km).\n"
"                 <dx>: sampling interval (km).\n"
"\n"
"    -Y<y1>/<y2>/<dy>\n"
"                 Set the equidistant points in the east direction.\n"
"                 <y1>: start coordinate (km).\n"
"                 <y2>: end coordinate (km).\n"
"                 <dy>: sampling interval (km).\n"
"\n"
"    -R<r1>,<r2>[,...]|<r1>/<r2>/<dr>|<file>\n"
"                 Multiple epicentral distances (km), support three ways:\n"
"                 + <r1>,<r2>[,...]: seperated by comma (strictly ascending).\n"
"                 + <r1>/<r2>/<dr>:  equal distance <dr> within [r1,r2].\n"
"                 + <file>: each line contains a distance value\n"
"                   (must be strictly ascending).\n"
"\n"
"    -O<outgrid>  Filepath to output nc grid.\n"
"\n"
"    -L[<length>][+l<Flength>][+a<Ftol>][+o<offset>]\n"
"                 Define the wavenumber integration interval\n"
"                 dk=(2*PI)/(<length>*rmax) and methods. \n"
"                 rmax is the maximum epicentral distance. \n"
"                 For DWM:\n"
"                 + (default) not set or set 0.0. \n"
"                   <length> will be %.1f.\n", GRT_GREENFN_L_LENGTH); printf(
"                 + manually set one POSITIVE <length>, e.g. -L20\n"
"                 For FIM or SAFIM (large epicentral distance only;\n"
"                 zero epicentral distance is not allowed):\n"
"                 + +l<Flength> defines the dk of the FIM.\n"
"                 + +a<Ftol> defines the tolerance of the SAFIM.\n"
"                   you can't set both.\n"
"                 + +o<offset> split the integration into two parts,\n"
"                   [0, k*] and [k*, kmax], in which k*=<offset>/rmax,\n"
"                   the former uses DWM and the latter uses FIM/SAFIM.\n"
"\n"
"    -Cd|p|n      Set global convergence method.\n"
"                 + d: Direct Convergence Method (DCM).\n"
"                 + p: Peak-Trough Averaging Method (PTAM).\n"
"                 + n: None.\n"
"                 DCM may still be applied internally when needed.\n"
"\n"
"    -Bf|F|r|R|h|H\n"
"                 Boundary condition of top layer (lowercase) and\n"
"                 bottom layer (uppercase).\n"
"                 f|F: Free boundary.\n"
"                 r|R: Rigid boundary.\n"
"                 h|H: Halfspace.\n"
"\n"
"    -K[+k<k0>][+f][+e<keps>]\n"
"                 Define the reference upper bound of wavenumber integration\n"
"                 kmax_ref = <k0> * PI / hs,\n"
"                 <k0>:   coefficient, default is %.1f, and multiply PI/hs \n", GRT_GREENFN_K_K0); printf(
"                         in program, where hs = max(fabs(depsrc-deprcv), %.1f).\n", GRT_MIN_DEPTH_GAP_SRC_RCV); printf(
"                         The program searches kmax in [dk, kmax_ref] based on\n"
"                         kernel amplitude. If the search reaches kmax_ref without\n"
"                         convergence, or source and receiver are at the same depth,\n"
"                         DCM will be applied (in Auto mode).\n"
"                         If use +f, directly set kmax to kmax_ref.\n"
"                 <keps>: a threshold for break wavenumber \n"
"                         integration in advance. See \n"
"                         (Yao and Harkrider, 1983) for details.\n"
"                         Default 0.0 not use.\n"
"\n"
"    -S           Output statsfile in wavenumber integration.\n"
"                 Only available for a single source/receiver depth.\n"
"\n"
"    -e           Compute the spatial derivatives, ui_z and ui_r,\n"
"                 of displacement u. In columns, prefix \"r\" means \n"
"                 ui_r and \"z\" means ui_z. The units of derivatives\n"
"                 for different sources are: \n"
"                 + Explosion:     1e-25 /(dyne-cm)\n"
"                 + Single Force:  1e-20 /(dyne)\n"
"                 + Shear:         1e-25 /(dyne-cm)\n"
"\n"
"    -h           Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt static greenfn -Mmilrow -D2/0 -X-10/10/1 -Y-10/10/1 -Ostgrn.nc\n"
"    grt static greenfn -Mmilrow -D2/0 -R0/20/1 -Ostgrn.nc\n"
"    grt static greenfn -Mmilrow -Ds0.2/50/0.5 -Dr0 -R0/500/0.5 -Ostgrn.nc -e\n"
"\n\n\n"
);
}


/**
 * 解析深度列表（语法同 -R），结果升序去重写入 *zs / *nz
 *
 * 三种输入形式按优先级依次尝试：
 *   1. 仅含数字与分隔符时：按逗号拆成离散列表 z1,z2,...
 *   2. 可扫成 z1/z2/dz：按等间距生成 [z1, z1+dz, ..., <=z2]
 *   3. 否则当作文件路径：逐行读入数值
 *
 * 随后转为 real_t、禁止负深度、升序排序并按 1e-8 容差去重
 *
 * @param[in]   optarg    -Ds/-Dr 后的字符串（不含前缀 s/r）
 * @param[out]  zs        新分配的深度数组，调用方释放
 * @param[out]  nz        去重后的点数
 * @param[in]   optname   用于报错的选项名（'s' 或 'r'）
 */
static void parse_depth_spec(const char *optarg, real_t **zs, size_t *nz, char optname)
{
    real_t a1, a2, delta;
    char **s_vals = NULL;
    size_t n = 0;

    // 形式 1：逗号分隔列表（字符集与 -R 一致）
    if(grt_string_composed_of(optarg, GRT_NUM_STR "eE+-" ".,")){
        s_vals = grt_string_split(optarg, ",", &n);
    }
    // 形式 2：等间距 z1/z2/dz
    else if(3 == sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &delta)){
        if(delta <= 0){
            GRTRaiseError("-%c: nonpositive spacing (%f).", optname, delta);
        }
        if(a1 > a2){
            GRTRaiseError("-%c: start (%f) > end (%f).", optname, a1, a2);
        }
        n = (size_t)floor((a2 - a1) / delta) + 1;
        s_vals = (char **)calloc(n, sizeof(char *));
        for(size_t i = 0; i < n; ++i){
            GRT_SAFE_ASPRINTF(&s_vals[i], "%.*f", 8, a1 + delta * i);
        }
    }
    // 形式 3：从文件逐行读取
    else {
        FILE *fp = GRTCheckOpenFile(optarg, "r");
        s_vals = grt_string_from_file(fp, &n);
        fclose(fp);
    }

    if(n == 0){
        GRTRaiseError("-%c: empty depth list.", optname);
    }

    // 字符串 -> 数值，并检查非负
    real_t *raw = (real_t *)calloc(n, sizeof(real_t));
    for(size_t i = 0; i < n; ++i){
        raw[i] = atof(s_vals[i]);
        if(raw[i] < 0.0){
            GRTRaiseError("-%c: negative depth (%f) is not supported.", optname, raw[i]);
        }
    }
    GRT_SAFE_FREE_PTR_ARRAY(s_vals, n);

    // 升序排序（argsort 写索引，再按索引取数）
    size_t *order = (size_t *)calloc(n, sizeof(size_t));
    for(size_t i = 0; i < n; ++i) order[i] = i;
    if(n > 1 && grt_argsort(raw, n, sizeof(*raw), grt_compare_real_t, order) != 0){
        GRTRaiseError("-%c: unable to sort depths.", optname);
    }

    // 容差去重，保留升序唯一深度
    real_t *uniq = (real_t *)calloc(n, sizeof(real_t));
    size_t nu = 0;
    for(size_t i = 0; i < n; ++i){
        real_t v = raw[order[i]];
        if(nu == 0 || fabs(v - uniq[nu - 1]) > 1e-8){
            uniq[nu++] = v;
        }
    }
    GRT_SAFE_FREE_PTR(raw);
    GRT_SAFE_FREE_PTR(order);

    *zs = uniq;
    *nz = nu;
}





/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    // 先为个别参数设置非0初始值
    Ctrl->B.topbound = GRT_BOUND_FREE;
    Ctrl->B.botbound = GRT_BOUND_HALFSPACE;

    Ctrl->K.k0 = GRT_GREENFN_K_K0;

    int opt;

    while ((opt = getopt(argc, argv, ":M:D:B:L:C:K:X:Y:R:O:Seh")) != -1) {
        switch (opt) {
            // 模型路径，其中每行分别为 
            //      厚度(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  Qp   Qs
            // 互相用空格隔开即可
            case 'M':
                Ctrl->M.active = true;
                Ctrl->M.s_modelpath = strdup(optarg);
                break;

            // -Dsrc/rcv 或 -Ds<spec> / -Dr<spec>
            case 'D':
                if(optarg[0] == 's'){
                    Ctrl->D.s_active = true;
                    parse_depth_spec(optarg + 1, &Ctrl->D.depsrcs, &Ctrl->D.ndepsrc, 's');
                } else if(optarg[0] == 'r'){
                    Ctrl->D.r_active = true;
                    parse_depth_spec(optarg + 1, &Ctrl->D.deprcvs, &Ctrl->D.ndeprcv, 'r');
                } else {
                    Ctrl->D.active = true;
                    real_t depsrc, deprcv;
                    if(2 != sscanf(optarg, "%lf/%lf", &depsrc, &deprcv)){
                        GRTBadOptionError(D, "");
                    }
                    if(depsrc < 0.0 || deprcv < 0.0){
                        GRTBadOptionError(D, "Negative value in -D is not supported.");
                    }
                    Ctrl->D.ndepsrc = 1;
                    Ctrl->D.ndeprcv = 1;
                    Ctrl->D.depsrcs = (real_t *)calloc(1, sizeof(real_t));
                    Ctrl->D.deprcvs = (real_t *)calloc(1, sizeof(real_t));
                    Ctrl->D.depsrcs[0] = depsrc;
                    Ctrl->D.deprcvs[0] = deprcv;
                }
                break;

            // 顶层和底层的边界条件  -Bf|F|r|R|h|H
            case 'B':
                Ctrl->B.active = true;
                if(strlen(optarg) == 0 || strlen(optarg) > 2)  GRTBadOptionError(B, "");
                for(size_t i = 0; i < strlen(optarg); ++i) {
                    switch(optarg[i]) {
                        case 'f': Ctrl->B.topbound = GRT_BOUND_FREE; break;
                        case 'r': Ctrl->B.topbound = GRT_BOUND_RIGID; break;
                        case 'h': Ctrl->B.topbound = GRT_BOUND_HALFSPACE; break;
                        case 'F': Ctrl->B.botbound = GRT_BOUND_FREE; break;
                        case 'R': Ctrl->B.botbound = GRT_BOUND_RIGID; break;
                        case 'H': Ctrl->B.botbound = GRT_BOUND_HALFSPACE; break;
                        default:
                            GRTBadOptionError(B, "unsupported -B%s.", optarg);
                    }
                }
                break;

            // 波数积分间隔 -L[<length>][+l<Flength>][+a<Ftol>][+o<offset>]
            case 'L':
                Ctrl->L.active = true;
                {
                    char *string = strdup(optarg);
                    char *token = strtok(string, "+");
                    // 如果首先不是加号，则先读取DWM的length
                    if(optarg[0] != '+'){
                        if(1 != sscanf(optarg, "%lf", &Ctrl->L.Length)){
                            GRTBadOptionError(L, "");
                        }
                        token = strtok(NULL, "+");
                    }

                    while(token != NULL){
                        switch(token[0]) {
                            case 'l':
                                Ctrl->L.FIM.active = true;
                                if(1 != sscanf(token+1, "%lf", &Ctrl->L.FIM.Length)){
                                    GRTBadOptionError(L+l, "");
                                }
                                break;
                            
                            case 'a':
                                Ctrl->L.SAFIM.active = true;
                                if(1 != sscanf(token+1, "%lf", &Ctrl->L.SAFIM.tol)){
                                    GRTBadOptionError(L+a, "");
                                }
                                break;
                            
                            case 'o':
                                if(1 != sscanf(token+1, "%lf", &Ctrl->L.kcut)){
                                    GRTBadOptionError(L+o, "");
                                }
                                break;

                            default:
                                GRTBadOptionError(L, "-L+%s is not supported.", token);
                                break;
                        }

                        token = strtok(NULL, "+");
                    }

                    if(Ctrl->L.FIM.active && Ctrl->L.SAFIM.active){
                        GRTBadOptionError(L, "You can't set -L+a and -L+l both.");
                    }

                    GRT_SAFE_FREE_PTR(string);
                }
                break;

            // 波数积分收敛方法  -Cd|p|n
            case 'C':
                Ctrl->C.active = true;
                if(strlen(optarg) == 0){
                    GRTBadOptionError(C, "");
                }
                switch (optarg[0]){
                    case 'p':
                        Ctrl->C.convmet = K_INTEG_CONVERG_PTAM;
                        break;
                    case 'd':
                        Ctrl->C.convmet = K_INTEG_CONVERG_DCM;
                        break;
                    case 'n':
                        Ctrl->C.convmet = K_INTEG_CONVERG_REFUSE;
                        break;
                    default:
                        GRTBadOptionError(C, "-C+%s is not supported.", optarg);
                        break;
                }
                break;

            // 波数积分相关变量 -K[+k<k0>][+f][+e<keps>]
            case 'K':
                Ctrl->K.active = true;
                {
                char *line = strdup(optarg);
                char *token = strtok(line, "+");
                while(token != NULL){
                    switch(token[0]) {
                        case 'k':
                            if(1 != sscanf(token+1, "%lf", &Ctrl->K.k0)){
                                GRTBadOptionError(K+k, "");
                            }
                            if(Ctrl->K.k0 < 0.0){
                                GRTBadOptionError(K, "Can't set negative k0(%f).", Ctrl->K.k0);
                            }
                            break;

                        case 'e':
                            if(1 != sscanf(token+1, "%lf", &Ctrl->K.keps)){
                                GRTBadOptionError(K+e, "");
                            }
                            break;

                        case 'f':
                            Ctrl->K.use_kmax_ref = true;
                            break;

                        default:
                            GRTBadOptionError(K, "-K+%s is not supported.", token);
                            break;
                    }

                    token = strtok(NULL, "+");
                }

                GRT_SAFE_FREE_PTR(line);
                }
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

            // -R：一维震中距序列，存成 nnorth=1, norths=[0], easts=R
            // -R<r1>,<r2>[,...]|<r1>/<r2>/<dr>|<file>
            case 'R':
                Ctrl->X.active = Ctrl->Y.active = true;
                {
                    real_t a1, a2, delta;
                    char **s_easts = NULL;
                    
                    // 如果输入仅由数字、小数点和间隔符组成，则直接读取
                    if(grt_string_composed_of(optarg, GRT_NUM_STR "eE+-" ".,")){
                        s_easts = grt_string_split(optarg, ",", &Ctrl->Y.neast);
                    }
                    // 尝试按照 <r1>/<r2>/<dr> 读取
                    else if(3 == sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &delta)){
                        if(delta <= 0){
                            GRTBadOptionError(R, "Can't set nonpositive dr(%f)", delta);
                        }
                        if(a1 > a2){
                            GRTBadOptionError(R, "r1(%f) > r2(%f).", a1, a2);
                        }

                        Ctrl->Y.neast = floor((a2-a1)/delta) + 1;
                        s_easts = (char **)calloc(Ctrl->Y.neast, sizeof(char*) * Ctrl->Y.neast);
                        for(size_t ir = 0; ir < Ctrl->Y.neast; ++ir){
                            GRT_SAFE_ASPRINTF(&s_easts[ir], "%.*f", 8, a1 + delta*ir);
                        }
                    }
                    // 否则从文件读取
                    else {
                        FILE *fp = GRTCheckOpenFile(optarg, "r");
                        s_easts = grt_string_from_file(fp, &Ctrl->Y.neast);
                        fclose(fp);
                    }

                    // 转为浮点数
                    Ctrl->Y.easts = (real_t*)realloc(Ctrl->Y.easts, sizeof(real_t)*(Ctrl->Y.neast));
                    for(size_t i=0; i < Ctrl->Y.neast; ++i){
                        Ctrl->Y.easts[i] = atof(s_easts[i]);
                        if(Ctrl->Y.easts[i] < 0.0){
                            GRTBadOptionError(R, "Can't set negative epicentral distance(%f).", Ctrl->Y.easts[i]);
                        }
                    }
                    GRT_SAFE_FREE_PTR_ARRAY(s_easts, Ctrl->Y.neast);

                    for(size_t i = 1; i < Ctrl->Y.neast; ++i){
                        if(!(Ctrl->Y.easts[i] > Ctrl->Y.easts[i - 1])){
                            GRTBadOptionError(R, "Epicentral distances must be strictly ascending.");
                        }
                    }

                    Ctrl->X.nnorth = 1;
                    Ctrl->X.norths = (real_t*)calloc(Ctrl->X.nnorth, sizeof(real_t));
                    Ctrl->X.norths[0] = 0.0;
                }
                break;

            // 输出 nc 文件名
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.s_outgrid = strdup(optarg);
                break;

            // 输出波数积分中间文件
            case 'S':
                Ctrl->S.active = true;
                break;

            // 是否计算位移空间导数
            case 'e':
                Ctrl->e.active = true;
                break;
            
            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, M);
    GRTCheckOptionActive(Ctrl, X);
    GRTCheckOptionActive(Ctrl, Y);
    GRTCheckOptionActive(Ctrl, O);

    // 深度选项：-D 与 -Ds/-Dr 互斥；-Ds/-Dr 必须成对
    if(Ctrl->D.active && (Ctrl->D.s_active || Ctrl->D.r_active)){
        GRTRaiseError("Options -D and -Ds/-Dr are mutually exclusive.");
    } else if(Ctrl->D.s_active != Ctrl->D.r_active){
        GRTRaiseError("Options -Ds and -Dr must be set together.");
    } else if(!Ctrl->D.active && !Ctrl->D.s_active){
        GRTRaiseError("Depth option required: -D<depsrc>/<deprcv> or -Ds... -Dr...");
    }

    // 设置震中距数组
    Ctrl->nr = Ctrl->X.nnorth*Ctrl->Y.neast;
    Ctrl->rs = (real_t*)calloc(Ctrl->nr, sizeof(real_t));
    for(size_t inorth=0; inorth<Ctrl->X.nnorth; ++inorth){
        for(size_t ieast=0; ieast<Ctrl->Y.neast; ++ieast){
            Ctrl->rs[ieast + inorth*Ctrl->Y.neast] = hypot(Ctrl->X.norths[inorth], Ctrl->Y.easts[ieast]);
        }
    }

}


/**
 * 静态积分前准备：按震中距与用户参数填充深度无关的 K_INTEG_PROCESS 字段
 * 不写入 Kproc->k0；调用方按 hs 自行缩放
 */
static void prepare_static_grn(
    size_t nr, real_t *rs,
    real_t Length,
    real_t filonLength, real_t safilonTol, real_t filonCut,
    real_t keps, bool use_kmax_ref,
    int convmet,
    K_INTEG_PROCESS *Kproc)
{
    // FIM/SAFIM 面向远震中距，公式含 1/r、1/√r，不能用于 r=0（含 kcut 分段）
    if(filonLength > 0.0 || safilonTol > 0.0){
        for(size_t ir = 0; ir < nr; ++ir){
            if(GRT_IS_ZERO(rs[ir])){
                GRTRaiseError("FIM/SAFIM cannot be used with zero epicentral distance.");
            }
        }
    }

    if(Length == 0.0){
        Length = GRT_GREENFN_L_LENGTH;
    }

    memset(Kproc, 0, sizeof(*Kproc));
    {
        Kproc->use_kmax_ref = use_kmax_ref;
        // 显式收敛方法时不使用 keps
        Kproc->keps = (convmet != K_INTEG_CONVERG_AUTO) ? 0.0 : keps;

        real_t rmax = rs[grt_findMax_real_t(rs, nr)];
        Kproc->kcut = filonCut / rmax;

        // rmax=0 时用阈值防止除零，此时 dk 偏大，后续由 GRT_MIN_NK 收紧
        Kproc->dk = PI2 / (Length * GRT_MAX(rmax, GRT_ZERO_DISTANCE));

        Kproc->applyFIM = filonLength > 0.0;
        Kproc->filondk = (filonLength > 0.0) ? PI2 / (filonLength * rmax) : 0.0;

        Kproc->applySAFIM = safilonTol > 0.0;
        Kproc->sa_tol = safilonTol;

        Kproc->cvgmet = convmet;
    }
}


/** 将积分结果按符号约定写入 STGRNLIB 的一层 */
static void copy_grn_slice_with_sign(
    STGRNLIB *lib, size_t is, size_t ir,
    size_t nr, bool calc_upar,
    const realChnlGrid *grn, const realChnlGrid *grn_uiz, const realChnlGrid *grn_uir)
{
    GRT_LOOP_ChnlGrid(im, c){
        int modr = GRT_SRC_M_ORDERS[im];
        if(modr == 0 && GRT_ZRT_CODES[c] == 'T') continue;
        int sgn0 = (GRT_ZRT_CODES[c] == 'Z') ? -1 : 1;
        for(size_t ipt = 0; ipt < nr; ++ipt){
            lib->u[is][ir][ipt][im][c] = sgn0 * grn[ipt][im][c];
            if(calc_upar){
                lib->uiz[is][ir][ipt][im][c] = (-1) * sgn0 * grn_uiz[ipt][im][c];
                lib->uir[is][ir][ipt][im][c] = sgn0 * grn_uir[ipt][im][c];
            }
        }
    }
}


/**
 * 循环计算多震源/接收深度静态格林函数并写入单个四维 nc
 *
 * Kproc 须已由 prepare_static_grn 填好深度无关字段；
 * 循环内按各 (depsrc, deprcv) 的 hs 更新局部拷贝的 k0
 */
static void compute_stgrnlib_to_nc(
    const char *modelpath,
    size_t ndepsrc, const real_t *depsrcs,
    size_t ndeprcv, const real_t *deprcvs,
    size_t nnorth,  const real_t *norths,
    size_t neast,   const real_t *easts,
    real_t k0,
    K_INTEG_PROCESS *Kproc,
    GRT_BOUND_TYPE topbound, GRT_BOUND_TYPE botbound,
    bool calc_upar,
    const char *outpath,
    const char *statsstr)
{
    if(modelpath == NULL || outpath == NULL || Kproc == NULL
       || depsrcs == NULL || deprcvs == NULL || norths == NULL || easts == NULL){
        GRTRaiseError("NULL argument.");
    }
    if(ndepsrc == 0 || ndeprcv == 0 || nnorth == 0 || neast == 0){
        GRTRaiseError("empty dimension.");
    }
    if(statsstr != NULL && (ndepsrc > 1 || ndeprcv > 1)){
        GRTRaiseError("-S / statsstr is only available for a single source/receiver depth.");
    }

    STGRNLIB *lib = grt_stgrnlib_alloc(
        ndepsrc, depsrcs, ndeprcv, deprcvs, nnorth, norths, neast, easts, calc_upar);
    size_t nr = lib->nr;
    real_t *rs = lib->rs;

    realChnlGrid *grn = (realChnlGrid *)calloc(nr, sizeof(*grn));
    realChnlGrid *grn_uiz = calc_upar ? (realChnlGrid *)calloc(nr, sizeof(*grn_uiz)) : NULL;
    realChnlGrid *grn_uir = calc_upar ? (realChnlGrid *)calloc(nr, sizeof(*grn_uir)) : NULL;

    size_t ntot = ndepsrc * ndeprcv;
    size_t idone = 0;
    const char *modelname = grt_get_basename(modelpath);

    for(size_t is = 0; is < ndepsrc; ++is){
        for(size_t ir = 0; ir < ndeprcv; ++ir){
            real_t zs = depsrcs[is];
            real_t zr = deprcvs[ir];

            MODEL1D *mod1d = NULL;
            if((mod1d = grt_read_mod1d_from_file(modelpath, zs, zr, false)) == NULL){
                exit(EXIT_FAILURE);
            }
            grt_set_mod1d_boundary(mod1d, topbound, botbound);

            // 拷贝模板，避免 integ 改写 cvgmet/dk 等影响后续深度
            K_INTEG_PROCESS local_K = *Kproc;
            real_t hs = GRT_MAX(fabs(zs - zr), GRT_MIN_DEPTH_GAP_SRC_RCV);
            local_K.k0 = k0 * PI / hs;

            memset(grn, 0, nr * sizeof(*grn));
            if(calc_upar){
                memset(grn_uiz, 0, nr * sizeof(*grn_uiz));
                memset(grn_uir, 0, nr * sizeof(*grn_uir));
            }

            grt_integ_static_grn(
                mod1d, nr, rs, &local_K,
                calc_upar, grn, grn_uiz, grn_uir, statsstr);

            lib->src_va[is] = mod1d->Va[mod1d->isrc];
            lib->src_vb[is] = mod1d->Vb[mod1d->isrc];
            lib->src_rho[is] = mod1d->Rho[mod1d->isrc];
            lib->rcv_va[ir] = mod1d->Va[mod1d->ircv];
            lib->rcv_vb[ir] = mod1d->Vb[mod1d->ircv];
            lib->rcv_rho[ir] = mod1d->Rho[mod1d->ircv];

            copy_grn_slice_with_sign(lib, is, ir, nr, calc_upar, grn, grn_uiz, grn_uir);

            idone++;
            GRTRaiseInfo("[%zu/%zu] depsrc=%.6g deprcv=%.6g (%s) done.",
                idone, ntot, zs, zr, modelname);

            grt_free_mod1d(mod1d);
        }
    }

    // 将模型矩阵写入库，供后续 syn/应力按深度查层
    {
        size_t nlayer = 0;
        real_t (*modarr)[GRT_MODARR_NCOL] = grt_read_modarr_from_file(modelpath, &nlayer, false);
        grt_stgrnlib_set_modarr(lib, nlayer, (const real_t (*)[GRT_MODARR_NCOL])modarr);
        GRT_SAFE_FREE_PTR(modarr);
    }

    grt_stgrnlib_save_nc(lib, outpath);
    GRTRaiseInfo("Static Green's function library saved in \"%s\".", outpath);

    GRT_SAFE_FREE_PTR(grn);
    GRT_SAFE_FREE_PTR(grn_uiz);
    GRT_SAFE_FREE_PTR(grn_uir);
    grt_stgrnlib_free(lib);
}


/** 子模块主函数 */
int static_greenfn_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    getopt_from_command(Ctrl, argc, argv);

    bool multi_depth = (Ctrl->D.ndepsrc > 1) || (Ctrl->D.ndeprcv > 1);
    if(Ctrl->S.active){
        if(multi_depth){
            GRTRaiseWarning("-S is ignored for multi-depth STGRNLIB computation.");
        } else {
            // 单深度：stgrtstats/<model>_<depsrc>_<deprcv>
            GRT_SAFE_ASPRINTF(&Ctrl->S.s_statsdir, "stgrtstats");
            GRTCheckMakeDir(Ctrl->S.s_statsdir);
            GRT_SAFE_ASPRINTF(
                &Ctrl->S.s_statsdir, "%s/%s_%g_%g",
                Ctrl->S.s_statsdir,
                grt_get_basename(Ctrl->M.s_modelpath),
                Ctrl->D.depsrcs[0], Ctrl->D.deprcvs[0]);
            GRTCheckMakeDir(Ctrl->S.s_statsdir);
        }
    }

    K_INTEG_PROCESS KPROC = {0};
    prepare_static_grn(
        Ctrl->nr, Ctrl->rs,
        Ctrl->L.Length,
        Ctrl->L.FIM.active ? Ctrl->L.FIM.Length : 0.0,
        Ctrl->L.SAFIM.active ? Ctrl->L.SAFIM.tol : 0.0,
        Ctrl->L.kcut,
        Ctrl->K.keps, Ctrl->K.use_kmax_ref,
        Ctrl->C.convmet,
        &KPROC);

    compute_stgrnlib_to_nc(
        Ctrl->M.s_modelpath,
        Ctrl->D.ndepsrc, Ctrl->D.depsrcs,
        Ctrl->D.ndeprcv, Ctrl->D.deprcvs,
        Ctrl->X.nnorth, Ctrl->X.norths,
        Ctrl->Y.neast, Ctrl->Y.easts,
        Ctrl->K.k0,
        &KPROC,
        Ctrl->B.topbound, Ctrl->B.botbound,
        Ctrl->e.active,
        Ctrl->O.s_outgrid,
        Ctrl->S.s_statsdir);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}

