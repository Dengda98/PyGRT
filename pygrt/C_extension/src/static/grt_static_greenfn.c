/**
 * @file   grt_static_greenfn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 *    计算静态位移
 * 
 */


#include "static/static_grn.h"
#include "common/const.h"
#include "common/model.h"
#include "common/integral.h"
#include "common/iostats.h"
#include "common/search.h"
#include "common/util.h"

#include "grt.h"

// 一些变量的非零默认值
#define GRT_GREENFN_V_VMIN_REF    0.1
#define GRT_GREENFN_K_K0          5.0
#define GRT_GREENFN_L_LENGTH     15.0


/** 该子模块的参数控制结构体 */
typedef struct {
    char *name;
    /** 输入模型 */
    struct {
        bool active;
        char *s_modelpath;        ///< 模型路径
        const char *s_modelname;  ///< 模型名称
        PYMODEL1D *pymod;         ///< 模型PYMODEL1D结构体指针
    } M;
    /** 震源和接收器深度 */
    struct {
        bool active;
        MYREAL depsrc;
        MYREAL deprcv;
        char *s_depsrc;
        char *s_deprcv;
    } D;
    /** 波数积分间隔 */
    struct {
        bool active;
        MYREAL Length;
        MYREAL filonLength;
        MYREAL safilonTol;
        MYREAL filonCut;
    } L;
    /** 波数积分上限 */
    struct {
        bool active;
        MYREAL keps;
        MYREAL k0;
    } K;
    /** 参考速度 */
    struct {
        bool active;
        MYREAL vmin_ref;
    } V;
    /** 波数积分过程的核函数文件 */
    struct {
        bool active;
        char *s_statsdir;  ///< 保存目录，和当前目录同级
    } S;
    /** X 坐标 */
    struct {
        bool active;
        MYINT nx;
        MYREAL *xs;
    } X;
    /** Y 坐标 */
    struct {
        bool active;
        MYINT ny;
        MYREAL *ys;
    } Y;
    /** 是否计算空间导数 */
    struct {
        bool active;
    } e;

    MYINT nr;
    MYREAL *rs;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    free(Ctrl->name);

    // M
    free(Ctrl->M.s_modelpath);
    free_pymod(Ctrl->M.pymod);
    
    // D
    free(Ctrl->D.s_depsrc);
    free(Ctrl->D.s_deprcv);

    // X
    free(Ctrl->X.xs);

    // Y
    free(Ctrl->Y.ys);

    free(Ctrl->rs);

    // S
    if(Ctrl->S.active){
        free(Ctrl->S.s_statsdir);
    }

    free(Ctrl);
}


/**
 * 打印使用说明
 */
static void print_help(){
printf("\n"
"[grt static greenfn] %s\n\n", GRT_VERSION);printf(
"    Compute static Green's Functions, output to stdout. \n"
"    The units and components are consistent with the dynamics, \n"
"    check \"grt greenfn -h\" for details.\n"
"\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt static greenfn -M<model> -D<depsrc>/<deprcv> -X<x1>/<x2>/<dx> \n"
"          -Y<y1>/<y2>/<dy>  [-L<length>] [-V<vmin_ref>] \n" 
"          [-K<k0>[/<keps>]] [-S]  [-e]\n"
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
"                 and the number of layers are unlimited.\n"
"\n"
"    -D<depsrc>/<deprcv>\n"
"                 <depsrc>: source depth (km).\n"
"                 <deprcv>: receiver depth (km).\n"
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
"    -L[a]<length>[/<Flength>/<Fcut>]\n"
"                 Define the wavenumber integration interval\n"
"                 dk=(2*PI)/(<length>*rmax). rmax is the maximum \n"
"                 epicentral distance. \n"
"                 There are 3 cases:\n"
"                 + (default) not set or set 0.0.\n"); printf(
"                   <length> will be %.1f.\n", GRT_GREENFN_L_LENGTH); printf(
"                 + manually set one POSITIVE value, e.g. -L20\n"
"                 + manually set three POSITIVE values, \n"
"                   e.g. -L20/5/10, means split the integration \n"
"                   into two parts, [0, k*] and [k*, kmax], \n"
"                   in which k*=<Fcut>/rmax, and use DWM with\n"
"                   <length> and FIM with <Flength>, respectively.\n"
"                 + manually set three POSITIVE values, with -La,\n"
"                   in this case, <Flength> will be <Ftol> for Self-\n"
"                   Adaptive FIM.\n"
"\n"
"    -V<vmin_ref> \n"
"                 (Inherited from the dynamic case, and the numerical\n"
"                 value will not be used in here, except its sign.)\n"
"                 + (default) not set or set 0.0.\n"); printf(
"                   <vmin_ref> will be the minimum velocity\n"
"                   of model, but limited to %.1f. and if the \n", GRT_GREENFN_V_VMIN_REF); printf(
"                   depth gap between source and receiver is \n"
"                   thinner than %.1f km, PTAM will be appled\n", MIN_DEPTH_GAP_SRC_RCV); printf(
"                   automatically.\n"
"                 + manually set POSITIVE value. \n"
"                 + manually set NEGATIVE value, \n"
"                   and PTAM will be appled.\n"
"\n"
"    -K<k0>[/<keps>]\n"
"                 Several parameters designed to define the\n"
"                 behavior in wavenumber integration. The upper\n"
"                 bound is k0,\n"
"                 <k0>:   default is %.1f, and \n", GRT_GREENFN_K_K0); printf(
"                         multiply PI/hs in program, \n"
"                         where hs = max(fabs(depsrc-deprcv), %.1f).\n", MIN_DEPTH_GAP_SRC_RCV); printf(
"                 <keps>: a threshold for break wavenumber \n"
"                         integration in advance. See \n"
"                         (Yao and Harkrider, 1983) for details.\n"
"                         Default 0.0 not use.\n"); printf(
"\n"
"    -S           Output statsfile in wavenumber integration.\n"
"\n"
"    -e           Compute the spatial derivatives, ui_z and ui_r,\n"
"                 of displacement u. In columns, prefix \"r\" means \n"
"                 ui_r and \"z\" means ui_z. The units of derivatives\n"
"                 for different sources are: \n"
"                 + Explosion:     1e-25 /(dyne-cm)\n"
"                 + Single Force:  1e-20 /(dyne)\n"
"                 + Shear:         1e-25 /(dyne-cm)\n"
"                 + Moment Tensor: 1e-25 /(dyne-cm)\n"
"\n"
"    -h           Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt static greenfn -Mmilrow -D2/0 -X-10/10/20 -Y-10/10/20 > grn\n"
"\n\n\n"
);
}





/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    char* command = Ctrl->name;

    // 先为个别参数设置非0初始值
    Ctrl->V.vmin_ref = GRT_GREENFN_V_VMIN_REF;
    Ctrl->K.k0 = GRT_GREENFN_K_K0;

    int opt;
    while ((opt = getopt(argc, argv, ":M:D:L:K:X:Y:V:Seh")) != -1) {
        switch (opt) {
            // 模型路径，其中每行分别为 
            //      厚度(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  Qp   Qs
            // 互相用空格隔开即可
            case 'M':
                Ctrl->M.active = true;
                Ctrl->M.s_modelpath = strdup(optarg);
                Ctrl->M.s_modelname = get_basename(Ctrl->M.s_modelpath);
                break;

            // 震源和场点深度， -Ddepsrc/deprcv
            case 'D':
                Ctrl->D.active = true;
                Ctrl->D.s_depsrc = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                Ctrl->D.s_deprcv = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                if(2 != sscanf(optarg, "%[^/]/%s", Ctrl->D.s_depsrc, Ctrl->D.s_deprcv)){
                    GRTBadOptionError(command, D, "");
                };
                if(1 != sscanf(Ctrl->D.s_depsrc, "%lf", &Ctrl->D.depsrc)){
                    GRTBadOptionError(command, D, "");
                }
                if(1 != sscanf(Ctrl->D.s_deprcv, "%lf", &Ctrl->D.deprcv)){
                    GRTBadOptionError(command, D, "");
                }
                if(Ctrl->D.depsrc < 0.0 || Ctrl->D.deprcv < 0.0){
                    GRTBadOptionError(command, D, "Negative value in -D is not supported.");
                }
                break;

            // 波数积分间隔 -L[a]<length>[/<Flength>/<Fcut>]
            case 'L':
                Ctrl->L.active = true;
                {
                    // 检查首字母是否为a，表明使用自适应Filon积分
                    int pos=0;
                    bool useSAFIM = false;
                    if(optarg[0] == 'a'){
                        pos++;
                        useSAFIM = true;
                    }
                    double filona = 0.0;
                    int n = sscanf(optarg+pos, "%lf/%lf/%lf", &Ctrl->L.Length, &filona, &Ctrl->L.filonCut);
                    if(n != 1 && n != 3){
                        GRTBadOptionError(command, L, "");
                    };
                    if(n == 1 && Ctrl->L.Length <= 0){
                        GRTBadOptionError(command, L, "Length should be positive.");
                    }
                    if(n == 3 && (filona <= 0 || Ctrl->L.filonCut < 0)){
                        GRTBadOptionError(command, L, "Flength/Ftol should be positive, Fcut should be nonnegative.");
                    }
                    if(n == 3){
                        useSAFIM ? (Ctrl->L.safilonTol = filona) : (Ctrl->L.filonLength = filona);
                    }
                }
                break;

            // 参考最小速度 -Vvmin_ref
            case 'V':
                Ctrl->V.active = true;
                if(0 == sscanf(optarg, "%lf", &Ctrl->V.vmin_ref)){
                    GRTBadOptionError(command, V, "");
                };
                break;

            // 波数积分相关变量 -Kk0/keps
            case 'K':
                Ctrl->K.active = true;
                if(0 == sscanf(optarg, "%lf/%lf", &Ctrl->K.k0, &Ctrl->K.keps)){
                    GRTBadOptionError(command, K, "");
                };
                if(Ctrl->K.k0 < 0.0){
                    GRTBadOptionError(command, K, "Can't set negative k0(%f).", Ctrl->K.k0);
                }
                break;

            // X坐标数组，-Xx1/x2/dx
            case 'X':
                Ctrl->X.active = true;
                {
                    MYREAL a1, a2, delta;
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &delta)){
                        GRTBadOptionError(command, X, "");
                    };
                    if(delta <= 0){
                        GRTBadOptionError(command, X, "Can't set nonpositive dx(%f)", delta);
                    }
                    if(a1 > a2){
                        GRTBadOptionError(command, X, "x1(%f) > x2(%f).", a1, a2);
                    }

                    Ctrl->X.nx = floor((a2-a1)/delta) + 1;
                    Ctrl->X.xs = (MYREAL*)calloc(Ctrl->X.nx, sizeof(MYREAL));
                    for(int i=0; i<Ctrl->X.nx; ++i){
                        Ctrl->X.xs[i] = a1 + delta*i;
                    }
                }
                break;

            // Y坐标数组，-Yy1/y2/dy
            case 'Y':
                Ctrl->Y.active = true;
                {
                    MYREAL a1, a2, delta;
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &delta)){
                        GRTBadOptionError(command, Y, "");
                    };
                    if(delta <= 0){
                        GRTBadOptionError(command, Y, "Can't set nonpositive dy(%f)", delta);
                    }
                    if(a1 > a2){
                        GRTBadOptionError(command, Y, "y1(%f) > y2(%f).", a1, a2);
                    }

                    Ctrl->Y.ny = floor((a2-a1)/delta) + 1;
                    Ctrl->Y.ys = (MYREAL*)calloc(Ctrl->Y.ny, sizeof(MYREAL));
                    for(int i=0; i<Ctrl->Y.ny; ++i){
                        Ctrl->Y.ys[i] = a1 + delta*i;
                    }
                }
                break;

            // 输出波数积分中间文件
            case 'S':
                Ctrl->S.active = true;
                break;

            // 是否计算位移空间导数
            case 'e':
                Ctrl->e.active = true;
                break;
            
            GRT_Common_Options_in_Switch(command, (char)(optopt));
        }
    }

    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(command, argc > 1);
    GRTCheckOptionActive(command, Ctrl, M);
    GRTCheckOptionActive(command, Ctrl, D);
    GRTCheckOptionActive(command, Ctrl, X);
    GRTCheckOptionActive(command, Ctrl, Y);

    // 设置震中距数组
    Ctrl->nr = Ctrl->X.nx*Ctrl->Y.ny;
    Ctrl->rs = (MYREAL*)calloc(Ctrl->nr, sizeof(MYREAL));
    for(int iy=0; iy<Ctrl->Y.ny; ++iy){
        for(int ix=0; ix<Ctrl->X.nx; ++ix){
            Ctrl->rs[ix + iy*Ctrl->X.nx] = GRT_MAX(sqrt(GRT_SQUARE(Ctrl->X.xs[ix]) + GRT_SQUARE(Ctrl->Y.ys[iy])), 1e-5);  // 避免0震中距
        }
    }

}



/**
 * 打印各分量的名称
 * 
 * @param[in]   prefix    前缀字符串
 */
static void print_grn_title(const char *prefix){
    for(int i=0; i<SRC_M_NUM; ++i){
        int modr = SRC_M_ORDERS[i];
        char s_title[10+strlen(prefix)];
        for(int c=0; c<CHANNEL_NUM; ++c){
            if(modr==0 && ZRTchs[c]=='T')  continue;

            snprintf(s_title, sizeof(s_title), "%s%s%c", prefix, SRC_M_NAME_ABBR[i], ZRTchs[c]);
            fprintf(stdout, GRT_STRING_FMT, s_title);
        }
    }
}

/**
 * 打印各分量的值
 * 
 * @param      grn       静态格林函数结果
 * @param      sgn0      全局符号
 */
static void print_grn_value(const MYREAL grn[SRC_M_NUM][CHANNEL_NUM], const int sgn0){
    for(int i=0; i<SRC_M_NUM; ++i){
        int modr = SRC_M_ORDERS[i];
        int sgn = 1;
        for(int c=0; c<CHANNEL_NUM; ++c){
            if(modr==0 && ZRTchs[c]=='T')  continue;

            sgn = (ZRTchs[c]=='Z') ? -sgn0 : sgn0;

            fprintf(stdout, GRT_REAL_FMT, sgn * grn[i][c]);
        }
    }
}


/** 子模块主函数 */
int static_greenfn_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    Ctrl->name = strdup(argv[0]);
    
    const char *command = Ctrl->name;

    // 传入参数 
    getopt_from_command(Ctrl, argc, argv);

    // 读入模型文件（暂先不考虑液体层）
    if((Ctrl->M.pymod = read_pymod_from_file(command, Ctrl->M.s_modelpath, Ctrl->D.depsrc, Ctrl->D.deprcv, false)) == NULL){
        exit(EXIT_FAILURE);
    }
    PYMODEL1D *pymod = Ctrl->M.pymod;

    // 最大最小速度
    MYREAL vmin, vmax;
    get_pymod_vmin_vmax(pymod, &vmin, &vmax);

    // 参考最小速度
    if(!Ctrl->V.active){
        Ctrl->V.vmin_ref = GRT_MAX(vmin, GRT_GREENFN_V_VMIN_REF);
    } 

    // 如果没有主动设置vmin_ref，则判断是否要自动使用PTAM
    if( !Ctrl->V.active && fabs(Ctrl->D.deprcv - Ctrl->D.depsrc) <= MIN_DEPTH_GAP_SRC_RCV) {
        Ctrl->V.vmin_ref = - fabs(Ctrl->V.vmin_ref);
    }
    
    // 设置积分间隔默认值
    if(Ctrl->L.Length == 0.0)  Ctrl->L.Length = GRT_GREENFN_L_LENGTH;

    // 波数积分输出目录
    if(Ctrl->S.active){
        Ctrl->S.s_statsdir = (char*)malloc(sizeof(char)*(strlen(Ctrl->M.s_modelpath)+strlen(Ctrl->D.s_depsrc)+strlen(Ctrl->D.s_deprcv)+100));
        sprintf(Ctrl->S.s_statsdir, "stgrtstats");
        // 建立保存目录
        GRTCheckMakeDir(command, Ctrl->S.s_statsdir);
        sprintf(Ctrl->S.s_statsdir, "%s/%s_%s_%s", Ctrl->S.s_statsdir, Ctrl->M.s_modelname, Ctrl->D.s_depsrc, Ctrl->D.s_deprcv);
        GRTCheckMakeDir(command, Ctrl->S.s_statsdir);
    }

    // 建立格林函数的浮点数
    MYREAL (*grn)[SRC_M_NUM][CHANNEL_NUM] = (MYREAL (*)[SRC_M_NUM][CHANNEL_NUM]) calloc(Ctrl->nr, sizeof(*grn));
    MYREAL (*grn_uiz)[SRC_M_NUM][CHANNEL_NUM] = (Ctrl->e.active)? (MYREAL (*)[SRC_M_NUM][CHANNEL_NUM]) calloc(Ctrl->nr, sizeof(*grn_uiz)) : NULL;
    MYREAL (*grn_uir)[SRC_M_NUM][CHANNEL_NUM] = (Ctrl->e.active)? (MYREAL (*)[SRC_M_NUM][CHANNEL_NUM]) calloc(Ctrl->nr, sizeof(*grn_uir)) : NULL;


    //==============================================================================
    // 计算静态格林函数
    integ_static_grn(
        pymod, Ctrl->nr, Ctrl->rs, Ctrl->V.vmin_ref, Ctrl->K.keps, Ctrl->K.k0, Ctrl->L.Length, Ctrl->L.filonLength, Ctrl->L.safilonTol, Ctrl->L.filonCut, 
        grn, Ctrl->e.active, grn_uiz, grn_uir,
        Ctrl->S.s_statsdir
    );
    //==============================================================================

    MYREAL src_va = pymod->Va[pymod->isrc];
    MYREAL src_vb = pymod->Vb[pymod->isrc];
    MYREAL src_rho = pymod->Rho[pymod->isrc];
    MYREAL rcv_va = pymod->Va[pymod->ircv];
    MYREAL rcv_vb = pymod->Vb[pymod->ircv];
    MYREAL rcv_rho = pymod->Rho[pymod->ircv];

    // 输出物性参数
    fprintf(stdout, "# "GRT_REAL_FMT" "GRT_REAL_FMT" "GRT_REAL_FMT"\n", src_va, src_vb, src_rho);
    fprintf(stdout, "# "GRT_REAL_FMT" "GRT_REAL_FMT" "GRT_REAL_FMT"\n", rcv_va, rcv_vb, rcv_rho);


    // 输出标题
    char XX[20];
    sprintf(XX, GRT_STRING_FMT, "X(km)"); XX[0]='#';
    fprintf(stdout, "%s", XX);
    fprintf(stdout, GRT_STRING_FMT, "Y(km)");
    print_grn_title("");

    if(Ctrl->e.active) {
        print_grn_title("z");
        print_grn_title("r");
    }
    fprintf(stdout, "\n");

    // 写结果
    for(int iy=0; iy<Ctrl->Y.ny; ++iy) {
        for(int ix=0; ix<Ctrl->X.nx; ++ix) {
            int ir = ix + iy * Ctrl->X.nx;
            fprintf(stdout, GRT_REAL_FMT GRT_REAL_FMT, Ctrl->X.xs[ix], Ctrl->Y.ys[iy]);

            print_grn_value(grn[ir], 1);

            if(Ctrl->e.active) {
                print_grn_value(grn_uiz[ir], -1);
                print_grn_value(grn_uir[ir], 1);
            }
            fprintf(stdout, "\n");
        }
    }

    // 释放内存
    free(grn);
    if(grn_uiz) free(grn_uiz);
    if(grn_uir) free(grn_uir);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}

