/**
 * @file   grt_static_greenfn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 *    计算静态位移
 * 
 */


#include "grt/static/static_grn.h"
#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/common/integral.h"
#include "grt/common/iostats.h"
#include "grt/common/search.h"
#include "grt/common/util.h"
#include "grt/common/mynetcdf.h"

#include "grt.h"

// 一些变量的非零默认值
#define GRT_GREENFN_K_VMIN        0.1
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
        GRT_MODEL1D *mod1d;         ///< 模型结构体指针
    } M;
    /** 震源和接收器深度 */
    struct {
        bool active;
        real_t depsrc;
        real_t deprcv;
        char *s_depsrc;
        char *s_deprcv;
    } D;
    /** 波数积分间隔 */
    struct {
        bool active;
        enum GRT_K_INTEG_METHOD method;
        real_t Length;
        real_t filonLength;
        real_t safilonTol;
        real_t filonCut;
    } L;
    /** 波数积分上限 */
    struct {
        bool active;
        real_t keps;
        real_t k0;
        real_t vmin;
        bool v_active;
    } K;
    /** 波数积分过程的核函数文件 */
    struct {
        bool active;
        char *s_statsdir;  ///< 保存目录，和当前目录同级
    } S;
    /** X 坐标 */
    struct {
        bool active;
        size_t nx;
        real_t *xs;
    } X;
    /** Y 坐标 */
    struct {
        bool active;
        size_t ny;
        real_t *ys;
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
    GRT_SAFE_FREE_PTR(Ctrl->name);

    // M
    GRT_SAFE_FREE_PTR(Ctrl->M.s_modelpath);
    grt_free_mod1d(Ctrl->M.mod1d);
    
    // D
    GRT_SAFE_FREE_PTR(Ctrl->D.s_depsrc);
    GRT_SAFE_FREE_PTR(Ctrl->D.s_deprcv);

    // X
    GRT_SAFE_FREE_PTR(Ctrl->X.xs);

    // Y
    GRT_SAFE_FREE_PTR(Ctrl->Y.ys);

    // O
    GRT_SAFE_FREE_PTR(Ctrl->O.s_outgrid);

    GRT_SAFE_FREE_PTR(Ctrl->rs);

    // S
    if(Ctrl->S.active){
        GRT_SAFE_FREE_PTR(Ctrl->S.s_statsdir);
    }

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
"    grt static greenfn -M<model> -D<depsrc>/<deprcv> -X<x1>/<x2>/<dx> \n"
"          -Y<y1>/<y2>/<dy>  -O<outgrid>  [-L<length>]  \n" 
"          [-K[+k<k0>][+e<keps>][+v<vmin>]] [-S]  [-e]\n"
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
"    -O<outgrid>  Filepath to output nc grid.\n"
"\n"
"    -L[a|l]<length>[/<Flength>|<Ftol>/<Fcut>]\n"
"                 Define the wavenumber integration interval\n"
"                 dk=(2*PI)/(<length>*rmax). rmax is the maximum \n"
"                 epicentral distance. \n"
"                 There are 4 cases:\n"
"                 + (default) not set or set 0.0.\n"); printf(
"                   <length> will be %.1f.\n", GRT_GREENFN_L_LENGTH); printf(
"                 + manually set one POSITIVE value, e.g. -Ll20\n"
"                 + manually set three POSITIVE values, \n"
"                   e.g. -Ll20/5/10, means split the integration \n"
"                   into two parts, [0, k*] and [k*, kmax], \n"
"                   in which k*=<Fcut>/rmax, and use DWM with\n"
"                   <length> and FIM with <Flength>, respectively.\n"
"                 + manually set three POSITIVE values, with -La,\n"
"                   in this case, <Flength> will be <Ftol> for Self-\n"
"                   Adaptive FIM.\n"
"\n"
"    -K[+k<k0>][+e<keps>][+v<vmin>]\n"
"                 Several parameters designed to define the\n"
"                 behavior in wavenumber integration. The upper\n"
"                 bound is k0,\n"
"                 <k0>:   default is %.1f, and \n", GRT_GREENFN_K_K0); printf(
"                         multiply PI/hs in program, \n"
"                         where hs = max(fabs(depsrc-deprcv), %.1f).\n", GRT_MIN_DEPTH_GAP_SRC_RCV); printf(
"                 <keps>: a threshold for break wavenumber \n"
"                         integration in advance. See \n"
"                         (Yao and Harkrider, 1983) for details.\n"
"                         Default 0.0 not use.\n"
"                 <vmin>: Minimum velocity (km/s) for reference. This\n"
"                         is designed to define the upper bound \n"
"                         of wavenumber integration.\n"
"                         There are 3 cases:\n"
"                         + (default) not set or set 0.0.\n"); printf(
"                           <vmin> will be the minimum velocity\n"
"                           of model, but limited to %.1f. and if \n", GRT_GREENFN_K_VMIN); printf(
"                           hs is thinner than %.1f km, PTAM will be appled\n", GRT_MIN_DEPTH_GAP_SRC_RCV); printf(
"                           automatically.\n"
"                         + manually set POSITIVE value. \n"
"                         + manually set NEGATIVE value, \n"
"                           and PTAM will be appled.\n"
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
"\n"
"    -h           Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt static greenfn -Mmilrow -D2/0 -X-10/10/20 -Y-10/10/20 -Ostgrn.nc\n"
"\n\n\n"
);
}





/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    char* command = Ctrl->name;

    // 先为个别参数设置非0初始值
    Ctrl->K.vmin = GRT_GREENFN_K_VMIN;
    Ctrl->K.k0 = GRT_GREENFN_K_K0;

    int opt;

    while ((opt = getopt(argc, argv, ":M:D:L:K:X:Y:O:Seh")) != -1) {
        switch (opt) {
            // 模型路径，其中每行分别为 
            //      厚度(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  Qp   Qs
            // 互相用空格隔开即可
            case 'M':
                Ctrl->M.active = true;
                Ctrl->M.s_modelpath = strdup(optarg);
                Ctrl->M.s_modelname = grt_get_basename(Ctrl->M.s_modelpath);
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

            // 波数积分间隔 -L[a|l]<length>[/<Flength>|<Ftol>/<Fcut>]
            case 'L':
                Ctrl->L.active = true;
                {   
                    // 若是纯数字，即未指定子选项，则直接使用DWM，并指定步长
                    if(isdigit(optarg[0])){
                        Ctrl->L.method = GRT_K_INTEG_METHOD_DWM;
                        // 仅接受一个值，若有多个值的分隔符则报错
                        if(strchr(optarg, '/') != NULL){
                            GRTBadOptionError(command, L, "single -L accept only 1 argument, but found %s.\n", optarg);
                        }
                        if(1 != sscanf(optarg, "%lf", &Ctrl->L.Length)){
                            GRTBadOptionError(command, L, "");
                        }
                    }
                    // 指定了子选项
                    else {
                        // 固定间隔 Filon 积分
                        if(optarg[0] == 'l'){
                            Ctrl->L.method = GRT_K_INTEG_METHOD_FIM;
                            if(3 != sscanf(optarg+1, "%lf/%lf/%lf", &Ctrl->L.Length, &Ctrl->L.filonLength, &Ctrl->L.filonCut)){
                                GRTBadOptionError(command, L, "");
                            }
                            if(Ctrl->L.filonLength <= 0.0){
                                GRTBadOptionError(command, L, "Flength should be positive.");
                            }
                        }
                        //  自适应采样
                        else if(optarg[0] == 'a'){
                            Ctrl->L.method = GRT_K_INTEG_METHOD_SAFIM;
                            if(3 != sscanf(optarg+1, "%lf/%lf/%lf", &Ctrl->L.Length, &Ctrl->L.safilonTol, &Ctrl->L.filonCut)){
                                GRTBadOptionError(command, L, "");
                            }
                            if(Ctrl->L.safilonTol <= 0.0){
                                GRTBadOptionError(command, L, "safilonTol should be positive.");
                            }
                        }

                        // 检查共有参数
                        if(Ctrl->L.Length <= 0.0){
                            GRTBadOptionError(command, L, "Length should be positive.");
                        }
                        if(Ctrl->L.filonCut < 0.0){
                            GRTBadOptionError(command, L, "Fcut should be nonnegative.");
                        }
                    }
                }
                break;

            // 波数积分相关变量 -K[+k<k0>][+e<keps>][+v<vmin>]
            case 'K':
                Ctrl->K.active = true;
                {
                char *line = strdup(optarg);
                char *token = strtok(line, "+");
                while(token != NULL){
                    switch(token[0]) {
                        case 'k':
                            if(1 != sscanf(token+1, "%lf", &Ctrl->K.k0)){
                                GRTBadOptionError(command, K+k, "");
                            }
                            if(Ctrl->K.k0 < 0.0){
                                GRTBadOptionError(command, K, "Can't set negative k0(%f).", Ctrl->K.k0);
                            }
                            break;

                        case 'e':
                            if(1 != sscanf(token+1, "%lf", &Ctrl->K.keps)){
                                GRTBadOptionError(command, K+e, "");
                            }
                            break;

                        case 'v':
                            Ctrl->K.v_active = true;
                            if(1 != sscanf(token+1, "%lf", &Ctrl->K.vmin)){
                                GRTBadOptionError(command, K+v, "");
                            }
                            break;

                        default:
                            GRTBadOptionError(command, K, "-K+%s is not supported.", token);
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
                        GRTBadOptionError(command, X, "");
                    };
                    if(delta <= 0){
                        GRTBadOptionError(command, X, "Can't set nonpositive dx(%f)", delta);
                    }
                    if(a1 > a2){
                        GRTBadOptionError(command, X, "x1(%f) > x2(%f).", a1, a2);
                    }

                    Ctrl->X.nx = floor((a2-a1)/delta) + 1;
                    Ctrl->X.xs = (real_t*)calloc(Ctrl->X.nx, sizeof(real_t));
                    for(size_t i=0; i<Ctrl->X.nx; ++i){
                        Ctrl->X.xs[i] = a1 + delta*i;
                    }
                }
                break;

            // Y坐标数组，-Yy1/y2/dy
            case 'Y':
                Ctrl->Y.active = true;
                {
                    real_t a1, a2, delta;
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
                    Ctrl->Y.ys = (real_t*)calloc(Ctrl->Y.ny, sizeof(real_t));
                    for(size_t i=0; i<Ctrl->Y.ny; ++i){
                        Ctrl->Y.ys[i] = a1 + delta*i;
                    }
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
            
            GRT_Common_Options_in_Switch(command, (char)(optopt));
        }
    }

    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(command, argc > 1);
    GRTCheckOptionActive(command, Ctrl, M);
    GRTCheckOptionActive(command, Ctrl, D);
    GRTCheckOptionActive(command, Ctrl, X);
    GRTCheckOptionActive(command, Ctrl, Y);
    GRTCheckOptionActive(command, Ctrl, O);

    // 设置震中距数组
    Ctrl->nr = Ctrl->X.nx*Ctrl->Y.ny;
    Ctrl->rs = (real_t*)calloc(Ctrl->nr, sizeof(real_t));
    for(size_t ix=0; ix<Ctrl->X.nx; ++ix){
        for(size_t iy=0; iy<Ctrl->Y.ny; ++iy){
            Ctrl->rs[iy + ix*Ctrl->Y.ny] = GRT_MAX(sqrt(GRT_SQUARE(Ctrl->X.xs[ix]) + GRT_SQUARE(Ctrl->Y.ys[iy])), GRT_MIN_DISTANCE);  // 避免0震中距
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
    if((Ctrl->M.mod1d = grt_read_mod1d_from_file(command, Ctrl->M.s_modelpath, Ctrl->D.depsrc, Ctrl->D.deprcv, false)) == NULL){
        exit(EXIT_FAILURE);
    }
    GRT_MODEL1D *mod1d = Ctrl->M.mod1d;

    // 最大最小速度
    real_t vmin, vmax;
    grt_get_mod1d_vmin_vmax(mod1d, &vmin, &vmax);

    // 参考最小速度
    if(!Ctrl->K.v_active){
        Ctrl->K.vmin = GRT_MAX(vmin, GRT_GREENFN_K_VMIN);
    } 

    // 如果没有主动设置vmin_ref，则判断是否要自动使用PTAM
    if( !Ctrl->K.v_active && fabs(Ctrl->D.deprcv - Ctrl->D.depsrc) <= GRT_MIN_DEPTH_GAP_SRC_RCV) {
        Ctrl->K.vmin = - fabs(Ctrl->K.vmin);
    }
    
    // 设置积分间隔默认值
    if(Ctrl->L.Length == 0.0)  Ctrl->L.Length = GRT_GREENFN_L_LENGTH;

    // 波数积分输出目录
    if(Ctrl->S.active){
        Ctrl->S.s_statsdir = NULL;
        GRT_SAFE_ASPRINTF(&Ctrl->S.s_statsdir, "stgrtstats");
        // 建立保存目录
        GRTCheckMakeDir(command, Ctrl->S.s_statsdir);
        GRT_SAFE_ASPRINTF(&Ctrl->S.s_statsdir, "%s/%s_%s_%s", Ctrl->S.s_statsdir, Ctrl->M.s_modelname, Ctrl->D.s_depsrc, Ctrl->D.s_deprcv);
        GRTCheckMakeDir(command, Ctrl->S.s_statsdir);
    }

    // 建立格林函数的浮点数
    realChnlGrid *grn = (realChnlGrid *) calloc(Ctrl->nr, sizeof(*grn));
    realChnlGrid *grn_uiz = (Ctrl->e.active)? (realChnlGrid *) calloc(Ctrl->nr, sizeof(*grn_uiz)) : NULL;
    realChnlGrid *grn_uir = (Ctrl->e.active)? (realChnlGrid *) calloc(Ctrl->nr, sizeof(*grn_uir)) : NULL;


    //==============================================================================
    // 计算静态格林函数
    grt_integ_static_grn(
        mod1d, Ctrl->nr, Ctrl->rs, Ctrl->K.vmin, Ctrl->K.keps, Ctrl->K.k0, Ctrl->L.Length, Ctrl->L.filonLength, Ctrl->L.safilonTol, Ctrl->L.filonCut, 
        grn, Ctrl->e.active, grn_uiz, grn_uir,
        Ctrl->S.s_statsdir
    );
    //==============================================================================

    real_t src_va = mod1d->Va[mod1d->isrc];
    real_t src_vb = mod1d->Vb[mod1d->isrc];
    real_t src_rho = mod1d->Rho[mod1d->isrc];
    real_t rcv_va = mod1d->Va[mod1d->ircv];
    real_t rcv_vb = mod1d->Vb[mod1d->ircv];
    real_t rcv_rho = mod1d->Rho[mod1d->ircv];


    // ==================================================================================
    // 将结果保存为 nc 格式
    // ==================================================================================
    int ncid, x_dimid, y_dimid;
    const int ndims = 2;
    int dimids[ndims];
    int x_varid, y_varid;
    intChnlGrid u_varids;
    intChnlGrid uiz_varids;
    intChnlGrid uir_varids;

    // 创建 NC 文件
    NC_CHECK(nc_create(Ctrl->O.s_outgrid, NC_CLOBBER, &ncid));

    // 写入全局属性
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, NC_GLOBAL, "src_va", NC_REAL, 1, &src_va));
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, NC_GLOBAL, "src_vb", NC_REAL, 1, &src_vb));
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, NC_GLOBAL, "src_rho", NC_REAL, 1, &src_rho));
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, NC_GLOBAL, "rcv_va", NC_REAL, 1, &rcv_va));
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, NC_GLOBAL, "rcv_vb", NC_REAL, 1, &rcv_vb));
    NC_CHECK(NC_FUNC_REAL(nc_put_att) (ncid, NC_GLOBAL, "rcv_rho", NC_REAL, 1, &rcv_rho));
    // 是否计算了位移偏导也直接写到全局属性
    {
        int tmp = Ctrl->e.active;
        NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "calc_upar", NC_INT, 1, &tmp));
    }

    // 定义维度
    NC_CHECK(nc_def_dim(ncid, "north", Ctrl->X.nx, &x_dimid));
    NC_CHECK(nc_def_dim(ncid, "east", Ctrl->Y.ny, &y_dimid));
    dimids[0] = x_dimid;
    dimids[1] = y_dimid;

    // 定义维度数组
    NC_CHECK(nc_def_var(ncid, "north", NC_REAL, 1, &x_dimid, &x_varid));
    NC_CHECK(nc_def_var(ncid, "east", NC_REAL, 1, &y_dimid, &y_varid));

    // 定义不同震源不同分量的格林函数数组
    GRT_LOOP_ChnlGrid(im, c){
        int modr = GRT_SRC_M_ORDERS[im];
        char *s_title = NULL;

        if(modr==0 && GRT_ZRT_CODES[c]=='T')  continue;

        GRT_SAFE_ASPRINTF(&s_title, "%s%c", GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
        NC_CHECK(nc_def_var(ncid, s_title, NC_REAL, ndims, dimids, &u_varids[im][c]));

        // 位移偏导
        if(Ctrl->e.active){
            GRT_SAFE_ASPRINTF(&s_title, "z%s%c", GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
            NC_CHECK(nc_def_var(ncid, s_title, NC_REAL, ndims, dimids, &uiz_varids[im][c]));
            GRT_SAFE_ASPRINTF(&s_title, "r%s%c", GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
            NC_CHECK(nc_def_var(ncid, s_title, NC_REAL, ndims, dimids, &uir_varids[im][c]));
        }
        GRT_SAFE_FREE_PTR(s_title);
    }

    // 结束定义模式
    NC_CHECK(nc_enddef(ncid));

    // 写入数据
    NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, x_varid, Ctrl->X.xs));
    NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, y_varid, Ctrl->Y.ys));
    real_t *tmpdata = (real_t *)calloc(Ctrl->nr, sizeof(real_t));
    GRT_LOOP_ChnlGrid(im, c){
        int modr = GRT_SRC_M_ORDERS[im];

        if(modr==0 && GRT_ZRT_CODES[c]=='T')  continue;

        int sgn0 = 1;
        sgn0 = (GRT_ZRT_CODES[c]=='Z')? -1 : 1;
        for(size_t ir=0; ir < Ctrl->nr; ++ir){
            tmpdata[ir] = sgn0 * grn[ir][im][c];
        }

        NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, u_varids[im][c], tmpdata));

        // 位移偏导
        if(Ctrl->e.active){
            for(size_t ir=0; ir < Ctrl->nr; ++ir){
                tmpdata[ir] = (-1) * sgn0 * grn_uiz[ir][im][c];  // 这里多乘的(-1)是因为对z的偏导，z需反向
            }
            NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, uiz_varids[im][c], tmpdata));
            for(size_t ir=0; ir < Ctrl->nr; ++ir){
                tmpdata[ir] = sgn0 * grn_uir[ir][im][c];  // 这里多乘的(-1)是因为对z的偏导，z需反向
            }
            NC_CHECK(NC_FUNC_REAL(nc_put_var) (ncid, uir_varids[im][c], tmpdata));
        }
    }
    GRT_SAFE_FREE_PTR(tmpdata);

    // 关闭文件
    NC_CHECK(nc_close(ncid));


    // 释放内存
    GRT_SAFE_FREE_PTR(grn);
    GRT_SAFE_FREE_PTR(grn_uiz);
    GRT_SAFE_FREE_PTR(grn_uir);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}

