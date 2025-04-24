/**
 * @file   stgrt_main.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 *    计算静态位移
 * 
 */

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>

#include "static/stgrt.h"
#include "common/const.h"
#include "common/model.h"
#include "common/colorstr.h"
#include "common/logo.h"
#include "common/integral.h"
#include "common/iostats.h"
#include "common/search.h"

static char *command;
static PYMODEL1D *pymod;
static double depsrc, deprcv;

//****************** 在该文件以内的全局变量 ***********************//
// 命令名称
static char *command = NULL;
// 模型路径，模型PYMODEL1D指针，全局最大最小速度
static char *s_modelpath = NULL;
static char *s_modelname = NULL;
static PYMODEL1D *pymod;
static double vmax, vmin;
// 震源和场点深度
static double depsrc, deprcv;
static char *s_depsrc = NULL, *s_deprcv = NULL;
// 波数积分间隔, Filon积分间隔，Filon积分起始点
static double Length=0.0, filonLength=0.0, filonCut=0.0;
static double Length0=15.0; // 默认Length
// 波数积分相关变量
static double keps=-1.0, k0=5.0;
// 参考最小速度，小于0表示使用峰谷平均法;
static double vmin_ref=0.0;
static const double min_vmin_ref=0.1;
// 自动使用峰谷平均法的最小厚度差
static const double hs_ptam = MIN_DEPTH_GAP_SRC_RCV;
// 接收点位置数组
static MYREAL *rs = NULL;
static MYREAL *xs = NULL;
static MYREAL *ys = NULL;
static int nr=0, nx=0, ny=0;

// 输出波数积分过程文件
static char *s_statsdir = NULL;

// 是否计算位移空间导数
static bool calc_upar=false;

// 各选项的标志变量，初始化为0，定义了则为1
static int M_flag=0, D_flag=0, 
            L_flag=0, V_flag=0,
            K_flag=0, S_flag=0,
            X_flag=0, Y_flag=0, 
            e_flag=0;

// 三分量代号
const char chs[3] = {'Z', 'R', 'T'};

/**
 * 打印使用说明
 */
static void print_help(){
print_logo();
printf("\n"
"[stgrt]\n\n"
"    Compute static Green's Functions, output to stdout. \n"
"    The units and components are consistent with the dynamics, \n"
"    check \"grt -h\" for details.\n"
"\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    stgrt -M<model> -D<depsrc>/<deprcv> -X<x1>/<x2>/<nx> \n"
"          -Y<y1>/<y2>/<ny>  [-L<length>] [-V<vmin_ref>] \n" 
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
"    -X<x1>/<x2>/<nx>\n"
"                 Set the equidistant points in the north direction.\n"
"                 <x1>: start coordinate (km).\n"
"                 <x2>: end coordinate (km).\n"
"                 <nx>: number of points.\n"
"\n"
"    -Y<y1>/<y2>/<ny>\n"
"                 Set the equidistant points in the east direction.\n"
"                 <y1>: start coordinate (km).\n"
"                 <y2>: end coordinate (km).\n"
"                 <ny>: number of points.\n"
"\n"
"    -L<length>[/<Flength>/<Fcut>]\n"
"                 Define the wavenumber integration interval\n"
"                 dk=(2*PI)/(<length>*rmax). rmax is the maximum \n"
"                 epicentral distance. \n"
"                 There are 3 cases:\n"
"                 + (default) not set or set %.1f.\n", Length); printf(
"                   <length> will be %.1f.\n", Length0); printf(
"                 + manually set one POSITIVE value, e.g. -L20\n"
"                 + manually set three POSITIVE values, \n"
"                   e.g. -L20/5/10, means split the integration \n"
"                   into two parts, [0, k*] and [k*, kmax], \n"
"                   in which k*=<Fcut>/rmax, and use DWM with\n"
"                   <length> and FIM with <Flength>, respectively.\n"
"\n"
"    -V<vmin_ref> \n"
"                 (Inherited from the dynamic case, and the numerical\n"
"                 value will not be used in here, except its sign.)\n"
"                 + (default) not set or set %.1f.\n", vmin_ref); printf(
"                   <vmin_ref> will be the minimum velocity\n"
"                   of model, but limited to %.1f. and if the \n", min_vmin_ref); printf(
"                   depth gap between source and receiver is \n"
"                   thinner than %.1f km, PTAM will be appled\n", hs_ptam); printf(
"                   automatically.\n"
"                 + manually set POSITIVE value. \n"
"                 + manually set NEGATIVE value, \n"
"                   and PTAM will be appled.\n"
"\n"
"    -K<k0>[/<keps>]\n"
"                 Several parameters designed to define the\n"
"                 behavior in wavenumber integration. The upper\n"
"                 bound is k0,\n"
"                 <k0>:   default is %.1f, and \n", k0); printf(
"                         multiply PI/hs in program, \n"
"                         where hs = max(fabs(depsrc-deprcv), %.1f).\n", MIN_DEPTH_GAP_SRC_RCV); printf(
"                 <keps>: a threshold for break wavenumber \n"
"                         integration in advance. See \n"
"                         (Yao and Harkrider, 1983) for details.\n"
"                         Default %.1f not use.\n", keps); printf(
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
"    stgrt -Mmilrow -D2/0 -X-10/10/20 -Y-10/10/20 > grn\n"
"\n\n\n"
);
}


/**
 * 从路径字符串中找到用/或\\分隔的最后一项
 * 
 * @param    path     路径字符串指针
 * 
 * @return   指向最后一项字符串的指针
 */
static char* get_basename(char* path) {
    // 找到最后一个 '/'
    char* last_slash = strrchr(path, '/'); 
    
#ifdef _WIN32
    char* last_backslash = strrchr(path, '\\');
    if (last_backslash && (!last_slash || last_backslash > last_slash)) {
        last_slash = last_backslash;
    }
#endif
    if (last_slash) {
        // 返回最后一个 '/' 之后的部分
        return last_slash + 1; 
    }
    // 如果没有 '/'，整个路径就是最后一项
    return path; 
}


/**
 * 从命令行中读取选项，处理后记录到全局变量中
 * 
 * @param     argc      命令行的参数个数
 * @param     argv      多个参数字符串指针
 */
static void getopt_from_command(int argc, char **argv){
    int opt;
    while ((opt = getopt(argc, argv, ":M:D:L:K:X:Y:V:Seh")) != -1) {
        switch (opt) {
            // 模型路径，其中每行分别为 
            //      厚度(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  Qp   Qs
            // 互相用空格隔开即可
            case 'M':
                M_flag = 1;
                s_modelpath = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                strcpy(s_modelpath, optarg);
                if(access(s_modelpath, F_OK) == -1){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! File \"%s\" set by -M not exists.\n" DEFAULT_RESTORE, command, s_modelpath);
                    exit(EXIT_FAILURE);
                }
            
                s_modelname = get_basename(s_modelpath);
                break;

            // 震源和场点深度， -Ddepsrc/deprcv
            case 'D':
                D_flag = 1;
                s_depsrc = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                s_deprcv = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                if(2 != sscanf(optarg, "%[^/]/%s", s_depsrc, s_deprcv)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -D.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                if(1 != sscanf(s_depsrc, "%lf", &depsrc)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -D.\n"  DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                if(1 != sscanf(s_deprcv, "%lf", &deprcv)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -D.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                if(depsrc < 0.0 || deprcv < 0.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Negative value in -D is not supported.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                break;

            // 波数积分间隔 -L<length>[/<Flength>/<Fcut>]
            case 'L':
                L_flag = 1;
                {
                    int n = sscanf(optarg, "%lf/%lf/%lf", &Length, &filonLength, &filonCut);
                    if(n != 1 && n != 3){
                        fprintf(stderr, "[%s] " BOLD_RED "Error in -L.\n" DEFAULT_RESTORE, command);
                        exit(EXIT_FAILURE);
                    };
                    if(n == 1 && Length <= 0){
                        fprintf(stderr, "[%s] " BOLD_RED "Error! In -L, length should be positive.\n" DEFAULT_RESTORE, command);
                        exit(EXIT_FAILURE);
                    }
                    if(n == 3 && (filonLength <= 0 || filonCut < 0)){
                        fprintf(stderr, "[%s] " BOLD_RED "Error! In -L, Flength should be positive, Fcut should be nonnegative.\n" DEFAULT_RESTORE, command);
                        exit(EXIT_FAILURE);
                    }
                }
                
                break;

            // 参考最小速度 -Vvmin_ref
            case 'V':
                V_flag = 1;
                if(0 == sscanf(optarg, "%lf", &vmin_ref)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -V.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                break;

            // 波数积分相关变量 -Kk0/keps
            case 'K':
                K_flag = 1;
                if(0 == sscanf(optarg, "%lf/%lf", &k0, &keps)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -K.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                
                if(k0 < 0.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Can't set negative k0(%f) in -K.\n" DEFAULT_RESTORE, command, k0);
                    exit(EXIT_FAILURE);
                }
                break;

            // X坐标数组，-Xx1/x2/nx
            case 'X':
                X_flag = 1;
                {
                    MYREAL a1, a2;
                    if(3 != sscanf(optarg, "%lf/%lf/%d", &a1, &a2, &nx)){
                        fprintf(stderr, "[%s] " BOLD_RED "Error in -X.\n" DEFAULT_RESTORE, command);
                        exit(EXIT_FAILURE);
                    };
                    if(nx <= 0){
                        fprintf(stderr, "[%s] " BOLD_RED "Error! Can't set nonpositive nx(%d) in -X.\n" DEFAULT_RESTORE, command, nx);
                        exit(EXIT_FAILURE);
                    }
                    if(a1 > a2){
                        fprintf(stderr, "[%s] " BOLD_RED "Error! x1(%f) > x2(%f) in -X.\n" DEFAULT_RESTORE, command, a1, a2);
                        exit(EXIT_FAILURE);
                    }

                    xs = (MYREAL*)calloc(nx, sizeof(MYREAL));
                    MYREAL delta = (a2 - a1)/((nx>1)? nx-1 : 1);
                    for(int i=0; i<nx; ++i){
                        xs[i] = a1 + delta*i;
                    }
                }
                break;

            // Y坐标数组，-Yy1/y2/ny
            case 'Y':
                Y_flag = 1;
                {
                    MYREAL a1, a2;
                    if(3 != sscanf(optarg, "%lf/%lf/%d", &a1, &a2, &ny)){
                        fprintf(stderr, "[%s] " BOLD_RED "Error in -Y.\n" DEFAULT_RESTORE, command);
                        exit(EXIT_FAILURE);
                    };
                    if(ny <= 0){
                        fprintf(stderr, "[%s] " BOLD_RED "Error! Can't set nonpositive ny(%d) in -Y.\n" DEFAULT_RESTORE, command, ny);
                        exit(EXIT_FAILURE);
                    }
                    if(a1 > a2){
                        fprintf(stderr, "[%s] " BOLD_RED "Error! y1(%f) > y2(%f) in -Y.\n" DEFAULT_RESTORE, command, a1, a2);
                        exit(EXIT_FAILURE);
                    }

                    ys = (MYREAL*)calloc(ny, sizeof(MYREAL));
                    MYREAL delta = (a2 - a1)/((ny>1)? ny-1 : 1);
                    for(int i=0; i<ny; ++i){
                        ys[i] = a1 + delta*i;
                    }
                }
                break;

            // 输出波数积分中间文件
            case 'S':
                S_flag = 1;
                break;

            // 是否计算位移空间导数
            case 'e':
                e_flag = 1;
                calc_upar = true;
                break;
            
            // 帮助
            case 'h':
                print_help();
                exit(EXIT_SUCCESS);
                break;

            // 参数缺失
            case ':':
                fprintf(stderr, "[%s] " BOLD_RED "Error! Option '-%c' requires an argument. Use '-h' for help.\n" DEFAULT_RESTORE, command, optopt);
                exit(EXIT_FAILURE);
                break;

            // 非法选项
            case '?':
            default:
                fprintf(stderr, "[%s] " BOLD_RED "Error! Option '-%c' is invalid. Use '-h' for help.\n" DEFAULT_RESTORE, command, optopt);
                exit(EXIT_FAILURE);
                break;
        }
    }

    // 检查必须设置的参数是否有设置
    if(argc == 1){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set options. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(M_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -M. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(D_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -D. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(X_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -X. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(Y_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -Y. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }

    // 设置震中距数组
    nr = nx*ny;
    rs = (MYREAL*)calloc(nr, sizeof(MYREAL));
    for(int iy=0; iy<ny; ++iy){
        for(int ix=0; ix<nx; ++ix){
            rs[ix + iy*nx] = SQRT(xs[ix]*xs[ix] + ys[iy]*ys[iy]);
            if(rs[ix + iy*nx] < 1e-5)  rs[ix + iy*nx] = 1e-5;  // 避免0震中距
        }
    }


}





int main(int argc, char **argv){
    command = argv[0];

    // 传入参数 
    getopt_from_command(argc, argv);

    // 读入模型文件
    if((pymod = read_pymod_from_file(command, s_modelpath, depsrc, deprcv)) ==NULL){
        exit(EXIT_FAILURE);
    }

    // 最大最小速度
    vmax = pymod->Va[findMinMax_MYREAL(pymod->Va, pymod->n, true)];
    vmin = pymod->Vb[findMinMax_MYREAL(pymod->Vb, pymod->n, false)];
    if(vmin > vmax) {
        double tmp;
        tmp = vmin; vmin = vmax; vmax = tmp;
    }

    // 参考最小速度
    if(vmin_ref == 0.0){
        vmin_ref = vmin;
        if(vmin_ref < min_vmin_ref) vmin_ref = min_vmin_ref;
    } 

    // 如果没有主动设置vmin_ref，则判断是否要自动使用PTAM
    if(V_flag == 0 && fabs(deprcv - depsrc) <= hs_ptam) {
        vmin_ref = - fabs(vmin_ref);
    }
    
    // 设置积分间隔默认值
    if(Length == 0.0)  Length = Length0;

    // 波数积分输出目录
    if(S_flag==1){
        s_statsdir = (char*)malloc(sizeof(char)*(strlen(s_modelpath)+strlen(s_depsrc)+strlen(s_deprcv)+100));
        sprintf(s_statsdir, "stgrtstats");
        // 建立保存目录
        if(mkdir(s_statsdir, 0777) != 0){
            if(errno != EEXIST){
                fprintf(stderr, "[%s] " BOLD_RED "Error! Unable to create folder %s. Error code: %d\n" DEFAULT_RESTORE, command, s_statsdir, errno);
                exit(EXIT_FAILURE);
            }
        }
        sprintf(s_statsdir, "%s/%s_%s_%s", s_statsdir, s_modelname, s_depsrc, s_deprcv);
        if(mkdir(s_statsdir, 0777) != 0){
            if(errno != EEXIST){
                fprintf(stderr, "[%s] " BOLD_RED "Error! Unable to create folder %s. Error code: %d\n" DEFAULT_RESTORE, command, s_statsdir, errno);
                exit(EXIT_FAILURE);
            }
        }
    }


    // 建立格林函数的complex数组
    MYREAL (*EXPgrn)[2] = (MYREAL(*)[2])calloc(nr, sizeof(*EXPgrn));
    MYREAL (*VFgrn)[2]  = (MYREAL(*)[2])calloc(nr, sizeof(*VFgrn));
    MYREAL (*HFgrn)[3]  = (MYREAL(*)[3])calloc(nr, sizeof(*HFgrn));
    MYREAL (*DDgrn)[2]  = (MYREAL(*)[2])calloc(nr, sizeof(*DDgrn));
    MYREAL (*DSgrn)[3]  = (MYREAL(*)[3])calloc(nr, sizeof(*DSgrn));
    MYREAL (*SSgrn)[3]  = (MYREAL(*)[3])calloc(nr, sizeof(*SSgrn));

    MYREAL (*EXPgrn_uiz)[2] = (calc_upar)? (MYREAL(*)[2])calloc(nr, sizeof(*EXPgrn)) : NULL;
    MYREAL (*VFgrn_uiz)[2]  = (calc_upar)? (MYREAL(*)[2])calloc(nr, sizeof(*VFgrn)) : NULL;
    MYREAL (*HFgrn_uiz)[3]  = (calc_upar)? (MYREAL(*)[3])calloc(nr, sizeof(*HFgrn)) : NULL;
    MYREAL (*DDgrn_uiz)[2]  = (calc_upar)? (MYREAL(*)[2])calloc(nr, sizeof(*DDgrn)) : NULL;
    MYREAL (*DSgrn_uiz)[3]  = (calc_upar)? (MYREAL(*)[3])calloc(nr, sizeof(*DSgrn)) : NULL;
    MYREAL (*SSgrn_uiz)[3]  = (calc_upar)? (MYREAL(*)[3])calloc(nr, sizeof(*SSgrn)) : NULL;

    MYREAL (*EXPgrn_uir)[2] = (calc_upar)? (MYREAL(*)[2])calloc(nr, sizeof(*EXPgrn)) : NULL;
    MYREAL (*VFgrn_uir)[2]  = (calc_upar)? (MYREAL(*)[2])calloc(nr, sizeof(*VFgrn)) : NULL;
    MYREAL (*HFgrn_uir)[3]  = (calc_upar)? (MYREAL(*)[3])calloc(nr, sizeof(*HFgrn)) : NULL;
    MYREAL (*DDgrn_uir)[2]  = (calc_upar)? (MYREAL(*)[2])calloc(nr, sizeof(*DDgrn)) : NULL;
    MYREAL (*DSgrn_uir)[3]  = (calc_upar)? (MYREAL(*)[3])calloc(nr, sizeof(*DSgrn)) : NULL;
    MYREAL (*SSgrn_uir)[3]  = (calc_upar)? (MYREAL(*)[3])calloc(nr, sizeof(*SSgrn)) : NULL;

    //==============================================================================
    // 计算静态格林函数
    integ_static_grn(
        pymod, nr, rs, vmin_ref, keps, k0, Length, filonLength, filonCut, 
        EXPgrn, VFgrn, HFgrn, DDgrn, DSgrn, SSgrn,
        calc_upar, 
        EXPgrn_uiz, VFgrn_uiz, HFgrn_uiz, DDgrn_uiz, DSgrn_uiz, SSgrn_uiz, 
        EXPgrn_uir, VFgrn_uir, HFgrn_uir, DDgrn_uir, DSgrn_uir, SSgrn_uir, 
        s_statsdir
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

    // 定义标题数组
    const char *titles[15] = {
        "EXZ", "EXR", "VFZ", "VFR", "HFZ", "HFR", "HFT",
        "DDZ", "DDR", "DSZ", "DSR", "DST", "SSZ", "SSR", "SST"
    };
    const char *upar_titles[30] = {
        "zEXZ", "zEXR", "zVFZ", "zVFR", "zHFZ", "zHFR", "zHFT",
        "zDDZ", "zDDR", "zDSZ", "zDSR", "zDST", "zSSZ", "zSSR", "zSST",
        "rEXZ", "rEXR", "rVFZ", "rVFR", "rHFZ", "rHFR", "rHFT",
        "rDDZ", "rDDR", "rDSZ", "rDSR", "rDST", "rSSZ", "rSSR", "rSST"
    };

    // 输出标题
    char XX[20];
    sprintf(XX, GRT_STRING_FMT, "X(km)"); XX[0]='#';
    fprintf(stdout, "%s", XX);
    fprintf(stdout, GRT_STRING_FMT, "Y(km)");
    for(int i=0; i<15; ++i) fprintf(stdout, GRT_STRING_FMT, titles[i]);
    if(calc_upar) {
        for(int i=0; i<30; ++i)  fprintf(stdout, GRT_STRING_FMT, upar_titles[i]);
    }
    fprintf(stdout, "\n");

    // 写结果
    for(int iy=0; iy<ny; ++iy) {
        for(int ix=0; ix<nx; ++ix) {
            int ir = ix + iy * nx;
            fprintf(stdout, GRT_REAL_FMT GRT_REAL_FMT, xs[ix], ys[iy]);
            MYREAL *grns[] = {
                EXPgrn[ir], VFgrn[ir], HFgrn[ir], DDgrn[ir], DSgrn[ir], SSgrn[ir]
            };
            int grn_sizes[] = {2, 2, 3, 2, 3, 3};
            // 对Z分量反向
            for(int i=0; i<6; ++i) {
                for (int j=0; j<grn_sizes[i]; ++j)
                    fprintf(stdout, GRT_REAL_FMT, (j == 0 ? -1.0 : 1.0) * grns[i][j]);
            }
            if(calc_upar) {
                // 前6个是对z的偏导，注意后续对z符号的判断
                MYREAL *upar_grns[] = {
                    EXPgrn_uiz[ir], VFgrn_uiz[ir], HFgrn_uiz[ir], DDgrn_uiz[ir], DSgrn_uiz[ir], SSgrn_uiz[ir],
                    EXPgrn_uir[ir], VFgrn_uir[ir], HFgrn_uir[ir], DDgrn_uir[ir], DSgrn_uir[ir], SSgrn_uir[ir]
                };
                for(int i=0; i<12; ++i) {
                    for(int j=0; j<grn_sizes[i % 6]; ++j)
                        fprintf(stdout, GRT_REAL_FMT, (j == 0 ? -1.0 : 1.0) * (i < 6 ? -1.0 : 1.0) * upar_grns[i][j]);
                }
            }
            fprintf(stdout, "\n");
        }
    }

}

