/**
 * @file   stgrt_syn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 *    根据计算好的静态格林函数，定义震源机制以及方位角等，生成合成的静态三分量位移场
 * 
 */



#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>

#include "common/const.h"
#include "common/logo.h"
#include "common/colorstr.h"
#include "common/radiation.h"
#include "common/coord.h"

extern char *optarg;
extern int optind;
extern int optopt;

//****************** 在该文件以内的全局变量 ***********************//
// 命令名称
static char *command = NULL;
// 放大系数，对于剪切源、爆炸源、张量震源，M0是标量地震矩；对于单力源，M0是放大系数
static double M0 = 0.0;
// 在放大系数上是否需要乘上震源处的剪切模量
static bool mult_src_mu = false;
// 存储不同震源的震源机制相关参数的数组
static double mchn[MECHANISM_NUM] = {0};
// 最终要计算的震源类型
static int computeType=GRT_SYN_COMPUTE_EX;
static char s_computeType[3] = "EX";

// 是否计算位移空间导数
static bool calc_upar=false;

// 输出分量格式，即是否需要旋转到ZNE
static bool rot2ZNE = false;

// 各选项的标志变量，初始化为0，定义了则为1
static int S_flag=0, M_flag=0, F_flag=0,
           T_flag=0, N_flag=0,
           e_flag=0;

// 计算和位移相关量的种类（1-位移，2-ui_z，3-ui_r，4-ui_t）
static int calcUTypes=1;

/**
 * 打印使用说明
 */
static void print_help(){
print_logo();
printf("\n"
"[stgrt.syn]\n\n"
"    Compute static displacement with the outputs of \n"
"    command `stgrt` (reading from stdin).\n"
"    Three components are:\n"
"       + Up (Z),\n"
"       + Radial Outward (R),\n"
"       + Transverse Clockwise (T),\n"
"    and the units are cm. You can add -N to rotate ZRT to ZNE.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    stgrt.syn -S[u]<scale> \n"
"              [-M<strike>/<dip>/<rake>]\n"
"              [-T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>]\n"
"              [-F<fn>/<fe>/<fz>] \n"
"              [-N] [-e]\n"
"              < <grn>\n"
"\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -S[u]<scale>  Scale factor to all kinds of source. \n"
"                  + For Explosion, Shear and Moment Tensor,\n"
"                    unit of <scale> is dyne-cm. \n"
"                  + For Single Force, unit of <scale> is dyne.\n"
"                  + Since \"\\mu\" exists in scalar seismic moment\n"
"                    (\\mu*A*D), you can simply set -Su<scale>, <scale>\n"
"                    equals A*D (Area*Slip, [cm^3]), and <scale> will \n"
"                    multiply \\mu automatically in program.\n"
"\n"
"    For source type, you can only set at most one of\n"
"    '-M', '-T' and '-F'. If none, an Explosion is used.\n"
"\n"
"    -M<strike>/<dip>/<rake>\n"
"                  Three angles to define a fault. \n"
"                  The angles are in degree.\n"
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
"    -N            Components of results will be Z, N, E.\n"
"\n"
"    -e            Compute the spatial derivatives, ui_z and ui_r,\n"
"                  of displacement u. In filenames, prefix \"r\" means \n"
"                  ui_r and \"z\" means ui_z. \n"
"\n"
"    -h            Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    Say you have computed Static Green's functions with following command:\n"
"        stgrt -Mmilrow -D2/0 -X-5/5/10 -Y-5/5/10 > grn\n"
"\n"
"    Then you can get static displacement of Explosion\n"
"        stgrt.syn -Su1e16 < grn > syn_ex\n"
"\n"
"    or Shear\n"
"        stgrt.syn -Su1e16 -M100/20/80 < grn > syn_dc\n"
"\n"
"    or Single Force\n"
"        stgrt.syn -S1e20 -F0.5/-1.2/3.3 < grn > syn_sf\n"
"\n"
"    or Moment Tensor\n"
"        stgrt.syn -Su1e16 -T2.3/0.2/-4.0/0.3/0.5/1.2 < grn > syn_mt\n"
"\n\n\n"
"\n"
);
}


/**
 * 从命令行中读取选项，处理后记录到全局变量中
 * 
 * @param     argc      命令行的参数个数
 * @param     argv      多个参数字符串指针
 */
static void getopt_from_command(int argc, char **argv){
    int opt;
    while ((opt = getopt(argc, argv, ":S:M:F:T:Neh")) != -1) {
        switch (opt) {
            // 放大系数
            case 'S':
                S_flag = 1;
                {   
                    // 检查是否存在字符u，若存在表明需要乘上震源处的剪切模量
                    char *upos=NULL;
                    if((upos=strchr(optarg, 'u')) != NULL){
                        mult_src_mu = true;
                        *upos = ' ';
                    }
                }
                
                if(0 == sscanf(optarg, "%lf", &M0)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -S.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                break;

            // 剪切震源
            case 'M':
                M_flag = 1; 
                computeType = GRT_SYN_COMPUTE_DC;
                double strike, dip, rake;
                sprintf(s_computeType, "%s", "DC");
                if(3 != sscanf(optarg, "%lf/%lf/%lf", &strike, &dip, &rake)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -M.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                if(strike < 0.0 || strike > 360.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Strike in -M must be in [0, 360].\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                if(dip < 0.0 || dip > 90.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Dip in -M must be in [0, 90].\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                if(rake < -180.0 || rake > 180.0){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! Rake in -M must be in [-180, 180].\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }
                mchn[0] = strike;
                mchn[1] = dip;
                mchn[2] = rake;
                break;

            // 单力源
            case 'F':
                F_flag = 1;
                computeType = GRT_SYN_COMPUTE_SF;
                double fn, fe, fz;
                sprintf(s_computeType, "%s", "SF");
                if(3 != sscanf(optarg, "%lf/%lf/%lf", &fn, &fe, &fz)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -F.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                mchn[0] = fn;
                mchn[1] = fe;
                mchn[2] = fz;
                break;

            // 张量震源
            case 'T':
                T_flag = 1;
                computeType = GRT_SYN_COMPUTE_MT;
                double Mxx, Mxy, Mxz, Myy, Myz, Mzz;
                sprintf(s_computeType, "%s", "MT");
                if(6 != sscanf(optarg, "%lf/%lf/%lf/%lf/%lf/%lf", &Mxx, &Mxy, &Mxz, &Myy, &Myz, &Mzz)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -T.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                };
                mchn[0] = Mxx;
                mchn[1] = Mxy;
                mchn[2] = Mxz;
                mchn[3] = Myy;
                mchn[4] = Myz;
                mchn[5] = Mzz;
                break;

            // 是否计算位移空间导数
            case 'e':
                e_flag = 1;
                calc_upar = true;
                calcUTypes = 4;
                break;

            // 是否旋转到ZNE
            case 'N':
                N_flag = 1;
                rot2ZNE = true;
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

    // 检查必选项有没有设置
    if(argc == 1){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set options. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
    if(S_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -S. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }

    // 只能使用一种震源
    if(M_flag + F_flag + T_flag > 1){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Only support at most one of '-M', '-F' and '-T'. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
}




//====================================================================================
//====================================================================================
//====================================================================================
int main(int argc, char **argv){
    command = argv[0];
    getopt_from_command(argc, argv);

    // 辐射因子
    double srcRadi[SRC_M_NUM][CHANNEL_NUM]={0};

    // 从标准输入中读取静态格林函数表
    double x0, y0, grn[SRC_M_NUM][CHANNEL_NUM]={0}, syn[CHANNEL_NUM]={0}, syn_upar[CHANNEL_NUM][CHANNEL_NUM]={0};
    double grn_uiz[SRC_M_NUM][CHANNEL_NUM]={0}, grn_uir[SRC_M_NUM][CHANNEL_NUM]={0};

    // 根据参数设置，选择分量名
    const char *chs = (rot2ZNE)? ZNEchs : ZRTchs;


    // 建立一个指针数组，方便读取多列数据
    const int max_ncol = 47;
    double *pt_grn[max_ncol];
    // 按照特定顺序
    {
        double **pt = &pt_grn[0];
        *(pt++) = &x0;
        *(pt++) = &y0;
        for(int m=0; m<3; ++m){
            for(int k=0; k<SRC_M_NUM; ++k){
                for(int c=0; c<CHANNEL_NUM; ++c){
                    if(SRC_M_ORDERS[k]==0 && ZRTchs[c] == 'T')  continue;

                    if(m==0){
                        *pt = &grn[k][c];
                    } else if(m==1){
                        *pt = &grn_uiz[k][c];
                    } else if(m==2){
                        *pt = &grn_uir[k][c];
                    }
                    pt++;
                }
            }
        }
    }

    // 是否已打印输出的列名
    bool printHead = false;

    // 输入列数
    int ncols = 0;

    // 方位角
    double azrad = 0.0;

    // 物性参数
    double src_va=0.0, src_vb=0.0, src_rho=0.0, src_mu=0.0;
    double rcv_va=0.0, rcv_vb=0.0, rcv_rho=0.0;

    // 用于计算位移空间导数的比例系数
    double upar_scale=1.0; 

    // 震中距
    double dist = 0.0;

    char line[1024];
    int iline = 0;
    while(fgets(line, sizeof(line), stdin)){
        iline++;
        if(iline == 1){
            // 读取震源物性参数
            if(3 != sscanf(line, "# %lf %lf %lf", &src_va, &src_vb, &src_rho)){
                fprintf(stderr, "[%s] " BOLD_RED "Error! Unable to read src property from \"%s\". \n" DEFAULT_RESTORE, command, line);
                exit(EXIT_FAILURE);
            }
            if(src_va <= 0.0 || src_vb < 0.0 || src_rho <= 0.0){
                fprintf(stderr, "[%s] " BOLD_RED "Error! Bad src_va, src_vb or src_rho from \"%s\".\n" DEFAULT_RESTORE, command, line);
                exit(EXIT_FAILURE);
            }
            if(src_vb == 0.0 && mult_src_mu){
                fprintf(stderr, "[%s] " BOLD_RED 
                    "Error! Zero src_vb from \"%s\". "
                    "Maybe you try to use -Su<scale> but the source is in the liquid. "
                    "Use -S<scale> instead.\n" 
                    DEFAULT_RESTORE, command, line);
                exit(EXIT_FAILURE);
                src_mu = src_vb*src_vb*src_rho*1e10;
            }

            if(mult_src_mu)  M0 *= src_mu;
        }
        else if(iline == 2){
            // 读取场点物性参数
            if(3 != sscanf(line, "# %lf %lf %lf", &rcv_va, &rcv_vb, &rcv_rho)){
                fprintf(stderr, "[%s] " BOLD_RED "Error! Unable to read rcv property from \"%s\". \n" DEFAULT_RESTORE, command, line);
                exit(EXIT_FAILURE);
            }
            if(rcv_va <= 0.0 || rcv_vb < 0.0 || rcv_rho <= 0.0){
                fprintf(stderr, "[%s] " BOLD_RED "Error! Bad rcv_va, rcv_vb or rcv_rho in line %d from \"%s\".\n" DEFAULT_RESTORE, command, iline, line);
                exit(EXIT_FAILURE);
            }
        }
        else if(iline == 3){
            // 根据列长度判断是否有位移空间导数
            char *copyline = strdup(line+1);  // +1去除首个#字符
            char *token = strtok(copyline, " ");
            while (token != NULL) {
                ncols++;
                token = strtok(NULL, " ");
            }
            free(copyline);

            // 想合成位移空间导数但输入的格林函数没有
            if(calc_upar && ncols < max_ncol){
                fprintf(stderr, "[%s] " BOLD_RED "Error! The input has no spatial derivatives. \n" DEFAULT_RESTORE, command);
                exit(EXIT_FAILURE);
            }
        }
        if(line[0] == '#')  continue;

        // 读取该行数据
        char *copyline = strdup(line);
        char *token = strtok(copyline, " ");
        for(int i=0; i<ncols; ++i){
            sscanf(token, "%lf", pt_grn[i]);  token = strtok(NULL, " ");
        }
        free(copyline);

        // 计算方位角
        azrad = atan2(y0, x0);

        // 计算震中距
        dist = sqrt(x0*x0 + y0*y0);
        if(dist < 1e-5)  dist=1e-5;

        // 先输出列名
        if(!printHead){
            // 打印物性参数
            fprintf(stdout, "# "GRT_REAL_FMT" "GRT_REAL_FMT" "GRT_REAL_FMT"\n", src_va, src_vb, src_rho);
            fprintf(stdout, "# "GRT_REAL_FMT" "GRT_REAL_FMT" "GRT_REAL_FMT"\n", rcv_va, rcv_vb, rcv_rho);
            
            char XX[20];
            sprintf(XX, GRT_STRING_FMT, "X(km)"); XX[0]='#';
            fprintf(stdout, "%s", XX);
            fprintf(stdout, GRT_STRING_FMT, "Y(km)");
            char s_channel[5];
            for(int i=0; i<CHANNEL_NUM; ++i){
                sprintf(s_channel, "%s%c", s_computeType, toupper(chs[i])); 
                fprintf(stdout, GRT_STRING_FMT, s_channel);
            }

            if(calc_upar){
                for(int k=0; k<CHANNEL_NUM; ++k){
                    for(int i=0; i<CHANNEL_NUM; ++i){
                        sprintf(s_channel, "%c%s%c", tolower(chs[k]), s_computeType, toupper(chs[i])); 
                        fprintf(stdout, GRT_STRING_FMT, s_channel);
                    }
                }
            }

            fprintf(stdout, "\n");
            printHead = true;
        }

        double (*grn3)[CHANNEL_NUM];  // 使用对应类型的格林函数
        double tmpsyn[CHANNEL_NUM];
        for(int ityp=0; ityp<calcUTypes; ++ityp){
            // 求位移空间导数时，需调整比例系数
            switch (ityp){
                // 合成位移
                case 0:
                    upar_scale=1.0;
                    break;

                // 合成ui_z
                case 1:
                // 合成ui_r
                case 2:
                    upar_scale=1e-5;
                    break;

                // 合成ui_t
                case 3:
                    upar_scale=1e-5 / dist;
                    break;
                    
                default:
                    break;
            }

            if(ityp==1){
                grn3 = grn_uiz;
            } else if(ityp==2){
                grn3 = grn_uir;
            } else {
                grn3 = grn;
            }

            tmpsyn[0] = tmpsyn[1] = tmpsyn[2] = 0.0;
            // 计算震源辐射因子
            set_source_radiation(srcRadi, computeType, ityp==3, M0, upar_scale, azrad, mchn);

            for(int i=0; i<CHANNEL_NUM; ++i){
                for(int k=0; k<SRC_M_NUM; ++k){
                    tmpsyn[i] += grn3[k][i] * srcRadi[k][i];
                }
            }

            // 保存数据
            for(int i=0; i<CHANNEL_NUM; ++i){
                if(ityp == 0){
                    syn[i] = tmpsyn[i];
                } else {
                    syn_upar[ityp-1][i] = tmpsyn[i];
                }
            }
        }

        // 是否要转到ZNE
        if(rot2ZNE){
            if(calc_upar){
                rot_zrt2zxy_upar(azrad, syn, syn_upar, dist*1e5);
            } else {
                rot_zxy2zrt_vec(-azrad, syn);
            }
        }

        // 输出数据
        fprintf(stdout, GRT_REAL_FMT GRT_REAL_FMT, x0, y0);
        for(int i=0; i<CHANNEL_NUM; ++i){
            fprintf(stdout, GRT_REAL_FMT, syn[i]);
        }
        if(calc_upar){
            for(int i=0; i<CHANNEL_NUM; ++i){
                for(int k=0; k<CHANNEL_NUM; ++k){
                    fprintf(stdout, GRT_REAL_FMT, syn_upar[i][k]);
                }
            }
        }
        

        fprintf(stdout, "\n");
        
    }

    if(iline==0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Empty input. \n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
}