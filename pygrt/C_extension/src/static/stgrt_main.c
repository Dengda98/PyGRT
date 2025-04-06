/**
 * @file   grt_static.c
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
// 输出目录
static char *s_output_dir = NULL;
// 波数积分间隔
static double Length=0.0;
// 波数积分相关变量
static double keps=-1.0, k0=5.0;
// 震中距数组以及保存对应初至波走时的数组
static char **s_rs = NULL;
static MYREAL *rs = NULL;
static int nr=0;
// 是否silence整个输出
static bool silenceInput=false;
// 计算哪些格林函数，确定震源类型, 默认计算全部
static bool doEXP=true, doVF=true, doHF=true, doDC=true;

// 是否计算位移空间导数
static bool calc_upar=false;

// 各选项的标志变量，初始化为0，定义了则为1
static int M_flag=0, D_flag=0, 
            O_flag=0, 
            L_flag=0, 
            K_flag=0, s_flag=0, 
            R_flag=0, 
            G_flag=0, e_flag=0;

// 三分量代号
const char chs[3] = {'Z', 'R', 'T'};

/**
 * 打印使用说明
 */
static void print_help(){
print_logo();
printf("\n"
"[grt.syn]\n\n"
"    Compute static Green's Functions.\n"
"\n"
"The units of output Green's Functions for different sources are: \n"
"    + Explosion:     1e-20 cm/(dyne-cm)\n"
"    + Single Force:  1e-15 cm/(dyne)\n"
"    + Double Couple: 1e-20 cm/(dyne-cm)\n"
"    + Moment Tensor: 1e-20 cm/(dyne-cm)\n" 
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt.static -M<model> -D<depsrc>/<deprcv> \n"
"        -R<r1>,<r2>[,...]    [-O<outdir>] \n"
"        [-L<length>]  \n" 
"        [-K<k0>[/<iwk0>/<ampk>/<keps>]]  \n"
"        [-G<b1>[/<b2>/<b3>/<b4>]] [-S<i1>,<i2>[,...]] [-e]\n"
"        [-s]\n"
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
    while ((opt = getopt(argc, argv, ":M:D:O:L:E:K:shR:G:e")) != -1) {
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

            // 输出路径 -Ooutput_dir
            case 'O':
                O_flag = 1;
                s_output_dir = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                strcpy(s_output_dir, optarg);
                break;

            // 波数积分间隔 -LLength
            case 'L':
                L_flag = 1;
                if(0 == sscanf(optarg, "%lf", &Length)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -L.\n" DEFAULT_RESTORE, command);
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

            // 震中距数组，-Rr1,r2,r3,r4 ...
            case 'R':
                R_flag = 1;
                {
                    char *token;
                    char *str_copy = strdup(optarg);  // 创建字符串副本，以免修改原始字符串
                    token = strtok(str_copy, ",");

                    while(token != NULL){
                        s_rs = (char**)realloc(s_rs, sizeof(char*)*(nr+1));
                        s_rs[nr] = NULL;
                        s_rs[nr] = (char*)realloc(s_rs[nr], sizeof(char)*(strlen(token)+1));
                        rs = (MYREAL*)realloc(rs, sizeof(MYREAL)*(nr+1));
                        strcpy(s_rs[nr], token);
                        rs[nr] = atof(token);
                        if(rs[nr] == 0.0){
                            fprintf(stderr, "[%s] " BOLD_RED "Warning! Add 1e-5 to Zero epicentral distance in -R.\n" DEFAULT_RESTORE, command);
                            rs[nr] += 1e-5;
                        }
                        if(rs[nr] < 0.0){
                            fprintf(stderr, "[%s] " BOLD_RED "Error! Can't set negative epicentral distance(%f) in -R.\n" DEFAULT_RESTORE, command, rs[nr]);
                            exit(EXIT_FAILURE);
                        }


                        token = strtok(NULL, ",");
                        nr++;
                    }
                    free(str_copy);
                }
                break;

            // 选择要计算的格林函数 -G1/1/1/1
            case 'G': 
                G_flag = 1;
                doEXP = doVF = doHF = doDC = false;
                {
                    int i1, i2, i3, i4;
                    i1 = i2 = i3 = i4 = 0;
                    if(0 == sscanf(optarg, "%d/%d/%d/%d", &i1, &i2, &i3, &i4)){
                        fprintf(stderr, "[%s] " BOLD_RED "Error in -G.\n" DEFAULT_RESTORE, command);
                        exit(EXIT_FAILURE);
                    };
                    doEXP = (i1!=0);
                    doVF  = (i2!=0);
                    doHF  = (i3!=0);
                    doDC  = (i4!=0);
                }
                
                // 至少要有一个真
                if(!(doEXP || doVF || doHF || doDC)){
                    fprintf(stderr, "[%s] " BOLD_RED "Error! At least set one true value in -G.\n" DEFAULT_RESTORE, command);
                    exit(EXIT_FAILURE);
                }

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
    if(R_flag == 0){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set -R. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }

    // 建立保存目录
    if(mkdir(s_output_dir, 0777) != 0){
        if(errno != EEXIST){
            fprintf(stderr, "[%s] " BOLD_RED "Error! Unable to create folder %s. Error code: %d\n" DEFAULT_RESTORE, command, s_output_dir, errno);
            exit(EXIT_FAILURE);
        }
    }

    // 在目录中保留命令
    char *dummy = (char*)malloc(sizeof(char)*(strlen(s_output_dir)+100));
    sprintf(dummy, "%s/command", s_output_dir);
    FILE *fp = fopen(dummy, "a");
    for(int i=0; i<argc; ++i){
        fprintf(fp, "%s ", argv[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
    free(dummy);
}





int main(int argc, char **argv){
    command = argv[0];

    // 传入参数 
    getopt_from_command(argc, argv);

    // 读入模型文件
    if((pymod = read_pymod_from_file(command, s_modelpath, depsrc, deprcv)) ==NULL){
        exit(EXIT_FAILURE);
    }
    
    // 设置积分间隔默认值
    if(Length == 0.0)  Length = 15.0;

    // 建立格林函数的complex数组
    MYREAL (*EXPgrn)[2] = (doEXP) ? (MYREAL(*)[2])calloc(nr, sizeof(*EXPgrn)) : NULL;
    MYREAL (*VFgrn)[2]  = (doVF)  ? (MYREAL(*)[2])calloc(nr, sizeof(*VFgrn))  : NULL;
    MYREAL (*HFgrn)[3]  = (doHF)  ? (MYREAL(*)[3])calloc(nr, sizeof(*HFgrn))  : NULL;
    MYREAL (*DDgrn)[2]  = (doDC)  ? (MYREAL(*)[2])calloc(nr, sizeof(*DDgrn))  : NULL;
    MYREAL (*DSgrn)[3]  = (doDC)  ? (MYREAL(*)[3])calloc(nr, sizeof(*DSgrn))  : NULL;
    MYREAL (*SSgrn)[3]  = (doDC)  ? (MYREAL(*)[3])calloc(nr, sizeof(*SSgrn))  : NULL;


    //==============================================================================
    // 计算静态格林函数
    integ_static_grn(
        pymod, nr, rs, keps, k0, Length, 
        EXPgrn, VFgrn, HFgrn, DDgrn, DSgrn, SSgrn
    );
    //==============================================================================


    // 输出
    for(int ir=0; ir<nr; ++ir){
        if(doEXP){
            fprintf(stdout, "%15.7e", EXPgrn[ir][0]);
            fprintf(stdout, "%15.7e", EXPgrn[ir][1]);
        }
        if(doVF){
            fprintf(stdout, "%15.7e", VFgrn[ir][0]);
            fprintf(stdout, "%15.7e", VFgrn[ir][1]);
        }
        if(doHF){
            fprintf(stdout, "%15.7e", HFgrn[ir][0]);
            fprintf(stdout, "%15.7e", HFgrn[ir][1]);
            fprintf(stdout, "%15.7e", HFgrn[ir][2]);
        }
        if(doDC){
            fprintf(stdout, "%15.7e", DDgrn[ir][0]);
            fprintf(stdout, "%15.7e", DDgrn[ir][1]);
            fprintf(stdout, "%15.7e", DSgrn[ir][0]);
            fprintf(stdout, "%15.7e", DSgrn[ir][1]);
            fprintf(stdout, "%15.7e", DSgrn[ir][2]);
            fprintf(stdout, "%15.7e", SSgrn[ir][0]);
            fprintf(stdout, "%15.7e", SSgrn[ir][1]);
            fprintf(stdout, "%15.7e", SSgrn[ir][2]);
        }

        fprintf(stdout, "\n");
    }
}

