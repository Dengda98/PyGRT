/**
 * @file   grt_travt.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12-2
 * 
 *    主程序，计算一维均匀半无限层状介质的初至走时
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


#include "travt/travt.h"
#include "common/const.h"
#include "common/model.h"
#include "common/logo.h"
#include "common/colorstr.h"

extern char *optarg;
extern int optind;
extern int optopt;


//****************** 在该文件以内的全局变量 ***********************//
static int M_flag=0, D_flag=0, R_flag=0;
// 命令名称
static char *command = NULL;
// 模型路径，模型PYMODEL1D指针，全局最大最小速度
static char *s_modelpath = NULL;
static PYMODEL1D *pymod;
// 震源和场点深度
static double depsrc, deprcv;
static char *s_depsrc = NULL, *s_deprcv = NULL;
// 震中距数组
static char **s_rs = NULL;
static MYREAL *rs = NULL;
static int nr=0;



/**
 * 打印使用说明
 */
static void print_help(){
print_logo();
printf("\n"
"[grt.travt]\n\n"
"    A Supplementary Tool of GRT to Compute First Arrival Traveltime\n"
"    of P-wave and S-wave in Horizontally Layerd Halfspace Model. \n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt.travt -M<model> -D<depsrc>/<deprcv> -R<r1>,<r2>[,...]\n"
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
"    -R<r1>,<r2>[,...]\n"
"                 Multiple epicentral distance (km), \n"
"                 seperated by comma.\n"
"\n"
"    -h           Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt.travt -Mmilrow -D2/0 -R10,20,30,40,50\n"
"\n\n\n"
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
    while ((opt = getopt(argc, argv, ":M:D:R:h")) != -1) {
        switch (opt) {
            // 模型路径，其中每行分别为 
            //      厚度(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  Qp   Qs
            // 互相用空格隔开即可
            case 'M':
                M_flag = 1;
                s_modelpath = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                strcpy(s_modelpath, optarg);
                // s_modelname = get_basename(s_modelpath);
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
                    fprintf(stderr, "[%s] " BOLD_RED "Error in -D.\n" DEFAULT_RESTORE, command);
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

}


int main(int argc, char **argv){
    command = argv[0];

    getopt_from_command(argc, argv);
    
    // 读入模型文件
    if((pymod = read_pymod_from_file(command, s_modelpath, depsrc, deprcv, true)) ==NULL){
        exit(EXIT_FAILURE);
    }
    // print_pymod(pymod);

    printf("------------------------------------------------\n");
    printf(" Distance(km)     Tp(secs)         Ts(secs)     \n");
    double travtP=-1, travtS=-1;
    for(int i=0; i<nr; ++i){
        travtP = compute_travt1d(
        pymod->Thk, pymod->Va, pymod->n, pymod->isrc, pymod->ircv, rs[i]);
        travtS = compute_travt1d(
        pymod->Thk, pymod->Vb, pymod->n, pymod->isrc, pymod->ircv, rs[i]);
        
        printf(" %-15s  %-15.3f  %-15.3f\n", s_rs[i], travtP, travtS);
    }
    printf("------------------------------------------------\n");

    
}