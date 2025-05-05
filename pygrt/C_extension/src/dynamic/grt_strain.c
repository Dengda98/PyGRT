/**
 * @file   grt_strain.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-28
 * 
 *    根据预先合成的位移空间导数，组合成应变张量
 * 
 */


#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>

#include "common/sacio2.h"
#include "common/const.h"
#include "common/logo.h"
#include "common/colorstr.h"


//****************** 在该文件以内的全局变量 ***********************//
// 命令名称
static char *command = NULL;

// 输出分量格式，即是否需要旋转到ZNE
static bool rot2ZNE = false;

// 三分量
const char *chs = NULL;


/**
 * 打印使用说明
 */
static void print_help(){
print_logo();
printf("\n"
"[grt.strain]\n\n"
"    Conbine spatial derivatives of displacements into strain tensor.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt.strain <syn_dir>/<name>\n"
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
    while ((opt = getopt(argc, argv, ":h")) != -1) {
        switch (opt) {

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
    if(argc != 2){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set options. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }
}



int main(int argc, char **argv){
    command = argv[0];

    getopt_from_command(argc, argv);

    
    // 合成地震图目录路径
    char *s_synpath = (char*)malloc(sizeof(char)*(strlen(argv[1])+1));
    // 保存文件前缀 
    char *s_prefix = (char*)malloc(sizeof(char)*(strlen(argv[1])+1));
    if(2 != sscanf(argv[1], "%[^/]/%s", s_synpath, s_prefix)){
        fprintf(stderr, "[%s] " BOLD_RED "Error format in \"%s\".\n" DEFAULT_RESTORE, command, argv[1]);
        exit(EXIT_FAILURE);
    }

    // 检查是否存在该目录
    DIR *dir = opendir(s_synpath);
    if (dir == NULL) {
        fprintf(stderr, "[%s] " BOLD_RED "Error! Directory \"%s\" not exists.\n" DEFAULT_RESTORE, command, s_synpath);
        exit(EXIT_FAILURE);
    } 


    // ----------------------------------------------------------------------------------
    // 开始读取计算，输出6个量
    float *arrin = NULL;
    char c1, c2;
    char *s_filepath = (char*)malloc(sizeof(char) * (strlen(s_synpath)+strlen(s_prefix)+100));

    // 判断标志性文件是否存在，来判断输出使用ZNE还是ZRT
    sprintf(s_filepath, "%s/n%sN.sac", s_synpath, s_prefix);
    rot2ZNE = (access(s_filepath, F_OK) == 0);

    // 指示特定的通道名
    chs = (rot2ZNE)? ZNEchs : ZRTchs;


    // 读取一个头段变量，获得基本参数，分配数组内存
    SACHEAD hd;
    sprintf(s_filepath, "%s/%c%s%c.sac", s_synpath, tolower(chs[0]), s_prefix, chs[0]);
    read_SAC_HEAD(command, s_filepath, &hd);
    int npts=hd.npts;
    float dist=hd.dist;
    float *arrout = (float*)calloc(npts, sizeof(float));

    // ----------------------------------------------------------------------------------
    // 循环6个分量
    for(int i1=0; i1<3; ++i1){
        c1 = chs[i1];
        for(int i2=i1; i2<3; ++i2){
            c2 = chs[i2];

            // 读取数据 u_{i,j}
            sprintf(s_filepath, "%s/%c%s%c.sac", s_synpath, tolower(c2), s_prefix, c1);
            arrin = read_SAC(command, s_filepath, &hd, arrin);

            // 累加
            for(int i=0; i<npts; ++i)  arrout[i] += arrin[i];

            // 读取数据 u_{j,i}
            sprintf(s_filepath, "%s/%c%s%c.sac", s_synpath, tolower(c1), s_prefix, c2);
            arrin = read_SAC(command, s_filepath, &hd, arrin);

            // 累加
            for(int i=0; i<npts; ++i)  arrout[i] = (arrout[i] + arrin[i]) * 0.5f;

            // 特殊情况需加上协变导数，1e-5是因为km->cm
            if(c1=='R' && c2=='T'){
                // 读取数据 u_T
                sprintf(s_filepath, "%s/%sT.sac", s_synpath, s_prefix);
                arrin = read_SAC(command, s_filepath, &hd, arrin);
                for(int i=0; i<npts; ++i)  arrout[i] -= 0.5f * arrin[i] / dist * 1e-5;
            }
            else if(c1=='T' && c2=='T'){
                // 读取数据 u_R
                sprintf(s_filepath, "%s/%sR.sac", s_synpath, s_prefix);
                arrin = read_SAC(command, s_filepath, &hd, arrin);
                for(int i=0; i<npts; ++i)  arrout[i] += arrin[i] / dist * 1e-5;
            }

            // 保存到SAC
            sprintf(hd.kcmpnm, "%c%c", c1, c2);
            sprintf(s_filepath, "%s/%s.strain.%c%c.sac", s_synpath, s_prefix, c1, c2);
            write_sac(s_filepath, hd, arrout);

            // 置零
            for(int i=0; i<npts; ++i)  arrout[i] = 0.0f;
        }
    }

    if(arrin)   free(arrin);
    if(arrout)  free(arrout);
    free(s_filepath);
    free(s_synpath);
    free(s_prefix);
}
