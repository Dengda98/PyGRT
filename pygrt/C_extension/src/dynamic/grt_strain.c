/**
 * @file   grt_strain.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-28
 * 
 *    根据预先合成的位移空间导数，组合成应变
 * 
 */


#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <ctype.h>
#include <string.h>

#include "common/sacio2.h"
#include "common/const.h"
#include "common/logo.h"
#include "common/colorstr.h"



/**
 * 打印使用说明
 */
static void print_help(){
print_logo();
printf("\n"
"[grt.strain]\n\n"
"    Conbine spatial derivatives of displacements into strain.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt.strain <syn_dir>/<name>\n"
"\n\n\n"
);
}



int main(int argc, char **argv){
    const char *command = argv[0];

    // 输入不够
    if(argc < 2){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set an input. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }

    // 输入过多
    if(argc > 2){
        fprintf(stderr, "[%s] " BOLD_RED "Error! You should set only one input. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }

    // 使用-h查看帮助
    if(strcmp(argv[1], "-h") == 0){
        print_help();
        exit(EXIT_SUCCESS);
    }

    
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
    const char chs[3] = {'R', 'T', 'Z'};

    // 读取一个头段变量，获得基本参数，分配数组内存
    SACHEAD hd;
    sprintf(s_filepath, "%s/r%sR.sac", s_synpath, s_prefix);
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
            sprintf(s_filepath, "%s/%c%s%c.sac", s_synpath, tolower(c1), s_prefix, c2);
            arrin = read_SAC(command, s_filepath, &hd, arrin);

            // 累加
            for(int i=0; i<npts; ++i)  arrout[i] += arrin[i];

            // 读取数据 u_{j,i}
            sprintf(s_filepath, "%s/%c%s%c.sac", s_synpath, tolower(c2), s_prefix, c1);
            arrin = read_SAC(command, s_filepath, &hd, arrin);

            // 累加
            for(int i=0; i<npts; ++i)  arrout[i] = (arrout[i] + arrin[i]) * 0.5f;

            // 特殊情况需加上协变导数
            if(c1=='R' && c2=='T'){
                // 读取数据 u_T
                sprintf(s_filepath, "%s/%sT.sac", s_synpath, s_prefix);
                arrin = read_SAC(command, s_filepath, &hd, arrin);
                for(int i=0; i<npts; ++i)  arrout[i] -= 0.5f * arrin[i] / dist;
            }
            else if(c1=='T' && c2=='T'){
                // 读取数据 u_R
                sprintf(s_filepath, "%s/%sR.sac", s_synpath, s_prefix);
                arrin = read_SAC(command, s_filepath, &hd, arrin);
                for(int i=0; i<npts; ++i)  arrout[i] += arrin[i] / dist;
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
