/**
 * @file   grt_b2a.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-27
 * 
 *    一个简单的小程序，将二进制SAC文件中的波形文件转为方便可读的文本文件，
 *    可供没有安装SAC程序和不使用Python的用户临时使用。
 * 
 */


#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "common/sacio2.h"
#include "common/logo.h"
#include "common/colorstr.h"

extern char *optarg;
extern int optind;
extern int optopt;

/**
 * 打印使用说明
 */
static void print_help(){
print_logo();
printf("\n"
"[grt.b2a]\n\n"
"    Convert a binary SAC file into an ASCII file, \n"
"    write to standard output (ignore header vars).\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt.b2a <sacfile>\n"
"\n\n\n"
);
}


int main(int argc, char **argv){
    const char *command = argv[0];

    // 输入不够
    if(argc < 2){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Need set a SAC file. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }

    // 输入过多
    if(argc > 2){
        fprintf(stderr, "[%s] " BOLD_RED "Error! You should set only one SAC file. Use '-h' for help.\n" DEFAULT_RESTORE, command);
        exit(EXIT_FAILURE);
    }

    // 使用-h查看帮助
    if(strcmp(argv[1], "-h") == 0){
        print_help();
        exit(EXIT_SUCCESS);
    }

    const char *filepath = argv[1];
    // 检查文件名是否存在
    if(access(filepath, F_OK) == -1){
        fprintf(stderr, "[%s] " BOLD_RED "Error! %s not exists.\n" DEFAULT_RESTORE, command, filepath);
        exit(EXIT_FAILURE);
    }


    // 读入SAC文件
    SACHEAD hd;
    float *arr = read_SAC(command, filepath, &hd, NULL);

    // 将波形写入标准输出，第一列时间，第二列振幅
    float begt = hd.b;
    float dt = hd.delta;
    int npts = hd.npts;
    for(int i=0; i<npts; ++i){
        printf("%13.7e  %13.7e\n", begt+dt*i, arr[i]);
    }

    free(arr);
}