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

//****************** 在该文件以内的全局变量 ***********************//
// 命令名称
static char *command = NULL;

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