/**
 * @file   grt_k2a.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-27
 * 
 *    一个简单的小程序，将波数积分过程中输出的二进制过程文件转为方便可读的文本文件，
 *    这可以作为临时查看，但更推荐使用Python读取
 * 
 */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "common/const.h"
#include "common/logo.h"
#include "common/colorstr.h"
#include "common/iostats.h"


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
"[grt.k2a]\n\n"
"    Convert a binary stats file generated during wavenumber integration\n"
"    into an ASCII file, write to standard output.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt.k2a <statsfile>\n"
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
static const char* get_basename(const char* path) {
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


/**
 * 处理传统离散波数积分以及Filon积分的过程文件
 * 
 * @param     fp       文件指针
 */
static void print_K(FILE *fp){
    // 打印标题
    extract_stats(NULL, stdout);
    fprintf(stdout, "\n");
    
    // 读取数据    
    while (true) {
        if(0 != extract_stats(fp, stdout))  break;

        fprintf(stdout, "\n");
    }
}

/**
 * 处理峰谷平均法的过程文件
 * 
 * @param     fp       文件指针
 */
static void print_PTAM(FILE *fp){
    // 打印标题
    extract_stats_ptam(NULL, stdout);
    fprintf(stdout, "\n");
    
    // 读取数据    
    while (true) {
        if(0 != extract_stats_ptam(fp, stdout))  break;

        fprintf(stdout, "\n");
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


    // 打开stats
    FILE *fp=NULL;
    if((fp = fopen(filepath, "rb")) == NULL){
        fprintf(stderr, "[%s] " BOLD_RED "Error! Can't read %s.\n" DEFAULT_RESTORE, command, filepath);
        exit(EXIT_FAILURE);
    }

    // 根据文件名确定函数
    const char *basename = get_basename(filepath);
    if(strncmp(basename, "PTAM", 4) == 0) {
        print_PTAM(fp);
    } else if(strncmp(basename, "K", 1) == 0) {
        print_K(fp);
    } else {
        fprintf(stderr, "[%s] " BOLD_RED "Error! Can't read %s.\n" DEFAULT_RESTORE, command, filepath);
        exit(EXIT_FAILURE);
    }

    // 检查是否是因为文件结束而退出
    if (ferror(fp)) {
        fprintf(stderr, "[%s] " BOLD_RED "Error reading file %s.\n" DEFAULT_RESTORE, command, filepath);
        exit(EXIT_FAILURE);
    }

    fclose(fp);
}