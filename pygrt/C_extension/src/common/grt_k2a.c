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
    MYREAL k;
    MYCOMPLEX res[4];
    const int nsrc = 6, ncols = 15;
    int nums[] = {2,2,3,2,3,3};
    const char *colnames[] = {
        "EXP_q0", "EXP_w0", "VF_q0", "VF_w0",
        "HF_q1", "HF_w1", "HF_v1",
        "DC_q0", "DC_w0",
        "DC_q1", "DC_w1", "DC_v1",
        "DC_q2", "DC_w2", "DC_v2"
    };

    // 先输出列名
    {
        char K[20];
        sprintf(K, GRT_STRING_FMT, "k");  K[0]='#';
        fprintf(stdout, "%s", K);
        for(int i=0; i<ncols; ++i){
            fprintf(stdout, GRT_STR_CMPLX_FMT, colnames[i]);
        }
        fprintf(stdout, "\n");
    }
    
    bool fullrow=false;
    // 读取数据    
    while (true) {
        fullrow=false;

        if(1 != fread(&k, sizeof(MYREAL), 1, fp)) break;
        fprintf(stdout, GRT_REAL_FMT, k);

        for(int i=0; i<nsrc; ++i){
            if(nums[i] != fread(res, sizeof(MYCOMPLEX), nums[i], fp)) break;
            for(int m=0; m<nums[i]; ++m){
                fprintf(stdout, GRT_CMPLX_FMT, CREAL(res[m]), CIMAG(res[m]));
            }
            if(i==nsrc-1)  fullrow=true;
        }

        fprintf(stdout, "\n");
        
        // 是否读完一整行
        if(! fullrow)  break;
    }

}

/**
 * 处理峰谷平均法的过程文件
 * 
 * @param     fp       文件指针
 */
static void print_PTAM(FILE *fp){
    MYREAL k;
    MYCOMPLEX res;
    const int ntyps = 18;
    const char *colnames[] = {
        "EXP_00", "VF_00", "DC_00", 
        "EXP_02", "VF_02", "DC_02", 
        "HF_10", "DC_10", 
        "HF_11", "DC_11", 
        "HF_12", "DC_12", 
        "HF_13", "DC_13",
        "DC_20", "DC_21", "DC_22", "DC_23"
    };

    // 先输出列名
    {
        char K[20], K2[20];
        sprintf(K2, "sum_%s_k", colnames[0]);
        sprintf(K, GRT_STRING_FMT, K2);  K[0]='#';
        fprintf(stdout, "%s", K);
        sprintf(K2, "sum_%s", colnames[0]);
        fprintf(stdout, GRT_STR_CMPLX_FMT, K2);
        for(int i=1; i<ntyps; ++i){
            sprintf(K2, "sum_%s_k", colnames[i]);
            fprintf(stdout, GRT_STRING_FMT, K2);
            sprintf(K2, "sum_%s", colnames[i]);
            fprintf(stdout, GRT_STR_CMPLX_FMT, K2);
        }
        fprintf(stdout, "\n");
    }
    
    bool fullrow=false;
    // 读取数据    
    while (true) {
        fullrow=false;

        for(int i=0; i<ntyps; ++i){
            if(1 != fread(&k, sizeof(MYREAL), 1, fp)) break;
            fprintf(stdout, GRT_REAL_FMT, k);

            if(1 != fread(&res, sizeof(MYCOMPLEX), 1, fp)) break;
            fprintf(stdout, GRT_CMPLX_FMT, CREAL(res), CIMAG(res));
            
            if(i==ntyps-1)  fullrow=true;
        }

        fprintf(stdout, "\n");
        
        // 是否读完一整行
        if(! fullrow)  break;
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