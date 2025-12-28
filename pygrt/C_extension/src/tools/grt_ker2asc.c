/**
 * @file   grt_ker2asc.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-27
 * 
 *    一个简单的小程序，将波数积分过程中输出的二进制过程文件转为方便可读的文本文件，
 *    这可以作为临时查看，但更推荐使用Python读取
 * 
 */

#include "grt/common/const.h"
#include "grt/integral/iostats.h"
#include "grt/common/util.h"

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    /** 输入文件路径 */
    char *s_filepath;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl->s_filepath);
    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt ker2asc] %s\n\n", GRT_VERSION);printf(
"    Convert a binary stats file generated during wavenumber integration\n"
"    into an ASCII file, write to standard output.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt ker2asc <statsfile>\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    (void)Ctrl;
    int opt;
    while ((opt = getopt(argc, argv, ":h")) != -1) {
        switch (opt) {
            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 检查必选项有没有设置
    GRTCheckOptionSet(argc > 1);
}


/**
 * 处理传统离散波数积分以及Filon积分的过程文件
 * 
 * @param     fp       文件指针
 */
static void print_K(FILE *fp, const char *col0_name){
    // 打印标题
    grt_extract_stats(NULL, stdout, col0_name);
    fprintf(stdout, "\n");
    
    // 读取数据    
    while (true) {
        if(0 != grt_extract_stats(fp, stdout, col0_name))  break;

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
    grt_extract_stats_ptam(NULL, stdout);
    fprintf(stdout, "\n");
    
    // 读取数据    
    while (true) {
        if(0 != grt_extract_stats_ptam(fp, stdout))  break;

        fprintf(stdout, "\n");
    }

}


/** 子模块主函数 */
int ker2asc_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    getopt_from_command(Ctrl, argc, argv);

    Ctrl->s_filepath = strdup(argv[1]);

    // 检查文件名是否存在
    GRTCheckFileExist(Ctrl->s_filepath);

    // 打开stats
    FILE *fp = GRTCheckOpenFile(Ctrl->s_filepath, "rb");

    // 根据文件名确定函数
    const char *basename = grt_get_basename(Ctrl->s_filepath);
    if(strncmp(basename, "PTAM", 4) == 0) {
        print_PTAM(fp);
    } else if(strncmp(basename, "K", 1) == 0) {
        print_K(fp, "k");
    } else if(strncmp(basename, "C", 1) == 0) {
        print_K(fp, "c");
    } else {
        // 文件名不符合要求
        GRTRaiseError("Unsupported File \"%s\".\n", Ctrl->s_filepath);
    }

    // 未知错误导致文件提前结束
    if (ferror(fp)) {
        GRTRaiseError("An unknown error caused the premature termination "
                      "of reading the file \"%s\".\n", Ctrl->s_filepath);
    }
    fclose(fp);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}