/**
 * @file   grt_sac2asc.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-27
 * 
 *    一个简单的小程序，将二进制SAC文件中的波形文件转为方便可读的文本文件，
 *    可供没有安装SAC程序和不使用Python的用户临时使用。
 * 
 */


#include "grt/common/sacio2.h"

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
"[grt sac2asc] %s\n\n", GRT_VERSION);printf(
"    Convert a binary SAC file into an ASCII file, \n"
"    write to standard output (ignore header vars).\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt sac2asc <sacfile>\n"
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


/** 子模块主函数 */
int sac2asc_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    getopt_from_command(Ctrl, argc, argv);

    Ctrl->s_filepath = strdup(argv[1]);
    
    // 检查文件名是否存在
    GRTCheckFileExist(Ctrl->s_filepath);

    // 读入SAC文件
    SACTRACE *insac = grt_read_SACTRACE(Ctrl->s_filepath, false);

    // 将波形写入标准输出，第一列时间，第二列振幅
    float begt = insac->hd.b;
    float dt = insac->hd.delta;
    int npts = insac->hd.npts;
    for(int i=0; i<npts; ++i){
        printf("%13.7e  %13.7e\n", begt+dt*i, insac->data[i]);
    }

    grt_free_SACTRACE(insac);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}