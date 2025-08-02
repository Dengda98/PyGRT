/**
 * @file   grt.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * C 程序的主函数，由此发起各个子模块的任务
 * 
 */


#include "grt.h"


typedef struct {
    char *name;
    
} GRT_MAIN_CTRL;


/**
 * 打印使用说明
 */
static void print_help(){
print_logo();
printf("\n"
"Usage: \n"
"----------------------------------------------------------------\n"
"    grt <submodule name> [options] ...\n\n\n");
printf("GRT supports the following submodules:\n"
"----------------------------------------------------------------\n");
for (MYINT n = 0; GRT_Submodule_Names[n] != NULL; ++n) {
    const char *name = GRT_Submodule_Names[n];
    printf("    %-s\n", name);
}
printf("\n\n");
printf("For each submodule, you can use -h to see the help message.\n\n");
}

/**
 * 从命令行中读取选项，处理后记录到参数控制结构体
 * 
 * @param     argc      命令行的参数个数
 * @param     argv      多个参数字符串指针
 */
static void getopt_from_command(GRT_MAIN_CTRL *Ctrl, int argc, char **argv){
    int opt;
    while ((opt = getopt(argc, argv, ":h")) != -1) {
        switch (opt) {
            GRT_Common_Options_in_Switch(Ctrl, optopt);
        }
    }

    // 必须有输入
    GRTCheckOptionEmpty(Ctrl, argc-1);
}


// 查找并执行子命令
int dispatch_command(int argc, char **argv) {
    char *entry_name = argv[1];
    
    for (GRT_SUBMODULE_ENTRY *entry = GRT_Submodules_Entry; entry->name != NULL; entry++) {
        if (strcmp(entry_name, entry->name) == 0) {
            return entry->func(argc - 1, argv + 1);
        }
    }
    

    // 未知子模块
    GRTRaiseError("Error! Unknown submodule %s. Use \"-h\" for help.\n", entry_name);
}


// 主函数
int main(int argc, char **argv) {
    GRT_MAIN_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    Ctrl->name = argv[0];

    if(argc == 2){
        getopt_from_command(Ctrl, argc, argv);
    }

    return dispatch_command(argc, argv);
}