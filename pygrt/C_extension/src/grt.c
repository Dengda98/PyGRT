/**
 * @file   grt.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * C 程序的主函数，由此发起各个子模块的任务
 * 
 */


#include "grt.h"

/** 注册所有子模块命令 */
const GRT_SUBMODULE_ENTRY GRT_Submodules_Entry[] = {
    #define X(name) {#name, name##_main},
        GRT_Submodule_List
    #undef X
    {NULL, NULL} // 结束标记
};

/** 定义包含子模块名称的字符串数组 */
const char *GRT_Submodule_Names[] = {
    #define X(name) #name ,
        GRT_Submodule_List
    #undef X
    NULL
};


/** 参数控制结构体 */
typedef struct {
    char *name;
    
} GRT_MAIN_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MAIN_CTRL *Ctrl){
    free(Ctrl->name);
    free(Ctrl);
}

/** 打印使用说明 */
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

/** 从命令行中读取选项，处理后记录到参数控制结构体 */
static void getopt_from_command(GRT_MAIN_CTRL *Ctrl, int argc, char **argv){
    char* command = Ctrl->name;
    int opt;
    while ((opt = getopt(argc, argv, ":h")) != -1) {
        switch (opt) {
            GRT_Common_Options_in_Switch(command, optopt);
        }
    }

    // 必须有输入
    GRTCheckOptionSet(command, argc > 1);
}


/** 查找并执行子模块 */
int dispatch_command(GRT_MAIN_CTRL *Ctrl, int argc, char **argv) {
    char *entry_name = (char*)malloc((strlen(argv[1]) + ((argc > 2)? strlen(argv[2]) : 0) + 100)*sizeof(char));
    strcpy(entry_name, argv[1]);

    // 计算静态解
    bool is_static = false;
    if(strcmp(argv[1], "static") == 0 && argc > 2){
        is_static = true;
        sprintf(entry_name, "static_%s", argv[2]);
    }
    
    for (const GRT_SUBMODULE_ENTRY *entry = GRT_Submodules_Entry; entry->name != NULL; entry++) {
        if (strcmp(entry_name, entry->name) == 0) {
            free(entry_name);
            return entry->func(argc - 1 - (int)is_static, argv + 1 + (int)is_static);
        }
    }

    // 未知子模块
    GRTRaiseError("[%s] Error! Unknown submodule %s. Use \"-h\" for help.\n", Ctrl->name, entry_name);
}


/** 主函数 */
int main(int argc, char **argv) {
    GRT_MAIN_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    Ctrl->name = strdup(argv[0]);

    if(argc <= 2){
        getopt_from_command(Ctrl, argc, argv);
    }
    
    dispatch_command(Ctrl, argc, argv);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}