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
const GRT_MODULE_ENTRY GRT_Modules_Entry[] = {
    #define X(name) {#name, name##_main},
        GRT_Module_List
    #undef X
    {NULL, NULL} // 结束标记
};

/** 定义包含子模块名称的字符串数组 */
const char *GRT_Module_Names[] = {
    #define X(name) #name ,
        GRT_Module_List
    #undef X
    NULL
};


/** 参数控制结构体 */
typedef struct {
    int dummy;
} GRT_MAIN_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MAIN_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
grt_print_logo();
printf("\n"
"Usage: \n"
"----------------------------------------------------------------\n"
"    grt [options]\n"
"    grt <module-name> [<module-options>] ...\n\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -v            Display the program version.\n"
"\n"
"    -h            Display this help message.\n"
"\n\n");
printf("GRT supports the following modules:\n"
"----------------------------------------------------------------\n");
for (size_t n = 0; GRT_Module_Names[n] != NULL; ++n) {
    const char *name = GRT_Module_Names[n];
    printf("    %-s\n", name);
}
printf("\n"
"For each module, you can use -h to see its help message, e.g.\n"
"    grt greenfn -h \n"
"\n");
}

/** 从命令行中读取选项，处理后记录到参数控制结构体 */
static void getopt_from_command(GRT_MAIN_CTRL *Ctrl, int argc, char **argv){
    (void)Ctrl;
    int opt;
    while ((opt = getopt(argc, argv, ":vh")) != -1) {
        switch (opt) {
            // 打印版本
            case 'v':
                printf(GRT_MAIN_COMMAND" %s\n", GRT_VERSION);
                exit(EXIT_SUCCESS);
                break;

            GRT_Common_Options_in_Switch(optopt);
        }
    }

    // 必须有输入
    GRTCheckOptionSet(argc > 1);
}


/** 查找并执行子模块 */
int dispatch_command(GRT_MAIN_CTRL *Ctrl, int argc, char **argv) {
    (void)Ctrl;
    char *entry_name = strdup(argv[1]);

    // 是否单独传入“static”以计算静态解
    bool is_single_static = false;
    if(strcmp(argv[1], "static") == 0 && argc > 2){
        is_single_static = true;
        GRT_SAFE_ASPRINTF(&entry_name, "static_%s", argv[2]);

        // 同理也“修改” argv[2] 参数
        char *newarg = NULL;
        GRT_SAFE_ASPRINTF(&newarg, "static_%s", argv[2]);
        
        // 原本指向系统分配的只读内存的指针被覆盖，故后续需要手动free
        argv[2] = newarg;
    }
    
    int return_code = EXIT_SUCCESS;
    bool valid_entry_name = false;
    for (const GRT_MODULE_ENTRY *entry = GRT_Modules_Entry; entry->name != NULL; entry++) {
        if (strcmp(entry_name, entry->name) == 0) {
            GRT_MODULE_NAME = entry_name;
            return_code = entry->func(argc - 1 - (int)is_single_static, argv + 1 + (int)is_single_static);
            valid_entry_name = true;
            break;
        }
    }
    
    // 未知子模块
    if( ! valid_entry_name){
        GRTRaiseError("Unknown module %s. Use \"-h\" for help.\n", entry_name);
    }

    if(is_single_static)  GRT_SAFE_FREE_PTR(argv[2]);
    GRT_SAFE_FREE_PTR(entry_name);
    return return_code;
}


/** 主函数 */
int main(int argc, char **argv) {
    GRT_MODULE_NAME = GRT_MAIN_COMMAND;
    GRT_MAIN_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    if(argc <= 2){
        getopt_from_command(Ctrl, argc, argv);
    }
    
    dispatch_command(Ctrl, argc, argv);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}