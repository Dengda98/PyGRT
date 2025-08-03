/**
 * @file   grt.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * C 程序的主函数，由此发起各个子模块的任务
 * 
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>

#include "common/const.h"
#include "common/logo.h"

// ------------------------------------------------------
/** X 宏，用于定义子模块命令。后续的命令名称和函数名称均与此匹配 */
#define GRT_Submodule_List   \
    /* dynamic */          \
    X(greenfn)             \
    X(syn)                 \
    X(rotation)            \
    X(strain)              \
    X(stress)              \
    /* static */           \
    X(static_greenfn)      \
    X(static_syn)          \
    X(static_rotation)     \
    X(static_strain)       \
    X(static_stress)       \
    /* other */            \
    X(k2a)                 \
    X(b2a)                 \
    X(travt)               \
// ------------------------------------------------------

/** 子模块的函数指针类型 */
typedef int (*GRT_SUBMODULE_FUNC)(int argc, char **argv);

/** 子模块命令注册结构 */
typedef struct {
    const char *name;
    GRT_SUBMODULE_FUNC func;
} GRT_SUBMODULE_ENTRY;


/** 声明所有子模块函数 */
#define X(name) int name##_main(int argc, char **argv);
    GRT_Submodule_List
#undef X

/** 注册所有子模块命令 */
extern const GRT_SUBMODULE_ENTRY GRT_Submodules_Entry[];

/** 定义包含子模块名称的字符串数组 */
extern const char *GRT_Submodule_Names[];





// GRT自定义报错信息
#define GRTRaiseError(ErrorMessage, ...) ({\
    fprintf(stderr, BOLD_RED ErrorMessage DEFAULT_RESTORE, ##__VA_ARGS__);\
    exit(EXIT_FAILURE);\
})

// GRT报错：选项设置不符要求
#define GRTBadOptionError(Ctrl, X, MoreErrorMessage, ...) ({\
    GRTRaiseError("[%s] Error in -"#X" . "MoreErrorMessage" Use \"-h\" for help.\n", Ctrl->name, ##__VA_ARGS__);\
})

// GRT报错：选项未设置参数    注意这里使用的是 %c 和 运行时变量X
#define GRTMissArgsError(Ctrl, X, MoreErrorMessage, ...) ({\
    GRTRaiseError("[%s] Error! Option \"-%c\" requires an argument. "MoreErrorMessage" Use \"-h\" for help.\n", Ctrl->name, X, ##__VA_ARGS__);\
})

// GRT报错：非法选项    注意这里使用的是 %c 和 运行时变量X
#define GRTInvalidOptionError(Ctrl, X, MoreErrorMessage, ...) ({\
    GRTRaiseError("[%s] Error! Option \"-%c\" is invalid. "MoreErrorMessage" Use \"-h\" for help.\n", Ctrl->name, X, ##__VA_ARGS__);\
})

// GRT报错：检查某个选项是否启用
#define GRTCheckOptionActive(Ctrl, X) ({\
    if( ! Ctrl->X.active){\
        GRTRaiseError("[%s] Error! Need set options -\""#X"\". Use \"-h\" for help.\n", Ctrl->name);\
    }\
})

// GRT报错：检查是否有使用选项
#define GRTCheckOptionEmpty(Ctrl, N) ({\
    if((N) == 0){\
        GRTRaiseError("[%s] Error! Need set options. Use \"-h\" for help.\n", Ctrl->name);\
    }\
})






/** 共有的命令行处理语句 */ 
#define GRT_Common_Options_in_Switch(Ctrl, X) \
    /** 帮助 */  \
    case 'h': \
        print_help(); \
        exit(EXIT_SUCCESS); \
        break; \
    /** 参数缺失 */  \
    case ':': \
        GRTMissArgsError(Ctrl, X, ""); \
        break; \
    /** 非法选项 */  \
    case '?': \
    default: \
        GRTInvalidOptionError(Ctrl, X, ""); \
        break; \

