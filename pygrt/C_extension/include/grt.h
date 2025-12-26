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
#include <sys/stat.h>
#include <dirent.h>
#include <ctype.h>
#include <stdbool.h>

#include "grt/common/const.h"
#include "grt/common/logo.h"
#include "grt/common/checkerror.h"


#define GRT_MAIN_COMMAND   "grt"   ///< 主程序名

// ------------------------------------------------------
/** X 宏，用于定义子模块命令。后续的命令名称和函数名称均与此匹配 */
#define GRT_Module_List   \
    /* dynamic */          \
    X(greenfn)             \
    X(syn)                 \
    X(rotation)            \
    X(strain)              \
    X(stress)              \
    X(kernel)              \
    /* static */           \
    X(static_greenfn)      \
    X(static_syn)          \
    X(static_rotation)     \
    X(static_strain)       \
    X(static_stress)       \
    /* other */            \
    X(ker2asc)             \
    X(sac2asc)             \
    X(travt)               \
    X(lamb1)               \
// ------------------------------------------------------

/** 子模块的函数指针类型 */
typedef int (*GRT_MODULE_FUNC)(int argc, char **argv);

/** 子模块命令注册结构 */
typedef struct {
    const char *name;
    GRT_MODULE_FUNC func;
} GRT_MODULE_ENTRY;


/** 声明所有子模块函数 */
#define X(name) int name##_main(int argc, char **argv);
    GRT_Module_List
#undef X

/** 注册所有子模块命令 */
extern const GRT_MODULE_ENTRY GRT_Modules_Entry[];

/** 定义包含子模块名称的字符串数组 */
extern const char *GRT_Module_Names[];





/** 共有的命令行处理语句 */ 
#define GRT_Common_Options_in_Switch(X) \
    /** 帮助 */  \
    case 'h': \
        print_help(); \
        exit(EXIT_SUCCESS); \
        break; \
    /** 参数缺失 */  \
    case ':': \
        GRTMissArgsError(X, ""); \
        break; \
    /** 非法选项 */  \
    case '?': \
    default: \
        GRTInvalidOptionError(X, ""); \
        break; \

