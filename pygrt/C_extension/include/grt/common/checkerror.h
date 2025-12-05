/**
 * @file   checkerror.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 一些检查和报错的宏
 * 
 */

#pragma once 

#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <stdlib.h>

#include "grt/common/colorstr.h"

// GRT自定义报错信息
#define GRTRaiseError(ErrorMessage, ...) ({\
    fprintf(stderr, BOLD_RED ErrorMessage "\n" DEFAULT_RESTORE, ##__VA_ARGS__);\
    fflush(stderr);\
    exit(EXIT_FAILURE);\
})

// GRT自定义警告信息，不结束程序
#define GRTRaiseWarning(WarnMessage, ...) ({\
    fprintf(stderr, BOLD_YELLOW WarnMessage "\n" DEFAULT_RESTORE, ##__VA_ARGS__);\
    fflush(stderr);\
})

// GRT报错：选项设置不符要求
#define GRTBadOptionError(name, X, MoreErrorMessage, ...) ({\
    GRTRaiseError("[%s] Error in \"-"#X"\". "MoreErrorMessage" Use \"-h\" for help.\n", name, ##__VA_ARGS__);\
})

// GRT报错：选项未设置参数    注意这里使用的是 %c 和 运行时变量X
#define GRTMissArgsError(name, X, MoreErrorMessage, ...) ({\
    GRTRaiseError("[%s] Error! Option \"-%c\" requires an argument. "MoreErrorMessage" Use \"-h\" for help.\n", name, X, ##__VA_ARGS__);\
})

// GRT报错：非法选项    注意这里使用的是 %c 和 运行时变量X
#define GRTInvalidOptionError(name, X, MoreErrorMessage, ...) ({\
    GRTRaiseError("[%s] Error! Option \"-%c\" is invalid. "MoreErrorMessage" Use \"-h\" for help.\n", name, X, ##__VA_ARGS__);\
})

// GRT报错：文件不存在
#define GRTFileNotFoundError(name, filepath) ({\
    GRTRaiseError("[%s] Error! File \"%s\" not found. Please check.\n", name, filepath);\
})

// GRT报错：文件打开失败
#define GRTFileOpenError(name, filepath) ({\
    GRTRaiseError("[%s] Error! Cannot open File \"%s\". Please check.\n", name, filepath);\
})

// GRT报错：目录创建失败
#define GRTMakeDirError(name, dirpath, errno) ({\
    GRTRaiseError("[%s] Error! Unable to create folder %s. Error code: %d\n", name, dirpath, errno);\
})

// GRT报错：目录不存在
#define GRTDirNotFoundError(name, dirpath) ({\
    GRTRaiseError("[%s] Error! Directory \"%s\" not found. Please check.\n", name, dirpath);\
})


// ============================================================================================================================================

// GRT检查：某个选项是否启用
#define GRTCheckOptionActive(name, Ctrl, X) ({\
    if(!(Ctrl->X.active)){\
        GRTRaiseError("[%s] Error! Need set options \"-"#X"\". Use \"-h\" for help.\n", name);\
    }\
})

// GRT检查：是否有传递任何参数
#define GRTCheckOptionSet(name, condition) ({\
    if(!(condition)){\
        GRTRaiseError("[%s] Error! Need set options. Use \"-h\" for help.\n", name);\
    }\
})

// GRT检查：文件是否存在
#define GRTCheckFileExist(name, filepath) ({\
    if(access(filepath, F_OK) == -1){\
        GRTFileNotFoundError(name, filepath);\
    }\
})

// GRT检查：读取文件+检查文件指针+返回文件指针
#define GRTCheckOpenFile(name, filepath, mode) ({\
    FILE *_fp_ = fopen(filepath, mode);\
    if(_fp_ == NULL) {\
        GRTFileOpenError(name, filepath);\
    }\
    /** 返回文件指针 */ \
    _fp_;\
})

// GRT检查：创建目录
#define GRTCheckMakeDir(name, dirpath) ({\
    if(mkdir(dirpath, 0777) != 0){\
        if(errno != EEXIST){\
            GRTMakeDirError(name, dirpath, errno);\
        }\
    }\
})

// GRT检查：目录是否存在
#define GRTCheckDirExist(name, dirpath) ({\
    DIR *dir = opendir(dirpath);\
    if(dir == NULL) {\
        GRTDirNotFoundError(name, dirpath);\
    }\
    closedir(dir);\
})
