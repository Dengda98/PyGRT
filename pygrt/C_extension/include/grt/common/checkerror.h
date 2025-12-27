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
    fprintf(stderr, BOLD_WHITE "%s:%d: In function ‘%s’: \n" BOLD_RED "[%s] Error! " ErrorMessage "\n" DEFAULT_RESTORE, __FILE__, __LINE__, __func__, GRT_MODULE_NAME, ##__VA_ARGS__);\
    fflush(stderr);\
    exit(EXIT_FAILURE);\
})

// GRT自定义一般信息
#define GRTRaiseInfo(ErrorMessage, ...) ({\
    fprintf(stdout, REGULAR_GREEN "[%s] " ErrorMessage "\n" DEFAULT_RESTORE, GRT_MODULE_NAME, ##__VA_ARGS__);\
    fflush(stdout);\
})

// GRT自定义警告信息，不结束程序
#define GRTRaiseWarning(WarnMessage, ...) ({\
    fprintf(stdout, BOLD_WHITE "%s:%d: In function ‘%s’: \n" BOLD_YELLOW "[%s] Warning! " WarnMessage "\n" DEFAULT_RESTORE, __FILE__, __LINE__, __func__, GRT_MODULE_NAME, ##__VA_ARGS__);\
    fflush(stdout);\
})

// GRT报错：选项设置不符要求
#define GRTBadOptionError(X, MoreErrorMessage, ...) ({\
    GRTRaiseError("Error in \"-"#X"\". "MoreErrorMessage" Use \"-h\" for help.\n", ##__VA_ARGS__);\
})

// GRT报错：选项未设置参数    注意这里使用的是 %c 和 运行时变量X
#define GRTMissArgsError(X, MoreErrorMessage, ...) ({\
    GRTRaiseError("Option \"-%c\" requires an argument. "MoreErrorMessage" Use \"-h\" for help.\n", X, ##__VA_ARGS__);\
})

// GRT报错：非法选项    注意这里使用的是 %c 和 运行时变量X
#define GRTInvalidOptionError(X, MoreErrorMessage, ...) ({\
    GRTRaiseError("Option \"-%c\" is invalid. "MoreErrorMessage" Use \"-h\" for help.\n", X, ##__VA_ARGS__);\
})

// GRT报错：文件不存在
#define GRTFileNotFoundError(filepath) ({\
    GRTRaiseError("File \"%s\" not found. Please check.\n", filepath);\
})

// GRT报错：文件打开失败
#define GRTFileOpenError(filepath) ({\
    GRTRaiseError("Cannot open File \"%s\". Please check.\n", filepath);\
})

// GRT报错：目录创建失败
#define GRTMakeDirError(dirpath, errno) ({\
    GRTRaiseError("Unable to create folder %s. Error code: %d\n", dirpath, errno);\
})

// GRT报错：目录不存在
#define GRTDirNotFoundError(dirpath) ({\
    GRTRaiseError("Directory \"%s\" not found. Please check.\n", dirpath);\
})


// ============================================================================================================================================

// GRT检查：某个选项是否启用
#define GRTCheckOptionActive(Ctrl, X) ({\
    if(!(Ctrl->X.active)){\
        GRTRaiseError("Need set options \"-"#X"\". Use \"-h\" for help.\n");\
    }\
})

// GRT检查：是否有传递任何参数
#define GRTCheckOptionSet(condition) ({\
    if(!(condition)){\
        GRTRaiseError("Need set options. Use \"-h\" for help.\n");\
    }\
})

// GRT检查：文件是否存在
#define GRTCheckFileExist(filepath) ({\
    if(access(filepath, F_OK) == -1){\
        GRTFileNotFoundError(filepath);\
    }\
})

// GRT检查：读取文件+检查文件指针+返回文件指针
#define GRTCheckOpenFile(filepath, mode) ({\
    FILE *_fp_ = fopen(filepath, mode);\
    if(_fp_ == NULL) {\
        GRTFileOpenError(filepath);\
    }\
    /** 返回文件指针 */ \
    _fp_;\
})

// GRT检查：创建目录
#define GRTCheckMakeDir(dirpath) ({\
    if(mkdir(dirpath, 0777) != 0){\
        if(errno != EEXIST){\
            GRTMakeDirError(dirpath, errno);\
        }\
    }\
})

// GRT检查：目录是否存在
#define GRTCheckDirExist(dirpath) ({\
    DIR *dir = opendir(dirpath);\
    if(dir == NULL) {\
        GRTDirNotFoundError(dirpath);\
    }\
    closedir(dir);\
})
