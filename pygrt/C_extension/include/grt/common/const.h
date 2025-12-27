/**
 * @file   const.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 *                   
 */

#pragma once

#include <complex.h> 
#include <tgmath.h>
#include <omp.h>

#include "grt/common/checkerror.h"

// CMPLX macro not exist on MacOS
#ifndef CMPLX
    #define CMPLX(real, imag) ((double)(real) + (double)(imag) * I)  ///< 复数扩展宏，添加此指令以适配MacOS
#endif


#if defined(_WIN32) || defined(__WIN32__)
    #define _TEST_WHETHER_WIN32_ 1
#else
    #define _TEST_WHETHER_WIN32_ 0  ///< 测试是否是windows系统
#endif

#if _TEST_WHETHER_WIN32_
    #include <direct.h> /* _mkdir */
    #define mkdir(x, y) _mkdir(x)  ///< 为windows系统修改mkdir函数
#endif


// 使用双精度
typedef double real_t;
typedef double complex cplx_t;
#define NC_REAL  NC_DOUBLE
#define NC_FUNC_REAL(func)  func##_double


// 常数
#define RTWOTHIRD 0.6666666666666667  ///< 2/3
#define PI  3.141592653589793        ///< \f$ \pi \f$
#define PI2 6.283185307179586        ///< \f$ 2\pi \f$
#define HALFPI 1.5707963267948966   ///< \f$ \frac{\pi}{2} \f$
#define QUARTERPI 0.7853981633974483   ///< \f$ \frac{\pi}{4} \f$
#define THREEQUARTERPI 2.356194490192345  ///< \f$ \frac{3\pi}{4} \f$
#define FIVEQUARTERPI  3.9269908169872414  ///< \f$ \frac{5\pi}{4} \f$
#define SEVENQUARTERPI  5.497787143782138   ///< \f$ \frac{7\pi}{4} \f$
#define INV_SQRT_TWO 0.7071067811865475   ///< \f$ \frac{1}{\sqrt{2}} \f$ 
#define DEG1 0.017453292519943295  ///< \f$ \frac{\pi}{180} \f$ 
#define GOLDEN_RATIO 0.6180339887498949  ///< \f$ \frac{\sqrt{5}-1}{2} \f$

#define GRT_MIN_DEPTH_GAP_SRC_RCV  1.0  ///< 震源和台站的最小深度差（不做绝对限制，仅用于参考波数积分上限，以及判断是否需要其它收敛方法）
#define GCC_ALWAYS_INLINE __attribute__((always_inline))  ///< gcc编译器不改动内联函数

#define GRT_SWAP(type, a, b) { type temp = a; a = b; b = temp; } ///< 交换两个变量的值
#define GRT_MIN_DISTANCE    1e-5   ///< 最小震中距，用于限制

#define GRT_STRING_FMT "%18s"  ///< 字符串输出格式
#define GRT_REAL_FMT "%18.8e"  ///< 浮点数输出格式
#define GRT_CMPLX_FMT "%18.8e%16.8e"   ///< 复数输出格式
#define GRT_STR_CMPLX_FMT "%34s"    ///< 与复数格式同长度的字符串输出格式
#define GRT_STRING_LONG_FMT "%25s"  ///< 字符串输出格式（更长）
#define GRT_REAL_LONG_FMT "%25.16e"  ///< 浮点数输出格式（更长）
#define GRT_CMPLX_LONG_FMT "%25.16e %25.16e"   ///< 复数输出格式（更长）
#define GRT_STR_CMPLX_LONG_FMT "%51s"    ///< 与复数格式同长度的字符串输出格式（更长）
#define GRT_CMPLX_SPLIT(x)  creal(x), cimag(x)   ///< 用于打印复数时将实部虚部分开

#define GRT_COMMENT_HEAD  '#'   ///< # 号， 作为注释的字符

#define GRT_NUM_STR "0123456789"   ///< 包括所有数字的字符串宏

#define GRT_MAX(a, b) ((a) > (b) ? (a) : (b))  ///< 求两者较大值
#define GRT_MIN(a, b) ((a) < (b) ? (a) : (b))  ///< 求两者较小值

#define GRT_SQUARE(x) ((x) * (x))  ///< 计算一个数的平方

// 内存管理
// 释放单个指针
#define GRT_SAFE_FREE_PTR(ptr) ({\
    if(ptr!=NULL) {\
        free(ptr);\
        ptr=NULL;\
    }\
})

// 释放指针数组
#define GRT_SAFE_FREE_PTR_ARRAY(ptr, count) ({\
    if(ptr!=NULL){\
        for(size_t i=0; i<count; ++i){\
            GRT_SAFE_FREE_PTR((ptr)[i]);\
            (ptr)[i] = NULL;\
        }\
        GRT_SAFE_FREE_PTR(ptr);\
        ptr=NULL;\
    }\
})

#define GRT_SAFE_ASPRINTF(ptr, fmt, ...) ({\
    int res;\
    if((res = asprintf(ptr, fmt, ##__VA_ARGS__)) == -1){\
        GRTRaiseError("Abnormal Error in function ‘asprintf’.\n");\
    };\
})

// -----------------------------------------------------------------------------
#define GRT_CHANNEL_NUM    3     ///< 3, 代码中分量个数（ZRT，ZNE），也是核函数类型个数(q, w, v)

#define GRT_INTEG_NUM   4    ///< 4, 代码中积分类型个数
#define GRT_MORDER_MAX   2    ///< 2, 代码中阶数m的最大值
#define GRT_SRC_M_NUM    6   ///< 6, 代码中不同震源、不同阶数的个数
#define GRT_MECHANISM_NUM   6   ///<  6, 描述震源机制的最多参数

#define GRT_PTAM_PT_MAX   36         ///< 36， 最后统计波峰波谷的目标数量
#define GRT_PTAM_WINDOW_SIZE  3      ///< 3，  使用连续点数判断是否为波峰或波谷
#define GRT_PTAM_WAITS_MAX    9      ///< 9,   判断波峰或波谷的最大等待次数，不能太小

#define GRT_INVERSE_SUCCESS   0      ///< 求逆或除法没有遇到除0错误
#define GRT_INVERSE_FAILURE   -1     ///< 求逆或除法遇到除0错误

#define GRT_GTYPES_MAX   2      ///< 2, 所有震源根据是否使用格林函数导数分为两类

typedef cplx_t     cplxChnlGrid[GRT_SRC_M_NUM][GRT_CHANNEL_NUM];
typedef cplx_t *pt_cplxChnlGrid[GRT_SRC_M_NUM][GRT_CHANNEL_NUM];
typedef real_t     realChnlGrid[GRT_SRC_M_NUM][GRT_CHANNEL_NUM];
typedef real_t *pt_realChnlGrid[GRT_SRC_M_NUM][GRT_CHANNEL_NUM];
typedef int        intChnlGrid [GRT_SRC_M_NUM][GRT_CHANNEL_NUM];

typedef cplx_t cplxIntegGrid[GRT_SRC_M_NUM][GRT_INTEG_NUM];
typedef real_t realIntegGrid[GRT_SRC_M_NUM][GRT_INTEG_NUM];
typedef size_t sizeIntegGrid[GRT_SRC_M_NUM][GRT_INTEG_NUM];

// 为以上静态多维数组的循环创建宏，谨慎使用，留意作用域
#define GRT_LOOP_ChnlGrid(im, c) \
    for(int im = 0; im < GRT_SRC_M_NUM; ++im) \
        for(int c = 0; c < GRT_CHANNEL_NUM; ++c)

#define GRT_LOOP_IntegGrid(im, v) \
    for(int im = 0; im < GRT_SRC_M_NUM; ++im) \
        for(int v = 0; v < GRT_INTEG_NUM; ++v)


/** 不同震源类型在大小为 GRT_SRC_M_NUM 的数组中的索引 */
enum {
    GRT_SRC_M_EX_INDEX = 0,
    GRT_SRC_M_VF_INDEX = 1,
    GRT_SRC_M_HF_INDEX = 2,
    GRT_SRC_M_DD_INDEX = 3,
    GRT_SRC_M_DS_INDEX = 4,
    GRT_SRC_M_SS_INDEX = 5,
};

/** 当前模块名，根据所调用模块进行切换 */
extern const char *GRT_MODULE_NAME;


/** 分别对应爆炸源(0阶)，垂直力源(0阶)，水平力源(1阶)，剪切源(0,1,2阶) */ 
extern const int GRT_SRC_M_ORDERS[GRT_SRC_M_NUM];

/** 不同震源类型使用的格林函数类型，0为Gij，1为格林函数导数Gij,k */
extern const int GRT_SRC_M_GTYPES[GRT_SRC_M_NUM];

/** 不同震源，不同阶数的名称简写，用于命名 */
extern const char *GRT_SRC_M_NAME_ABBR[GRT_SRC_M_NUM];

/** q, w, v 名称代号 */
extern const char GRT_QWV_CODES[];

/** ZRT三分量代号 */
extern const char GRT_ZRT_CODES[];

/** ZNE三分量代号 */
extern const char GRT_ZNE_CODES[];

/**
 * 设置OpenMP多线程数
 * 
 * @param[in]   num_threads        线程数
 */
void grt_set_num_threads(int num_threads);