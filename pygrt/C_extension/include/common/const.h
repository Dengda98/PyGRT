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



// #define GRT_USE_FLOAT  ///< 是否使用单精度浮点数

typedef int MYINT;  ///< 整数 


#ifdef GRT_USE_FLOAT 
    typedef float _Complex  MYCOMPLEX;   ///< 复数
    typedef float MYREAL;   ///< 浮点数
#else 
    typedef double _Complex MYCOMPLEX;
    typedef double MYREAL;
#endif

// 常数
#define RZERO 0.0            ///< 0.0
#define RQUART 0.25          ///< 0.25
#define RHALF 0.5             ///< 0.5
#define RONE 1.0            ///< 1.0
#define RTWO 2.0            ///< 2.0
#define RTHREE 3.0            ///< 3.0
#define RFOUR 4.0            ///< 4.0
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

#define CZERO CMPLX(RZERO, RZERO)  ///< 0.0 + j0.0
#define CONE  CMPLX(RONE,  RZERO)  ///< 1.0 + j0.0
#define INIT_C_ZERO_2x2_MATRIX {{CZERO, CZERO}, {CZERO, CZERO}}   ///< 初始化复数0矩阵
#define INIT_C_IDENTITY_2x2_MATRIX {{CONE, CZERO}, {CZERO, CONE}}  ///< 初始化复数单位阵

#define MIN_DEPTH_GAP_SRC_RCV  1.0  ///< 震源和台站的最小深度差（不做绝对限制，仅用于参考波数积分上限）
#define GCC_ALWAYS_INLINE __attribute__((always_inline))  ///< gcc编译器不改动内联函数

#define GRT_SWAP(type, a, b) { type temp = a; a = b; b = temp; } ///< 交换两个变量的值

#define GRT_STRING_FMT "%18s"  ///< 字符串输出格式
#define GRT_REAL_FMT "%18.8e"  ///< 浮点数输出格式
#define GRT_CMPLX_FMT "%18.8e%-+14.8ej"   ///< 复数输出格式
#define GRT_STR_CMPLX_FMT "%34s"    ///< 与复数格式同长度的字符串输出格式


// -----------------------------------------------------------------------------
#define CHANNEL_NUM    3     ///< 3, 代码中分量个数（ZRT，ZNE）

#define QWV_NUM     3   ///< 3, 代码中核函数类型个数(q, w, v)
#define INTEG_NUM   4    ///< 4, 代码中积分类型个数
#define MORDER_MAX   2    ///< 2, 代码中阶数m的最大值
#define SRC_M_NUM    6   ///< 6, 代码中不同震源、不同阶数的个数
#define MECHANISM_NUM   6   ///<  6, 描述震源机制的最多参数

#define PTAM_MAX_PT   36         ///< 36， 最后统计波峰波谷的目标数量
#define PTAM_WINDOW_SIZE  3      ///< 3，  使用连续点数判断是否为波峰或波谷
#define PTAM_MAX_WAITS    9      ///< 9,   判断波峰或波谷的最大等待次数，不能太小

#define INVERSE_SUCCESS   0      ///< 求逆或除法没有遇到除0错误
#define INVERSE_FAILURE   -1     ///< 求逆或除法遇到除0错误

#define GTYPES_MAX   2      ///< 2, 所有震源根据是否使用格林函数导数分为两类

/** 不同震源类型在全局数组中的索引 */
enum {
    SRC_M_EX_INDEX = 0,
    SRC_M_VF_INDEX = 1,
    SRC_M_HF_INDEX = 2,
    SRC_M_DD_INDEX = 3,
    SRC_M_DS_INDEX = 4,
    SRC_M_SS_INDEX = 5,
};

/** 分别对应爆炸源(0阶)，垂直力源(0阶)，水平力源(1阶)，剪切源(0,1,2阶) */ 
extern const MYINT SRC_M_ORDERS[SRC_M_NUM];

/** 不同震源类型使用的格林函数类型，0为Gij，1为格林函数导数Gij,k */
extern const MYINT SRC_M_GTYPES[SRC_M_NUM];

/** 不同震源，不同阶数的名称简写，用于命名 */
extern const char *SRC_M_NAME_ABBR[SRC_M_NUM];

/** q, w, v 名称代号 */
extern const char qwvchs[];

/** ZRT三分量代号 */
extern const char ZRTchs[];

/** ZNE三分量代号 */
extern const char ZNEchs[];