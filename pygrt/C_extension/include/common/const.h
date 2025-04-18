/**
 * @file   const.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 *                   
 */

#pragma once

#include <complex.h> 
#include <math.h>
// #include <tgmath.h>

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
    // 数学函数
    #define EXP(x)    (expf(x))      ///< \f$ \exp(x) \f$
    #define CEXP(x)   (cexpf(x))      ///< 复数域 \f$ \exp(x) \f$
    #define CLOG(x)   (clogf(x))      ///< 复数域 \f$ \ln(x) \f$
    #define CSQRT(x)  (csqrtf(x))     ///< 复数域 \f$ x^{1/2} \f$
    #define SQRT(x)   (sqrtf(x))     ///< \f$ \sqrt{x} \f$
    #define SIN(x)    (sinf(x))      ///< \f$ \sin(x) \f$
    #define COS(x)    (cosf(x))      ///< \f$ \cos(x) \f$
    #define FABS(x)   (fabsf(x))     ///< \f$ |x| \f$
    #define CABS(x)   (cabsf(x))     ///< 复数域 \f$ |x| \f$
    #define J0(x)     (j0f(x))       ///< \f$ J_0(x) \f$
    #define J1(x)     (j1f(x))       ///< \f$ J_1(x) \f$
    #define JN(n,x)   (jnf((n),(x))) ///< \f$ J_n(x) \f$

    #define CREAL(x)  (crealf(x))    ///< 复数实部 \f$ \Re[x] \f$
    #define CIMAG(x)  (cimagf(x))    ///< 复数虚部 \f$ \Im[x] \f$

    // 常数
    #define RZERO 0.0f            ///< 0.0
    #define RQUART 0.25f          ///< 0.25
    #define RHALF 0.5f            ///< 0.5
    #define RONE 1.0f            ///< 1.0
    #define RTWO 2.0f            ///< 2.0
    #define RTHREE 3.0f            ///< 3.0
    #define RFOUR 4.0f            ///< 4.0
    #define PI  3.1415926f         ///< \f$ \pi \f$
    #define PI2 6.2831853f        ///< \f$ 2\pi \f$
    #define HALFPI 1.5707963f   ///< \f$ \frac{\pi}{2} \f$
    #define QUARTERPI 0.78539816f   ///< \f$ \frac{\pi}{4} \f$
    #define THREEQUARTERPI 2.35619450f  ///< \f$ \frac{3\pi}{4} \f$
    #define FIVEQUARTERPI  3.92699082f  ///< \f$ \frac{5\pi}{4} \f$
    #define SEVENQUARTERPI  5.49778714f   ///< \f$ \frac{7\pi}{4} \f$
    #define INV_SQRT_TWO 0.70710678f   ///< \f$ \frac{1}{\sqrt{2}} \f$ 
    #define DEG1 0.017453293f  ///< \f$ \frac{\pi}{180} \f$ 
    
#else 
    typedef double _Complex MYCOMPLEX;
    typedef double MYREAL;
    // 数学函数
    #define EXP(x)    (exp(x))      ///< \f$ \exp(x) \f$
    #define CEXP(x)   (cexp(x))      ///< 复数域 \f$ \exp(x) \f$
    #define CLOG(x)   (clog(x))      ///< 复数域 \f$ \ln(x) \f$
    #define CSQRT(x)  (csqrt(x))     ///< 复数域 \f$ x^{1/2} \f$
    #define SQRT(x)   (sqrt(x))     ///< \f$ \sqrt{x} \f$
    #define SIN(x)    (sin(x))      ///< \f$ \sin(x) \f$
    #define COS(x)    (cos(x))      ///< \f$ \cos(x) \f$
    #define FABS(x)   (fabs(x))     ///< \f$ |x| \f$
    #define CABS(x)   (cabs(x))     ///< 复数域 \f$ |x| \f$
    #define J0(x)     (j0(x))       ///< \f$ J_0(x) \f$
    #define J1(x)     (j1(x))       ///< \f$ J_1(x) \f$
    #define JN(n,x)   (jn((n),(x))) ///< \f$ J_n(x) \f$

    #define CREAL(x)  (creal(x))    ///< 复数实部 \f$ \Re[x] \f$
    #define CIMAG(x)  (cimag(x))    ///< 复数虚部 \f$ \Im[x] \f$

    // 常数
    #define RZERO 0.0            ///< 0.0
    #define RQUART 0.25          ///< 0.25
    #define RHALF 0.5             ///< 0.5
    #define RONE 1.0            ///< 1.0
    #define RTWO 2.0            ///< 2.0
    #define RTHREE 3.0            ///< 3.0
    #define RFOUR 4.0            ///< 4.0
    #define PI  3.141592653589793        ///< \f$ \pi \f$
    #define PI2 6.283185307179586        ///< \f$ 2\pi \f$
    #define HALFPI 1.5707963267948966   ///< \f$ \frac{\pi}{2} \f$
    #define QUARTERPI 0.7853981633974483   ///< \f$ \frac{\pi}{4} \f$
    #define THREEQUARTERPI 2.356194490192345  ///< \f$ \frac{3\pi}{4} \f$
    #define FIVEQUARTERPI  3.9269908169872414  ///< \f$ \frac{5\pi}{4} \f$
    #define SEVENQUARTERPI  5.497787143782138   ///< \f$ \frac{7\pi}{4} \f$
    #define INV_SQRT_TWO 0.7071067811865475   ///< \f$ \frac{1}{\sqrt{2}} \f$ 
    #define DEG1 0.017453292519943295  ///< \f$ \frac{\pi}{180} \f$ 
    
#endif

#define CZERO CMPLX(RZERO, RZERO)  ///< 0.0 + j0.0
#define CONE  CMPLX(RONE,  RZERO)  ///< 1.0 + j0.0
#define INIT_C_ZERO_2x2_MATRIX {{CZERO, CZERO}, {CZERO, CZERO}}   ///< 初始化复数0矩阵
#define INIT_C_IDENTITY_2x2_MATRIX {{CONE, CZERO}, {CZERO, CONE}}  ///< 初始化复数单位阵

#define MIN_DEPTH_GAP_SRC_RCV  1.0  ///< 震源和台站的最小深度差（不做绝对限制，仅用于参考波数积分上限）
#define GCC_ALWAYS_INLINE __attribute__((always_inline))  ///< gcc编译器不改动内联函数

#define GRT_STRING_FMT "%18s"  ///< 字符串输出格式
#define GRT_REAL_FMT "%18.8e"  ///< 浮点数输出格式
#define GRT_CMPLX_FMT "%18.8e%-+14.8eJ"   ///< 复数输出格式
#define GRT_STR_CMPLX_FMT "%34s"    ///< 与复数格式同长度的字符串输出格式

