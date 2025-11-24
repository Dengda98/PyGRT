/**
 * @file   myfftw.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 使用FFTW的辅助结构体和函数
 * 
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

// 二者保持以下顺序，使fftw_complex 变为C标准的双精度复数
#include <complex.h>
#include <fftw3.h>

#include "grt/common/const.h"


#define FOR_EACH_FFTW_TYPE   \
    X(double, , )       \
    X(float, f, F)


#define GRT_SAFE_FFTW_FREE_PTR(ptr, s) ({\
    if(ptr!=NULL) {\
        fftw##s##_free(ptr);\
    }\
})


#define X(T, s, S) \
typedef struct {\
    /*时间序列长度*/ \
    size_t nt; \
    /* 时间间隔 */ \
    real_t dt; \
    /*有效频谱长度*/\
    size_t nf_valid; \
    /*总频谱长度 nf = nt/2+1*/\
    size_t nf; \
    /*频率间隔*/\
    real_t df; \
    /*时间域实序列*/\
    T *w_t; \
    /*频率域复序列*/\
    fftw##s##_complex *W_f; \
    /*FFTW执行计划*/\
    fftw##s##_plan plan; \
    /*为一些特殊情况而准备的变量*/ \
    /*发生频移*/\
    real_t f0; \
    /*最终不执行FFTW而使用最朴素的傅里叶逆变换*/\
    bool naive_inv; \
} GRT_FFTW##S##_HOLDER;\
\
/**
 * 初始化 FFTW_HOLDER 结构体指针，进行复数频谱到实数序列的逆变换
 * 
 * @param[in]    nt           时间序列点数
 * @param[in]    dt           时间间隔
 * @param[in]    nf_valid     有效频谱点数
 * @param[in]    df           频率间隔
 * 
 * @return    FFTW_HOLDER 结构体指针
 * 
 */\
GRT_FFTW##S##_HOLDER * grt_create_fftw##s##_holder_C2R_1D(const size_t nt, const real_t dt, const size_t nf_valid, const real_t df); \
/** 初始化 FFTW_HOLDER 结构体指针，进行实数序列到复数频谱的正变换 */\
GRT_FFTW##S##_HOLDER * grt_create_fftw##s##_holder_R2C_1D(const size_t nt, const real_t dt, const size_t nf_valid, const real_t df); \
\
/** 将内部数据全部置零  */\
void grt_reset_fftw##s##_holder_zero(GRT_FFTW##S##_HOLDER *fh);\
\
/** 清理函数：释放结构体内存，防止内存泄漏  */\
void grt_destroy_fftw##s##_holder(GRT_FFTW##S##_HOLDER *fh);\
\
\
/**
 * 最朴素的非均匀频域采样逆变换
 * 严格按照傅里叶逆变换定义：g(t) = \int G(f) * e^(i2πft) df
 */ \
void grt_naive_inverse_transform_##T(GRT_FFTW##S##_HOLDER *fh);\

FOR_EACH_FFTW_TYPE
#undef X

