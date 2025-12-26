/**
 * @file   myfftw.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 使用FFTW的辅助结构体和函数
 * 
 */

#include <string.h>
#include "grt/common/myfftw.h"
#include "grt/common/checkerror.h"

#define X(T, s, S) \
static GRT_FFTW##S##_HOLDER *  grt_init_fftw##s##_holder(const size_t nt, const real_t dt, const size_t nf_valid, const real_t df)\
{\
    GRT_FFTW##S##_HOLDER *fh = (GRT_FFTW##S##_HOLDER*)calloc(1, sizeof(GRT_FFTW##S##_HOLDER));\
    if (!fh) {\
        GRTRaiseError("Failed to allocate memory.\n");\
    }\
    fh->nt = nt;\
    fh->dt = dt;\
    fh->nf_valid = nf_valid;\
    fh->nf = nt/2+1;\
    fh->df = df;\
    fh->w_t = (T*)calloc(nt, sizeof(T));\
    fh->W_f = (fftw##s##_complex*)fftw##s##_malloc(sizeof(fftw##s##_complex) * fh->nf);\
    memset(fh->w_t, 0, sizeof(T)*nt);\
    memset(fh->W_f, 0, sizeof(fftw##s##_complex)*fh->nf);\
    if (!fh->w_t || !fh->W_f) {\
        GRTRaiseError("Failed to allocate arrays.\n");\
    }\
    return fh;\
}\
\
\
void grt_reset_fftw##s##_holder_zero(GRT_FFTW##S##_HOLDER *fh)\
{\
    memset(fh->w_t, 0, sizeof(T)*fh->nt);\
    memset(fh->W_f, 0, sizeof(fftw##s##_complex)*fh->nf);\
}\
\
\
GRT_FFTW##S##_HOLDER * grt_create_fftw##s##_holder_C2R_1D(const size_t nt, const real_t dt, const size_t nf_valid, const real_t df)\
{\
    GRT_FFTW##S##_HOLDER * fh = grt_init_fftw##s##_holder(nt, dt, nf_valid, df);\
    fh->plan = fftw##s##_plan_dft_c2r_1d(nt, fh->W_f, fh->w_t, FFTW_ESTIMATE);\
    return fh;\
}\
\
\
GRT_FFTW##S##_HOLDER * grt_create_fftw##s##_holder_R2C_1D(const size_t nt, const real_t dt, const size_t nf_valid, const real_t df)\
{\
    GRT_FFTW##S##_HOLDER * fh = grt_init_fftw##s##_holder(nt, dt, nf_valid, df);\
    fh->plan = fftw##s##_plan_dft_r2c_1d(nt, fh->w_t, fh->W_f, FFTW_ESTIMATE);\
    return fh;\
}\
\
\
void grt_destroy_fftw##s##_holder(GRT_FFTW##S##_HOLDER *fh)\
{\
    if (fh) {\
        fftw##s##_destroy_plan(fh->plan);\
        GRT_SAFE_FREE_PTR(fh->w_t);\
        GRT_SAFE_FFTW_FREE_PTR(fh->W_f, s);\
        GRT_SAFE_FREE_PTR(fh);\
    }\
}\
\
\
void grt_naive_inverse_transform_##T(GRT_FFTW##S##_HOLDER *fh)\
{\
    size_t nt = fh->nt;\
    size_t nf = fh->nf;\
    size_t nf_valid = fh->nf_valid;\
    T f;\
    for(size_t i = 0; i < nt; i++) {\
        T t = i*fh->dt;\
        fh->w_t[i] = 0.0;\
        /*单独处理零频*/\
        if(fh->f0 == 0.0)  fh->w_t[i] += creal(fh->W_f[0]); \
        /*对频率点进行求和（数值积分）*/\
        for(size_t k = 1; k < nf_valid-1; k++) {\
            f = k*fh->df + fh->f0;\
            fh->w_t[i] += 2.0*creal(fh->W_f[k] * exp(I * PI2 * f * t));\
        }\
        /*单独讨论最后一个点*/\
        f = fh->df*(nf_valid-1) + fh->f0;\
        if(nf == nf_valid){\
            if(nt%2==0){\
                fh->w_t[i] += creal(creal(fh->W_f[nf_valid-1]) * exp(I * PI2 * f * t));\
            } else {\
                fh->w_t[i] += 2.0*creal(fh->W_f[nf_valid-1] * exp(I * PI2 * f * t));\
            }\
        }\
    }\
}\

FOR_EACH_FFTW_TYPE
#undef X