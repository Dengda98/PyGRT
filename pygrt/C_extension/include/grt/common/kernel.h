/**
 * @file   kernel.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-06
 * 
 *    动态或静态下计算核函数的函数指针
 * 
 */


#pragma once 

#include "grt/common/model.h"
#include "grt/common/matrix.h"

/**
 * 计算核函数的函数指针，动态与静态的接口一致
 */
typedef void (*GRT_KernelFunc) (
    const GRT_MODEL1D *mod1d, cplx_t QWV[GRT_SRC_M_NUM][GRT_QWV_NUM],
    bool calc_uiz, cplx_t QWV_uiz[GRT_SRC_M_NUM][GRT_QWV_NUM]);


#if defined(__DYNAMIC_KERNEL__) || defined(__STATIC_KERNEL__)

    // ================ 动态解 =====================
    #ifdef __DYNAMIC_KERNEL__
        #define CONCAT_PREFIX(base)  grt_##base

        #include "grt/dynamic/layer.h"
        #include "grt/dynamic/source.h"
    #endif

    // ================ 静态解 =====================
    #ifdef __STATIC_KERNEL__
        #define CONCAT_PREFIX(base)  grt_static_##base

        #include "grt/static/static_layer.h"
        #include "grt/static/static_source.h"
    #endif


    #define __KERNEL_FUNC__           CONCAT_PREFIX(kernel)
    #define __grt_RT_matrix_PSV       CONCAT_PREFIX(RT_matrix_PSV)
    #define __grt_RT_matrix_SH        CONCAT_PREFIX(RT_matrix_SH)
    #define __grt_delay_RT_matrix     CONCAT_PREFIX(delay_RT_matrix)
    #define __grt_source_coef         CONCAT_PREFIX(source_coef)
    #define __grt_topfree_RU          CONCAT_PREFIX(topfree_RU)
    #define __grt_wave2qwv_REV_PSV    CONCAT_PREFIX(wave2qwv_REV_PSV)
    #define __grt_wave2qwv_REV_SH     CONCAT_PREFIX(wave2qwv_REV_SH)
    #define __grt_wave2qwv_z_REV_PSV  CONCAT_PREFIX(wave2qwv_z_REV_PSV)
    #define __grt_wave2qwv_z_REV_SH   CONCAT_PREFIX(wave2qwv_z_REV_SH)

    #include "grt/common/kernel_template.h"

#endif