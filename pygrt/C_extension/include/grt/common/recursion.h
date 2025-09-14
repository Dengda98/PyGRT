/**
 * @file   recursion.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 * 以下代码通过递推公式计算两层的广义反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *
 */


#pragma once 


#include "grt/common/const.h"


// 对4个矩阵赋值
#define GRT_RT_PSV_ASSIGN(suffix) ({\
    for(MYINT kk=0; kk<2; ++kk){\
        for(MYINT pp=0; pp<2; ++pp){\
            RD_##suffix[kk][pp] = RD[kk][pp];\
            RU_##suffix[kk][pp] = RU[kk][pp];\
            TD_##suffix[kk][pp] = TD[kk][pp];\
            TU_##suffix[kk][pp] = TU[kk][pp];\
        }\
    }\
})

#define GRT_RT_SH_ASSIGN(suffix) ({\
    RDL_##suffix = RDL;\
    RUL_##suffix = RUL;\
    TDL_##suffix = TDL;\
    TUL_##suffix = TUL;\
})


/** 合并 recursion_RD_PSV(SH) */
void grt_recursion_RD(
    const MYCOMPLEX RD1[2][2], MYCOMPLEX RDL1, const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats);

/**
 * 根据公式(5.5.3(1))进行递推 P-SV
 * 
 * @param[in]      RD1             1层 P-SV 下传反射系数矩阵
 * @param[in]      RU1             1层 P-SV 上传反射系数矩阵
 * @param[in]      TD1             1层 P-SV 下传透射系数矩阵
 * @param[in]      TU1             1层 P-SV 上传透射系数矩阵
 * @param[in]      RD2             2层 P-SV 下传反射系数矩阵
 * @param[out]     RD              1+2层 P-SV 下传反射系数矩阵
 * @param[out]     inv_2x2T        非NULL时，返回公式中的 \f$ (\mathbf{I} - \mathbf{R}_U^1 \mathbf{R}_D^2)^{-1} \mathbf{T}_D^1 \f$ 一项   
 * @param[out]     stats           状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_recursion_RD_PSV(
    const MYCOMPLEX RD1[2][2], const MYCOMPLEX RU1[2][2],
    const MYCOMPLEX TD1[2][2], const MYCOMPLEX TU1[2][2],
    const MYCOMPLEX RD2[2][2],
    MYCOMPLEX RD[2][2], MYCOMPLEX inv_2x2T[2][2], MYINT *stats);

/**
 * 根据公式(5.5.3(1))进行递推 SH
 * 
 * @param[in]      RDL1            1层 SH 下传反射系数
 * @param[in]      RUL1            1层 SH 上传反射系数
 * @param[in]      TDL1            1层 SH 下传透射系数
 * @param[in]      TUL1            1层 SH 上传透射系数
 * @param[in]      RDL2            2层 SH 下传反射系数
 * @param[out]     RDL             1+2层 SH 下传反射系数
 * @param[out]     invT            非NULL时，返回上面inv_2x2T的标量形式    
 * @param[out]     stats           状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_recursion_RD_SH(
    MYCOMPLEX RDL1, MYCOMPLEX RUL1,
    MYCOMPLEX TDL1, MYCOMPLEX TUL1,
    MYCOMPLEX RDL2, MYCOMPLEX *RDL, MYCOMPLEX *invT, MYINT *stats);


/** 合并 recursion_TD_PSV(SH) */
void grt_recursion_TD(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, 
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, 
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats);

/**
 * 根据公式(5.5.3(2))进行递推 P-SV
 * 
 * @param[in]     RU1                 1层 P-SV 上传反射系数矩阵
 * @param[in]     TD1                 1层 P-SV 下传透射系数矩阵
 * @param[in]     RD2                 2层 P-SV 下传反射系数矩阵
 * @param[in]     TD2                 2层 P-SV 下传透射系数矩阵
 * @param[out]    TD                  1+2层 P-SV 下传透射系数矩阵
 * @param[out]    inv_2x2T            非NULL时，返回公式中的 \f$ (\mathbf{I} - \mathbf{R}_U^1 \mathbf{R}_D^2)^{-1} \mathbf{T}_D^1 \f$ 一项   
 * @param[out]    stats               状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_recursion_TD_PSV(
    const MYCOMPLEX RU1[2][2], const MYCOMPLEX TD1[2][2],
    const MYCOMPLEX RD2[2][2], const MYCOMPLEX TD2[2][2],
    MYCOMPLEX TD[2][2], MYCOMPLEX inv_2x2T[2][2], MYINT *stats);

/**
 * 根据公式(5.5.3(2))进行递推 SH
 * 
 * @param[in]     RUL1                1层 SH 上传反射系数
 * @param[in]     TDL1                1层 SH 下传透射系数
 * @param[in]     RDL2                2层 SH 下传反射系数
 * @param[in]     TDL2                2层 SH 下传透射系数
 * @param[out]    TDL                 1+2层 SH 下传透射系数
 * @param[out]    invT                非NULL时，返回上面inv_2x2T的标量形式      
 * @param[out]    stats               状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_recursion_TD_SH(
    MYCOMPLEX RUL1, MYCOMPLEX TDL1, 
    MYCOMPLEX RDL2, MYCOMPLEX TDL2, 
    MYCOMPLEX *TDL, MYCOMPLEX *invT, MYINT *stats);


/** 合并 recursion_RU_PSV(SH) */
void grt_recursion_RU(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, const MYCOMPLEX RU2[2][2], MYCOMPLEX RUL2,
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats);

/**
 * 根据公式(5.5.3(3))进行递推 P-SV
 * 
 * @param[in]     RU1                 1层 P-SV 上传反射系数矩阵
 * @param[in]     RD2                 2层 P-SV 下传反射系数矩阵
 * @param[in]     RU2                 2层 P-SV 上传反射系数矩阵
 * @param[in]     TD2                 2层 P-SV 下传透射系数矩阵
 * @param[in]     TU2                 2层 P-SV 上传透射系数矩阵
 * @param[out]    RU                  1+2层 P-SV 上传反射系数矩阵
 * @param[out]    inv_2x2T            非NULL时，返回公式中的 \f$ (\mathbf{I} - \mathbf{R}_D^2 \mathbf{R}_U^1)^{-1} \mathbf{T}_U^2 \f$ 一项   
 * @param[out]    stats               状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_recursion_RU_PSV(
    const MYCOMPLEX RU1[2][2], const MYCOMPLEX RD2[2][2], 
    const MYCOMPLEX RU2[2][2], const MYCOMPLEX TD2[2][2], 
    const MYCOMPLEX TU2[2][2], 
    MYCOMPLEX RU[2][2], MYCOMPLEX inv_2x2T[2][2], MYINT *stats);

/**
 * 根据公式(5.5.3(3))进行递推 SH
 * 
 * @param[in]     RUL1                1层 SH 上传反射系数
 * @param[in]     RDL2                2层 SH 下传反射系数
 * @param[in]     RUL2                2层 SH 上传反射系数
 * @param[in]     TDL2                2层 SH 下传透射系数
 * @param[in]     TUL2                2层 SH 上传透射系数
 * @param[out]    RUL                 1+2层 SH 上传反射系数
 * @param[out]    invT                非NULL时，返回上面inv_2x2T的标量形式      
 * @param[out]    stats               状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_recursion_RU_SH(
    MYCOMPLEX RUL1, MYCOMPLEX RDL2,
    MYCOMPLEX RUL2, MYCOMPLEX TDL2, 
    MYCOMPLEX TUL2,
    MYCOMPLEX *RUL, MYCOMPLEX *invT, MYINT *stats);


/** 合并 recursion_TU_PSV(SH) */
void grt_recursion_TU(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2,
    const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats);

/**
 * 根据公式(5.5.3(4))进行递推 P-SV
 * 
 * @param[in]     RU1                 1层 P-SV 上传反射系数矩阵
 * @param[in]     RD2                 2层 P-SV 下传反射系数矩阵
 * @param[in]     RD2                 2层 P-SV 下传反射系数矩阵
 * @param[in]     TU2                 2层 P-SV 上传透射系数矩阵
 * @param[out]    TU                  1+2层 P-SV 上传透射系数矩阵
 * @param[out]    inv_2x2T            非NULL时，返回公式中的 \f$ (\mathbf{I} - \mathbf{R}_D^2 \mathbf{R}_U^1)^{-1} \mathbf{T}_U^2 \f$ 一项   
 * @param[out]    stats               状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_recursion_TU_PSV(
    const MYCOMPLEX RU1[2][2], const MYCOMPLEX TU1[2][2],
    const MYCOMPLEX RD2[2][2], const MYCOMPLEX TU2[2][2],
    MYCOMPLEX TU[2][2], MYCOMPLEX inv_2x2T[2][2], MYINT *stats);

/**
 * 根据公式(5.5.3(4))进行递推 SH
 * 
 * @param[in]     RUL1                1层 SH 上传反射系数
 * @param[in]     RDL2                2层 SH 下传反射系数
 * @param[in]     RDL2                2层 SH 下传反射系数
 * @param[in]     TUL2                2层 SH 上传透射系数
 * @param[out]    TUL                 1+2层 SH 上传透射系数
 * @param[out]    invT                非NULL时，返回上面inv_2x2T的标量形式      
 * @param[out]    stats               状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_recursion_TU_SH(
    MYCOMPLEX RUL1, MYCOMPLEX TUL1,
    MYCOMPLEX RDL2, MYCOMPLEX TUL2,
    MYCOMPLEX *TUL, MYCOMPLEX *invT, MYINT *stats);




/** 合并 recursion_RT_matrix_PSV(SH) */
void grt_recursion_RT_matrix(
    const MYCOMPLEX RD1[2][2], MYCOMPLEX RDL1, const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, const MYCOMPLEX RU2[2][2], MYCOMPLEX RUL2,
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL,
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats);

/**
 * 根据公式(5.5.3)进行递推 P-SV ，相当于对应四个函数合并，
 * 内部使用了共有变量防止重复计算
 * 
 * @param[in]     RD1                 1层 P-SV 下传反射系数矩阵
 * @param[in]     RU1                 1层 P-SV 上传反射系数矩阵
 * @param[in]     TD1                 1层 P-SV 下传透射系数矩阵
 * @param[in]     TU1                 1层 P-SV 上传透射系数矩阵
 * @param[in]     RD2                 2层 P-SV 下传反射系数矩阵
 * @param[in]     RU2                 2层 P-SV 上传反射系数矩阵
 * @param[in]     TD2                 2层 P-SV 下传透射系数矩阵
 * @param[in]     TU2                 2层 P-SV 上传透射系数矩阵
 * @param[out]    RD                  1+2层 P-SV 下传反射系数矩阵
 * @param[out]    RU                  1+2层 P-SV 上传反射系数矩阵
 * @param[out]    TD                  1+2层 P-SV 下传透射系数矩阵
 * @param[out]    TU                  1+2层 P-SV 上传透射系数矩阵
 * @param[out]    stats               状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_recursion_RT_matrix_PSV(
    const MYCOMPLEX RD1[2][2], const MYCOMPLEX RU1[2][2],
    const MYCOMPLEX TD1[2][2], const MYCOMPLEX TU1[2][2],
    const MYCOMPLEX RD2[2][2], const MYCOMPLEX RU2[2][2],
    const MYCOMPLEX TD2[2][2], const MYCOMPLEX TU2[2][2],
    MYCOMPLEX RD[2][2], MYCOMPLEX RU[2][2],
    MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2], MYINT *stats);

/**
 * 根据公式(5.5.3)进行递推 SH ，相当于对应四个函数合并，
 * 内部使用了共有变量防止重复计算
 * 
 * @param[in]     RDL1                1层 SH 下传反射系数
 * @param[in]     RUL1                1层 SH 上传反射系数
 * @param[in]     TDL1                1层 SH 下传透射系数
 * @param[in]     TUL1                1层 SH 上传透射系数
 * @param[in]     RDL2                2层 SH 下传反射系数
 * @param[in]     RUL2                2层 SH 上传反射系数
 * @param[in]     TDL2                2层 SH 下传透射系数
 * @param[in]     TUL2                2层 SH 上传透射系数
 * @param[out]    RDL                 1+2层 SH 下传反射系数
 * @param[out]    RUL                 1+2层 SH 上传反射系数
 * @param[out]    TDL                 1+2层 SH 下传透射系数
 * @param[out]    TUL                 1+2层 SH 上传透射系数
 * @param[out]    stats               状态代码，是否有除零错误，非0为异常值
 * 
 */
void grt_recursion_RT_matrix_SH(
    MYCOMPLEX RDL1, MYCOMPLEX RUL1,
    MYCOMPLEX TDL1, MYCOMPLEX TUL1,
    MYCOMPLEX RDL2, MYCOMPLEX RUL2,
    MYCOMPLEX TDL2, MYCOMPLEX TUL2,
    MYCOMPLEX *RDL, MYCOMPLEX *RUL,
    MYCOMPLEX *TDL, MYCOMPLEX *TUL, MYINT *stats);



/** 合并 recursion_RT_matrix_PSV(SH)_virtual */
void grt_recursion_RT_matrix_virtual(
    MYCOMPLEX xa1, MYCOMPLEX xb1, MYREAL thk, MYREAL k, // 使用上层的厚度
    MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL);

/**
 * 递推虚拟层位的 P-SV 矩阵，即上下层是相同的物性参数，对公式(5.5.3)进行简化，只剩下时间延迟因子
 * 
 * @param[in]         xa1            P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]         xb1            S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]         thk            厚度
 * @param[in]         k              波数
 * @param[in,out]     RU             上层 P-SV 上传反射系数矩阵
 * @param[in,out]     TD             上层 P-SV 下传透射系数矩阵
 * @param[in,out]     TU             上层 P-SV 上传透射系数矩阵
 */
void grt_recursion_RT_matrix_PSV_virtual(
    MYCOMPLEX xa1, MYCOMPLEX xb1, MYREAL thk, MYREAL k, // 使用上层的厚度
    MYCOMPLEX RU[2][2], MYCOMPLEX TD[2][2], MYCOMPLEX TU[2][2]);

/**
 * 递推虚拟层位的 SH 矩阵，即上下层是相同的物性参数，对公式(5.5.3)进行简化，只剩下时间延迟因子
 * 
 * @param[in]         xb1            S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]         thk            厚度
 * @param[in]         k              波数
 * @param[in,out]     RUL            上层 SH 上传反射系数
 * @param[in,out]     TDL            上层 SH 下传透射系数
 * @param[in,out]     TUL            上层 SH 上传透射系数
 */
void grt_recursion_RT_matrix_SH_virtual(
    MYCOMPLEX xb1, MYREAL thk, MYREAL k, // 使用上层的厚度
    MYCOMPLEX *RUL, MYCOMPLEX *TDL, MYCOMPLEX *TUL);



