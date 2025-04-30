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


#include "common/const.h"


/**
 * 根据公式(5.5.3(1))进行递推  
 * 
 * @param[in]      RD1             1层 P-SV 下传反射系数矩阵
 * @param[in]      RDL1            1层 SH 下传反射系数
 * @param[in]      RU1             1层 P-SV 上传反射系数矩阵
 * @param[in]      RUL1            1层 SH 上传反射系数
 * @param[in]      TD1             1层 P-SV 下传透射系数矩阵
 * @param[in]      TDL1            1层 SH 下传透射系数
 * @param[in]      TU1             1层 P-SV 上传透射系数矩阵
 * @param[in]      TUL1            1层 SH 上传透射系数
 * @param[in]      RD2             2层 P-SV 下传反射系数矩阵
 * @param[in]      RDL2            2层 SH 下传反射系数
 * @param[out]     RD              1+2层 P-SV 下传反射系数矩阵
 * @param[out]     RDL             1+2层 SH 下传反射系数
 * @param[out]     inv_2x2T        非NULL时，返回公式中的 \f$ (\mathbf{I} - \mathbf{R}_U^1 \mathbf{R}_D^2)^{-1} \mathbf{T}_D^1 \f$ 一项   
 * @param[out]     invT            非NULL时，返回上面inv_2x2T的标量形式    
 * @param[out]     stats           状态代码，是否有除零错误，非0为异常值
 * 
 */
void recursion_RD(
    const MYCOMPLEX RD1[2][2], MYCOMPLEX RDL1, const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats);


/**
 * 根据公式(5.5.3(2))进行递推 
 * 
 * @param[in]     RU1                 1层 P-SV 上传反射系数矩阵
 * @param[in]     RUL1                1层 SH 上传反射系数
 * @param[in]     TD1                 1层 P-SV 下传透射系数矩阵
 * @param[in]     TDL1                1层 SH 下传透射系数
 * @param[in]     RD2                 2层 P-SV 下传反射系数矩阵
 * @param[in]     RDL2                2层 SH 下传反射系数
 * @param[in]     TD2                 2层 P-SV 下传透射系数矩阵
 * @param[in]     TDL2                2层 SH 下传透射系数
 * @param[out]    TD                  1+2层 P-SV 下传透射系数矩阵
 * @param[out]    TDL                 1+2层 SH 下传透射系数
 * @param[out]    inv_2x2T            非NULL时，返回公式中的 \f$ (\mathbf{I} - \mathbf{R}_U^1 \mathbf{R}_D^2)^{-1} \mathbf{T}_D^1 \f$ 一项   
 * @param[out]    invT                非NULL时，返回上面inv_2x2T的标量形式      
 * @param[out]    stats               状态代码，是否有除零错误，非0为异常值
 * 
 */
void recursion_TD(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, 
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, 
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats);




/**
 * 根据公式(5.5.3(3))进行递推  
 * 
 * @param[in]     RU1                 1层 P-SV 上传反射系数矩阵
 * @param[in]     RUL1                1层 SH 上传反射系数
 * @param[in]     RD2                 2层 P-SV 下传反射系数矩阵
 * @param[in]     RDL2                2层 SH 下传反射系数
 * @param[in]     RU2                 2层 P-SV 上传反射系数矩阵
 * @param[in]     RUL2                2层 SH 上传反射系数
 * @param[in]     TD2                 2层 P-SV 下传透射系数矩阵
 * @param[in]     TDL2                2层 SH 下传透射系数
 * @param[in]     TU2                 2层 P-SV 上传透射系数矩阵
 * @param[in]     TUL2                2层 SH 上传透射系数
 * @param[out]    RU                  1+2层 P-SV 上传反射系数矩阵
 * @param[out]    RUL                 1+2层 SH 上传反射系数
 * @param[out]    inv_2x2T            非NULL时，返回公式中的 \f$ (\mathbf{I} - \mathbf{R}_D^2 \mathbf{R}_U^1)^{-1} \mathbf{T}_U^2 \f$ 一项   
 * @param[out]    invT                非NULL时，返回上面inv_2x2T的标量形式      
 * @param[out]    stats               状态代码，是否有除零错误，非0为异常值
 * 
 */
void recursion_RU(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, const MYCOMPLEX RU2[2][2], MYCOMPLEX RUL2,
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats);

/**
 * 根据公式(5.5.3(4))进行递推
 * 
 * @param[in]     RU1                 1层 P-SV 上传反射系数矩阵
 * @param[in]     RUL1                1层 SH 上传反射系数
 * @param[in]     RD2                 2层 P-SV 下传反射系数矩阵
 * @param[in]     RDL2                2层 SH 下传反射系数
 * @param[in]     RD2                 2层 P-SV 下传反射系数矩阵
 * @param[in]     RDL2                2层 SH 下传反射系数
 * @param[in]     TU2                 2层 P-SV 上传透射系数矩阵
 * @param[in]     TUL2                2层 SH 上传透射系数
 * @param[out]    TU                  1+2层 P-SV 上传透射系数矩阵
 * @param[out]    TUL                 1+2层 SH 上传透射系数
 * @param[out]    inv_2x2T            非NULL时，返回公式中的 \f$ (\mathbf{I} - \mathbf{R}_D^2 \mathbf{R}_U^1)^{-1} \mathbf{T}_U^2 \f$ 一项   
 * @param[out]    invT                非NULL时，返回上面inv_2x2T的标量形式      
 * @param[out]    stats               状态代码，是否有除零错误，非0为异常值
 * 
 * 
 */
void recursion_TU(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2,
    const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT, MYINT *stats);



/**
 * 根据公式(5.5.3)进行递推，相当于上面四个函数合并
 * 
 * @param[in]     RD1                 1层 P-SV 下传反射系数矩阵
 * @param[in]     RDL1                1层 SH 下传反射系数
 * @param[in]     RU1                 1层 P-SV 上传反射系数矩阵
 * @param[in]     RUL1                1层 SH 上传反射系数
 * @param[in]     TD1                 1层 P-SV 下传透射系数矩阵
 * @param[in]     TDL1                1层 SH 下传透射系数
 * @param[in]     TU1                 1层 P-SV 上传透射系数矩阵
 * @param[in]     TUL1                1层 SH 上传透射系数
 * @param[in]     RD2                 2层 P-SV 下传反射系数矩阵
 * @param[in]     RDL2                2层 SH 下传反射系数
 * @param[in]     RU2                 2层 P-SV 上传反射系数矩阵
 * @param[in]     RUL2                2层 SH 上传反射系数
 * @param[in]     TD2                 2层 P-SV 下传透射系数矩阵
 * @param[in]     TDL2                2层 SH 下传透射系数
 * @param[in]     TU2                 2层 P-SV 上传透射系数矩阵
 * @param[in]     TUL2                2层 SH 上传透射系数
 * @param[out]    RD                  1+2层 P-SV 下传反射系数矩阵
 * @param[out]    RDL                 1+2层 SH 下传反射系数
 * @param[out]    RU                  1+2层 P-SV 上传反射系数矩阵
 * @param[out]    RUL                 1+2层 SH 上传反射系数
 * @param[out]    TD                  1+2层 P-SV 下传透射系数矩阵
 * @param[out]    TDL                 1+2层 SH 下传透射系数
 * @param[out]    TU                  1+2层 P-SV 上传透射系数矩阵
 * @param[out]    TUL                 1+2层 SH 上传透射系数
 * @param[out]    stats               状态代码，是否有除零错误，非0为异常值
 * 
 */
void recursion_RT_2x2(
    const MYCOMPLEX RD1[2][2], MYCOMPLEX RDL1, const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, const MYCOMPLEX RU2[2][2], MYCOMPLEX RUL2,
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL,
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYINT *stats);


/**
 * 对于虚拟层位，即上下层是相同的物性参数，对公式(5.5.3)进行简化，只剩下时间延迟因子
 * 
 * @param[in]         xa1            P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param[in]         xb1            S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param[in]         thk            厚度
 * @param[in]         k              波数
 * @param[in,out]     RU             上层 P-SV 上传反射系数矩阵
 * @param[in,out]     RUL            上层 SH 上传反射系数
 * @param[in,out]     TD             上层 P-SV 下传透射系数矩阵
 * @param[in,out]     TDL            上层 SH 下传透射系数
 * @param[in,out]     TU             上层 P-SV 上传透射系数矩阵
 * @param[in,out]     TUL            上层 SH 上传透射系数
 */
void recursion_RT_2x2_imaginary(
    MYCOMPLEX xa1, MYCOMPLEX xb1, MYREAL thk, MYREAL k, // 使用上层的厚度
    MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL);



/**
 * 最终公式(5.7.12,13,26,27)简化为 (P-SV波) :
 * + 当台站在震源上方时：
 * 
 * \f[ 
 * \begin{pmatrix} q_m \\ w_m  \end{pmatrix} = \mathbf{R_1} 
 * \left[ 
 * \mathbf{R_2} \begin{pmatrix}  P_m^+ \\ SV_m^+  \end{pmatrix}
 * + \begin{pmatrix}  P_m^- \\ SV_m^- \end{pmatrix}
 * \right]
 * \f]
 * 
 * + 当台站在震源下方时：
 * 
 * \f[
 * \begin{pmatrix} q_m \\ w_m  \end{pmatrix} = \mathbf{R_1}
 * \left[
 * \begin{pmatrix} P_m^+ \\ SV_m^+ \end{pmatrix}
 * + \mathbf{R_2} \begin{pmatrix} P_m^- \\ SV_m^- \end{pmatrix}
 * \right]
 * \f]
 * 
 * SH波类似，但是是标量形式。 
 * 
 * @param[in]     ircvup        接收层是否浅于震源层
 * @param[in]     R1            P-SV波，\f$\mathbf{R_1}\f$矩阵
 * @param[in]     RL1           SH波，  \f$ R_1\f$
 * @param[in]     R2            P-SV波，\f$\mathbf{R_2}\f$矩阵
 * @param[in]     RL2           SH波，  \f$ R_2\f$
 * @param[in]     coef          震源系数，\f$ P_m, SV_m, SH_m \f$ ，维度2表示下行波(p=0)和上行波(p=1)
 * @param[out]    qwv           最终通过矩阵传播计算出的在台站位置的\f$ q_m,w_m,v_m\f$
 */
void get_qwv(
    bool ircvup, 
    const MYCOMPLEX R1[2][2], MYCOMPLEX RL1, 
    const MYCOMPLEX R2[2][2], MYCOMPLEX RL2, 
    const MYCOMPLEX coef[QWV_NUM][2], MYCOMPLEX qwv[QWV_NUM]);
