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
 * @param     RD1[2][2]       (in)1层 P-SV 下传反射系数矩阵
 * @param     RDL1            (in)1层 SH 下传反射系数
 * @param     RU1[2][2]       (in)1层 P-SV 上传反射系数矩阵
 * @param     RUL1            (in)1层 SH 上传反射系数
 * @param     TD1[2][2]       (in)1层 P-SV 下传透射系数矩阵
 * @param     TDL1            (in)1层 SH 下传透射系数
 * @param     TU1[2][2]       (in)1层 P-SV 上传透射系数矩阵
 * @param     TUL1            (in)1层 SH 上传透射系数
 * @param     RD2[2][2]       (in)2层 P-SV 下传反射系数矩阵
 * @param     RDL2            (in)2层 SH 下传反射系数
 * @param     RD[2][2]        (out)1+2层 P-SV 下传反射系数矩阵
 * @param     RDL             (out)1+2层 SH 下传反射系数
 * @param     inv_2x2T[2][2]  (out) 非NULL时，返回公式中的 \f$ (\mathbf{I} - \mathbf{R}_U^1 \mathbf{R}_D^2)^{-1} \mathbf{T}_D^1 \f$ 一项   
 * @param     invT            (out) 非NULL时，返回上面inv_2x2T的标量形式      
 * 
 */
void recursion_RD(
    const MYCOMPLEX RD1[2][2], MYCOMPLEX RDL1, const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, 
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT);


/**
 * 根据公式(5.5.3(2))进行递推 
 * 
 * @param     RU1[2][2]       (in)1层 P-SV 上传反射系数矩阵
 * @param     RUL1            (in)1层 SH 上传反射系数
 * @param     TD1[2][2]       (in)1层 P-SV 下传透射系数矩阵
 * @param     TDL1            (in)1层 SH 下传透射系数
 * @param     RD2[2][2]       (in)2层 P-SV 下传反射系数矩阵
 * @param     RDL2            (in)2层 SH 下传反射系数
 * @param     TD2[2][2]       (in)2层 P-SV 下传透射系数矩阵
 * @param     TDL2            (in)2层 SH 下传透射系数
 * @param     TD[2][2]        (out)1+2层 P-SV 下传透射系数矩阵
 * @param     TDL             (out)1+2层 SH 下传透射系数
 * @param     inv_2x2T[2][2]  (out) 非NULL时，返回公式中的 \f$ (\mathbf{I} - \mathbf{R}_U^1 \mathbf{R}_D^2)^{-1} \mathbf{T}_D^1 \f$ 一项   
 * @param     invT            (out) 非NULL时，返回上面inv_2x2T的标量形式      
 * 
 */
void recursion_TD(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, 
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, 
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, 
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT);




/**
 * 根据公式(5.5.3(3))进行递推  
 * 
 * @param     RU1[2][2]       (in)1层 P-SV 上传反射系数矩阵
 * @param     RUL1            (in)1层 SH 上传反射系数
 * @param     RD2[2][2]       (in)2层 P-SV 下传反射系数矩阵
 * @param     RDL2            (in)2层 SH 下传反射系数
 * @param     RU2[2][2]       (in)2层 P-SV 上传反射系数矩阵
 * @param     RUL2            (in)2层 SH 上传反射系数
 * @param     TD2[2][2]       (in)2层 P-SV 下传透射系数矩阵
 * @param     TDL2            (in)2层 SH 下传透射系数
 * @param     TU2[2][2]       (in)2层 P-SV 上传透射系数矩阵
 * @param     TUL2            (in)2层 SH 上传透射系数
 * @param     RU[2][2]        (out)1+2层 P-SV 上传反射系数矩阵
 * @param     RUL             (out)1+2层 SH 上传反射系数
 * @param     inv_2x2T[2][2]  (out) 非NULL时，返回公式中的 \f$ (\mathbf{I} - \mathbf{R}_D^2 \mathbf{R}_U^1)^{-1} \mathbf{T}_U^2 \f$ 一项   
 * @param     invT            (out) 非NULL时，返回上面inv_2x2T的标量形式      
 * 
 */
void recursion_RU(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, const MYCOMPLEX RU2[2][2], MYCOMPLEX RUL2,
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX RU[2][2], MYCOMPLEX *RUL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT);

/**
 * 根据公式(5.5.3(4))进行递推
 * 
 * @param     RU1[2][2]       (in)1层 P-SV 上传反射系数矩阵
 * @param     RUL1            (in)1层 SH 上传反射系数
 * @param     RD2[2][2]       (in)2层 P-SV 下传反射系数矩阵
 * @param     RDL2            (in)2层 SH 下传反射系数
 * @param     RD2[2][2]       (in)2层 P-SV 下传反射系数矩阵
 * @param     RDL2            (in)2层 SH 下传反射系数
 * @param     TU2[2][2]       (in)2层 P-SV 上传透射系数矩阵
 * @param     TUL2            (in)2层 SH 上传透射系数
 * @param     TU[2][2]        (out)1+2层 P-SV 上传透射系数矩阵
 * @param     TUL             (out)1+2层 SH 上传透射系数
 * @param     inv_2x2T[2][2]  (out) 非NULL时，返回公式中的 \f$ (\mathbf{I} - \mathbf{R}_D^2 \mathbf{R}_U^1)^{-1} \mathbf{T}_U^2 \f$ 一项   
 * @param     invT            (out) 非NULL时，返回上面inv_2x2T的标量形式      
 * 
 * 
 */
void recursion_TU(
    const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2,
    const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX TU[2][2], MYCOMPLEX *TUL, MYCOMPLEX inv_2x2T[2][2], MYCOMPLEX *invT);



/**
 * 根据公式(5.5.3)进行递推，相当于上面四个函数合并
 * 
 * @param     RD1[2][2]       (in)1层 P-SV 下传反射系数矩阵
 * @param     RDL1            (in)1层 SH 下传反射系数
 * @param     RU1[2][2]       (in)1层 P-SV 上传反射系数矩阵
 * @param     RUL1            (in)1层 SH 上传反射系数
 * @param     TD1[2][2]       (in)1层 P-SV 下传透射系数矩阵
 * @param     TDL1            (in)1层 SH 下传透射系数
 * @param     TU1[2][2]       (in)1层 P-SV 上传透射系数矩阵
 * @param     TUL1            (in)1层 SH 上传透射系数
 * @param     RD2[2][2]       (in)2层 P-SV 下传反射系数矩阵
 * @param     RDL2            (in)2层 SH 下传反射系数
 * @param     RU2[2][2]       (in)2层 P-SV 上传反射系数矩阵
 * @param     RUL2            (in)2层 SH 上传反射系数
 * @param     TD2[2][2]       (in)2层 P-SV 下传透射系数矩阵
 * @param     TDL2            (in)2层 SH 下传透射系数
 * @param     TU2[2][2]       (in)2层 P-SV 上传透射系数矩阵
 * @param     TUL2            (in)2层 SH 上传透射系数
 * @param     RD[2][2]        (out)1+2层 P-SV 下传反射系数矩阵
 * @param     RDL             (out)1+2层 SH 下传反射系数
 * @param     RU[2][2]        (out)1+2层 P-SV 上传反射系数矩阵
 * @param     RUL             (out)1+2层 SH 上传反射系数
 * @param     TD[2][2]        (out)1+2层 P-SV 下传透射系数矩阵
 * @param     TDL             (out)1+2层 SH 下传透射系数
 * @param     TU[2][2]        (out)1+2层 P-SV 上传透射系数矩阵
 * @param     TUL             (out)1+2层 SH 上传透射系数
 * 
 */
void recursion_RT_2x2(
    const MYCOMPLEX RD1[2][2], MYCOMPLEX RDL1, const MYCOMPLEX RU1[2][2], MYCOMPLEX RUL1,
    const MYCOMPLEX TD1[2][2], MYCOMPLEX TDL1, const MYCOMPLEX TU1[2][2], MYCOMPLEX TUL1,
    const MYCOMPLEX RD2[2][2], MYCOMPLEX RDL2, const MYCOMPLEX RU2[2][2], MYCOMPLEX RUL2,
    const MYCOMPLEX TD2[2][2], MYCOMPLEX TDL2, const MYCOMPLEX TU2[2][2], MYCOMPLEX TUL2,
    MYCOMPLEX RD[2][2], MYCOMPLEX *RDL, MYCOMPLEX RU[2][2], MYCOMPLEX *RUL,
    MYCOMPLEX TD[2][2], MYCOMPLEX *TDL, MYCOMPLEX TU[2][2], MYCOMPLEX *TUL);


/**
 * 对于虚拟层位，即上下层是相同的物性参数，对公式(5.5.3)进行简化，只剩下时间延迟因子
 * 
 * @param     xa1      (in)P波归一化垂直波数 \f$ \sqrt{1 - (k_a/k)^2} \f$
 * @param     xb1      (in)S波归一化垂直波数 \f$ \sqrt{1 - (k_b/k)^2} \f$
 * @param     thk      (in)厚度
 * @param     k         (in)波数
 * @param     RU[2][2]       (inout)上层 P-SV 上传反射系数矩阵
 * @param     RUL            (inout)上层 SH 上传反射系数
 * @param     TD[2][2]       (inout)上层 P-SV 下传透射系数矩阵
 * @param     TDL            (inout)上层 SH 下传透射系数
 * @param     TU[2][2]       (inout)上层 P-SV 上传透射系数矩阵
 * @param     TUL            (inout)上层 SH 上传透射系数
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
 * \begin{pmatrix}
 * q_m \\
 * w_m 
 * \end{pmatrix}
 * =
 * \mathbf{R_1}
 * 
 * \left[
 * \mathbf{R_2}
 * \begin{pmatrix}
 * P_m^+ \\
 * SV_m^+ 
 * \end{pmatrix}
 * +
 * \begin{pmatrix}
 * P_m^- \\
 * SV_m^- 
 * \end{pmatrix}
 * 
 * \right]
 * 
 * \f]
 * 
 * + 当台站在震源下方时：
 * 
 * \f[
 * \begin{pmatrix}
 * q_m \\
 * w_m 
 * \end{pmatrix}
 * =
 * \mathbf{R_1}
 * 
 * \left[
 * \begin{pmatrix}
 * P_m^+ \\
 * SV_m^+ 
 * \end{pmatrix}
 * +
 * \mathbf{R_2}
 * \begin{pmatrix}
 * P_m^- \\
 * SV_m^- 
 * \end{pmatrix}
 * 
 * \right]
 * 
 * \f]
 * 
 * SH波类似，但是是标量形式。 
 * 
 * @param     ircvup    (in)接收层是否浅于震源层
 * @param     R1[2][2]  (in)P-SV波，\f$\mathbf{R_1}\f$矩阵
 * @param     RL1       (in)SH波，  \f$ R_1\f$
 * @param     R2[2][2]  (in)P-SV波，\f$\mathbf{R_2}\f$矩阵
 * @param     RL2       (in)SH波，  \f$ R_2\f$
 * @param     coef[3][2]  (in)震源系数，维度3表示震源附近的\f$ q_m,w_m,v_m\f$  ，维度2表示下行波(p=0)和上行波(p=1)
 * @param     qwv[3]      (out)最终通过矩阵传播计算出的在台站位置的\f$ q_m,w_m,v_m\f$
 */
void get_qwv(
    bool ircvup, 
    const MYCOMPLEX R1[2][2], MYCOMPLEX RL1, 
    const MYCOMPLEX R2[2][2], MYCOMPLEX RL2, 
    const MYCOMPLEX coef[3][2], MYCOMPLEX qwv[3]);
