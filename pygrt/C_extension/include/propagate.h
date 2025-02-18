/**
 * @file   propagate.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码通过递推公式实现 广义反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient 
 *               matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699
 * 
 *                   
 */

#pragma once 

#include "const.h"
#include "model.h"




/**
 * kernel函数根据(5.5.3)式递推计算广义反射透射矩阵， 再根据公式得到
 * 
 *      1.EXP 爆炸源， (P0)   
 *      2.VF  垂直力源, (P0, SV0)  
 *      3.HF  水平力源, (P1, SV1, SH1)  
 *      4.DC  剪切位错源, (P0, SV0), (P1, SV1, SH1), (P2, SV2, SH2)  
 *
 *  的 \f$ q_m, w_m, v_m \f$ 系数(\f$ m=0,1,2 \f$), 
 *
 *  eg. DC_qwv[i][j]表示 \f$ m=i \f$ 阶时的 \f$ q_m(j=0), w_m(j=1), v_m(j=2) \f$ 系数
 *
 * 在递推得到广义反射透射矩阵后，计算位移系数的公式本质是类似的，但根据震源和接受点的相对位置，
 * 空间划分为多个层，公式也使用不同的矩阵，具体为
 *
 *
 * \f[
 * \begin{array}{c}
 * \\\\  \hline
 * \hspace{5cm}\text{Free Surface(自由表面)}\hspace{5cm} \\\\ 
 * \text{...} \\\\  \hline
 * \text{Source/Receiver interface(震源/接收虚界面) (A面)} \\\\ 
 * \text{...} \\\\  \hline
 * \text{Receiver/Source interface(接收/震源虚界面) (B面)} \\\\ 
 * \text{...} \\\\  \hline
 * \text{Lower interface(底界面)} \\\\ 
 * \text{...} \\
 * \text{(无穷深)} \\
 * \text{...} \\ 
 * 
 * 
 * \end{array}
 * \f]
 *
 *  界面之间构成一个广义层，每个层都对应2个反射系数矩阵RD/RU和2个透射系数矩阵TD/TU,
 *  根据公式的整理结果，但实际需要的矩阵为：
 *  
 * |  广义层   | **台站在震源上方** | **台站在震源下方** |
 * |----------|-------------------|-------------------|
 * | FS (震源 <-> 表面) | RU             | RD, RU, TD, TU |
 * | FR (接收 <-> 表面) | RD, RU, TD, TU |       /        |
 * | RS (震源 <-> 接收) | RD, RU, TD, TU | RD, RU, TD, TU |
 * | SL (震源 <-> 底面) | RD             | RD             |
 * | RL (接收 <-> 底面) |       /        | RD             |
 * 
 * 
 *
 * 
 *  @note 关于与自由表面相关的系数矩阵要注意，FS表示(z1, zR+)之间的效应，但通常我们
 *        定义KP表示(zK+, zP+)之间的效应，所以这里F表示已经加入了自由表面的作用，
 *        对应的我们使用ZR表示(z1+, zR+)的效应，FR和ZR也满足类似的递推关系。
 *  @note  从公式推导上，例如RD_RS，描述的是(zR+, zS-)的效应，但由于我们假定
 *         震源位于介质层内，则z=zS并不是介质的物理分界面，此时 \f$ D_{j-1}^{-1} * D_j = I \f$，
 *         故在程序可更方便的编写。
 *  @note  接收点位于自由表面的情况 不再单独考虑，合并在接受点浅于震源的情况
 *
 *
 *  为了尽量减少冗余的计算，且保证程序的可读性，可将震源层和接收层抽象为A,B层，
 *  即空间划分为FA,AB,BL, 计算这三个广义层的系数矩阵，再讨论震源层和接收层的深浅，
 *  计算相应的矩阵。  
 *
 *  @param  mod1d     (in)`MODEL1D` 结构体指针
 *  @param  omega     (in)复数频率
 *  @param   k        (in)波数
 *  @param    EXP_qwv[3][3]    (out)爆炸源核函数
 *  @param    VF_qwv[3][3]     (out)垂直力源核函数
 *  @param    HF_qwv[3][3]     (out)水平力源核函数
 *  @param    DC_qwv[3][3]     (out)双力偶源核函数 
 * 
 */
void kernel(
    const MODEL1D *mod1d, MYCOMPLEX omega, MYREAL k,
    MYCOMPLEX EXP_qwv[3][3], MYCOMPLEX VF_qwv[3][3], MYCOMPLEX HF_qwv[3][3], MYCOMPLEX DC_qwv[3][3]);



/**
 * 计算核函数和Bessel函数的乘积，相当于计算了一个小积分区间内的值。参数中涉及两种数组形状：
 *    + [3][3]. 存储的是核函数，第一个维度3代表阶数m=0,1,2，第二个维度3代表三类系数qm,wm,vm  
 *    + [3][4]. 存储的是该dk区间内的积分值，维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型
 * 
 * 
 * @param     k     (in)波数
 * @param     r     (in)震中距 
 * @param    EXP_qwv[3][3]    (in)爆炸源核函数
 * @param    VF_qwv[3][3]     (in)垂直力源核函数
 * @param    HF_qwv[3][3]     (in)水平力源核函数
 * @param    DC_qwv[3][3]     (in)双力偶源核函数
 * @param    EXP_J[3][4]      (out)爆炸源，该dk区间内的积分值，下同
 * @param    VF_J[3][4]       (out)垂直力源
 * @param    HF_J[3][4]       (out)水平力源
 * @param    DC_J[3][4]       (out)双力偶源
 * 
 */
void int_Pk(
    MYREAL k, MYREAL r, 
    const MYCOMPLEX EXP_qwv[3][3], const MYCOMPLEX VF_qwv[3][3], 
    const MYCOMPLEX HF_qwv[3][3],  const MYCOMPLEX DC_qwv[3][3], 
    MYCOMPLEX EXP_J[3][4], MYCOMPLEX VF_J[3][4], 
    MYCOMPLEX HF_J[3][4],  MYCOMPLEX DC_J[3][4] );


/**
 * 将最终计算好的多个积分值，按照公式(5.6.22)组装成3分量。数组形状[3][4]，\
 * 存储的是最终的积分值，维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型
 * 
 * @param    sum_EXP_J[3][4]      (in)爆炸源，最终的积分值，下同
 * @param    sum_VF_J[3][4]       (in)垂直力源
 * @param    sum_HF_J[3][4]       (in)水平力源
 * @param    sum_DC_J[3][4]       (in)双力偶源
 * @param    tol_EXP[2]           (out)爆炸源的Z、R分量频谱结果
 * @param    tol_VF[2]            (out)垂直力源的Z、R分量频谱结果
 * @param    tol_HF[3]            (out)水平力源的Z、R、T分量频谱结果
 * @param    tol_DD[2]            (out)45度倾滑的Z、R分量频谱结果
 * @param    tol_DS[3]            (out)90度倾滑的Z、R、T分量频谱结果
 * @param    tol_SS[3]            (out)90度走滑的Z、R、T分量频谱结果
 */
void merge_Pk(
    const MYCOMPLEX sum_EXP_J[3][4], const MYCOMPLEX sum_VF_J[3][4], 
    const MYCOMPLEX sum_HF_J[3][4],  const MYCOMPLEX sum_DC_J[3][4], 
    MYCOMPLEX tol_EXP[2], MYCOMPLEX tol_VF[2], MYCOMPLEX tol_HF[3],
    MYCOMPLEX tol_DD[2],  MYCOMPLEX tol_DS[3], MYCOMPLEX tol_SS[3]);


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
