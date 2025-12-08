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

/**
 * 计算核函数的函数指针，动态与静态的接口一致
 */
typedef void (*GRT_KernelFunc) (
    GRT_MODEL1D *mod1d, const real_t k, cplxChnlGrid QWV,
    bool calc_uiz, cplxChnlGrid QWV_uiz);


/**
 * kernel函数根据(5.5.3)式递推计算广义反射透射矩阵， 再根据公式得到
 * 
 *  1.EX 爆炸源， (P0)   
 *  2.VF  垂直力源, (P0, SV0)  
 *  3.HF  水平力源, (P1, SV1, SH1)  
 *  4.DC  剪切源, (P0, SV0), (P1, SV1, SH1), (P2, SV2, SH2)  
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
 *         故在程序可更方便的编写。（这个在静态情况下不成立，不能以此优化）
 *  @note  接收点位于自由表面的情况 不再单独考虑，合并在接受点浅于震源的情况
 *
 *
 *  为了尽量减少冗余的计算，且保证程序的可读性，可将震源层和接收层抽象为A,B层，
 *  即空间划分为FA,AB,BL, 计算这三个广义层的系数矩阵，再讨论震源层和接收层的深浅，
 *  计算相应的矩阵。  
 *
 *  @param[in,out]     mod1d           `MODEL1D` 结构体指针
 *  @param[in]     k               波数
 *  @param[out]    QWV             不同震源，不同阶数的核函数 \f$ q_m, w_m, v_m \f$
 *  @param[in]     calc_uiz        是否计算ui_z（位移u对坐标z的偏导）
 *  @param[out]    QWV_uiz         不同震源，不同阶数的核函数对z的偏导 \f$ \frac{\partial q_m}{\partial z}, \frac{\partial w_m}{\partial z}, \frac{\partial v_m}{\partial z} \f$
 * 
 */
void grt_kernel(
    GRT_MODEL1D *mod1d, const real_t k, cplxChnlGrid QWV,
    bool calc_uiz, cplxChnlGrid QWV_uiz);

/** 构建广义反射透射系数矩阵。作为 kernel 函数中的第一部分 */
void grt_GRT_matrix(GRT_MODEL1D *mod1d, const real_t k);

/** 从广义 R/T 矩阵出发，计算每个震源对应的核函数 QWV。 作为 kernel 函数中的第二部分 */
void grt_GRT_build_QWV(
    GRT_MODEL1D *mod1d, cplxChnlGrid QWV,
    bool calc_uiz, cplxChnlGrid QWV_uiz);

/** 静态解的核函数 */
void grt_static_kernel(
    GRT_MODEL1D *mod1d, const real_t k, cplxChnlGrid QWV,
    bool calc_uiz, cplxChnlGrid QWV_uiz);

/** 静态广义反射透射系数矩阵 */
void grt_static_GRT_matrix(GRT_MODEL1D *mod1d, const real_t k);

/** 静态 QWV */
void grt_static_GRT_build_QWV(
    GRT_MODEL1D *mod1d, cplxChnlGrid QWV,
    bool calc_uiz, cplxChnlGrid QWV_uiz);