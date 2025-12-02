/**
 * @file   kernel_template.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    广义 R/T 矩阵计算函数模板，以下代码通过递推公式实现 广义反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient 
 *               matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699
 * 
 * 
 */


// 注意这里不要使用 #pragma once


// 只有定义了特定宏才可以被 include
#if defined(__DYNAMIC_KERNEL__) || defined(__STATIC_KERNEL__)

// 不可以同时定义
#if defined(__DYNAMIC_KERNEL__) && defined(__STATIC_KERNEL__)
#error "Cannot define both __DYNAMIC_KERNEL__ and __STATIC_KERNEL__ simultaneously"
#endif

#include <stdio.h>
#include <complex.h>
#include <string.h>

#include "grt/common/RT_matrix.h"
#include "grt/common/model.h"
#include "grt/common/matrix.h"
#include "grt/common/prtdbg.h"



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
 * @param[in]     coef          震源系数，\f$ P_m, SV_m，SH_m\f$ ，维度2表示下行波(p=0)和上行波(p=1)
 * @param[out]    qwv           最终通过矩阵传播计算出的在台站位置的\f$ q_m,w_m,v_m\f$
 */
inline GCC_ALWAYS_INLINE void grt_construct_qwv(
    bool ircvup, 
    const cplx_t R1[2][2], cplx_t RL1, 
    const cplx_t R2[2][2], cplx_t RL2, 
    const cplx_t coef[GRT_QWV_NUM][2], cplx_t qwv[GRT_QWV_NUM])
{
    cplx_t qw0[2], qw1[2], v0;
    cplx_t coefD[2] = {coef[0][0], coef[1][0]};
    cplx_t coefU[2] = {coef[0][1], coef[1][1]};
    if(ircvup){
        grt_cmat2x1_mul(R2, coefD, qw0);
        qw0[0] += coefU[0]; qw0[1] += coefU[1]; 
        v0 = RL1 * (RL2*coef[2][0] + coef[2][1]);
    } else {
        grt_cmat2x1_mul(R2, coefU, qw0);
        qw0[0] += coefD[0]; qw0[1] += coefD[1]; 
        v0 = RL1 * (coef[2][0] + RL2*coef[2][1]);
    }
    grt_cmat2x1_mul(R1, qw0, qw1);

    qwv[0] = qw1[0];
    qwv[1] = qw1[1];
    qwv[2] = v0;
}



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
static void __KERNEL_FUNC__(
    GRT_MODEL1D *mod1d, const real_t k, cplx_t QWV[GRT_SRC_M_NUM][GRT_QWV_NUM],
    bool calc_uiz, cplx_t QWV_uiz[GRT_SRC_M_NUM][GRT_QWV_NUM])
{
    // 初始化qwv为0
    memset(QWV, 0, sizeof(cplx_t)*GRT_SRC_M_NUM*GRT_QWV_NUM);
    if(calc_uiz)  memset(QWV_uiz, 0, sizeof(cplx_t)*GRT_SRC_M_NUM*GRT_QWV_NUM);

    mod1d->k = k;
    // 为动态解计算xa,xb,caca,cbcb
    #if defined(__DYNAMIC_KERNEL__)
        grt_mod1d_xa_xb(mod1d, k);
    #endif

    bool ircvup = mod1d->ircvup;
    size_t isrc = mod1d->isrc; // 震源所在虚拟层位, isrc>=1
    size_t ircv = mod1d->ircv; // 接收点所在虚拟层位, ircv>=1, ircv != isrc
    size_t imin, imax; // 相对浅层深层层位
    imin = GRT_MIN(isrc, ircv);
    imax = GRT_MAX(isrc, ircv);

    // 广义层的反射透射系数
    RT_MATRIX *M_BL = &mod1d->M_BL;
    RT_MATRIX *M_AL = &mod1d->M_AL;
    RT_MATRIX *M_RS = &mod1d->M_RS;
    RT_MATRIX *M_FA = &mod1d->M_FA;
    RT_MATRIX *M_FB = &mod1d->M_FB;
    // 自由表面
    RT_MATRIX *M_top = &mod1d->M_top;

    // 定义物理层内的反射透射系数矩阵
    grt_init_RT_matrix(M);

    // 从顶到底进行矩阵递推, 公式(5.5.3)
    for(size_t iy=1; iy<mod1d->n; ++iy){ // 因为n>=3, 故一定会进入该循环

        // 只有动态解才可以使用这个 if 简化，
        // 对于静态解，即使是震源层、接收层这种虚拟层位也需要显式计算R/T矩阵
        #ifdef __DYNAMIC_KERNEL__
        if(iy != isrc && iy != ircv){
        #endif
            // 对第iy层的系数矩阵赋值
            __grt_RT_matrix_PSV(mod1d, iy, M);
            __grt_RT_matrix_SH(mod1d, iy, M);
            if(M->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

            // 加入时间延迟因子(第iy-1界面与第iy界面之间)
            __grt_delay_RT_matrix(mod1d, iy, M);
        #ifdef __DYNAMIC_KERNEL__
        }
        #endif

        // FA
        if(iy <= imin){ 
            if(iy == 1){ // 初始化FA
                memcpy(M_FA, M, sizeof(*M));
            }
            #ifdef __DYNAMIC_KERNEL__
            else if(iy == imin){
                grt_delay_GRT_matrix(mod1d, iy, M_FA);
            }
            #endif
            else { // 递推FA
                grt_recursion_RT_matrix(M_FA, M, M_FA);  // 写入原地址
                if(M_FA->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
            }
        } 
        // RS
        else if(iy <= imax){
            if(iy == imin+1){// 初始化RS
                memcpy(M_RS, M, sizeof(*M));
            }
            #ifdef __DYNAMIC_KERNEL__
            else if(iy == imax){
                grt_delay_GRT_matrix(mod1d, iy, M_RS);
            }
            #endif
            else { // 递推RS
                grt_recursion_RT_matrix(M_RS, M, M_RS);  // 写入原地址
                if(M_RS->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
            }
        } 
        // BL
        else {
            if(iy == imax+1){// 初始化BL
                memcpy(M_BL, M, sizeof(*M));
            } else { // 递推BL
                // 只有 RD 矩阵最终会被使用到
                grt_recursion_RT_matrix(M_BL, M, M_BL);  // 写入原地址
                if(M_BL->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
            }
        } // END if


    } // END for loop 
    //===================================================================================

    // return;


    // 计算震源系数
    __grt_source_coef(mod1d);

    // 临时中转矩阵 (temperary)
    cplx_t tmpR2[2][2], tmp2x2[2][2], tmpRL, tmp2x2_uiz[2][2], tmpRL2;

    // 递推RU_FA
    __grt_topfree_RU(mod1d);
    if(M_top->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
    grt_recursion_RU(M_top, M_FA, M_FA);
    if(M_FA->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
    
    // 根据震源和台站相对位置，计算最终的系数
    if(ircvup){ // A接收  B震源

        // 计算R_EV
        __grt_wave2qwv_REV_PSV(mod1d);
        __grt_wave2qwv_REV_SH(mod1d);

        // 递推RU_FS
        grt_recursion_RU(M_FA, M_RS, M_FB); // 已从ZR变为FR，加入了自由表面的效应
        if(M_FB->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
        
        // 公式(5.7.12-14)
        grt_cmat2x2_mul(M_BL->RD, M_FB->RU, tmpR2);
        grt_cmat2x2_one_sub(tmpR2);
        if((M_FB->stats = grt_cmat2x2_inv(tmpR2, tmpR2)) == GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;// (I - xx)^-1
        grt_cmat2x2_mul(M_FB->invT, tmpR2, tmp2x2);

        if(calc_uiz) grt_cmat2x2_assign(tmp2x2, tmp2x2_uiz); // 为后续计算空间导数备份

        grt_cmat2x2_mul(mod1d->R_EV, tmp2x2, tmp2x2);
        tmpRL = M_FB->invTL  / (1.0 - M_BL->RDL * M_FB->RUL);
        tmpRL2 = mod1d->R_EVL * tmpRL;

        for(int i=0; i<GRT_SRC_M_NUM; ++i){
            grt_construct_qwv(ircvup, tmp2x2, tmpRL2, M_BL->RD, M_BL->RDL, mod1d->src_coef[i], QWV[i]);
        }


        if(calc_uiz){
            __grt_wave2qwv_z_REV_PSV(mod1d);
            __grt_wave2qwv_z_REV_SH(mod1d);
            grt_cmat2x2_mul(mod1d->uiz_R_EV, tmp2x2_uiz, tmp2x2_uiz);
            tmpRL2 = mod1d->uiz_R_EVL * tmpRL;

            for(int i=0; i<GRT_SRC_M_NUM; ++i){
                grt_construct_qwv(ircvup, tmp2x2_uiz, tmpRL2, M_BL->RD, M_BL->RDL, mod1d->src_coef[i], QWV_uiz[i]);
            }    
        }
    } 
    else { // A震源  B接收

        // 计算R_EV
        __grt_wave2qwv_REV_PSV(mod1d);    
        __grt_wave2qwv_REV_SH(mod1d);    

        // 递推RD_SL
        grt_recursion_RD(M_RS, M_BL, M_AL);
        if(M_AL->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
        
        // 公式(5.7.26-27)
        grt_cmat2x2_mul(M_FA->RU, M_AL->RD, tmpR2);
        grt_cmat2x2_one_sub(tmpR2);
        if((M_AL->stats = grt_cmat2x2_inv(tmpR2, tmpR2)) == GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;// (I - xx)^-1
        grt_cmat2x2_mul(M_AL->invT, tmpR2, tmp2x2);

        if(calc_uiz) grt_cmat2x2_assign(tmp2x2, tmp2x2_uiz); // 为后续计算空间导数备份

        grt_cmat2x2_mul(mod1d->R_EV, tmp2x2, tmp2x2);
        tmpRL = M_AL->invTL / (1.0 - M_FA->RUL * M_AL->RDL);
        tmpRL2 = mod1d->R_EVL * tmpRL;

        for(int i=0; i<GRT_SRC_M_NUM; ++i){
            grt_construct_qwv(ircvup, tmp2x2, tmpRL2, M_FA->RU, M_FA->RUL, mod1d->src_coef[i], QWV[i]);
        }


        if(calc_uiz){
            __grt_wave2qwv_z_REV_PSV(mod1d);    
            __grt_wave2qwv_z_REV_SH(mod1d);    
            grt_cmat2x2_mul(mod1d->uiz_R_EV, tmp2x2_uiz, tmp2x2_uiz);
            tmpRL2 = mod1d->uiz_R_EVL * tmpRL;
            
            for(int i=0; i<GRT_SRC_M_NUM; ++i){
                grt_construct_qwv(ircvup, tmp2x2_uiz, tmpRL2, M_FA->RU, M_FA->RUL, mod1d->src_coef[i], QWV_uiz[i]);
            }
        }

    } // END if



    BEFORE_RETURN:

    // 对一些特殊情况的修正
    // 当震源和场点均位于地表时，可理论验证DS分量恒为0，这里直接赋0以避免后续的精度干扰
    if(mod1d->depsrc == 0.0 && mod1d->deprcv == 0.0)
    {
        for(int c=0; c<GRT_QWV_NUM; ++c){
            QWV[GRT_SRC_M_DS_INDEX][c] = 0.0;
            if(calc_uiz)  QWV_uiz[GRT_SRC_M_DS_INDEX][c] = 0.0;
        }
    }

}



#endif