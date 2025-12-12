/**
 * @file   kernel.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-12
 * 
 *    动态或静态下计算核函数
 * 
 */

#include <stdbool.h>
#include "grt/common/matrix.h"
#include "grt/integral/kernel.h"

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
 * @param[in]     coefD         下行震源系数，\f$ P_m, SV_m，SH_m\f$
 * @param[in]     coefU         上行震源系数，\f$ P_m, SV_m，SH_m\f$
 * @param[out]    QWV           最终通过矩阵传播计算出的在台站位置的\f$ q_m,w_m,v_m\f$
 */
inline GCC_ALWAYS_INLINE void grt_construct_qwv(
    bool ircvup, 
    const cplx_t R1[2][2], cplx_t RL1, 
    const cplx_t R2[2][2], cplx_t RL2, 
    const cplxChnlGrid coefD, const cplxChnlGrid coefU, cplxChnlGrid QWV)
{
    for(int i = 0; i < GRT_SRC_M_NUM; ++i)
    {
        if(ircvup){
            grt_cmat2x1_mul(R2, coefD[i], QWV[i]);
            QWV[i][0] += coefU[i][0];
            QWV[i][1] += coefU[i][1]; 
            QWV[i][2] = RL1 * (RL2*coefD[i][2] + coefU[i][2]);
        } else {
            grt_cmat2x1_mul(R2, coefU[i], QWV[i]);
            QWV[i][0] += coefD[i][0];
            QWV[i][1] += coefD[i][1]; 
            QWV[i][2] = RL1 * (coefD[i][2] + RL2*coefU[i][2]);
        }
        grt_cmat2x1_mul(R1, QWV[i], QWV[i]);
    }
}




// ================ 动态解 =====================
#define __DYNAMIC_KERNEL__ 
    #include "kernel_template.c_"
#undef __DYNAMIC_KERNEL__

// ================ 静态解 =====================
#define __STATIC_KERNEL__ 
    #include "kernel_template.c_"
#undef __STATIC_KERNEL__
