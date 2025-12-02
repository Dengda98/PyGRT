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
#include "grt/common/kernel.h"

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


// ================ 动态解 =====================
#define __DYNAMIC_KERNEL__ 
    #include "kernel_template.c_"
#undef __DYNAMIC_KERNEL__

// ================ 静态解 =====================
#define __STATIC_KERNEL__ 
    #include "kernel_template.c_"
#undef __STATIC_KERNEL__
