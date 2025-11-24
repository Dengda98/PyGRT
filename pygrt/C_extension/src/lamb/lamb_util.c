/**
 * @file   lamb_util.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    一些使用广义闭合解求解 Lamb 问题过程中可能用到的辅助函数
 */

#include "grt/lamb/lamb_util.h"


/**
 * 求解一元三次方程的根， \f$ x^3 + ax^2 + bx + c = 0 \f$
 * 
 * @param[in]      a     系数 a
 * @param[in]      b     系数 b
 * @param[in]      c     系数 c
 * @param[out]     y3    三个根，其中 y3[2] 为正根
 */
static void _roots3(real_t a, real_t b, real_t c, cplx_t y3[3])
{
    real_t Q, R;
    Q = (a*a - 3.0*b) / 9.0;
    R = (2.0*a*a*a - 9.0*a*b + 27.0*c) / 54.0;

    real_t Q3, R2;
    R2 = R*R;
    Q3 = Q*Q*Q;

    y3[0] = y3[1] = y3[2] = 0.0;

    if (Q > 0.0 && Q3 > R2){
        real_t theta;
        theta = acos(R / sqrt(Q3));
        y3[0] = - 2.0 * sqrt(Q) * cos(theta/3.0) - a/3.0;
        y3[1] = - 2.0 * sqrt(Q) * cos((theta - M_PI*2.0)/3.0) - a/3.0;
        y3[2] = - 2.0 * sqrt(Q) * cos((theta + M_PI*2.0)/3.0) - a/3.0;
    } else {
        real_t A = pow(fabs(R) + sqrt(R2 - Q3), 1/3.0);
        A = (R > 0.0) ? - A : A;
        real_t B = (A != 0.0)? Q/A : 0.0;

        y3[0] = - 0.5 * (A+B) -  a/3.0 + I*sqrt(3.0)/2.0*(A - B);
        y3[1] = - 0.5 * (A+B) -  a/3.0 - I*sqrt(3.0)/2.0*(A - B);
        y3[2] = A + B - a/3.0;
    }
}

void grt_rayleigh1_roots(real_t nu, cplx_t y3[3])
{
    real_t a, b, c;
    real_t nu2, nu3, nu4;
    nu2 = nu*nu;
    nu3 = nu2*nu;
    nu4 = nu3*nu;
    real_t snu = 1.0 - nu;
    real_t snu2 = snu*snu;
    real_t snu3 = snu2*snu;
    a = -0.5 * (2.0*nu2 + 1.0)/snu;
    b = 0.25 * (4.0*nu3 - 4.0*nu2 + 4.0*nu - 1.0)/snu2;
    c = -0.125*nu4/snu3;
    _roots3(a, b, c, y3);
}

void grt_rayleigh2_roots(real_t m, cplx_t y3[3])
{
    real_t a, b, c;
    a = 0.5*(2.0*m - 3.0)/(1.0 - m);
    b = 0.5/(1.0 - m);
    c = - 0.0625/(1 - m);
    _roots3(a, b, c, y3);
}

cplx_t grt_evalpoly2(const cplx_t *C, const int n, const cplx_t y, const int offset)
{
    cplx_t res = 0.0;
    cplx_t p = 1.0;
    for(int i=0; i<=n; ++i){
        res += C[2*i+offset] * p;
        p *= y;
    }
    return res;
}