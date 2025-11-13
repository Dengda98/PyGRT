/**
 * @file   lamb_util.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    一些使用广义闭合解求解 Lamb 问题过程中可能用到的辅助函数
 */

#include "grt/lamb/lamb_util.h"

void grt_roots3(double a, double b, double c, MYCOMPLEX y3[3])
{
    double Q, R;
    Q = (a*a - 3.0*b) / 9.0;
    R = (2.0*a*a*a - 9.0*a*b + 27.0*c) / 54.0;

    double Q3, R2;
    R2 = R*R;
    Q3 = Q*Q*Q;

    y3[0] = y3[1] = y3[2] = 0.0;

    if (Q > 0.0 && Q3 > R2){
        double theta;
        theta = acos(R / sqrt(Q3));
        y3[0] = - 2.0 * sqrt(Q) * cos(theta/3.0) - a/3.0;
        y3[1] = - 2.0 * sqrt(Q) * cos((theta - M_PI*2.0)/3.0) - a/3.0;
        y3[2] = - 2.0 * sqrt(Q) * cos((theta + M_PI*2.0)/3.0) - a/3.0;
    } else {
        double A = pow(fabs(R) + sqrt(R2 - Q3), 1/3.0);
        A = (R > 0.0) ? - A : A;
        double B = (A != 0.0)? Q/A : 0.0;

        y3[0] = - 0.5 * (A+B) -  a/3.0 + I*sqrt(3.0)/2.0*(A - B);
        y3[1] = - 0.5 * (A+B) -  a/3.0 - I*sqrt(3.0)/2.0*(A - B);
        y3[2] = A + B - a/3.0;
    }
}

MYCOMPLEX grt_evalpoly2(const MYCOMPLEX *C, const int n, const MYCOMPLEX y, const int offset)
{
    MYCOMPLEX res = 0.0;
    MYCOMPLEX p = 1.0;
    for(int i=0; i<=n; ++i){
        res += C[2*i+offset] * p;
        p *= y;
    }
    return res;
}