/**
 * @file   quadratic.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 已知x1,x2,x3三点以及对应的函数值(复数)，确定这三点构成的二次函数的系数， 
 * \f[
 * f(x) = ax^2 + bx + c
 * \f]
 *       
 */

#include <stdio.h>

#include "grt/integral/quadratic.h"
#include "grt/common/const.h"



void grt_quad_term(const real_t x[3], const cplx_t f[3], cplx_t *pa, cplx_t *pb, cplx_t *pc)
{
    real_t x1, x2, x3, w1, w2, w3;
    x1 = x[0];
    x2 = x[1];
    x3 = x[2];

    w1 = x1*x1;
    w2 = x2*x2;
    w3 = x3*x3;

    cplx_t f1, f2, f3;
    f1 = f[0];
    f2 = f[1];
    f3 = f[2];

    cplx_t d;
    d = (x1-x2)*(x2-x3)*(x3-x1);
    *pa = - (f1*(x2-x3) + f2*(x3-x1) + f3*(x1-x2)) / d;
    *pb =   (f1*(w2-w3) + f2*(w3-w1) + f3*(w1-w2)) / d;
    *pc = - (f1*x2*x3*(x2-x3) + f2*x3*x1*(x3-x1) + f3*x1*x2*(x1-x2)) / d;
}



cplx_t grt_quad_eval(real_t x, cplx_t a, cplx_t b, cplx_t c)
{
    return a*x*x + b*x + c;
}


cplx_t grt_quad_integral(real_t x1, real_t x2, cplx_t a, cplx_t b, cplx_t c)
{
    real_t xx1, xx2, xxx1, xxx2;
    xx1 = x1*x1;    xx2 = x2*x2;
    xxx1 = xx1*x1;  xxx2 = xx2*x2; 
    return a/3.0*(xxx2 - xxx1) + b/2.0*(xx2 - xx1) + c*(x2 - x1);
}


