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

#include "common/quadratic.h"
#include "common/const.h"



void quad_term(const MYREAL x[3], const MYCOMPLEX f[3], MYCOMPLEX *pa, MYCOMPLEX *pb, MYCOMPLEX *pc)
{
    MYREAL x1, x2, x3, w1, w2, w3;
    x1 = x[0];
    x2 = x[1];
    x3 = x[2];

    w1 = x1*x1;
    w2 = x2*x2;
    w3 = x3*x3;

    MYCOMPLEX f1, f2, f3;
    f1 = f[0];
    f2 = f[1];
    f3 = f[2];

    MYCOMPLEX d;
    d = (x1-x2)*(x2-x3)*(x3-x1);
    *pa = - (f1*(x2-x3) + f2*(x3-x1) + f3*(x1-x2)) / d;
    *pb =   (f1*(w2-w3) + f2*(w3-w1) + f3*(w1-w2)) / d;
    *pc = - (f1*x2*x3*(x2-x3) + f2*x3*x1*(x3-x1) + f3*x1*x2*(x1-x2)) / d;
}



MYCOMPLEX quad_eval(MYREAL x, MYCOMPLEX a, MYCOMPLEX b, MYCOMPLEX c)
{
    return a*x*x + b*x + c;
}


MYCOMPLEX quad_integral(MYREAL x1, MYREAL x2, MYCOMPLEX a, MYCOMPLEX b, MYCOMPLEX c)
{
    MYREAL xx1, xx2, xxx1, xxx2;
    xx1 = x1*x1;    xx2 = x2*x2;
    xxx1 = xx1*x1;  xxx2 = xx2*x2; 
    return a/RTHREE*(xxx2 - xxx1) + b/RTWO*(xx2 - xx1) + c*(x2 - x1);
}


