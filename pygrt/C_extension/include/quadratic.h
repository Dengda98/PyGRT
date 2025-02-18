/**
 * @file   quadratic.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 已知x1,x2,x3三点以及对应的函数值(复数)，确定这三点构成的二次函数的系数， 
 * \f[
 * f(x) = ax^2 + bx + c
 * \f]
 *                   
 */

#pragma once 

#include "const.h"


/**
 * 已知三个点x1,x2,x3,以及对应的函数值f1,f2,f3, 拟合函数 
 * 
 * 
 * @param    x[3]      (in)自变量 
 * @param    f[3]      (in)因变量 
 * @param    pa       (out)拟合a值
 * @param    pb       (out)拟合b值
 * @param    pc       (out)拟合c值
 */
void quad_term(const MYREAL x[3], const MYCOMPLEX f[3], MYCOMPLEX *pa, MYCOMPLEX *pb, MYCOMPLEX *pc);


/**
 * 给定x，根据a,b,c值，估计 \f$ f(x) \f$
 * 
 * @param     x    (in)自变量 
 * @param     a    (in)a值
 * @param     b    (in)b值
 * @param     c    (in)c值
 * @return    \f$ f(x) = ax^2 + bx + c \f$
 * 
 */
MYCOMPLEX quad_eval(MYREAL x, MYCOMPLEX a, MYCOMPLEX b, MYCOMPLEX c);


/**
 * 给定x，根据a,b,c值，估计 \f$ \MYINT_{x_1}^{x_2} f(s)ds \f$
 * 
 * @param     x1   (in)积分下限
 * @param     x2   (in)积分上限
 * @param     a    (in)a值
 * @param     b    (in)b值
 * @param     c    (in)c值
 * @return    \f$  \frac{1}{3}a(x_2^3-x_1^3) + \frac{1}{2}b(x_2^2-x_1^2) + c(x_2-x_1) \f$ 
 * 
 */
MYCOMPLEX quad_integral(MYREAL x1, MYREAL x2, MYCOMPLEX a, MYCOMPLEX b, MYCOMPLEX c);