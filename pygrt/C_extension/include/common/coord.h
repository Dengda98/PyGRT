/**
 * @file   coord.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-10
 * 
 * 关于坐标变换的一些函数
 * 
 */

#pragma once 

#include <stdbool.h>

/**
 * 直角坐标zxy到柱坐标zrt的矢量旋转
 * 
 * @param   theta        (in)r轴相对x轴的旋转弧度(负数表示逆变换，即zrt->zxy)
 * @param   A[3]         (inout)待旋转的矢量(s1, s2, s3)
 */
void rot_zxy2zrt_vec(double theta, double A[3]);



/**
 * 直角坐标zxy到柱坐标zrt的二阶对称张量旋转
 * 
 * @param   theta        (in)r轴相对x轴的旋转弧度(负数表示逆变换，即zrt->zxy)
 * @param   A[6]         (inout)待旋转的二阶对称张量(s11, s12, s13, s22, s23, s33)
 */
void rot_zxy2zrt_symtensor2odr(double theta, double A[6]);


/**
 * 柱坐标下的位移偏导 ∂u(z,r,t)/∂(z,r,t) 转到 直角坐标 ∂u(z,x,y)/∂(z,x,y)
 * 
 *            uz       ur       ut                 uz       ux       uy
 *  ∂z                                       ∂z
 *  ∂r                                 --->  ∂x
 *  1/r*∂t                                   ∂y
 * 
 * 
 * @param   theta        (in)r轴相对x轴的旋转弧度
 * @param   u[3]         (inout)柱坐标下的位移矢量
 * @param   upar[3][3]   (inout)柱坐标下的位移空间偏导
 * @param   r            (in)r轴坐标
 */
void rot_zrt2zxy_upar(const double theta, double u[3], double upar[3][3], const double r);