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

#include "grt/common/const.h"

/**
 * 直角坐标zxy到柱坐标zrt的矢量旋转
 * 
 * @param[in]    theta        r轴相对x轴的旋转弧度(负数表示逆变换，即zrt->zxy)
 * @param[out]   A            待旋转的矢量(s1, s2, s3)
 */
void grt_rot_zxy2zrt_vec(real_t theta, real_t A[3]);



/**
 * 直角坐标zxy到柱坐标zrt的二阶对称张量旋转
 * 
 * @param[in]    theta       r轴相对x轴的旋转弧度(负数表示逆变换，即zrt->zxy)
 * @param[out]   A           待旋转的二阶对称张量(s11, s12, s13, s22, s23, s33)
 */
void grt_rot_zxy2zrt_symtensor2odr(real_t theta, real_t A[6]);


/**
 * 柱坐标下的位移偏导 ∂u(z,r,t)/∂(z,r,t) 转到 直角坐标 ∂u(z,x,y)/∂(z,x,y)
 * 
 * |          |    uz     |     ur    |     ut    |
 * |----------|-----------|-----------|-----------|
 * |    ∂z    |           |           |           |
 * |    ∂r    |           |           |           |
 * |  1/r*∂t  |           |           |           |
 * 
 * 
 * |          |    uz     |     ux    |     uy    |
 * |----------|-----------|-----------|-----------|
 * |    ∂z    |           |           |           |
 * |    ∂x    |           |           |           |
 * |    ∂y    |           |           |           |
 * 
 * 
 * 
 * @param[in]       theta      r轴相对x轴的旋转弧度
 * @param[in,out]   u          柱坐标下的位移矢量
 * @param[in,out]   upar       柱坐标下的位移空间偏导
 * @param[in]       r          r轴坐标
 */
void grt_rot_zrt2zxy_upar(const real_t theta, real_t u[3], real_t upar[3][3], const real_t r);