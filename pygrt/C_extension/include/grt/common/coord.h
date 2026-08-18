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
 * @param[in,out]   upar       柱坐标下的位移空间偏导（第三行已是 (1/r)∂_θ 有限部分）
 * @param[in]       r          r 坐标 (cm)；r=0 时联络项 u/r 改用 ∂_r u
 */
void grt_rot_zrt2zxy_upar(const real_t theta, real_t u[3], real_t upar[3][3], const real_t r);


/**
 * 直角坐标 zxy 到柱坐标 zrt 的位移及位移偏导旋转
 *
 * 输入偏导矩阵为直角坐标分量对 z、x、y 的偏导
 * 输出偏导矩阵为柱坐标分量对 z、r、theta 的偏导，第三行为 (1/r)∂_theta 对柱坐标位移分量的偏导
 * 变换过程中包含位移基矢变化产生的联络项
 *
 * @param[in]       theta      r 轴相对 x 轴的旋转弧度
 * @param[in,out]   u          待旋转的位移矢量
 * @param[in,out]   upar       待旋转的位移偏导矩阵
 * @param[in]       r          r 坐标，单位为 cm；r=0 时使用轴线上有限极限
 */
void grt_rot_zxy2zrt_upar(const real_t theta, real_t u[3], real_t upar[3][3], const real_t r);
