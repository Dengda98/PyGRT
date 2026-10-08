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
 * 使用参考点处的局部切平面近似，将局部北向/东向坐标转换为纬度/经度
 *
 * 地球半径取 6371 km，输出经度规范化到 [-180, 180)
 * 转换后的纬度超出 [-90, 90] 时终止程序并报错
 *
 * @param[in]  north 北向坐标，单位为 km
 * @param[in]  east  东向坐标，单位为 km
 * @param[in]  lat0  参考点纬度，单位为度，必须在 (-90, 90) 内
 * @param[in]  lon0  参考点经度，单位为度，必须在 [-180, 180] 内
 * @param[out] lat   转换后的纬度，单位为度
 * @param[out] lon   转换后的经度，单位为度
 */
void grt_xy2geo(real_t north, real_t east, real_t lat0, real_t lon0, real_t *lat, real_t *lon);


/**
 * 使用参考点处的局部切平面近似，将纬度/经度转换为局部北向/东向坐标
 *
 * 地球半径取 6371 km，经度差规范化到 [-180, 180)
 * 输入纬度超出 [-90, 90] 时终止程序并报错
 *
 * @param[in]  lat   纬度，单位为度
 * @param[in]  lon   经度，单位为度
 * @param[in]  lat0  参考点纬度，单位为度，必须在 (-90, 90) 内
 * @param[in]  lon0  参考点经度，单位为度，必须在 [-180, 180] 内
 * @param[out] north 转换后的北向坐标，单位为 km
 * @param[out] east  转换后的东向坐标，单位为 km
 */
void grt_geo2xy(real_t lat, real_t lon, real_t lat0, real_t lon0, real_t *north, real_t *east);


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

/**
 * 根据断层三要素构造 N、E、Z-up 顺序的平面法向和面内切向单位矢量
 * @param[in]  strike  走向，度
 * @param[in]  dip     倾角，度
 * @param[in]  rake    滑动角，度
 * @param[out] nvec    平面法向单位矢量
 * @param[out] tvec    沿滑动角方向的面内切向单位矢量
 */
void grt_fault_plane_vectors(real_t strike, real_t dip, real_t rake, real_t nvec[3], real_t tvec[3]);

/**
 * 将 ZNE 应力张量投影为接收断层面上的法向应力和沿滑动方向的剪应力
 * @param[in]  stress   ZZ、ZN、ZE、NN、NE、EE
 * @param[in]  nvec     N、E、Z 顺序的法向量
 * @param[in]  tvec     N、E、Z 顺序的滑动方向
 * @param[out] sigma_n  法向应力，张为正
 * @param[out] tau_s    滑动方向剪应力
 */
void grt_project_stress_to_fault_plane(const real_t stress[6], const real_t nvec[3], const real_t tvec[3], real_t *sigma_n, real_t *tau_s);
