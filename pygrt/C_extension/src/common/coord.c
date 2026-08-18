/**
 * @file   coord.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-10
 * 
 * 关于坐标变换的一些函数
 * 
 */

#include <stdbool.h>
#include <tgmath.h>

#include "grt/common/coord.h"

void grt_rot_zxy2zrt_vec(real_t theta, real_t A[3]){
    real_t s1, s2, s3;
    s1 = A[0];  s2 = A[1];  s3 = A[2];
    real_t st = sin(theta);
    real_t ct = cos(theta);
    A[0] = s1;
    A[1] = s2*ct + s3*st;
    A[2] = -s2*st + s3*ct;
}



void grt_rot_zxy2zrt_symtensor2odr(real_t theta, real_t A[6]) {
    real_t s11, s12, s13, s22, s23, s33;
    s11 = A[0];   s12 = A[1];   s13 = A[2];
                  s22 = A[3];   s23 = A[4];
                                s33 = A[5];
    real_t st = sin(theta);
    real_t ct = cos(theta);
    real_t sst = st*st;
    real_t cct = ct*ct;
    real_t sct = st*ct;
    A[0] = s11;
    A[1] = s12*ct + s13*st;
    A[2] = -s12*st + s13*ct;
    A[3] = s22*cct + s33*sst + 2.0*s23*sct;
    A[4] = (s33 - s22)*sct + s23*(cct - sst);
    A[5] = s22*sst + s33*cct - 2.0*s23*sct;
    
}



void grt_rot_zrt2zxy_upar(const real_t theta, real_t u[3], real_t upar[3][3], const real_t r){
    real_t s00, s01, s02;
    real_t s10, s11, s12;
    real_t s20, s21, s22;
    //           uz       ur       ut
    //  ∂z
    //  ∂r
    //  1/r*∂t
    s00 = upar[0][0]; s01 = upar[0][1]; s02 = upar[0][2];
    s10 = upar[1][0]; s11 = upar[1][1]; s12 = upar[1][2];
    s20 = upar[2][0]; s21 = upar[2][1]; s22 = upar[2][2];

    real_t u0, u1, u2;
    u0 = u[0];  u1 = u[1];  u2 = u[2];

    real_t st = sin(theta);
    real_t ct = cos(theta);
    real_t sst = st*st;
    real_t cct = ct*ct;
    real_t sct = st*ct;

    // 变换含联络项 u_r/r、u_θ/r（r 单位 cm）。
    // r=0: u/r 联络项改用 ∂_r u_r、∂_r u_θ（s11,s12），与 syn 中 (1/r)∂_θ 有限部分配套。
    real_t u1_over_r, u2_over_r;
    if(GRT_IS_ZERO(r * 1e-5)){  // cm → km 后再判零
        u1_over_r = s11;
        u2_over_r = s12;
    } else {
        u1_over_r = u1/r;
        u2_over_r = u2/r;
    }

    //           uz       ux       uy
    //  ∂z
    //  ∂x
    //  ∂y

    // ∂ uz / ∂ z
    upar[0][0] = s00;
    // ∂ ux / ∂ z
    upar[0][1] = s01*ct - s02*st;
    // ∂ uy / ∂ z
    upar[0][2] = s01*st + s02*ct;


    // ∂ uz / ∂ x
    upar[1][0] = s10*ct - s20*st;
    // ∂ ux / ∂ x
    upar[1][1] = s11*cct + s22*sst - (s12+s21)*sct + u1_over_r*sst + u2_over_r*sct;
    // ∂ uy / ∂ x
    upar[1][2] = s12*cct - s21*sst + (s11-s22)*sct - u1_over_r*sct + u2_over_r*sst;


    // ∂ uz / ∂ y
    upar[2][0] = s10*st + s20*ct;
    // ∂ ux / ∂ y
    upar[2][1] = s21*cct - s12*sst + (s11-s22)*sct - u1_over_r*sct - u2_over_r*cct;
    // ∂ uy / ∂ y
    upar[2][2] = s22*cct + s11*sst + (s12+s21)*sct + u1_over_r*cct - u2_over_r*sct;


    // 转矢量
    u[0] = u0;
    u[1] = u1*ct - u2*st;
    u[2] = u1*st + u2*ct;
}



/**
 * 直角坐标 zxy 到柱坐标 zrt 的位移及位移偏导旋转
 *
 * 该函数与 grt_rot_zrt2zxy_upar 互为相反方向的坐标变换
 */
void grt_rot_zxy2zrt_upar(const real_t theta, real_t u[3], real_t upar[3][3], const real_t r){
    real_t s00, s01, s02;
    real_t s10, s11, s12;
    real_t s20, s21, s22;

    // 保存直角坐标下的位移偏导，行表示偏导坐标，列表示位移分量
    s00 = upar[0][0];  s01 = upar[0][1];  s02 = upar[0][2];
    s10 = upar[1][0];  s11 = upar[1][1];  s12 = upar[1][2];
    s20 = upar[2][0];  s21 = upar[2][1];  s22 = upar[2][2];

    real_t st = sin(theta);
    real_t ct = cos(theta);
    real_t sst = st * st;
    real_t cct = ct * ct;
    real_t sct = st * ct;

    // 先将水平偏导的偏导方向和位移分量同时旋转到 r、theta
    real_t hrr = s11 * cct + (s12 + s21) * sct + s22 * sst;
    real_t hrt = -s11 * sct + s12 * cct - s21 * sst + s22 * sct;
    real_t htr = -s11 * sct - s12 * sst + s21 * cct + s22 * sct;
    real_t htt = s11 * sst - (s12 + s21) * sct + s22 * cct;

    // 先旋转位移矢量，后续联络项需要使用 ur、ut 分量
    grt_rot_zxy2zrt_vec(theta, u);

    // 联络项使用 cm 作为 r 的长度单位，与位移的 cm 单位匹配
    real_t ur_over_r, ut_over_r;
    if(GRT_IS_ZERO(r * 1e-5)){  // cm → km 后再判零
        ur_over_r = hrr;
        ut_over_r = hrt;
    } else {
        ur_over_r = u[1] / r;
        ut_over_r = u[2] / r;
    }

    // z 方向偏导只需要旋转位移分量
    upar[0][0] = s00;
    upar[0][1] = s01 * ct + s02 * st;
    upar[0][2] = -s01 * st + s02 * ct;

    // r、theta 方向对 z 分量的偏导只需要旋转偏导方向
    upar[1][0] = s10 * ct + s20 * st;
    upar[2][0] = -s10 * st + s20 * ct;

    // 柱坐标偏导中补充位移基矢变化产生的联络项
    upar[1][1] = hrr;
    upar[1][2] = hrt;
    upar[2][1] = htr + ut_over_r;
    upar[2][2] = htt - ur_over_r;

}