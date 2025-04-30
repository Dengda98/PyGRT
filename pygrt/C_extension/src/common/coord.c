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



void rot_zxy2zrt_vec(double theta, double A[3]){
    double s1, s2, s3;
    s1 = A[0];  s2 = A[1];  s3 = A[2];
    double st = sin(theta);
    double ct = cos(theta);
    A[0] = s1;
    A[1] = s2*ct + s3*st;
    A[2] = -s2*st + s3*ct;
}



void rot_zxy2zrt_symtensor2odr(double theta, double A[6]) {
    double s11, s12, s13, s22, s23, s33;
    s11 = A[0];   s12 = A[1];   s13 = A[2];
                  s22 = A[3];   s23 = A[4];
                                s33 = A[5];
    double st = sin(theta);
    double ct = cos(theta);
    double sst = st*st;
    double cct = ct*ct;
    double sct = st*ct;
    A[0] = s11;
    A[1] = s12*ct + s13*st;
    A[2] = -s12*st + s13*ct;
    A[3] = s22*cct + s33*sst + 2.0*s23*sct;
    A[4] = (s33 - s22)*sct + s23*(cct - sst);
    A[5] = s22*sst + s33*cct - 2.0*s23*sct;
    
}



void rot_zrt2zxy_upar(const double theta, double u[3], double upar[3][3], const double r){
    double s00, s01, s02;
    double s10, s11, s12;
    double s20, s21, s22;
    //           uz       ur       ut
    //  ∂z
    //  ∂r
    //  1/r*∂t
    s00 = upar[0][0]; s01 = upar[0][1]; s02 = upar[0][2];
    s10 = upar[1][0]; s11 = upar[1][1]; s12 = upar[1][2];
    s20 = upar[2][0]; s21 = upar[2][1]; s22 = upar[2][2];

    double u0, u1, u2;
    u0 = u[0];  u1 = u[1];  u2 = u[2];

    double st = sin(theta);
    double ct = cos(theta);
    double sst = st*st;
    double cct = ct*ct;
    double sct = st*ct;

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
    upar[1][1] = s11*cct + s22*sst - (s12+s21)*sct + u1*sst/r + u2*sct/r;
    // ∂ uy / ∂ x
    upar[1][2] = s12*cct - s21*sst + (s11-s22)*sct - u1*sct/r + u2*sst/r;


    // ∂ uz / ∂ y
    upar[2][0] = s10*st + s20*ct;
    // ∂ ux / ∂ y
    upar[2][1] = s21*cct - s12*sst + (s11-s22)*sct - u1*sct/r - u2*cct/r;
    // ∂ uy / ∂ y
    upar[2][2] = s22*cct + s11*sst + (s12+s21)*sct + u1*cct/r - u2*sct/r;


    // 转矢量
    u[0] = u0;
    u[1] = u1*ct - u2*st;
    u[2] = u1*st + u2*ct;
}