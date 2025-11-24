/**
 * @file   elliptic.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    使用渐近解和 Carlson 算法计算三类完全椭圆积分，参考：
 *    
 *        1. Carlson, B. C. 1995. “Numerical Computation of Real or Complex Elliptic Integrals.” 
 *           Numerical Algorithms 10 (1): 13–26. https://doi.org/10.1007/BF02198293.
 *        2. 张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社
 */


#include <stdio.h>
#include <tgmath.h>
#include <stdbool.h>
#include <complex.h>

#include "grt/lamb/elliptic.h"
#include "grt/common/checkerror.h"


real_t grt_ellipticK(const real_t m)
{   
    if(m <= 0.0 || m>= 1.0){
        GRTRaiseError("For the first complete elliptic integral, m should be in (0,1), but get %f.\n", m);
    }

    #define N 5

    real_t m0 = 1 - m;

    static const real_t a[N] = {
        1.38629436112,
        0.09666344259,
        0.03590092383,
        0.03742563713,
        0.01451196212,
    };

    static const real_t b[N] = {
        0.50000000000,
        0.12498593597,
        0.06880248576,
        0.03328355346,
        0.00441787012,
    };

    real_t K1 = 0.0, K2 = 0.0;
    real_t p = 1.0;
    K1 += a[0];
    K2 += b[0];
    for(int i = 1; i < N; ++i){
        p *= m0;
        K1 += a[i] * p;
        K2 += b[i] * p;
    }

    return K1 + K2*log(1.0/m0);
    #undef N
}

real_t grt_ellipticE(const real_t m)
{   
    if(m <= 0.0 || m>= 1.0){
        GRTRaiseError("For the second complete elliptic integral, m should be in (0,1), but get %f.\n", m);
    }

    #define N 5

    real_t m0 = 1 - m;

    static const real_t a[N] = {
        1.00000000000,
        0.44325141463,
        0.06260601220,
        0.04757383546,
        0.01736506451,
    };

    static const real_t b[N] = {
        0.00000000000,
        0.24998368310,
        0.09200180037,
        0.04069697526,
        0.00526449639
    };

    real_t K1 = 0.0, K2 = 0.0;
    real_t p = 1.0;
    K1 += a[0];
    K2 += b[0];
    for(int i = 1; i < N; ++i){
        p *= m0;
        K1 += a[i] * p;
        K2 += b[i] * p;
    }

    return K1 + K2*log(1.0/m0);
    #undef N
}

static real_t MAX3(const real_t x, const real_t y, const real_t z)
{
    real_t xy = ((x) > (y) ? (x) : (y));
    return ((xy) > (z) ? (xy) : (z));
}

static real_t MAX4(const real_t x, const real_t y, const real_t z, const real_t p)
{
    real_t xy = ((x) > (y) ? (x) : (y));
    real_t xyz = ((xy) > (z) ? (xy) : (z));
    return ((xyz) > (p) ? (xyz) : (p));
}

static cplx_t grt_ellipticRF(const cplx_t x, const cplx_t y, const cplx_t z)
{
    real_t errtol = 0.001;
    int nmax = 10000;

    real_t eps;
    cplx_t lam;
    cplx_t X, Y, Z, An;
    cplx_t xn, yn, zn;
    cplx_t xr, yr, zr;
    int n = 0;
    xn = x;
    yn = y;
    zn = z;
    while (true) {
        An = (xn + yn + zn) / 3.0;
        X = 2.0 - (An + xn) / An;
        Y = 2.0 - (An + yn) / An;
        Z = 2.0 - (An + zn) / An;
        eps = MAX3(fabs(X), fabs(Y), fabs(Z));
        if (eps < errtol)  break;

        xr = sqrt(xn);
        yr = sqrt(yn);
        zr = sqrt(zn);
        lam = xr*yr + xr*zr + yr*zr;
        xn = (xn + lam) * 0.25;
        yn = (yn + lam) * 0.25;
        zn = (zn + lam) * 0.25;
        n++;
        if(n >= nmax)  break;
    }

    cplx_t E2, E3;
    E2 = X*Y - Z*Z;
    E3 = X*Y*Z;

    cplx_t R;
    R = 1.0 - 0.1*E2 + 1/14.0 * E3 + 1/24.0 * E2*E2 - 3/44.0 *E2*E3;
    R /= sqrt(An);

    return R;
}


static cplx_t grt_ellipticRJ(const cplx_t x, const cplx_t y, const cplx_t z, const cplx_t p)
{
    real_t errtol = 0.001;
    int nmax = 10000;

    real_t eps;
    cplx_t lam, dm, em;
    cplx_t X, Y, Z, P, An;
    cplx_t xn, yn, zn, pn;
    cplx_t xr, yr, zr, pr;
    cplx_t sum1 = 0.0;
    real_t pow4 = 1.0;
    int n = 0;
    xn = x;
    yn = y;
    zn = z;
    pn = p;

    while (true) {
        An = (xn + yn + zn + 2.0*pn) * 0.2;
        X = 2.0 - (An + xn) / An;
        Y = 2.0 - (An + yn) / An;
        Z = 2.0 - (An + zn) / An;
        P = 2.0 - (An + zn) / An;
        eps = MAX4(fabs(X), fabs(Y), fabs(Z), fabs(P));
        if (eps < errtol)  break;

        xr = sqrt(xn);
        yr = sqrt(yn);
        zr = sqrt(zn);
        pr = sqrt(pn);
        lam = xr*yr + xr*zr + yr*zr;
        dm = (pr + xr)*(pr + yr)*(pr + zr);
        em = (pn - xn)*(pn - yn)*(pn - zn) / (dm*dm);
        sum1 += grt_ellipticRF(1.0, 1.0+em, 1.0+em)/dm * pow4;
        pow4 *= 0.25;

        xn = (xn + lam) * 0.25;
        yn = (yn + lam) * 0.25;
        zn = (zn + lam) * 0.25;
        pn = (pn + lam) * 0.25;
        n++;
        if(n >= nmax)  break;
    }

    cplx_t E2, E3, E4, E5;
    cplx_t PPP = P*P*P;
    E2 = X*Y + Y*Z + X*Z - 3.0*P*P;
    E3 = X*Y*Z + 3.0*E2*P + 4.0*PPP;
    E4 = (2.0*X*Y*Z + E2*P + 3.0*PPP)*P;
    E5 = X*Y*Z*P*P;

    cplx_t res = 0.0;
    res = 1.0 - 3.0/14 * E2 + 1.0/6 * E3 + 9.0/88*E2*E2 - 3.0/22*E4 - 9.0/52*E2*E3 + 3.0/26*E5;
    res /= An * sqrt(An);
    res *= pow4;

    res += 6.0 * sum1;

    return res;
}


cplx_t grt_ellipticPi(const cplx_t n, const real_t m)
{
    if(m <= 0.0 || m>= 1.0){
        GRTRaiseError("For the third complete elliptic integral, m should be in (0,1), but get %f.\n", m);
    }

    if(cimag(n) == 0.0 && creal(n) > 1.0){
        return grt_ellipticK(m) - grt_ellipticPi(m/n, m);
    } else {
        return grt_ellipticK(m) + n/3.0 * grt_ellipticRJ(0.0, 1.0 - m, 1.0, 1.0 - n);
    }
}

