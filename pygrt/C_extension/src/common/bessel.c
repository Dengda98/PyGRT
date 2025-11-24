/**
 * @file   bessel.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 
 */


#include "grt/common/bessel.h"
#include "grt/common/const.h"

void grt_bessel012(real_t x, real_t *bj0, real_t *bj1, real_t *bj2){
    *bj0 = j0(x);
    *bj1 = j1(x);
    *bj2 = jn(2, x);
}

void grt_besselp012(real_t x, real_t *bj0, real_t *bj1, real_t *bj2){
    real_t j0=*bj0;
    real_t j1=*bj1;
    real_t j2=*bj2;
    *bj0 = -j1;
    *bj1 = j0 - 1.0/x * j1;
    *bj2 = j1 - 2.0/x * j2;
}
