/**
 * @file   bessel.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 
 */


#include "common/bessel.h"
#include "common/const.h"

void bessel012(MYREAL x, MYREAL *bj0, MYREAL *bj1, MYREAL *bj2){
    *bj0 = j0(x);
    *bj1 = j1(x);
    *bj2 = jn(2, x);
}

void besselp012(MYREAL x, MYREAL *bj0, MYREAL *bj1, MYREAL *bj2){
    MYREAL j0=*bj0;
    MYREAL j1=*bj1;
    MYREAL j2=*bj2;
    *bj0 = -j1;
    *bj1 = j0 - RONE/x * j1;
    *bj2 = j1 - RTWO/x * j2;
}
