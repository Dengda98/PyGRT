/**
 * @file   bessel.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 
 */


#include "bessel.h"
#include "const.h"

void bessel012(MYREAL x, MYREAL *bj0, MYREAL *bj1, MYREAL *bj2){
    *bj0 = J0(x);
    *bj1 = J1(x);
    *bj2 = JN(2, x);
}

