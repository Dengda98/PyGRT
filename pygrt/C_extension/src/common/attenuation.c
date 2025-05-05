/**
 * @file   attenuation.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 
 */


#include "common/attenuation.h"
#include "common/const.h"



MYCOMPLEX attenuation_law(MYREAL Qinv, MYCOMPLEX omega){
    return RONE + Qinv/PI * log(omega/PI2) + RHALF*Qinv*I;
    // return RONE;
}

void py_attenuation_law(MYREAL Qinv, MYREAL omg[2], MYREAL atte[2]){
    // 用于在python中调用attenuation_law
    MYCOMPLEX omega = omg[0] + I*omg[1];
    MYCOMPLEX atte0 = attenuation_law(Qinv, omega);
    atte[0] = creal(atte0);
    atte[1] = cimag(atte0);
}