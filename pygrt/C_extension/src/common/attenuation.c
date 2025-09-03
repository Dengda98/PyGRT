/**
 * @file   attenuation.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 
 */


#include "grt/common/attenuation.h"
#include "grt/common/const.h"



MYCOMPLEX grt_attenuation_law(MYREAL Qinv, MYCOMPLEX omega){
    return 1.0 + Qinv/PI * log(omega/PI2) + 0.5*Qinv*I;
    // return 1.0;
}

void grt_py_attenuation_law(MYREAL Qinv, MYREAL omg[2], MYREAL atte[2]){
    // 用于在python中调用attenuation_law
    MYCOMPLEX omega = omg[0] + I*omg[1];
    MYCOMPLEX atte0 = grt_attenuation_law(Qinv, omega);
    atte[0] = creal(atte0);
    atte[1] = cimag(atte0);
}