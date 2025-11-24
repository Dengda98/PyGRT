/**
 * @file   attenuation.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 
 */


#include "grt/common/attenuation.h"
#include "grt/common/const.h"



cplx_t grt_attenuation_law(real_t Qinv, cplx_t omega){
    return 1.0 + Qinv/PI * log(omega/PI2) + 0.5*Qinv*I;
    // return 1.0;
}

void grt_py_attenuation_law(real_t Qinv, real_t omg[2], real_t atte[2]){
    // 用于在python中调用attenuation_law
    cplx_t omega = omg[0] + I*omg[1];
    cplx_t atte0 = grt_attenuation_law(Qinv, omega);
    atte[0] = creal(atte0);
    atte[1] = cimag(atte0);
}