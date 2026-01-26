/**
 * @file   attenuation.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 
 */


#include "grt/common/attenuation.h"
#include "grt/common/const.h"



cplx_t grt_attenuation_law(real_t Qinv, cplx_t omgref, cplx_t omega){
    return 1.0 + Qinv/PI * log(omega/omgref) + 0.5*Qinv*I;
}