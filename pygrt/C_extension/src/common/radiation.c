/**
 * @file   radiation.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-06
 * 
 *    计算不同震源的辐射因子
 * 
 */

#include <stdbool.h>

#include "common/radiation.h"
#include "common/const.h"

void set_source_radiation(
    double srcRadi[SRC_M_NUM][CHANNEL_NUM], const int computeType, const bool par_theta,
    const double M0, const double coef, const double azrad, const double mchn[MECHANISM_NUM]
){
    double mult;
    if(computeType == GRT_SYN_COMPUTE_SF){
        mult = 1e-15*M0*coef;
    } else {
        mult = 1e-20*M0*coef;
    }

    double saz, caz;
    saz = sin(azrad);
    caz = cos(azrad);

    if(computeType == GRT_SYN_COMPUTE_EX){
        srcRadi[0][0] = srcRadi[0][1] = (par_theta)? 0.0 : mult; // Z/R
        srcRadi[0][2] = 0.0; // T
    }  
    else if(computeType == GRT_SYN_COMPUTE_SF){
        double A0, A1, A4;
        double fn, fe, fz;
        fn=mchn[0];   fe=mchn[1];  fz=mchn[2];   
        A0 = fz*mult;
        A1 = (fn*caz + fe*saz)*mult;
        A4 = (- fn*saz + fe*caz)*mult;

        // 公式(4.6.20)
        srcRadi[1][0] = srcRadi[1][1] = (par_theta)? 0.0 : A0; // VF, Z/R
        srcRadi[2][0] = srcRadi[2][1] = (par_theta)? A4 : A1; // HF, Z/R
        srcRadi[1][2] = 0.0; // VF, T
        srcRadi[2][2] = (par_theta)? -A1 : A4; // HF, T
    }
    else if(computeType == GRT_SYN_COMPUTE_DC){
        double strike, dip, rake;
        strike=mchn[0];   dip=mchn[1];   rake=mchn[2];
        // 公式(4.8.35)
        double stkrad = strike*DEG1;
        double diprad = dip*DEG1;
        double rakrad = rake*DEG1;
        double therad = azrad - stkrad;
        double srak, crak, sdip, cdip, sdip2, cdip2, sthe, cthe, sthe2, cthe2;
        srak = sin(rakrad);     crak = cos(rakrad);
        sdip = sin(diprad);     cdip = cos(diprad);
        sdip2 = 2.0*sdip*cdip;  cdip2 = 2.0*cdip*cdip - 1.0;
        sthe = sin(therad);     cthe = cos(therad);
        sthe2 = 2.0*sthe*cthe;  cthe2 = 2.0*cthe*cthe - 1.0;

        double A0, A1, A2, A4, A5;
        A0 = mult * (0.5*sdip2*srak);
        A1 = mult * (cdip*crak*cthe - cdip2*srak*sthe);
        A2 = mult * (0.5*sdip2*srak*cthe2 + sdip*crak*sthe2);
        A4 = mult * (- cdip2*srak*cthe - cdip*crak*sthe);
        A5 = mult * (sdip*crak*cthe2 - 0.5*sdip2*srak*sthe2);

        srcRadi[3][0] = srcRadi[3][1] = (par_theta)? 0.0 : A0; // DD, Z/R
        srcRadi[4][0] = srcRadi[4][1] = (par_theta)? A4 : A1; // DS, Z/R
        srcRadi[5][0] = srcRadi[5][1] = (par_theta)? 2.0*A5 : A2; // SS, Z/R
        srcRadi[3][2] = 0.0; // DD, T
        srcRadi[4][2] = (par_theta)? -A1 : A4;  // DS, T
        srcRadi[5][2] = (par_theta)? -2.0*A2 : A5;  // DS, T
    }
    else if(computeType == GRT_SYN_COMPUTE_MT){
        // 公式(4.9.7)但修改了各向同性的量
        double M11, M12, M13, M22, M23, M33;
        M11 = mchn[0];   M12 = mchn[1];   M13 = mchn[2];
                         M22 = mchn[3];   M23 = mchn[4];
                                          M33 = mchn[5];
        double Mexp = (M11 + M22 + M33)/3.0;
        M11 -= Mexp;
        M22 -= Mexp;
        M33 -= Mexp;

        double saz2, caz2;
        saz2 = 2.0*saz*caz;
        caz2 = 2.0*caz*caz - 1.0;

        double A0, A1, A2, A4, A5;
        A0 = mult * ((2.0*M33 - M11 - M22)/6.0 );
        A1 = mult * (- (M13*caz + M23*saz));
        A2 = mult * (0.5*(M11 - M22)*caz2+ M12*saz2);
        A4 = mult * (M13*saz - M23*caz);
        A5 = mult * (-0.5*(M11 - M22)*saz2 + M12*caz2);

        srcRadi[0][0] = srcRadi[0][1] = (par_theta)? 0.0 : mult*Mexp; // EX, Z/R
        srcRadi[3][0] = srcRadi[3][1] = (par_theta)? 0.0 : A0; // DD, Z/R
        srcRadi[4][0] = srcRadi[4][1] = (par_theta)? A4 : A1; // DS, Z/R
        srcRadi[5][0] = srcRadi[5][1] = (par_theta)? 2.0*A5 : A2; // SS, Z/R
        srcRadi[0][2] = 0.0; // EX, T
        srcRadi[3][2] = 0.0; // DD, T
        srcRadi[4][2] = (par_theta)? -A1 : A4;  // DS, T
        srcRadi[5][2] = (par_theta)? -2.0*A2 : A5;  // DS, T
    }
}