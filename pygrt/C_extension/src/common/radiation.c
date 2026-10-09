/**
 * @file   radiation.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-06
 * 
 *    计算不同震源的辐射因子
 * 
 */

#include <stdbool.h>

#include "grt/common/radiation.h"
#include "grt/common/const.h"
#include "grt/common/checkerror.h"

void grt_set_source_radiation(
    realChnlGrid srcRadi, const GRT_SYN_TYPE source_type, const bool par_theta,
    const real_t scale, const real_t coef, const real_t VpVs_ratio, const real_t azrad, const real_t mchn[GRT_MECHANISM_NUM]
){
    real_t mult;
    if(source_type == GRT_SYN_SF){
        mult = 1e-15*scale*coef;
    } else {
        mult = 1e-20*scale*coef;
    }

    real_t saz, caz;
    saz = sin(azrad);
    caz = cos(azrad);

    if(source_type == GRT_SYN_EX){
        srcRadi[0][0] = srcRadi[0][1] = (par_theta)? 0.0 : mult; // Z/R
        srcRadi[0][2] = 0.0; // T
    }  
    else if(source_type == GRT_SYN_SF){
        real_t A0, A1, A4;
        real_t fn, fe, fz;
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
    else if(source_type == GRT_SYN_DC){
        real_t strike, dip, rake;
        strike=mchn[0];   dip=mchn[1];   rake=mchn[2];
        // 公式(4.8.35)
        real_t stkrad = strike*DEG1;
        real_t diprad = dip*DEG1;
        real_t rakrad = rake*DEG1;
        real_t therad = azrad - stkrad;
        real_t srak, crak, sdip, cdip, sdip2, cdip2, sthe, cthe, sthe2, cthe2;
        srak = sin(rakrad);     crak = cos(rakrad);
        sdip = sin(diprad);     cdip = cos(diprad);
        sdip2 = 2.0*sdip*cdip;  cdip2 = 2.0*cdip*cdip - 1.0;
        sthe = sin(therad);     cthe = cos(therad);
        sthe2 = 2.0*sthe*cthe;  cthe2 = 2.0*cthe*cthe - 1.0;

        real_t A0, A1, A2, A4, A5;
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
    else if(source_type == GRT_SYN_TS){
        real_t strike, dip;
        strike=mchn[0];   dip=mchn[1];
        // 公式(4.8.40)但修改了各向同性的量
        real_t stkrad = strike*DEG1;
        real_t diprad = dip*DEG1;
        real_t therad = azrad - stkrad;
        real_t sdip, cdip, sdip2, sthe, cthe, sthe2, cthe2;
        sdip = sin(diprad);     cdip = cos(diprad);
        sdip2 = 2.0*sdip*cdip;
        sthe = sin(therad);     cthe = cos(therad);
        sthe2 = 2.0*sthe*cthe;  cthe2 = 2.0*cthe*cthe - 1.0;

        real_t A0, A1, A2, A4, A5;
        A0 = mult * (2.0/3.0 - sdip*sdip);
        A1 = mult * sdip2*sthe;
        A2 = mult * ( - sdip*sdip*cthe2);
        A4 = mult * sdip2*cthe;
        A5 = mult * sdip*sdip*sthe2;

        srcRadi[0][0] = srcRadi[0][1] = (par_theta)? 0.0 : mult * (GRT_SQUARE(VpVs_ratio) - 4.0/3.0); // EX, Z/R
        srcRadi[3][0] = srcRadi[3][1] = (par_theta)? 0.0 : A0; // DD, Z/R
        srcRadi[4][0] = srcRadi[4][1] = (par_theta)? A4 : A1; // DS, Z/R
        srcRadi[5][0] = srcRadi[5][1] = (par_theta)? 2.0*A5 : A2; // SS, Z/R
        srcRadi[0][2] = 0.0; // EX, T
        srcRadi[3][2] = 0.0; // DD, T
        srcRadi[4][2] = (par_theta)? -A1 : A4;  // DS, T
        srcRadi[5][2] = (par_theta)? -2.0*A2 : A5;  // DS, T
    }
    else if(source_type == GRT_SYN_MT){
        // 公式(4.9.7)但修改了各向同性的量
        real_t M11, M12, M13, M22, M23, M33;
        M11 = mchn[0];   M12 = mchn[1];   M13 = mchn[2];
                         M22 = mchn[3];   M23 = mchn[4];
                                          M33 = mchn[5];
        real_t Mexp = (M11 + M22 + M33)/3.0;
        M11 -= Mexp;
        M22 -= Mexp;
        M33 -= Mexp;

        real_t saz2, caz2;
        saz2 = 2.0*saz*caz;
        caz2 = 2.0*caz*caz - 1.0;

        real_t A0, A1, A2, A4, A5;
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
    else {
        GRTRaiseError("Unsupported source type.");
    }
}


void grt_source_moment_tensor(GRT_SYN_TYPE type, real_t scale, real_t VpVs_ratio, const real_t mchn[GRT_MECHANISM_NUM], real_t tensor[6])
{
    // 在方位角为零时，从 EX、DD、DS、SS 辐射系数恢复六个 NED 张量分量
    realChnlGrid rad = {{0}};
    // 用 coef 抵消辐射系数中的 1e-20 量纲换算
    grt_set_source_radiation(rad, type, false, scale, 1e20, VpVs_ratio, 0, mchn);

    // 以下关系式从上面的 “else if(source_type == GRT_SYN_MT)” 中反解得到
    real_t e = rad[GRT_SRC_M_EX_INDEX][0], dd = rad[GRT_SRC_M_DD_INDEX][0], ss = rad[GRT_SRC_M_SS_INDEX][0];
    tensor[0] = e - dd + ss;
    tensor[3] = e - dd - ss;
    tensor[5] = e + 2 * dd;
    tensor[1] = rad[GRT_SRC_M_SS_INDEX][2];
    tensor[2] = -rad[GRT_SRC_M_DS_INDEX][0];
    tensor[4] = -rad[GRT_SRC_M_DS_INDEX][2];
}


bool grt_source_has_component(GRT_SYN_TYPE source_type, int im)
{
    // 按震源类型筛选所需基本源型，供格林函数按需加载
    if(source_type == GRT_SYN_EX) {
        return im == GRT_SRC_M_EX_INDEX;
    } else if(source_type == GRT_SYN_SF) {
        return im == GRT_SRC_M_VF_INDEX || im == GRT_SRC_M_HF_INDEX;
    } else if(source_type == GRT_SYN_DC) {
        return im >= GRT_SRC_M_DD_INDEX;
    } else if(source_type == GRT_SYN_TS || source_type == GRT_SYN_MT) {
        return im >= GRT_SRC_M_DD_INDEX || im == GRT_SRC_M_EX_INDEX;
    }
    return false;
}

void grt_check_source_medium(GRT_SYN_TYPE source_type, bool with_mu, real_t vs)
{
    if(vs == 0 && (source_type != GRT_SYN_EX || with_mu)) {
        GRTRaiseError("Liquid source media support only explosion sources without -Su.");
    }
}
