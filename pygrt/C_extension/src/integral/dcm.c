/**
 * @file   dcm.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-12
 * 
 *     DCM 的校正系数
 *                   
 */

#include "grt/common/const.h"
#include "grt/integral/k_integ.h"

/** DCM 的校正项，其中积分定义与 grt_int_Pk() 函数保持一致 */
static void _correct_Wm(real_t coefs[GRT_MORDER_MAX+1], bool keep_nearfield, cplxChnlGrid QWV_kmax, cplxIntegGrid SUM)
{
    for(int i = 0; i < GRT_SRC_M_NUM; ++i){
        int modr = GRT_SRC_M_ORDERS[i];  // 对应m阶数
        if(modr == 0){
            SUM[i][0] += - QWV_kmax[i][0] * coefs[1];   // - q0*J1*k
            SUM[i][2] +=   QWV_kmax[i][1] * coefs[0];   //   w0*J0*k
        }
        else{
            // 近场项这里沿用coefs是公式恰好巧合
            SUM[i][0]  +=   QWV_kmax[i][0] * coefs[modr-1];         // qm*Jm-1*k
            if(keep_nearfield) {
                SUM[i][1]  += - (QWV_kmax[i][0] + QWV_kmax[i][2]) * coefs[modr];    // - m*(qm+vm)*Jm*k/kr
            }
            SUM[i][2]  +=   QWV_kmax[i][1] * coefs[modr];           // wm*Jm*k
            SUM[i][3]  += - QWV_kmax[i][2] * coefs[modr-1];         // -vm*Jm-1*k
        }
    }
}

void grt_dcm_correction(size_t nr, real_t *rs, K_INTEG *Kint, bool keep_nearfield)
{
    for(size_t ir = 0; ir < nr; ++ir)    {
        real_t r = rs[ir];
        real_t cc = 1.0 / r;
        
        // 1 / r^2
        cc = cc * cc;
        real_t coefs[] = {0.0, cc, 2.0*cc};

        // 1 / r^3
        cc = cc * cc;
        real_t coefs_r[] = {0.0, - 2.0*cc, - 4.0*cc};

        _correct_Wm(coefs, keep_nearfield, Kint->QWV_kmax, Kint->sumJ[ir]);
        if(Kint->calc_upar){
            _correct_Wm(coefs,   keep_nearfield, Kint->QWVz_kmax, Kint->sumJz[ir]);
            _correct_Wm(coefs_r, keep_nearfield, Kint->QWV_kmax,  Kint->sumJr[ir]);
        }
    }
}