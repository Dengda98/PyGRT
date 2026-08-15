/**
 * @file   dcm.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-12
 * 
 *     DCM 的校正系数
 *                   
 */

#include "grt/common/const.h"
#include "grt/integral/dcm.h"
#include "grt/integral/k_integ.h"

/** DCM 的校正项，其中积分定义与 grt_int_Pk() 函数保持一致 */
static void _correct_Wm(
    size_t i, real_t coefs[GRT_MORDER_MAX+1], real_t coefs_near[GRT_MORDER_MAX+1], 
    bool keep_nearfield, cplxChnlGrid QWV_kmax, real_t coef, cplxIntegGrid SUM)
{
    int modr = GRT_SRC_M_ORDERS[i];  // 对应m阶数
    if(modr == 0){
        SUM[i][0] += - coef * QWV_kmax[i][0] * coefs[1];   // - q0*J1*k
        SUM[i][2] +=   coef * QWV_kmax[i][1] * coefs[0];   //   w0*J0*k
    }
    else{
        SUM[i][0]  +=   coef * QWV_kmax[i][0] * coefs[modr-1];         // qm*Jm-1*k
        if(keep_nearfield) {
            SUM[i][1]  += - modr * coef * (QWV_kmax[i][0] + QWV_kmax[i][2]) * coefs_near[modr];    // - m*(qm+vm)*Jm*k/kr
        }
        SUM[i][2]  +=   coef * QWV_kmax[i][1] * coefs[modr];           // wm*Jm*k
        SUM[i][3]  += - coef * QWV_kmax[i][2] * coefs[modr-1];         // -vm*Jm-1*k
    }
}

void grt_dcm_correction(size_t nr, real_t *rs, real_t dk, real_t kcut, K_INTEG *Kint, bool keep_nearfield)
{
    for(size_t ir = 0; ir < nr; ++ir){
        real_t r = rs[ir];
        // 修正系数含 1/r；r=0 跳过（近场极限已在 k_integ 的 Bessel 极限中处理）
        if(GRT_IS_ZERO(r)) continue;

        real_t c = 1.0 / r;

        // m = 0 时没有近场项，数组下标 0 是未使用的哨兵值
        
        // =====================================================
        real_t c2 = c * c;  // 1 / r^2
        // m=0 时 J_0(kr) 在 k=0 时为 1，这导致数值积分需要考虑 k=0
        real_t force_coefs[]      = {c - 0.5*dk,   c,     c};            // \int_0^\infty J_m(kr) dk
        // m=1 时 J_1(kr)/k 在 k=0 的极限为 r/2；DWM 从 dk 开始，
        // 因而右端点矩形和的 Euler--Maclaurin 修正为 -dk/4。
        real_t force_coefs_near[] = {0.0, c - 0.25*dk, 0.5*c};  // 1/r \int_0^\infty J_m(kr) 1/k dk

        real_t force_coefs_z[]      = {0.0,   c2, 2.0*c2};       // \int_0^\infty J_m(kr) k dk
        real_t force_coefs_z_near[] = {0.0,   c2,     c2};       // 1/r \int_0^\infty J_m(kr) dk

        real_t force_coefs_r[]      = {-c2,    -c2,     -c2};    // \int_0^\infty J_m^'(kr) k dk
        real_t force_coefs_r_near[] = {0.0,    -c2, -0.5*c2};    // \int_0^\infty d/dr[1/r J_m(kr)] 1/k dk = 1/r \int_0^\infty J_m^'(kr) dk - 1/r^2 \int_0^\infty J_m(kr) 1/k dk
        
        // =====================================================
        real_t c3 = c2 * c;    // 1 / r^3
        real_t couple_coefs[]      = {0.0, c2, 2.0*c2};          // \int_0^\infty J_m(kr) k dk
        real_t couple_coefs_near[] = {0.0, c2,     c2};          // 1/r \int_0^\infty J_m(kr) dk

        real_t couple_coefs_z[]      = {-1.0*c3, 0.0, 3.0*c3};   // \int_0^\infty J_m(kr) k^2 dk
        real_t couple_coefs_z_near[] = {0.0,      c3, 2.0*c3};   // 1/r \int_0^\infty J_m(kr) k dk
        
        real_t couple_coefs_r[]      = {0.0, -2.0*c3, -4.0*c3};  // \int_0^\infty J_m^'(kr) k^2 dk
        real_t couple_coefs_r_near[] = {0.0, -2.0*c3, -2.0*c3};  // \int_0^\infty d/dr[1/r J_m(kr)] dk = 1/r \int_0^\infty J_m^'(kr) k dk - 1/r^2 \int_0^\infty J_m(kr) dk

        // 这里 dk 只是先把数值积分的 dk 补上
        for(int i = 0; i < GRT_SRC_M_NUM; ++i){
            if(GRT_SRC_M_INDEX_IS_FORCE(i)){
                _correct_Wm(i, force_coefs, force_coefs_near, keep_nearfield, Kint->QWV_kmax, kcut / dk, Kint->sumJ[ir]);
                if(Kint->calc_upar){
                    _correct_Wm(i, force_coefs_z, force_coefs_z_near, keep_nearfield, Kint->QWVz_kmax, 1.0 / dk, Kint->sumJz[ir]);
                    _correct_Wm(i, force_coefs_r, force_coefs_r_near, keep_nearfield, Kint->QWV_kmax,  kcut / dk, Kint->sumJr[ir]);
                }
            } else {
                _correct_Wm(i, couple_coefs, couple_coefs_near, keep_nearfield, Kint->QWV_kmax, 1.0 / dk, Kint->sumJ[ir]);
                if(Kint->calc_upar){
                    _correct_Wm(i, couple_coefs_z, couple_coefs_z_near, keep_nearfield, Kint->QWVz_kmax, 1.0 / kcut / dk, Kint->sumJz[ir]);
                    _correct_Wm(i, couple_coefs_r, couple_coefs_r_near, keep_nearfield, Kint->QWV_kmax,  1.0 / dk, Kint->sumJr[ir]);
                }
            }
        }
        
    }
}
