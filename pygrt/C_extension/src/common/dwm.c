/**
 * @file   dwm.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 使用离散波数法求积分，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient 
 *               matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699
 * 
 */


#include <stdio.h> 
#include <stdlib.h>
#include <string.h>

#include "grt/common/dwm.h"
#include "grt/common/kernel.h"
#include "grt/common/integral.h"
#include "grt/common/iostats.h"
#include "grt/common/model.h"
#include "grt/common/const.h"


real_t grt_discrete_integ(
    GRT_MODEL1D *mod1d, real_t dk, real_t kmax, real_t keps,
    size_t nr, real_t *rs,
    cplxIntegGrid sum_J[nr],
    bool calc_upar,
    cplxIntegGrid sum_uiz_J[nr],
    cplxIntegGrid sum_uir_J[nr],
    FILE *fstats, GRT_KernelFunc kerfunc)
{
    cplxIntegGrid SUM = {0};

    // 不同震源不同阶数的核函数 F(k, w) 
    cplxChnlGrid QWV = {0};
    cplxChnlGrid QWV_uiz = {0};
    
    real_t k = 0.0;
    size_t ik = 0;

    // 所有震中距的k循环是否结束
    bool iendk = true;

    // 每个震中距的k循环是否结束
    bool *iendkrs = (bool *)calloc(nr, sizeof(bool)); // 自动初始化为 false
    bool iendk0 = false;

    // 波数k循环 (5.9.2)
    while(true){
        
        if(k > kmax && ik > 2)  break;
        k += dk; 

        // 计算核函数 F(k, w)
        kerfunc(mod1d, k, QWV, calc_upar, QWV_uiz); 
        if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
        
        // 记录积分核函数
        if(fstats!=NULL)  grt_write_stats(fstats, k, QWV);

        // 震中距rs循环
        iendk = true;
        for(size_t ir=0; ir<nr; ++ir){
            if(iendkrs[ir]) continue; // 该震中距下的波数k积分已收敛

            memset(SUM, 0, sizeof(cplxIntegGrid));
            
            // 计算被积函数一项 F(k,w)Jm(kr)k
            grt_int_Pk(k, rs[ir], QWV, false, SUM);
            
            iendk0 = true;

            GRT_LOOP_IntegGrid(im, v){
                int modr = GRT_SRC_M_ORDERS[im];
                sum_J[ir][im][v] += SUM[im][v];
                    
                // 是否提前判断达到收敛
                if(keps <= 0.0 || (modr==0 && v!=0 && v!=2))  continue;
                
                iendk0 = iendk0 && (fabs(SUM[im][v])/ fabs(sum_J[ir][im][v]) <= keps);
            }
            
            if(keps > 0.0){
                iendkrs[ir] = iendk0;
                iendk = iendk && iendkrs[ir];
            } else {
                iendk = iendkrs[ir] = false;
            }
            

            // ---------------- 位移空间导数，SUM数组重复利用 --------------------------
            if(calc_upar){
                // ------------------------------- ui_z -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                grt_int_Pk(k, rs[ir], QWV_uiz, false, SUM);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uiz也收敛
                GRT_LOOP_IntegGrid(im, v){
                    sum_uiz_J[ir][im][v] += SUM[im][v];
                }

                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                grt_int_Pk(k, rs[ir], QWV, true, SUM);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uiz也收敛
                GRT_LOOP_IntegGrid(im, v){
                    sum_uir_J[ir][im][v] += SUM[im][v];
                }
            } // END if calc_upar

        } // END rs loop

        ++ik;

        // 所有震中距的格林函数都已收敛
        if(iendk) break;

    } // END k loop

    BEFORE_RETURN:
    GRT_SAFE_FREE_PTR(iendkrs);

    return k;

}

