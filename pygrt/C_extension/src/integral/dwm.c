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

#include "grt/integral/dwm.h"
#include "grt/integral/kernel.h"
#include "grt/integral/k_integ.h"
#include "grt/integral/iostats.h"
#include "grt/common/model.h"
#include "grt/common/const.h"


real_t grt_discrete_integ(
    GRT_MODEL1D *mod1d, real_t dk, real_t kmax, real_t keps,
    size_t nr, real_t *rs, K_INTEG *K, FILE *fstats, GRT_KernelFunc kerfunc)
{
    if(kmax == 0.0)  return 0.0;
    
    real_t k = 0.0;

    // 所有震中距的k循环是否结束
    bool iendk = true;

    // 每个震中距的k循环是否结束
    bool *iendkrs = (bool *)calloc(nr, sizeof(bool)); // 自动初始化为 false
    bool iendk0 = false;

    size_t nk = floor(kmax / dk) + 1L;

    // 波数k循环 (5.9.2)
    for(size_t ik = 0; ik < nk; ++ik){
        
        k += dk; 

        // 计算核函数 F(k, w)
        kerfunc(mod1d, k, K->QWV, K->calc_upar, K->QWVz); 
        if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

        if(K->applyDCM){
            GRT_LOOP_ChnlGrid(im, c){
                K->QWV[im][c] -= K->QWV_kmax[im][c];
                if(K->calc_upar) K->QWVz[im][c] -= K->QWVz_kmax[im][c];
            }
        }

        // 记录积分核函数
        if(fstats!=NULL)  grt_write_stats(fstats, k, (K->calc_upar)? K->QWVz : K->QWV);

        // 震中距rs循环
        iendk = true;
        for(size_t ir=0; ir<nr; ++ir){
            if(iendkrs[ir]) continue; // 该震中距下的波数k积分已收敛

            memset(K->SUM, 0, sizeof(cplxIntegGrid));
            
            // 计算被积函数一项 F(k,w)Jm(kr)k
            grt_int_Pk(k, rs[ir], K->QWV, false, K->SUM);
            
            iendk0 = true;

            GRT_LOOP_IntegGrid(im, v){
                int modr = GRT_SRC_M_ORDERS[im];
                K->sumJ[ir][im][v] += K->SUM[im][v];
                    
                // 是否提前判断达到收敛
                if(keps <= 0.0 || (modr==0 && v!=0 && v!=2))  continue;
                
                iendk0 = iendk0 && (fabs(K->SUM[im][v])/ fabs(K->sumJ[ir][im][v]) <= keps);
            }
            
            if(keps > 0.0){
                iendkrs[ir] = iendk0;
                iendk = iendk && iendkrs[ir];
            } else {
                iendk = iendkrs[ir] = false;
            }
            

            // ---------------- 位移空间导数，SUM数组重复利用 --------------------------
            if(K->calc_upar){
                // ------------------------------- ui_z -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                grt_int_Pk(k, rs[ir], K->QWVz, false, K->SUM);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uiz也收敛
                GRT_LOOP_IntegGrid(im, v){
                    K->sumJz[ir][im][v] += K->SUM[im][v];
                }

                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                grt_int_Pk(k, rs[ir], K->QWV, true, K->SUM);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uiz也收敛
                GRT_LOOP_IntegGrid(im, v){
                    K->sumJr[ir][im][v] += K->SUM[im][v];
                }
            } // END if calc_upar

        } // END rs loop

        // 所有震中距的格林函数都已收敛
        if(iendk) break;

    } // END k loop

    BEFORE_RETURN:
    GRT_SAFE_FREE_PTR(iendkrs);

    return k;

}

