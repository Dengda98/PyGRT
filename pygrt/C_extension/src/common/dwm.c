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

#include "common/dwm.h"
#include "common/kernel.h"
#include "common/integral.h"
#include "common/iostats.h"
#include "common/model.h"
#include "common/const.h"


MYREAL discrete_integ(
    const MODEL1D *mod1d, MYREAL dk, MYREAL kmax, MYREAL keps, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_J[nr][SRC_M_NUM][INTEG_NUM],
    bool calc_upar,
    MYCOMPLEX sum_uiz_J[nr][SRC_M_NUM][INTEG_NUM],
    MYCOMPLEX sum_uir_J[nr][SRC_M_NUM][INTEG_NUM],
    FILE *fstats, KernelFunc kerfunc, MYINT *stats)
{
    MYCOMPLEX SUM[SRC_M_NUM][INTEG_NUM];

    // 不同震源不同阶数的核函数 F(k, w) 
    MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM];
    MYCOMPLEX QWV_uiz[SRC_M_NUM][QWV_NUM];
    
    MYREAL k = 0.0;
    MYINT ik = 0;

    // 所有震中距的k循环是否结束
    bool iendk = true;

    // 每个震中距的k循环是否结束
    bool *iendkrs = (bool *)malloc(nr * sizeof(bool));
    bool iendk0 = false;
    for(MYINT ir=0; ir<nr; ++ir) iendkrs[ir] = false;
    

    // 波数k循环 (5.9.2)
    while(true){
        
        if(k > kmax && ik > 2)  break;
        k += dk; 

        // printf("w=%15.5e, ik=%d\n", creal(omega), ik);
        // 计算核函数 F(k, w)
        kerfunc(mod1d, omega, k, QWV, calc_upar, QWV_uiz, stats); 
        if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN;
        
        // 记录积分核函数
        if(fstats!=NULL)  write_stats(fstats, k, QWV);

        // 震中距rs循环
        iendk = true;
        for(MYINT ir=0; ir<nr; ++ir){
            if(iendkrs[ir]) continue; // 该震中距下的波数k积分已收敛

            for(MYINT i=0; i<SRC_M_NUM; ++i){
                for(MYINT v=0; v<INTEG_NUM; ++v){
                    SUM[i][v] = CZERO;
                }
            }
            
            // 计算被积函数一项 F(k,w)Jm(kr)k
            int_Pk(k, rs[ir], QWV, false, SUM);
            
            iendk0 = true;
            for(MYINT i=0; i<SRC_M_NUM; ++i){
                MYINT modr = SRC_M_ORDERS[i];

                for(MYINT v=0; v<INTEG_NUM; ++v){
                    sum_J[ir][i][v] += SUM[i][v];
                    
                    // 是否提前判断达到收敛
                    if(keps < RZERO || (modr==0 && v!=0 && v!=2))  continue;
                    
                    iendk0 = iendk0 && (fabs(SUM[i][v])/ fabs(sum_J[ir][i][v]) <= keps);
                }
            }
            
            if(keps > RZERO){
                iendkrs[ir] = iendk0;
                iendk = iendk && iendkrs[ir];
            } else {
                iendk = iendkrs[ir] = false;
            }
            

            // ---------------- 位移空间导数，SUM数组重复利用 --------------------------
            if(calc_upar){
                // ------------------------------- ui_z -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                int_Pk(k, rs[ir], QWV_uiz, false, SUM);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uiz也收敛
                for(MYINT i=0; i<SRC_M_NUM; ++i){
                    for(MYINT v=0; v<INTEG_NUM; ++v){
                        sum_uiz_J[ir][i][v] += SUM[i][v];
                    }
                }

                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                int_Pk(k, rs[ir], QWV, true, SUM);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uiz也收敛
                for(MYINT i=0; i<SRC_M_NUM; ++i){
                    for(MYINT v=0; v<INTEG_NUM; ++v){
                        sum_uir_J[ir][i][v] += SUM[i][v];
                    }
                }
            } // END if calc_upar

        } // END rs loop

        ++ik;

        // 所有震中距的格林函数都已收敛
        if(iendk) break;

    } // END k loop

    // printf("w=%15.5e, ik=%d\n", creal(omega), ik);

    BEFORE_RETURN:
    free(iendkrs);

    return k;

}

