/**
 * @file   fim.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是基于线性插值的Filon积分，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. 纪晨, 姚振兴. 1995. 区域地震范围的宽频带理论地震图算法研究. 地球物理学报. 38(4)
 * 
 */

#include <stdio.h> 
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "grt/integral/fim.h"
#include "grt/integral/k_integ.h"
#include "grt/integral/iostats.h"
#include "grt/common/const.h"
#include "grt/common/model.h"



real_t grt_linear_filon_integ(
    GRT_MODEL1D *mod1d, real_t k0, real_t dk0, real_t dk, real_t kmax, real_t keps, bool applyDCM,
    size_t nr, real_t *rs, K_INTEG *K0, FILE *fstats, GRT_KernelFunc kerfunc)
{   
    if(k0 + dk0 >= kmax)  return k0;
    
    // 从0开始，存储第二部分Filon积分的结果
    K_INTEG *K = grt_init_K_INTEG(K0->calc_upar, nr);

    // 不同震源不同阶数的核函数 F(k, w) 
    cplxChnlGrid QWV = {0};
    cplxChnlGrid QWVz = {0};

    cplxChnlGrid QWV_kmax = {0};
    cplxChnlGrid QWVz_kmax = {0};

    real_t k=k0; 
    
    // 所有震中距的k循环是否结束
    bool iendk = true;

    // 每个震中距的k循环是否结束
    bool *iendkrs = (bool *)calloc(nr, sizeof(bool)); // 自动初始化为 false
    bool iendk0 = false;

    size_t nk = floor((kmax - k0) / dk) + 1L;

    if(applyDCM){
        kerfunc(mod1d, nk*dk, QWV_kmax, K->calc_upar, QWVz_kmax); 
    }

    // k循环 
    for(size_t ik = 0; ik < nk; ++ik){
        
        k += dk; 

        // 计算核函数 F(k, w)
        kerfunc(mod1d, k, QWV, K->calc_upar, QWVz); 
        if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

        if(applyDCM){
            GRT_LOOP_ChnlGrid(im, c){
                QWV[im][c] -= QWV_kmax[im][c];
                if(K->calc_upar) QWVz[im][c] -= QWVz_kmax[im][c];
            }
        }

        // 记录积分结果
        if(fstats!=NULL)  grt_write_stats(fstats, k, QWV);

        // 震中距rs循环
        iendk = true;
        for(size_t ir=0; ir<nr; ++ir){
            if(iendkrs[ir]) continue; // 该震中距下的波数k积分已收敛

            memset(K->SUM, 0, sizeof(cplxIntegGrid));
            
            // F(k, w)*Jm(kr)k 的近似公式, sqrt(k) * F(k,w) * cos
            grt_int_Pk_filon(k, rs[ir], true, QWV, false, K->SUM);

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
                grt_int_Pk_filon(k, rs[ir], true, QWVz, false, K->SUM);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uiz也收敛
                GRT_LOOP_IntegGrid(im, v){
                    K->sumJz[ir][im][v] += K->SUM[im][v];
                }

                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                grt_int_Pk_filon(k, rs[ir], true, QWV, true, K->SUM);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uir也收敛
                GRT_LOOP_IntegGrid(im, v){
                    K->sumJr[ir][im][v] += K->SUM[im][v];
                }
            } // END if calc_upar

            
        }  // end rs loop 
        
        // 所有震中距的格林函数都已收敛
        if(iendk) break;

    } // end k loop



    // ------------------------------------------------------------------------------
    // 为累计项乘系数
    for(size_t ir=0; ir<nr; ++ir){
        real_t tmp = 2.0*(1.0 - cos(dk*rs[ir])) / (rs[ir]*rs[ir]*dk);

        GRT_LOOP_IntegGrid(im, v){
            K->sumJ[ir][im][v] *= tmp;

            if(K->calc_upar){
                K->sumJz[ir][im][v] *= tmp;
                K->sumJr[ir][im][v] *= tmp;
            }
        }
    }


    // -------------------------------------------------------------------------------
    // 计算余项, [2]表示k积分的第一个点和最后一个点
    cplxIntegGrid SUM_Gc[2] = {0};
    cplxIntegGrid SUM_Gs[2] = {0};


    // 计算来自第一个点和最后一个点的余项
    for(int iik=0; iik<2; ++iik){ 
        real_t k0N;
        int sgn;
        if(0==iik)       {k0N = k0+dk; sgn =  1.0;}
        else if(1==iik)  {k0N = k;     sgn = -1.0;}
        else {
            fprintf(stderr, "Filon error.\n");
            exit(EXIT_FAILURE);
        }

        // 计算核函数 F(k, w)
        kerfunc(mod1d, k0N, QWV, K->calc_upar, QWVz);
        if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN; 

        if(applyDCM){
            GRT_LOOP_ChnlGrid(im, c){
                QWV[im][c] -= QWV_kmax[im][c];
                if(K->calc_upar) QWVz[im][c] -= QWVz_kmax[im][c];
            }
        }

        for(size_t ir=0; ir<nr; ++ir){
            // Gc
            grt_int_Pk_filon(k0N, rs[ir], true, QWV, false, SUM_Gc[iik]);
            
            // Gs
            grt_int_Pk_filon(k0N, rs[ir], false, QWV, false, SUM_Gs[iik]);

            
            real_t tmp = 1.0 / (rs[ir]*rs[ir]*dk);
            real_t tmpc = tmp * (1.0 - cos(dk*rs[ir]));
            real_t tmps = sgn * tmp * sin(dk*rs[ir]);

            GRT_LOOP_IntegGrid(im, v){
                K->sumJ[ir][im][v] += (- tmpc*SUM_Gc[iik][im][v] + tmps*SUM_Gs[iik][im][v] - sgn*SUM_Gs[iik][im][v]/rs[ir]);
            }

            // ---------------- 位移空间导数，SUM_Gc/s数组重复利用 --------------------------
            if(K->calc_upar){
                // ------------------------------- ui_z -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                // Gc
                grt_int_Pk_filon(k0N, rs[ir], true, QWVz, false, SUM_Gc[iik]);
                
                // Gs
                grt_int_Pk_filon(k0N, rs[ir], false, QWVz, false, SUM_Gs[iik]);

                GRT_LOOP_IntegGrid(im, v){
                    K->sumJz[ir][im][v] += (- tmpc*SUM_Gc[iik][im][v] + tmps*SUM_Gs[iik][im][v] - sgn*SUM_Gs[iik][im][v]/rs[ir]);
                }

                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                // Gc
                grt_int_Pk_filon(k0N, rs[ir], true, QWV, true, SUM_Gc[iik]);
                
                // Gs
                grt_int_Pk_filon(k0N, rs[ir], false, QWV, true, SUM_Gs[iik]);

                GRT_LOOP_IntegGrid(im, v){
                    K->sumJr[ir][im][v] += (- tmpc*SUM_Gc[iik][im][v] + tmps*SUM_Gs[iik][im][v] - sgn*SUM_Gs[iik][im][v]/rs[ir]);
                }
            } // END if calc_upar
          
        }  // END rs loop
    
    }  // END k 2-points loop

    // 乘上总系数 sqrt(2.0/(PI*r)) / dk0,  除dks0是在该函数外还会再乘dk0, 并将结果加到原数组中
    for(size_t ir=0; ir<nr; ++ir){
        real_t tmp = sqrt(2.0/(PI*rs[ir])) / dk0;

        GRT_LOOP_IntegGrid(im, v){
            K0->sumJ[ir][im][v] += K->sumJ[ir][im][v] * tmp;

            if(K->calc_upar){
                K0->sumJz[ir][im][v] += K->sumJz[ir][im][v] * tmp;
                K0->sumJr[ir][im][v] += K->sumJr[ir][im][v] * tmp;
            }
        }
    }


    BEFORE_RETURN:
    grt_free_K_INTEG(K);

    GRT_SAFE_FREE_PTR(iendkrs);

    return k;
}

