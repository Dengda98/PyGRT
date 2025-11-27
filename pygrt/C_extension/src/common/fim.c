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

#include "grt/common/fim.h"
#include "grt/common/integral.h"
#include "grt/common/iostats.h"
#include "grt/common/const.h"
#include "grt/common/model.h"



real_t grt_linear_filon_integ(
    GRT_MODEL1D *mod1d, real_t k0, real_t dk0, real_t dk, real_t kmax, real_t keps,
    size_t nr, real_t *rs,
    cplx_t sum_J0[nr][GRT_SRC_M_NUM][GRT_INTEG_NUM],
    bool calc_upar,
    cplx_t sum_uiz_J0[nr][GRT_SRC_M_NUM][GRT_INTEG_NUM],
    cplx_t sum_uir_J0[nr][GRT_SRC_M_NUM][GRT_INTEG_NUM],
    FILE *fstats, GRT_KernelFunc kerfunc)
{   
    // 从0开始，存储第二部分Filon积分的结果
    cplx_t (*sum_J)[GRT_SRC_M_NUM][GRT_INTEG_NUM] = (cplx_t(*)[GRT_SRC_M_NUM][GRT_INTEG_NUM])calloc(nr, sizeof(*sum_J));
    cplx_t (*sum_uiz_J)[GRT_SRC_M_NUM][GRT_INTEG_NUM] = (calc_upar)? (cplx_t(*)[GRT_SRC_M_NUM][GRT_INTEG_NUM])calloc(nr, sizeof(*sum_uiz_J)) : NULL;
    cplx_t (*sum_uir_J)[GRT_SRC_M_NUM][GRT_INTEG_NUM] = (calc_upar)? (cplx_t(*)[GRT_SRC_M_NUM][GRT_INTEG_NUM])calloc(nr, sizeof(*sum_uir_J)) : NULL;

    cplx_t SUM[GRT_SRC_M_NUM][GRT_INTEG_NUM];

    // 不同震源不同阶数的核函数 F(k, w) 
    cplx_t QWV[GRT_SRC_M_NUM][GRT_QWV_NUM];
    cplx_t QWV_uiz[GRT_SRC_M_NUM][GRT_QWV_NUM];

    real_t k=k0; 
    size_t ik=0;
    
    // 所有震中距的k循环是否结束
    bool iendk = true;

    // 每个震中距的k循环是否结束
    bool *iendkrs = (bool *)calloc(nr, sizeof(bool)); // 自动初始化为 false
    bool iendk0 = false;

    // k循环 
    ik = 0;
    while(true){
        
        if(k > kmax && ik > 2) break;
        k += dk; 

        // 计算核函数 F(k, w)
        kerfunc(mod1d, k, QWV, calc_upar, QWV_uiz); 
        if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

        // 记录积分结果
        if(fstats!=NULL)  grt_write_stats(fstats, k, QWV);

        // 震中距rs循环
        iendk = true;
        for(size_t ir=0; ir<nr; ++ir){
            if(iendkrs[ir]) continue; // 该震中距下的波数k积分已收敛

            for(int i=0; i<GRT_SRC_M_NUM; ++i){
                for(int v=0; v<GRT_INTEG_NUM; ++v){
                    SUM[i][v] = 0.0;
                }
            }
            
            // F(k, w)*Jm(kr)k 的近似公式, sqrt(k) * F(k,w) * cos
            grt_int_Pk_filon(k, rs[ir], true, QWV, false, SUM);

            iendk0 = true;
            for(int i=0; i<GRT_SRC_M_NUM; ++i){
                int modr = GRT_SRC_M_ORDERS[i];

                for(int v=0; v<GRT_INTEG_NUM; ++v){
                    sum_J[ir][i][v] += SUM[i][v];
                    
                    // 是否提前判断达到收敛
                    if(keps <= 0.0 || (modr==0 && v!=0 && v!=2))  continue;
                    
                    iendk0 = iendk0 && (fabs(SUM[i][v])/ fabs(sum_J[ir][i][v]) <= keps);
                }
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
                grt_int_Pk_filon(k, rs[ir], true, QWV_uiz, false, SUM);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uiz也收敛
                for(int i=0; i<GRT_SRC_M_NUM; ++i){
                    for(int v=0; v<GRT_INTEG_NUM; ++v){
                        sum_uiz_J[ir][i][v] += SUM[i][v];
                    }
                }

                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                grt_int_Pk_filon(k, rs[ir], true, QWV, true, SUM);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uir也收敛
                for(int i=0; i<GRT_SRC_M_NUM; ++i){
                    for(int v=0; v<GRT_INTEG_NUM; ++v){
                        sum_uir_J[ir][i][v] += SUM[i][v];
                    }
                }
            } // END if calc_upar

            
        }  // end rs loop 
        
        ++ik;
        // 所有震中距的格林函数都已收敛
        if(iendk) break;

    } // end k loop



    // ------------------------------------------------------------------------------
    // 为累计项乘系数
    for(size_t ir=0; ir<nr; ++ir){
        real_t tmp = 2.0*(1.0 - cos(dk*rs[ir])) / (rs[ir]*rs[ir]*dk);

        for(int i=0; i<GRT_SRC_M_NUM; ++i){
            for(int v=0; v<GRT_INTEG_NUM; ++v){
                sum_J[ir][i][v] *= tmp;

                if(calc_upar){
                    sum_uiz_J[ir][i][v] *= tmp;
                    sum_uir_J[ir][i][v] *= tmp;
                }
            }
        }
    }


    // -------------------------------------------------------------------------------
    // 计算余项, [2]表示k积分的第一个点和最后一个点
    cplx_t SUM_Gc[2][GRT_SRC_M_NUM][GRT_INTEG_NUM] = {0};
    cplx_t SUM_Gs[2][GRT_SRC_M_NUM][GRT_INTEG_NUM] = {0};


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
        kerfunc(mod1d, k0N, QWV, calc_upar, QWV_uiz);
        if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN; 

        for(size_t ir=0; ir<nr; ++ir){
            // Gc
            grt_int_Pk_filon(k0N, rs[ir], true, QWV, false, SUM_Gc[iik]);
            
            // Gs
            grt_int_Pk_filon(k0N, rs[ir], false, QWV, false, SUM_Gs[iik]);

            
            real_t tmp = 1.0 / (rs[ir]*rs[ir]*dk);
            real_t tmpc = tmp * (1.0 - cos(dk*rs[ir]));
            real_t tmps = sgn * tmp * sin(dk*rs[ir]);

            for(int i=0; i<GRT_SRC_M_NUM; ++i){
                for(int v=0; v<GRT_INTEG_NUM; ++v){
                    sum_J[ir][i][v] += (- tmpc*SUM_Gc[iik][i][v] + tmps*SUM_Gs[iik][i][v] - sgn*SUM_Gs[iik][i][v]/rs[ir]);
                }
            }

            // ---------------- 位移空间导数，SUM_Gc/s数组重复利用 --------------------------
            if(calc_upar){
                // ------------------------------- ui_z -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                // Gc
                grt_int_Pk_filon(k0N, rs[ir], true, QWV_uiz, false, SUM_Gc[iik]);
                
                // Gs
                grt_int_Pk_filon(k0N, rs[ir], false, QWV_uiz, false, SUM_Gs[iik]);

                for(int i=0; i<GRT_SRC_M_NUM; ++i){
                    for(int v=0; v<GRT_INTEG_NUM; ++v){
                        sum_uiz_J[ir][i][v] += (- tmpc*SUM_Gc[iik][i][v] + tmps*SUM_Gs[iik][i][v] - sgn*SUM_Gs[iik][i][v]/rs[ir]);
                    }
                }


                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                // Gc
                grt_int_Pk_filon(k0N, rs[ir], true, QWV, true, SUM_Gc[iik]);
                
                // Gs
                grt_int_Pk_filon(k0N, rs[ir], false, QWV, true, SUM_Gs[iik]);

                for(int i=0; i<GRT_SRC_M_NUM; ++i){
                    for(int v=0; v<GRT_INTEG_NUM; ++v){
                        sum_uir_J[ir][i][v] += (- tmpc*SUM_Gc[iik][i][v] + tmps*SUM_Gs[iik][i][v] - sgn*SUM_Gs[iik][i][v]/rs[ir]);
                    }
                }
            } // END if calc_upar
          
        }  // END rs loop
    
    }  // END k 2-points loop

    // 乘上总系数 sqrt(2.0/(PI*r)) / dk0,  除dks0是在该函数外还会再乘dk0, 并将结果加到原数组中
    for(size_t ir=0; ir<nr; ++ir){
        real_t tmp = sqrt(2.0/(PI*rs[ir])) / dk0;
        for(int i=0; i<GRT_SRC_M_NUM; ++i){
            for(int v=0; v<GRT_INTEG_NUM; ++v){
                sum_J0[ir][i][v] += sum_J[ir][i][v] * tmp;

                if(calc_upar){
                    sum_uiz_J0[ir][i][v] += sum_uiz_J[ir][i][v] * tmp;
                    sum_uir_J0[ir][i][v] += sum_uir_J[ir][i][v] * tmp;
                }
            }
        }
    }


    BEFORE_RETURN:
    GRT_SAFE_FREE_PTR(sum_J);
    GRT_SAFE_FREE_PTR(sum_uiz_J);
    GRT_SAFE_FREE_PTR(sum_uir_J);


    GRT_SAFE_FREE_PTR(iendkrs);

    return k;
}

