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

#include "common/fim.h"
#include "common/integral.h"
#include "common/iostats.h"
#include "common/const.h"
#include "common/model.h"



MYREAL linear_filon_integ(
    const MODEL1D *mod1d, MYREAL k0, MYREAL dk0, MYREAL dk, MYREAL kmax, MYREAL keps, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_J0[nr][SRC_M_NUM][INTEG_NUM],
    bool calc_upar,
    MYCOMPLEX sum_uiz_J0[nr][SRC_M_NUM][INTEG_NUM],
    MYCOMPLEX sum_uir_J0[nr][SRC_M_NUM][INTEG_NUM],
    FILE *fstats, KernelFunc kerfunc, MYINT *stats)
{   
    // 从0开始，存储第二部分Filon积分的结果
    MYCOMPLEX (*sum_J)[SRC_M_NUM][INTEG_NUM] = (MYCOMPLEX(*)[SRC_M_NUM][INTEG_NUM])calloc(nr, sizeof(*sum_J));
    MYCOMPLEX (*sum_uiz_J)[SRC_M_NUM][INTEG_NUM] = (calc_upar)? (MYCOMPLEX(*)[SRC_M_NUM][INTEG_NUM])calloc(nr, sizeof(*sum_uiz_J)) : NULL;
    MYCOMPLEX (*sum_uir_J)[SRC_M_NUM][INTEG_NUM] = (calc_upar)? (MYCOMPLEX(*)[SRC_M_NUM][INTEG_NUM])calloc(nr, sizeof(*sum_uir_J)) : NULL;

    MYCOMPLEX SUM[SRC_M_NUM][INTEG_NUM];

    // 不同震源不同阶数的核函数 F(k, w) 
    MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM];
    MYCOMPLEX QWV_uiz[SRC_M_NUM][QWV_NUM];

    MYREAL k=k0; 
    MYINT ik=0;
    
    bool iendk, iendk0;

    // 每个震中距的k循环是否结束
    bool *iendkrs = (bool *)malloc(nr * sizeof(bool));
    for(MYINT ir=0; ir<nr; ++ir) iendkrs[ir] = false;

    // k循环 
    ik = 0;
    while(true){
        
        if(k > kmax && ik > 2) break;
        k += dk; 

        // 计算核函数 F(k, w)
        kerfunc(mod1d, omega, k, QWV, calc_upar, QWV_uiz, stats); 
        if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN;

        // 记录积分结果
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
            
            // F(k, w)*Jm(kr)k 的近似公式, sqrt(k) * F(k,w) * cos
            int_Pk_filon(k, rs[ir], true, QWV, false, SUM);

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
                int_Pk_filon(k, rs[ir], true, QWV_uiz, false, SUM);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uiz也收敛
                for(MYINT i=0; i<SRC_M_NUM; ++i){
                    for(MYINT v=0; v<INTEG_NUM; ++v){
                        sum_uiz_J[ir][i][v] += SUM[i][v];
                    }
                }

                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                int_Pk_filon(k, rs[ir], true, QWV, true, SUM);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uir也收敛
                for(MYINT i=0; i<SRC_M_NUM; ++i){
                    for(MYINT v=0; v<INTEG_NUM; ++v){
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
    for(MYINT ir=0; ir<nr; ++ir){
        MYREAL tmp = RTWO*(RONE - cos(dk*rs[ir])) / (rs[ir]*rs[ir]*dk);

        for(MYINT i=0; i<SRC_M_NUM; ++i){
            for(MYINT v=0; v<INTEG_NUM; ++v){
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
    MYCOMPLEX SUM_Gc[2][SRC_M_NUM][INTEG_NUM] = {0};
    MYCOMPLEX SUM_Gs[2][SRC_M_NUM][INTEG_NUM] = {0};


    // 计算来自第一个点和最后一个点的余项
    for(MYINT iik=0; iik<2; ++iik){ 
        MYREAL k0N;
        MYINT sgn;
        if(0==iik)       {k0N = k0+dk; sgn =  RONE;}
        else if(1==iik)  {k0N = k;     sgn = -RONE;}
        else {
            fprintf(stderr, "Filon error.\n");
            exit(EXIT_FAILURE);
        }

        // 计算核函数 F(k, w)
        kerfunc(mod1d, omega, k0N, QWV, calc_upar, QWV_uiz, stats);
        if(*stats==INVERSE_FAILURE)  goto BEFORE_RETURN; 

        for(MYINT ir=0; ir<nr; ++ir){
            // Gc
            int_Pk_filon(k0N, rs[ir], true, QWV, false, SUM_Gc[iik]);
            
            // Gs
            int_Pk_filon(k0N, rs[ir], false, QWV, false, SUM_Gs[iik]);

            
            MYREAL tmp = RONE / (rs[ir]*rs[ir]*dk);
            MYREAL tmpc = tmp * (RONE - cos(dk*rs[ir]));
            MYREAL tmps = sgn * tmp * sin(dk*rs[ir]);

            for(MYINT i=0; i<SRC_M_NUM; ++i){
                for(MYINT v=0; v<INTEG_NUM; ++v){
                    sum_J[ir][i][v] += (- tmpc*SUM_Gc[iik][i][v] + tmps*SUM_Gs[iik][i][v] - sgn*SUM_Gs[iik][i][v]/rs[ir]);
                }
            }

            // ---------------- 位移空间导数，SUM_Gc/s数组重复利用 --------------------------
            if(calc_upar){
                // ------------------------------- ui_z -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                // Gc
                int_Pk_filon(k0N, rs[ir], true, QWV_uiz, false, SUM_Gc[iik]);
                
                // Gs
                int_Pk_filon(k0N, rs[ir], false, QWV_uiz, false, SUM_Gs[iik]);

                for(MYINT i=0; i<SRC_M_NUM; ++i){
                    for(MYINT v=0; v<INTEG_NUM; ++v){
                        sum_uiz_J[ir][i][v] += (- tmpc*SUM_Gc[iik][i][v] + tmps*SUM_Gs[iik][i][v] - sgn*SUM_Gs[iik][i][v]/rs[ir]);
                    }
                }


                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                // Gc
                int_Pk_filon(k0N, rs[ir], true, QWV, true, SUM_Gc[iik]);
                
                // Gs
                int_Pk_filon(k0N, rs[ir], false, QWV, true, SUM_Gs[iik]);

                for(MYINT i=0; i<SRC_M_NUM; ++i){
                    for(MYINT v=0; v<INTEG_NUM; ++v){
                        sum_uir_J[ir][i][v] += (- tmpc*SUM_Gc[iik][i][v] + tmps*SUM_Gs[iik][i][v] - sgn*SUM_Gs[iik][i][v]/rs[ir]);
                    }
                }
            } // END if calc_upar
          
        }  // END rs loop
    
    }  // END k 2-points loop

    // 乘上总系数 sqrt(RTWO/(PI*r)) / dk0,  除dks0是在该函数外还会再乘dk0, 并将结果加到原数组中
    for(MYINT ir=0; ir<nr; ++ir){
        MYREAL tmp = sqrt(RTWO/(PI*rs[ir])) / dk0;
        for(MYINT i=0; i<SRC_M_NUM; ++i){
            for(MYINT v=0; v<INTEG_NUM; ++v){
                sum_J0[ir][i][v] += sum_J[ir][i][v] * tmp;

                if(calc_upar){
                    sum_uiz_J0[ir][i][v] += sum_uiz_J[ir][i][v] * tmp;
                    sum_uir_J0[ir][i][v] += sum_uir_J[ir][i][v] * tmp;
                }
            }
        }
    }


    BEFORE_RETURN:
    free(sum_J);
    if(sum_uiz_J) free(sum_uiz_J);
    if(sum_uir_J) free(sum_uir_J);


    free(iendkrs);

    return k;
}

