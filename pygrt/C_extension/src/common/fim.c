/**
 * @file   filon.c
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
    MYCOMPLEX sum_EXP_J0[nr][3][4], MYCOMPLEX sum_VF_J0[nr][3][4],  
    MYCOMPLEX sum_HF_J0[nr][3][4],  MYCOMPLEX sum_DC_J0[nr][3][4],  
    bool calc_upar,
    MYCOMPLEX sum_EXP_uiz_J0[nr][3][4], MYCOMPLEX sum_VF_uiz_J0[nr][3][4],  
    MYCOMPLEX sum_HF_uiz_J0[nr][3][4],  MYCOMPLEX sum_DC_uiz_J0[nr][3][4],  
    MYCOMPLEX sum_EXP_uir_J0[nr][3][4], MYCOMPLEX sum_VF_uir_J0[nr][3][4],  
    MYCOMPLEX sum_HF_uir_J0[nr][3][4],  MYCOMPLEX sum_DC_uir_J0[nr][3][4],  
    FILE *fstats, KernelFunc kerfunc)
{   
    // 从0开始，存储第二部分Filon积分的结果
    MYCOMPLEX (*sum_EXP_J)[3][4] = (sum_EXP_J0!=NULL)? (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_EXP_J)) : NULL;
    MYCOMPLEX (*sum_VF_J)[3][4] =  (sum_VF_J0!=NULL)?  (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_VF_J)) : NULL;
    MYCOMPLEX (*sum_HF_J)[3][4] =  (sum_HF_J0!=NULL)?  (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_HF_J)) : NULL;
    MYCOMPLEX (*sum_DC_J)[3][4] =  (sum_DC_J0!=NULL)?  (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_DC_J)) : NULL;

    MYCOMPLEX (*sum_EXP_uiz_J)[3][4] = (sum_EXP_uiz_J0!=NULL)? (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_EXP_uiz_J)) : NULL;
    MYCOMPLEX (*sum_VF_uiz_J)[3][4] =  (sum_VF_uiz_J0!=NULL)?  (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_VF_uiz_J)) : NULL;
    MYCOMPLEX (*sum_HF_uiz_J)[3][4] =  (sum_HF_uiz_J0!=NULL)?  (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_HF_uiz_J)) : NULL;
    MYCOMPLEX (*sum_DC_uiz_J)[3][4] =  (sum_DC_uiz_J0!=NULL)?  (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_DC_uiz_J)) : NULL;

    MYCOMPLEX (*sum_EXP_uir_J)[3][4] = (sum_EXP_uir_J0!=NULL)? (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_EXP_uir_J)) : NULL;
    MYCOMPLEX (*sum_VF_uir_J)[3][4] =  (sum_VF_uir_J0!=NULL)?  (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_VF_uir_J)) : NULL;
    MYCOMPLEX (*sum_HF_uir_J)[3][4] =  (sum_HF_uir_J0!=NULL)?  (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_HF_uir_J)) : NULL;
    MYCOMPLEX (*sum_DC_uir_J)[3][4] =  (sum_DC_uir_J0!=NULL)?  (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_DC_uir_J)) : NULL;


    MYCOMPLEX EXP_J[3][4], VF_J[3][4], HF_J[3][4],  DC_J[3][4];
    for(MYINT m=0; m<3; ++m){
        for(MYINT v=0; v<4; ++v){
            EXP_J[m][v] = VF_J[m][v] = HF_J[m][v] = DC_J[m][v] = 0.0;
        }
    }


    MYCOMPLEX EXP_qwv[3][3], VF_qwv[3][3], HF_qwv[3][3], DC_qwv[3][3]; // 不同震源的核函数
    MYCOMPLEX (*pEXP_qwv)[3] = (sum_EXP_J!=NULL)? EXP_qwv : NULL;
    MYCOMPLEX (*pVF_qwv)[3]  = (sum_VF_J!=NULL)?  VF_qwv  : NULL;
    MYCOMPLEX (*pHF_qwv)[3]  = (sum_HF_J!=NULL)?  HF_qwv  : NULL;
    MYCOMPLEX (*pDC_qwv)[3]  = (sum_DC_J!=NULL)?  DC_qwv  : NULL;

    MYCOMPLEX EXP_uiz_qwv[3][3], VF_uiz_qwv[3][3], HF_uiz_qwv[3][3], DC_uiz_qwv[3][3]; 
    MYCOMPLEX (*pEXP_uiz_qwv)[3] = (sum_EXP_uiz_J!=NULL)? EXP_uiz_qwv : NULL;
    MYCOMPLEX (*pVF_uiz_qwv)[3]  = (sum_VF_uiz_J!=NULL)?  VF_uiz_qwv  : NULL;
    MYCOMPLEX (*pHF_uiz_qwv)[3]  = (sum_HF_uiz_J!=NULL)?  HF_uiz_qwv  : NULL;
    MYCOMPLEX (*pDC_uiz_qwv)[3]  = (sum_DC_uiz_J!=NULL)?  DC_uiz_qwv  : NULL;

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
        kerfunc(mod1d, omega, k, pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv,
                calc_upar, pEXP_uiz_qwv, pVF_uiz_qwv, pHF_uiz_qwv, pDC_uiz_qwv); 

        // 记录积分结果
        if(fstats!=NULL){
            write_stats(
                fstats, k, 
                pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv);
        }

        // 震中距rs循环
        iendk = true;
        for(MYINT ir=0; ir<nr; ++ir){
            if(iendkrs[ir]) continue; // 该震中距下的波数k积分已收敛
            
            // F(k, w)*Jm(kr)k 的近似公式, sqrt(k) * F(k,w) * cos
            int_Pk_filon(
                k, rs[ir], true,
                pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv, false,
                EXP_J, VF_J, HF_J, DC_J);


            iendk0 = true;
            for(MYINT m=0; m<3; ++m){
                for(MYINT v=0; v<4; ++v){
                    if(sum_EXP_J!=NULL) sum_EXP_J[ir][m][v] += EXP_J[m][v];
                    if(sum_VF_J!=NULL)  sum_VF_J[ir][m][v]  += VF_J[m][v];
                    if(sum_HF_J!=NULL)  sum_HF_J[ir][m][v]  += HF_J[m][v];
                    if(sum_DC_J!=NULL)  sum_DC_J[ir][m][v]  += DC_J[m][v];

                    if(keps > 0.0){
                        // 判断是否达到收敛条件
                        if(sum_EXP_J!=NULL && m==0 && (v==0||v==2)) iendk0 = iendk0 && (CABS(EXP_J[m][v])/ CABS(sum_EXP_J[ir][m][v]) <= keps);
                        if(sum_VF_J!=NULL  && m==0 && (v==0||v==2)) iendk0 = iendk0 && (CABS(VF_J[m][v]) / CABS(sum_VF_J[ir][m][v])  <= keps);
                        if(sum_HF_J!=NULL  && m==1) iendk0 = iendk0 && (CABS(HF_J[m][v]) / CABS(sum_HF_J[ir][m][v])  <= keps);
                        if(sum_DC_J!=NULL  && ((m==0 && (v==0||v==2)) || m!=0)) iendk0 = iendk0 && (CABS(DC_J[m][v]) / CABS(sum_DC_J[ir][m][v])  <= keps);
                    } 
                }
            }
            
            if(keps > 0.0){
                iendkrs[ir] = iendk0;
                iendk = iendk && iendkrs[ir];
            } else {
                iendk = iendkrs[ir] = false;
            }
            

            // ---------------- 位移空间导数，EXP_J, VF_J, HF_J, DC_J数组重复利用 --------------------------
            if(calc_upar){
                // ------------------------------- ui_z -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                int_Pk_filon(k, rs[ir], true,
                       pEXP_uiz_qwv, pVF_uiz_qwv, pHF_uiz_qwv, pDC_uiz_qwv, false,
                       EXP_J, VF_J, HF_J, DC_J);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uiz也收敛
                for(MYINT m=0; m<3; ++m){
                    for(MYINT v=0; v<4; ++v){
                        if(sum_EXP_uiz_J!=NULL) sum_EXP_uiz_J[ir][m][v] += EXP_J[m][v];
                        if(sum_VF_uiz_J!=NULL)  sum_VF_uiz_J[ir][m][v]  += VF_J[m][v];
                        if(sum_HF_uiz_J!=NULL)  sum_HF_uiz_J[ir][m][v]  += HF_J[m][v];
                        if(sum_DC_uiz_J!=NULL)  sum_DC_uiz_J[ir][m][v]  += DC_J[m][v];
                    }
                }


                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                int_Pk_filon(k, rs[ir], true,
                       pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv, true,
                       EXP_J, VF_J, HF_J, DC_J);
                
                // keps不参与计算位移空间导数的积分，背后逻辑认为u收敛，则uir也收敛
                for(MYINT m=0; m<3; ++m){
                    for(MYINT v=0; v<4; ++v){
                        if(sum_EXP_uir_J!=NULL) sum_EXP_uir_J[ir][m][v] += EXP_J[m][v];
                        if(sum_VF_uir_J!=NULL)  sum_VF_uir_J[ir][m][v]  += VF_J[m][v];
                        if(sum_HF_uir_J!=NULL)  sum_HF_uir_J[ir][m][v]  += HF_J[m][v];
                        if(sum_DC_uir_J!=NULL)  sum_DC_uir_J[ir][m][v]  += DC_J[m][v];
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
        MYREAL tmp = RTWO*(RONE - COS(dk*rs[ir])) / (rs[ir]*rs[ir]*dk);
        for(MYINT m=0; m<3; ++m){
            for(MYINT v=0; v<4; ++v){
                if(sum_EXP_J!=NULL) sum_EXP_J[ir][m][v] *= tmp;
                if(sum_VF_J!=NULL)  sum_VF_J[ir][m][v]  *= tmp;
                if(sum_HF_J!=NULL)  sum_HF_J[ir][m][v]  *= tmp;
                if(sum_DC_J!=NULL)  sum_DC_J[ir][m][v]  *= tmp;

                if(calc_upar){
                    if(sum_EXP_uiz_J!=NULL) sum_EXP_uiz_J[ir][m][v] *= tmp;
                    if(sum_VF_uiz_J!=NULL)  sum_VF_uiz_J[ir][m][v]  *= tmp;
                    if(sum_HF_uiz_J!=NULL)  sum_HF_uiz_J[ir][m][v]  *= tmp;
                    if(sum_DC_uiz_J!=NULL)  sum_DC_uiz_J[ir][m][v]  *= tmp;

                    if(sum_EXP_uir_J!=NULL) sum_EXP_uir_J[ir][m][v] *= tmp;
                    if(sum_VF_uir_J!=NULL)  sum_VF_uir_J[ir][m][v]  *= tmp;
                    if(sum_HF_uir_J!=NULL)  sum_HF_uir_J[ir][m][v]  *= tmp;
                    if(sum_DC_uir_J!=NULL)  sum_DC_uir_J[ir][m][v]  *= tmp;
                }
            }
        }
    }

    // -------------------------------------------------------------------------------
    // 计算余项, [2]表示k积分的第一个点和最后一个点
    MYCOMPLEX EXP_Gc[2][3][4] = {0}, EXP_Gs[2][3][4] = {0};
    MYCOMPLEX VF_Gc[2][3][4] = {0}, VF_Gs[2][3][4] = {0};
    MYCOMPLEX HF_Gc[2][3][4] = {0}, HF_Gs[2][3][4] = {0};
    MYCOMPLEX DC_Gc[2][3][4] = {0}, DC_Gs[2][3][4] = {0};

    // for(MYINT s=0; s<2; ++s){
    //     for(MYINT m=0; m<3; ++m){
    //         for(MYINT v=0; v<4; ++v){
    //             EXP_Gc[s][m][v] = VF_Gc[s][m][v] = HF_Gc[s][m][v] = DC_Gc[s][m][v] = 0.0;
    //             EXP_Gs[s][m][v] = VF_Gs[s][m][v] = HF_Gs[s][m][v] = DC_Gs[s][m][v] = 0.0;
    //         }
    //     }
    // }

    // 计算来自第一个点和最后一个点的余项
    for(MYINT iik=0; iik<2; ++iik){ 
        MYREAL k0N;
        MYINT sgn;
        if(0==iik)       {k0N = k0+dk; sgn = RONE;}
        else if(1==iik)  {k0N = k;  sgn = -RONE;}
        else {
            fprintf(stderr, "Filon error.\n");
            exit(EXIT_FAILURE);
        }

        // 计算核函数 F(k, w)
        kerfunc(mod1d, omega, k0N, pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv,
            calc_upar, pEXP_uiz_qwv, pVF_uiz_qwv, pHF_uiz_qwv, pDC_uiz_qwv); 

        for(MYINT ir=0; ir<nr; ++ir){
            // Gc
            int_Pk_filon(
                k0N, rs[ir], true,
                pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv, false,
                EXP_Gc[iik], VF_Gc[iik], HF_Gc[iik], DC_Gc[iik]);
            
            // Gs
            int_Pk_filon(
                k0N, rs[ir], false,
                pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv, false,
                EXP_Gs[iik], VF_Gs[iik], HF_Gs[iik], DC_Gs[iik]);
            
            MYREAL tmp = RONE / (rs[ir]*rs[ir]*dk);
            MYREAL tmpc = tmp * (RONE - COS(dk*rs[ir]));
            MYREAL tmps = sgn * tmp * SIN(dk*rs[ir]);

            for(MYINT m=0; m<3; ++m){
                for(MYINT v=0; v<4; ++v){
                    if(sum_EXP_J!=NULL) sum_EXP_J[ir][m][v] += (- tmpc*EXP_Gc[iik][m][v] + tmps*EXP_Gs[iik][m][v] - sgn*EXP_Gs[iik][m][v]/rs[ir]);
                    if(sum_VF_J!=NULL)  sum_VF_J[ir][m][v]  += (- tmpc*VF_Gc[iik][m][v] + tmps*VF_Gs[iik][m][v] - sgn*VF_Gs[iik][m][v]/rs[ir]);
                    if(sum_HF_J!=NULL)  sum_HF_J[ir][m][v]  += (- tmpc*HF_Gc[iik][m][v] + tmps*HF_Gs[iik][m][v] - sgn*HF_Gs[iik][m][v]/rs[ir]);
                    if(sum_DC_J!=NULL)  sum_DC_J[ir][m][v]  += (- tmpc*DC_Gc[iik][m][v] + tmps*DC_Gs[iik][m][v] - sgn*DC_Gs[iik][m][v]/rs[ir]);
                }
            }


            // ---------------- 位移空间导数，EXP_Gc/s, VF_Gc/s, HF_Gc/s, DC_Gc/s数组重复利用 --------------------------
            if(calc_upar){
                // ------------------------------- ui_z -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                // Gc
                int_Pk_filon(
                    k0N, rs[ir], true,
                    pEXP_uiz_qwv, pVF_uiz_qwv, pHF_uiz_qwv, pDC_uiz_qwv, false,
                    EXP_Gc[iik], VF_Gc[iik], HF_Gc[iik], DC_Gc[iik]);
                
                // Gs
                int_Pk_filon(
                    k0N, rs[ir], false,
                    pEXP_uiz_qwv, pVF_uiz_qwv, pHF_uiz_qwv, pDC_uiz_qwv, false,
                    EXP_Gs[iik], VF_Gs[iik], HF_Gs[iik], DC_Gs[iik]);

                for(MYINT m=0; m<3; ++m){
                    for(MYINT v=0; v<4; ++v){
                        if(sum_EXP_uiz_J!=NULL) sum_EXP_uiz_J[ir][m][v] += (- tmpc*EXP_Gc[iik][m][v] + tmps*EXP_Gs[iik][m][v] - sgn*EXP_Gs[iik][m][v]/rs[ir]);
                        if(sum_VF_uiz_J!=NULL)  sum_VF_uiz_J[ir][m][v]  += (- tmpc*VF_Gc[iik][m][v] + tmps*VF_Gs[iik][m][v] - sgn*VF_Gs[iik][m][v]/rs[ir]);
                        if(sum_HF_uiz_J!=NULL)  sum_HF_uiz_J[ir][m][v]  += (- tmpc*HF_Gc[iik][m][v] + tmps*HF_Gs[iik][m][v] - sgn*HF_Gs[iik][m][v]/rs[ir]);
                        if(sum_DC_uiz_J!=NULL)  sum_DC_uiz_J[ir][m][v]  += (- tmpc*DC_Gc[iik][m][v] + tmps*DC_Gs[iik][m][v] - sgn*DC_Gs[iik][m][v]/rs[ir]);
                    }
                }


                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                // Gc
                int_Pk_filon(
                    k0N, rs[ir], true,
                    pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv, true,
                    EXP_Gc[iik], VF_Gc[iik], HF_Gc[iik], DC_Gc[iik]);
                
                // Gs
                int_Pk_filon(
                    k0N, rs[ir], false,
                    pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv, true,
                    EXP_Gs[iik], VF_Gs[iik], HF_Gs[iik], DC_Gs[iik]);
                
                for(MYINT m=0; m<3; ++m){
                    for(MYINT v=0; v<4; ++v){
                        if(sum_EXP_uir_J!=NULL) sum_EXP_uir_J[ir][m][v] += (- tmpc*EXP_Gc[iik][m][v] + tmps*EXP_Gs[iik][m][v] - sgn*EXP_Gs[iik][m][v]/rs[ir]);
                        if(sum_VF_uir_J!=NULL)  sum_VF_uir_J[ir][m][v]  += (- tmpc*VF_Gc[iik][m][v] + tmps*VF_Gs[iik][m][v] - sgn*VF_Gs[iik][m][v]/rs[ir]);
                        if(sum_HF_uir_J!=NULL)  sum_HF_uir_J[ir][m][v]  += (- tmpc*HF_Gc[iik][m][v] + tmps*HF_Gs[iik][m][v] - sgn*HF_Gs[iik][m][v]/rs[ir]);
                        if(sum_DC_uir_J!=NULL)  sum_DC_uir_J[ir][m][v]  += (- tmpc*DC_Gc[iik][m][v] + tmps*DC_Gs[iik][m][v] - sgn*DC_Gs[iik][m][v]/rs[ir]);
                    }
                }
            } // END if calc_upar
          
        }  // END rs loop
    
    }  // END k 2-points loop

    // 乘上总系数 SQRT(RTWO/(PI*r)) / dk0,  除dks0是在该函数外还会再乘dk0
    for(MYINT ir=0; ir<nr; ++ir){
        MYREAL tmp = SQRT(RTWO/(PI*rs[ir])) / dk0;
        for(MYINT m=0; m<3; ++m){
            for(MYINT v=0; v<4; ++v){
                if(sum_EXP_J!=NULL) sum_EXP_J[ir][m][v] *= tmp;
                if(sum_VF_J!=NULL)  sum_VF_J[ir][m][v]  *= tmp;
                if(sum_HF_J!=NULL)  sum_HF_J[ir][m][v]  *= tmp;
                if(sum_DC_J!=NULL)  sum_DC_J[ir][m][v]  *= tmp;

                if(calc_upar){
                    if(sum_EXP_uiz_J!=NULL) sum_EXP_uiz_J[ir][m][v] *= tmp;
                    if(sum_VF_uiz_J!=NULL)  sum_VF_uiz_J[ir][m][v]  *= tmp;
                    if(sum_HF_uiz_J!=NULL)  sum_HF_uiz_J[ir][m][v]  *= tmp;
                    if(sum_DC_uiz_J!=NULL)  sum_DC_uiz_J[ir][m][v]  *= tmp;

                    if(sum_EXP_uir_J!=NULL) sum_EXP_uir_J[ir][m][v] *= tmp;
                    if(sum_VF_uir_J!=NULL)  sum_VF_uir_J[ir][m][v]  *= tmp;
                    if(sum_HF_uir_J!=NULL)  sum_HF_uir_J[ir][m][v]  *= tmp;
                    if(sum_DC_uir_J!=NULL)  sum_DC_uir_J[ir][m][v]  *= tmp;
                }
            }
        }
    }


    // 将结果加到原数组中
    for(MYINT ir=0; ir<nr; ++ir){
        for(MYINT m=0; m<3; ++m){
            for(MYINT v=0; v<4; ++v){
                if(sum_EXP_J!=NULL) sum_EXP_J0[ir][m][v] += sum_EXP_J[ir][m][v];
                if(sum_VF_J!=NULL)  sum_VF_J0[ir][m][v]  += sum_VF_J[ir][m][v];
                if(sum_HF_J!=NULL)  sum_HF_J0[ir][m][v]  += sum_HF_J[ir][m][v];
                if(sum_DC_J!=NULL)  sum_DC_J0[ir][m][v]  += sum_DC_J[ir][m][v] ;

                if(calc_upar){
                    if(sum_EXP_uiz_J!=NULL) sum_EXP_uiz_J0[ir][m][v] += sum_EXP_uiz_J[ir][m][v];
                    if(sum_VF_uiz_J!=NULL)  sum_VF_uiz_J0[ir][m][v]  += sum_VF_uiz_J[ir][m][v];
                    if(sum_HF_uiz_J!=NULL)  sum_HF_uiz_J0[ir][m][v]  += sum_HF_uiz_J[ir][m][v];
                    if(sum_DC_uiz_J!=NULL)  sum_DC_uiz_J0[ir][m][v]  += sum_DC_uiz_J[ir][m][v];

                    if(sum_EXP_uir_J!=NULL) sum_EXP_uir_J0[ir][m][v] += sum_EXP_uir_J[ir][m][v];
                    if(sum_VF_uir_J!=NULL)  sum_VF_uir_J0[ir][m][v]  += sum_VF_uir_J[ir][m][v];
                    if(sum_HF_uir_J!=NULL)  sum_HF_uir_J0[ir][m][v]  += sum_HF_uir_J[ir][m][v];
                    if(sum_DC_uir_J!=NULL)  sum_DC_uir_J0[ir][m][v]  += sum_DC_uir_J[ir][m][v];
                }
            }
        }
    }

    
    if(sum_EXP_J) free(sum_EXP_J); 
    if(sum_VF_J) free(sum_VF_J); 
    if(sum_HF_J) free(sum_HF_J); 
    if(sum_DC_J) free(sum_DC_J);
    if(sum_EXP_uiz_J) free(sum_EXP_uiz_J); 
    if(sum_VF_uiz_J) free(sum_VF_uiz_J); 
    if(sum_HF_uiz_J) free(sum_HF_uiz_J); 
    if(sum_DC_uiz_J) free(sum_DC_uiz_J);
    if(sum_EXP_uir_J) free(sum_EXP_uir_J); 
    if(sum_VF_uir_J) free(sum_VF_uir_J); 
    if(sum_HF_uir_J) free(sum_HF_uir_J); 
    if(sum_DC_uir_J) free(sum_DC_uir_J);

    free(iendkrs);

    return k;
}

