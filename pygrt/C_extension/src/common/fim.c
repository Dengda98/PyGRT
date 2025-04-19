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
    const MODEL1D *mod1d, MYREAL dk, MYREAL kmax, MYREAL keps, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_EXP_J[nr][3][4], MYCOMPLEX sum_VF_J[nr][3][4],  
    MYCOMPLEX sum_HF_J[nr][3][4],  MYCOMPLEX sum_DC_J[nr][3][4],  
    bool calc_upar,
    MYCOMPLEX sum_EXP_uiz_J[nr][3][4], MYCOMPLEX sum_VF_uiz_J[nr][3][4],  
    MYCOMPLEX sum_HF_uiz_J[nr][3][4],  MYCOMPLEX sum_DC_uiz_J[nr][3][4],  
    MYCOMPLEX sum_EXP_uir_J[nr][3][4], MYCOMPLEX sum_VF_uir_J[nr][3][4],  
    MYCOMPLEX sum_HF_uir_J[nr][3][4],  MYCOMPLEX sum_DC_uir_J[nr][3][4],  
    FILE *fstats, KernelFunc kerfunc)
{   
    for(MYINT ir=0; ir<nr; ++ir){
        for(MYINT m=0; m<3; ++m){
            for(MYINT v=0; v<4; ++v){
                if(sum_EXP_J!=NULL) sum_EXP_J[ir][m][v] = 0.0;
                if(sum_EXP_J!=NULL) sum_VF_J[ir][m][v]  = 0.0;
                if(sum_HF_J!=NULL)  sum_HF_J[ir][m][v]  = 0.0;
                if(sum_DC_J!=NULL)  sum_DC_J[ir][m][v]  = 0.0;
            }
        }
    }

    MYCOMPLEX EXP_J[3][4], VF_J[3][4], HF_J[3][4],  DC_J[3][4];
    for(MYINT ir=0; ir<nr; ++ir){
        for(MYINT m=0; m<3; ++m){
            for(MYINT v=0; v<4; ++v){
                EXP_J[m][v] = VF_J[m][v] = HF_J[m][v] = DC_J[m][v] = 0.0;
            }
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

    MYREAL k=0.0, r; 
    MYINT ik=0;

    MYCOMPLEX coef[nr];
    for(MYINT ir=0; ir<nr; ++ir){
        r = rs[ir];
        for(MYINT m=0; m<3; ++m){
            // NOTICE: 这里对参数进行了设计（基于我个人理解，需进一步讨论）
            // 
            // 在(5.9.11)式中以及纪晨等(1995)的文章中是 2* (1 - cos(dk*r))， 
            // 推导过程是基于向外传播的Hankel函数Hm(x) = Jm(x) + i*Ym(x)，
            // 但由于推导中只保留(kr)的零阶项，导致Hm(x)只剩下实部Jm(x), 引入了误差，
            // 
            // 第一类Bessel函数的近似公式为 
            //             Jm(x) = sqrt(2/(pi*x)) * cos(x - m*pi/2 - pi/4)
            // 对cos函数运用欧拉公式:
            //             cos(x) = 0.5 * ( exp(j*x) + exp(-j*x) )
            // 此时带入待求积分式中，发现和(5.9.11)式相比多了0.5的系数，故这里系数2被"抵消"了,
            // 本质原因其实是bessel函数被替换为复数形式的exp(-j*x)，如果使用欧拉公式进一步调整
            // 则bessel函数替换为余弦函数，系数2可以保留。
            // 
            // 另外提出了dk系数，故分母dk为二次
            coef[ir] = SQRT(RTWO/(PI*r)) * (RONE - COS(dk*r)) / (r*r*dk*dk);
        }
    }
    
    bool iendk, iendk0;

    // 每个震中距的k循环是否结束
    bool *iendkrs = (bool *)malloc(nr * sizeof(bool));
    for(MYINT ir=0; ir<nr; ++ir) iendkrs[ir] = false;

    // k循环 
    ik = 0;
    while(true){
        
        if(k > kmax) break;
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
            
            // F(k, w)*Jm(kr)k 的近似公式
            int_Pk_filon(
                k, rs[ir], 
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
                int_Pk_filon(k, rs[ir], 
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
                int_Pk_filon(k, rs[ir], 
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

    // 乘上系数
    for(MYINT ir=0; ir<nr; ++ir){
        for(MYINT m=0; m<3; ++m){
            for(MYINT v=0; v<4; ++v){
                if(sum_EXP_J!=NULL) sum_EXP_J[ir][m][v] *= coef[ir];
                if(sum_VF_J!=NULL)  sum_VF_J[ir][m][v]  *= coef[ir];
                if(sum_HF_J!=NULL)  sum_HF_J[ir][m][v]  *= coef[ir];
                if(sum_DC_J!=NULL)  sum_DC_J[ir][m][v]  *= coef[ir];

                if(calc_upar){
                    if(sum_EXP_uiz_J!=NULL) sum_EXP_uiz_J[ir][m][v] *= coef[ir];
                    if(sum_VF_uiz_J!=NULL)  sum_VF_uiz_J[ir][m][v]  *= coef[ir];
                    if(sum_HF_uiz_J!=NULL)  sum_HF_uiz_J[ir][m][v]  *= coef[ir];
                    if(sum_DC_uiz_J!=NULL)  sum_DC_uiz_J[ir][m][v]  *= coef[ir];

                    if(sum_EXP_uir_J!=NULL) sum_EXP_uir_J[ir][m][v] *= coef[ir];
                    if(sum_VF_uir_J!=NULL)  sum_VF_uir_J[ir][m][v]  *= coef[ir];
                    if(sum_HF_uir_J!=NULL)  sum_HF_uir_J[ir][m][v]  *= coef[ir];
                    if(sum_DC_uir_J!=NULL)  sum_DC_uir_J[ir][m][v]  *= coef[ir];
                }
            }
        }
    }
    

    free(iendkrs);

    return k;
}



void int_Pk_filon(
    MYREAL k, MYREAL r, 
    const MYCOMPLEX EXP_qwv[3][3], const MYCOMPLEX VF_qwv[3][3], 
    const MYCOMPLEX HF_qwv[3][3],  const MYCOMPLEX DC_qwv[3][3], 
    bool calc_uir,
    MYCOMPLEX EXP_J[3][4], MYCOMPLEX VF_J[3][4], 
    MYCOMPLEX HF_J[3][4],  MYCOMPLEX DC_J[3][4] )
{
    MYREAL kr = k*r;
    MYREAL kr_inv = RONE/kr;
    MYREAL kcoef = SQRT(k);
    MYCOMPLEX bj0k, bj1k, bj2k;

    MYCOMPLEX J1coef, J2coef;

    if(calc_uir){
        kcoef *= k;

        bj0k = - CEXP(-I*(kr - THREEQUARTERPI));
        bj1k = - CEXP(-I*(kr - FIVEQUARTERPI));
        bj2k = - CEXP(-I*(kr - SEVENQUARTERPI));
    } else {
        bj0k = CEXP(-I*(kr - QUARTERPI));
        bj1k = CEXP(-I*(kr - THREEQUARTERPI));
        bj2k = CEXP(-I*(kr - FIVEQUARTERPI));
    }
    J1coef = bj1k*kr_inv;
    J2coef = bj2k*kr_inv;

    J1coef *= kcoef;
    J2coef *= kcoef;

    bj0k *= kcoef;
    bj1k *= kcoef;
    bj2k *= kcoef;

    
    if(EXP_qwv!=NULL){
    // 公式(5.6.22), 将公式分解为F(k,w)Jm(kr)k的形式
    // m=0 爆炸源
    EXP_J[0][0] = - EXP_qwv[0][0]*bj1k;
    EXP_J[0][2] =   EXP_qwv[0][1]*bj0k;
    }

    if(VF_qwv!=NULL){
    // m=0 垂直力源
    VF_J[0][0] = - VF_qwv[0][0]*bj1k;
    VF_J[0][2] =   VF_qwv[0][1]*bj0k;
    }

    if(HF_qwv!=NULL){
    // m=1 水平力源
    HF_J[1][0]  =   HF_qwv[1][0]*bj0k;         // q1*J0*k
    HF_J[1][1]  = - (HF_qwv[1][0] + HF_qwv[1][2])*J1coef;    // - (q1+v1)*J1*k/kr
    HF_J[1][2]  =   HF_qwv[1][1]*bj1k;         // w1*J1*k
    HF_J[1][3]  = - HF_qwv[1][2]*bj0k;         // -v1*J0*k
    }

    if(DC_qwv!=NULL){
    // m=0 剪切源
    DC_J[0][0] = - DC_qwv[0][0]*bj1k;
    DC_J[0][2] =   DC_qwv[0][1]*bj0k;

    // m=1 剪切源
    DC_J[1][0]  =   DC_qwv[1][0]*bj0k;         // q1*J0*k
    DC_J[1][1]  = - (DC_qwv[1][0] + DC_qwv[1][2])*J1coef;    // - (q1+v1)*J1*k/kr
    DC_J[1][2]  =   DC_qwv[1][1]*bj1k;         // w1*J1*k
    DC_J[1][3]  = - DC_qwv[1][2]*bj0k;         // -v1*J0*k

    // m=2 剪切源
    DC_J[2][0]  =   DC_qwv[2][0]*bj1k;         // q2*J1*k
    DC_J[2][1]  = - RTWO*(DC_qwv[2][0] + DC_qwv[2][2])*J2coef;    // - (q2+v2)*J2*k/kr
    DC_J[2][2]  =   DC_qwv[2][1]*bj2k;         // w2*J2*k
    DC_J[2][3]  = - DC_qwv[2][2]*bj1k;         // -v2*J1*k
    }
}
