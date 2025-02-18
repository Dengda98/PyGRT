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

#include "fim.h"
#include "iostats.h"
#include "const.h"
#include "model.h"
#include "propagate.h"






MYREAL linear_filon_integ(
    const MODEL1D *mod1d, MYREAL dk, MYREAL kmax, MYREAL keps, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_EXP_J[nr][3][4], MYCOMPLEX sum_VF_J[nr][3][4],  
    MYCOMPLEX sum_HF_J[nr][3][4],  MYCOMPLEX sum_DC_J[nr][3][4],  
    FILE *(fstats[nr]))
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

    MYREAL k=0.0, r; 
    MYINT ik=0;

    // exp( i * (2*m+1)*pi / 4 )
    MYCOMPLEX ecoef[3];
    ecoef[0] =   INV_SQRT_TWO +  INV_SQRT_TWO*I;
    ecoef[1] = - INV_SQRT_TWO +  INV_SQRT_TWO*I;
    ecoef[2] = - INV_SQRT_TWO -  INV_SQRT_TWO*I;

    MYCOMPLEX coef[nr][3];
    for(MYINT ir=0; ir<nr; ++ir){
        r = rs[ir];
        for(MYINT m=0; m<3; ++m){
            // NOTICE: 这里对参数进行了设计（基于我个人理解，需进一步讨论）
            // 
            // 在(5.9.11)式中以及纪晨等(1995)的文章中是 2* (1 - cos(dk*r))， 
            // 推导过程是基于向外传播的Hankel函数Hm(x) = Jm(x) - i*Ym(x)，
            // 但由于推导中只保留(kr)的零阶项，导致Hm(x)只剩下实部Jm(x), 引入了误差，
            // 
            // 第一类Bessel函数的近似公式为 
            //             Jm(x) = sqrt(2/(pi*x)) * cos(x - m*pi/2 - pi/4)
            // 对cos函数运用欧拉公式:
            //             cos(x) = 0.5 * ( exp(j*x) + exp(-j*x) )
            // 此时带入待求积分式中，发现和(5.9.11)式相比多了0.5的系数，故这里系数2被"抵消"了
            // 
            // 另外提出了dk系数，故分母dk为二次
            coef[ir][m] = SQRT(RTWO/(PI*r)) * (RONE - COS(dk*r)) / (r*r*dk*dk) * ecoef[m];
        }
    }
    
    bool iendk, iendk0;

    // 每个震中距的k循环是否结束
    bool iendkrs[nr];
    for(MYINT ir=0; ir<nr; ++ir) iendkrs[ir] = false;

    // k循环 
    ik = 0;
    while(true){
        k += dk; 

        if(k > kmax) break;

        // 计算核函数 F(k, w)
        kernel(mod1d, omega, k, pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv); 

        // 震中距rs循环
        iendk = true;
        for(MYINT ir=0; ir<nr; ++ir){
            if(iendkrs[ir]) continue; // 该震中距下的波数k积分已收敛
            
            // F(k, w)*Jm(kr)k 的近似公式
            int_Pk_filon(
                k, rs[ir], 
                pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv,
                EXP_J, VF_J, HF_J, DC_J);

            // 记录积分结果
            if(fstats[ir]!=NULL){
                write_stats(
                    fstats[ir], k, 
                    pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv,
                    EXP_J, VF_J, HF_J, DC_J);
            }


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
            

            
        }  // end rs loop 
        
        ++ik;
        // 所有震中距的格林函数都已收敛
        if(iendk) break;

    } // end k loop

    // 乘上系数
    for(MYINT ir=0; ir<nr; ++ir){
        if(sum_EXP_J!=NULL){
            sum_EXP_J[ir][0][0] *= coef[ir][1];
            sum_EXP_J[ir][0][2] *= coef[ir][0];
        }
        
        if(sum_VF_J!=NULL){
            sum_VF_J[ir][0][0] *= coef[ir][1];
            sum_VF_J[ir][0][2] *= coef[ir][0];
        }
        
        if(sum_HF_J!=NULL){
            sum_HF_J[ir][1][0] *= coef[ir][0];
            sum_HF_J[ir][1][1] *= coef[ir][1];
            sum_HF_J[ir][1][2] *= coef[ir][1];
            sum_HF_J[ir][1][3] *= coef[ir][0];
        }

        if(sum_DC_J!=NULL){
            sum_DC_J[ir][0][0] *= coef[ir][1];
            sum_DC_J[ir][0][2] *= coef[ir][0];

            sum_DC_J[ir][1][0] *= coef[ir][0];
            sum_DC_J[ir][1][1] *= coef[ir][1];
            sum_DC_J[ir][1][2] *= coef[ir][1];
            sum_DC_J[ir][1][3] *= coef[ir][0];

            sum_DC_J[ir][2][0] *= coef[ir][1];
            sum_DC_J[ir][2][1] *= coef[ir][2];
            sum_DC_J[ir][2][2] *= coef[ir][2];
            sum_DC_J[ir][2][3] *= coef[ir][1];
        }
        
    }

    return k;
}



void int_Pk_filon(
    MYREAL k, MYREAL r, 
    const MYCOMPLEX EXP_qwv[3][3], const MYCOMPLEX VF_qwv[3][3], 
    const MYCOMPLEX HF_qwv[3][3],  const MYCOMPLEX DC_qwv[3][3], 
    MYCOMPLEX EXP_J[3][4], MYCOMPLEX VF_J[3][4], 
    MYCOMPLEX HF_J[3][4],  MYCOMPLEX DC_J[3][4] )
{
    MYREAL kr = k*r;
    MYCOMPLEX ekr = CEXP(- I*kr) * SQRT(k);
    
    if(EXP_qwv!=NULL){
    // 公式(5.6.22), 将公式分解为F(k,w)Jm(kr)k的形式
    // m=0 爆炸源
    EXP_J[0][0] = - EXP_qwv[0][0]*ekr;
    EXP_J[0][2] =   EXP_qwv[0][1]*ekr;
    }

    if(VF_qwv!=NULL){
    // m=0 垂直力源
    VF_J[0][0] = - VF_qwv[0][0]*ekr;
    VF_J[0][2] =   VF_qwv[0][1]*ekr;
    }

    if(HF_qwv!=NULL){
    // m=1 水平力源
    HF_J[1][0]  =   HF_qwv[1][0]*ekr;         // q1*J0
    HF_J[1][1]  = - (HF_qwv[1][0] + HF_qwv[1][2])*ekr/(kr);    // - (q1+v1)*J1/kr
    HF_J[1][2]  =   HF_qwv[1][1]*ekr;         // w1*J1
    HF_J[1][3]  = - HF_qwv[1][2]*ekr;         // -v1*J0
    }

    if(DC_qwv!=NULL){
    // m=0 双力偶源
    DC_J[0][0] = - DC_qwv[0][0]*ekr;
    DC_J[0][2] =   DC_qwv[0][1]*ekr;

    // m=1 双力偶源
    DC_J[1][0]  =   DC_qwv[1][0]*ekr;         // q1*J0
    DC_J[1][1]  = - (DC_qwv[1][0] + DC_qwv[1][2])*ekr/(kr);    // - (q1+v1)*J1/kr
    DC_J[1][2]  =   DC_qwv[1][1]*ekr;         // w1*J1
    DC_J[1][3]  = - DC_qwv[1][2]*ekr;         // -v1*J0

    // m=2 双力偶源
    DC_J[2][0]  =   DC_qwv[2][0]*ekr;         // q2*J1
    DC_J[2][1]  = - RTWO*(DC_qwv[2][0] + DC_qwv[2][2])*ekr/(kr);    // - (q2+v2)*J2/kr
    DC_J[2][2]  =   DC_qwv[2][1]*ekr;         // w2*J2
    DC_J[2][3]  = - DC_qwv[2][2]*ekr;         // -v2*J1
    }
}
