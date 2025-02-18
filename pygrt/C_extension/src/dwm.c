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

#include "dwm.h"
#include "const.h"
#include "model.h"
#include "propagate.h"
#include "iostats.h"


MYREAL discrete_integ(
    const MODEL1D *mod1d, MYREAL dk, MYREAL kmax, MYREAL keps, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_EXP_J[nr][3][4], MYCOMPLEX sum_VF_J[nr][3][4],  
    MYCOMPLEX sum_HF_J[nr][3][4],  MYCOMPLEX sum_DC_J[nr][3][4],  
    FILE *(fstats[nr]))
{
    MYCOMPLEX EXP_J[3][4], VF_J[3][4], HF_J[3][4],  DC_J[3][4];

    // 不同震源的核函数 F(k, w) 
    // 第一个维度3代表阶数m=0,1,2，第二个维度3代表三类系数qm,wm,vm 
    // 实际上对于不同震源只有特定阶数/系数才有值，不需要建立3x3的小矩阵，
    // 但这里还是为了方便可读性，牺牲了部分性能 
    MYCOMPLEX EXP_qwv[3][3], VF_qwv[3][3], HF_qwv[3][3], DC_qwv[3][3]; 
    MYCOMPLEX (*pEXP_qwv)[3] = (sum_EXP_J!=NULL)? EXP_qwv : NULL;
    MYCOMPLEX (*pVF_qwv)[3]  = (sum_VF_J!=NULL)?  VF_qwv  : NULL;
    MYCOMPLEX (*pHF_qwv)[3]  = (sum_HF_J!=NULL)?  HF_qwv  : NULL;
    MYCOMPLEX (*pDC_qwv)[3]  = (sum_DC_J!=NULL)?  DC_qwv  : NULL;
    
    MYREAL k = 0.0;
    MYINT ik = 0;

    // 所有震中距的k循环是否结束
    bool iendk = true;

    // 每个震中距的k循环是否结束
    bool iendkrs[nr], iendk0=false;
    for(MYINT ir=0; ir<nr; ++ir) iendkrs[ir] = false;
    

    // 波数k循环 (5.9.2)
    while(true){
        k += dk; 

        if(k > kmax)  break;

        // printf("w=%15.5e, ik=%d\n", CREAL(omega), ik);
        // 计算核函数 F(k, w)
        kernel(mod1d, omega, k, pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv); 


        // 震中距rs循环
        iendk = true;
        for(MYINT ir=0; ir<nr; ++ir){
            if(iendkrs[ir]) continue; // 该震中距下的波数k积分已收敛

            for(MYINT m=0; m<3; ++m){
                for(MYINT v=0; v<4; ++v){
                    EXP_J[m][v] = VF_J[m][v] = HF_J[m][v] = DC_J[m][v] = CZERO;
                }
            }
            
            // 计算被积函数一项 F(k,w)Jm(kr)k
            int_Pk(k, rs[ir], 
                    pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv,
                    EXP_J, VF_J, HF_J, DC_J);

            // 记录积分结果
            if(fstats[ir]!=NULL){
                write_stats(
                    fstats[ir], k, 
                    EXP_qwv, VF_qwv, HF_qwv, DC_qwv,
                    EXP_J, VF_J, HF_J, DC_J);
            }

            
            iendk0 = true;
            for(MYINT m=0; m<3; ++m){
                for(MYINT v=0; v<4; ++v){
                    if(sum_EXP_J!=NULL) sum_EXP_J[ir][m][v] += EXP_J[m][v];
                    if(sum_VF_J!=NULL)  sum_VF_J[ir][m][v]  += VF_J[m][v];
                    if(sum_HF_J!=NULL)  sum_HF_J[ir][m][v]  += HF_J[m][v];
                    if(sum_DC_J!=NULL)  sum_DC_J[ir][m][v]  += DC_J[m][v];

                    if(keps > RZERO){
                        // 判断是否达到收敛条件
                        if(sum_EXP_J!=NULL && m==0 && (v==0||v==2)) iendk0 = iendk0 && (CABS(EXP_J[m][v])/ CABS(sum_EXP_J[ir][m][v]) <= keps);
                        if(sum_VF_J!=NULL  && m==0 && (v==0||v==2)) iendk0 = iendk0 && (CABS(VF_J[m][v]) / CABS(sum_VF_J[ir][m][v])  <= keps);
                        if(sum_HF_J!=NULL && m==1)                  iendk0 = iendk0 && (CABS(HF_J[m][v]) / CABS(sum_HF_J[ir][m][v])  <= keps);
                        if(sum_DC_J!=NULL && ((m==0 && (v==0||v==2)) || m!=0)) iendk0 = iendk0 && (CABS(DC_J[m][v]) / CABS(sum_DC_J[ir][m][v])  <= keps);
                    } 
                }
            }
            
            if(keps > RZERO){
                iendkrs[ir] = iendk0;
                iendk = iendk && iendkrs[ir];
            } else {
                iendk = iendkrs[ir] = false;
            }
            
            

        } // END rs loop

        ++ik;

        // 所有震中距的格林函数都已收敛
        if(iendk) break;

    } // END k loop

    // printf("w=%15.5e, ik=%d\n", CREAL(omega), ik);

    return k;

}

