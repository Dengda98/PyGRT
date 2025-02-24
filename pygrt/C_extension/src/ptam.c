/**
 * @file   ptam.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 峰谷平均法 ，参考：
 * 
 *         1. 张海明. 2021. 地震学中的Lamb问题（上）. 科学出版社
 *         2. Zhang, H. M., Chen, X. F., & Chang, S. (2003). 
 *               An efficient numerical method for computing synthetic seismograms 
 *               for a layered half-space with sources and receivers at close or same depths. 
 *               Seismic motion, lithospheric structures, earthquake and volcanic sources: 
 *               The Keiiti Aki volume, 467-486.
 * 
 */

#include <stdio.h> 
#include <complex.h>
#include <stdlib.h>

#include "ptam.h"
#include "iostats.h"
#include "const.h"
#include "model.h"
#include "propagate.h"
#include "quadratic.h"


/**
 * 处理并确定波峰或波谷                                    
 * 
 * @param ir        震中距索引                          
 * @param m         Bessel函数阶                          
 * @param v         积分形式索引                          
 * @param maxNpt    最大峰谷数                                
 * @param maxnwait  最大等待次数                        
 * @param k         波数                             
 * @param dk        波数步长                              
 * @param J3        存储的采样幅值数组                  
 * @param kpt       存储的采样值对应的波数数组             
 * @param pt        用于存储波峰/波谷点的幅值数组      
 * @param ipt       用于存储波峰/波谷点的个数数组         
 * @param gpt       用于存储等待迭次数的数组      
 * @param iendk0    一个布尔指针，用于指示是否满足结束条件 
 */
static void process_peak_or_trough(
    MYINT ir, MYINT m, MYINT v, MYINT maxNpt, MYINT maxnwait, 
    MYREAL k, MYREAL dk, MYCOMPLEX (*J3)[3][3][4], MYREAL (*kpt)[3][4][maxNpt], 
    MYCOMPLEX (*pt)[3][4][maxNpt], MYINT (*ipt)[3][4], MYINT (*gpt)[3][4], bool *iendk0)
{
    MYCOMPLEX tmp0;
    if (gpt[ir][m][v] >= 2 && ipt[ir][m][v] < maxNpt) {
        if (cplx_peak_or_trough(m, v, J3[ir], k, dk, &kpt[ir][m][v][ipt[ir][m][v]], &tmp0) != 0) {
            pt[ir][m][v][ipt[ir][m][v]++] = tmp0;
            gpt[ir][m][v] = 0;
        } else if (gpt[ir][m][v] >= maxnwait) {
            kpt[ir][m][v][ipt[ir][m][v]] = k - dk;
            pt[ir][m][v][ipt[ir][m][v]++] = J3[ir][1][m][v];
            gpt[ir][m][v] = 0;
        }
    }
    *iendk0 = *iendk0 && (ipt[ir][m][v] == maxNpt);
}


void PTA_method(
    const MODEL1D *mod1d, MYREAL k0, MYREAL predk, MYREAL rmin, MYREAL rmax, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs, 
    MYCOMPLEX sum_EXP_J0[nr][3][4], MYCOMPLEX sum_VF_J0[nr][3][4],  
    MYCOMPLEX sum_HF_J0[nr][3][4],  MYCOMPLEX sum_DC_J0[nr][3][4],  
    FILE *(fstats[nr]), FILE *(ptam_fstats[nr]))
{   
    // 需要兼容对正常收敛而不具有规律波峰波谷的序列
    // 有时序列收敛比较好，不表现为规律的波峰波谷，
    // 此时设置最大等待次数，超过直接设置为中间值

    MYINT ik=0;
    const MYINT maxnwait = 9;     // 最大等待次数，不能太小
    const MYREAL dk=PI/((maxnwait-1)*rmax); 
    MYREAL k=0.0;
    const MYREAL precoef = dk/predk; // 提前乘dk系数，以抵消格林函数主函数计算时最后乘dk

    MYCOMPLEX EXP_qwv[3][3], VF_qwv[3][3], HF_qwv[3][3], DC_qwv[3][3]; // 不同震源的核函数
    MYCOMPLEX (*pEXP_qwv)[3] = (sum_EXP_J0!=NULL)? EXP_qwv : NULL;
    MYCOMPLEX (*pVF_qwv)[3]  = (sum_VF_J0!=NULL)?  VF_qwv  : NULL;
    MYCOMPLEX (*pHF_qwv)[3]  = (sum_HF_J0!=NULL)?  HF_qwv  : NULL;
    MYCOMPLEX (*pDC_qwv)[3]  = (sum_DC_J0!=NULL)?  DC_qwv  : NULL;

    static const MYINT maxNpt=PTAM_MAX_PEAK_TROUGH; // 波峰波谷的目标

    // 根据波峰波谷的目标也给出一个kmax，+5以防万一 
    const MYREAL kmax = k0 + (maxNpt+5)*PI/rmin;

    // 每个震中距是否已找齐慢收敛序列
    bool iendk = true, iendk0 = false;
    bool *iendkrs = (bool *)calloc(nr, sizeof(bool));
    for(MYINT ir=0; ir<nr; ++ir) iendkrs[ir] = false;

    // 用于接收F(ki,w)Jm(ki*r)ki
    // 存储采样的值，维度3表示通过连续3个点来判断波峰或波谷
    // 既用于存储被积函数，也最后用于存储求和的结果
    MYCOMPLEX (*EXP_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*EXP_J3));
    MYCOMPLEX (*VF_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*VF_J3));
    MYCOMPLEX (*HF_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*HF_J3));
    MYCOMPLEX (*DC_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*DC_J3));

    // 之前求和的值
    MYCOMPLEX (*sum_EXP_J)[3][4] = (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_EXP_J));
    MYCOMPLEX (*sum_VF_J)[3][4] = (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_VF_J));
    MYCOMPLEX (*sum_HF_J)[3][4] = (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_HF_J));
    MYCOMPLEX (*sum_DC_J)[3][4] = (MYCOMPLEX (*)[3][4])calloc(nr, sizeof(*sum_DC_J));

    // 存储波峰波谷的位置和值
    MYREAL (*kEXPpt)[3][4][maxNpt] = (MYREAL (*)[3][4][maxNpt])calloc(nr, sizeof(*kEXPpt));
    MYREAL (*kVFpt)[3][4][maxNpt] = (MYREAL (*)[3][4][maxNpt])calloc(nr, sizeof(*kVFpt));
    MYREAL (*kHFpt)[3][4][maxNpt] = (MYREAL (*)[3][4][maxNpt])calloc(nr, sizeof(*kHFpt));
    MYREAL (*kDCpt)[3][4][maxNpt] = (MYREAL (*)[3][4][maxNpt])calloc(nr, sizeof(*kDCpt));
    MYCOMPLEX (*EXPpt)[3][4][maxNpt] = (MYCOMPLEX (*)[3][4][maxNpt])calloc(nr, sizeof(*EXPpt));
    MYCOMPLEX (*VFpt)[3][4][maxNpt] = (MYCOMPLEX (*)[3][4][maxNpt])calloc(nr, sizeof(*VFpt));
    MYCOMPLEX (*HFpt)[3][4][maxNpt] = (MYCOMPLEX (*)[3][4][maxNpt])calloc(nr, sizeof(*HFpt));
    MYCOMPLEX (*DCpt)[3][4][maxNpt] = (MYCOMPLEX (*)[3][4][maxNpt])calloc(nr, sizeof(*DCpt));
    MYINT (*iEXPpt)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*iEXPpt));
    MYINT (*iVFpt)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*iVFpt));
    MYINT (*iHFpt)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*iHFpt));
    MYINT (*iDCpt)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*iDCpt));

    // 记录点数，当峰谷找到后，清零
    MYINT (*gEXPpt)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gEXPpt));
    MYINT (*gVFpt)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gVFpt));
    MYINT (*gHFpt)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gHFpt));
    MYINT (*gDCpt)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gDCpt));
    
    
    for(MYINT ir=0; ir<nr; ++ir){
        for(MYINT m=0; m<3; ++m){
            for(MYINT v=0; v<4; ++v){
                if(sum_EXP_J0!=NULL) sum_EXP_J[ir][m][v] = sum_EXP_J0[ir][m][v];
                if(sum_VF_J0!=NULL)  sum_VF_J[ir][m][v]  = sum_VF_J0[ir][m][v];
                if(sum_HF_J0!=NULL)  sum_HF_J[ir][m][v]  = sum_HF_J0[ir][m][v];
                if(sum_DC_J0!=NULL)  sum_DC_J[ir][m][v]  = sum_DC_J0[ir][m][v];

                iEXPpt[ir][m][v] = gEXPpt[ir][m][v] = 0;
                iVFpt[ir][m][v]  = gVFpt[ir][m][v]  = 0;
                iHFpt[ir][m][v]  = gHFpt[ir][m][v]  = 0;
                iDCpt[ir][m][v]  = gDCpt[ir][m][v]  = 0;

            }
        }
        iendkrs[ir] = false;
    }


    k = k0 - dk;
    while(true){
        k += dk;
        if(k > kmax) break;

        // 计算核函数 F(k, w)
        kernel(mod1d, omega, k, pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv); 

        iendk=true;
        for(MYINT ir=0; ir<nr; ++ir){
            if(iendkrs[ir]) continue;        // 该震中距下的慢收敛序列已找齐

            // 计算被积函数一项 F(k,w)Jm(kr)k
            int_Pk(k, rs[ir],
                    pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv,
                    EXP_J3[ir][2], VF_J3[ir][2], HF_J3[ir][2], DC_J3[ir][2]);  // [2]表示把新点值放在最后
        

            // 记录积分结果
            if(fstats[ir]!=NULL){
                write_stats(
                    fstats[ir], k, 
                    pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv,
                    EXP_J3[ir][2], VF_J3[ir][2], HF_J3[ir][2], DC_J3[ir][2]);
            }

            // 赋更新量
            for(MYINT m=0; m<3; ++m){
                for(MYINT v=0; v<4; ++v){
                    // 转为求和结果
                    if(sum_EXP_J0!=NULL)  {
                        sum_EXP_J[ir][m][v] += EXP_J3[ir][2][m][v] * precoef;
                        EXP_J3[ir][2][m][v] = sum_EXP_J[ir][m][v];
                    }
                    if(sum_VF_J0!=NULL){
                        sum_VF_J[ir][m][v]  += VF_J3[ir][2][m][v] * precoef;
                        VF_J3[ir][2][m][v]  = sum_VF_J[ir][m][v];
                    }
                    if(sum_HF_J0!=NULL){
                        sum_HF_J[ir][m][v]  += HF_J3[ir][2][m][v] * precoef;
                        HF_J3[ir][2][m][v]  = sum_HF_J[ir][m][v];
                    }
                    if(sum_DC_J0!=NULL){
                        sum_DC_J[ir][m][v]  += DC_J3[ir][2][m][v] * precoef;
                        DC_J3[ir][2][m][v]  = sum_DC_J[ir][m][v];
                    }
                    
                }
            } 


            // 3点以上，判断波峰波谷 
            iendk0 = true;
            for (MYINT m = 0; m < 3; ++m) {
                for (MYINT v = 0; v < 4; ++v) {
                    if (sum_EXP_J0 != NULL && m == 0 && (v == 0 || v == 2)) {
                        process_peak_or_trough(ir, m, v, maxNpt, maxnwait, k, dk, EXP_J3, kEXPpt, EXPpt, iEXPpt, gEXPpt, &iendk0);
                    }
                    if (sum_VF_J0 != NULL && m == 0 && (v == 0 || v == 2)) {
                        process_peak_or_trough(ir, m, v, maxNpt, maxnwait, k, dk, VF_J3, kVFpt, VFpt, iVFpt, gVFpt, &iendk0);
                    }
                    if (sum_HF_J0 != NULL && m == 1) {
                        process_peak_or_trough(ir, m, v, maxNpt, maxnwait, k, dk, HF_J3, kHFpt, HFpt, iHFpt, gHFpt, &iendk0);
                    }
                    if (sum_DC_J0 != NULL && ((m == 0 && (v == 0 || v == 2)) || m != 0)) {
                        process_peak_or_trough(ir, m, v, maxNpt, maxnwait, k, dk, DC_J3, kDCpt, DCpt, iDCpt, gDCpt, &iendk0);
                    }
                }
            }
            iendkrs[ir] = iendk0;
            iendk = iendk && iendkrs[ir];

            // 左移动点, 
            for(MYINT m=0; m<3; ++m){
                for(MYINT v=0; v<4; ++v){
                    for(MYINT jj=0; jj<2; ++jj){
                        if(sum_EXP_J0!=NULL) EXP_J3[ir][jj][m][v] = EXP_J3[ir][jj+1][m][v];
                        if(sum_VF_J0!=NULL)  VF_J3[ir][jj][m][v]  = VF_J3[ir][jj+1][m][v];
                        if(sum_HF_J0!=NULL)  HF_J3[ir][jj][m][v]  = HF_J3[ir][jj+1][m][v];
                        if(sum_DC_J0!=NULL)  DC_J3[ir][jj][m][v]  = DC_J3[ir][jj+1][m][v];
                    }

                    // 未找到峰谷---点数+1
                    if(sum_EXP_J0!=NULL) gEXPpt[ir][m][v]++;
                    if(sum_VF_J0!=NULL)  gVFpt[ir][m][v]++;
                    if(sum_HF_J0!=NULL)  gHFpt[ir][m][v]++;
                    if(sum_DC_J0!=NULL)  gDCpt[ir][m][v]++;
                }
            }

        } // end rs loop

        ++ik;

        // 所有震中距的慢收敛序列都已找到
        if(iendk) break;

    } // end k loop

    // printf("w=%f, ik=%d\n", CREAL(omega), ik);


    // 做缩减序列，赋值最终解
    for(MYINT ir=0; ir<nr; ++ir){
        // 记录到文件
        if(ptam_fstats[ir]!=NULL){
            write_stats_ptam(
                ptam_fstats[ir], k, maxNpt, 
                EXPpt[ir], VFpt[ir], HFpt[ir], DCpt[ir],
                // iEXPpt[ir], iVFpt[ir], iHFpt[ir], iDCpt[ir],
                kEXPpt[ir], kVFpt[ir], kHFpt[ir], kDCpt[ir]);
        }

        for(MYINT m=0; m<3; ++m){
            for(MYINT v=0; v<4; ++v){
                if(sum_EXP_J0!=NULL)  {cplx_shrink(iEXPpt[ir][m][v], EXPpt[ir][m][v]);  sum_EXP_J0[ir][m][v] = EXPpt[ir][m][v][0];}
                if(sum_VF_J0!=NULL)   {cplx_shrink(iVFpt[ir][m][v],  VFpt[ir][m][v]);   sum_VF_J0[ir][m][v]  = VFpt[ir][m][v][0];}
                if(sum_HF_J0!=NULL)   {cplx_shrink(iHFpt[ir][m][v],  HFpt[ir][m][v]);   sum_HF_J0[ir][m][v]  = HFpt[ir][m][v][0];}
                if(sum_DC_J0!=NULL)   {cplx_shrink(iDCpt[ir][m][v],  DCpt[ir][m][v]);   sum_DC_J0[ir][m][v]  = DCpt[ir][m][v][0];}
            }
        }
    }


    free(iendkrs);
    free(EXP_J3); free(VF_J3); free(HF_J3); free(DC_J3);
    free(sum_EXP_J); free(sum_VF_J); free(sum_HF_J); free(sum_DC_J);
    free(kEXPpt); free(kVFpt); free(kHFpt); free(kDCpt);
    free(EXPpt); free(VFpt); free(HFpt); free(DCpt);
    free(iEXPpt); free(iVFpt); free(iHFpt); free(iDCpt);
    free(gEXPpt); free(gVFpt); free(gHFpt); free(gDCpt);

}




MYINT cplx_peak_or_trough(MYINT idx1, MYINT idx2, const MYCOMPLEX arr[3][3][4], MYREAL k, MYREAL dk, MYREAL *pk, MYCOMPLEX *value){
    MYCOMPLEX f1, f2, f3;
    MYREAL rf1, rf2, rf3;
    MYINT stat=0;

    f1 = arr[0][idx1][idx2];
    f2 = arr[1][idx1][idx2];
    f3 = arr[2][idx1][idx2];

    rf1 = CREAL(f1);
    rf2 = CREAL(f2);
    rf3 = CREAL(f3);
    if     ( (rf1 <= rf2) && (rf2 >= rf3) )  stat = 1;
    else if( (rf1 >= rf2) && (rf2 <= rf3) )  stat = -1;
    else                                     stat =  0;

    if(stat==0)  return stat;

    MYREAL x1, x2, x3; 
    x3 = k;
    x2 = x3-dk;
    x1 = x2-dk;

    MYREAL xarr[3] = {x1, x2, x3};
    MYCOMPLEX farr[3] = {f1, f2, f3};

    // 二次多项式
    MYCOMPLEX a, b, c;
    quad_term(xarr, farr, &a, &b, &c);

    MYREAL k0 = x2;
    *pk = k0;
    *value = 0.0;
    if(a != 0.0+0.0*I){
        k0 = - b / (2*a);

        // 拟合二次多项式可能会有各种潜在问题，例如f1,f2,f3几乎相同，此时a,b很小，k0值非常不稳定
        // 这里暂且使用范围来框定，如果在范围外，就直接使用x2的值
        if(k0 < x3 && k0 > x1){
            // printf("a=%f%+fI, b=%f%+fI, c=%f%+fI, xarr=(%f,%f,%f), yarr=(%f%+fI, %f%+fI, %f%+fI)\n", 
            //         CREAL(a),CIMAG(a),CREAL(b),CIMAG(b),CREAL(c),CIMAG(c),x1,x2,x3,CREAL(f1),CIMAG(f1),CREAL(f2),CIMAG(f2),CREAL(f3),CIMAG(f3));
            *pk = k0;
            *value = a*k0*k0 + b*k0;
        }
    } 
    *value += c;
    
    return stat;
}


void cplx_shrink(MYINT n1, MYCOMPLEX *arr){
    for(MYINT n=n1; n>1; --n){
        for(MYINT i=0; i<n-1; ++i){
            arr[i] = 0.5*(arr[i] + arr[i+1]);
        }
    }
}