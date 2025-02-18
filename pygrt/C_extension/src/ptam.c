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

#include "ptam.h"
#include "iostats.h"
#include "const.h"
#include "model.h"
#include "propagate.h"
#include "quadratic.h"






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
    bool iendkrs[nr], iendk=true, iendk0=false;
    for(MYINT ir=0; ir<nr; ++ir) iendkrs[ir] = false;

    // 用于接收F(ki,w)Jm(ki*r)ki
    // 存储采样的值，维度3表示通过连续3个点来判断波峰或波谷
    // 既用于存储被积函数，也最后用于存储求和的结果
    MYCOMPLEX EXP_J3[nr][3][3][4], VF_J3[nr][3][3][4], HF_J3[nr][3][3][4],  DC_J3[nr][3][3][4];

    // 之前求和的值
    MYCOMPLEX sum_EXP_J[nr][3][4], sum_VF_J[nr][3][4];
    MYCOMPLEX sum_HF_J[nr][3][4],  sum_DC_J[nr][3][4];

    // 存储波峰波谷的位置和值
    MYREAL kEXPpt[nr][3][4][maxNpt], kVFpt[nr][3][4][maxNpt];
    MYREAL kHFpt[nr][3][4][maxNpt],  kDCpt[nr][3][4][maxNpt];
    MYCOMPLEX EXPpt[nr][3][4][maxNpt], VFpt[nr][3][4][maxNpt];
    MYCOMPLEX HFpt[nr][3][4][maxNpt],  DCpt[nr][3][4][maxNpt];
    MYINT iEXPpt[nr][3][4], iVFpt[nr][3][4], iHFpt[nr][3][4], iDCpt[nr][3][4];

    // 记录点数，当峰谷找到后，清零
    MYINT gEXPpt[nr][3][4], gVFpt[nr][3][4], gHFpt[nr][3][4], gDCpt[nr][3][4];
    
    
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
            MYCOMPLEX tmp0;
            iendk0 = true;
            for(MYINT m=0; m<3; ++m){
                for(MYINT v=0; v<4; ++v){
                    if(sum_EXP_J0!=NULL && m==0 && (v==0||v==2)){
                        if(gEXPpt[ir][m][v] >= 2 && iEXPpt[ir][m][v] < maxNpt){
                            if(cplx_peak_or_trough(m, v, EXP_J3[ir], k, dk, &kEXPpt[ir][m][v][iEXPpt[ir][m][v]], &tmp0) != 0){
                                // 成功找到峰谷
                                EXPpt[ir][m][v][iEXPpt[ir][m][v]++] = tmp0;
                                gEXPpt[ir][m][v] = 0;
                            }
                            else if(gEXPpt[ir][m][v] >= maxnwait){
                                // 等待过多，直接取中间值
                                kEXPpt[ir][m][v][iEXPpt[ir][m][v]] = k-dk;
                                EXPpt[ir][m][v][iEXPpt[ir][m][v]++] = EXP_J3[ir][1][m][v];
                                gEXPpt[ir][m][v] = 0;
                            }
                        }
                        iendk0 = iendk0 && (iEXPpt[ir][m][v]==maxNpt);
                    }
                    if(sum_VF_J0!=NULL && m==0 && (v==0||v==2) ){
                        if(gVFpt[ir][m][v] >= 2 && iVFpt[ir][m][v] < maxNpt){ 
                            if(cplx_peak_or_trough(m, v, VF_J3[ir], k, dk, &kVFpt[ir][m][v][iVFpt[ir][m][v]], &tmp0) != 0){
                                VFpt[ir][m][v][iVFpt[ir][m][v]++] = tmp0;
                                gVFpt[ir][m][v] = 0;
                            }
                            else if(gVFpt[ir][m][v] >= maxnwait){
                                kVFpt[ir][m][v][iVFpt[ir][m][v]] = k-dk;
                                VFpt[ir][m][v][iVFpt[ir][m][v]++] = VF_J3[ir][1][m][v];
                                gVFpt[ir][m][v] = 0;
                            }
                        }
                        iendk0 = iendk0 && (iVFpt[ir][m][v]==maxNpt);
                    }
                    if(sum_HF_J0!=NULL && m==1){
                        if(gHFpt[ir][m][v] >= 2 && iHFpt[ir][m][v] < maxNpt){ 
                            if(cplx_peak_or_trough(m, v, HF_J3[ir], k, dk, &kHFpt[ir][m][v][iHFpt[ir][m][v]], &tmp0) != 0){
                                HFpt[ir][m][v][iHFpt[ir][m][v]++] = tmp0;
                                gHFpt[ir][m][v] = 0;
                            }
                            else if(gHFpt[ir][m][v] >= maxnwait){
                                kHFpt[ir][m][v][iHFpt[ir][m][v]] = k-dk;
                                HFpt[ir][m][v][iHFpt[ir][m][v]++] = HF_J3[ir][1][m][v];
                                gHFpt[ir][m][v] = 0;
                            }
                        }
                        iendk0 = iendk0 && (iHFpt[ir][m][v]==maxNpt);
                    }
                    if(sum_DC_J0!=NULL && ((m==0 && (v==0||v==2)) || m!=0)){
                        if(gDCpt[ir][m][v] >= 2 && iDCpt[ir][m][v] < maxNpt){ 
                            if(cplx_peak_or_trough(m, v, DC_J3[ir], k, dk, &kDCpt[ir][m][v][iDCpt[ir][m][v]], &tmp0) != 0){
                                DCpt[ir][m][v][iDCpt[ir][m][v]++] = tmp0;
                                gDCpt[ir][m][v] = 0;
                            }
                            else if(gDCpt[ir][m][v] >= maxnwait){
                                kDCpt[ir][m][v][iDCpt[ir][m][v]] = k-dk;
                                DCpt[ir][m][v][iDCpt[ir][m][v]++] = DC_J3[ir][1][m][v];
                                gDCpt[ir][m][v] = 0;
                            }
                        }
                        iendk0 = iendk0 && (iDCpt[ir][m][v]==maxNpt);
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