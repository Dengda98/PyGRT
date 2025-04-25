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

#include "common/ptam.h"
#include "common/quadratic.h"
#include "common/integral.h"
#include "common/iostats.h"
#include "common/const.h"
#include "common/model.h"


/**
 * 处理并确定波峰或波谷                                    
 * 
 * @param ir        震中距索引                          
 * @param im        不同震源不同阶数的索引              
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
    MYINT ir, MYINT im, MYINT v, MYINT maxNpt, MYINT maxnwait, 
    MYREAL k, MYREAL dk, MYCOMPLEX (*J3)[3][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS], MYREAL (*kpt)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt], 
    MYCOMPLEX (*pt)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt], MYINT (*ipt)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS], MYINT (*gpt)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS], bool *iendk0)
{
    MYCOMPLEX tmp0;
    if (gpt[ir][im][v] >= 2 && ipt[ir][im][v] < maxNpt) {
        if (cplx_peak_or_trough(im, v, J3[ir], k, dk, &kpt[ir][im][v][ipt[ir][im][v]], &tmp0) != 0) {
            pt[ir][im][v][ipt[ir][im][v]++] = tmp0;
            gpt[ir][im][v] = 0;
        } else if (gpt[ir][im][v] >= maxnwait) {
            kpt[ir][im][v][ipt[ir][im][v]] = k - dk;
            pt[ir][im][v][ipt[ir][im][v]++] = J3[ir][1][im][v];
            gpt[ir][im][v] = 0;
        }
    }
    *iendk0 = *iendk0 && (ipt[ir][im][v] == maxNpt);
}


/**
 * 在输入被积函数的情况下，对不同震源使用峰谷平均法
 * 
 * @param       ir                  震中距索引
 * @param       nr                  震中距个数
 * @param       precoef             积分值系数
 * @param       maxNpt              最大峰谷数  
 * @param       maxnwait            最大等待次数      
 * @param       k                   波数                             
 * @param       dk                  波数步长       
 * @param       SUM3                被积函数的幅值数组 
 * @param       sum_J               的积分值数组 
 * 
 * @param       Kpt                 积分值峰谷的波数数组     
 * @param       Fpt                 用于存储波峰/波谷点的幅值数组 
 * @param       Ipt                 用于存储波峰/波谷点的个数数组 
 * @param       Gpt                 用于存储等待迭次数的数组 
 * 
 * @param       iendk0              是否收集足够峰谷
 * 
 */
static void ptam_once(
    const MYINT ir, const MYINT nr, const MYREAL precoef, 
    MYINT maxNpt, MYINT maxnwait, MYREAL k, MYREAL dk, 
    MYCOMPLEX SUM3[nr][3][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS],
    MYCOMPLEX sum_J[nr][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS],
    MYREAL Kpt[nr][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt],
    MYCOMPLEX Fpt[nr][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt],
    MYINT Ipt[nr][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS],
    MYINT Gpt[nr][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS],
    bool *iendk0)
{
    *iendk0 = true;
    for(MYINT i=0; i<GRT_SRC_M_COUNTS; ++i){
        MYINT modr = GRT_SRC_M_ORDERS[i];
        for(MYINT v=0; v<GRT_SRC_P_COUNTS; ++v){
            if(modr == 0 && v!=0 && v!= 2)  continue;

            // 赋更新量
            // SUM3转为求和结果
            sum_J[ir][i][v] += SUM3[ir][2][i][v] * precoef;
            SUM3[ir][2][i][v] = sum_J[ir][i][v];         
            
            // 3点以上，判断波峰波谷 
            process_peak_or_trough(ir, i, v, maxNpt, maxnwait, k, dk, SUM3, Kpt, Fpt, Ipt, Gpt, iendk0);

            // 左移动点, 
            for(MYINT jj=0; jj<2; ++jj){
                SUM3[ir][jj][i][v] = SUM3[ir][jj+1][i][v];
            }

            // 点数+1
            Gpt[ir][i][v]++;
        }
    } 
}


void PTA_method(
    const MODEL1D *mod1d, MYREAL k0, MYREAL predk, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_J0[nr][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS],
    bool calc_upar,
    MYCOMPLEX sum_uiz_J0[nr][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS],
    MYCOMPLEX sum_uir_J0[nr][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS],
    FILE *ptam_fstatsnr[nr][2], KernelFunc kerfunc)
{   
    // 需要兼容对正常收敛而不具有规律波峰波谷的序列
    // 有时序列收敛比较好，不表现为规律的波峰波谷，
    // 此时设置最大等待次数，超过直接设置为中间值

    const MYINT maxnwait = 9;     // 最大等待次数，不能太小
    MYREAL k=0.0;

    // 不同震源不同阶数的核函数 F(k, w) 
    MYCOMPLEX QWV[GRT_SRC_M_COUNTS][GRT_SRC_QWV_COUNTS];
    MYCOMPLEX QWV_uiz[GRT_SRC_M_COUNTS][GRT_SRC_QWV_COUNTS];

    static const MYINT maxNpt=PTAM_MAX_PEAK_TROUGH; // 波峰波谷的目标


    // 用于接收F(ki,w)Jm(ki*r)ki
    // 存储采样的值，维度3表示通过连续3个点来判断波峰或波谷
    // 既用于存储被积函数，也最后用于存储求和的结果
    MYCOMPLEX (*SUM3)[3][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS] = (MYCOMPLEX (*)[3][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS])calloc(nr, sizeof(*SUM3));
    MYCOMPLEX (*SUM3_uiz)[3][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS] = (MYCOMPLEX (*)[3][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS])calloc(nr, sizeof(*SUM3_uiz));
    MYCOMPLEX (*SUM3_uir)[3][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS] = (MYCOMPLEX (*)[3][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS])calloc(nr, sizeof(*SUM3_uir));

    // 之前求和的值
    MYCOMPLEX (*sum_J)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS] = (MYCOMPLEX(*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS])calloc(nr, sizeof(*sum_J));
    MYCOMPLEX (*sum_uiz_J)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS] =  (MYCOMPLEX(*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS])calloc(nr, sizeof(*sum_uiz_J));
    MYCOMPLEX (*sum_uir_J)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS] =  (MYCOMPLEX(*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS])calloc(nr, sizeof(*sum_uir_J));

    // 存储波峰波谷的位置和值
    MYREAL (*Kpt)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt] = (MYREAL (*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt])calloc(nr, sizeof(*Kpt));
    MYCOMPLEX (*Fpt)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt] = (MYCOMPLEX (*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt])calloc(nr, sizeof(*Fpt));
    MYINT (*Ipt)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS] = (MYINT (*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS])calloc(nr, sizeof(*Ipt));

    MYREAL (*Kpt_uiz)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt] = (MYREAL (*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt])calloc(nr, sizeof(*Kpt_uiz));
    MYCOMPLEX (*Fpt_uiz)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt] = (MYCOMPLEX (*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt])calloc(nr, sizeof(*Fpt_uiz));
    MYINT (*Ipt_uiz)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS] = (MYINT (*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS])calloc(nr, sizeof(*Ipt_uiz));

    MYREAL (*Kpt_uir)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt] = (MYREAL (*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt])calloc(nr, sizeof(*Kpt_uir));
    MYCOMPLEX (*Fpt_uir)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt] = (MYCOMPLEX (*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt])calloc(nr, sizeof(*Fpt_uir));
    MYINT (*Ipt_uir)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS] = (MYINT (*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS])calloc(nr, sizeof(*Ipt_uir));

    // 记录点数，当峰谷找到后，清零
    MYINT (*Gpt)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS] = (MYINT (*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS])calloc(nr, sizeof(*Gpt));
    MYINT (*Gpt_uiz)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS] = (MYINT (*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS])calloc(nr, sizeof(*Gpt_uiz));
    MYINT (*Gpt_uir)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS] = (MYINT (*)[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS])calloc(nr, sizeof(*Gpt_uir));

    for(MYINT ir=0; ir<nr; ++ir){
        for(MYINT i=0; i<GRT_SRC_M_COUNTS; ++i){
            for(MYINT v=0; v<GRT_SRC_P_COUNTS; ++v){
                sum_J[ir][i][v] = sum_J0[ir][i][v];

                if(calc_upar){
                    sum_uiz_J[ir][i][v] = sum_uiz_J0[ir][i][v];
                    sum_uir_J[ir][i][v] = sum_uir_J0[ir][i][v];
                }

                Ipt[ir][i][v] = Gpt[ir][i][v] = 0;
                Ipt_uiz[ir][i][v] = Gpt_uiz[ir][i][v] = 0;
                Ipt_uir[ir][i][v] = Gpt_uir[ir][i][v] = 0;
            }
        }
    }


    // 对于PTAM，不同震中距使用不同dk
    for(MYINT ir=0; ir<nr; ++ir){
        MYREAL dk = PI/((maxnwait-1)*rs[ir]); 
        MYREAL precoef = dk/predk; // 提前乘dk系数，以抵消格林函数主函数计算时最后乘dk
        // 根据波峰波谷的目标也给出一个kmax，+5以防万一 
        MYREAL kmax = k0 + (maxNpt+5)*PI/rs[ir];

        bool iendk0=false;

        // 积分过程文件
        FILE *fstatsK = ptam_fstatsnr[ir][0];

        k = k0;
        while(true){
            if(k > kmax) break;
            k += dk;

            // 计算核函数 F(k, w)
            kerfunc(mod1d, omega, k, QWV, calc_upar, QWV_uiz); 

            // 记录核函数
            if(fstatsK!=NULL)  write_stats(fstatsK, k, QWV);

            // 计算被积函数一项 F(k,w)Jm(kr)k
            int_Pk(k, rs[ir], QWV, false, SUM3[ir][2]);  // [2]表示把新点值放在最后

            // 
            ptam_once(
                ir, nr, precoef, maxNpt, maxnwait, k, dk, 
                SUM3, sum_J, Kpt, Fpt, Ipt, Gpt, &iendk0
            );
            
            // -------------------------- 位移空间导数 ------------------------------------
            if(calc_upar){
                // ------------------------------- ui_z -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                int_Pk(k, rs[ir], QWV_uiz, false, SUM3_uiz[ir][2]);  // [2]表示把新点值放在最后
                
                ptam_once(
                    ir, nr, precoef, maxNpt, maxnwait, k, dk, 
                    SUM3_uiz, sum_uiz_J, Kpt_uiz, Fpt_uiz, Ipt_uiz, Gpt_uiz, &iendk0
                );

                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                int_Pk(k, rs[ir], QWV, true, SUM3_uir[ir][2]);  // [2]表示把新点值放在最后
                
                ptam_once(
                    ir, nr, precoef, maxNpt, maxnwait, k, dk, 
                    SUM3_uir, sum_uir_J, Kpt_uir, Fpt_uir, Ipt_uir, Gpt_uir, &iendk0
                );
            
            } // END if calc_upar


            if(iendk0) break;
        }// end k loop
    }

    // printf("w=%f, ik=%d\n", CREAL(omega), ik);


    // 做缩减序列，赋值最终解
    for(MYINT ir=0; ir<nr; ++ir){
        FILE *fstatsP = ptam_fstatsnr[ir][1];
        // 记录到文件
        if(fstatsP!=NULL)  write_stats_ptam(fstatsP, k, maxNpt, Kpt[ir], Fpt[ir]);

        for(MYINT i=0; i<GRT_SRC_M_COUNTS; ++i){
            for(MYINT v=0; v<GRT_SRC_P_COUNTS; ++v){
                cplx_shrink(Ipt[ir][i][v], Fpt[ir][i][v]);  
                sum_J0[ir][i][v] = Fpt[ir][i][v][0];

                if(calc_upar){
                    cplx_shrink(Ipt_uiz[ir][i][v], Fpt_uiz[ir][i][v]);  
                    sum_uiz_J0[ir][i][v] = Fpt_uiz[ir][i][v][0];
                
                    cplx_shrink(Ipt_uir[ir][i][v], Fpt_uir[ir][i][v]);  
                    sum_uir_J0[ir][i][v] = Fpt_uir[ir][i][v][0];
                }
            }
        }
    }


    free(SUM3); free(SUM3_uiz); free(SUM3_uir);
    free(sum_J); free(sum_uiz_J); free(sum_uir_J); 

    free(Kpt);  free(Fpt);  free(Ipt);  free(Gpt);
    free(Kpt_uiz);  free(Fpt_uiz);  free(Ipt_uiz);  free(Gpt_uiz);
    free(Kpt_uir);  free(Fpt_uir);  free(Ipt_uir);  free(Gpt_uir);

}




MYINT cplx_peak_or_trough(MYINT idx1, MYINT idx2, const MYCOMPLEX arr[3][GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS], MYREAL k, MYREAL dk, MYREAL *pk, MYCOMPLEX *value){
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