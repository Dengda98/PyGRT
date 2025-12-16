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
#include <string.h>

#include "grt/integral/ptam.h"
#include "grt/integral/quadratic.h"
#include "grt/integral/k_integ.h"
#include "grt/integral/iostats.h"
#include "grt/common/const.h"
#include "grt/common/model.h"




/**
 * 观察连续3个点的函数值的实部变化，判断是波峰(1)还是波谷(-1), 并计算对应值。
 * 
 * @param[in]     idx1        阶数索引
 * @param[in]     idx2        积分类型索引 
 * @param[in]     arr         存有连续三个点的函数值的数组 
 * @param[in]     k           三个点的起始波数
 * @param[in]     dk          三个点的波数间隔，这样使用k和dk定义了三个点的位置
 * @param[out]    pk          估计的波峰或波谷处的波数
 * @param[out]    value       估计的波峰或波谷处的函数值
 * 
 * @return    波峰(1)，波谷(-1)，其它(0)
 *  
 */
static int _cplx_peak_or_trough(
    int idx1, int idx2, const cplxIntegGrid arr[GRT_PTAM_WINDOW_SIZE], 
    real_t k, real_t dk, real_t *pk, cplx_t *value)
{
    cplx_t f1, f2, f3;
    real_t rf1, rf2, rf3;
    int stat=0;

    f1 = arr[0][idx1][idx2];
    f2 = arr[1][idx1][idx2];
    f3 = arr[2][idx1][idx2];

    rf1 = creal(f1);
    rf2 = creal(f2);
    rf3 = creal(f3);
    if     ( (rf1 <= rf2) && (rf2 >= rf3) )  stat = 1;
    else if( (rf1 >= rf2) && (rf2 <= rf3) )  stat = -1;
    else                                     stat =  0;

    if(stat==0)  return stat;

    real_t x1, x2, x3; 
    x3 = k;
    x2 = x3-dk;
    x1 = x2-dk;

    real_t xarr[3] = {x1, x2, x3};
    cplx_t farr[3] = {f1, f2, f3};

    // 二次多项式
    cplx_t a, b, c;
    grt_quad_term(xarr, farr, &a, &b, &c);

    real_t k0 = x2;
    *pk = k0;
    *value = 0.0;
    if(a != 0.0+0.0*I){
        k0 = - b / (2*a);

        // 拟合二次多项式可能会有各种潜在问题，例如f1,f2,f3几乎相同，此时a,b很小，k0值非常不稳定
        // 这里暂且使用范围来框定，如果在范围外，就直接使用x2的值
        if(k0 < x3 && k0 > x1){
            // printf("a=%f%+fI, b=%f%+fI, c=%f%+fI, xarr=(%f,%f,%f), yarr=(%f%+fI, %f%+fI, %f%+fI)\n", 
            //         creal(a),cimag(a),creal(b),cimag(b),creal(c),cimag(c),x1,x2,x3,creal(f1),cimag(f1),creal(f2),cimag(f2),creal(f3),cimag(f3));
            *pk = k0;
            *value = a*k0*k0 + b*k0;
        }
    } 
    *value += c;
    
    return stat;
}


/**
 * 处理并确定波峰或波谷                                    
 * 
 * @param[in]           ir        震中距索引                          
 * @param[in]           im        不同震源不同阶数的索引              
 * @param[in]           v         积分形式索引                          
 * @param[in]           k         波数                             
 * @param[in]           dk        波数步长                              
 * @param[in]           J3        存储的积分采样幅值数组                  
 * @param[in,out]       Kpt       积分值峰谷的波数数组     
 * @param[in,out]       Fpt       用于存储波峰/波谷点的幅值数组 
 * @param[in,out]       Ipt       用于存储波峰/波谷点的个数数组 
 * @param[in,out]       Gpt       用于存储等待迭次数的数组    
 * @param[in,out]       iendk0    一个布尔指针，用于指示是否满足结束条件 
 */
static void process_peak_or_trough(
    size_t ir, int im, int v, real_t k, real_t dk, 
    cplxIntegGrid (*J3)[GRT_PTAM_WINDOW_SIZE], realIntegGrid (*Kpt)[GRT_PTAM_PT_MAX], 
    cplxIntegGrid (*Fpt)[GRT_PTAM_PT_MAX], sizeIntegGrid (*Ipt), sizeIntegGrid (*Gpt), bool *iendk0)
{
    cplx_t tmp0;
    if (Gpt[ir][im][v] >= GRT_PTAM_WINDOW_SIZE-1 && Ipt[ir][im][v] < GRT_PTAM_PT_MAX) {
        if (_cplx_peak_or_trough(im, v, J3[ir], k, dk, &Kpt[ir][Ipt[ir][im][v]][im][v], &tmp0) != 0) {
            Fpt[ir][Ipt[ir][im][v]++][im][v] = tmp0;
            Gpt[ir][im][v] = 0;
        } else if (Gpt[ir][im][v] >= GRT_PTAM_WAITS_MAX) {  // 不再等待，直接取中点作为波峰波谷
            Kpt[ir][Ipt[ir][im][v]][im][v] = k - dk;
            Fpt[ir][Ipt[ir][im][v]++][im][v] = J3[ir][1][im][v];
            Gpt[ir][im][v] = 0;
        }
    }
    *iendk0 = *iendk0 && (Ipt[ir][im][v] == GRT_PTAM_PT_MAX);
}


/**
 * 在输入被积函数的情况下，对不同震源使用峰谷平均法
 * 
 * @param[in]           ir                  震中距索引
 * @param[in]           nr                  震中距个数
 * @param[in]           precoef             积分值系数
 * @param[in]           k                   波数                             
 * @param[in]           dk                  波数步长       
 * @param[in,out]       SUM3                被积函数的幅值数组 
 * @param[in,out]       sumJ                积分值数组 
 * 
 * @param[in,out]       Kpt                 积分值峰谷的波数数组     
 * @param[in,out]       Fpt                 用于存储波峰/波谷点的幅值数组 
 * @param[in,out]       Ipt                 用于存储波峰/波谷点的个数数组 
 * @param[in,out]       Gpt                 用于存储等待迭次数的数组 
 * 
 * @param[in,out]       iendk0              是否收集足够峰谷
 * 
 */
static void ptam_once(
    const size_t ir, const size_t nr, const real_t precoef, real_t k, real_t dk, 
    cplxIntegGrid SUM3[nr][GRT_PTAM_WINDOW_SIZE],
    cplxIntegGrid sumJ[nr],
    realIntegGrid Kpt[nr][GRT_PTAM_PT_MAX],
    cplxIntegGrid Fpt[nr][GRT_PTAM_PT_MAX],
    sizeIntegGrid Ipt[nr],
    sizeIntegGrid Gpt[nr],
    bool *iendk0)
{
    *iendk0 = true;

    GRT_LOOP_IntegGrid(im, v){
        int modr = GRT_SRC_M_ORDERS[im];
        if(modr == 0 && v!=0 && v!= 2)  continue;

        // 赋更新量
        // SUM3转为求和结果
        sumJ[ir][im][v] += SUM3[ir][GRT_PTAM_WINDOW_SIZE-1][im][v] * precoef;
        SUM3[ir][GRT_PTAM_WINDOW_SIZE-1][im][v] = sumJ[ir][im][v];         
        
        // 3点以上，判断波峰波谷 
        process_peak_or_trough(ir, im, v, k, dk, SUM3, Kpt, Fpt, Ipt, Gpt, iendk0);

        // 左移动点, 
        for(int jj=0; jj<GRT_PTAM_WINDOW_SIZE-1; ++jj){
            SUM3[ir][jj][im][v] = SUM3[ir][jj+1][im][v];
        }

        // 点数+1
        Gpt[ir][im][v]++;
    }
}




/**
 * 递归式地计算缩减序列的值，
 * \f[
 * M_i = 0.5\times (M_i + M_{i+1})
 * \f]
 * 
 * @param[in]         n1          数组长度 
 * @param[in]         ir          震中距索引                          
 * @param[in]         im          不同震源不同阶数的索引              
 * @param[in]         v           积分形式索引  
 * @param[in,out]     Fpt         用于存储波峰/波谷点的幅值数组，最终收敛值在第一个 
 * 
 */
static void _cplx_shrink(size_t n1, size_t ir,  int im, int v, cplxIntegGrid (*Fpt)[GRT_PTAM_PT_MAX]){
    for(size_t n=n1; n>1; --n){
        for(size_t i=0; i<n-1; ++i){
            Fpt[ir][i][im][v] = 0.5*(Fpt[ir][i][im][v] + Fpt[ir][i+1][im][v]);
        }
    }
}



void grt_PTA_method(
    GRT_MODEL1D *mod1d, real_t k0, real_t predk,
    size_t nr, real_t *rs, K_INTEG *K, FILE *ptam_fstatsnr[nr][2], GRT_KernelFunc kerfunc)
{   
    // 需要兼容对正常收敛而不具有规律波峰波谷的序列
    // 有时序列收敛比较好，不表现为规律的波峰波谷，
    // 此时设置最大等待次数，超过直接设置为中间值

    real_t k=0.0;

    // 使用宏函数，方便定义
    #define __CALLOC_ARRAY(VAR, TYP, __ARR) \
        TYP (*VAR)__ARR = (TYP (*)__ARR)calloc(nr, sizeof(*VAR));

    // 用于接收F(ki,w)Jm(ki*r)ki
    // 存储采样的值，维度3表示通过连续3个点来判断波峰或波谷
    // 既用于存储被积函数，也最后用于存储求和的结果
    #define __ARR [GRT_PTAM_WINDOW_SIZE]
        __CALLOC_ARRAY(SUM3, cplxIntegGrid, __ARR);
        __CALLOC_ARRAY(SUM3_uiz, cplxIntegGrid, __ARR);
        __CALLOC_ARRAY(SUM3_uir, cplxIntegGrid, __ARR);
    #undef __ARR

    // 存储波峰波谷的位置和值
    #define __ARR [GRT_PTAM_PT_MAX]
        __CALLOC_ARRAY(Kpt, realIntegGrid, __ARR);
        __CALLOC_ARRAY(Fpt, cplxIntegGrid, __ARR);

        __CALLOC_ARRAY(Kpt_uiz, realIntegGrid, __ARR);
        __CALLOC_ARRAY(Fpt_uiz, cplxIntegGrid, __ARR);

        __CALLOC_ARRAY(Kpt_uir, realIntegGrid, __ARR);
        __CALLOC_ARRAY(Fpt_uir, cplxIntegGrid, __ARR);
    #undef __ARR

    #define __ARR
        // 存储波峰波谷的总个数
        __CALLOC_ARRAY(Ipt,     sizeIntegGrid, __ARR);
        __CALLOC_ARRAY(Ipt_uiz, sizeIntegGrid, __ARR);
        __CALLOC_ARRAY(Ipt_uir, sizeIntegGrid, __ARR);

        // 记录点数，当峰谷找到后，清零
        __CALLOC_ARRAY(Gpt,     sizeIntegGrid, __ARR);
        __CALLOC_ARRAY(Gpt_uiz, sizeIntegGrid, __ARR);
        __CALLOC_ARRAY(Gpt_uir, sizeIntegGrid, __ARR);
    #undef __ARR
    #undef __CALLOC_ARRAY


    // 对于PTAM，不同震中距使用不同dk
    for(size_t ir=0; ir<nr; ++ir){
        real_t dk = PI/((GRT_PTAM_WAITS_MAX-1)*rs[ir]); 
        real_t precoef = dk/predk; // 提前乘dk系数，以抵消格林函数主函数计算时最后乘dk
        // 根据波峰波谷的目标也给出一个kmax，+5以防万一 
        real_t kmax = k0 + (GRT_PTAM_PT_MAX+5)*PI/rs[ir];

        bool iendk0=false;

        k = k0;
        while(true){
            if(k > kmax) break;
            k += dk;

            // 计算核函数 F(k, w)
            kerfunc(mod1d, k, K->QWV, K->calc_upar, K->QWVz); 
            if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

            // 记录核函数
            if(ptam_fstatsnr != NULL)  grt_write_stats(ptam_fstatsnr[ir][0], k, (K->calc_upar)? K->QWVz : K->QWV);

            // 计算被积函数一项 F(k,w)Jm(kr)k
            grt_int_Pk(k, rs[ir], K->QWV, false, SUM3[ir][GRT_PTAM_WINDOW_SIZE-1]);  // [GRT_PTAM_WINDOW_SIZE-1]表示把新点值放在最后
            // 判断和记录波峰波谷
            ptam_once(ir, nr, precoef, k, dk, SUM3, K->sumJ, Kpt, Fpt, Ipt, Gpt, &iendk0);
            
            // -------------------------- 位移空间导数 ------------------------------------
            if(K->calc_upar){
                // ------------------------------- ui_z -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                grt_int_Pk(k, rs[ir], K->QWVz, false, SUM3_uiz[ir][GRT_PTAM_WINDOW_SIZE-1]);  // [GRT_PTAM_WINDOW_SIZE-1]表示把新点值放在最后
                // 判断和记录波峰波谷
                ptam_once(ir, nr, precoef, k, dk, SUM3_uiz, K->sumJz, Kpt_uiz, Fpt_uiz, Ipt_uiz, Gpt_uiz, &iendk0);

                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                grt_int_Pk(k, rs[ir], K->QWV, true, SUM3_uir[ir][GRT_PTAM_WINDOW_SIZE-1]);  // [GRT_PTAM_WINDOW_SIZE-1]表示把新点值放在最后
                // 判断和记录波峰波谷
                ptam_once(ir, nr, precoef, k, dk, SUM3_uir, K->sumJr, Kpt_uir, Fpt_uir, Ipt_uir, Gpt_uir, &iendk0);
            
            } // END if calc_upar


            if(iendk0) break;
        }// end k loop
    }

    // 做缩减序列，赋值最终解
    for(size_t ir=0; ir<nr; ++ir){
        // 记录到文件
        if(ptam_fstatsnr != NULL)  grt_write_stats_ptam(ptam_fstatsnr[ir][1], Kpt[ir], (K->calc_upar)? Fpt_uiz[ir] : Fpt[ir]);

        GRT_LOOP_IntegGrid(im, v){
            _cplx_shrink(Ipt[ir][im][v], ir, im, v, Fpt);  
            K->sumJ[ir][im][v] = Fpt[ir][0][im][v];

            if(K->calc_upar){
                _cplx_shrink(Ipt_uiz[ir][im][v], ir, im, v, Fpt_uiz);  
                K->sumJz[ir][im][v] = Fpt_uiz[ir][0][im][v];
            
                _cplx_shrink(Ipt_uir[ir][im][v], ir, im, v, Fpt_uir);  
                K->sumJr[ir][im][v] = Fpt_uir[ir][0][im][v];
            }
        }
    }

    BEFORE_RETURN:

    #define __FREE_ALL_ARRAY \
        X(SUM3)    X(SUM3_uiz)    X(SUM3_uir) \
        X(Kpt)     X(Fpt)         X(Ipt)       X(Gpt) \
        X(Kpt_uiz)     X(Fpt_uiz)         X(Ipt_uiz)       X(Gpt_uiz) \
        X(Kpt_uir)     X(Fpt_uir)         X(Ipt_uir)       X(Gpt_uir) \

    #define X(A)  GRT_SAFE_FREE_PTR(A);
        __FREE_ALL_ARRAY
    #undef X
    #undef __FREE_ALL_ARRAY
}