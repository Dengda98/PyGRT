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
 * @param       EXP_J3              爆炸源对应的被积函数的幅值数组，下同
 * @param       VF_J3               垂直力源
 * @param       HF_J3               水平力源
 * @param       DC_J3               剪切源
 * @param       sum_EXP_J           爆炸源对应的积分值数组，下同
 * @param       sum_VF_J            垂直力源
 * @param       sum_HF_J            水平力源
 * @param       sum_DC_J            剪切源
 * 
 * @param       kEXPpt              爆炸源对应的积分值峰谷的波数数组，下同             
 * @param       EXPpt               爆炸源对应的用于存储波峰/波谷点的幅值数组，下同      
 * @param       iEXPpt              爆炸源对应的用于存储波峰/波谷点的个数数组，下同         
 * @param       gEXPpt              爆炸源对应的用于存储等待迭次数的数组，下同
 * @param       kVFpt               垂直力源
 * @param       VFpt                垂直力源
 * @param       iVFpt               垂直力源
 * @param       gVFpt               垂直力源
 * @param       kHFpt               水平力源
 * @param       HFpt                水平力源
 * @param       iHFpt               水平力源
 * @param       gHFpt               水平力源
 * @param       kDCpt               剪切源
 * @param       DCpt                剪切源
 * @param       iDCpt               剪切源
 * @param       gDCpt               剪切源
 * 
 * 
 */
static void ptam_once(
    const MYINT ir, const MYINT nr, const MYREAL precoef, 
    MYINT maxNpt, MYINT maxnwait, MYREAL k, MYREAL dk, 
    MYCOMPLEX EXP_J3[nr][3][3][4], MYCOMPLEX VF_J3[nr][3][3][4], 
    MYCOMPLEX HF_J3[nr][3][3][4], MYCOMPLEX DC_J3[nr][3][3][4], 
    MYCOMPLEX sum_EXP_J[nr][3][4], MYCOMPLEX sum_VF_J[nr][3][4],  
    MYCOMPLEX sum_HF_J[nr][3][4],  MYCOMPLEX sum_DC_J[nr][3][4],  
    MYREAL kEXPpt[nr][3][4][maxNpt], MYCOMPLEX EXPpt[nr][3][4][maxNpt], MYINT iEXPpt[nr][3][4], MYINT gEXPpt[nr][3][4],
    MYREAL kVFpt[nr][3][4][maxNpt], MYCOMPLEX VFpt[nr][3][4][maxNpt], MYINT iVFpt[nr][3][4], MYINT gVFpt[nr][3][4],
    MYREAL kHFpt[nr][3][4][maxNpt], MYCOMPLEX HFpt[nr][3][4][maxNpt], MYINT iHFpt[nr][3][4], MYINT gHFpt[nr][3][4],
    MYREAL kDCpt[nr][3][4][maxNpt], MYCOMPLEX DCpt[nr][3][4][maxNpt], MYINT iDCpt[nr][3][4], MYINT gDCpt[nr][3][4],
    bool *iendk0)
{
    // 赋更新量
    for(MYINT m=0; m<3; ++m){
        for(MYINT v=0; v<4; ++v){
            // EXP_J3, VF_J3, HF_J3, DC_J3转为求和结果
            if(sum_EXP_J!=NULL)  {
                sum_EXP_J[ir][m][v] += EXP_J3[ir][2][m][v] * precoef;
                EXP_J3[ir][2][m][v] = sum_EXP_J[ir][m][v];
            }
            if(sum_VF_J!=NULL){
                sum_VF_J[ir][m][v]  += VF_J3[ir][2][m][v] * precoef;
                VF_J3[ir][2][m][v]  = sum_VF_J[ir][m][v];
            }
            if(sum_HF_J!=NULL){
                sum_HF_J[ir][m][v]  += HF_J3[ir][2][m][v] * precoef;
                HF_J3[ir][2][m][v]  = sum_HF_J[ir][m][v];
            }
            if(sum_DC_J!=NULL){
                sum_DC_J[ir][m][v]  += DC_J3[ir][2][m][v] * precoef;
                DC_J3[ir][2][m][v]  = sum_DC_J[ir][m][v];
            }
            
        }
    } 

    // 3点以上，判断波峰波谷 
    *iendk0 = true;
    for (MYINT m = 0; m < 3; ++m) {
        for (MYINT v = 0; v < 4; ++v) {
            if (sum_EXP_J != NULL && m == 0 && (v == 0 || v == 2)) {
                process_peak_or_trough(ir, m, v, maxNpt, maxnwait, k, dk, EXP_J3, kEXPpt, EXPpt, iEXPpt, gEXPpt, iendk0);
            }
            if (sum_VF_J != NULL && m == 0 && (v == 0 || v == 2)) {
                process_peak_or_trough(ir, m, v, maxNpt, maxnwait, k, dk, VF_J3, kVFpt, VFpt, iVFpt, gVFpt, iendk0);
            }
            if (sum_HF_J != NULL && m == 1) {
                process_peak_or_trough(ir, m, v, maxNpt, maxnwait, k, dk, HF_J3, kHFpt, HFpt, iHFpt, gHFpt, iendk0);
            }
            if (sum_DC_J != NULL && ((m == 0 && (v == 0 || v == 2)) || m != 0)) {
                process_peak_or_trough(ir, m, v, maxNpt, maxnwait, k, dk, DC_J3, kDCpt, DCpt, iDCpt, gDCpt, iendk0);
            }
        }
    }
    

    // 左移动点, 
    for(MYINT m=0; m<3; ++m){
        for(MYINT v=0; v<4; ++v){
            for(MYINT jj=0; jj<2; ++jj){
                if(sum_EXP_J!=NULL) EXP_J3[ir][jj][m][v] = EXP_J3[ir][jj+1][m][v];
                if(sum_VF_J!=NULL)  VF_J3[ir][jj][m][v]  = VF_J3[ir][jj+1][m][v];
                if(sum_HF_J!=NULL)  HF_J3[ir][jj][m][v]  = HF_J3[ir][jj+1][m][v];
                if(sum_DC_J!=NULL)  DC_J3[ir][jj][m][v]  = DC_J3[ir][jj+1][m][v];
            }

            // 点数+1
            if(sum_EXP_J!=NULL) gEXPpt[ir][m][v]++;
            if(sum_VF_J!=NULL)  gVFpt[ir][m][v]++;
            if(sum_HF_J!=NULL)  gHFpt[ir][m][v]++;
            if(sum_DC_J!=NULL)  gDCpt[ir][m][v]++;
        }
    }
}


void PTA_method(
    const MODEL1D *mod1d, MYREAL k0, MYREAL predk, MYCOMPLEX omega, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_EXP_J0[nr][3][4], MYCOMPLEX sum_VF_J0[nr][3][4],  
    MYCOMPLEX sum_HF_J0[nr][3][4],  MYCOMPLEX sum_DC_J0[nr][3][4],  
    bool calc_upar,
    MYCOMPLEX sum_EXP_uiz_J0[nr][3][4], MYCOMPLEX sum_VF_uiz_J0[nr][3][4],  
    MYCOMPLEX sum_HF_uiz_J0[nr][3][4],  MYCOMPLEX sum_DC_uiz_J0[nr][3][4],  
    MYCOMPLEX sum_EXP_uir_J0[nr][3][4], MYCOMPLEX sum_VF_uir_J0[nr][3][4],  
    MYCOMPLEX sum_HF_uir_J0[nr][3][4],  MYCOMPLEX sum_DC_uir_J0[nr][3][4],  
    FILE *ptam_fstatsnr[nr][2], KernelFunc kerfunc)
{   
    // 需要兼容对正常收敛而不具有规律波峰波谷的序列
    // 有时序列收敛比较好，不表现为规律的波峰波谷，
    // 此时设置最大等待次数，超过直接设置为中间值

    const MYINT maxnwait = 9;     // 最大等待次数，不能太小
    MYREAL k=0.0;

    MYCOMPLEX EXP_qwv[3][3], VF_qwv[3][3], HF_qwv[3][3], DC_qwv[3][3]; // 不同震源的核函数
    MYCOMPLEX (*pEXP_qwv)[3] = (sum_EXP_J0!=NULL)? EXP_qwv : NULL;
    MYCOMPLEX (*pVF_qwv)[3]  = (sum_VF_J0!=NULL)?  VF_qwv  : NULL;
    MYCOMPLEX (*pHF_qwv)[3]  = (sum_HF_J0!=NULL)?  HF_qwv  : NULL;
    MYCOMPLEX (*pDC_qwv)[3]  = (sum_DC_J0!=NULL)?  DC_qwv  : NULL;

    MYCOMPLEX EXP_uiz_qwv[3][3], VF_uiz_qwv[3][3], HF_uiz_qwv[3][3], DC_uiz_qwv[3][3]; 
    MYCOMPLEX (*pEXP_uiz_qwv)[3] = (sum_EXP_uiz_J0!=NULL)? EXP_uiz_qwv : NULL;
    MYCOMPLEX (*pVF_uiz_qwv)[3]  = (sum_VF_uiz_J0!=NULL)?  VF_uiz_qwv  : NULL;
    MYCOMPLEX (*pHF_uiz_qwv)[3]  = (sum_HF_uiz_J0!=NULL)?  HF_uiz_qwv  : NULL;
    MYCOMPLEX (*pDC_uiz_qwv)[3]  = (sum_DC_uiz_J0!=NULL)?  DC_uiz_qwv  : NULL;
    

    static const MYINT maxNpt=PTAM_MAX_PEAK_TROUGH; // 波峰波谷的目标


    // 用于接收F(ki,w)Jm(ki*r)ki
    // 存储采样的值，维度3表示通过连续3个点来判断波峰或波谷
    // 既用于存储被积函数，也最后用于存储求和的结果
    MYCOMPLEX (*EXP_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*EXP_J3));
    MYCOMPLEX (*VF_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*VF_J3));
    MYCOMPLEX (*HF_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*HF_J3));
    MYCOMPLEX (*DC_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*DC_J3));

    MYCOMPLEX (*EXP_uiz_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*EXP_uiz_J3));
    MYCOMPLEX (*VF_uiz_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*VF_uiz_J3));
    MYCOMPLEX (*HF_uiz_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*HF_uiz_J3));
    MYCOMPLEX (*DC_uiz_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*DC_uiz_J3));

    MYCOMPLEX (*EXP_uir_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*EXP_uir_J3));
    MYCOMPLEX (*VF_uir_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*VF_uir_J3));
    MYCOMPLEX (*HF_uir_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*HF_uir_J3));
    MYCOMPLEX (*DC_uir_J3)[3][3][4] = (MYCOMPLEX (*)[3][3][4])calloc(nr, sizeof(*DC_uir_J3));

    // 之前求和的值
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

    MYREAL (*kEXPpt_uiz)[3][4][maxNpt] = (MYREAL (*)[3][4][maxNpt])calloc(nr, sizeof(*kEXPpt_uiz));
    MYREAL (*kVFpt_uiz)[3][4][maxNpt] = (MYREAL (*)[3][4][maxNpt])calloc(nr, sizeof(*kVFpt_uiz));
    MYREAL (*kHFpt_uiz)[3][4][maxNpt] = (MYREAL (*)[3][4][maxNpt])calloc(nr, sizeof(*kHFpt_uiz));
    MYREAL (*kDCpt_uiz)[3][4][maxNpt] = (MYREAL (*)[3][4][maxNpt])calloc(nr, sizeof(*kDCpt_uiz));
    MYCOMPLEX (*EXPpt_uiz)[3][4][maxNpt] = (MYCOMPLEX (*)[3][4][maxNpt])calloc(nr, sizeof(*EXPpt_uiz));
    MYCOMPLEX (*VFpt_uiz)[3][4][maxNpt] = (MYCOMPLEX (*)[3][4][maxNpt])calloc(nr, sizeof(*VFpt_uiz));
    MYCOMPLEX (*HFpt_uiz)[3][4][maxNpt] = (MYCOMPLEX (*)[3][4][maxNpt])calloc(nr, sizeof(*HFpt_uiz));
    MYCOMPLEX (*DCpt_uiz)[3][4][maxNpt] = (MYCOMPLEX (*)[3][4][maxNpt])calloc(nr, sizeof(*DCpt_uiz));
    MYINT (*iEXPpt_uiz)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*iEXPpt_uiz));
    MYINT (*iVFpt_uiz)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*iVFpt_uiz));
    MYINT (*iHFpt_uiz)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*iHFpt_uiz));
    MYINT (*iDCpt_uiz)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*iDCpt_uiz));

    MYREAL (*kEXPpt_uir)[3][4][maxNpt] = (MYREAL (*)[3][4][maxNpt])calloc(nr, sizeof(*kEXPpt_uir));
    MYREAL (*kVFpt_uir)[3][4][maxNpt] = (MYREAL (*)[3][4][maxNpt])calloc(nr, sizeof(*kVFpt_uir));
    MYREAL (*kHFpt_uir)[3][4][maxNpt] = (MYREAL (*)[3][4][maxNpt])calloc(nr, sizeof(*kHFpt_uir));
    MYREAL (*kDCpt_uir)[3][4][maxNpt] = (MYREAL (*)[3][4][maxNpt])calloc(nr, sizeof(*kDCpt_uir));
    MYCOMPLEX (*EXPpt_uir)[3][4][maxNpt] = (MYCOMPLEX (*)[3][4][maxNpt])calloc(nr, sizeof(*EXPpt_uir));
    MYCOMPLEX (*VFpt_uir)[3][4][maxNpt] = (MYCOMPLEX (*)[3][4][maxNpt])calloc(nr, sizeof(*VFpt_uir));
    MYCOMPLEX (*HFpt_uir)[3][4][maxNpt] = (MYCOMPLEX (*)[3][4][maxNpt])calloc(nr, sizeof(*HFpt_uir));
    MYCOMPLEX (*DCpt_uir)[3][4][maxNpt] = (MYCOMPLEX (*)[3][4][maxNpt])calloc(nr, sizeof(*DCpt_uir));
    MYINT (*iEXPpt_uir)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*iEXPpt_uir));
    MYINT (*iVFpt_uir)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*iVFpt_uir));
    MYINT (*iHFpt_uir)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*iHFpt_uir));
    MYINT (*iDCpt_uir)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*iDCpt_uir));

    // 记录点数，当峰谷找到后，清零
    MYINT (*gEXPpt)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gEXPpt));
    MYINT (*gVFpt)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gVFpt));
    MYINT (*gHFpt)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gHFpt));
    MYINT (*gDCpt)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gDCpt));

    MYINT (*gEXPpt_uiz)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gEXPpt_uiz));
    MYINT (*gVFpt_uiz)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gVFpt_uiz));
    MYINT (*gHFpt_uiz)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gHFpt_uiz));
    MYINT (*gDCpt_uiz)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gDCpt_uiz));

    MYINT (*gEXPpt_uir)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gEXPpt_uir));
    MYINT (*gVFpt_uir)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gVFpt_uir));
    MYINT (*gHFpt_uir)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gHFpt_uir));
    MYINT (*gDCpt_uir)[3][4] = (MYINT (*)[3][4])calloc(nr, sizeof(*gDCpt_uir));
    
    
    for(MYINT ir=0; ir<nr; ++ir){
        for(MYINT m=0; m<3; ++m){
            for(MYINT v=0; v<4; ++v){
                if(sum_EXP_J0!=NULL) sum_EXP_J[ir][m][v] = sum_EXP_J0[ir][m][v];
                if(sum_VF_J0!=NULL)  sum_VF_J[ir][m][v]  = sum_VF_J0[ir][m][v];
                if(sum_HF_J0!=NULL)  sum_HF_J[ir][m][v]  = sum_HF_J0[ir][m][v];
                if(sum_DC_J0!=NULL)  sum_DC_J[ir][m][v]  = sum_DC_J0[ir][m][v];

                if(calc_upar){
                    if(sum_EXP_uiz_J0!=NULL) sum_EXP_uiz_J[ir][m][v] = sum_EXP_uiz_J0[ir][m][v];
                    if(sum_VF_uiz_J0!=NULL)  sum_VF_uiz_J[ir][m][v]  = sum_VF_uiz_J0[ir][m][v];
                    if(sum_HF_uiz_J0!=NULL)  sum_HF_uiz_J[ir][m][v]  = sum_HF_uiz_J0[ir][m][v];
                    if(sum_DC_uiz_J0!=NULL)  sum_DC_uiz_J[ir][m][v]  = sum_DC_uiz_J0[ir][m][v];

                    if(sum_EXP_uir_J0!=NULL) sum_EXP_uir_J[ir][m][v] = sum_EXP_uir_J0[ir][m][v];
                    if(sum_VF_uir_J0!=NULL)  sum_VF_uir_J[ir][m][v]  = sum_VF_uir_J0[ir][m][v];
                    if(sum_HF_uir_J0!=NULL)  sum_HF_uir_J[ir][m][v]  = sum_HF_uir_J0[ir][m][v];
                    if(sum_DC_uir_J0!=NULL)  sum_DC_uir_J[ir][m][v]  = sum_DC_uir_J0[ir][m][v];
                }

                iEXPpt[ir][m][v] = gEXPpt[ir][m][v] = 0;
                iVFpt[ir][m][v]  = gVFpt[ir][m][v]  = 0;
                iHFpt[ir][m][v]  = gHFpt[ir][m][v]  = 0;
                iDCpt[ir][m][v]  = gDCpt[ir][m][v]  = 0;

                iEXPpt_uiz[ir][m][v] = gEXPpt_uiz[ir][m][v] = 0;
                iVFpt_uiz[ir][m][v]  = gVFpt_uiz[ir][m][v]  = 0;
                iHFpt_uiz[ir][m][v]  = gHFpt_uiz[ir][m][v]  = 0;
                iDCpt_uiz[ir][m][v]  = gDCpt_uiz[ir][m][v]  = 0;

                iEXPpt_uir[ir][m][v] = gEXPpt_uir[ir][m][v] = 0;
                iVFpt_uir[ir][m][v]  = gVFpt_uir[ir][m][v]  = 0;
                iHFpt_uir[ir][m][v]  = gHFpt_uir[ir][m][v]  = 0;
                iDCpt_uir[ir][m][v]  = gDCpt_uir[ir][m][v]  = 0;

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
            kerfunc(mod1d, omega, k, pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv,
                    calc_upar, pEXP_uiz_qwv, pVF_uiz_qwv, pHF_uiz_qwv, pDC_uiz_qwv); 

            // 记录核函数
            if(fstatsK!=NULL){
                write_stats(
                    fstatsK, k, 
                    EXP_qwv, VF_qwv, HF_qwv, DC_qwv);
            }

            // 计算被积函数一项 F(k,w)Jm(kr)k
            int_Pk(k, rs[ir],
                   pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv, false,
                   EXP_J3[ir][2], VF_J3[ir][2], HF_J3[ir][2], DC_J3[ir][2]);  // [2]表示把新点值放在最后

            // 
            ptam_once(
                ir, nr, precoef, maxNpt, maxnwait, k, dk, 
                EXP_J3, VF_J3, HF_J3, DC_J3, 
                sum_EXP_J, sum_VF_J, sum_HF_J, sum_DC_J, 
                kEXPpt, EXPpt, iEXPpt, gEXPpt, 
                kVFpt, VFpt, iVFpt, gVFpt, 
                kHFpt, HFpt, iHFpt, gHFpt, 
                kDCpt, DCpt, iDCpt, gDCpt, 
                &iendk0);
            
            // -------------------------- 位移空间导数 ------------------------------------
            if(calc_upar){
                // ------------------------------- ui_z -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                int_Pk(k, rs[ir],
                       pEXP_uiz_qwv, pVF_uiz_qwv, pHF_uiz_qwv, pDC_uiz_qwv, false,
                       EXP_uiz_J3[ir][2], VF_uiz_J3[ir][2], HF_uiz_J3[ir][2], DC_uiz_J3[ir][2]);  // [2]表示把新点值放在最后
                
                ptam_once(
                    ir, nr, precoef, maxNpt, maxnwait, k, dk, 
                    EXP_uiz_J3, VF_uiz_J3, HF_uiz_J3, DC_uiz_J3, 
                    sum_EXP_uiz_J, sum_VF_uiz_J, sum_HF_uiz_J, sum_DC_uiz_J, 
                    kEXPpt_uiz, EXPpt_uiz, iEXPpt_uiz, gEXPpt_uiz, 
                    kVFpt_uiz, VFpt_uiz, iVFpt_uiz, gVFpt_uiz, 
                    kHFpt_uiz, HFpt_uiz, iHFpt_uiz, gHFpt_uiz, 
                    kDCpt_uiz, DCpt_uiz, iDCpt_uiz, gDCpt_uiz, 
                    &iendk0);

                // ------------------------------- ui_r -----------------------------------
                // 计算被积函数一项 F(k,w)Jm(kr)k
                int_Pk(k, rs[ir], 
                       pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv, true,
                       EXP_uir_J3[ir][2], VF_uir_J3[ir][2], HF_uir_J3[ir][2], DC_uir_J3[ir][2]);  // [2]表示把新点值放在最后
                
                ptam_once(
                    ir, nr, precoef, maxNpt, maxnwait, k, dk, 
                    EXP_uir_J3, VF_uir_J3, HF_uir_J3, DC_uir_J3, 
                    sum_EXP_uir_J, sum_VF_uir_J, sum_HF_uir_J, sum_DC_uir_J, 
                    kEXPpt_uir, EXPpt_uir, iEXPpt_uir, gEXPpt_uir, 
                    kVFpt_uir, VFpt_uir, iVFpt_uir, gVFpt_uir, 
                    kHFpt_uir, HFpt_uir, iHFpt_uir, gHFpt_uir, 
                    kDCpt_uir, DCpt_uir, iDCpt_uir, gDCpt_uir, 
                    &iendk0);
            
            } // END if calc_upar


            if(iendk0) break;
        }// end k loop
    }

    // printf("w=%f, ik=%d\n", CREAL(omega), ik);


    // 做缩减序列，赋值最终解
    for(MYINT ir=0; ir<nr; ++ir){
        FILE *fstatsP = ptam_fstatsnr[ir][1];
        // 记录到文件
        if(fstatsP!=NULL){
            write_stats_ptam(
                fstatsP, k, maxNpt, 
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
            
                if(calc_upar){
                    if(sum_EXP_uiz_J0!=NULL)  {cplx_shrink(iEXPpt_uiz[ir][m][v], EXPpt_uiz[ir][m][v]);  sum_EXP_uiz_J0[ir][m][v] = EXPpt_uiz[ir][m][v][0];}
                    if(sum_VF_uiz_J0!=NULL)   {cplx_shrink(iVFpt_uiz[ir][m][v],  VFpt_uiz[ir][m][v]);   sum_VF_uiz_J0[ir][m][v]  = VFpt_uiz[ir][m][v][0];}
                    if(sum_HF_uiz_J0!=NULL)   {cplx_shrink(iHFpt_uiz[ir][m][v],  HFpt_uiz[ir][m][v]);   sum_HF_uiz_J0[ir][m][v]  = HFpt_uiz[ir][m][v][0];}
                    if(sum_DC_uiz_J0!=NULL)   {cplx_shrink(iDCpt_uiz[ir][m][v],  DCpt_uiz[ir][m][v]);   sum_DC_uiz_J0[ir][m][v]  = DCpt_uiz[ir][m][v][0];}
                
                    if(sum_EXP_uir_J0!=NULL)  {cplx_shrink(iEXPpt_uir[ir][m][v], EXPpt_uir[ir][m][v]);  sum_EXP_uir_J0[ir][m][v] = EXPpt_uir[ir][m][v][0];}
                    if(sum_VF_uir_J0!=NULL)   {cplx_shrink(iVFpt_uir[ir][m][v],  VFpt_uir[ir][m][v]);   sum_VF_uir_J0[ir][m][v]  = VFpt_uir[ir][m][v][0];}
                    if(sum_HF_uir_J0!=NULL)   {cplx_shrink(iHFpt_uir[ir][m][v],  HFpt_uir[ir][m][v]);   sum_HF_uir_J0[ir][m][v]  = HFpt_uir[ir][m][v][0];}
                    if(sum_DC_uir_J0!=NULL)   {cplx_shrink(iDCpt_uir[ir][m][v],  DCpt_uir[ir][m][v]);   sum_DC_uir_J0[ir][m][v]  = DCpt_uir[ir][m][v][0];}
                }
            }
        }
    }


    free(EXP_J3); free(VF_J3); free(HF_J3); free(DC_J3);
    free(EXP_uiz_J3); free(VF_uiz_J3); free(HF_uiz_J3); free(DC_uiz_J3);
    free(EXP_uir_J3); free(VF_uir_J3); free(HF_uir_J3); free(DC_uir_J3);
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

    free(kEXPpt); free(kVFpt); free(kHFpt); free(kDCpt);
    free(EXPpt);  free(VFpt);  free(HFpt);  free(DCpt);
    free(iEXPpt); free(iVFpt); free(iHFpt); free(iDCpt);
    free(gEXPpt); free(gVFpt); free(gHFpt); free(gDCpt);

    free(kEXPpt_uiz); free(kVFpt_uiz); free(kHFpt_uiz); free(kDCpt_uiz);
    free(EXPpt_uiz);  free(VFpt_uiz);  free(HFpt_uiz);  free(DCpt_uiz);
    free(iEXPpt_uiz); free(iVFpt_uiz); free(iHFpt_uiz); free(iDCpt_uiz);
    free(gEXPpt_uiz); free(gVFpt_uiz); free(gHFpt_uiz); free(gDCpt_uiz);

    free(kEXPpt_uir); free(kVFpt_uir); free(kHFpt_uir); free(kDCpt_uir);
    free(EXPpt_uir);  free(VFpt_uir);  free(HFpt_uir);  free(DCpt_uir);
    free(iEXPpt_uir); free(iVFpt_uir); free(iHFpt_uir); free(iDCpt_uir);
    free(gEXPpt_uir); free(gVFpt_uir); free(gHFpt_uir); free(gDCpt_uir);
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