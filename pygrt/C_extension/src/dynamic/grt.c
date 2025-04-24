/**
 * @file   grt.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 广义反射透射系数矩阵+离散波数法 计算理论地震图，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient 
 *               matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699
 * 
 */

#include <stdio.h> 
#include <sys/stat.h>
#include <errno.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>

#include "dynamic/grt.h"
#include "dynamic/propagate.h"
#include "common/ptam.h"
#include "common/fim.h"
#include "common/dwm.h"
#include "common/integral.h"
#include "common/iostats.h"
#include "common/const.h"
#include "common/model.h"
#include "common/prtdbg.h"
#include "common/search.h"
#include "common/progressbar.h"



void set_num_threads(int num_threads){
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif
}



void integ_grn_spec_in_C(
    PYMODEL1D *pymod1d, MYINT nf1, MYINT nf2, MYINT nf, MYREAL *freqs,  
    MYINT nr, MYREAL *rs, MYREAL wI, 
    MYREAL vmin_ref, MYREAL keps, MYREAL ampk, MYREAL k0, MYREAL Length, MYREAL filonLength, MYREAL filonCut,       
    bool print_progressbar, 

    // 返回值，维度2代表Z、R分量，维度3代表Z、R、T分量
    MYCOMPLEX *EXPcplx[nr][2], // EXZ, EXR 的实部和虚部
    MYCOMPLEX *VFcplx[nr][2],  // VFZ, VFR 的实部和虚部
    MYCOMPLEX *HFcplx[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    MYCOMPLEX *DDcplx[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    MYCOMPLEX *DScplx[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    MYCOMPLEX *SScplx[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]

    bool calc_upar,
    MYCOMPLEX *EXPcplx_uiz[nr][2], // EXZ, EXR 的实部和虚部
    MYCOMPLEX *VFcplx_uiz[nr][2],  // VFZ, VFR 的实部和虚部
    MYCOMPLEX *HFcplx_uiz[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    MYCOMPLEX *DDcplx_uiz[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    MYCOMPLEX *DScplx_uiz[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    MYCOMPLEX *SScplx_uiz[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]
    MYCOMPLEX *EXPcplx_uir[nr][2], // EXZ, EXR 的实部和虚部
    MYCOMPLEX *VFcplx_uir[nr][2],  // VFZ, VFR 的实部和虚部
    MYCOMPLEX *HFcplx_uir[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    MYCOMPLEX *DDcplx_uir[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    MYCOMPLEX *DScplx_uir[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    MYCOMPLEX *SScplx_uir[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]

    const char *statsstr, // 积分结果输出
    MYINT  nstatsidxs, // 仅输出特定频点
    MYINT *statsidxs
){
    // 定义接收结果的GRN结构体
    GRN *(*EXPgrn)[2] = (EXPcplx != NULL) ? (GRN*(*)[2])calloc(nr, sizeof(*EXPgrn)) : NULL;
    GRN *(*VFgrn)[2]  = (VFcplx != NULL) ? (GRN*(*)[2])calloc(nr, sizeof(*VFgrn)) : NULL;
    GRN *(*HFgrn)[3]  = (HFcplx != NULL) ? (GRN*(*)[3])calloc(nr, sizeof(*HFgrn)) : NULL;
    GRN *(*DDgrn)[2]  = (DDcplx != NULL) ? (GRN*(*)[2])calloc(nr, sizeof(*DDgrn)) : NULL;
    GRN *(*DSgrn)[3]  = (DScplx != NULL) ? (GRN*(*)[3])calloc(nr, sizeof(*DSgrn)) : NULL;
    GRN *(*SSgrn)[3]  = (SScplx != NULL) ? (GRN*(*)[3])calloc(nr, sizeof(*SSgrn)) : NULL;

    GRN *(*EXPgrn_uiz)[2] = (EXPcplx_uiz != NULL) ? (GRN*(*)[2])calloc(nr, sizeof(*EXPgrn_uiz)) : NULL;
    GRN *(*VFgrn_uiz)[2]  = (VFcplx_uiz != NULL) ? (GRN*(*)[2])calloc(nr, sizeof(*VFgrn_uiz)) : NULL;
    GRN *(*HFgrn_uiz)[3]  = (HFcplx_uiz != NULL) ? (GRN*(*)[3])calloc(nr, sizeof(*HFgrn_uiz)) : NULL;
    GRN *(*DDgrn_uiz)[2]  = (DDcplx_uiz != NULL) ? (GRN*(*)[2])calloc(nr, sizeof(*DDgrn_uiz)) : NULL;
    GRN *(*DSgrn_uiz)[3]  = (DScplx_uiz != NULL) ? (GRN*(*)[3])calloc(nr, sizeof(*DSgrn_uiz)) : NULL;
    GRN *(*SSgrn_uiz)[3]  = (SScplx_uiz != NULL) ? (GRN*(*)[3])calloc(nr, sizeof(*SSgrn_uiz)) : NULL;

    GRN *(*EXPgrn_uir)[2] = (EXPcplx_uir != NULL) ? (GRN*(*)[2])calloc(nr, sizeof(*EXPgrn_uir)) : NULL;
    GRN *(*VFgrn_uir)[2]  = (VFcplx_uir != NULL) ? (GRN*(*)[2])calloc(nr, sizeof(*VFgrn_uir)) : NULL;
    GRN *(*HFgrn_uir)[3]  = (HFcplx_uir != NULL) ? (GRN*(*)[3])calloc(nr, sizeof(*HFgrn_uir)) : NULL;
    GRN *(*DDgrn_uir)[2]  = (DDcplx_uir != NULL) ? (GRN*(*)[2])calloc(nr, sizeof(*DDgrn_uir)) : NULL;
    GRN *(*DSgrn_uir)[3]  = (DScplx_uir != NULL) ? (GRN*(*)[3])calloc(nr, sizeof(*DSgrn_uir)) : NULL;
    GRN *(*SSgrn_uir)[3]  = (SScplx_uir != NULL) ? (GRN*(*)[3])calloc(nr, sizeof(*SSgrn_uir)) : NULL;
    
    for(int ir=0; ir<nr; ++ir){
        for(int i=0; i<3; ++i){
            if(i<2){
                //
                if(EXPcplx) {
                    EXPgrn[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                    EXPgrn[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                    EXPgrn[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
                }
                if(EXPcplx_uiz) {
                    EXPgrn_uiz[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                    EXPgrn_uiz[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                    EXPgrn_uiz[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
                }
                if(EXPcplx_uir) {
                    EXPgrn_uir[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                    EXPgrn_uir[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                    EXPgrn_uir[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
                }
                //
                if(VFcplx) {
                    VFgrn[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                    VFgrn[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                    VFgrn[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
                }
                if(VFcplx_uiz) {
                    VFgrn_uiz[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                    VFgrn_uiz[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                    VFgrn_uiz[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
                }
                if(VFcplx_uir) {
                    VFgrn_uir[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                    VFgrn_uir[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                    VFgrn_uir[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
                }
                //
                if(DDcplx) {
                    DDgrn[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                    DDgrn[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                    DDgrn[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
                }
                if(DDcplx_uiz) {
                    DDgrn_uiz[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                    DDgrn_uiz[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                    DDgrn_uiz[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
                }
                if(DDcplx_uir) {
                    DDgrn_uir[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                    DDgrn_uir[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                    DDgrn_uir[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
                }
            }
            //
            if(HFcplx) {
                HFgrn[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                HFgrn[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                HFgrn[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
            }
            if(HFcplx_uiz) {
                HFgrn_uiz[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                HFgrn_uiz[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                HFgrn_uiz[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
            }
            if(HFcplx_uir) {
                HFgrn_uir[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                HFgrn_uir[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                HFgrn_uir[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
            }
            //
            if(DScplx) {
                DSgrn[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                DSgrn[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                DSgrn[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
            }
            if(DScplx_uiz) {
                DSgrn_uiz[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                DSgrn_uiz[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                DSgrn_uiz[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
            }
            if(DScplx_uir) {
                DSgrn_uir[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                DSgrn_uir[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                DSgrn_uir[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
            }
            //
            if(SScplx) {
                SSgrn[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                SSgrn[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                SSgrn[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
            }
            if(SScplx_uiz) {
                SSgrn_uiz[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                SSgrn_uiz[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                SSgrn_uiz[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
            }
            if(SScplx_uir) {
                SSgrn_uir[ir][i] = (GRN*)calloc(1, sizeof(GRN));
                SSgrn_uir[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                SSgrn_uir[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
            }
        }
    }


    //==============================================================================
    // 计算格林函数
    integ_grn_spec(
        pymod1d, nf1, nf2, nf, freqs, nr, rs, wI,
        vmin_ref, keps, ampk, k0, Length, filonLength, filonCut, print_progressbar,
        EXPgrn, VFgrn, HFgrn, DDgrn, DSgrn, SSgrn, 
        calc_upar,
        EXPgrn_uiz, VFgrn_uiz, HFgrn_uiz, DDgrn_uiz, DSgrn_uiz, SSgrn_uiz, 
        EXPgrn_uir, VFgrn_uir, HFgrn_uir, DDgrn_uir, DSgrn_uir, SSgrn_uir, 
        statsstr, nstatsidxs, statsidxs
    );
    //==============================================================================
    

    // 写入complex数组
    for(int ir=0; ir<nr; ++ir){
        for(int i=0; i<3; ++i){
            for(int n=nf1; n<=nf2; ++n){
                if(i<2){
                    //
                    if(EXPcplx) EXPcplx[ir][i][n] = CMPLX(EXPgrn[ir][i]->Re[n], EXPgrn[ir][i]->Im[n]);
                    if(EXPcplx_uiz) EXPcplx_uiz[ir][i][n] = CMPLX(EXPgrn_uiz[ir][i]->Re[n], EXPgrn_uiz[ir][i]->Im[n]);
                    if(EXPcplx_uir) EXPcplx_uir[ir][i][n] = CMPLX(EXPgrn_uir[ir][i]->Re[n], EXPgrn_uir[ir][i]->Im[n]);
                    //
                    if(VFcplx) VFcplx[ir][i][n] = CMPLX(VFgrn[ir][i]->Re[n], VFgrn[ir][i]->Im[n]);
                    if(VFcplx_uiz) VFcplx_uiz[ir][i][n] = CMPLX(VFgrn_uiz[ir][i]->Re[n], VFgrn_uiz[ir][i]->Im[n]);
                    if(VFcplx_uir) VFcplx_uir[ir][i][n] = CMPLX(VFgrn_uir[ir][i]->Re[n], VFgrn_uir[ir][i]->Im[n]);
                    if(DDcplx) DDcplx[ir][i][n] = CMPLX(DDgrn[ir][i]->Re[n], DDgrn[ir][i]->Im[n]);
                    //
                    if(DDcplx_uiz) DDcplx_uiz[ir][i][n] = CMPLX(DDgrn_uiz[ir][i]->Re[n], DDgrn_uiz[ir][i]->Im[n]);
                    if(DDcplx_uir) DDcplx_uir[ir][i][n] = CMPLX(DDgrn_uir[ir][i]->Re[n], DDgrn_uir[ir][i]->Im[n]);
                }
                //
                if(HFcplx) HFcplx[ir][i][n] = CMPLX(HFgrn[ir][i]->Re[n], HFgrn[ir][i]->Im[n]);
                if(HFcplx_uiz) HFcplx_uiz[ir][i][n] = CMPLX(HFgrn_uiz[ir][i]->Re[n], HFgrn_uiz[ir][i]->Im[n]);
                if(HFcplx_uir) HFcplx_uir[ir][i][n] = CMPLX(HFgrn_uir[ir][i]->Re[n], HFgrn_uir[ir][i]->Im[n]);
                //
                if(DScplx) DScplx[ir][i][n] = CMPLX(DSgrn[ir][i]->Re[n], DSgrn[ir][i]->Im[n]);
                if(DScplx_uiz) DScplx_uiz[ir][i][n] = CMPLX(DSgrn_uiz[ir][i]->Re[n], DSgrn_uiz[ir][i]->Im[n]);
                if(DScplx_uir) DScplx_uir[ir][i][n] = CMPLX(DSgrn_uir[ir][i]->Re[n], DSgrn_uir[ir][i]->Im[n]);
                //
                if(SScplx) SScplx[ir][i][n] = CMPLX(SSgrn[ir][i]->Re[n], SSgrn[ir][i]->Im[n]);
                if(SScplx_uiz) SScplx_uiz[ir][i][n] = CMPLX(SSgrn_uiz[ir][i]->Re[n], SSgrn_uiz[ir][i]->Im[n]);
                if(SScplx_uir) SScplx_uir[ir][i][n] = CMPLX(SSgrn_uir[ir][i]->Re[n], SSgrn_uir[ir][i]->Im[n]);
            }
        }
    }


    // Free allocated memory
    for(int ir=0; ir<nr; ++ir){
        for(int i=0; i<3; ++i){
            if(i<2){
                //
                if(EXPgrn) {
                    free(EXPgrn[ir][i]->Re);
                    free(EXPgrn[ir][i]->Im);
                    free(EXPgrn[ir][i]);
                }
                if(EXPgrn_uiz) {
                    free(EXPgrn_uiz[ir][i]->Re);
                    free(EXPgrn_uiz[ir][i]->Im);
                    free(EXPgrn_uiz[ir][i]);
                }
                if(EXPgrn_uir) {
                    free(EXPgrn_uir[ir][i]->Re);
                    free(EXPgrn_uir[ir][i]->Im);
                    free(EXPgrn_uir[ir][i]);
                }
                //
                if(VFgrn) {
                    free(VFgrn[ir][i]->Re);
                    free(VFgrn[ir][i]->Im);
                    free(VFgrn[ir][i]);
                }
                if(VFgrn_uiz) {
                    free(VFgrn_uiz[ir][i]->Re);
                    free(VFgrn_uiz[ir][i]->Im);
                    free(VFgrn_uiz[ir][i]);
                }
                if(VFgrn_uir) {
                    free(VFgrn_uir[ir][i]->Re);
                    free(VFgrn_uir[ir][i]->Im);
                    free(VFgrn_uir[ir][i]);
                }
                //
                if(DDgrn) {
                    free(DDgrn[ir][i]->Re);
                    free(DDgrn[ir][i]->Im);
                    free(DDgrn[ir][i]);
                }
                if(DDgrn_uiz) {
                    free(DDgrn_uiz[ir][i]->Re);
                    free(DDgrn_uiz[ir][i]->Im);
                    free(DDgrn_uiz[ir][i]);
                }
                if(DDgrn_uir) {
                    free(DDgrn_uir[ir][i]->Re);
                    free(DDgrn_uir[ir][i]->Im);
                    free(DDgrn_uir[ir][i]);
                }
            }
            //
            if(HFcplx) {
                free(HFgrn[ir][i]->Re);
                free(HFgrn[ir][i]->Im);
                free(HFgrn[ir][i]);
            }
            if(HFcplx_uiz) {
                free(HFgrn_uiz[ir][i]->Re);
                free(HFgrn_uiz[ir][i]->Im);
                free(HFgrn_uiz[ir][i]);
            }
            if(HFcplx_uir) {
                free(HFgrn_uir[ir][i]->Re);
                free(HFgrn_uir[ir][i]->Im);
                free(HFgrn_uir[ir][i]);
            }
            //
            if(DScplx) {
                free(DSgrn[ir][i]->Re);
                free(DSgrn[ir][i]->Im);
                free(DSgrn[ir][i]);
            }
            if(DScplx_uiz) {
                free(DSgrn_uiz[ir][i]->Re);
                free(DSgrn_uiz[ir][i]->Im);
                free(DSgrn_uiz[ir][i]);
            }
            if(DScplx_uir) {
                free(DSgrn_uir[ir][i]->Re);
                free(DSgrn_uir[ir][i]->Im);
                free(DSgrn_uir[ir][i]);
            }
            //
            if(SScplx) {
                free(SSgrn[ir][i]->Re);
                free(SSgrn[ir][i]->Im);
                free(SSgrn[ir][i]);
            }
            if(SScplx_uiz) {
                free(SSgrn_uiz[ir][i]->Re);
                free(SSgrn_uiz[ir][i]->Im);
                free(SSgrn_uiz[ir][i]);
            }
            if(SScplx_uir) {
                free(SSgrn_uir[ir][i]->Re);
                free(SSgrn_uir[ir][i]->Im);
                free(SSgrn_uir[ir][i]);
            }
        }
    }
    if(EXPgrn) free(EXPgrn);
    if(EXPgrn_uiz) free(EXPgrn_uiz);
    if(EXPgrn_uir) free(EXPgrn_uir);
    if(VFgrn) free(VFgrn);
    if(VFgrn_uiz) free(VFgrn_uiz);
    if(VFgrn_uir) free(VFgrn_uir);
    if(HFgrn) free(HFgrn);
    if(HFgrn_uiz) free(HFgrn_uiz);
    if(HFgrn_uir) free(HFgrn_uir);
    if(DDgrn) free(DDgrn);
    if(DDgrn_uiz) free(DDgrn_uiz);
    if(DDgrn_uir) free(DDgrn_uir);
    if(DSgrn) free(DSgrn);
    if(DSgrn_uiz) free(DSgrn_uiz);
    if(DSgrn_uir) free(DSgrn_uir);
    if(SSgrn) free(SSgrn);
    if(SSgrn_uiz) free(SSgrn_uiz);
    if(SSgrn_uir) free(SSgrn_uir);
}


/**
 * 将计算好的复数形式的积分结果记录到GRN结构体中
 * 
 * @param    iw      (in)当前频率索引值
 * @param    nr      (in)震中距个数
 * @param    coef    (in)统一系数
 * @param  sum_EXP_J[nr][3][4]  (in)爆炸源
 * @param  sum_VF_J[nr][3][4]   (in)垂直力源
 * @param  sum_HF_J[nr][3][4]   (in)水平力源
 * @param  sum_DC_J[nr][3][4]   (in)剪切源
 * @param      EXPgrn[nr][2]      (out)`GRN` 结构体指针，爆炸源的Z、R分量频谱结果
 * @param      VFgrn[nr][2]       (out)`GRN` 结构体指针，垂直力源的Z、R分量频谱结果
 * @param      HFgrn[nr][3]       (out)`GRN` 结构体指针，水平力源的Z、R、T分量频谱结果
 * @param      DDgrn[nr][2]       (out)`GRN` 结构体指针，45度倾滑的Z、R分量频谱结果
 * @param      DSgrn[nr][3]       (out)`GRN` 结构体指针，90度倾滑的Z、R、T分量频谱结果
 * @param      SSgrn[nr][3]       (out)`GRN` 结构体指针，90度走滑的Z、R、T分量频谱结果
 */
static void recordin_GRN(
    MYINT iw, MYINT nr, MYCOMPLEX coef, 
    MYCOMPLEX sum_EXP_J[nr][3][4], MYCOMPLEX sum_VF_J[nr][3][4],  
    MYCOMPLEX sum_HF_J[nr][3][4],  MYCOMPLEX sum_DC_J[nr][3][4],  
    GRN *EXPgrn[nr][2], GRN *VFgrn[nr][2], GRN *HFgrn[nr][3],
    GRN *DDgrn[nr][2], GRN *DSgrn[nr][3], GRN *SSgrn[nr][3]
)
{
    // 局部变量，将某个频点的格林函数谱临时存放
    MYCOMPLEX (*tmp_EXP)[2] = (MYCOMPLEX(*)[2])calloc(nr, sizeof(*tmp_EXP));
    MYCOMPLEX (*tmp_VF)[2] = (MYCOMPLEX(*)[2])calloc(nr, sizeof(*tmp_VF));
    MYCOMPLEX (*tmp_HF)[3] = (MYCOMPLEX(*)[3])calloc(nr, sizeof(*tmp_HF));
    MYCOMPLEX (*tmp_DD)[2] = (MYCOMPLEX(*)[2])calloc(nr, sizeof(*tmp_DD));
    MYCOMPLEX (*tmp_DS)[3] = (MYCOMPLEX(*)[3])calloc(nr, sizeof(*tmp_DS));
    MYCOMPLEX (*tmp_SS)[3] = (MYCOMPLEX(*)[3])calloc(nr, sizeof(*tmp_SS));

    for(MYINT ir=0; ir<nr; ++ir){
        for(MYINT ii=0; ii<3; ++ii){
            if(ii<2){
                tmp_EXP[ir][ii] = RZERO;
                tmp_VF[ir][ii] = RZERO; 
                tmp_DD[ir][ii] = RZERO;
            }
            tmp_HF[ir][ii] = RZERO;
            tmp_DS[ir][ii] = RZERO;
            tmp_SS[ir][ii] = RZERO;
        }
        merge_Pk(
            (sum_EXP_J)? sum_EXP_J[ir] : NULL, 
            (sum_VF_J)?  sum_VF_J[ir]  : NULL, 
            (sum_HF_J)?  sum_HF_J[ir]  : NULL, 
            (sum_DC_J)?  sum_DC_J[ir]  : NULL, 
            tmp_EXP[ir], tmp_VF[ir],  tmp_HF[ir], 
            tmp_DD[ir], tmp_DS[ir], tmp_SS[ir]);

        MYCOMPLEX mtmp;
        for(MYINT ii=0; ii<3; ++ii) {
            if(ii<2){
                if(EXPgrn!=NULL){
                    mtmp = coef*tmp_EXP[ir][ii]; // m=0 爆炸源
                    EXPgrn[ir][ii]->Re[iw] = CREAL(mtmp);
                    EXPgrn[ir][ii]->Im[iw] = CIMAG(mtmp);
                }
                if(VFgrn!=NULL){
                    mtmp = coef*tmp_VF[ir][ii]; // m=0 垂直力源
                    VFgrn[ir][ii]->Re[iw] = CREAL(mtmp);
                    VFgrn[ir][ii]->Im[iw] = CIMAG(mtmp);
                }
                if(DDgrn!=NULL){
                    mtmp = coef*tmp_DD[ir][ii]; // m=0 45-倾滑
                    DDgrn[ir][ii]->Re[iw] = CREAL(mtmp);
                    DDgrn[ir][ii]->Im[iw] = CIMAG(mtmp);
                }
            }
            if(HFgrn!=NULL){
                mtmp = coef*tmp_HF[ir][ii]; // m=1 水平力源
                HFgrn[ir][ii]->Re[iw] = CREAL(mtmp);
                HFgrn[ir][ii]->Im[iw] = CIMAG(mtmp);
            }
            if(DSgrn!=NULL){
                mtmp = coef*tmp_DS[ir][ii]; // m=1 90-倾滑
                DSgrn[ir][ii]->Re[iw] = CREAL(mtmp);
                DSgrn[ir][ii]->Im[iw] = CIMAG(mtmp);
            }
            if(SSgrn!=NULL){
                mtmp = coef*tmp_SS[ir][ii]; // m=2 走滑 
                SSgrn[ir][ii]->Re[iw] = CREAL(mtmp);
                SSgrn[ir][ii]->Im[iw] = CIMAG(mtmp);
            }

        }
    }

    free(tmp_EXP);
    free(tmp_VF);
    free(tmp_HF);
    free(tmp_DD);
    free(tmp_DS);
    free(tmp_SS);
}



void integ_grn_spec(
    PYMODEL1D *pymod1d, MYINT nf1, MYINT nf2, MYINT nf, MYREAL *freqs,  
    MYINT nr, MYREAL *rs, MYREAL wI, 
    MYREAL vmin_ref, MYREAL keps, MYREAL ampk, MYREAL k0, MYREAL Length, MYREAL filonLength, MYREAL filonCut,      
    bool print_progressbar, 

    // 返回值，维度2代表Z、R分量，维度3代表Z、R、T分量
    GRN *EXPgrn[nr][2], // EXZ, EXR 的实部和虚部
    GRN *VFgrn[nr][2],  // VFZ, VFR 的实部和虚部
    GRN *HFgrn[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    GRN *DDgrn[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    GRN *DSgrn[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    GRN *SSgrn[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]

    bool calc_upar,
    GRN *EXPgrn_uiz[nr][2], // EXZ, EXR 的实部和虚部
    GRN *VFgrn_uiz[nr][2],  // VFZ, VFR 的实部和虚部
    GRN *HFgrn_uiz[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    GRN *DDgrn_uiz[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    GRN *DSgrn_uiz[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    GRN *SSgrn_uiz[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]
    GRN *EXPgrn_uir[nr][2], // EXZ, EXR 的实部和虚部
    GRN *VFgrn_uir[nr][2],  // VFZ, VFR 的实部和虚部
    GRN *HFgrn_uir[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    GRN *DDgrn_uir[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    GRN *DSgrn_uir[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    GRN *SSgrn_uir[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]

    const char *statsstr, // 积分结果输出
    MYINT  nstatsidxs, // 仅输出特定频点
    MYINT *statsidxs
){
    // 程序运行开始时间
    struct timeval begin_t;
    gettimeofday(&begin_t, NULL);

    // 最大震中距
    MYINT irmax = findMinMax_MYREAL(rs, nr, true);
    MYREAL rmax=rs[irmax];   

    // pymod1d -> mod1d
    MODEL1D *main_mod1d = init_mod1d(pymod1d->n);
    get_mod1d(pymod1d, main_mod1d);

    const LAYER *src_lay = main_mod1d->lays + main_mod1d->isrc;
    const MYREAL Rho = src_lay->Rho; // 震源区密度
    const MYREAL fac = RONE/(RFOUR*PI*Rho);
    const MYREAL hs = (FABS(pymod1d->depsrc - pymod1d->deprcv) < MIN_DEPTH_GAP_SRC_RCV)? 
                      MIN_DEPTH_GAP_SRC_RCV : FABS(pymod1d->depsrc - pymod1d->deprcv); // hs=max(震源和台站深度差,1.0)

    // 乘相应系数
    k0 *= PI/hs;
    const MYREAL k02 = k0*k0;
    const MYREAL ampk2 = ampk*ampk;

    if(vmin_ref < RZERO)  keps = -RONE;  // 若使用峰谷平均法，则不使用keps进行收敛判断

    const MYREAL dk=PI2/(Length*rmax);     // 波数积分间隔
    const MYREAL filondk = (filonLength > RZERO) ? PI2/(filonLength*rmax) : RZERO;  // Filon积分间隔
    const MYREAL filonK = filonCut/rmax;  // 波数积分和Filon积分的分割点


    // PTAM的积分中间结果, 每个震中距两个文件，因为PTAM对不同震中距使用不同的dk
    // 在文件名后加后缀，区分不同震中距
    char *ptam_fstatsdir[nr];
    for(MYINT ir=0; ir<nr; ++ir) {ptam_fstatsdir[ir] = NULL;}
    if(statsstr!=NULL && nstatsidxs > 0 && vmin_ref < RZERO){
        for(MYINT ir=0; ir<nr; ++ir){
            ptam_fstatsdir[ir] = (char*)malloc((strlen(statsstr)+200)*sizeof(char));
            ptam_fstatsdir[ir][0] = '\0';
            // 新建文件夹目录 
            sprintf(ptam_fstatsdir[ir], "%s/PTAM_%04d_%.5e", statsstr, ir, rs[ir]);
            if(mkdir(ptam_fstatsdir[ir], 0777) != 0){
                if(errno != EEXIST){
                    printf("Unable to create folder %s. Error code: %d\n", ptam_fstatsdir[ir], errno);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }


    // 进度条变量 
    MYINT progress=0;

    // 频率omega循环
    // schedule语句可以动态调度任务，最大程度地使用计算资源
    #pragma omp parallel for schedule(guided) default(shared) 
    for(MYINT iw=nf1; iw<=nf2; ++iw){
        MYREAL k=RZERO;               // 波数
        MYREAL w = freqs[iw]*PI2;     // 实频率
        MYCOMPLEX omega = w - wI*I; // 复数频率 omega = w - i*wI
        MYCOMPLEX omega2_inv = RONE/omega; // 1/omega^2
        omega2_inv = omega2_inv*omega2_inv; 
        MYCOMPLEX coef = -dk*fac*omega2_inv; // 最终要乘上的系数

        // 局部变量，用于求和 sum F(ki,w)Jm(ki*r)ki 
        // 维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型，详见int_Pk()函数内的注释
        MYCOMPLEX (*sum_EXP_J)[3][4] = (EXPgrn != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_EXP_J)) : NULL;
        MYCOMPLEX (*sum_VF_J)[3][4] = (VFgrn != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_VF_J)) : NULL;
        MYCOMPLEX (*sum_HF_J)[3][4] = (HFgrn != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_HF_J)) : NULL;
        MYCOMPLEX (*sum_DC_J)[3][4] = (DDgrn != NULL || DSgrn != NULL || SSgrn != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_DC_J)) : NULL;

        MYCOMPLEX (*sum_EXP_uiz_J)[3][4] = (EXPgrn_uiz != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_EXP_uiz_J)) : NULL;
        MYCOMPLEX (*sum_VF_uiz_J)[3][4] = (VFgrn_uiz != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_VF_uiz_J)) : NULL;
        MYCOMPLEX (*sum_HF_uiz_J)[3][4] = (HFgrn_uiz != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_HF_uiz_J)) : NULL;
        MYCOMPLEX (*sum_DC_uiz_J)[3][4] = (DDgrn_uiz != NULL || DSgrn_uiz != NULL || SSgrn_uiz != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_DC_uiz_J)) : NULL;

        MYCOMPLEX (*sum_EXP_uir_J)[3][4] = (EXPgrn_uir != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_EXP_uir_J)) : NULL;
        MYCOMPLEX (*sum_VF_uir_J)[3][4] = (VFgrn_uir != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_VF_uir_J)) : NULL;
        MYCOMPLEX (*sum_HF_uir_J)[3][4] = (HFgrn_uir != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_HF_uir_J)) : NULL;
        MYCOMPLEX (*sum_DC_uir_J)[3][4] = (DDgrn_uir != NULL || DSgrn_uir != NULL || SSgrn_uir != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_DC_uir_J)) : NULL;

        

        MODEL1D *local_mod1d = NULL;
    #ifdef _OPENMP 
        // 定义局部模型对象
        local_mod1d = init_mod1d(main_mod1d->n);
        copy_mod1d(main_mod1d, local_mod1d);
    #else 
        local_mod1d = main_mod1d;
    #endif
        update_mod1d_omega(local_mod1d, omega);

        // 是否要输出积分过程文件
        bool needfstats = (statsstr!=NULL && ((findElement_MYINT(statsidxs, nstatsidxs, iw) >= 0) || (findElement_MYINT(statsidxs, nstatsidxs, -1) >= 0)));

        // 为当前频率创建波数积分记录文件
        FILE *fstats = NULL;
        // PTAM为每个震中距都创建波数积分记录文件
        FILE *(*ptam_fstatsnr)[2] = (FILE *(*)[2])malloc(nr * sizeof(*ptam_fstatsnr));
        {
            MYINT len0 = (statsstr!=NULL) ? strlen(statsstr) : 0;
            char *fname = (char *)malloc((len0+200)*sizeof(char));
            if(needfstats){
                sprintf(fname, "%s/K_%04d_%.5e", statsstr, iw, freqs[iw]);
                fstats = fopen(fname, "wb");
            }
            for(MYINT ir=0; ir<nr; ++ir){
                for(MYINT m=0; m<3; ++m){
                    for(MYINT v=0; v<4; ++v){
                        if(sum_EXP_J) sum_EXP_J[ir][m][v] = RZERO;
                        if(sum_VF_J) sum_VF_J[ir][m][v] = RZERO;
                        if(sum_HF_J) sum_HF_J[ir][m][v] = RZERO;
                        if(sum_DC_J) sum_DC_J[ir][m][v] = RZERO;
    
                        if(sum_EXP_uiz_J) sum_EXP_uiz_J[ir][m][v] = RZERO;
                        if(sum_VF_uiz_J) sum_VF_uiz_J[ir][m][v] = RZERO;
                        if(sum_HF_uiz_J) sum_HF_uiz_J[ir][m][v] = RZERO;
                        if(sum_DC_uiz_J) sum_DC_uiz_J[ir][m][v] = RZERO;
    
                        if(sum_EXP_uir_J) sum_EXP_uir_J[ir][m][v] = RZERO;
                        if(sum_VF_uir_J) sum_VF_uir_J[ir][m][v] = RZERO;
                        if(sum_HF_uir_J) sum_HF_uir_J[ir][m][v] = RZERO;
                        if(sum_DC_uir_J) sum_DC_uir_J[ir][m][v] = RZERO;
                    }
                }
    
                ptam_fstatsnr[ir][0] = ptam_fstatsnr[ir][1] = NULL;
                if(needfstats && vmin_ref < RZERO){
                    // 峰谷平均法
                    sprintf(fname, "%s/K_%04d_%.5e", ptam_fstatsdir[ir], iw, freqs[iw]);
                    ptam_fstatsnr[ir][0] = fopen(fname, "wb");
                    sprintf(fname, "%s/PTAM_%04d_%.5e", ptam_fstatsdir[ir], iw, freqs[iw]);
                    ptam_fstatsnr[ir][1] = fopen(fname, "wb");
                }
            } // end init rs loop
            free(fname);
        }

        


        MYREAL kmax;
        // vmin_ref的正负性在这里不影响
        kmax = SQRT(k02 + ampk2*(w/vmin_ref)*(w/vmin_ref));


        // 常规的波数积分
        k = discrete_integ(
            local_mod1d, dk, (filondk > RZERO)? filonK : kmax, keps, omega, nr, rs, 
            sum_EXP_J, sum_VF_J, sum_HF_J, sum_DC_J, 
            calc_upar,
            sum_EXP_uiz_J, sum_VF_uiz_J, sum_HF_uiz_J, sum_DC_uiz_J,
            sum_EXP_uir_J, sum_VF_uir_J, sum_HF_uir_J, sum_DC_uir_J,
            fstats, kernel);
            
        // 基于线性插值的Filon积分
        if(filondk > RZERO){
            k = linear_filon_integ(
                local_mod1d, k, dk, filondk, kmax, keps, omega, nr, rs, 
                sum_EXP_J, sum_VF_J, sum_HF_J, sum_DC_J, 
                calc_upar,
                sum_EXP_uiz_J, sum_VF_uiz_J, sum_HF_uiz_J, sum_DC_uiz_J,
                sum_EXP_uir_J, sum_VF_uir_J, sum_HF_uir_J, sum_DC_uir_J,
                fstats, kernel);
        }

        // k之后的部分使用峰谷平均法进行显式收敛，建议在浅源地震的时候使用   
        if(vmin_ref < RZERO){
            PTA_method(
                local_mod1d, k, dk, omega, nr, rs, 
                sum_EXP_J, sum_VF_J, sum_HF_J, sum_DC_J, 
                calc_upar,
                sum_EXP_uiz_J, sum_VF_uiz_J, sum_HF_uiz_J, sum_DC_uiz_J,
                sum_EXP_uir_J, sum_VF_uir_J, sum_HF_uir_J, sum_DC_uir_J,
                ptam_fstatsnr, kernel);
        }

        // printf("iw=%d, w=%.5e, k=%.5e, dk=%.5e, nk=%d\n", iw, w, k, dk, (int)(k/dk));



        // 记录到格林函数结构体内
        recordin_GRN(
            iw, nr, coef, 
            sum_EXP_J, sum_VF_J, sum_HF_J, sum_DC_J,
            EXPgrn, VFgrn, HFgrn, DDgrn, DSgrn, SSgrn);
        if(calc_upar){
            recordin_GRN(
                iw, nr, coef, 
                sum_EXP_uiz_J, sum_VF_uiz_J, sum_HF_uiz_J, sum_DC_uiz_J,
                EXPgrn_uiz, VFgrn_uiz, HFgrn_uiz, DDgrn_uiz, DSgrn_uiz, SSgrn_uiz);
            recordin_GRN(
                iw, nr, coef, 
                sum_EXP_uir_J, sum_VF_uir_J, sum_HF_uir_J, sum_DC_uir_J,
                EXPgrn_uir, VFgrn_uir, HFgrn_uir, DDgrn_uir, DSgrn_uir, SSgrn_uir);
        }
        

        if(fstats!=NULL) fclose(fstats);
        for(MYINT ir=0; ir<nr; ++ir){
            if(ptam_fstatsnr[ir][0]!=NULL){
                fclose(ptam_fstatsnr[ir][0]);
            }
            if(ptam_fstatsnr[ir][1]!=NULL){
                fclose(ptam_fstatsnr[ir][1]);
            }
        }
        free(ptam_fstatsnr);

    #ifdef _OPENMP
        free_mod1d(local_mod1d);
    #endif

        // 记录进度条变量 
        #pragma omp critical
        {
            progress++;
            if(print_progressbar) printprogressBar("Computing Green Functions: ", progress*100/(nf2-nf1+1));
        } 
        



        // Free allocated memory for temporary variables
        if (sum_EXP_J) free(sum_EXP_J);
        if (sum_VF_J) free(sum_VF_J);
        if (sum_HF_J) free(sum_HF_J);
        if (sum_DC_J) free(sum_DC_J);

        if (sum_EXP_uiz_J) free(sum_EXP_uiz_J);
        if (sum_VF_uiz_J) free(sum_VF_uiz_J);
        if (sum_HF_uiz_J) free(sum_HF_uiz_J);
        if (sum_DC_uiz_J) free(sum_DC_uiz_J);

        if (sum_EXP_uir_J) free(sum_EXP_uir_J);
        if (sum_VF_uir_J) free(sum_VF_uir_J);
        if (sum_HF_uir_J) free(sum_HF_uir_J);
        if (sum_DC_uir_J) free(sum_DC_uir_J);

    } // END omega loop



    free_mod1d(main_mod1d);

    for(MYINT ir=0; ir<nr; ++ir){
        if(ptam_fstatsdir[ir]!=NULL){
            free(ptam_fstatsdir[ir]);
        } 
    }

    // 程序运行结束时间
    struct timeval end_t;
    gettimeofday(&end_t, NULL);
    if(print_progressbar) printf("Runtime: %.3f s\n", (end_t.tv_sec - begin_t.tv_sec) + (end_t.tv_usec - begin_t.tv_usec) / 1e6);
    fflush(stdout);
}






