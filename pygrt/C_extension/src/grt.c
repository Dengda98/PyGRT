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

#include "grt.h"
#include "dwm.h"
#include "const.h"
#include "model.h"
#include "propagate.h"
#include "prtdbg.h"
#include "search.h"
#include "ptam.h"
#include "fim.h"
#include "iostats.h"
#include "progressbar.h"



void set_num_threads(int num_threads){
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif
}



void integ_grn_spec_in_C(
    PYMODEL1D *pymod1d, MYINT nf1, MYINT nf2, MYINT nf, MYREAL *freqs,  
    MYINT nr, MYREAL *rs, MYREAL wI, 
    MYREAL vmin_ref, MYREAL keps, MYREAL ampk, bool iwk0, MYREAL k0, MYREAL Length,       
    bool print_progressbar, 

    // 返回值，维度2代表Z、R分量，维度3代表Z、R、T分量
    MYCOMPLEX *EXPcplx[nr][2], // EXZ, EXR 的实部和虚部
    MYCOMPLEX *VFcplx[nr][2],  // VFZ, VFR 的实部和虚部
    MYCOMPLEX *HFcplx[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    MYCOMPLEX *DDcplx[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    MYCOMPLEX *DScplx[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    MYCOMPLEX *SScplx[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]

    const char *statsstr, // 积分结果输出
    MYINT  nstatsidxs, // 仅输出特定频点
    MYINT *statsidxs
){
    // 定义接收结果的GRN结构体
    GRN *EXPgrn[nr][2];
    GRN *VFgrn[nr][2];
    GRN *HFgrn[nr][3];
    GRN *DDgrn[nr][2];
    GRN *DSgrn[nr][3];
    GRN *SSgrn[nr][3];
    GRN *((*pEXPgrn)[2]) = (EXPcplx!=NULL)? EXPgrn : NULL;
    GRN *((*pVFgrn)[2])  = (VFcplx!=NULL) ?  VFgrn : NULL;
    GRN *((*pHFgrn)[3])  = (HFcplx!=NULL) ?  HFgrn : NULL;
    GRN *((*pDDgrn)[2])  = (DDcplx!=NULL) ?  DDgrn : NULL;
    GRN *((*pDSgrn)[3])  = (DScplx!=NULL) ?  DSgrn : NULL;
    GRN *((*pSSgrn)[3])  = (SScplx!=NULL) ?  SSgrn : NULL;
    for(int ir=0; ir<nr; ++ir){
        for(int i=0; i<3; ++i){
            if(i<2){
                if(EXPcplx!=NULL){
                    EXPgrn[ir][i] = (GRN*)malloc(sizeof(GRN));
                    EXPgrn[ir][i]->nf = nf;
                    EXPgrn[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                    EXPgrn[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
                } else {
                    EXPgrn[ir][i] = NULL;
                }
                
                if(VFcplx!=NULL){
                    VFgrn[ir][i] = (GRN*)malloc(sizeof(GRN));
                    VFgrn[ir][i]->nf = nf;
                    VFgrn[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                    VFgrn[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
                } else {
                    VFgrn[ir][i] = NULL;
                }
                
                if(DDcplx!=NULL){
                    DDgrn[ir][i] = (GRN*)malloc(sizeof(GRN));
                    DDgrn[ir][i]->nf = nf;
                    DDgrn[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                    DDgrn[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
                } else {
                    DDgrn[ir][i] = NULL;
                }
                
            }
            
            if(HFcplx!=NULL){
                HFgrn[ir][i] = (GRN*)malloc(sizeof(GRN));
                HFgrn[ir][i]->nf = nf;
                HFgrn[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                HFgrn[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
            } else {
                HFgrn[ir][i] = NULL;
            }
            
            if(DScplx!=NULL){
                DSgrn[ir][i] = (GRN*)malloc(sizeof(GRN));
                DSgrn[ir][i]->nf = nf;
                DSgrn[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                DSgrn[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
            } else {
                DSgrn[ir][i] = NULL;
            }
            
            if(SScplx!=NULL){
                SSgrn[ir][i] = (GRN*)malloc(sizeof(GRN));
                SSgrn[ir][i]->nf = nf;
                SSgrn[ir][i]->Re = (MYREAL*)calloc(nf, sizeof(MYREAL));
                SSgrn[ir][i]->Im = (MYREAL*)calloc(nf, sizeof(MYREAL));
            } else {
                SSgrn[ir][i] = NULL;
            }
            
        }
    }


    //==============================================================================
    // 计算格林函数
    integ_grn_spec(
        pymod1d, nf1, nf2, nf, freqs, nr, rs, wI,
        vmin_ref, keps, ampk, iwk0, k0, Length, print_progressbar,
        pEXPgrn, pVFgrn, pHFgrn, pDDgrn, pDSgrn, pSSgrn, 
        statsstr, nstatsidxs, statsidxs
    );
    //==============================================================================
    

    // 写入complex数组
    for(int ir=0; ir<nr; ++ir){
        for(int i=0; i<3; ++i){
            for(int n=nf1; n<=nf2; ++n){
                if(i<2){
                    if(EXPcplx!=NULL){
                        EXPcplx[ir][i][n] = CMPLX(EXPgrn[ir][i]->Re[n], EXPgrn[ir][i]->Im[n]);
                    } 
                    
                    if(VFcplx!=NULL){
                        VFcplx[ir][i][n] = CMPLX(VFgrn[ir][i]->Re[n], VFgrn[ir][i]->Im[n]);
                    } 
                    
                    if(DDcplx!=NULL){
                        DDcplx[ir][i][n] = CMPLX(DDgrn[ir][i]->Re[n], DDgrn[ir][i]->Im[n]);
                    } 
                }
                
                if(HFcplx!=NULL){
                    HFcplx[ir][i][n] = CMPLX(HFgrn[ir][i]->Re[n], HFgrn[ir][i]->Im[n]);
                } 
                
                if(DScplx!=NULL){
                    DScplx[ir][i][n] = CMPLX(DSgrn[ir][i]->Re[n], DSgrn[ir][i]->Im[n]);
                } 
                
                if(SScplx!=NULL){
                    SScplx[ir][i][n] = CMPLX(SSgrn[ir][i]->Re[n], SSgrn[ir][i]->Im[n]);
                } 
            }
        }
    }
}


void integ_grn_spec(
    PYMODEL1D *pymod1d, MYINT nf1, MYINT nf2, MYINT nf, MYREAL *freqs,  
    MYINT nr, MYREAL *rs, MYREAL wI, 
    MYREAL vmin_ref, MYREAL keps, MYREAL ampk, bool iwk0, MYREAL k0, MYREAL Length,       
    bool print_progressbar, 

    // 返回值，维度2代表Z、R分量，维度3代表Z、R、T分量
    GRN *EXPgrn[nr][2], // EXZ, EXR 的实部和虚部
    GRN *VFgrn[nr][2],  // VFZ, VFR 的实部和虚部
    GRN *HFgrn[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    GRN *DDgrn[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    GRN *DSgrn[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    GRN *SSgrn[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]

    const char *statsstr, // 积分结果输出
    MYINT  nstatsidxs, // 仅输出特定频点
    MYINT *statsidxs
){
    // 程序运行开始时间
    struct timeval begin_t;
    gettimeofday(&begin_t, NULL);

    MYREAL rmin=rs[findMinMax_MYREAL(rs, nr, false)];  // 最小震中距
    MYREAL rmax=rs[findMinMax_MYREAL(rs, nr, true)];   // 最大震中距

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


    const MYREAL wmax = freqs[nf-1]/PI2;  // 最大圆频率
    const MYREAL dk=FABS(PI2/(Length*rmax));     // 波数积分间隔


    // 输出波数积分中间结果, 每个震中距一个文件
    // 在文件名后加后缀，区分不同震中距
    char *fstatsdir[nr];
    for(MYINT ir=0; ir<nr; ++ir) {fstatsdir[ir] = NULL;}
    if(statsstr!=NULL && nstatsidxs > 0){
        for(MYINT ir=0; ir<nr; ++ir){
            fstatsdir[ir] = (char*)malloc((strlen(statsstr)+200)*sizeof(char));
            fstatsdir[ir][0] = '\0';
            // 新建文件夹目录 
            sprintf(fstatsdir[ir], "%s_%.3f_%.3f_%.3f", statsstr, pymod1d->depsrc, pymod1d->deprcv, rs[ir]);
            if(mkdir(fstatsdir[ir], 0777) != 0){
                if(errno != EEXIST){
                    printf("Unable to create folder %s. Error code: %d\n", fstatsdir[ir], errno);
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

        // 局部变量，将某个频点的格林函数谱临时存放
        MYCOMPLEX tmp_EXP[nr][2], tmp_VF[nr][2], tmp_HF[nr][3];
        MYCOMPLEX tmp_DD[nr][2], tmp_DS[nr][3], tmp_SS[nr][3];

        // 局部变量，用于求和 sum F(ki,w)Jm(ki*r)ki 
        // 维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型，详见int_Pk()函数内的注释
        MYCOMPLEX sum_EXP_J[nr][3][4], sum_VF_J[nr][3][4], sum_HF_J[nr][3][4], sum_DC_J[nr][3][4];
        MYCOMPLEX (*psum_EXP_J)[3][4] = (EXPgrn!=NULL)? sum_EXP_J : NULL;
        MYCOMPLEX (*psum_VF_J)[3][4]  = (VFgrn!=NULL)?  sum_VF_J  : NULL;
        MYCOMPLEX (*psum_HF_J)[3][4]  = (HFgrn!=NULL)?  sum_HF_J  : NULL;
        MYCOMPLEX (*psum_DC_J)[3][4]  = (DDgrn!=NULL || DSgrn!=NULL || SSgrn!=NULL )?  sum_DC_J  : NULL;

        MODEL1D *local_mod1d = NULL;
    #ifdef _OPENMP 
        // 定义局部模型对象
        local_mod1d = init_mod1d(main_mod1d->n);
        copy_mod1d(main_mod1d, local_mod1d);
    #else 
        local_mod1d = main_mod1d;
    #endif
        update_mod1d_omega(local_mod1d, omega);



        // 给每个频率创建波数积分记录文件
        FILE *(fstats[nr]);
        FILE *(ptam_fstats[nr]);

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

            for(MYINT m=0; m<3; ++m){
                for(MYINT v=0; v<4; ++v){
                    sum_EXP_J[ir][m][v] = sum_VF_J[ir][m][v] = 
                    sum_HF_J[ir][m][v]  = sum_DC_J[ir][m][v] = RZERO;
                }
            }


            fstats[ir] = NULL;
            ptam_fstats[ir] = NULL;
            if(statsstr!=NULL && ((findElement_MYINT(statsidxs, nstatsidxs, iw) >= 0) || (findElement_MYINT(statsidxs, nstatsidxs, -1) >= 0))){
                char *fname = (char *)malloc((strlen(fstatsdir[ir])+200)*sizeof(char));
                if(Length > 0){
                    // 常规的波数积分
                    sprintf(fname, "%s/K_%d_%.3e", fstatsdir[ir], iw, freqs[iw]);
                } else {
                    // Filon积分
                    sprintf(fname, "%s/Filon_%d_%.3e", fstatsdir[ir], iw, freqs[iw]);
                }
                
                fstats[ir] = fopen(fname, "wb");

                if(vmin_ref < 0){
                    // 峰谷平均法
                    sprintf(fname, "%s/PTAM_%d_%.3e", fstatsdir[ir], iw, freqs[iw]);
                    ptam_fstats[ir] = fopen(fname, "wb");
                }
                free(fname);
            }
            
        } // end init rs loop


        MYREAL kmax;
        MYREAL ampk02 = RONE;
        if(iwk0)  ampk02 = (w/wmax)*(w/wmax);
        // vmin_ref的正负性在这里不影响
        kmax = SQRT(k02*ampk02 + ampk2*(w/vmin_ref)*(w/vmin_ref));


        if(Length > RZERO){
            // 常规的波数积分
            k = discrete_integ(
                local_mod1d, dk, kmax, keps, omega, nr, rs, 
                psum_EXP_J, psum_VF_J, psum_HF_J, psum_DC_J, fstats);
        } 
        else {
            // 基于线性插值的Filon积分
            k = linear_filon_integ(
                local_mod1d, dk, kmax, keps, omega, nr, rs, 
                psum_EXP_J, psum_VF_J, psum_HF_J, psum_DC_J, fstats);
        }

        // k之后的部分使用峰谷平均法进行显式收敛，建议在浅源地震的时候使用   
        if(vmin_ref < RZERO){
            PTA_method(
                local_mod1d, k, dk, rmin, rmax, omega, nr, rs, 
                psum_EXP_J, psum_VF_J, psum_HF_J, psum_DC_J, fstats, ptam_fstats);
        }

        // printf("iw=%d, w=%.5e, k=%.5e, dk=%.5e, nk=%d\n", iw, w, k, dk, (int)(k/dk));



        // 记录到格林函数结构体内
        for(MYINT ir=0; ir<nr; ++ir){
            merge_Pk(
                (psum_EXP_J!=NULL)? psum_EXP_J[ir] : NULL, 
                (psum_VF_J!=NULL)?  psum_VF_J[ir]  : NULL, 
                (psum_HF_J!=NULL)?  psum_HF_J[ir]  : NULL, 
                (psum_DC_J!=NULL)?  psum_DC_J[ir]  : NULL, 
                tmp_EXP[ir], tmp_VF[ir],  tmp_HF[ir], 
                tmp_DD[ir], tmp_DS[ir], tmp_SS[ir]);

            MYCOMPLEX mtmp0 = dk * fac * omega2_inv * (-RONE); // 
            MYCOMPLEX mtmp;
            for(MYINT ii=0; ii<3; ++ii) {
                if(ii<2){
                    if(EXPgrn!=NULL){
                        mtmp = mtmp0*tmp_EXP[ir][ii]; // m=0 爆炸源
                        EXPgrn[ir][ii]->Re[iw] = CREAL(mtmp);
                        EXPgrn[ir][ii]->Im[iw] = CIMAG(mtmp);
                    }
                    if(VFgrn!=NULL){
                        mtmp = mtmp0*tmp_VF[ir][ii]; // m=0 垂直力源
                        VFgrn[ir][ii]->Re[iw] = CREAL(mtmp);
                        VFgrn[ir][ii]->Im[iw] = CIMAG(mtmp);
                    }
                    if(DDgrn!=NULL){
                        mtmp = mtmp0*tmp_DD[ir][ii]; // m=0 45-倾滑
                        DDgrn[ir][ii]->Re[iw] = CREAL(mtmp);
                        DDgrn[ir][ii]->Im[iw] = CIMAG(mtmp);
                    }
                }
                if(HFgrn!=NULL){
                    mtmp = mtmp0*tmp_HF[ir][ii]; // m=1 水平力源
                    HFgrn[ir][ii]->Re[iw] = CREAL(mtmp);
                    HFgrn[ir][ii]->Im[iw] = CIMAG(mtmp);
                }
                if(DSgrn!=NULL){
                    mtmp = mtmp0*tmp_DS[ir][ii]; // m=1 90-倾滑
                    DSgrn[ir][ii]->Re[iw] = CREAL(mtmp);
                    DSgrn[ir][ii]->Im[iw] = CIMAG(mtmp);
                }
                if(SSgrn!=NULL){
                    mtmp = mtmp0*tmp_SS[ir][ii]; // m=2 走滑 
                    SSgrn[ir][ii]->Re[iw] = CREAL(mtmp);
                    SSgrn[ir][ii]->Im[iw] = CIMAG(mtmp);
                }

            }
        }

        for(MYINT ir=0; ir<nr; ++ir){
            if(fstats[ir]!=NULL){
                fclose(fstats[ir]);
            }
            if(ptam_fstats[ir]!=NULL){
                fclose(ptam_fstats[ir]);
            }
        }

    #ifdef _OPENMP
        free_mod1d(local_mod1d);
    #endif

        // 记录进度条变量 
        #pragma omp critical
        {
            progress++;
            if(print_progressbar) printprogressBar("Computing Green Functions: ", progress*100/(nf2-nf1+1));
        } 


    } // END omega loop



    free_mod1d(main_mod1d);

    for(MYINT ir=0; ir<nr; ++ir){
        if(fstatsdir[ir]!=NULL){
            free(fstatsdir[ir]);
        } 
    }

    // 程序运行结束时间
    struct timeval end_t;
    gettimeofday(&end_t, NULL);
    if(print_progressbar) printf("Runtime: %.3f s\n", (end_t.tv_sec - begin_t.tv_sec) + (end_t.tv_usec - begin_t.tv_usec) / 1e6);
    fflush(stdout);
}






