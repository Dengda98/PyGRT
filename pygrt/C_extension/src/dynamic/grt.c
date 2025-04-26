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


/**
 * 将计算好的复数形式的积分结果记录到GRN结构体中
 * 
 * @param[in]    iw             当前频率索引值
 * @param[in]    nr             震中距个数
 * @param[in]    coef           统一系数
 * @param[in]    sum_J          积分结果
 * @param[out]   grn            三分量频谱
 */
static void recordin_GRN(
    MYINT iw, MYINT nr, MYCOMPLEX coef, MYCOMPLEX sum_J[nr][SRC_M_NUM][INTEG_NUM],
    MYCOMPLEX *grn[nr][SRC_M_NUM][CHANNEL_NUM]
)
{
    // 局部变量，将某个频点的格林函数谱临时存放
    MYCOMPLEX (*tmp_grn)[SRC_M_NUM][CHANNEL_NUM] = (MYCOMPLEX(*)[SRC_M_NUM][CHANNEL_NUM])calloc(nr, sizeof(*tmp_grn));

    for(MYINT ir=0; ir<nr; ++ir){
        merge_Pk(sum_J[ir], tmp_grn[ir]);

        for(MYINT i=0; i<SRC_M_NUM; ++i) {
            for(MYINT c=0; c<CHANNEL_NUM; ++c){
                grn[ir][i][c][iw] = coef * tmp_grn[ir][i][c];
            }

        }
    }

    free(tmp_grn);
}



void integ_grn_spec(
    PYMODEL1D *pymod1d, MYINT nf1, MYINT nf2, MYINT nf, MYREAL *freqs,  
    MYINT nr, MYREAL *rs, MYREAL wI, 
    MYREAL vmin_ref, MYREAL keps, MYREAL ampk, MYREAL k0, MYREAL Length, MYREAL filonLength, MYREAL filonCut,      
    bool print_progressbar, 

    // 返回值，代表Z、R、T分量
    MYCOMPLEX *grn[nr][SRC_M_NUM][CHANNEL_NUM],

    bool calc_upar,
    MYCOMPLEX *grn_uiz[nr][SRC_M_NUM][CHANNEL_NUM],
    MYCOMPLEX *grn_uir[nr][SRC_M_NUM][CHANNEL_NUM],

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
        // 关于形状详见int_Pk()函数内的注释
        MYCOMPLEX (*sum_J)[SRC_M_NUM][INTEG_NUM] = (MYCOMPLEX(*)[SRC_M_NUM][INTEG_NUM])calloc(nr, sizeof(*sum_J));
        MYCOMPLEX (*sum_uiz_J)[SRC_M_NUM][INTEG_NUM] = (calc_upar)? (MYCOMPLEX(*)[SRC_M_NUM][INTEG_NUM])calloc(nr, sizeof(*sum_uiz_J)) : NULL;
        MYCOMPLEX (*sum_uir_J)[SRC_M_NUM][INTEG_NUM] = (calc_upar)? (MYCOMPLEX(*)[SRC_M_NUM][INTEG_NUM])calloc(nr, sizeof(*sum_uir_J)) : NULL;

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
                for(MYINT i=0; i<SRC_M_NUM; ++i){
                    for(MYINT v=0; v<INTEG_NUM; ++v){
                        sum_J[ir][i][v] = CZERO;
                        if(calc_upar){
                            sum_uiz_J[ir][i][v] = CZERO;
                            sum_uir_J[ir][i][v] = CZERO;
                        }
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
            sum_J, calc_upar, sum_uiz_J, sum_uir_J,
            fstats, kernel);
            
        // 基于线性插值的Filon积分
        if(filondk > RZERO){
            k = linear_filon_integ(
                local_mod1d, k, dk, filondk, kmax, keps, omega, nr, rs, 
                sum_J, calc_upar, sum_uiz_J, sum_uir_J,
                fstats, kernel);
        }

        // k之后的部分使用峰谷平均法进行显式收敛，建议在浅源地震的时候使用   
        if(vmin_ref < RZERO){
            PTA_method(
                local_mod1d, k, dk, omega, nr, rs, 
                sum_J, calc_upar, sum_uiz_J, sum_uir_J,
                ptam_fstatsnr, kernel);
        }

        // printf("iw=%d, w=%.5e, k=%.5e, dk=%.5e, nk=%d\n", iw, w, k, dk, (int)(k/dk));



        // 记录到格林函数结构体内
        recordin_GRN(iw, nr, coef, sum_J, grn);
        if(calc_upar){
            recordin_GRN(iw, nr, coef, sum_uiz_J, grn_uiz);
            recordin_GRN(iw, nr, coef, sum_uir_J, grn_uir);
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
        free(sum_J);
        if(sum_uiz_J) free(sum_uiz_J);
        if(sum_uir_J) free(sum_uir_J);

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






