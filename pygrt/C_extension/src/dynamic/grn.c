/**
 * @file   grn.c
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

#include "grt/dynamic/grn.h"
#include "grt/dynamic/propagate.h"
#include "grt/common/ptam.h"
#include "grt/common/fim.h"
#include "grt/common/safim.h"
#include "grt/common/dwm.h"
#include "grt/common/integral.h"
#include "grt/common/iostats.h"
#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/common/prtdbg.h"
#include "grt/common/search.h"
#include "grt/common/progressbar.h"


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
    MYINT iw, MYINT nr, MYCOMPLEX coef, MYCOMPLEX sum_J[nr][GRT_SRC_M_NUM][GRT_INTEG_NUM],
    MYCOMPLEX *grn[nr][GRT_SRC_M_NUM][GRT_CHANNEL_NUM]
)
{
    // 局部变量，将某个频点的格林函数谱临时存放
    MYCOMPLEX (*tmp_grn)[GRT_SRC_M_NUM][GRT_CHANNEL_NUM] = (MYCOMPLEX(*)[GRT_SRC_M_NUM][GRT_CHANNEL_NUM])calloc(nr, sizeof(*tmp_grn));

    for(MYINT ir=0; ir<nr; ++ir){
        grt_merge_Pk(sum_J[ir], tmp_grn[ir]);

        for(MYINT i=0; i<GRT_SRC_M_NUM; ++i) {
            MYINT modr = GRT_SRC_M_ORDERS[i];
            for(MYINT c=0; c<GRT_CHANNEL_NUM; ++c){
                if(modr == 0 && GRT_ZRT_CODES[c] == 'T')  continue;

                grn[ir][i][c][iw] = coef * tmp_grn[ir][i][c];
            }

        }
    }

    GRT_SAFE_FREE_PTR(tmp_grn);
}



void grt_integ_grn_spec(
    GRT_MODEL1D *mod1d, MYINT nf1, MYINT nf2, MYREAL *freqs,  
    MYINT nr, MYREAL *rs, MYREAL wI, 
    MYREAL vmin_ref, MYREAL keps, MYREAL ampk, MYREAL k0, MYREAL Length, MYREAL filonLength, MYREAL safilonTol, MYREAL filonCut,      
    bool print_progressbar, 

    // 返回值，代表Z、R、T分量
    MYCOMPLEX *grn[nr][GRT_SRC_M_NUM][GRT_CHANNEL_NUM],

    bool calc_upar,
    MYCOMPLEX *grn_uiz[nr][GRT_SRC_M_NUM][GRT_CHANNEL_NUM],
    MYCOMPLEX *grn_uir[nr][GRT_SRC_M_NUM][GRT_CHANNEL_NUM],

    const char *statsstr, // 积分结果输出
    MYINT  nstatsidxs, // 仅输出特定频点
    MYINT *statsidxs
){
    // 程序运行开始时间
    struct timeval begin_t;
    gettimeofday(&begin_t, NULL);

    // 最大震中距
    MYINT irmax = grt_findMax_MYREAL(rs, nr);
    MYREAL rmax=rs[irmax];   

    const MYREAL Rho = mod1d->Rho[mod1d->isrc]; // 震源区密度
    const MYREAL fac = 1.0/(4.0*PI*Rho);
    const MYREAL hs = GRT_MAX(fabs(mod1d->depsrc - mod1d->deprcv), GRT_MIN_DEPTH_GAP_SRC_RCV); // hs=max(震源和台站深度差,1.0)

    // 乘相应系数
    k0 *= PI/hs;
    const MYREAL k02 = k0*k0;
    const MYREAL ampk2 = ampk*ampk;

    if(vmin_ref < 0.0)  keps = 0.0;  // 若使用峰谷平均法，则不使用keps进行收敛判断

    bool useFIM = (filonLength > 0.0) || (safilonTol > 0.0) ;    // 是否使用Filon积分（包括自适应Filon）
    const MYREAL dk=PI2/(Length*rmax);     // 波数积分间隔
    const MYREAL filondk = (filonLength > 0.0) ? PI2/(filonLength*rmax) : 0.0;  // Filon积分间隔
    const MYREAL filonK = filonCut/rmax;  // 波数积分和Filon积分的分割点


    // PTAM的积分中间结果, 每个震中距两个文件，因为PTAM对不同震中距使用不同的dk
    // 在文件名后加后缀，区分不同震中距
    char **ptam_fstatsdir = (char**)calloc(nr, sizeof(char*));
    if(statsstr!=NULL && nstatsidxs > 0 && vmin_ref < 0.0){
        for(MYINT ir=0; ir<nr; ++ir){
            // 新建文件夹目录 
            GRT_SAFE_ASPRINTF(&ptam_fstatsdir[ir], "%s/PTAM_%04d_%.5e", statsstr, ir, rs[ir]);
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

    // 每个频率的计算中是否有除0错误
    MYINT freq_invstats[nf2+1];

    // 频率omega循环
    // schedule语句可以动态调度任务，最大程度地使用计算资源
    #pragma omp parallel for schedule(guided) default(shared) 
    for(MYINT iw=nf1; iw<=nf2; ++iw){
        MYREAL k=0.0;               // 波数
        MYREAL w = freqs[iw]*PI2;     // 实频率
        MYCOMPLEX omega = w - wI*I; // 复数频率 omega = w - i*wI
        MYCOMPLEX omega2_inv = 1.0/omega; // 1/omega^2
        omega2_inv = omega2_inv*omega2_inv; 
        MYCOMPLEX coef = -dk*fac*omega2_inv; // 最终要乘上的系数

        // 局部变量，用于求和 sum F(ki,w)Jm(ki*r)ki 
        // 关于形状详见int_Pk()函数内的注释
        MYCOMPLEX (*sum_J)[GRT_SRC_M_NUM][GRT_INTEG_NUM] = (MYCOMPLEX(*)[GRT_SRC_M_NUM][GRT_INTEG_NUM])calloc(nr, sizeof(*sum_J));
        MYCOMPLEX (*sum_uiz_J)[GRT_SRC_M_NUM][GRT_INTEG_NUM] = (calc_upar)? (MYCOMPLEX(*)[GRT_SRC_M_NUM][GRT_INTEG_NUM])calloc(nr, sizeof(*sum_uiz_J)) : NULL;
        MYCOMPLEX (*sum_uir_J)[GRT_SRC_M_NUM][GRT_INTEG_NUM] = (calc_upar)? (MYCOMPLEX(*)[GRT_SRC_M_NUM][GRT_INTEG_NUM])calloc(nr, sizeof(*sum_uir_J)) : NULL;

        GRT_MODEL1D *local_mod1d = NULL;
    #ifdef _OPENMP 
        // 定义局部模型对象
        local_mod1d = grt_copy_mod1d(mod1d);
    #else 
        local_mod1d = mod1d;
    #endif
        grt_attenuate_mod1d(local_mod1d, omega);

        // 是否要输出积分过程文件
        bool needfstats = (statsstr!=NULL && ((grt_findElement_MYINT(statsidxs, nstatsidxs, iw) >= 0) || (grt_findElement_MYINT(statsidxs, nstatsidxs, -1) >= 0)));

        // 为当前频率创建波数积分记录文件
        FILE *fstats = NULL;
        // PTAM为每个震中距都创建波数积分记录文件
        FILE *(*ptam_fstatsnr)[2] = (FILE *(*)[2])malloc(nr * sizeof(*ptam_fstatsnr));
        {
            char *fname = NULL;
            if(needfstats){
                GRT_SAFE_ASPRINTF(&fname, "%s/K_%04d_%.5e", statsstr, iw, freqs[iw]);
                fstats = fopen(fname, "wb");
            }
            for(MYINT ir=0; ir<nr; ++ir){
                for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
                    for(MYINT v=0; v<GRT_INTEG_NUM; ++v){
                        sum_J[ir][i][v] = 0.0;
                        if(calc_upar){
                            sum_uiz_J[ir][i][v] = 0.0;
                            sum_uir_J[ir][i][v] = 0.0;
                        }
                    }
                }
    
                ptam_fstatsnr[ir][0] = ptam_fstatsnr[ir][1] = NULL;
                if(needfstats && vmin_ref < 0.0){
                    // 峰谷平均法
                    GRT_SAFE_ASPRINTF(&fname, "%s/K_%04d_%.5e", ptam_fstatsdir[ir], iw, freqs[iw]);
                    ptam_fstatsnr[ir][0] = fopen(fname, "wb");
                    GRT_SAFE_ASPRINTF(&fname, "%s/PTAM_%04d_%.5e", ptam_fstatsdir[ir], iw, freqs[iw]);
                    ptam_fstatsnr[ir][1] = fopen(fname, "wb");
                }
            } // end init rs loop
            GRT_SAFE_FREE_PTR(fname);
        }

        


        MYREAL kmax;
        // vmin_ref的正负性在这里不影响
        kmax = sqrt(k02 + ampk2*(w/vmin_ref)*(w/vmin_ref));

        // 计算核函数过程中是否有遇到除零错误
        freq_invstats[iw]=GRT_INVERSE_SUCCESS;

        // 常规的波数积分
        k = grt_discrete_integ(
            local_mod1d, dk, (useFIM)? filonK : kmax, keps, omega, nr, rs, 
            sum_J, calc_upar, sum_uiz_J, sum_uir_J,
            fstats, grt_kernel, &freq_invstats[iw]);
    
        // 使用Filon积分
        if(useFIM && freq_invstats[iw]==GRT_INVERSE_SUCCESS){
            if(filondk > 0.0){
                // 基于线性插值的Filon积分，固定采样间隔
                k = grt_linear_filon_integ(
                    local_mod1d, k, dk, filondk, kmax, keps, omega, nr, rs, 
                    sum_J, calc_upar, sum_uiz_J, sum_uir_J,
                    fstats, grt_kernel, &freq_invstats[iw]);
            }
            else if(safilonTol > 0.0){
                // 基于自适应采样的Filon积分
                k = grt_sa_filon_integ(
                    local_mod1d, fabs(vmin_ref)/ampk, k, dk, safilonTol, kmax, omega, nr, rs, 
                    sum_J, calc_upar, sum_uiz_J, sum_uir_J,
                    fstats, grt_kernel, &freq_invstats[iw]);
            }
        }

        // k之后的部分使用峰谷平均法进行显式收敛，建议在浅源地震的时候使用   
        if(vmin_ref < 0.0 && freq_invstats[iw]==GRT_INVERSE_SUCCESS){
            grt_PTA_method(
                local_mod1d, k, dk, omega, nr, rs, 
                sum_J, calc_upar, sum_uiz_J, sum_uir_J,
                ptam_fstatsnr, grt_kernel, &freq_invstats[iw]);
        }

        // fprintf(stderr, "iw=%d, w=%.5e, k=%.5e, dk=%.5e, nk=%d\n", iw, w, k, dk, (int)(k/dk));
        // fflush(stderr);

        // 记录到格林函数结构体内
        // 如果计算核函数过程中存在除零错误，则放弃该频率【通常在大震中距的低频段】
        if(freq_invstats[iw]==GRT_INVERSE_SUCCESS){ 
            recordin_GRN(iw, nr, coef, sum_J, grn);
            if(calc_upar){
                recordin_GRN(iw, nr, coef, sum_uiz_J, grn_uiz);
                recordin_GRN(iw, nr, coef, sum_uir_J, grn_uir);
            }
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
        GRT_SAFE_FREE_PTR(ptam_fstatsnr);

    #ifdef _OPENMP
        grt_free_mod1d(local_mod1d);
    #endif

        // 记录进度条变量 
        #pragma omp critical
        {
            progress++;
            if(print_progressbar) grt_printprogressBar("Computing Green Functions: ", progress*100/(nf2-nf1+1));
        } 
        

        // Free allocated memory for temporary variables
        GRT_SAFE_FREE_PTR(sum_J);
        GRT_SAFE_FREE_PTR(sum_uiz_J);
        GRT_SAFE_FREE_PTR(sum_uir_J);

    } // END omega loop


    GRT_SAFE_FREE_PTR_ARRAY(ptam_fstatsdir, nr);

    // 打印freq_invstats
    for(MYINT iw=nf1; iw<=nf2; ++iw){
        if(freq_invstats[iw]==GRT_INVERSE_FAILURE){
            fprintf(stderr, "iw=%d, freq=%e(Hz), meet Zero Divison Error, results are filled with 0.\n", iw, freqs[iw]);
        }
    }

    // 程序运行结束时间
    struct timeval end_t;
    gettimeofday(&end_t, NULL);
    if(print_progressbar) printf("Runtime: %.3f s\n", (end_t.tv_sec - begin_t.tv_sec) + (end_t.tv_usec - begin_t.tv_usec) / 1e6);
    fflush(stdout);

}






