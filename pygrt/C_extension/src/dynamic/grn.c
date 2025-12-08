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

#include "grt/common/kernel.h"
#include "grt/dynamic/grn.h"
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
    size_t iw, size_t nr, cplx_t coef, cplxIntegGrid sum_J[nr],
    pt_cplxChnlGrid grn[nr])
{
    // 局部变量，将某个频点的格林函数谱临时存放
    cplxChnlGrid *tmp_grn = (cplxChnlGrid *)calloc(nr, sizeof(*tmp_grn));

    for(size_t ir=0; ir<nr; ++ir){
        grt_merge_Pk(sum_J[ir], tmp_grn[ir]);

        GRT_LOOP_ChnlGrid(im, c){
            int modr = GRT_SRC_M_ORDERS[im];
            if(modr == 0 && GRT_ZRT_CODES[c] == 'T')  continue;

            grn[ir][im][c][iw] = coef * tmp_grn[ir][im][c];
        }
    }

    GRT_SAFE_FREE_PTR(tmp_grn);
}



void grt_integ_grn_spec(
    GRT_MODEL1D *mod1d, size_t nf1, size_t nf2, real_t *freqs,  
    size_t nr, real_t *rs, real_t wI, bool keepAllFreq,
    real_t vmin_ref, real_t keps, real_t ampk, real_t k0, real_t Length, real_t filonLength, real_t safilonTol, real_t filonCut,      
    bool print_progressbar, 

    // 返回值，代表Z、R、T分量
    pt_cplxChnlGrid grn[nr],

    bool calc_upar,
    pt_cplxChnlGrid grn_uiz[nr],
    pt_cplxChnlGrid grn_uir[nr],

    const char *statsstr, // 积分结果输出
    size_t  nstatsidxs, // 仅输出特定频点
    size_t *statsidxs
){
    // 程序运行开始时间
    struct timeval begin_t;
    gettimeofday(&begin_t, NULL);

    // 最大震中距
    size_t irmax = grt_findMax_real_t(rs, nr);
    real_t rmax=rs[irmax];   

    const real_t Rho = mod1d->Rho[mod1d->isrc]; // 震源区密度
    const real_t fac = 1.0/(4.0*PI*Rho);
    const real_t hs = GRT_MAX(fabs(mod1d->depsrc - mod1d->deprcv), GRT_MIN_DEPTH_GAP_SRC_RCV); // hs=max(震源和台站深度差,1.0)

    // 乘相应系数
    k0 *= PI/hs;
    const real_t k02 = k0*k0;
    const real_t ampk2 = ampk*ampk;

    if(vmin_ref < 0.0)  keps = 0.0;  // 若使用峰谷平均法，则不使用keps进行收敛判断

    bool useFIM = (filonLength > 0.0) || (safilonTol > 0.0) ;    // 是否使用Filon积分（包括自适应Filon）
    const real_t dk=PI2/(Length*rmax);     // 波数积分间隔
    const real_t filondk = (filonLength > 0.0) ? PI2/(filonLength*rmax) : 0.0;  // Filon积分间隔
    const real_t filonK = filonCut/rmax;  // 波数积分和Filon积分的分割点


    // PTAM的积分中间结果, 每个震中距两个文件，因为PTAM对不同震中距使用不同的dk
    // 在文件名后加后缀，区分不同震中距
    char **ptam_fstatsdir = (char**)calloc(nr, sizeof(char*));
    if(statsstr!=NULL && nstatsidxs > 0 && vmin_ref < 0.0){
        for(size_t ir=0; ir<nr; ++ir){
            // 新建文件夹目录 
            GRT_SAFE_ASPRINTF(&ptam_fstatsdir[ir], "%s/PTAM_%04zu_%.5e", statsstr, ir, rs[ir]);
            if(mkdir(ptam_fstatsdir[ir], 0777) != 0){
                if(errno != EEXIST){
                    printf("Unable to create folder %s. Error code: %d\n", ptam_fstatsdir[ir], errno);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }


    // 进度条变量 
    int progress=0;

    // 记录每个频率的计算中是否有除0错误
    int *freq_invstats = (int *)calloc(nf2+1, sizeof(int));

    // 实际计算的频点数
    size_t nf_valid = nf2 - nf1 + 1;

    // 频率omega循环
    // schedule语句可以动态调度任务，最大程度地使用计算资源
    #pragma omp parallel for schedule(guided) default(shared) 
    for(size_t iw=nf1; iw<=nf2; ++iw){
        real_t k=0.0;               // 波数
        real_t w = freqs[iw]*PI2;     // 实频率
        cplx_t omega = w - wI*I; // 复数频率 omega = w - i*wI

        // 如果在虚频率的帮助下，频率仍然距离原点太近，
        // 计算会有严重的数值问题，因此直接根据频率距离原点的距离，
        // 跳过该频率，没有必要再计算
        if( ! keepAllFreq )
        {
            real_t R = 0.1; // 完全经验性地设定，暂不必要暴露在用户可控层面
            if(fabs(omega) < R){
                #pragma omp critical
                {
                    GRTRaiseWarning("Skip low frequency (iw=%zu, freq=%.5e).", iw, freqs[iw]);
                    nf_valid--;
                }
                if(nf_valid == 0)  GRTRaiseError("NO VALID FREQUENCIES.");
                continue;
            }
        }

        cplx_t omega2_inv = 1.0/omega; // 1/omega^2
        omega2_inv = omega2_inv*omega2_inv; 
        cplx_t coef = -dk*fac*omega2_inv; // 最终要乘上的系数

        // 局部变量，用于求和 sum F(ki,w)Jm(ki*r)ki 
        // 关于形状详见int_Pk()函数内的注释
        cplxIntegGrid *sum_J = (cplxIntegGrid *)calloc(nr, sizeof(*sum_J));
        cplxIntegGrid *sum_uiz_J = (calc_upar)? (cplxIntegGrid *)calloc(nr, sizeof(*sum_uiz_J)) : NULL;
        cplxIntegGrid *sum_uir_J = (calc_upar)? (cplxIntegGrid *)calloc(nr, sizeof(*sum_uir_J)) : NULL;

        GRT_MODEL1D *local_mod1d = NULL;
    #ifdef _OPENMP 
        // 定义局部模型对象
        local_mod1d = grt_copy_mod1d(mod1d);
    #else 
        local_mod1d = mod1d;
    #endif

        // 将 omega 计入模型结构体
        local_mod1d->omega = omega;

        grt_attenuate_mod1d(local_mod1d, omega);

        // 是否要输出积分过程文件
        bool needfstats = (statsstr!=NULL && (grt_findElement_size_t(statsidxs, nstatsidxs, iw) >= 0));

        // 为当前频率创建波数积分记录文件
        FILE *fstats = NULL;
        // PTAM为每个震中距都创建波数积分记录文件
        FILE *(*ptam_fstatsnr)[2] = (FILE *(*)[2])malloc(nr * sizeof(*ptam_fstatsnr));
        {
            char *fname = NULL;
            if(needfstats){
                GRT_SAFE_ASPRINTF(&fname, "%s/K_%04zu_%.5e", statsstr, iw, freqs[iw]);
                fstats = fopen(fname, "wb");
            }
            for(size_t ir=0; ir<nr; ++ir){
                ptam_fstatsnr[ir][0] = ptam_fstatsnr[ir][1] = NULL;
                if(needfstats && vmin_ref < 0.0){
                    // 峰谷平均法
                    GRT_SAFE_ASPRINTF(&fname, "%s/K_%04zu_%.5e", ptam_fstatsdir[ir], iw, freqs[iw]);
                    ptam_fstatsnr[ir][0] = fopen(fname, "wb");
                    GRT_SAFE_ASPRINTF(&fname, "%s/PTAM_%04zu_%.5e", ptam_fstatsdir[ir], iw, freqs[iw]);
                    ptam_fstatsnr[ir][1] = fopen(fname, "wb");
                }
            }
            GRT_SAFE_FREE_PTR(fname);
        }

        


        real_t kmax;
        // vmin_ref的正负性在这里不影响
        kmax = sqrt(k02 + ampk2*(w/vmin_ref)*(w/vmin_ref));

        // 计算核函数过程中是否有遇到除零错误
        freq_invstats[iw]=GRT_INVERSE_SUCCESS;


        // ===================================================================================
        //                          Wavenumber Integration
        // 常规的波数积分
        k = grt_discrete_integ(
            local_mod1d, dk, (useFIM)? filonK : kmax, keps, nr, rs, 
            sum_J, calc_upar, sum_uiz_J, sum_uir_J,
            fstats, grt_kernel);
        if(local_mod1d->stats==GRT_INVERSE_FAILURE)  goto NEXT_FREQ;
    
        // 使用Filon积分
        if(useFIM){
            if(filondk > 0.0){
                // 基于线性插值的Filon积分，固定采样间隔
                k = grt_linear_filon_integ(
                    local_mod1d, k, dk, filondk, kmax, keps, nr, rs, 
                    sum_J, calc_upar, sum_uiz_J, sum_uir_J,
                    fstats, grt_kernel);
            }
            else if(safilonTol > 0.0){
                // 基于自适应采样的Filon积分
                k = grt_sa_filon_integ(
                    local_mod1d, k, dk, safilonTol, kmax, creal(omega)/fabs(vmin_ref)*ampk, nr, rs, 
                    sum_J, calc_upar, sum_uiz_J, sum_uir_J,
                    fstats, grt_kernel);
            }
            if(local_mod1d->stats==GRT_INVERSE_FAILURE)  goto NEXT_FREQ;
        }

        // k之后的部分使用峰谷平均法进行显式收敛，建议在浅源地震的时候使用   
        if(vmin_ref < 0.0){
            grt_PTA_method(
                local_mod1d, k, dk, nr, rs, 
                sum_J, calc_upar, sum_uiz_J, sum_uir_J,
                ptam_fstatsnr, grt_kernel);
            if(local_mod1d->stats==GRT_INVERSE_FAILURE)  goto NEXT_FREQ;
        }

        // fprintf(stderr, "iw=%d, w=%.5e, k=%.5e, dk=%.5e, nk=%d\n", iw, w, k, dk, (int)(k/dk));
        // fflush(stderr);

        // 记录到格林函数结构体内
        // 如果计算核函数过程中存在除零错误，则放弃该频率【通常在大震中距的低频段】
        recordin_GRN(iw, nr, coef, sum_J, grn);
        if(calc_upar){
            recordin_GRN(iw, nr, coef, sum_uiz_J, grn_uiz);
            recordin_GRN(iw, nr, coef, sum_uir_J, grn_uir);
        }
        // ===================================================================================

        // 如果有什么计算意外，从以上的波数积分部分跳至此处
        NEXT_FREQ:
        freq_invstats[iw] = local_mod1d->stats;


        if(fstats!=NULL) fclose(fstats);
        for(size_t ir=0; ir<nr; ++ir){
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
            if(print_progressbar) grt_printprogressBar("Computing Green Functions: ", progress*100/nf_valid);
        } 
        

        // Free allocated memory for temporary variables
        GRT_SAFE_FREE_PTR(sum_J);
        GRT_SAFE_FREE_PTR(sum_uiz_J);
        GRT_SAFE_FREE_PTR(sum_uir_J);

    } // END omega loop


    GRT_SAFE_FREE_PTR_ARRAY(ptam_fstatsdir, nr);

    // 打印 freq_invstats
    for(size_t iw=nf1; iw<=nf2; ++iw){
        if(freq_invstats[iw]==GRT_INVERSE_FAILURE){
            fprintf(stderr, "iw=%zu, freq=%e(Hz), meet Zero Divison Error, results are filled with 0.\n", iw, freqs[iw]);
        }
    }
    GRT_SAFE_FREE_PTR(freq_invstats);

    // 程序运行结束时间
    struct timeval end_t;
    gettimeofday(&end_t, NULL);
    if(print_progressbar) printf("Runtime: %.3f s\n", (end_t.tv_sec - begin_t.tv_sec) + (end_t.tv_usec - begin_t.tv_usec) / 1e6);
    fflush(stdout);

}






