/**
 * @file   grn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是利用多种波数积分类方法计算理论地震图，参考：
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
#include "grt/dynamic/grnspec.h"
#include "grt/integral/integ_method.h"

#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/common/search.h"
#include "grt/common/progressbar.h"


/**
 * 将计算好的复数形式的积分结果转化为频谱
 * 
 * @param[in]    iw             当前频率索引值
 * @param[in]    nr             震中距个数
 * @param[in]    coef           统一系数
 * @param[in]    sumJ            积分结果
 * @param[out]   u              三分量频谱
 */
static void recordin_GRN(size_t iw, size_t nr, cplx_t coef, cplxIntegGrid sumJ[nr], pcplxChnlGrid u[nr])
{
    // 局部变量，将某个频点的格林函数谱临时存放
    cplxChnlGrid *tmp_u = (cplxChnlGrid *)calloc(nr, sizeof(*tmp_u));

    for(size_t ir=0; ir<nr; ++ir){
        grt_merge_Pk(sumJ[ir], tmp_u[ir]);

        GRT_LOOP_ChnlGrid(im, c){
            int modr = GRT_SRC_M_ORDERS[im];
            if(modr == 0 && GRT_ZRT_CODES[c] == 'T')  continue;

            u[ir][im][c][iw] = coef * tmp_u[ir][im][c];
        }
    }

    GRT_SAFE_FREE_PTR(tmp_u);
}



void grt_integ_grn_spec(MODEL1D *mod1d, K_INTEG_METHOD *Kmet, GRNSPEC *grn, const bool print_progressbar)
{
    // 程序运行开始时间
    struct timeval begin_t;
    gettimeofday(&begin_t, NULL);

    const real_t Rho = mod1d->Rho[mod1d->isrc]; // 震源区密度
    const real_t fac = 1.0/(4.0*PI*Rho);
    const real_t dk = Kmet->dk;

    // 进度条变量 
    int progress=0;

    // 记录每个频率的计算中是否有除0错误
    int *freq_invstats = (int *)calloc(grn->nf2 + 1, sizeof(int));

    // 实际计算的频点数
    size_t nf_valid = grn->nf2 - grn->nf1 + 1;

    mod1d->omgref = PI2*grn->freqs[grn->nf2];

    // 频率omega循环
    // schedule语句可以动态调度任务，最大程度地使用计算资源
    #pragma omp parallel for schedule(guided) default(shared) 
    for(size_t iw = grn->nf1; iw <= grn->nf2; ++iw){
        real_t w = grn->freqs[iw]*PI2;     // 实频率
        cplx_t omega = w - grn->wI*I; // 复数频率 omega = w - i*wI

        // 如果在虚频率的帮助下，频率仍然距离原点太近，
        // 计算会有严重的数值问题，因此直接根据频率距离原点的距离，
        // 跳过该频率，没有必要再计算
        if( ! grn->keepAllFreq )
        {
            real_t R = 0.1; // 完全经验性地设定，暂不必要暴露在用户可控层面
            if(fabs(omega) < R){
                #pragma omp critical
                {
                    GRTRaiseWarning("Skip low frequency (iw=%zu, freq=%.5e).", iw, grn->freqs[iw]);
                    nf_valid--;
                }
                if(nf_valid == 0)  GRTRaiseError("NO VALID FREQUENCIES.");
                continue;
            }
        }

        cplx_t coef = - dk*fac / GRT_SQUARE(omega); // 最终要乘上的系数

        K_INTEG_METHOD *local_Kmet = NULL;
    #ifdef _OPENMP 
        // 定义局部对象
        K_INTEG_METHOD __KMET = *Kmet;  // 只需浅拷贝
        local_Kmet = &__KMET;
    #else 
        local_Kmet = Kmet;
    #endif

        MODEL1D_STATE *local_mstat = grt_init_mod1d_state(mod1d);

        // 将 omega 计入模型结构体
        grt_update_mod1d_state_omega(local_mstat, omega, false);

        // 是否要输出积分过程文件
        bool needfstats = (grn->statsstr!=NULL && (grt_findElement_size_t(grn->statsidxs, grn->nstatsidxs, iw) >= 0));

        // 为当前频率创建波数积分记录文件
        if(needfstats){
            char *suffix = NULL;
            GRT_SAFE_ASPRINTF(&suffix, "_%04zu_%.5e", iw, grn->freqs[iw]);
            grt_KMET_init_fstats(grn->nr, grn->rs, grn->statsstr, suffix, local_Kmet);
            GRT_SAFE_FREE_PTR(suffix);
        }
        
        // 计算核函数过程中是否有遇到除零错误
        freq_invstats[iw]=GRT_INVERSE_SUCCESS;


        // ===================================================================================
        //                          Wavenumber Integration
        // 波数积分上限
        local_Kmet->kmax = hypot(local_Kmet->k0, local_Kmet->ampk * w / local_Kmet->vmin);
        K_INTEG *Kint = grt_wavenumber_integral(local_mstat, grn->nr, grn->rs, local_Kmet, grn->calc_upar, grt_kernel);

        // 记录到格林函数结构体内
        // 如果计算核函数过程中存在除零错误，则放弃该频率【通常在大震中距的低频段】
        recordin_GRN(iw, grn->nr, coef, Kint->sumJ, grn->u);
        if(grn->calc_upar){
            recordin_GRN(iw, grn->nr, coef, Kint->sumJz, grn->uiz);
            recordin_GRN(iw, grn->nr, coef, Kint->sumJr, grn->uir);
        }
        // ===================================================================================

        freq_invstats[iw] = local_mstat->stats;

        if(needfstats)  grt_KMET_destroy_fstats(grn->nr, local_Kmet);

        grt_free_mod1d_state(local_mstat);

        // 记录进度条变量 
        #pragma omp critical
        {
            progress++;
            if(print_progressbar) grt_printprogressBar("Computing Green Functions: ", progress*100/nf_valid);
        } 
        
        // Free allocated memory for temporary variables
        grt_free_K_INTEG(Kint);

    } // END omega loop

    // 打印 freq_invstats
    for(size_t iw = grn->nf1; iw <= grn->nf2; ++iw){
        if(freq_invstats[iw]==GRT_INVERSE_FAILURE){
            GRTRaiseWarning("iw=%zu, freq=%e(Hz), meet Zero Divison Error, results are filled with 0.\n", iw, grn->freqs[iw]);
        }
    }
    GRT_SAFE_FREE_PTR(freq_invstats);

    // 程序运行结束时间
    struct timeval end_t;
    gettimeofday(&end_t, NULL);
    if(print_progressbar) printf("Runtime: %.3f s\n", (end_t.tv_sec - begin_t.tv_sec) + (end_t.tv_usec - begin_t.tv_usec) / 1e6);
    fflush(stdout);

}






