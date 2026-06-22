/**
 * @file   static_grn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 * 以下代码实现的是 广义反射透射系数矩阵 计算静态格林函数，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *              透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 * 
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>

#include "grt/static/static_grn.h"
#include "grt/static/static_util.h"
#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/integral/integ_process.h"
#include "grt/common/search.h"


/**
 * 将计算好的复数形式的积分结果取实部记录到浮点数中
 * 
 * @param[in]    nr             震中距个数
 * @param[in]    coef           统一系数
 * @param[in]    sumJ           积分结果
 * @param[out]   grn            三分量结果，浮点数数组
 */
static void recordin_GRN(
    size_t nr, cplx_t coef, cplxIntegGrid sumJ[nr],
    realChnlGrid grn[nr])
{
    // 局部变量，将某个频点的格林函数谱临时存放
    cplxChnlGrid *tmp_grn = (cplxChnlGrid *)calloc(nr, sizeof(*tmp_grn));

    for(size_t ir=0; ir<nr; ++ir){
        grt_merge_Pk(sumJ[ir], tmp_grn[ir]);

        GRT_LOOP_ChnlGrid(im, c){
            int modr = GRT_SRC_M_ORDERS[im];
            if(modr == 0 && GRT_ZRT_CODES[c] == 'T')  continue;

            grn[ir][im][c] = creal(coef * tmp_grn[ir][im][c]);
        }

    }

    GRT_SAFE_FREE_PTR(tmp_grn);
}



void grt_integ_static_grn(
    MODEL1D *mod1d, size_t nr, real_t *rs, K_INTEG_PROCESS *Kproc,
    bool calc_upar, 
    realChnlGrid grn[nr],
    realChnlGrid grn_uiz[nr],
    realChnlGrid grn_uir[nr],
    const char *statsstr) 
{
    // 根据静态解传入的震中距形式，有必要对震中距去重，从而减少不必要的计算量
    size_t uniq_nr = nr;
    real_t *uniq_rs = (real_t *)calloc(nr, sizeof(real_t));
    size_t *rs_refidx = (size_t *)calloc(nr, sizeof(size_t));
    memcpy(uniq_rs, rs, sizeof(real_t)*nr);
    uniq_nr = 1;
    for(size_t ir = 1; ir < nr; ++ir){
        size_t jr = 0;
        for( ; jr < uniq_nr; ++jr){
            // 找到重复，跳过
            if(fabs(uniq_rs[jr] - rs[ir]) < 1e-8) break;
        }
        // 无重复
        if(jr == uniq_nr){
            uniq_rs[uniq_nr++] = rs[ir];
        }
        rs_refidx[ir] = jr;
    }

    // 是否要输出积分过程文件
    bool needfstats = (statsstr!=NULL);

    // 输出积分过程文件
    if(needfstats) grt_KPROC_init_fstats(uniq_nr, uniq_rs, statsstr, "", Kproc);

    // ===================================================================================
    //                          Wavenumber Integration
    // 波数积分上限
    if(Kproc->fixed_k0){
        Kproc->kmax = Kproc->k0;
    } else {
        size_t ncount = 0, nk = 0;
        real_t static_kmax = 0.0;
        static_kmax = grt_predict_static_kmax(mod1d, Kproc->k0, &ncount);
        nk = (size_t)(static_kmax / Kproc->dk);
        GRTRaiseInfo(
            "for a proper k0, ncount = %zu, k0 = %.3e, k0_ref = %.3e, nk = %zu", 
            ncount, static_kmax, Kproc->k0, nk);
        
        // 若 nk 不够，适当调整 dk
        if(nk < GRT_MIN_STATIC_NK){
            real_t new_dk = static_kmax / GRT_MIN_STATIC_NK;
            GRTRaiseWarning("nk(%zu) < %d, adjust dk from %.3e to %.3e", nk, GRT_MIN_STATIC_NK, Kproc->dk, new_dk);
            Kproc->dk = new_dk;
        }
        
        if(Kproc->cvgmet == K_INTEG_CONVERG_AUTO){
            // 如果上限触发边界，且未指定收敛算法，则强制使用 DCM 进行收敛
            if(static_kmax >= Kproc->k0){
                Kproc->cvgmet = K_INTEG_CONVERG_DCM;
                Kproc->keps = 0.0;
                GRTRaiseWarning("k0 reaches k0_ref, apply %s.", GRT_EXPLAIN_CVGMETHOD(Kproc->cvgmet));
            }
        } else if(Kproc->cvgmet != K_INTEG_CONVERG_REFUSE) {
            // 正常打印手动选择的收敛方法
            GRTRaiseInfo("manually set the %s.", GRT_EXPLAIN_CVGMETHOD(Kproc->cvgmet));
        }

        Kproc->kmax = static_kmax;
    }


    // 模型状态
    MODEL1D_STATE *mstat = grt_init_mod1d_state(mod1d);
    grt_update_mod1d_state_omega(mstat, 1.0, true);
    
    K_INTEG *Kint = grt_wavenumber_integral(mstat, uniq_nr, uniq_rs, Kproc, calc_upar, grt_static_kernel);
    
    cplx_t src_mu = mstat->mu[mod1d->isrc];
    cplx_t fac = Kproc->dk * 1.0/(4.0*PI * src_mu);
    
    // 将积分结果记录到浮点数数组中
    recordin_GRN(uniq_nr, fac, Kint->sumJ, grn);
    if(calc_upar){
        recordin_GRN(uniq_nr, fac, Kint->sumJz, grn_uiz);
        recordin_GRN(uniq_nr, fac, Kint->sumJr, grn_uir);
    }
    // ===================================================================================

    // 从去重震中距恢复结果
    if(nr != uniq_nr){
        for(size_t ir = 0; ir < nr; ++ir){
            size_t ir_inv = nr-ir-1;
            size_t idx = rs_refidx[ir_inv];
            GRT_LOOP_ChnlGrid(im, c){
                grn[ir_inv][im][c] = grn[idx][im][c];
                if(calc_upar){
                    grn_uiz[ir_inv][im][c] = grn_uiz[idx][im][c];
                    grn_uir[ir_inv][im][c] = grn_uir[idx][im][c];
                }
            }
        }
    }

    // Free allocated memory for temporary variables
    grt_free_K_INTEG(Kint);
    grt_free_mod1d_state(mstat);

    if(needfstats)  grt_KPROC_destroy_fstats(uniq_nr, Kproc);

    GRT_SAFE_FREE_PTR(uniq_rs);
    GRT_SAFE_FREE_PTR(rs_refidx);
}