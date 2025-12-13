/**
 * @file   static_grn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 * 以下代码实现的是 广义反射透射系数矩阵+离散波数法 计算静态格林函数，参考：
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
#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/integral/integ_method.h"
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
    GRT_MODEL1D *mod1d, size_t nr, real_t *rs, K_INTEG_METHOD *Kmet,
    bool calc_upar, 
    realChnlGrid grn[nr],
    realChnlGrid grn_uiz[nr],
    realChnlGrid grn_uir[nr],
    const char *statsstr) 
{
    // 是否要输出积分过程文件
    bool needfstats = (statsstr!=NULL);

    // 输出积分过程文件
    if(needfstats) grt_KMET_init_fstats(nr, rs, statsstr, "", Kmet);

    // ===================================================================================
    //                          Wavenumber Integration
    // 波数积分上限
    Kmet->kmax = Kmet->k0;
    K_INTEG *Kint = grt_wavenumber_integral(mod1d, nr, rs, Kmet, calc_upar, grt_static_kernel);
    
    cplx_t src_mu = mod1d->mu[mod1d->isrc];
    cplx_t fac = Kmet->dk * 1.0/(4.0*PI * src_mu);
    
    // 将积分结果记录到浮点数数组中
    recordin_GRN(nr, fac, Kint->sumJ, grn);
    if(calc_upar){
        recordin_GRN(nr, fac, Kint->sumJz, grn_uiz);
        recordin_GRN(nr, fac, Kint->sumJr, grn_uir);
    }
    // ===================================================================================

    // Free allocated memory for temporary variables
    grt_free_K_INTEG(Kint);

    if(needfstats)  grt_KMET_destroy_fstats(nr, Kmet);

}