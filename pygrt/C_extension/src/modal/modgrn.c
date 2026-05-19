/**
 * @file   modgrn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 根据频散结果计算面波解
 *   
 *                 
 */

#include <string.h>

#include "grt/modal/modgrn.h"
#include "grt/dynamic/source.h"
#include "grt/modal/secular.h"
#include "grt/modal/energy.h"


void grt_modsum_grn_Rayl(
    const MODEL1D_STATE *mstat, const cplx_t (*mod_potRaylLove_Down)[GRT_RAYL_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_RAYL_DIM],
    const cplx_t *egyint, const size_t iw, GRNSPEC *grn)
{
    // 注意模型中并未插入震源层和接收层
    size_t isrc = mstat->mod1d->isrc;
    size_t ircv = mstat->mod1d->ircv;

    real_t eigenK = mstat->k;

    // 震源层物性参数
    cplx_t src_mu = mstat->mu[isrc];
    real_t src_rho = mstat->mod1d->Rho[isrc];
    cplx_t src_a = eigenK * mstat->xa[isrc];
    cplx_t src_b = eigenK * mstat->xb[isrc];
    cplx_t src_Omg = eigenK*eigenK * (1.0 - 0.5*mstat->cbcb[isrc]);

    // 计算震源系数
    cplxChnlGrid src_coefD = {0};
    cplxChnlGrid src_coefU = {0};
    grt_source_coef_PSV(mstat, src_coefD, src_coefU);

    // 构造所需的震源系数
    cplx_t fcoef[GRT_SRC_M_NUM][GRT_RAYL_DIM];
    for(int i = 0; i < GRT_SRC_M_NUM; ++i){
        fcoef[i][0] =   eigenK*(src_coefD[i][0] - src_coefU[i][0]) - src_b*(src_coefD[i][1] + src_coefU[i][1]);
        fcoef[i][1] = - src_a *(src_coefD[i][0] + src_coefU[i][0]) + eigenK*(src_coefD[i][1] - src_coefU[i][1]);
        fcoef[i][2] =  2.0*src_mu*(src_Omg*(src_coefD[i][0] - src_coefU[i][0]) - eigenK*src_b*(src_coefD[i][1] + src_coefU[i][1]));
        fcoef[i][3] =  2.0*src_mu*(- eigenK*src_a*(src_coefD[i][0] + src_coefU[i][0]) + src_Omg*(src_coefD[i][1] - src_coefU[i][1]));
    }

    // 计算震源处的本征函数
    cplx_t T0[GRT_RAYL_DIM][GRT_RAYL_DIM] = {0};
    cplx_t src_eigenfn[4] = {0};
    grt_get_eigenfn_single_depth_Rayl(mstat, mod_potRaylLove_Down, mod_potRaylLove_Up, T0, false, mstat->mod1d->depsrc, isrc, src_eigenfn);

    // 构造公式中的 Fm^R 项
    cplx_t F_item[GRT_SRC_M_NUM] = {0};
    for(int i = 0; i < GRT_SRC_M_NUM; ++i){
        F_item[i] =   fcoef[i][0]*src_eigenfn[3] + fcoef[i][1]*src_eigenfn[2]
                    - fcoef[i][2]*src_eigenfn[1] - fcoef[i][3]*src_eigenfn[0];
    }

    // 计算接收处的本征函数
    cplx_t rcv_eigenfn[4] = {0};
    cplx_t rcv_eigenfn_z[4] = {0};
    grt_get_eigenfn_single_depth_Rayl(mstat, mod_potRaylLove_Down, mod_potRaylLove_Up, T0, false, mstat->mod1d->deprcv, ircv, rcv_eigenfn);
    if(grn->calc_upar){
        cplx_t xa = mstat->xa[ircv];
        cplx_t xb = mstat->xb[ircv];
        bool   isLiquid = mstat->mod1d->isLiquid[ircv];
        real_t k = mstat->k;

        cplx_t ak = k*k*xa;
        cplx_t bk = k*k*xb;
        cplx_t bb = xb*bk;
        cplx_t aa = xa*ak;
        cplx_t T0_z[GRT_RAYL_DIM][GRT_RAYL_DIM] = {
            {ak, bb, -ak, bb},
            {aa, bk, aa, -bk}
        };

        if(! isLiquid){
            grt_get_eigenfn_single_depth_Rayl(mstat, mod_potRaylLove_Down, mod_potRaylLove_Up, T0_z, true, mstat->mod1d->deprcv, ircv, rcv_eigenfn_z);
        }
    }

    // 计算该频点、该本征值下的频谱
    cplx_t i_m[GRT_MORDER_MAX+1] = {1.0, I, -1.0};
    cplx_t dom = 4.0*egyint[1] - 2.0*egyint[2]/eigenK;
    cplx_t wcoef, qcoef;
    cplx_t wcoefz, qcoefz;
    for(int im = 0; im < GRT_SRC_M_NUM; ++im){
        int modr = GRT_SRC_M_ORDERS[im];
        // Rayleigh 波仅考虑 q, w 项
        
        cplx_t aa = - 1.0/(4.0*src_rho*GRT_SQUARE(mstat->omega));

        // Z  
        wcoef = rcv_eigenfn[1] * i_m[modr] * F_item[im] / dom * aa;
        wcoefz = rcv_eigenfn_z[1] * i_m[modr] * F_item[im] / dom * aa;
        // R
        qcoef = I * rcv_eigenfn[0] * i_m[modr] * F_item[im] / dom * aa;
        qcoefz = I * rcv_eigenfn_z[0] * i_m[modr] * F_item[im] / dom * aa;

        for(size_t ir = 0; ir < grn->nr; ++ir){
            // 草稿推导有误？似乎应该使用 H^(2)(kr)
            cplx_t ekr = sqrt(2.0 / (PI*eigenK*grn->rs[ir])) * exp( - I * (eigenK*grn->rs[ir] + QUARTERPI));
            
            grn->u[ir][im][0][iw] += wcoef * ekr;
            grn->u[ir][im][1][iw] += - qcoef * ekr;
            if(grn->calc_upar){
                grn->uiz[ir][im][0][iw] += wcoefz * ekr;
                grn->uiz[ir][im][1][iw] += - qcoefz * ekr;

                cplx_t ekrr = - eigenK * ekr * I;
                grn->uir[ir][im][0][iw] += wcoef * ekrr;
                grn->uir[ir][im][1][iw] += - qcoef * ekrr;
            }
        }
    }

}



void grt_modsum_grn_Love(
    const MODEL1D_STATE *mstat, const cplx_t (*mod_potRaylLove_Down)[GRT_LOVE_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_LOVE_DIM],
    const cplx_t *egyint, const size_t iw, GRNSPEC *grn)
{
    // 注意模型中并未插入震源层和接收层
    size_t isrc = mstat->mod1d->isrc;
    size_t ircv = mstat->mod1d->ircv;

    real_t eigenK = mstat->k;

    // 震源层物性参数
    cplx_t src_mu = mstat->mu[isrc];
    real_t src_rho = mstat->mod1d->Rho[isrc];
    cplx_t src_b = eigenK * mstat->xb[isrc];

    // 计算震源系数
    cplxChnlGrid src_coefD = {0};
    cplxChnlGrid src_coefU = {0};
    grt_source_coef_SH(mstat, src_coefD, src_coefU);

    // 构造所需的震源系数
    cplx_t fcoef[GRT_SRC_M_NUM][GRT_LOVE_DIM] = {0};
    for(int i = 0; i < GRT_SRC_M_NUM; ++i){
        fcoef[i][0] =   eigenK*(src_coefD[i][2] - src_coefU[i][2]);
        fcoef[i][1] = - src_mu*eigenK*src_b*(src_coefD[i][2] + src_coefU[i][2]);
    }

    // 计算震源处的本征函数
    cplx_t T0[GRT_LOVE_DIM][GRT_LOVE_DIM] = {0};
    cplx_t src_eigenfn[4] = {0};
    grt_get_eigenfn_single_depth_Love(mstat, mod_potRaylLove_Down, mod_potRaylLove_Up, T0, false, mstat->mod1d->depsrc, isrc, src_eigenfn);

    // 构造公式中的 Fm^L 项
    cplx_t F_item[GRT_SRC_M_NUM] = {0};
    for(int i = 0; i < GRT_SRC_M_NUM; ++i){
        F_item[i] = fcoef[i][0]*src_eigenfn[1] - fcoef[i][1]*src_eigenfn[0];
    }

    // 计算接收处的本征函数
    cplx_t rcv_eigenfn[4] = {0};
    cplx_t rcv_eigenfn_z[4] = {0};
    grt_get_eigenfn_single_depth_Love(mstat, mod_potRaylLove_Down, mod_potRaylLove_Up, T0, false, mstat->mod1d->deprcv, ircv, rcv_eigenfn);
    if(grn->calc_upar){
        cplx_t xb = mstat->xb[ircv];
        bool   isLiquid = mstat->mod1d->isLiquid[ircv];
        real_t k = mstat->k;

        cplx_t bk = k*k*xb;
        cplx_t T0_z[GRT_LOVE_DIM][GRT_LOVE_DIM] = {{bk, -bk}};

        if(! isLiquid){
            grt_get_eigenfn_single_depth_Love(mstat, mod_potRaylLove_Down, mod_potRaylLove_Up, T0_z, true, mstat->mod1d->deprcv, ircv, rcv_eigenfn_z);
        }
    }

    // 计算该频点、该本征值下的频谱
    cplx_t i_m[GRT_MORDER_MAX+1] = {1.0, I, -1.0};
    cplx_t dom = 4.0*egyint[1];
    cplx_t vcoef;
    cplx_t vcoefz;
    for(int im = 0; im < GRT_SRC_M_NUM; ++im){
        int modr = GRT_SRC_M_ORDERS[im];
        // Love 波仅考虑 v 项
        if(modr == 0)  continue;
        
        cplx_t aa = - 1.0/(4.0*src_rho*GRT_SQUARE(mstat->omega));

        // T
        vcoef = I * rcv_eigenfn[0] * i_m[modr] * F_item[im] / dom * aa;
        vcoefz = I * rcv_eigenfn_z[0] * i_m[modr] * F_item[im] / dom * aa;

        for(size_t ir = 0; ir < grn->nr; ++ir){
            // 草稿推导有误？似乎应该使用 H^(2)(kr)
            cplx_t ekr = sqrt(2.0 / (PI*eigenK*grn->rs[ir])) * exp( - I * (eigenK*grn->rs[ir] + QUARTERPI));
            grn->u[ir][im][2][iw] += vcoef * ekr;
            if(grn->calc_upar){
                grn->uiz[ir][im][2][iw] += vcoefz * ekr;
                
                cplx_t ekrr = - eigenK * ekr * I;
                grn->uir[ir][im][2][iw] += vcoef * ekrr;
            }
        }
    }

}

void grt_modsum_grn(
    const MODEL1D_STATE *mstat, const DISPER_TYPE wtype, const size_t ncols, 
    const cplx_t (*mod_potRaylLove_Down)[ncols], const cplx_t (*mod_potRaylLove_Up)[ncols],
    const cplx_t *egyint, const size_t iw, GRNSPEC *grn)
{
    if(wtype == GRT_DISPERSION_RAYL && ncols == GRT_RAYL_DIM){
        grt_modsum_grn_Rayl(mstat, mod_potRaylLove_Down, mod_potRaylLove_Up, egyint, iw, grn);
    } 
    else if(wtype == GRT_DISPERSION_LOVE && ncols == GRT_LOVE_DIM) {
        grt_modsum_grn_Love(mstat, mod_potRaylLove_Down, mod_potRaylLove_Up, egyint, iw, grn);
    }
    else {
        GRTRaiseError("Wrong execution.");
    }

    // 对一些特殊情况的修正
    // 当震源和场点均位于地表时，可理论验证DS分量恒为0，这里直接赋0以避免后续的精度干扰
    if(mstat->mod1d->depsrc == 0.0 && mstat->mod1d->deprcv == 0.0)
    {
        for(size_t ir = 0; ir < grn->nr; ++ir){
            for(int c = 0; c < GRT_CHANNEL_NUM; ++c){
                grn->u[ir][GRT_SRC_M_DS_INDEX][c][iw] = 0.0;
                if(grn->calc_upar){
                    grn->uiz[ir][GRT_SRC_M_DS_INDEX][c][iw] = 0.0;
                    grn->uir[ir][GRT_SRC_M_DS_INDEX][c][iw] = 0.0;
                }
            }
        }
    }
}


void grt_modsum_grn_spec(MODEL1D *mod1d, const DISPER_TYPE wtype, EIGENFN_INFO *eigfnmet, GRNSPEC *grn)
{
    real_t fft_df = grn->freqs[1];
    size_t fft_offset = round(eigfnmet->freqs[0] / fft_df);

    size_t nlay = mod1d->n;

    // Rayleigh or Love
    const int ncols = (wtype == GRT_DISPERSION_RAYL)? GRT_RAYL_DIM : GRT_LOVE_DIM;

    // 循环所需要的频率
    #pragma omp parallel for schedule(guided) default(shared)
    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        MODEL1D_STATE *local_mstat = grt_init_mod1d_state(mod1d);

        real_t omega = PI2*eigfnmet->freqs[iw];
        EIGENV *eigv = eigfnmet->eigv + iw;
        
        grt_update_mod1d_state_omega(local_mstat, omega, true);
        
        // 假设这些本征值都是从 iref 层的久期函数算出来的，
        // 1. 先根据 iref（z_i+）的垂直波函数，计算出每个物理分界面上侧（z_j-）和下侧（z_j+）的垂直波函数
        // 2. 循环每个深度采样点，使用上侧（z_{j+1}-）和下侧（z_j+）的垂直波函数推导采样点上的垂直波函数以及本征函数
        for(size_t ic = 0; ic < eigv->n; ++ic){
            real_t cphase = eigv->c_roots[ic];
            
            local_mstat->c_phase = cphase;
            local_mstat->k = creal(omega)/cphase;
            grt_update_mod1d_state_k(local_mstat, local_mstat->k);

            size_t iref = eigv->c_roots_iref[ic];

            // 模型每个物理层介质下方 z_j+ 的垂直波函数
            cplx_t (*mod_potRaylLove_Down)[ncols] = (cplx_t (*)[ncols])calloc(nlay, sizeof(cplx_t)*ncols);
            // 模型每个物理层介质上方 z_j- 的垂直波函数，申请 n+1 的内存，方便后续不必讨论“半空间中的上行波场”
            cplx_t (*mod_potRaylLove_Up)[ncols] = (cplx_t (*)[ncols])calloc(nlay + 1, sizeof(cplx_t)*ncols);

            cplx_t potRaylLove[ncols];  memset(potRaylLove, 0, sizeof(cplx_t)*ncols);
            cplx_t potRaylLoveUp[ncols];  memset(potRaylLoveUp, 0, sizeof(cplx_t)*ncols);
            cplx_t secRaylLove = 0.0;
            grt_secular_function_potential(local_mstat, cphase, iref, eigfnmet->wtype, &secRaylLove, potRaylLove, potRaylLoveUp);
 
            // 计算出每个物理分界面上侧（z_{j+1}-）和 下侧（z_j+）的垂直波函数
            grt_get_mod_potential_Up_Down(local_mstat, eigfnmet->wtype, ncols, iref, potRaylLove, potRaylLoveUp, mod_potRaylLove_Down, mod_potRaylLove_Up);

            // 计算能量积分
            grt_energy_integrals(
                local_mstat, eigfnmet->wtype, ncols, 
                mod_potRaylLove_Down, mod_potRaylLove_Up, eigfnmet, &eigfnmet->eigfn[iw][ic]);

            // 计算面波格林函数
            grt_modsum_grn(
                local_mstat, eigfnmet->wtype, ncols, 
                mod_potRaylLove_Down, mod_potRaylLove_Up, eigfnmet->eigfn[iw][ic].egyint,
                iw + fft_offset, grn);
            
            GRT_SAFE_FREE_PTR(mod_potRaylLove_Down);
            GRT_SAFE_FREE_PTR(mod_potRaylLove_Up);
        }

        grt_free_mod1d_state(local_mstat);
    }
}