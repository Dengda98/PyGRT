/**
 * @file   grnspec.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-04
 * 
 * 将格林函数频谱用一个结构体包裹
 *   
 *                 
 */

#include <string.h>

#include "grt/dynamic/grnspec.h"
#include "grt/common/myfftw.h"
#include "grt/common/sacio.h"

void grt_grnspec_allocate_u(GRNSPEC *grn)
{
    grn->u = (pcplxChnlGrid *) calloc(grn->nr, sizeof(*grn->u));
    grn->uiz = (grn->calc_upar)? (pcplxChnlGrid *) calloc(grn->nr, sizeof(*grn->uiz)) : NULL;
    grn->uir = (grn->calc_upar)? (pcplxChnlGrid *) calloc(grn->nr, sizeof(*grn->uir)) : NULL;

    for(size_t ir = 0; ir < grn->nr; ++ir){
        GRT_LOOP_ChnlGrid(im, c){
            grn->u[ir][im][c] = (cplx_t*) calloc(grn->nf, sizeof(cplx_t));
            if(grn->uiz)  grn->uiz[ir][im][c] = (cplx_t*)calloc(grn->nf, sizeof(cplx_t));
            if(grn->uir)  grn->uir[ir][im][c] = (cplx_t*)calloc(grn->nf, sizeof(cplx_t));
        }
    }
}

void grt_grnspec_free_u(GRNSPEC *grn)
{
    for(size_t ir = 0; ir < grn->nr; ++ir){
        GRT_LOOP_ChnlGrid(im, c){
            GRT_SAFE_FREE_PTR(grn->u[ir][im][c]);
            if(grn->uiz) GRT_SAFE_FREE_PTR(grn->uiz[ir][im][c]);
            if(grn->uir) GRT_SAFE_FREE_PTR(grn->uir[ir][im][c]);
        }
    }
    GRT_SAFE_FREE_PTR(grn->u);
    GRT_SAFE_FREE_PTR(grn->uiz);
    GRT_SAFE_FREE_PTR(grn->uir);
}



/** 将一条数据反变换回时间域再进行处理，保存到SAC文件 */
static void write_one_to_sac(
    const char *srcname, const char ch, GRT_FFTW_HOLDER *fh, const real_t wI, 
    SACTRACE *sac, const char *s_output_subdir, const char *s_prefix,
    const int sgn, bool skipImagComps, const cplx_t *grncplx)
{
    snprintf(sac->hd.kcmpnm, sizeof(sac->hd.kcmpnm), "%s%s%c", s_prefix, srcname, ch);
    
    char *s_outpath = NULL;
    GRT_SAFE_ASPRINTF(&s_outpath, "%s/%s.sac", s_output_subdir, sac->hd.kcmpnm);

    // 执行fft任务会修改数组，需重新置零
    grt_reset_fftw_holder_zero(fh);
    
    // 赋值复数，包括时移
    cplx_t cfac, ccoef;
    cfac = exp(I*PI2*fh->df*sac->hd.b);
    ccoef = sgn;
    // 只赋值有效长度，其余均为0
    for(size_t i = 0; i < fh->nf_valid; ++i){
        fh->W_f[i] = grncplx[i] * ccoef;
        ccoef *= cfac;
    }


    if(! fh->naive_inv){
        // 发起fft任务 
        fftw_execute(fh->plan);
    } else {
        grt_naive_inverse_transform_double(fh);
    }

    // 归一化，并处理虚频
    // 并转为 SAC 需要的单精度类型
    real_t fac, coef;
    coef = fh->df;
    fac = 1.0;
    if (! skipImagComps) {
        coef *= exp(sac->hd.b * wI);
        fac = exp(wI*fh->dt);
    }
    for(size_t i = 0; i < fh->nt; ++i){
        sac->data[i] = fh->w_t[i] * coef;
        coef *= fac;
    }
    

    // 以sac文件保存到本地
    grt_write_SACTRACE(s_outpath, sac);

    GRT_SAFE_FREE_PTR(s_outpath);
}


void grt_grnspec_write_sac(
    const GRNSPEC *grn, const real_t (*travtPS)[2], const real_t *begintimes, char **outputdirs, GRT_FFTW_HOLDER *fh, SACTRACE *sac, 
    const char *validChnls, const bool skipImagComps, const bool saveEX, const bool saveVF, const bool saveHF, const bool saveDC)
{
    // 做反傅里叶变换，保存SAC文件
    for(size_t ir = 0; ir < grn->nr; ++ir){
        real_t dist = grn->rs[ir];
        sac->hd.dist = dist;

        // 文件保存子目录
        const char *dirpath = outputdirs[ir];
        GRTCheckMakeDir(dirpath);

        // 计算理论走时
        sac->hd.t0 = travtPS[ir][0];
        strcpy(sac->hd.kt0, "P");
        sac->hd.t1 = travtPS[ir][1];
        strcpy(sac->hd.kt1, "S");

        // 时间延迟
        sac->hd.b = begintimes[ir];

        GRT_LOOP_ChnlGrid(im, c){
            if(strchr(validChnls, GRT_ZRT_CODES[c]) == NULL)  continue;

            if(! saveEX && im==GRT_SRC_M_EX_INDEX)  continue;
            if(! saveVF && im==GRT_SRC_M_VF_INDEX)  continue;
            if(! saveHF && im==GRT_SRC_M_HF_INDEX)  continue;
            if(! saveDC && im>=GRT_SRC_M_DD_INDEX)  continue;

            int modr = GRT_SRC_M_ORDERS[im];
            int sgn=1;  // 用于反转Z分量

            if(modr==0 && GRT_ZRT_CODES[c]=='T')  continue;  // 跳过输出0阶的T分量

            // Z分量反转
            sgn = (GRT_ZRT_CODES[c]=='Z') ? -1 : 1;

            char ch = GRT_ZRT_CODES[c];

            write_one_to_sac(GRT_SRC_M_NAME_ABBR[im], ch, fh, grn->wI, sac, dirpath, "", sgn, skipImagComps, grn->u[ir][im][c]);

            if(grn->calc_upar){
                write_one_to_sac(GRT_SRC_M_NAME_ABBR[im], ch, fh, grn->wI, sac, dirpath, "z", sgn*(-1), skipImagComps, grn->uiz[ir][im][c]);
                write_one_to_sac(GRT_SRC_M_NAME_ABBR[im], ch, fh, grn->wI, sac, dirpath, "r", sgn     , skipImagComps, grn->uir[ir][im][c]);
            }
        }
    }
}