/**
 * @file   integ_method.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-12
 * 
 * 使用结构体来管理与波数积分相关的参数和方法
 * 
 */

#include "grt/common/model.h"
#include "grt/integral/ptam.h"
#include "grt/integral/fim.h"
#include "grt/integral/safim.h"
#include "grt/integral/dwm.h"
#include "grt/integral/dcm.h"
#include "grt/integral/integ_method.h"



void grt_KMET_init_fstats(
    const size_t nr, const real_t *rs, 
    const char *statsstr, const char *suffix, K_INTEG_METHOD *Kmet)
{
    if(statsstr == NULL)  return;

    if(Kmet->fstats != NULL){
        GRTRaiseError("fstats != NULL in K_INTEG_METHOD.");
    }

    // 为当前频率创建波数积分记录文件
    // PTAM为每个震中距都创建波数积分记录文件
    Kmet->ptam_fstatsnr = (FILE *(*)[2])calloc(nr, sizeof(*Kmet->ptam_fstatsnr));
    char *fname = NULL;
    GRT_SAFE_ASPRINTF(&fname, "%s/K%s", statsstr, suffix);
    Kmet->fstats = fopen(fname, "wb");

    // PTAM的积分中间结果, 每个震中距两个文件，因为PTAM对不同震中距使用不同的dk
    // 在文件名后加后缀，区分不同震中距
    char *ptam_dirname = NULL;
    if(Kmet->applyPTAM){
        for(size_t ir = 0; ir < nr; ++ir){
            // 新建文件夹目录 
            GRT_SAFE_ASPRINTF(&ptam_dirname, "%s/PTAM_%04zu_%.5e", statsstr, ir, rs[ir]);
            GRTCheckMakeDir(ptam_dirname);

            Kmet->ptam_fstatsnr[ir][0] = Kmet->ptam_fstatsnr[ir][1] = NULL;
            GRT_SAFE_ASPRINTF(&fname, "%s/K%s", ptam_dirname, suffix);
            Kmet->ptam_fstatsnr[ir][0] = fopen(fname, "wb");
            GRT_SAFE_ASPRINTF(&fname, "%s/PTAM%s", ptam_dirname, suffix);
            Kmet->ptam_fstatsnr[ir][1] = fopen(fname, "wb");
        }
    }
    
    GRT_SAFE_FREE_PTR(fname);
    GRT_SAFE_FREE_PTR(ptam_dirname);
}


void grt_KMET_destroy_fstats(const size_t nr, K_INTEG_METHOD *Kmet)
{
    if(Kmet->fstats!=NULL) fclose(Kmet->fstats);

    if(Kmet->ptam_fstatsnr != NULL){
        for(size_t ir=0; ir<nr; ++ir){
            if(Kmet->ptam_fstatsnr[ir][0]!=NULL){
                fclose(Kmet->ptam_fstatsnr[ir][0]);
            }
            if(Kmet->ptam_fstatsnr[ir][1]!=NULL){
                fclose(Kmet->ptam_fstatsnr[ir][1]);
            }
        }
        GRT_SAFE_FREE_PTR(Kmet->ptam_fstatsnr);
    }

}




K_INTEG * grt_wavenumber_integral(
    GRT_MODEL1D *mod1d, size_t nr, real_t *rs, K_INTEG_METHOD *Kmet, bool calc_upar, GRT_KernelFunc kerfunc)
{
    real_t k = 0.0;

    real_t kcut = Kmet->kmax;

    bool isFilon = Kmet->applyFIM || Kmet->applySAFIM;

    if(isFilon){
        kcut = GRT_MIN(Kmet->kcut, Kmet->kmax);
    }

    // 求和 sum F(ki,w)Jm(ki*r)ki 
    K_INTEG *Kint = grt_init_K_INTEG(calc_upar, nr);
    Kint->applyDCM = Kmet->applyDCM;
    Kint->kmax = Kmet->kmax;

    // 准备 DCM，计算波数上限处的核函数
    if(Kint->applyDCM){
        kerfunc(mod1d, Kint->kmax, Kint->QWV_kmax, Kint->calc_upar, Kint->QWVz_kmax);
    }
    
    // DWM
    k = grt_discrete_integ(
        mod1d, Kmet->dk, kcut, Kmet->keps, nr, rs, 
        Kint, Kmet->fstats, kerfunc);
    if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

    // 基于线性插值的Filon积分，固定采样间隔
    if(Kmet->applyFIM){
        k = grt_linear_filon_integ(
            mod1d, k, Kmet->dk, Kmet->filondk, Kmet->kmax, Kmet->keps, nr, rs, 
            Kint, Kmet->fstats, kerfunc);
        if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
    }
    // 基于自适应采样的Filon积分
    else if(Kmet->applySAFIM){
        real_t kref = (creal(mod1d->omega) > 0.0) ? creal(mod1d->omega) / Kmet->vmin * Kmet->ampk : 0.0; // 静态解中频率位置给了负值
        k = grt_sa_filon_integ(
            mod1d, k, Kmet->dk, Kmet->sa_tol, Kmet->kmax, kref, nr, rs, 
            Kint, Kmet->fstats, kerfunc);
        if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
    }

    // 显式收敛
    if(Kmet->applyPTAM){
        grt_PTA_method(
            mod1d, k, Kmet->dk, nr, rs, 
            Kint, Kmet->ptam_fstatsnr, kerfunc);
        if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
    }
    else if(Kmet->applyDCM){
        grt_dcm_correction(nr, rs, Kint, !isFilon);
    }


    BEFORE_RETURN:

    return Kint;
}
