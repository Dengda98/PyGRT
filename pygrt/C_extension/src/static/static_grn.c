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
#include "grt/static/static_propagate.h"
#include "grt/common/dwm.h"
#include "grt/common/ptam.h"
#include "grt/common/fim.h"
#include "grt/common/safim.h"
#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/common/integral.h"
#include "grt/common/search.h"



/**
 * 将计算好的复数形式的积分结果取实部记录到浮点数中
 * 
 * @param[in]    nr             震中距个数
 * @param[in]    coef           统一系数
 * @param[in]    sum_J          积分结果
 * @param[out]   grn            三分量结果，浮点数数组
 */
static void recordin_GRN(
    MYINT nr, MYCOMPLEX coef, MYCOMPLEX sum_J[nr][GRT_SRC_M_NUM][GRT_INTEG_NUM],
    MYREAL grn[nr][GRT_SRC_M_NUM][GRT_CHANNEL_NUM]
){
    // 局部变量，将某个频点的格林函数谱临时存放
    MYCOMPLEX (*tmp_grn)[GRT_SRC_M_NUM][GRT_CHANNEL_NUM] = (MYCOMPLEX(*)[GRT_SRC_M_NUM][GRT_CHANNEL_NUM])calloc(nr, sizeof(*tmp_grn));

    for(MYINT ir=0; ir<nr; ++ir){
        grt_merge_Pk(sum_J[ir], tmp_grn[ir]);

        for(MYINT i=0; i<GRT_SRC_M_NUM; ++i) {
            for(MYINT c=0; c<GRT_CHANNEL_NUM; ++c){
                grn[ir][i][c] = creal(coef * tmp_grn[ir][i][c]);
            }

        }
    }

    GRT_SAFE_FREE_PTR(tmp_grn);
}



void grt_integ_static_grn(
    GRT_MODEL1D *mod1d, MYINT nr, MYREAL *rs, MYREAL vmin_ref, MYREAL keps, MYREAL k0, MYREAL Length,
    MYREAL filonLength, MYREAL safilonTol, MYREAL filonCut, 

    // 返回值，代表Z、R、T分量
    MYREAL grn[nr][GRT_SRC_M_NUM][GRT_CHANNEL_NUM],

    bool calc_upar,
    MYREAL grn_uiz[nr][GRT_SRC_M_NUM][GRT_CHANNEL_NUM],
    MYREAL grn_uir[nr][GRT_SRC_M_NUM][GRT_CHANNEL_NUM],

    const char *statsstr // 积分结果输出
){
    MYREAL rmax=rs[grt_findMax_MYREAL(rs, nr)];   // 最大震中距

    const MYREAL hs = GRT_MAX(fabs(mod1d->depsrc - mod1d->deprcv), GRT_MIN_DEPTH_GAP_SRC_RCV); // hs=max(震源和台站深度差,1.0)

    // 乘相应系数
    k0 *= PI/hs;

    if(vmin_ref < 0.0)  keps = -1.0;  // 若使用峰谷平均法，则不使用keps进行收敛判断

    MYREAL k=0.0;
    bool useFIM = (filonLength > 0.0) || (safilonTol > 0.0) ;    // 是否使用Filon积分（包括自适应Filon）
    const MYREAL dk=fabs(PI2/(Length*rmax));     // 波数积分间隔
    const MYREAL filondk = (filonLength > 0.0) ? PI2/(filonLength*rmax) : 0.0;  // Filon积分间隔
    const MYREAL filonK = filonCut/rmax;  // 波数积分和Filon积分的分割点

    const MYREAL kmax = k0;
    // 求和 sum F(ki,w)Jm(ki*r)ki 
    // 关于形状详见int_Pk()函数内的注释
    MYCOMPLEX (*sum_J)[GRT_SRC_M_NUM][GRT_INTEG_NUM] = (MYCOMPLEX(*)[GRT_SRC_M_NUM][GRT_INTEG_NUM])calloc(nr, sizeof(*sum_J));
    MYCOMPLEX (*sum_uiz_J)[GRT_SRC_M_NUM][GRT_INTEG_NUM] = (calc_upar)? (MYCOMPLEX(*)[GRT_SRC_M_NUM][GRT_INTEG_NUM])calloc(nr, sizeof(*sum_uiz_J)) : NULL;
    MYCOMPLEX (*sum_uir_J)[GRT_SRC_M_NUM][GRT_INTEG_NUM] = (calc_upar)? (MYCOMPLEX(*)[GRT_SRC_M_NUM][GRT_INTEG_NUM])calloc(nr, sizeof(*sum_uir_J)) : NULL;

    // 是否要输出积分过程文件
    bool needfstats = (statsstr!=NULL);

    // PTAM的积分中间结果, 每个震中距两个文件，因为PTAM对不同震中距使用不同的dk
    // 在文件名后加后缀，区分不同震中距
    char **ptam_fstatsdir = (char**)calloc(nr, sizeof(char*));
    if(needfstats && vmin_ref < 0.0){
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
    
    // 创建波数积分记录文件
    FILE *fstats = NULL;
    // PTAM为每个震中距都创建波数积分记录文件
    FILE *(*ptam_fstatsnr)[2] = (FILE *(*)[2])malloc(nr * sizeof(*ptam_fstatsnr));
    {   
        char *fname = NULL;
        if(needfstats){
            GRT_SAFE_ASPRINTF(&fname, "%s/K", statsstr);
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
                GRT_SAFE_ASPRINTF(&fname, "%s/K", ptam_fstatsdir[ir]);
                ptam_fstatsnr[ir][0] = fopen(fname, "wb");
                GRT_SAFE_ASPRINTF(&fname, "%s/PTAM", ptam_fstatsdir[ir]);
                ptam_fstatsnr[ir][1] = fopen(fname, "wb");
            }
        }  
        GRT_SAFE_FREE_PTR(fname);
    }

    // 计算核函数过程中是否有遇到除零错误
    //【静态解理论上不会有除零错误，这里是对应动态解的函数接口，作为一个占位符】
    MYINT inv_stats=GRT_INVERSE_SUCCESS;

    // 常规的波数积分
    k = grt_discrete_integ(
        mod1d, dk, (useFIM)? filonK : kmax, keps, 0.0, nr, rs, 
        sum_J, calc_upar, sum_uiz_J, sum_uir_J,
        fstats, grt_static_kernel, &inv_stats);
    
    // 基于线性插值的Filon积分
    if(useFIM){
        if(filondk > 0.0){
            // 基于线性插值的Filon积分，固定采样间隔
            k = grt_linear_filon_integ(
                mod1d, k, dk, filondk, kmax, keps, 0.0, nr, rs, 
                sum_J, calc_upar, sum_uiz_J, sum_uir_J,
                fstats, grt_static_kernel, &inv_stats);
        }
        else if(safilonTol > 0.0){
            // 基于自适应采样的Filon积分
            k = grt_sa_filon_integ(
                mod1d, kmax, k, dk, safilonTol, kmax, 0.0, nr, rs, 
                sum_J, calc_upar, sum_uiz_J, sum_uir_J,
                fstats, grt_static_kernel, &inv_stats);
        }
    }

    // k之后的部分使用峰谷平均法进行显式收敛，建议在浅源地震的时候使用   
    if(vmin_ref < 0.0){
        grt_PTA_method(
            mod1d, k, dk, 0.0, nr, rs, 
            sum_J, calc_upar, sum_uiz_J, sum_uir_J,
            ptam_fstatsnr, grt_static_kernel, &inv_stats);
    }


    
    MYCOMPLEX src_mu = mod1d->mu[mod1d->isrc];
    MYCOMPLEX fac = dk * 1.0/(4.0*PI * src_mu);
    
    // 将积分结果记录到浮点数数组中
    recordin_GRN(nr, fac, sum_J, grn);
    if(calc_upar){
        recordin_GRN(nr, fac, sum_uiz_J, grn_uiz);
        recordin_GRN(nr, fac, sum_uir_J, grn_uir);
    }


    // Free allocated memory for temporary variables
    GRT_SAFE_FREE_PTR(sum_J);
    GRT_SAFE_FREE_PTR(sum_uiz_J);
    GRT_SAFE_FREE_PTR(sum_uir_J);

    GRT_SAFE_FREE_PTR_ARRAY(ptam_fstatsdir, nr);

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
}