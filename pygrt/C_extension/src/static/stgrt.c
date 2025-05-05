/**
 * @file   stgrt.c
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

#include "static/stgrt.h"
#include "static/static_propagate.h"
#include "common/dwm.h"
#include "common/ptam.h"
#include "common/fim.h"
#include "common/safim.h"
#include "common/const.h"
#include "common/model.h"
#include "common/integral.h"
#include "common/search.h"



/**
 * 将计算好的复数形式的积分结果取实部记录到浮点数中
 * 
 * @param[in]    nr             震中距个数
 * @param[in]    coef           统一系数
 * @param[in]    sum_J          积分结果
 * @param[out]   grn            三分量结果，浮点数数组
 */
static void recordin_GRN(
    MYINT nr, MYCOMPLEX coef, MYCOMPLEX sum_J[nr][SRC_M_NUM][INTEG_NUM],
    MYREAL grn[nr][SRC_M_NUM][CHANNEL_NUM]
){
    // 局部变量，将某个频点的格林函数谱临时存放
    MYCOMPLEX (*tmp_grn)[SRC_M_NUM][CHANNEL_NUM] = (MYCOMPLEX(*)[SRC_M_NUM][CHANNEL_NUM])calloc(nr, sizeof(*tmp_grn));

    for(MYINT ir=0; ir<nr; ++ir){
        merge_Pk(sum_J[ir], tmp_grn[ir]);

        for(MYINT i=0; i<SRC_M_NUM; ++i) {
            for(MYINT c=0; c<CHANNEL_NUM; ++c){
                grn[ir][i][c] = creal(coef * tmp_grn[ir][i][c]);
            }

        }
    }

    free(tmp_grn);
}



void integ_static_grn(
    PYMODEL1D *pymod1d, MYINT nr, MYREAL *rs, MYREAL vmin_ref, MYREAL keps, MYREAL k0, MYREAL Length,
    MYREAL filonLength, MYREAL safilonTol, MYREAL filonCut, 

    // 返回值，代表Z、R、T分量
    MYREAL grn[nr][SRC_M_NUM][CHANNEL_NUM],

    bool calc_upar,
    MYREAL grn_uiz[nr][SRC_M_NUM][CHANNEL_NUM],
    MYREAL grn_uir[nr][SRC_M_NUM][CHANNEL_NUM],

    const char *statsstr // 积分结果输出
){
    MYREAL rmax=rs[findMinMax_MYREAL(rs, nr, true)];   // 最大震中距

    // pymod1d -> mod1d
    MODEL1D *mod1d = init_mod1d(pymod1d->n);
    get_mod1d(pymod1d, mod1d);

    const MYREAL hs = (fabs(pymod1d->depsrc - pymod1d->deprcv) < MIN_DEPTH_GAP_SRC_RCV)? 
                      MIN_DEPTH_GAP_SRC_RCV : fabs(pymod1d->depsrc - pymod1d->deprcv); // hs=max(震源和台站深度差,1.0)
    // 乘相应系数
    k0 *= PI/hs;

    if(vmin_ref < RZERO)  keps = -RONE;  // 若使用峰谷平均法，则不使用keps进行收敛判断

    MYREAL k=0.0;
    bool useFIM = (filonLength > RZERO) || (safilonTol > RZERO) ;    // 是否使用Filon积分（包括自适应Filon）
    const MYREAL dk=fabs(PI2/(Length*rmax));     // 波数积分间隔
    const MYREAL filondk = (filonLength > RZERO) ? PI2/(filonLength*rmax) : RZERO;  // Filon积分间隔
    const MYREAL filonK = filonCut/rmax;  // 波数积分和Filon积分的分割点

    const MYREAL kmax = k0;
    // 求和 sum F(ki,w)Jm(ki*r)ki 
    // 关于形状详见int_Pk()函数内的注释
    MYCOMPLEX (*sum_J)[SRC_M_NUM][INTEG_NUM] = (MYCOMPLEX(*)[SRC_M_NUM][INTEG_NUM])calloc(nr, sizeof(*sum_J));
    MYCOMPLEX (*sum_uiz_J)[SRC_M_NUM][INTEG_NUM] = (calc_upar)? (MYCOMPLEX(*)[SRC_M_NUM][INTEG_NUM])calloc(nr, sizeof(*sum_uiz_J)) : NULL;
    MYCOMPLEX (*sum_uir_J)[SRC_M_NUM][INTEG_NUM] = (calc_upar)? (MYCOMPLEX(*)[SRC_M_NUM][INTEG_NUM])calloc(nr, sizeof(*sum_uir_J)) : NULL;

    // 是否要输出积分过程文件
    bool needfstats = (statsstr!=NULL);

    // PTAM的积分中间结果, 每个震中距两个文件，因为PTAM对不同震中距使用不同的dk
    // 在文件名后加后缀，区分不同震中距
    char *ptam_fstatsdir[nr];
    for(MYINT ir=0; ir<nr; ++ir) {ptam_fstatsdir[ir] = NULL;}
    if(needfstats && vmin_ref < RZERO){
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
    
    // 创建波数积分记录文件
    FILE *fstats = NULL;
    // PTAM为每个震中距都创建波数积分记录文件
    FILE *(*ptam_fstatsnr)[2] = (FILE *(*)[2])malloc(nr * sizeof(*ptam_fstatsnr));
    {   
        MYINT len0 = (statsstr!=NULL) ? strlen(statsstr) : 0;
        char *fname = (char *)malloc((len0+200)*sizeof(char));
        if(needfstats){
            sprintf(fname, "%s/K", statsstr);
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
                sprintf(fname, "%s/K", ptam_fstatsdir[ir]);
                ptam_fstatsnr[ir][0] = fopen(fname, "wb");
                sprintf(fname, "%s/PTAM", ptam_fstatsdir[ir]);
                ptam_fstatsnr[ir][1] = fopen(fname, "wb");
            }
        }  
        free(fname);
    }

    // 计算核函数过程中是否有遇到除零错误
    //【静态解理论上不会有除零错误，这里是对应动态解的函数接口，作为一个占位符】
    MYINT inv_stats=INVERSE_SUCCESS;

    // 常规的波数积分
    k = discrete_integ(
        mod1d, dk, (useFIM)? filonK : kmax, keps, 0.0, nr, rs, 
        sum_J, calc_upar, sum_uiz_J, sum_uir_J,
        fstats, static_kernel, &inv_stats);
    
    // 基于线性插值的Filon积分
    if(useFIM){
        if(filondk > RZERO){
            // 基于线性插值的Filon积分，固定采样间隔
            k = linear_filon_integ(
                mod1d, k, dk, filondk, kmax, keps, 0.0, nr, rs, 
                sum_J, calc_upar, sum_uiz_J, sum_uir_J,
                fstats, static_kernel, &inv_stats);
        }
        else if(safilonTol > RZERO){
            // 基于自适应采样的Filon积分
            k = sa_filon_integ(
                mod1d, kmax, k, dk, safilonTol, kmax, 0.0, nr, rs, 
                sum_J, calc_upar, sum_uiz_J, sum_uir_J,
                fstats, static_kernel, &inv_stats);
        }
    }

    // k之后的部分使用峰谷平均法进行显式收敛，建议在浅源地震的时候使用   
    if(vmin_ref < RZERO){
        PTA_method(
            mod1d, k, dk, 0.0, nr, rs, 
            sum_J, calc_upar, sum_uiz_J, sum_uir_J,
            ptam_fstatsnr, static_kernel, &inv_stats);
    }


    
    MYCOMPLEX src_mu = (mod1d->lays + mod1d->isrc)->mu;
    MYCOMPLEX fac = dk * RONE/(RFOUR*PI * src_mu);
    
    // 将积分结果记录到浮点数数组中
    recordin_GRN(nr, fac, sum_J, grn);
    if(calc_upar){
        recordin_GRN(nr, fac, sum_uiz_J, grn_uiz);
        recordin_GRN(nr, fac, sum_uir_J, grn_uir);
    }


    // Free allocated memory for temporary variables
    free(sum_J);
    if(sum_uiz_J) free(sum_uiz_J);
    if(sum_uir_J) free(sum_uir_J);

    free_mod1d(mod1d);

    for(MYINT ir=0; ir<nr; ++ir){
        if(ptam_fstatsdir[ir]!=NULL){
            free(ptam_fstatsdir[ir]);
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

    free(ptam_fstatsnr);
}