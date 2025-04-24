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
#include "common/const.h"
#include "common/model.h"
#include "common/integral.h"
#include "common/search.h"



/**
 * 将计算好的复数形式的积分结果取实部记录到浮点数中
 * 
 * @param    iw      (in)当前频率索引值
 * @param    nr      (in)震中距个数
 * @param    coef    (in)统一系数
 * @param  sum_EXP_J[nr][3][4]  (in)爆炸源
 * @param  sum_VF_J[nr][3][4]   (in)垂直力源
 * @param  sum_HF_J[nr][3][4]   (in)水平力源
 * @param  sum_DC_J[nr][3][4]   (in)剪切源
 * @param      EXPgrn[nr][2]      (out)浮点数数组，爆炸源的Z、R分量频谱结果
 * @param      VFgrn[nr][2]       (out)浮点数数组，垂直力源的Z、R分量频谱结果
 * @param      HFgrn[nr][3]       (out)浮点数数组，水平力源的Z、R、T分量频谱结果
 * @param      DDgrn[nr][2]       (out)浮点数数组，45度倾滑的Z、R分量频谱结果
 * @param      DSgrn[nr][3]       (out)浮点数数组，90度倾滑的Z、R、T分量频谱结果
 * @param      SSgrn[nr][3]       (out)浮点数数组，90度走滑的Z、R、T分量频谱结果
 */
static void recordin_GRN(
    MYINT nr, MYCOMPLEX coef, 
    MYCOMPLEX sum_EXP_J[nr][3][4], MYCOMPLEX sum_VF_J[nr][3][4],  
    MYCOMPLEX sum_HF_J[nr][3][4],  MYCOMPLEX sum_DC_J[nr][3][4],  
    MYREAL EXPgrn[nr][2], MYREAL VFgrn[nr][2], MYREAL HFgrn[nr][3],
    MYREAL DDgrn[nr][2],  MYREAL DSgrn[nr][3], MYREAL SSgrn[nr][3]
){
    // 局部变量，将某个频点的格林函数谱临时存放
    MYCOMPLEX (*tmp_EXP)[2] = (MYCOMPLEX(*)[2])calloc(nr, sizeof(*tmp_EXP));
    MYCOMPLEX (*tmp_VF)[2] = (MYCOMPLEX(*)[2])calloc(nr, sizeof(*tmp_VF));
    MYCOMPLEX (*tmp_HF)[3] = (MYCOMPLEX(*)[3])calloc(nr, sizeof(*tmp_HF));
    MYCOMPLEX (*tmp_DD)[2] = (MYCOMPLEX(*)[2])calloc(nr, sizeof(*tmp_DD));
    MYCOMPLEX (*tmp_DS)[3] = (MYCOMPLEX(*)[3])calloc(nr, sizeof(*tmp_DS));
    MYCOMPLEX (*tmp_SS)[3] = (MYCOMPLEX(*)[3])calloc(nr, sizeof(*tmp_SS));

    for(MYINT ir=0; ir<nr; ++ir){
        for(MYINT ii=0; ii<3; ++ii){
            if(ii<2){
                tmp_EXP[ir][ii] = RZERO;
                tmp_VF[ir][ii] = RZERO; 
                tmp_DD[ir][ii] = RZERO;
            }
            tmp_HF[ir][ii] = RZERO;
            tmp_DS[ir][ii] = RZERO;
            tmp_SS[ir][ii] = RZERO;
        }
        merge_Pk(
            (sum_EXP_J)? sum_EXP_J[ir] : NULL, 
            (sum_VF_J)?  sum_VF_J[ir]  : NULL, 
            (sum_HF_J)?  sum_HF_J[ir]  : NULL, 
            (sum_DC_J)?  sum_DC_J[ir]  : NULL, 
            tmp_EXP[ir], tmp_VF[ir],  tmp_HF[ir], 
            tmp_DD[ir], tmp_DS[ir], tmp_SS[ir]);

        MYCOMPLEX mtmp;
        for(MYINT ii=0; ii<3; ++ii) {
            if(ii<2){
                if(EXPgrn!=NULL){
                    mtmp = coef*tmp_EXP[ir][ii]; // m=0 爆炸源
                    EXPgrn[ir][ii] = CREAL(mtmp);
                }
                if(VFgrn!=NULL){
                    mtmp = coef*tmp_VF[ir][ii]; // m=0 垂直力源
                    VFgrn[ir][ii] = CREAL(mtmp);
                }
                if(DDgrn!=NULL){
                    mtmp = coef*tmp_DD[ir][ii]; // m=0 45-倾滑
                    DDgrn[ir][ii] = CREAL(mtmp);
                }
            }
            if(HFgrn!=NULL){
                mtmp = coef*tmp_HF[ir][ii]; // m=1 水平力源
                HFgrn[ir][ii] = CREAL(mtmp);
            }
            if(DSgrn!=NULL){
                mtmp = coef*tmp_DS[ir][ii]; // m=1 90-倾滑
                DSgrn[ir][ii] = CREAL(mtmp);
            }
            if(SSgrn!=NULL){
                mtmp = coef*tmp_SS[ir][ii]; // m=2 走滑 
                SSgrn[ir][ii] = CREAL(mtmp);
            }

        }
    }

    free(tmp_EXP);
    free(tmp_VF);
    free(tmp_HF);
    free(tmp_DD);
    free(tmp_DS);
    free(tmp_SS);
}



void integ_static_grn(
    PYMODEL1D *pymod1d, MYINT nr, MYREAL *rs, MYREAL vmin_ref, MYREAL keps, MYREAL k0, MYREAL Length,
    MYREAL filonLength, MYREAL filonCut, 

    // 返回值，维度2代表Z、R分量，维度3代表Z、R、T分量
    MYREAL EXPgrn[nr][2], // EXZ, EXR 的实部和虚部
    MYREAL VFgrn[nr][2],  // VFZ, VFR 的实部和虚部
    MYREAL HFgrn[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    MYREAL DDgrn[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    MYREAL DSgrn[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    MYREAL SSgrn[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]

    bool calc_upar,
    MYREAL EXPgrn_uiz[nr][2], // EXZ, EXR 的实部和虚部
    MYREAL VFgrn_uiz[nr][2],  // VFZ, VFR 的实部和虚部
    MYREAL HFgrn_uiz[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    MYREAL DDgrn_uiz[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    MYREAL DSgrn_uiz[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    MYREAL SSgrn_uiz[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]
    MYREAL EXPgrn_uir[nr][2], // EXZ, EXR 的实部和虚部
    MYREAL VFgrn_uir[nr][2],  // VFZ, VFR 的实部和虚部
    MYREAL HFgrn_uir[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    MYREAL DDgrn_uir[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    MYREAL DSgrn_uir[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    MYREAL SSgrn_uir[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]

    const char *statsstr // 积分结果输出
){
    MYREAL rmax=rs[findMinMax_MYREAL(rs, nr, true)];   // 最大震中距

    // pymod1d -> mod1d
    MODEL1D *mod1d = init_mod1d(pymod1d->n);
    get_mod1d(pymod1d, mod1d);

    const MYREAL hs = (FABS(pymod1d->depsrc - pymod1d->deprcv) < MIN_DEPTH_GAP_SRC_RCV)? 
                      MIN_DEPTH_GAP_SRC_RCV : FABS(pymod1d->depsrc - pymod1d->deprcv); // hs=max(震源和台站深度差,1.0)
    // 乘相应系数
    k0 *= PI/hs;

    if(vmin_ref < RZERO)  keps = -RONE;  // 若使用峰谷平均法，则不使用keps进行收敛判断

    MYREAL k=0.0;
    const MYREAL dk=FABS(PI2/(Length*rmax));     // 波数积分间隔
    const MYREAL filondk = (filonLength > RZERO) ? PI2/(filonLength*rmax) : RZERO;  // Filon积分间隔
    const MYREAL filonK = filonCut/rmax;  // 波数积分和Filon积分的分割点

    const MYREAL kmax = k0;
    // 局部变量，用于求和 sum F(ki,w)Jm(ki*r)ki 
    // 维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型，详见int_Pk()函数内的注释
    MYCOMPLEX (*sum_EXP_J)[3][4] = (EXPgrn != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_EXP_J)) : NULL;
    MYCOMPLEX (*sum_VF_J)[3][4] = (VFgrn != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_VF_J)) : NULL;
    MYCOMPLEX (*sum_HF_J)[3][4] = (HFgrn != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_HF_J)) : NULL;
    MYCOMPLEX (*sum_DC_J)[3][4] = (DDgrn != NULL || DSgrn != NULL || SSgrn != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_DC_J)) : NULL;
    
    MYCOMPLEX (*sum_EXP_uiz_J)[3][4] = (EXPgrn_uiz != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_EXP_uiz_J)) : NULL;
    MYCOMPLEX (*sum_VF_uiz_J)[3][4] = (VFgrn_uiz != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_VF_uiz_J)) : NULL;
    MYCOMPLEX (*sum_HF_uiz_J)[3][4] = (HFgrn_uiz != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_HF_uiz_J)) : NULL;
    MYCOMPLEX (*sum_DC_uiz_J)[3][4] = (DDgrn_uiz != NULL || DSgrn_uiz != NULL || SSgrn_uiz != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_DC_uiz_J)) : NULL;

    MYCOMPLEX (*sum_EXP_uir_J)[3][4] = (EXPgrn_uir != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_EXP_uir_J)) : NULL;
    MYCOMPLEX (*sum_VF_uir_J)[3][4] = (VFgrn_uir != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_VF_uir_J)) : NULL;
    MYCOMPLEX (*sum_HF_uir_J)[3][4] = (HFgrn_uir != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_HF_uir_J)) : NULL;
    MYCOMPLEX (*sum_DC_uir_J)[3][4] = (DDgrn_uir != NULL || DSgrn_uir != NULL || SSgrn_uir != NULL) ? (MYCOMPLEX(*)[3][4])calloc(nr, sizeof(*sum_DC_uir_J)) : NULL;

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
            for(MYINT m=0; m<3; ++m){
                for(MYINT v=0; v<4; ++v){
                    if(sum_EXP_J) sum_EXP_J[ir][m][v] = RZERO;
                    if(sum_VF_J) sum_VF_J[ir][m][v] = RZERO;
                    if(sum_HF_J) sum_HF_J[ir][m][v] = RZERO;
                    if(sum_DC_J) sum_DC_J[ir][m][v] = RZERO;
    
                    if(sum_EXP_uiz_J) sum_EXP_uiz_J[ir][m][v] = RZERO;
                    if(sum_VF_uiz_J) sum_VF_uiz_J[ir][m][v] = RZERO;
                    if(sum_HF_uiz_J) sum_HF_uiz_J[ir][m][v] = RZERO;
                    if(sum_DC_uiz_J) sum_DC_uiz_J[ir][m][v] = RZERO;
    
                    if(sum_EXP_uir_J) sum_EXP_uir_J[ir][m][v] = RZERO;
                    if(sum_VF_uir_J) sum_VF_uir_J[ir][m][v] = RZERO;
                    if(sum_HF_uir_J) sum_HF_uir_J[ir][m][v] = RZERO;
                    if(sum_DC_uir_J) sum_DC_uir_J[ir][m][v] = RZERO;
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



    // 常规的波数积分
    k = discrete_integ(
        mod1d, dk, (filondk > RZERO)? filonK : kmax, keps, 0.0, nr, rs, 
        sum_EXP_J, sum_VF_J, sum_HF_J, sum_DC_J, 
        calc_upar,
        sum_EXP_uiz_J, sum_VF_uiz_J, sum_HF_uiz_J, sum_DC_uiz_J,
        sum_EXP_uir_J, sum_VF_uir_J, sum_HF_uir_J, sum_DC_uir_J,
        fstats, static_kernel);
    
    // 基于线性插值的Filon积分
    if(filondk > RZERO){
        k = linear_filon_integ(
            mod1d, k, dk, filondk, kmax, keps, 0.0, nr, rs, 
            sum_EXP_J, sum_VF_J, sum_HF_J, sum_DC_J, 
            calc_upar,
            sum_EXP_uiz_J, sum_VF_uiz_J, sum_HF_uiz_J, sum_DC_uiz_J,
            sum_EXP_uir_J, sum_VF_uir_J, sum_HF_uir_J, sum_DC_uir_J,
            fstats, static_kernel);
    }

    // k之后的部分使用峰谷平均法进行显式收敛，建议在浅源地震的时候使用   
    if(vmin_ref < RZERO){
        PTA_method(
            mod1d, k, dk, 0.0, nr, rs, 
            sum_EXP_J, sum_VF_J, sum_HF_J, sum_DC_J, 
            calc_upar,
            sum_EXP_uiz_J, sum_VF_uiz_J, sum_HF_uiz_J, sum_DC_uiz_J,
            sum_EXP_uir_J, sum_VF_uir_J, sum_HF_uir_J, sum_DC_uir_J,
            ptam_fstatsnr, static_kernel);
    }


    
    MYCOMPLEX src_mu = (mod1d->lays + mod1d->isrc)->mu;
    MYCOMPLEX fac = dk * RONE/(RFOUR*PI * src_mu);
    
    // 将积分结果记录到浮点数数组中
    recordin_GRN(
        nr, fac, 
        sum_EXP_J, sum_VF_J, sum_HF_J, sum_DC_J,
        EXPgrn, VFgrn, HFgrn, DDgrn, DSgrn, SSgrn);
    if(calc_upar){
        recordin_GRN(
            nr, fac, 
            sum_EXP_uiz_J, sum_VF_uiz_J, sum_HF_uiz_J, sum_DC_uiz_J,
            EXPgrn_uiz, VFgrn_uiz, HFgrn_uiz, DDgrn_uiz, DSgrn_uiz, SSgrn_uiz);
        recordin_GRN(
            nr, fac, 
            sum_EXP_uir_J, sum_VF_uir_J, sum_HF_uir_J, sum_DC_uir_J,
            EXPgrn_uir, VFgrn_uir, HFgrn_uir, DDgrn_uir, DSgrn_uir, SSgrn_uir);
    }


    // Free allocated memory for temporary variables
    if (sum_EXP_J) free(sum_EXP_J);
    if (sum_VF_J) free(sum_VF_J);
    if (sum_HF_J) free(sum_HF_J);
    if (sum_DC_J) free(sum_DC_J);

    if (sum_EXP_uiz_J) free(sum_EXP_uiz_J);
    if (sum_VF_uiz_J) free(sum_VF_uiz_J);
    if (sum_HF_uiz_J) free(sum_HF_uiz_J);
    if (sum_DC_uiz_J) free(sum_DC_uiz_J);

    if (sum_EXP_uir_J) free(sum_EXP_uir_J);
    if (sum_VF_uir_J) free(sum_VF_uir_J);
    if (sum_HF_uir_J) free(sum_HF_uir_J);
    if (sum_DC_uir_J) free(sum_DC_uir_J);

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