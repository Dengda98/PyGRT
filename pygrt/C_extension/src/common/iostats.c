/**
 * @file   iostats.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 将波数积分过程中的核函数F(k,w)以及F(k,w)Jm(kr)k的值记录在文件中
 * 
 */

#include <stdio.h> 
#include <string.h>
#include <stdbool.h>
#include <complex.h>

#include "grt/common/iostats.h"
#include "grt/common/const.h"



void grt_write_stats(
    FILE *f0, MYREAL k, const MYCOMPLEX QWV[GRT_SRC_M_NUM][GRT_QWV_NUM])
{
    fwrite(&k, sizeof(MYREAL), 1, f0);

    for(MYINT im=0; im<GRT_SRC_M_NUM; ++im){
        MYINT modr = GRT_SRC_M_ORDERS[im];
        for(MYINT c=0; c<GRT_QWV_NUM; ++c){
            if(modr == 0 && GRT_QWV_CODES[c] == 'v')   continue;

            fwrite(&QWV[im][c], sizeof(MYCOMPLEX), 1, f0);
        }
    }
}


MYINT grt_extract_stats(FILE *bf0, FILE *af0){
    // 打印标题
    if(bf0 == NULL){
        char K[20];
        snprintf(K, sizeof(K), GRT_STRING_FMT, "k");  K[0]=GRT_COMMENT_HEAD;
        fprintf(af0, "%s", K);

        for(MYINT im=0; im<GRT_SRC_M_NUM; ++im){
            MYINT modr = GRT_SRC_M_ORDERS[im];
            for(MYINT c=0; c<GRT_QWV_NUM; ++c){
                if(modr == 0 && GRT_QWV_CODES[c] == 'v')   continue;

                snprintf(K, sizeof(K), "%s_%c", GRT_SRC_M_NAME_ABBR[im], GRT_QWV_CODES[c]);
                fprintf(af0, GRT_STR_CMPLX_FMT, K);
            }
        }

        return 0;
    }

    MYREAL k;
    MYCOMPLEX val;

    if(1 != fread(&k, sizeof(MYREAL), 1, bf0))  return -1;
    fprintf(af0, GRT_REAL_FMT, k);

    for(MYINT im=0; im<GRT_SRC_M_NUM; ++im){
        MYINT modr = GRT_SRC_M_ORDERS[im];
        for(MYINT c=0; c<GRT_QWV_NUM; ++c){
            if(modr == 0 && GRT_QWV_CODES[c] == 'v')   continue;

            if(1 != fread(&val, sizeof(MYCOMPLEX), 1, bf0))  return -1;
            fprintf(af0, GRT_CMPLX_FMT, creal(val), cimag(val));
        }
    }

    return 0;
}


void grt_write_stats_ptam(
    FILE *f0, 
    MYREAL Kpt[GRT_SRC_M_NUM][GRT_INTEG_NUM][GRT_PTAM_PT_MAX],
    MYCOMPLEX Fpt[GRT_SRC_M_NUM][GRT_INTEG_NUM][GRT_PTAM_PT_MAX]
){

    for(MYINT i=0; i<GRT_PTAM_PT_MAX; ++i){
        for(MYINT im=0; im<GRT_SRC_M_NUM; ++im){
            MYINT modr = GRT_SRC_M_ORDERS[im];
            for(MYINT v=0; v<GRT_INTEG_NUM; ++v){
                if(modr == 0 && v!=0 && v!=2)  continue;
                
                fwrite(&Kpt[im][v][i], sizeof(MYREAL),  1, f0);
                fwrite(&Fpt[im][v][i], sizeof(MYCOMPLEX), 1, f0);
            }
        }
    }
    
}


MYINT grt_extract_stats_ptam(FILE *bf0, FILE *af0){
    // 打印标题
    if(bf0 == NULL){
        char K[20], K2[20];
        MYINT icol=0;

        for(MYINT im=0; im<GRT_SRC_M_NUM; ++im){
            MYINT modr = GRT_SRC_M_ORDERS[im];
            for(MYINT v=0; v<GRT_INTEG_NUM; ++v){
                if(modr == 0 && v!=0 && v!=2)  continue;

                snprintf(K2, sizeof(K2), "sum_%s_%d_k", GRT_SRC_M_NAME_ABBR[im], v);
                if(icol==0){
                    snprintf(K, sizeof(K), GRT_STRING_FMT, K2);  K2[0]=GRT_COMMENT_HEAD;
                    fprintf(af0, "%s", K);
                } else {
                    fprintf(af0, GRT_STRING_FMT, K2);
                }
                snprintf(K2, sizeof(K2), "sum_%s_%d", GRT_SRC_M_NAME_ABBR[im], v);
                fprintf(af0, GRT_STR_CMPLX_FMT, K2);
                
                icol++;
            }
        }

        return 0;
    }


    MYREAL k;
    MYCOMPLEX val;

    for(MYINT im=0; im<GRT_SRC_M_NUM; ++im){
        MYINT modr = GRT_SRC_M_ORDERS[im];
        for(MYINT v=0; v<GRT_INTEG_NUM; ++v){
            if(modr == 0 && v!=0 && v!=2)  continue;

            if(1 != fread(&k, sizeof(MYREAL), 1, bf0))  return -1;
            fprintf(af0, GRT_REAL_FMT, k);
            if(1 != fread(&val, sizeof(MYCOMPLEX), 1, bf0))  return -1;
            fprintf(af0, GRT_CMPLX_FMT, creal(val), cimag(val));
        }
    }

    return 0;
}