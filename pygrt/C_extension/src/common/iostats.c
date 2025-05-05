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

#include "common/iostats.h"
#include "common/const.h"



void write_stats(
    FILE *f0, MYREAL k, const MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM])
{
    fwrite(&k, sizeof(MYREAL), 1, f0);

    for(MYINT im=0; im<SRC_M_NUM; ++im){
        MYINT modr = SRC_M_ORDERS[im];
        for(MYINT c=0; c<QWV_NUM; ++c){
            if(modr == 0 && qwvchs[c] == 'v')   continue;

            fwrite(&QWV[im][c], sizeof(MYCOMPLEX), 1, f0);
        }
    }
}


MYINT extract_stats(FILE *bf0, FILE *af0){
    // 打印标题
    if(bf0 == NULL){
        char K[20];
        snprintf(K, sizeof(K), GRT_STRING_FMT, "k");  K[0]='#';
        fprintf(af0, "%s", K);

        for(MYINT im=0; im<SRC_M_NUM; ++im){
            MYINT modr = SRC_M_ORDERS[im];
            for(MYINT c=0; c<QWV_NUM; ++c){
                if(modr == 0 && qwvchs[c] == 'v')   continue;

                snprintf(K, sizeof(K), "%s_%c", SRC_M_NAME_ABBR[im], qwvchs[c]);
                fprintf(af0, GRT_STR_CMPLX_FMT, K);
            }
        }

        return 0;
    }

    MYREAL k;
    MYCOMPLEX val;

    if(1 != fread(&k, sizeof(MYREAL), 1, bf0))  return -1;
    fprintf(af0, GRT_REAL_FMT, k);

    for(MYINT im=0; im<SRC_M_NUM; ++im){
        MYINT modr = SRC_M_ORDERS[im];
        for(MYINT c=0; c<QWV_NUM; ++c){
            if(modr == 0 && qwvchs[c] == 'v')   continue;

            if(1 != fread(&val, sizeof(MYCOMPLEX), 1, bf0))  return -1;
            fprintf(af0, GRT_CMPLX_FMT, creal(val), cimag(val));
        }
    }

    return 0;
}


void write_stats_ptam(
    FILE *f0, 
    MYREAL Kpt[SRC_M_NUM][INTEG_NUM][PTAM_MAX_PT],
    MYCOMPLEX Fpt[SRC_M_NUM][INTEG_NUM][PTAM_MAX_PT]
){

    for(MYINT i=0; i<PTAM_MAX_PT; ++i){
        for(MYINT im=0; im<SRC_M_NUM; ++im){
            MYINT modr = SRC_M_ORDERS[im];
            for(MYINT v=0; v<INTEG_NUM; ++v){
                if(modr == 0 && v!=0 && v!=2)  continue;
                
                fwrite(&Kpt[im][v][i], sizeof(MYREAL),  1, f0);
                fwrite(&Fpt[im][v][i], sizeof(MYCOMPLEX), 1, f0);
            }
        }
    }
    
}


MYINT extract_stats_ptam(FILE *bf0, FILE *af0){
    // 打印标题
    if(bf0 == NULL){
        char K[20], K2[20];
        MYINT icol=0;

        for(MYINT im=0; im<SRC_M_NUM; ++im){
            MYINT modr = SRC_M_ORDERS[im];
            for(MYINT v=0; v<INTEG_NUM; ++v){
                if(modr == 0 && v!=0 && v!=2)  continue;

                snprintf(K2, sizeof(K2), "sum_%s_%d_k", SRC_M_NAME_ABBR[im], v);
                if(icol==0){
                    snprintf(K, sizeof(K), GRT_STRING_FMT, K2);  K2[0]='#';
                    fprintf(af0, "%s", K);
                } else {
                    fprintf(af0, GRT_STRING_FMT, K2);
                }
                snprintf(K2, sizeof(K2), "sum_%s_%d", SRC_M_NAME_ABBR[im], v);
                fprintf(af0, GRT_STR_CMPLX_FMT, K2);
                
                icol++;
            }
        }

        return 0;
    }


    MYREAL k;
    MYCOMPLEX val;

    for(MYINT im=0; im<SRC_M_NUM; ++im){
        MYINT modr = SRC_M_ORDERS[im];
        for(MYINT v=0; v<INTEG_NUM; ++v){
            if(modr == 0 && v!=0 && v!=2)  continue;

            if(1 != fread(&k, sizeof(MYREAL), 1, bf0))  return -1;
            fprintf(af0, GRT_REAL_FMT, k);
            if(1 != fread(&val, sizeof(MYCOMPLEX), 1, bf0))  return -1;
            fprintf(af0, GRT_CMPLX_FMT, creal(val), cimag(val));
        }
    }

    return 0;
}