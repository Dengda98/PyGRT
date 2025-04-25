/**
 * @file   iostats.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 将波数积分过程中的核函数F(k,w)以及F(k,w)Jm(kr)k的值记录在文件中
 * 
 */

#include <stdio.h> 
#include <stdbool.h>
#include <complex.h>

#include "common/iostats.h"
#include "common/const.h"



void write_stats(
    FILE *f0, MYREAL k, const MYCOMPLEX QWV[GRT_SRC_M_COUNTS][GRT_SRC_QWV_COUNTS])
{
    fwrite(&k, sizeof(MYREAL), 1, f0);

    fwrite(&QWV[0][0], sizeof(MYCOMPLEX), 2, f0);

    fwrite(&QWV[1][0], sizeof(MYCOMPLEX), 2, f0);
    
    fwrite(&QWV[2][0], sizeof(MYCOMPLEX), 3, f0);

    fwrite(&QWV[3][0], sizeof(MYCOMPLEX), 2, f0);

    fwrite(&QWV[4][0], sizeof(MYCOMPLEX), 3, f0);

    fwrite(&QWV[5][0], sizeof(MYCOMPLEX), 3, f0);

}




void write_stats_ptam(
    FILE *f0, MYREAL k, MYINT maxNpt, 
    MYREAL Kpt[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt],
    MYCOMPLEX Fpt[GRT_SRC_M_COUNTS][GRT_SRC_P_COUNTS][maxNpt]
){

    for(MYINT i=0; i<maxNpt; ++i){
        for(MYINT im=0; im<GRT_SRC_M_COUNTS; ++im){
            MYINT modr = GRT_SRC_M_ORDERS[im];
            for(MYINT v=0; v<GRT_SRC_P_COUNTS; ++v){
                if(modr == 0 && v!=0 && v!=2)  continue;
                
                fwrite(&Kpt[im][v][i], sizeof(MYREAL),  1, f0);
                fwrite(&Fpt[im][v][i], sizeof(MYCOMPLEX), 1, f0);
            }
        }
    }
    
}