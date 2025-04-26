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
    FILE *f0, MYREAL k, const MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM])
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
    FILE *f0, MYREAL k, 
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