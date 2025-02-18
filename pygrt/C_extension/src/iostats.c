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

#include "iostats.h"
#include "const.h"



void write_stats(
    FILE *f0, MYREAL k, 
    const MYCOMPLEX EXP_qwv[3][3], const MYCOMPLEX VF_qwv[3][3], 
    const MYCOMPLEX HF_qwv[3][3],  const MYCOMPLEX DC_qwv[3][3], 
    const MYCOMPLEX EXP_J[3][4], const MYCOMPLEX VF_J[3][4], 
    const MYCOMPLEX HF_J[3][4],  const MYCOMPLEX DC_J[3][4])
{
    fwrite(&k, sizeof(MYREAL), 1, f0);

    fwrite(&EXP_qwv[0][0], sizeof(MYCOMPLEX), 2, f0);

    fwrite(&VF_qwv[0][0], sizeof(MYCOMPLEX), 2, f0);
    
    fwrite(&HF_qwv[1][0], sizeof(MYCOMPLEX), 3, f0);

    fwrite(&DC_qwv[0][0], sizeof(MYCOMPLEX), 2, f0);

    fwrite(&DC_qwv[1][0], sizeof(MYCOMPLEX), 3, f0);

    fwrite(&DC_qwv[2][0], sizeof(MYCOMPLEX), 3, f0);

    fwrite(&EXP_J[0][0], sizeof(MYCOMPLEX), 1, f0);
    fwrite(&EXP_J[0][2], sizeof(MYCOMPLEX), 1, f0);

    fwrite(&VF_J[0][0], sizeof(MYCOMPLEX), 1, f0);
    fwrite(&VF_J[0][2], sizeof(MYCOMPLEX), 1, f0);

    fwrite(&HF_J[1], sizeof(MYCOMPLEX), 4, f0);

    fwrite(&DC_J[0][0], sizeof(MYCOMPLEX), 1, f0);
    fwrite(&DC_J[0][2], sizeof(MYCOMPLEX), 1, f0);

    fwrite(&DC_J[1], sizeof(MYCOMPLEX), 4, f0);
    fwrite(&DC_J[2], sizeof(MYCOMPLEX), 4, f0);

}




void write_stats_ptam(
    FILE *f0, MYREAL k, MYINT maxNpt, 
    const MYCOMPLEX EXPpt[3][4][maxNpt], const MYCOMPLEX VFpt[3][4][maxNpt],
    const MYCOMPLEX HFpt[3][4][maxNpt],  const MYCOMPLEX DCpt[3][4][maxNpt],
    const MYREAL kEXPpt[3][4][maxNpt], const MYREAL kVFpt[3][4][maxNpt],
    const MYREAL kHFpt[3][4][maxNpt],  const MYREAL kDCpt[3][4][maxNpt])
{
    
    MYINT i, m, v;

    for(i=0; i<maxNpt; ++i){
        for(m=0; m<3; ++m){
            for(v=0; v<4; ++v){
                if(m==0 && (v==0||v==2)){
                    fwrite(&kEXPpt[m][v][i], sizeof(MYREAL),  1, f0);
                    fwrite( &EXPpt[m][v][i], sizeof(MYCOMPLEX), 1, f0);

                    fwrite(&kVFpt[m][v][i], sizeof(MYREAL),  1, f0);
                    fwrite( &VFpt[m][v][i], sizeof(MYCOMPLEX), 1, f0);
                }

                if(m==1){
                    fwrite(&kHFpt[m][v][i], sizeof(MYREAL),  1, f0);
                    fwrite( &HFpt[m][v][i], sizeof(MYCOMPLEX), 1, f0);
                }

                if(((m==0 && (v==0||v==2)) || m!=0)){
                    fwrite(&kDCpt[m][v][i], sizeof(MYREAL),  1, f0);
                    fwrite( &DCpt[m][v][i], sizeof(MYCOMPLEX), 1, f0);
                }

            }
        }
    }
    
}