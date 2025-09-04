/**
 * @file   sacio2.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-31
 * 
 */



#include <stdio.h>
#include <stdlib.h>

#include "grt/common/sacio2.h"
#include "grt/common/sacio.h"

#include "grt/common/checkerror.h"


void grt_read_SAC_HEAD(const char *command, const char *name, SACHEAD *hd){
    int lswap = read_sac_head(name, hd);
    if(lswap == -1){
        GRTRaiseError("[%s] read %s head failed.\n", command, name);
    }
}


float * grt_read_SAC(const char *command, const char *name, SACHEAD *hd, float *arrout){
    float *arrin=NULL;
    if((arrin = read_sac(name, hd)) == NULL){
        GRTRaiseError("[%s] read %s failed.\n", command, name);
    }

    if(arrout!=NULL){
        for(int i=0; i<hd->npts; ++i)  arrout[i] = arrin[i];
        free(arrin);
        arrin = arrout;
    }
    
    return arrin;
}