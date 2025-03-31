/**
 * @file   sacio2.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-31
 * 
 */



#include <stdio.h>
#include <stdlib.h>

#include "common/sacio2.h"
#include "common/sacio.h"
#include "common/colorstr.h"


void read_SAC_HEAD(const char *command, const char *name, SACHEAD *hd){
    int lswap = read_sac_head(name, hd);
    if(lswap == -1){
        fprintf(stderr, "[%s] " BOLD_RED "read %s head failed.\n" DEFAULT_RESTORE, command, name);
        exit(EXIT_FAILURE);
    }
}


float * read_SAC(const char *command, const char *name, SACHEAD *hd, float *arrout){
    float *arrin=NULL;
    if((arrin = read_sac(name, hd)) == NULL){
        fprintf(stderr, "[%s] " BOLD_RED "read %s failed.\n" DEFAULT_RESTORE, command, name);
        exit(EXIT_FAILURE);
    }

    if(arrout!=NULL){
        for(int i=0; i<hd->npts; ++i)  arrout[i] = arrin[i];
        free(arrin);
        arrin = arrout;
    }
    
    return arrin;
}