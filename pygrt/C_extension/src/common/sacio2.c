/**
 * @file   sacio2.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-31
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grt/common/sacio2.h"
#include "grt/common/const.h"
#include "grt/common/checkerror.h"


SACTRACE * grt_read_SACTRACE(const char *path, const bool headonly)
{
    GRTCheckFileExist(path);

    SACTRACE *sac = (SACTRACE *)calloc(1, sizeof(SACTRACE));

    if (headonly) {
        if(read_sac_head(path, &sac->hd) != 0){
            GRTRaiseError("read %s head failed.\n", path);
        }
        return sac;
    }
    
    if ((sac->data = read_sac(path, &sac->hd)) == NULL){
        GRTRaiseError("read %s failed.\n", path);
    }
    
    return sac;
}

SACTRACE * grt_copy_SACTRACE(SACTRACE *sac, bool zero_value)
{
    SACTRACE *sac2 = (SACTRACE *)calloc(1, sizeof(SACTRACE));
    *sac2 = *sac;
    sac2->data = (float *)calloc(sac->hd.npts, sizeof(float));
    if(!zero_value) memcpy(sac2->data, sac->data, sizeof(float)*sac->hd.npts);
    return sac2;
}

SACTRACE * grt_new_SACTRACE(float dt, int nt, float b0)
{
    SACTRACE *sac = (SACTRACE *)calloc(1, sizeof(SACTRACE));
    sac->hd = new_sac_head(dt, nt, b0);
    sac->data = (float *)calloc(nt, sizeof(float));
    return sac;
}

int grt_write_SACTRACE(const char *path, SACTRACE *sac)
{
    return write_sac(path, sac->hd, sac->data);
}


void grt_free_SACTRACE(SACTRACE *sac)
{
    if (sac == NULL)  return;
    GRT_SAFE_FREE_PTR(sac->data);
    GRT_SAFE_FREE_PTR(sac);
}