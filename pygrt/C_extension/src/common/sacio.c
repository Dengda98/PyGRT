/**
 * @file   sacio.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-31
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grt/common/sacio.h"
#include "grt/common/const.h"
#include "grt/common/checkerror.h"

#include "sacio.c_"

SACTRACE * grt_read_SACTRACE(const char *path, const bool headonly)
{
    GRTCheckFileExist(path);

    SACTRACE *sac = GRT_SAFE_CALLOC(1, sizeof(SACTRACE));

    if (headonly) {
        if(read_sac_head(path, &sac->hd) != 0){
            GRTRaiseError("read %s head failed.\n", path);
        }
        return sac;
    }
    
    // 文件数据为 float，读入后立即转换为工作精度
    float *data = read_sac(path, &sac->hd);
    if (data == NULL){
        GRTRaiseError("read %s failed.\n", path);
    }
    size_t count = (size_t)sac->hd.npts * (sac->hd.iftype == IXY ? 2 : 1);
    sac->data = GRT_SAFE_MALLOC(sizeof(*sac->data)*count);
    for(size_t n = 0; n < count; ++n) {
        sac->data[n] = data[n];
    }
    GRT_SAFE_FREE_PTR(data);
    
    return sac;
}

SACTRACE * grt_copy_SACTRACE(SACTRACE *sac, bool zero_value)
{
    SACTRACE *sac2 = GRT_SAFE_CALLOC(1, sizeof(SACTRACE));
    *sac2 = *sac;
    size_t count = (size_t)sac->hd.npts * (sac->hd.iftype == IXY ? 2 : 1);
    sac2->data = GRT_SAFE_CALLOC(count, sizeof(*sac2->data));
    if(!zero_value) {
        memcpy(sac2->data, sac->data, sizeof(*sac2->data)*count);
    }
    return sac2;
}

SACTRACE * grt_new_SACTRACE(real_t dt, int nt, real_t b0)
{
    SACTRACE *sac = GRT_SAFE_CALLOC(1, sizeof(SACTRACE));
    sac->hd = new_sac_head((float)dt, nt, (float)b0);
    sac->data = GRT_SAFE_CALLOC(nt, sizeof(*sac->data));
    return sac;
}

int grt_write_SACTRACE(const char *path, SACTRACE *sac)
{
    // 仅写出时量化到 SAC 的存储精度，保留调用方的 real_t 波形
    size_t count = (size_t)sac->hd.npts * (sac->hd.iftype == IXY ? 2 : 1);
    float *data = GRT_SAFE_MALLOC(sizeof(*data)*count);
    for(size_t n = 0; n < count; ++n) {
        data[n] = (float)sac->data[n];
    }
    int status = write_sac(path, sac->hd, data);
    GRT_SAFE_FREE_PTR(data);
    return status;
}


void grt_free_SACTRACE(SACTRACE *sac)
{
    if (sac == NULL)  return;
    GRT_SAFE_FREE_PTR(sac->data);
    GRT_SAFE_FREE_PTR(sac);
}
