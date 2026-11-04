/**
 * @file   sacio.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-31
 * 
 */

#define _FILE_OFFSET_BITS 64  ///< 统一使用 64 位文件偏移，需在系统头文件之前定义

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "grt/common/sacio.h"
#include "grt/common/const.h"
#include "grt/common/checkerror.h"

#include "sacio.c_"

/** 按磁盘格式和大文件偏移定位 SAC 样本 */
static void seek_SAC_sample(FILE *fp, const char *path, long long first)
{
    off_t offset = GRT_SAC_HEADER_BYTES + (off_t)SAC_DATA_SIZEOF * first;
    if(fseeko(fp, offset, SEEK_SET) != 0) {
        GRTRaiseError("Cannot seek to sample %lld in %s.", first, path);
    }
}

void grt_create_zero_SAC(const char *path, const SACHEAD *hd)
{
    FILE *fp = GRTCheckOpenFile(path, "wb");
    if(write_head_out(path, *hd, fp) != 0) {
        GRTRaiseError("Cannot write SAC header: %s.", path);
    }

    // 真实写入零数据，内存占用与记录长度无关
    float *buffer = GRT_SAFE_CALLOC(GRT_SAC_BLOCK_SAMPLES, sizeof(*buffer));
    size_t remaining = hd->npts;
    while(remaining > 0) {
        size_t count = GRT_MIN(remaining, GRT_SAC_BLOCK_SAMPLES);
        if(fwrite(buffer, sizeof(*buffer), count, fp) != count) {
            GRTRaiseError("Cannot initialize SAC data: %s.", path);
        }
        remaining -= count;
    }
    GRT_SAFE_FREE_PTR(buffer);
    if(fflush(fp) != 0) {
        GRTRaiseError("Cannot flush SAC data: %s.", path);
    }
    off_t actual = ftello(fp);
    off_t expected = GRT_SAC_HEADER_BYTES + (off_t)SAC_DATA_SIZEOF * hd->npts;
    if(actual != expected) {
        GRTRaiseError("Unexpected SAC file size: %s.", path);
    }
    if(fclose(fp) != 0) {
        GRTRaiseError("Cannot close SAC file: %s.", path);
    }
}

void grt_add_SAC_window(const char *path, long long first, size_t count, const float *data)
{
    FILE *fp = GRTCheckOpenFile(path, "r+b");
    SACHEAD hd;
    int swap = read_head_in(path, &hd, fp);
    if(swap < 0 || hd.iftype != ITIME || !hd.leven || first < 0 || first + (long long)count > hd.npts) {
        GRTRaiseError("Invalid SAC update window: %s.", path);
    }

    float *buffer = GRT_SAFE_MALLOC(GRT_SAC_BLOCK_SAMPLES * sizeof(*buffer));
    size_t done = 0;
    while(done < count) {
        size_t n = GRT_MIN(count - done, GRT_SAC_BLOCK_SAMPLES);
        seek_SAC_sample(fp, path, first + (long long)done);
        if(fread(buffer, sizeof(*buffer), n, fp) != n) {
            GRTRaiseError("Cannot read SAC window: %s.", path);
        }
        if(swap) {
            byte_swap((char *)buffer, n * sizeof(*buffer));
        }
        for(size_t i = 0; i < n; ++i) {
            buffer[i] += data[done + i];
        }
        if(swap) {
            byte_swap((char *)buffer, n * sizeof(*buffer));
        }
        // 更新模式的读写切换必须重新定位，写回原来的字节位置
        seek_SAC_sample(fp, path, first + (long long)done);
        if(fwrite(buffer, sizeof(*buffer), n, fp) != n) {
            GRTRaiseError("Cannot write SAC window: %s.", path);
        }
        done += n;
    }
    GRT_SAFE_FREE_PTR(buffer);
    if(fclose(fp) != 0) {
        GRTRaiseError("Cannot close SAC file: %s.", path);
    }
}

void grt_sachead_set_begin(SACHEAD *hd, real_t begin)
{
    hd->b = begin;
    hd->e = hd->b + (hd->npts - 1) * hd->delta;
}

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

void grt_write_SACTRACE(const char *path, SACTRACE *sac)
{
    // 仅写出时量化到 SAC 的存储精度，保留调用方的 real_t 波形
    size_t count = (size_t)sac->hd.npts * (sac->hd.iftype == IXY ? 2 : 1);
    float *data = GRT_SAFE_MALLOC(sizeof(*data)*count);
    for(size_t n = 0; n < count; ++n) {
        data[n] = (float)sac->data[n];
    }
    if(write_sac(path, sac->hd, data) != 0) {
        GRTRaiseError("Failed to write SAC file %s.", path);
    }
    GRT_SAFE_FREE_PTR(data);
}


void grt_free_SACTRACE(SACTRACE *sac)
{
    if (sac == NULL)  return;
    GRT_SAFE_FREE_PTR(sac->data);
    GRT_SAFE_FREE_PTR(sac);
}
