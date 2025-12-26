/**
 * @file   syn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-12
 * 
 *    将合成地震图部分单独用一个源文件和若干函数来管理
 * 
 */

#include <stdio.h>
#include <string.h>
#include "grt/dynamic/syn.h"
#include "grt/common/radiation.h"

void grt_syn(const realChnlGrid srcRadi, const int computeType, const char *dirpath, const char *prefix, SACTRACE *synsac[GRT_CHANNEL_NUM])
{
    GRT_LOOP_ChnlGrid(im, c) {
        int modr = GRT_SRC_M_ORDERS[im];
        if(modr == 0 && GRT_QWV_CODES[c] == 'v')   continue;

        if (computeType == GRT_SYN_COMPUTE_EX) {
            if (im != GRT_SRC_M_EX_INDEX) continue;
        } else if (computeType == GRT_SYN_COMPUTE_SF) {
            if (im != GRT_SRC_M_VF_INDEX && im != GRT_SRC_M_HF_INDEX) continue;
        } else if (computeType == GRT_SYN_COMPUTE_DC) {
            if (im < GRT_SRC_M_DD_INDEX) continue;
        } else if (computeType == GRT_SYN_COMPUTE_MT) {
            if (im < GRT_SRC_M_DD_INDEX && im != GRT_SRC_M_EX_INDEX) continue;
        } else {
            GRTRaiseError("Not Supported.");
        }

        const real_t coef = srcRadi[im][c];
        const char ch = GRT_ZRT_CODES[c];

        // 读取格林函数
        SACTRACE *grnsac = NULL;
        {
            char *grnpath = NULL;
            GRT_SAFE_ASPRINTF(&grnpath, "%s/%s%s%c.sac", dirpath, prefix, GRT_SRC_M_NAME_ABBR[im], ch);
            grnsac = grt_read_SACTRACE(grnpath, false);
            GRT_SAFE_FREE_PTR(grnpath);
        }

        // 如果是第一次读取，先初始化结果指针
        if (synsac[0] == NULL) {
            for(int ii = 0; ii < GRT_CHANNEL_NUM; ++ii){
                synsac[ii] = grt_copy_SACTRACE(grnsac, true);
            }
        }

        // 线性叠加
        for(int n = 0; n < grnsac->hd.npts; ++n){
            synsac[c]->data[n] += grnsac->data[n] * coef;
        }

        grt_free_SACTRACE(grnsac);
    }
}