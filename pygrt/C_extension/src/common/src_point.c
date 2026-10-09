/**
 * @file   src_point.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-10-03
 */

#include "grt/common/src_point.h"

SRC_POINT *grt_src_points_from_faults(size_t nfault, const FINITE_FAULT *faults, size_t *npts)
{
    // 先统计全部剖分点，再按断层及断层内顺序一次性展开
    *npts = 0;
    for(size_t i = 0; i < nfault; ++i) {
        *npts += faults[i].nW * faults[i].nL;
    }
    SRC_POINT *points = GRT_SAFE_MALLOC(*npts * sizeof(*points));
    size_t ipt = 0;
    for(size_t i = 0; i < nfault; ++i) {
        const FINITE_FAULT *fault = &faults[i];
        for(size_t isub = 0; isub < fault->nW * fault->nL; ++isub) {
            points[ipt++] = (SRC_POINT){.north = fault->north[isub], .east = fault->east[isub], .depth = fault->depth[isub],
                                       .fault = fault, .isub = isub};
        }
    }
    return points;
}
