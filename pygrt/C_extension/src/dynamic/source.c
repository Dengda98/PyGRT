/**
 * @file   source.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现的是 震源系数————爆炸源，垂直力源，水平力源，剪切源， 参考：
 *             1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *
 */


#include <stdio.h>
#include <complex.h>
#include <string.h>

#include "grt/dynamic/source.h"


inline GCC_ALWAYS_INLINE void _source_PSV(
    const cplx_t xa, const cplx_t caca, 
    const cplx_t xb, const cplx_t cbcb, const real_t k, cplxChnlGrid coefD, cplxChnlGrid coefU)
{
    cplx_t tmp;

    // 爆炸源， 通过(4.9.8)的矩张量源公式，提取各向同性的量(M11+M22+M33)，-a+k^2/a -> ka^2/a
    coefD[0][0] = tmp = (caca / xa) * k;   coefU[0][0] = tmp;    
    
    // 垂直力源 (4.6.15)
    coefD[1][0] = tmp = -1.0;              coefU[1][0] = - tmp;
    coefD[1][1] = tmp = -1.0 / xb;         coefU[1][1] = tmp;

    // 水平力源 (4.6.21,26)
    coefD[2][0] = tmp = -1.0 / xa;       coefU[2][0] = tmp;
    coefD[2][1] = tmp = -1.0;            coefU[2][1] = - tmp;

    // 剪切位错 (4.8.34)
    // m=0
    coefD[3][0] = tmp = ((2.0*caca - 3.0) / xa) * k;    coefU[3][0] = tmp;
    coefD[3][1] = tmp = -3.0*k;                         coefU[3][1] = - tmp;
    // m=1
    coefD[4][0] = tmp = 2.0*k;                      coefU[4][0] = - tmp;
    coefD[4][1] = tmp = ((2.0 - cbcb) / xb) * k;    coefU[4][1] = tmp;

    // m=2
    coefD[5][0] = tmp = - (1.0 / xa) * k;            coefU[5][0] = tmp;
    coefD[5][1] = tmp = - k;                         coefU[5][1] = - tmp;

}

inline GCC_ALWAYS_INLINE void _source_SH(const cplx_t xb, const cplx_t cbcb, const real_t k, cplxChnlGrid coefD, cplxChnlGrid coefU)
{
    cplx_t tmp;

    // 水平力源 (4.6.21,26)
    coefD[2][2] = tmp = cbcb / xb;    coefU[2][2] = tmp;

    // 剪切位错 (4.8.34)
    // m=1
    coefD[4][2] = tmp = - cbcb * k;              coefU[4][2] = - tmp;

    // m=2
    coefD[5][2] = tmp = (cbcb / xb) * k;         coefU[5][2] = tmp;
}


void grt_source_coef(GRT_MODEL1D *mod1d)
{
    // 先全部赋0 
    memset(mod1d->src_coefD, 0, sizeof(cplxChnlGrid));
    memset(mod1d->src_coefU, 0, sizeof(cplxChnlGrid));

    grt_source_coef_PSV(mod1d);
    grt_source_coef_SH(mod1d);
}


void grt_source_coef_PSV(GRT_MODEL1D *mod1d)
{
    size_t isrc = mod1d->isrc;
    cplx_t xa = mod1d->xa[isrc];
    cplx_t caca = mod1d->caca[isrc];
    cplx_t xb = mod1d->xb[isrc];
    cplx_t cbcb = mod1d->cbcb[isrc];
    real_t k = mod1d->k;

    _source_PSV(xa, caca, xb, cbcb, k, mod1d->src_coefD, mod1d->src_coefU);
}


void grt_source_coef_SH(GRT_MODEL1D *mod1d)
{
    size_t isrc = mod1d->isrc;
    cplx_t xb = mod1d->xb[isrc];
    cplx_t cbcb = mod1d->cbcb[isrc];
    real_t k = mod1d->k;

    _source_SH(xb, cbcb, k, mod1d->src_coefD, mod1d->src_coefU);
}


