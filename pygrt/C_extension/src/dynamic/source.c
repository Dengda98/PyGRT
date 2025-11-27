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


inline void _source_PSV(
    const cplx_t xa, const cplx_t caca, 
    const cplx_t xb, const cplx_t cbcb, const real_t k, cplx_t coef[GRT_SRC_M_NUM][GRT_QWV_NUM][2])
{
    cplx_t tmp;

    // 爆炸源， 通过(4.9.8)的矩张量源公式，提取各向同性的量(M11+M22+M33)，-a+k^2/a -> ka^2/a
    coef[0][0][0] = tmp = (caca / xa) * k;   coef[0][0][1] = tmp;    
    
    // 垂直力源 (4.6.15)
    coef[1][0][0] = tmp = -1.0;              coef[1][0][1] = - tmp;
    coef[1][1][0] = tmp = -1.0 / xb;         coef[1][1][1] = tmp;

    // 水平力源 (4.6.21,26)
    coef[2][0][0] = tmp = -1.0 / xa;       coef[2][0][1] = tmp;
    coef[2][1][0] = tmp = -1.0;            coef[2][1][1] = - tmp;

    // 剪切位错 (4.8.34)
    // m=0
    coef[3][0][0] = tmp = ((2.0*caca - 3.0) / xa) * k;    coef[3][0][1] = tmp;
    coef[3][1][0] = tmp = -3.0*k;                         coef[3][1][1] = - tmp;
    // m=1
    coef[4][0][0] = tmp = 2.0*k;                      coef[4][0][1] = - tmp;
    coef[4][1][0] = tmp = ((2.0 - cbcb) / xb) * k;    coef[4][1][1] = tmp;

    // m=2
    coef[5][0][0] = tmp = - (1.0 / xa) * k;            coef[5][0][1] = tmp;
    coef[5][1][0] = tmp = - k;                         coef[5][1][1] = - tmp;

}

inline void _source_SH(const cplx_t xb, const cplx_t cbcb, const real_t k, cplx_t coef[GRT_SRC_M_NUM][GRT_QWV_NUM][2])
{
    cplx_t tmp;

    // 水平力源 (4.6.21,26)
    coef[2][2][0] = tmp = cbcb / xb;    coef[2][2][1] = tmp;

    // 剪切位错 (4.8.34)
    // m=1
    coef[4][2][0] = tmp = - cbcb * k;              coef[4][2][1] = - tmp;

    // m=2
    coef[5][2][0] = tmp = (cbcb / xb) * k;         coef[5][2][1] = tmp;
}


void grt_source_coef(GRT_MODEL1D *mod1d)
{
    // 先全部赋0 
    memset(mod1d->src_coef, 0, sizeof(cplx_t)*GRT_SRC_M_NUM*GRT_QWV_NUM*2);

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

    _source_PSV(xa, caca, xb, cbcb, k, mod1d->src_coef);
}


void grt_source_coef_SH(GRT_MODEL1D *mod1d)
{
    size_t isrc = mod1d->isrc;
    cplx_t xb = mod1d->xb[isrc];
    cplx_t cbcb = mod1d->cbcb[isrc];
    real_t k = mod1d->k;

    _source_SH(xb, cbcb, k, mod1d->src_coef);
}


