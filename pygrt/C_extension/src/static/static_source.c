/**
 * @file   static_source.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 * 以下代码实现的是 静态震源系数————剪切源， 参考：
 *             1. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *                透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 *
 */


#include <stdio.h>
#include <complex.h>
#include <string.h>

#include "grt/static/static_source.h"


inline GCC_ALWAYS_INLINE void _source_PSV(const cplx_t delta, const real_t k, cplx_t coef[GRT_SRC_M_NUM][GRT_QWV_NUM][2])
{
    cplx_t tmp;
    cplx_t A = 1.0+delta;

    // 爆炸源
    coef[0][0][0] = tmp = (delta-1.0)/A;         coef[0][0][1] = tmp;    

    // 垂直力源
    coef[1][0][0] = tmp = -1.0/(2.0*A*k);        coef[1][0][1] = - tmp;   
    coef[1][1][0] = tmp;                         coef[1][1][1] = - tmp;

    // 水平力源
    coef[2][0][0] = tmp = 1.0/(2.0*A*k);          coef[2][0][1] = tmp;   
    coef[2][1][0] = - tmp;                        coef[2][1][1] = - tmp;

    // 剪切位错
    // m=0
    coef[3][0][0] = tmp = (-1.0+4.0*delta)/(2.0*A);    coef[3][0][1] = tmp;
    coef[3][1][0] = tmp = -3.0/(2.0*A);                coef[3][1][1] = tmp;
    // m=1
    coef[4][0][0] = tmp = -delta/A;                       coef[4][0][1] = -tmp;
    coef[4][1][0] = tmp = 1.0/A;                          coef[4][1][1] = -tmp;
    // m=2
    coef[5][0][0] = tmp = 1.0/(2.0*A);                   coef[5][0][1] = tmp;
    coef[5][1][0] = tmp = -1.0/(2.0*A);                  coef[5][1][1] = tmp;
}


inline GCC_ALWAYS_INLINE void _source_SH(const real_t k, cplx_t coef[GRT_SRC_M_NUM][GRT_QWV_NUM][2])
{
    cplx_t tmp;

    // 水平力源
    coef[2][2][0] = tmp = - 1.0/k;                 coef[2][2][1] = tmp;

    // 剪切位错
    // m=1
    coef[4][2][0] = tmp = 1.0;                     coef[4][2][1] = -tmp;
    // m=2
    coef[5][2][0] = tmp = - 1.0;                   coef[5][2][1] = tmp;
}


void grt_static_source_coef(GRT_MODEL1D *mod1d)
{
    // 先全部赋0 
    memset(mod1d->src_coef, 0, sizeof(cplx_t)*GRT_SRC_M_NUM*GRT_QWV_NUM*2);
    
    grt_static_source_coef_PSV(mod1d);
    grt_static_source_coef_SH(mod1d);
}


void grt_static_source_coef_PSV(GRT_MODEL1D *mod1d)
{
    size_t isrc = mod1d->isrc;
    cplx_t delta = mod1d->delta[isrc];
    real_t k = mod1d->k;

    _source_PSV(delta, k, mod1d->src_coef);
}


void grt_static_source_coef_SH(GRT_MODEL1D *mod1d)
{
    real_t k = mod1d->k;
    
    _source_SH(k, mod1d->src_coef);
}


