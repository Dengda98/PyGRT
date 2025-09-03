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

#include "grt/static/static_source.h"
#include "grt/common/const.h"

void grt_static_source_coef_PSV(MYCOMPLEX delta, MYREAL k, MYCOMPLEX coef[GRT_SRC_M_NUM][GRT_QWV_NUM-1][2])
{
    // 先全部赋0 
    for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
        for(MYINT j=0; j<GRT_QWV_NUM-1; ++j){
            for(MYINT p=0; p<2; ++p){
                coef[i][j][p] = 0.0;
            }
        }
    }

    MYCOMPLEX tmp;
    MYCOMPLEX A = 1.0+delta;

    // 爆炸源
    coef[0][0][0] = tmp = (delta-1.0)/A;         coef[0][0][1] = tmp;    

    // 垂直力源
    coef[1][0][0] = tmp = -1.0/(2.0*A*k);        coef[1][0][1] = - tmp;   
    coef[1][1][0] = tmp;                           coef[1][1][1] = - tmp;

    // 水平力源
    coef[2][0][0] = tmp = 1.0/(2.0*A*k);        coef[2][0][1] = tmp;   
    coef[2][1][0] = - tmp;                        coef[2][1][1] = - tmp;

    // 剪切位错
    // m=0
    coef[3][0][0] = tmp = (-1.0+4.0*delta)/(2.0*A);    coef[3][0][1] = tmp;
    coef[3][1][0] = tmp = -3.0/(2.0*A);                coef[3][1][1] = tmp;
    // m=1
    coef[4][0][0] = tmp = -delta/A;                        coef[4][0][1] = -tmp;
    coef[4][1][0] = tmp = 1.0/A;                          coef[4][1][1] = -tmp;
    // m=2
    coef[5][0][0] = tmp = 1.0/(2.0*A);                   coef[5][0][1] = tmp;
    coef[5][1][0] = tmp = -1.0/(2.0*A);                  coef[5][1][1] = tmp;
}


void grt_static_source_coef_SH(MYREAL k, MYCOMPLEX coef[GRT_SRC_M_NUM][2])
{
    // 先全部赋0 
    for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
        for(MYINT p=0; p<2; ++p){
            coef[i][p] = 0.0;
        }
    }

    MYCOMPLEX tmp;

    // 水平力源
    coef[2][0] = tmp = -1.0/k;                coef[2][1] = tmp;

    // 剪切位错
    // m=1
    coef[4][0] = tmp = 1.0;                            coef[4][1] = -tmp;
    // m=2
    coef[5][0] = tmp = -1.0;                           coef[5][1] = tmp;
}


