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

#include "dynamic/source.h"
#include "common/model.h"
#include "common/matrix.h"
#include "common/prtdbg.h"




void source_coef(
    MYCOMPLEX src_xa, MYCOMPLEX src_xb, MYCOMPLEX src_kaka, MYCOMPLEX src_kbkb, 
    MYREAL k,
    MYCOMPLEX coef[SRC_M_NUM][QWV_NUM][2])
{
    // 先全部赋0 
    for(MYINT i=0; i<SRC_M_NUM; ++i){
        for(MYINT j=0; j<QWV_NUM; ++j){
            for(MYINT p=0; p<2; ++p){
                coef[i][j][p] = CZERO;
            }
        }
    }


    MYCOMPLEX a = k*src_xa;
    MYCOMPLEX b = k*src_xb;
    MYREAL kk = k*k;
    MYCOMPLEX tmp;


    // 爆炸源， 通过(4.9.8)的矩张量源公式，提取各向同性的量(M11+M22+M33)，-a+k^2/a -> ka^2/a
    coef[0][0][0] = tmp = src_kaka / a;         coef[0][0][1] = tmp;    
    
    // 垂直力源 (4.6.15)
    coef[1][0][0] = tmp = -RONE;                 coef[1][0][1] = - tmp;
    coef[1][1][0] = tmp = -k / b;                coef[1][1][1] = tmp;

    // 水平力源 (4.6.21,26), 这里可以把x1,x2方向的力转到r,theta方向
    // 推导可发现，r方向的力形成P,SV波, theta方向的力形成SH波
    // 方向性因子包含水平力方向与震源台站连线方向的夹角
    coef[2][0][0] = tmp = -k / a;              coef[2][0][1] = tmp;
    coef[2][1][0] = tmp = -RONE;               coef[2][1][1] = - tmp;
    coef[2][2][0] = tmp = src_kbkb / k / b;    coef[2][2][1] = tmp;

    // 剪切位错 (4.8.34)
    // m=0
    coef[3][0][0] = tmp = (RTWO*src_kaka - RTHREE*kk) / a;    coef[3][0][1] = tmp;
    coef[3][1][0] = tmp = -RTHREE*k;                          coef[3][1][1] = - tmp;
    // m=1
    coef[4][0][0] = tmp = RTWO*k;                      coef[4][0][1] = - tmp;
    coef[4][1][0] = tmp = (RTWO*kk - src_kbkb) / b;    coef[4][1][1] = tmp;
    coef[4][2][0] = tmp = - src_kbkb / k;              coef[4][2][1] = - tmp;

    // m=2
    coef[5][0][0] = tmp = - kk / a;                    coef[5][0][1] = tmp;
    coef[5][1][0] = tmp = - k;                         coef[5][1][1] = - tmp;
    coef[5][2][0] = tmp = src_kbkb / b;                coef[5][2][1] = tmp;


}




