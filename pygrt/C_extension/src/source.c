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

#include "source.h"
#include "model.h"
#include "matrix.h"
#include "prtdbg.h"




void source_coef(
    MYCOMPLEX src_xa, MYCOMPLEX src_xb, MYCOMPLEX src_kaka, MYCOMPLEX src_kbkb, 
    MYCOMPLEX omega, MYREAL k,
    MYCOMPLEX EXP[3][3][2], MYCOMPLEX VF[3][3][2], MYCOMPLEX HF[3][3][2], MYCOMPLEX DC[3][3][2])
{
    // 先全部赋0 
    for(MYINT i=0; i<3; ++i){
        for(MYINT j=0; j<3; ++j){
            for(MYINT p=0; p<2; ++p){
                if(EXP!=NULL) EXP[i][j][p] = RZERO;
                if(VF!=NULL)  VF[i][j][p] = RZERO;
                if(HF!=NULL)  HF[i][j][p] = RZERO;
                if(DC!=NULL)  DC[i][j][p] = RZERO;
            }
        }
    }


    MYCOMPLEX src_a_inv = RONE / (k*src_xa);
    MYCOMPLEX src_b_inv = RONE / (k*src_xb);
    MYREAL kk = k*k;
    MYREAL k_inv = RONE / k;
    MYCOMPLEX tmp;

    if(EXP!=NULL){
    // 爆炸源， 通过(4.9.8)的矩张量源公式，提取各向同性的量(M11+M22+M33)，-a+k^2/a -> ka^2/a
    EXP[0][0][0] = tmp = src_kaka * src_a_inv;         EXP[0][0][1] = tmp;    
    }
    
    if(VF!=NULL){
    // 垂直力源 (4.6.15)
    VF[0][0][0] = tmp = -RONE;                         VF[0][0][1] = - tmp;
    VF[0][1][0] = tmp = -k * src_b_inv;                VF[0][1][1] = tmp;
    }

    if(HF!=NULL){
    // 水平力源 (4.6.21,26), 这里可以把x1,x2方向的力转到r,theta方向
    // 推导可发现，r方向的力形成P,SV波, theta方向的力形成SH波
    // 方向性因子包含水平力方向与震源台站连线方向的夹角
    HF[1][0][0] = tmp = -k * src_a_inv;                   HF[1][0][1] = tmp;
    HF[1][1][0] = tmp = -RONE;                            HF[1][1][1] = - tmp;
    HF[1][2][0] = tmp = src_kbkb * k_inv * src_b_inv;     HF[1][2][1] = tmp;
    }

    if(DC!=NULL){
    // 剪切位错 (4.8.34)
    // m=0
    DC[0][0][0] = tmp = (RTWO*src_kaka - RTHREE*kk) * src_a_inv;    DC[0][0][1] = tmp;
    DC[0][1][0] = tmp = -RTHREE*k;                                  DC[0][1][1] = - tmp;
    // m=1
    DC[1][0][0] = tmp = RTWO*k;                              DC[1][0][1] = - tmp;
    DC[1][1][0] = tmp = (RTWO*kk - src_kbkb) * src_b_inv;    DC[1][1][1] = tmp;
    DC[1][2][0] = tmp = - src_kbkb * k_inv;                  DC[1][2][1] = - tmp;

    // m=2
    DC[2][0][0] = tmp = - kk * src_a_inv;                    DC[2][0][1] = tmp;
    DC[2][1][0] = tmp = - k;                                 DC[2][1][1] = - tmp;
    DC[2][2][0] = tmp = src_kbkb * src_b_inv;                DC[2][2][1] = tmp;
    }


}




