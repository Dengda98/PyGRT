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

#include "static/static_source.h"
#include "common/const.h"


void static_source_coef(
    MYCOMPLEX delta, MYREAL k,
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

    MYCOMPLEX tmp;
    MYCOMPLEX A = RONE+delta;

    if(EXP!=NULL){
    EXP[0][0][0] = tmp = (delta-RONE)/A;         EXP[0][0][1] = tmp;    
    }

    if(VF!=NULL){
    VF[0][0][0] = tmp = -RONE/(RTWO*A*k);        VF[0][0][1] = - tmp;   
    VF[0][1][0] = tmp;                           VF[0][1][1] = - tmp;
    }

    if(HF!=NULL){
    HF[1][0][0] = tmp = RONE/(RTWO*A*k);        HF[1][0][1] = tmp;   
    HF[1][1][0] = - tmp;                        HF[1][1][1] = - tmp;
    HF[1][2][0] = tmp = -RONE/k;                HF[1][2][1] = tmp;
    }


    if(DC!=NULL){
    // m=0
    DC[0][0][0] = tmp = (-RONE+RFOUR*delta)/(RTWO*A);    DC[0][0][1] = tmp;
    DC[0][1][0] = tmp = -RTHREE/(RTWO*A);                DC[0][1][1] = tmp;
    // m=1
    DC[1][0][0] = tmp = -delta/A;                        DC[1][0][1] = -tmp;
    DC[1][1][0] = tmp = RONE/A;                          DC[1][1][1] = -tmp;
    DC[1][2][0] = tmp = RONE;                            DC[1][2][1] = -tmp;
    // m=2
    DC[2][0][0] = tmp = RONE/(RTWO*A);                   DC[2][0][1] = tmp;
    DC[2][1][0] = tmp = -RONE/(RTWO*A);                  DC[2][1][1] = tmp;
    DC[2][2][0] = tmp = -RONE;                           DC[2][2][1] = tmp;
    }
}


