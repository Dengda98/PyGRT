/**
 * @file   integral.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 *     将被积函数的逐点值累加成积分值
 *                   
 */


#include <stdio.h>
#include <stdbool.h>

#include "common/integral.h"
#include "common/const.h"
#include "common/bessel.h"



void int_Pk(
    MYREAL k, MYREAL r, 
    // F(ki,w)， 第一个维度表示不同震源，不同阶数，第二个维度3代表三类系数qm,wm,vm 
    const MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM],
    // F(ki,w)Jm(ki*r)ki，第一个维度表示不同震源，不同阶数，第二个维度代表4种类型的F(k,w)Jm(kr)k的类型
    bool calc_uir,
    MYCOMPLEX SUM[SRC_M_NUM][INTEG_NUM])
{
    MYREAL bjmk[MORDER_MAX+1] = {0};
    MYREAL kr = k*r;
    MYREAL kr_inv = RONE/kr;
    MYREAL kcoef = k;

    MYREAL Jmcoef[MORDER_MAX+1] = {0};

    bessel012(kr, &bjmk[0], &bjmk[1], &bjmk[2]); 
    if(calc_uir){
        MYREAL bjmk0[MORDER_MAX+1] = {0};
        for(MYINT i=0; i<=MORDER_MAX; ++i)  bjmk0[i] = bjmk[i];

        besselp012(kr, &bjmk[0], &bjmk[1], &bjmk[2]); 
        kcoef = k*k;

        for(MYINT i=1; i<=MORDER_MAX; ++i)  Jmcoef[i] = kr_inv * (-kr_inv * bjmk0[i] + bjmk[i]);
    } 
    else {
        for(MYINT i=1; i<=MORDER_MAX; ++i)  Jmcoef[i] = bjmk[i]*kr_inv;
    }

    for(MYINT i=1; i<=MORDER_MAX; ++i)  Jmcoef[i] *= kcoef;
    for(MYINT i=0; i<=MORDER_MAX; ++i)  bjmk[i] *= kcoef;


    // 公式(5.6.22), 将公式分解为F(k,w)Jm(kr)k的形式
    for(MYINT i=0; i<SRC_M_NUM; ++i){
        MYINT modr = SRC_M_ORDERS[i];  // 对应m阶数
        if(modr == 0){
            SUM[i][0] = - QWV[i][0] * bjmk[1];   // - q0*J1*k
            SUM[i][2] =   QWV[i][1] * bjmk[0];   //   w0*J0*k
        }
        else{
            SUM[i][0]  =   QWV[i][0] * bjmk[modr-1];         // qm*Jm-1*k
            SUM[i][1]  = - modr*(QWV[i][0] + QWV[i][2]) * Jmcoef[modr];    // - m*(qm+vm)*Jm*k/kr
            SUM[i][2]  =   QWV[i][1] * bjmk[modr];           // wm*Jm*k
            SUM[i][3]  = - QWV[i][2] * bjmk[modr-1];         // -vm*Jm-1*k
        }
    }

}


void int_Pk_filon(
    MYREAL k, MYREAL r, bool iscos,
    // F(ki,w)， 第一个维度表示不同震源，不同阶数，第二个维度3代表三类系数qm,wm,vm 
    const MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM],
    // F(ki,w)Jm(ki*r)ki，第一个维度表示不同震源，不同阶数，第二个维度代表4种类型的F(k,w)Jm(kr)k的类型
    bool calc_uir,
    MYCOMPLEX SUM[SRC_M_NUM][INTEG_NUM])
{
    MYREAL phi0 = 0.0;
    if(! iscos)  phi0 = - HALFPI;  // 在cos函数中添加的相位差，用于计算sin函数

    MYREAL kr = k*r;
    MYREAL kr_inv = RONE/kr;
    MYREAL kcoef = SQRT(k);
    MYREAL bjmk[MORDER_MAX+1] = {0};

    MYREAL Jmcoef[MORDER_MAX+1] = {0};

    if(calc_uir){
        kcoef *= k;
        // 使用bessel递推公式 Jm'(x) = m/x * Jm(x) - J_{m+1}(x)
        // 考虑大震中距，忽略第一项，再使用bessel渐近公式
        bjmk[0] = - COS(kr - THREEQUARTERPI - phi0);
        bjmk[1] = - COS(kr - FIVEQUARTERPI - phi0);
        bjmk[2] = - COS(kr - SEVENQUARTERPI - phi0);
    } else {
        bjmk[0] = COS(kr - QUARTERPI - phi0);
        bjmk[1] = COS(kr - THREEQUARTERPI - phi0);
        bjmk[2] = COS(kr - FIVEQUARTERPI - phi0);
    }

    for(MYINT i=1; i<=MORDER_MAX; ++i)  Jmcoef[i] = bjmk[i]*kr_inv*kcoef;
    for(MYINT i=0; i<=MORDER_MAX; ++i)  bjmk[i] *= kcoef;
    
    // 公式(5.6.22), 将公式分解为F(k,w)Jm(kr)k的形式
    for(MYINT i=0; i<SRC_M_NUM; ++i){
        MYINT modr = SRC_M_ORDERS[i];  // 对应m阶数
        if(modr == 0){
            SUM[i][0] = - QWV[i][0] * bjmk[1];   // - q0*J1*k
            SUM[i][2] =   QWV[i][1] * bjmk[0];   //   w0*J0*k
        }
        else{
            SUM[i][0]  =   QWV[i][0] * bjmk[modr-1];         // qm*Jm-1*k
            SUM[i][1]  = - modr*(QWV[i][0] + QWV[i][2]) * Jmcoef[modr];    // - m*(qm+vm)*Jm*k/kr
            SUM[i][2]  =   QWV[i][1] * bjmk[modr];           // wm*Jm*k
            SUM[i][3]  = - QWV[i][2] * bjmk[modr-1];         // -vm*Jm-1*k
        }
    }
}





void merge_Pk(
    // F(ki,w)Jm(ki*r)ki，
    const MYCOMPLEX sum_J[SRC_M_NUM][INTEG_NUM], 
    // 累积求和，Z、R、T分量 
    MYCOMPLEX tol[SRC_M_NUM][CHANNEL_NUM])
{   
    for(MYINT i=0; i<SRC_M_NUM; ++i){
        MYINT modr = SRC_M_ORDERS[i];
        if(modr==0){
            tol[i][0] = sum_J[i][2];
            tol[i][1] = sum_J[i][0];
        }
        else {
            tol[i][0] = sum_J[i][2];
            tol[i][1] = sum_J[i][0] + sum_J[i][1];
            tol[i][2] = - sum_J[i][1] + sum_J[i][3];
        }
    }
}
