/**
 * @file   grt_static.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 *    计算静态位移
 * 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "const.h"
#include "model.h"
#include "colorstr.h"
#include "iostats.h"
#include "propagate.h"
#include "static_propagate.h"

#define DEG1 0.017453292519943295

static char *command;
static PYMODEL1D *pymod;
static double depsrc, deprcv;

MYREAL static_discrete_integ(
    const MODEL1D *mod1d, MYREAL dk, MYREAL kmax, MYREAL keps, 
    MYINT nr, MYREAL *rs,
    MYCOMPLEX sum_EXP_J[nr][3][4], MYCOMPLEX sum_VF_J[nr][3][4],  
    MYCOMPLEX sum_HF_J[nr][3][4],  MYCOMPLEX sum_DC_J[nr][3][4],  
    FILE *(fstats[nr]))
{
    MYCOMPLEX EXP_J[3][4], VF_J[3][4], HF_J[3][4],  DC_J[3][4];

    // 不同震源的核函数 F(k, w) 
    // 第一个维度3代表阶数m=0,1,2，第二个维度3代表三类系数qm,wm,vm 
    // 实际上对于不同震源只有特定阶数/系数才有值，不需要建立3x3的小矩阵，
    // 但这里还是为了方便可读性，牺牲了部分性能 
    MYCOMPLEX EXP_qwv[3][3], VF_qwv[3][3], HF_qwv[3][3], DC_qwv[3][3]; 
    MYCOMPLEX (*pEXP_qwv)[3] = (sum_EXP_J!=NULL)? EXP_qwv : NULL;
    MYCOMPLEX (*pVF_qwv)[3]  = (sum_VF_J!=NULL)?  VF_qwv  : NULL;
    MYCOMPLEX (*pHF_qwv)[3]  = (sum_HF_J!=NULL)?  HF_qwv  : NULL;
    MYCOMPLEX (*pDC_qwv)[3]  = (sum_DC_J!=NULL)?  DC_qwv  : NULL;
    
    MYREAL k = 0.0;
    MYINT ik = 0;

    // 所有震中距的k循环是否结束
    bool iendk = true;

    // 每个震中距的k循环是否结束
    bool iendkrs[nr], iendk0=false;
    for(MYINT ir=0; ir<nr; ++ir) iendkrs[ir] = false;
    

    // 波数k循环 (5.9.2)
    while(true){
        k += dk; 

        if(k > kmax)  break;

        // printf("w=%15.5e, ik=%d\n", CREAL(omega), ik);
        // 计算核函数 F(k, w)
        static_kernel(mod1d, k, pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv); 


        // 震中距rs循环
        iendk = true;
        for(MYINT ir=0; ir<nr; ++ir){
            if(iendkrs[ir]) continue; // 该震中距下的波数k积分已收敛

            for(MYINT m=0; m<3; ++m){
                for(MYINT v=0; v<4; ++v){
                    EXP_J[m][v] = VF_J[m][v] = HF_J[m][v] = DC_J[m][v] = CZERO;
                }
            }
            
            // 计算被积函数一项 F(k,w)Jm(kr)k
            int_Pk(k, rs[ir], 
                    pEXP_qwv, pVF_qwv, pHF_qwv, pDC_qwv,
                    EXP_J, VF_J, HF_J, DC_J);
            
            // printf("k=%f, r=%f, ", k, rs[ir]);
            // for(int i=0; i<3; ++i){
            //     for(int j=0; j<4; ++j){
            //         printf("%.5e%+.5ej ", CREAL(sum_DC_J[ir][i][j]), CIMAG(sum_DC_J[ir][i][j]));
            //     }
            // }
            // printf("\n");

            // 记录积分结果
            if(fstats[ir]!=NULL){
                write_stats(
                    fstats[ir], k, 
                    EXP_qwv, VF_qwv, HF_qwv, DC_qwv,
                    EXP_J, VF_J, HF_J, DC_J);
            }

            
            iendk0 = true;
            for(MYINT m=0; m<3; ++m){
                for(MYINT v=0; v<4; ++v){
                    if(sum_EXP_J!=NULL) sum_EXP_J[ir][m][v] += EXP_J[m][v];
                    if(sum_VF_J!=NULL)  sum_VF_J[ir][m][v]  += VF_J[m][v];
                    if(sum_HF_J!=NULL)  sum_HF_J[ir][m][v]  += HF_J[m][v];
                    if(sum_DC_J!=NULL)  sum_DC_J[ir][m][v]  += DC_J[m][v];

                    if(keps > RZERO){
                        // 判断是否达到收敛条件
                        if(sum_EXP_J!=NULL && m==0 && (v==0||v==2)) iendk0 = iendk0 && (CABS(EXP_J[m][v])/ CABS(sum_EXP_J[ir][m][v]) <= keps);
                        if(sum_VF_J!=NULL  && m==0 && (v==0||v==2)) iendk0 = iendk0 && (CABS(VF_J[m][v]) / CABS(sum_VF_J[ir][m][v])  <= keps);
                        if(sum_HF_J!=NULL && m==1)                  iendk0 = iendk0 && (CABS(HF_J[m][v]) / CABS(sum_HF_J[ir][m][v])  <= keps);
                        if(sum_DC_J!=NULL && ((m==0 && (v==0||v==2)) || m!=0)) iendk0 = iendk0 && (CABS(DC_J[m][v]) / CABS(sum_DC_J[ir][m][v])  <= keps);
                    } 
                }
            }
            
            if(keps > RZERO){
                iendkrs[ir] = iendk0;
                iendk = iendk && iendkrs[ir];
            } else {
                iendk = iendkrs[ir] = false;
            }
            
            

        } // END rs loop

        ++ik;

        // 所有震中距的格林函数都已收敛
        if(iendk) break;

    } // END k loop

    // printf("w=%15.5e, ik=%d\n", CREAL(omega), ik);

    return k;

}


int run(){

    depsrc = 2.0;
    deprcv = 0.0;

    // 读入模型文件
    if((pymod = read_pymod_from_file(command, "mmm", depsrc, deprcv)) ==NULL){
        exit(EXIT_FAILURE);
    }

    // pymod1d -> mod1d
    MODEL1D *mod1d = init_mod1d(pymod->n);
    get_mod1d(pymod, mod1d);
    // print_mod1d(mod1d);

    int nr=1;
    MYREAL *rs = (MYREAL*)malloc(sizeof(MYREAL)*nr);
    MYREAL *ys = (MYREAL*)malloc(sizeof(MYREAL)*nr);
    double x = 1.0, y;
    for(int ir=0; ir<nr; ++ir){
        y = 4.0 + ir*0.5;
        rs[ir] = sqrt(x*x + y*y);
        ys[ir] = y;
    }
        

    // for(int i=0; i<nr; ++i){
    //     rs[i] = 0.0 + 0.05*i;
    // }
    // rs[0] += 1e-5;

    MYREAL Sigma = 0.1*0.1; // 0.1*0.1km^2
    MYREAL fac = RONE/(RFOUR*PI) * Sigma;

    MYREAL dk = PI2/(100*20);
    MYREAL kmax = 15*PI/2.0;
    MYCOMPLEX sum_EXP_J[nr][3][4], sum_VF_J[nr][3][4], sum_HF_J[nr][3][4], sum_DC_J[nr][3][4];
    MYCOMPLEX (*psum_EXP_J)[3][4] = sum_EXP_J;
    MYCOMPLEX (*psum_VF_J)[3][4]  = sum_VF_J;
    MYCOMPLEX (*psum_HF_J)[3][4]  = sum_HF_J;
    MYCOMPLEX (*psum_DC_J)[3][4]  = sum_DC_J;
    FILE *fstats[nr];
    for(int i=0; i<nr; ++i)  fstats[i] = NULL;

    static_discrete_integ(mod1d, dk, kmax, -1, nr, rs, psum_EXP_J, psum_VF_J, psum_HF_J, psum_DC_J, fstats);


    // 局部变量，将某个频点的格林函数谱临时存放
    MYCOMPLEX tmp_EXP[nr][2], tmp_VF[nr][2], tmp_HF[nr][3];
    MYCOMPLEX tmp_DD[nr][2], tmp_DS[nr][3], tmp_SS[nr][3];

    MYCOMPLEX mtmp0 = dk * fac; // 

    // 记录到格林函数结构体内
    for(MYINT ir=0; ir<nr; ++ir){
        merge_Pk(
            (psum_EXP_J!=NULL)? psum_EXP_J[ir] : NULL, 
            (psum_VF_J!=NULL)?  psum_VF_J[ir]  : NULL, 
            (psum_HF_J!=NULL)?  psum_HF_J[ir]  : NULL, 
            (psum_DC_J!=NULL)?  psum_DC_J[ir]  : NULL, 
            tmp_EXP[ir], tmp_VF[ir],  tmp_HF[ir], 
            tmp_DD[ir], tmp_DS[ir], tmp_SS[ir]);
        
        // printf("%f ", rs[ir]);
        MYCOMPLEX mtmp;
        for(MYINT ii=0; ii<3; ++ii) {
            if(ii<2){
                mtmp = mtmp0*tmp_DD[ir][ii]; // m=0 45-倾滑
                // printf("%.8e ", CREAL(mtmp));
                // printf("%.8e%+.8ej ", CREAL(mtmp), CIMAG(mtmp));
            }
            else{
                // printf("%.8e ", 0.0);
            }
            mtmp = mtmp0*tmp_DS[ir][ii]; // m=1 90-倾滑
            // printf("%.8e ", CREAL(mtmp));
                // printf("%.8e%+.8ej ", CREAL(mtmp), CIMAG(mtmp));
            mtmp = mtmp0*tmp_SS[ir][ii]; // m=2 走滑 
            // printf("%.8e ", CREAL(mtmp));
            // printf("%.8e%+.8ej ", CREAL(mtmp), CIMAG(mtmp));
        }
        // printf("\n");

        printf("%f ", ys[ir]);
        double strike = 0.0;
        double dip = 70.0;
        double rake = 0.0;
        double azrad = atan2(ys[ir], x);
        // printf("%f ", azrad/DEG1);
        if(azrad < 0.0) azrad+=PI2;

        double stkrad = strike*DEG1;
        double diprad = dip*DEG1;
        double rakrad = rake*DEG1;
        double therad = azrad - stkrad;
        double srak, crak, sdip, cdip, sdip2, cdip2, sthe, cthe, sthe2, cthe2;
        srak = sin(rakrad);     crak = cos(rakrad);
        sdip = sin(diprad);     cdip = cos(diprad);
        sdip2 = 2.0*sdip*cdip;  cdip2 = 2.0*cdip*cdip - 1.0;
        sthe = sin(therad);     cthe = cos(therad);
        sthe2 = 2.0*sthe*cthe;  cthe2 = 2.0*cthe*cthe - 1.0;

        double A0, A1, A2, A4, A5;
        double mult = 1.0;
        A0 = mult * (0.5*sdip2*srak);
        A1 = mult * (cdip*crak*cthe - cdip2*srak*sthe);
        A2 = mult * (0.5*sdip2*srak*cthe2 + sdip*crak*sthe2);
        A4 = mult * (- cdip2*srak*cthe - cdip*crak*sthe);
        A5 = mult * (sdip*crak*cthe2 - 0.5*sdip2*srak*sthe2);
        // printf("A0=%.5e, A1=%.5e, A2=%.5e, A4=%.5e, A5=%.5e\n", A0,A1,A2,A4,A5);

        MYCOMPLEX Z, R, T;
        // Z
        Z = A0*tmp_DD[ir][0] + A1*tmp_DS[ir][0] + A2*tmp_SS[ir][0];
        Z *= mtmp0;
        // R
        R = A0*tmp_DD[ir][1] + A1*tmp_DS[ir][1] + A2*tmp_SS[ir][1];
        R *= mtmp0;
        // T
        T =                   A4*tmp_DS[ir][2] + A5*tmp_SS[ir][2];
        T *= mtmp0;
        
        MYCOMPLEX X, Y;
        X = R*cthe - T*sthe;
        Y = R*sthe + T*cthe;

        printf("%.8e ", CREAL(X));
        printf("%.8e ", CREAL(Y));
        printf("%.8e ", CREAL(-Z));
        printf("\n");
    }

}





int main(int argc, char **argv){
    command = argv[0];

    run();
}