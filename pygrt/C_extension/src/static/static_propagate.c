/**
 * @file   static_propagate.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 * 以下代码实现的是 静态广义反射透射系数矩阵 ，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *              透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 *
 */


#include <stdio.h>
#include <complex.h>
#include <string.h>

#include "grt/static/static_propagate.h"
#include "grt/static/static_layer.h"
#include "grt/static/static_source.h"
#include "grt/common/RT_matrix.h"
#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/common/matrix.h"


void grt_static_kernel(
    const GRT_MODEL1D *mod1d, MYCOMPLEX QWV[GRT_SRC_M_NUM][GRT_QWV_NUM],
    bool calc_uiz, MYCOMPLEX QWV_uiz[GRT_SRC_M_NUM][GRT_QWV_NUM])
{   
    // 初始化qwv为0
    for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
        for(MYINT j=0; j<GRT_QWV_NUM; ++j){
            QWV[i][j] = 0.0;
            if(calc_uiz)  QWV_uiz[i][j] = 0.0;
        }
    }


    bool ircvup = mod1d->ircvup;
    MYINT isrc = mod1d->isrc; // 震源所在虚拟层位, isrc>=1
    MYINT ircv = mod1d->ircv; // 接收点所在虚拟层位, ircv>=1, ircv != isrc
    MYINT imin, imax; // 相对浅层深层层位
    imin = GRT_MIN(isrc, ircv);
    imax = GRT_MAX(isrc, ircv);

    // 初始化广义反射透射系数矩阵
    // BL
    grt_init_RT_matrix(M_BL);
    // AL
    grt_init_RT_matrix(M_AL);
    // RS
    grt_init_RT_matrix(M_RS);
    // FA (实际先计算ZA，再递推到FA)
    grt_init_RT_matrix(M_FA);
    // FB (实际先计算ZB，再递推到FB)
    grt_init_RT_matrix(M_FB);

    // 自由表面的反射系数
    grt_init_RT_matrix(M_top);

    // 接收点处的接收矩阵
    MYCOMPLEX R_EV[2][2], R_EVL;
    
    // 接收点处的接收矩阵(转为位移导数ui_z的(B_m, C_m, P_m)系分量)
    MYCOMPLEX uiz_R_EV[2][2], uiz_R_EVL;

    // 从顶到底进行矩阵递推, 公式(5.5.3)
    for(MYINT iy=1; iy<mod1d->n; ++iy){ // 因为n>=3, 故一定会进入该循环

        // 定义物理层内的反射透射系数矩阵
        grt_init_RT_matrix(M);

        // 这里和动态解情况不同，即使是震源层、接收层这种虚拟层位也需要显式计算R/T矩阵
        // 对第iy层的系数矩阵赋值
        grt_static_RT_matrix_PSV(mod1d, iy, M);
        grt_static_RT_matrix_SH(mod1d, iy, M);

        // 加入时间延迟因子(第iy-1界面与第iy界面之间)
        grt_static_delay_RT_matrix(mod1d, iy, M);

        // FA
        if(iy <= imin){
            if(iy == 1){ // 初始化FA
                memcpy(M_FA, M, sizeof(*M));
            } else { // 递推FA
                grt_recursion_RT_matrix(M_FA, M, M_FA);  // 写入原地址 
            }
        }
        // RS
        else if(iy <= imax){
            if(iy == imin+1){// 初始化RS
                memcpy(M_RS, M, sizeof(*M));
            } else { // 递推RS
                grt_recursion_RT_matrix(M_RS, M, M_RS);  // 写入原地址
            }
        } 
        // BL
        else {
            if(iy == imax+1){// 初始化BL
                memcpy(M_BL, M, sizeof(*M));
            } else { // 递推BL
                // 只有 RD 矩阵最终会被使用到
                grt_recursion_RT_matrix(M_BL, M, M_BL);  // 写入原地址
            }
        } // END if

    } // END for loop
    //===================================================================================


    // 计算震源系数
    MYCOMPLEX src_coef_PSV[GRT_SRC_M_NUM][GRT_QWV_NUM-1][2] = {0};
    MYCOMPLEX src_coef_SH[GRT_SRC_M_NUM][2] = {0};
    grt_static_source_coef_PSV(mod1d, src_coef_PSV);
    grt_static_source_coef_SH(mod1d, src_coef_SH);

    // 临时中转矩阵 (temperary)
    MYCOMPLEX tmpR2[2][2], tmp2x2[2][2], tmpRL, tmp2x2_uiz[2][2], tmpRL_uiz;

    // 递推RU_FA
    grt_static_topfree_RU(mod1d, M_top);
    grt_recursion_RU(M_top, M_FA, M_FA);

    // 根据震源和台站相对位置，计算最终的系数
    if(ircvup){ // A接收  B震源
        // 计算R_EV
        grt_static_wave2qwv_REV_PSV(mod1d, M_FA->RU, R_EV);
        grt_static_wave2qwv_REV_SH(M_FA->RUL, &R_EVL);

        // 递推RU_FS
        grt_recursion_RU(M_FA, M_RS, M_FB); // 已从ZR变为FR，加入了自由表面的效应
        
        // 公式(5.7.12-14)
        grt_cmat2x2_mul(M_BL->RD, M_FB->RU, tmpR2);
        grt_cmat2x2_one_sub(tmpR2);
        grt_cmat2x2_inv(tmpR2, tmpR2, &M_FB->stats);// (I - xx)^-1
        grt_cmat2x2_mul(M_FB->invT, tmpR2, tmp2x2);

        if(calc_uiz) grt_cmat2x2_assign(tmp2x2, tmp2x2_uiz); // 为后续计算空间导数备份

        grt_cmat2x2_mul(R_EV, tmp2x2, tmp2x2);
        tmpRL = R_EVL * M_FB->invTL  / (1.0 - M_BL->RDL * M_FB->RUL);

        for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
            grt_construct_qwv(ircvup, tmp2x2, tmpRL, M_BL->RD, M_BL->RDL, src_coef_PSV[i], src_coef_SH[i], QWV[i]);
        }
        

        if(calc_uiz){
            grt_static_wave2qwv_z_REV_PSV(mod1d, M_FA->RU, uiz_R_EV);
            grt_static_wave2qwv_z_REV_SH(mod1d, M_FA->RUL, &uiz_R_EVL);
            grt_cmat2x2_mul(uiz_R_EV, tmp2x2_uiz, tmp2x2_uiz);
            tmpRL_uiz = tmpRL / R_EVL * uiz_R_EVL;
            
            for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
                grt_construct_qwv(ircvup, tmp2x2_uiz, tmpRL_uiz, M_BL->RD, M_BL->RDL, src_coef_PSV[i], src_coef_SH[i], QWV_uiz[i]);
            }    
        }
    }
    else { // A震源  B接收

        // 计算R_EV
        grt_static_wave2qwv_REV_PSV(mod1d, M_BL->RD, R_EV);    
        grt_static_wave2qwv_REV_SH(M_BL->RDL, &R_EVL);    

        // 递推RD_SL
        grt_recursion_RD(M_RS, M_BL, M_AL);
        
        // 公式(5.7.26-27)
        grt_cmat2x2_mul(M_FA->RU, M_AL->RD, tmpR2);
        grt_cmat2x2_one_sub(tmpR2);
        grt_cmat2x2_inv(tmpR2, tmpR2, &M_AL->stats);// (I - xx)^-1
        grt_cmat2x2_mul(M_AL->invT, tmpR2, tmp2x2);
        
        if(calc_uiz) grt_cmat2x2_assign(tmp2x2, tmp2x2_uiz); // 为后续计算空间导数备份

        grt_cmat2x2_mul(R_EV, tmp2x2, tmp2x2);
        tmpRL = R_EVL * M_AL->invTL / (1.0 - M_FA->RUL * M_AL->RDL);
        
        for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
            grt_construct_qwv(ircvup, tmp2x2, tmpRL, M_FA->RU, M_FA->RUL, src_coef_PSV[i], src_coef_SH[i], QWV[i]);
        }

        if(calc_uiz){
            grt_static_wave2qwv_z_REV_PSV(mod1d, M_BL->RD, uiz_R_EV);    
            grt_static_wave2qwv_z_REV_SH(mod1d, M_BL->RDL, &uiz_R_EVL);    
            grt_cmat2x2_mul(uiz_R_EV, tmp2x2_uiz, tmp2x2_uiz);
            tmpRL_uiz = tmpRL / R_EVL * uiz_R_EVL;
            
            for(MYINT i=0; i<GRT_SRC_M_NUM; ++i){
                grt_construct_qwv(ircvup, tmp2x2_uiz, tmpRL_uiz, M_FA->RU, M_FA->RUL, src_coef_PSV[i], src_coef_SH[i], QWV_uiz[i]);
            }
        }
    } // END if
}
