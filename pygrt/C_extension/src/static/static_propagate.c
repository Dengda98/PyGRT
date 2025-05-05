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

#include "static/static_propagate.h"
#include "static/static_layer.h"
#include "static/static_source.h"
#include "common/recursion.h"
#include "common/model.h"
#include "common/const.h"
#include "common/matrix.h"

#define CMAT_ASSIGN_SPLIT 0  // 2x2的小矩阵赋值合并为1个循环，程序速度提升微小


void static_kernel(
    const MODEL1D *mod1d, MYCOMPLEX omega, MYREAL k, MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM],
    bool calc_uiz, MYCOMPLEX QWV_uiz[SRC_M_NUM][QWV_NUM], MYINT *stats)
{   
    // 初始化qwv为0
    for(MYINT i=0; i<SRC_M_NUM; ++i){
        for(MYINT j=0; j<QWV_NUM; ++j){
            QWV[i][j] = CZERO;
            if(calc_uiz)  QWV_uiz[i][j] = CZERO;
        }
    }


    bool ircvup = mod1d->ircvup;
    MYINT isrc = mod1d->isrc; // 震源所在虚拟层位, isrc>=1
    MYINT ircv = mod1d->ircv; // 接收点所在虚拟层位, ircv>=1, ircv != isrc
    MYINT imin, imax; // 相对浅层深层层位
    imin = mod1d->imin;
    imax = mod1d->imax;

    // 初始化广义反射透射系数矩阵
    // BL
    MYCOMPLEX RD_BL[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL_BL = CZERO;
    MYCOMPLEX RU_BL[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL_BL = CZERO;
    MYCOMPLEX TD_BL[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TDL_BL = CONE;
    MYCOMPLEX TU_BL[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TUL_BL = CONE;
    // AL
    MYCOMPLEX RD_AL[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL_AL = CZERO;
    // RS
    MYCOMPLEX RD_RS[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL_RS = CZERO;
    MYCOMPLEX RU_RS[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL_RS = CZERO;
    MYCOMPLEX TD_RS[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TDL_RS = CONE;
    MYCOMPLEX TU_RS[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TUL_RS = CONE;
    // FA (实际先计算ZA，再递推到FA)
    MYCOMPLEX RD_FA[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RDL_FA = CZERO;
    MYCOMPLEX RU_FA[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL_FA = CZERO;
    MYCOMPLEX TD_FA[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TDL_FA = CONE;
    MYCOMPLEX TU_FA[2][2] = INIT_C_IDENTITY_2x2_MATRIX;
    MYCOMPLEX TUL_FA = CONE;
    // FB (实际先计算ZB，再递推到FB)
    MYCOMPLEX RU_FB[2][2] = INIT_C_ZERO_2x2_MATRIX;
    MYCOMPLEX RUL_FB = CZERO;

    // 抽象指针 
    // BL
    MYCOMPLEX *const pRDL_BL = &RDL_BL;
    MYCOMPLEX *const pRUL_BL = &RUL_BL;
    MYCOMPLEX *const pTDL_BL = &TDL_BL;
    MYCOMPLEX *const pTUL_BL = &TUL_BL;
    // AL
    MYCOMPLEX *const pRDL_AL = &RDL_AL;
    // RS
    MYCOMPLEX *const pRDL_RS = &RDL_RS;
    MYCOMPLEX *const pRUL_RS = &RUL_RS;
    MYCOMPLEX *const pTDL_RS = &TDL_RS;
    MYCOMPLEX *const pTUL_RS = &TUL_RS;
    // FA
    MYCOMPLEX *const pRDL_FA = &RDL_FA;
    MYCOMPLEX *const pRUL_FA = &RUL_FA;
    MYCOMPLEX *const pTDL_FA = &TDL_FA;
    MYCOMPLEX *const pTUL_FA = &TUL_FA;
    // FB 
    MYCOMPLEX *const pRUL_FB = &RUL_FB;


    // 定义物理层内的反射透射系数矩阵，相对于界面上的系数矩阵增加了时间延迟因子
    MYCOMPLEX RD[2][2], RDL, TD[2][2], TDL;
    MYCOMPLEX RU[2][2], RUL, TU[2][2], TUL;
    MYCOMPLEX *const pRDL = &RDL;
    MYCOMPLEX *const pTDL = &TDL;
    MYCOMPLEX *const pRUL = &RUL;
    MYCOMPLEX *const pTUL = &TUL;


    // 自由表面的反射系数
    MYCOMPLEX R_tilt[2][2] = INIT_C_ZERO_2x2_MATRIX; // SH波在自由表面的反射系数为1，不必定义变量

    // 接收点处的接收矩阵
    MYCOMPLEX R_EV[2][2], R_EVL;
    MYCOMPLEX *const pR_EVL = &R_EVL;
    
    // 接收点处的接收矩阵(转为位移导数ui_z的(B_m, C_m, P_m)系分量)
    MYCOMPLEX uiz_R_EV[2][2], uiz_R_EVL;
    MYCOMPLEX *const puiz_R_EVL = &uiz_R_EVL;


    // 模型参数
    // 后缀0，1分别代表上层和下层
    LAYER *lay = NULL;
    MYREAL mod1d_thk0, mod1d_thk1;
    MYCOMPLEX mod1d_mu0, mod1d_mu1;
    MYCOMPLEX mod1d_delta0, mod1d_delta1;
    MYCOMPLEX top_delta = CZERO;
    MYCOMPLEX src_delta = CZERO;
    MYCOMPLEX rcv_delta = CZERO;
    

    // 从顶到底进行矩阵递推, 公式(5.5.3)
    for(MYINT iy=0; iy<mod1d->n; ++iy){ // 因为n>=3, 故一定会进入该循环
        lay = mod1d->lays + iy;

        // 赋值上层 
        mod1d_thk0 = mod1d_thk1;
        mod1d_mu0 = mod1d_mu1;
        mod1d_delta0 = mod1d_delta1;

        // 更新模型参数
        mod1d_thk1 = lay->thk;
        mod1d_mu1 = lay->mu;
        mod1d_delta1 = lay->delta;

        if(0==iy){
            top_delta = mod1d_delta1;
            continue;
        }

        // 确定上下层的物性参数
        if(ircv==iy){
            rcv_delta = mod1d_delta1;
        } else if(isrc==iy){
            src_delta = mod1d_delta1;
        }

        // 对第iy层的系数矩阵赋值，加入时间延迟因子(第iy-1界面与第iy界面之间)
        calc_static_RT_2x2(
            mod1d_delta0, mod1d_mu0,
            mod1d_delta1, mod1d_mu1,
            mod1d_thk0, k, // 使用iy-1层的厚度
            RD, pRDL, RU, pRUL, 
            TD, pTDL, TU, pTUL);

        // FA
        if(iy <= imin){
            if(iy == 1){ // 初始化FA
#if CMAT_ASSIGN_SPLIT == 1
                cmat2x2_assign(RD, RD_FA);  RDL_FA = RDL;
                cmat2x2_assign(RU, RU_FA);  RUL_FA = RUL;
                cmat2x2_assign(TD, TD_FA);  TDL_FA = TDL;
                cmat2x2_assign(TU, TU_FA);  TUL_FA = TUL;
#else 
                for(MYINT kk=0; kk<2; ++kk){
                    for(MYINT pp=0; pp<2; ++pp){
                        RD_FA[kk][pp] = RD[kk][pp];
                        RU_FA[kk][pp] = RU[kk][pp];
                        TD_FA[kk][pp] = TD[kk][pp];
                        TU_FA[kk][pp] = TU[kk][pp];
                    }
                }
                RDL_FA = RDL;
                RUL_FA = RUL;
                TDL_FA = TDL;
                TUL_FA = TUL;
#endif
            } else { // 递推FA

                recursion_RT_2x2(
                    RD_FA, RDL_FA, RU_FA, RUL_FA, 
                    TD_FA, TDL_FA, TU_FA, TUL_FA,
                    RD, RDL, RU, RUL, 
                    TD, TDL, TU, TUL,
                    RD_FA, pRDL_FA, RU_FA, pRUL_FA, 
                    TD_FA, pTDL_FA, TU_FA, pTUL_FA, stats);  
            }

        }
        // RS
        else if(iy <= imax){
            if(iy == imin+1){// 初始化RS
#if CMAT_ASSIGN_SPLIT == 1
                cmat2x2_assign(RD, RD_RS);  RDL_RS = RDL;
                cmat2x2_assign(RU, RU_RS);  RUL_RS = RUL;
                cmat2x2_assign(TD, TD_RS);  TDL_RS = TDL;
                cmat2x2_assign(TU, TU_RS);  TUL_RS = TUL;
#else
                for(MYINT kk=0; kk<2; ++kk){
                    for(MYINT pp=0; pp<2; ++pp){
                        RD_RS[kk][pp] = RD[kk][pp];
                        RU_RS[kk][pp] = RU[kk][pp];
                        TD_RS[kk][pp] = TD[kk][pp];
                        TU_RS[kk][pp] = TU[kk][pp];
                    }
                }
                RDL_RS = RDL;
                RUL_RS = RUL;
                TDL_RS = TDL;
                TUL_RS = TUL;
#endif
            } else { // 递推RS
                recursion_RT_2x2(
                    RD_RS, RDL_RS, RU_RS, RUL_RS, 
                    TD_RS, TDL_RS, TU_RS, TUL_RS,
                    RD, RDL, RU, RUL, 
                    TD, TDL, TU, TUL,
                    RD_RS, pRDL_RS, RU_RS, pRUL_RS, 
                    TD_RS, pTDL_RS, TU_RS, pTUL_RS, stats);  // 写入原地址
            }
        } 
        // BL
        else {
            if(iy == imax+1){// 初始化BL
#if CMAT_ASSIGN_SPLIT == 1
                cmat2x2_assign(RD, RD_BL);  RDL_BL = RDL;
                cmat2x2_assign(RU, RU_BL);  RUL_BL = RUL;
                cmat2x2_assign(TD, TD_BL);  TDL_BL = TDL;
                cmat2x2_assign(TU, TU_BL);  TUL_BL = TUL;
#else 
                for(MYINT kk=0; kk<2; ++kk){
                    for(MYINT pp=0; pp<2; ++pp){
                        RD_BL[kk][pp] = RD[kk][pp];
                        RU_BL[kk][pp] = RU[kk][pp];
                        TD_BL[kk][pp] = TD[kk][pp];
                        TU_BL[kk][pp] = TU[kk][pp];
                    }
                }
                RDL_BL = RDL;
                RUL_BL = RUL;
                TDL_BL = TDL;
                TUL_BL = TUL;
#endif
            } else { // 递推BL

                // 这个IF纯粹是为了优化，因为不论是SL还是RL，只有RD矩阵最终会被使用到
                // 这里最终只把RD矩阵的值记录下来，其它的舍去，以减少部分运算
                if(iy < mod1d->n - 1){
                    recursion_RT_2x2(
                        RD_BL, RDL_BL, RU_BL, RUL_BL, 
                        TD_BL, TDL_BL, TU_BL, TUL_BL,
                        RD, RDL, RU, RUL, 
                        TD, TDL, TU, TUL,
                        RD_BL, pRDL_BL, RU_BL, pRUL_BL, 
                        TD_BL, pTDL_BL, TU_BL, pTUL_BL, stats);  // 写入原地址
                } else {
                    recursion_RT_2x2(
                        RD_BL, RDL_BL, RU_BL, RUL_BL, 
                        TD_BL, TDL_BL, TU_BL, TUL_BL,
                        RD, RDL, RU, RUL, 
                        TD, TDL, TU, TUL,
                        RD_BL, pRDL_BL, NULL, NULL, 
                        NULL, NULL, NULL, NULL, stats);  // 写入原地址
                }
                
            }

        } // END if

    } // END for loop
    //===================================================================================


    // 计算震源系数
    MYCOMPLEX src_coef[SRC_M_NUM][QWV_NUM][2] = {0};
    static_source_coef(src_delta, k, src_coef);

    // 临时中转矩阵 (temperary)
    MYCOMPLEX tmpR2[2][2], tmp2x2[2][2], tmpRL, tmp2x2_uiz[2][2], tmpRL_uiz;
    MYCOMPLEX inv_2x2T[2][2], invT;

    // 递推RU_FA
    calc_static_R_tilt(top_delta, R_tilt);
    recursion_RU(
        R_tilt, RONE, 
        RD_FA, RDL_FA,
        RU_FA, RUL_FA, 
        TD_FA, TDL_FA,
        TU_FA, TUL_FA,
        RU_FA, pRUL_FA, NULL, NULL, stats);

    // 根据震源和台站相对位置，计算最终的系数
    if(ircvup){ // A接收  B震源
        // 计算R_EV
        calc_static_R_EV(ircvup, RU_FA, RUL_FA, R_EV, pR_EVL);

        // 递推RU_FS
        recursion_RU(
            RU_FA, RUL_FA, // 已从ZR变为FR，加入了自由表面的效应
            RD_RS, RDL_RS,
            RU_RS, RUL_RS, 
            TD_RS, TDL_RS,
            TU_RS, TUL_RS,
            RU_FB, pRUL_FB, inv_2x2T, &invT, stats);
        
        // 公式(5.7.12-14)
        cmat2x2_mul(RD_BL, RU_FB, tmpR2);
        cmat2x2_one_sub(tmpR2);
        cmat2x2_inv(tmpR2, tmpR2, stats);// (I - xx)^-1
        cmat2x2_mul(inv_2x2T, tmpR2, tmp2x2);

        if(calc_uiz) cmat2x2_assign(tmp2x2, tmp2x2_uiz); // 为后续计算空间导数备份

        cmat2x2_mul(R_EV, tmp2x2, tmp2x2);
        tmpRL = R_EVL * invT  / (RONE - RDL_BL * RUL_FB);

        for(MYINT i=0; i<SRC_M_NUM; ++i){
            get_qwv(ircvup, tmp2x2, tmpRL, RD_BL, RDL_BL, src_coef[i], QWV[i]);
        }
        

        if(calc_uiz){
            calc_static_uiz_R_EV(rcv_delta, ircvup, k, RU_FA, RUL_FA, uiz_R_EV, puiz_R_EVL);
            cmat2x2_mul(uiz_R_EV, tmp2x2_uiz, tmp2x2_uiz);
            tmpRL_uiz = tmpRL / R_EVL * uiz_R_EVL;
            
            for(MYINT i=0; i<SRC_M_NUM; ++i){
                get_qwv(ircvup, tmp2x2_uiz, tmpRL_uiz, RD_BL, RDL_BL, src_coef[i], QWV_uiz[i]);
            }    
        }
    }
    else { // A震源  B接收

        // 计算R_EV
        calc_static_R_EV(ircvup, RD_BL, RDL_BL, R_EV, pR_EVL);    

        // 递推RD_SL
        recursion_RD(
            RD_RS, RDL_RS,
            RU_RS, RUL_RS,
            TD_RS, TDL_RS,
            TU_RS, TUL_RS,
            RD_BL, RDL_BL,
            RD_AL, pRDL_AL, inv_2x2T, &invT, stats);
        
        // 公式(5.7.26-27)
        cmat2x2_mul(RU_FA, RD_AL, tmpR2);
        cmat2x2_one_sub(tmpR2);
        cmat2x2_inv(tmpR2, tmpR2, stats);// (I - xx)^-1
        cmat2x2_mul(inv_2x2T, tmpR2, tmp2x2);
        
        if(calc_uiz) cmat2x2_assign(tmp2x2, tmp2x2_uiz); // 为后续计算空间导数备份

        cmat2x2_mul(R_EV, tmp2x2, tmp2x2);
        tmpRL = R_EVL * invT / (RONE - RUL_FA * RDL_AL);
        
        for(MYINT i=0; i<SRC_M_NUM; ++i){
            get_qwv(ircvup, tmp2x2, tmpRL, RU_FA, RUL_FA, src_coef[i], QWV[i]);
        }

        if(calc_uiz){
            calc_static_uiz_R_EV(rcv_delta, ircvup, k, RD_BL, RDL_BL, uiz_R_EV, puiz_R_EVL);    
            cmat2x2_mul(uiz_R_EV, tmp2x2_uiz, tmp2x2_uiz);
            tmpRL_uiz = tmpRL / R_EVL * uiz_R_EVL;
            
            for(MYINT i=0; i<SRC_M_NUM; ++i){
                get_qwv(ircvup, tmp2x2_uiz, tmpRL_uiz, RU_FA, RUL_FA, src_coef[i], QWV_uiz[i]);
            }
        }
    } // END if
}
