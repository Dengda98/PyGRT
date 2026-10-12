/**
 * @file   safim.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-4-27
 * 
 * 以下代码实现的是自适应Filon积分，参考：
 * 
 *         Chen, X., and H. Zhang (2001). An Efficient Method for Computing Green’s Functions 
 *         for a Layered Half-Space at Large Epicentral Distances, Bulletin of the Seismological 
 *         Society of America 91, no. 4, 858–869, doi: 10.1785/0120000113.
 * 
 */


#include <stdio.h> 
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "grt/integral/safim.h"
#include "grt/integral/k_integ.h"
#include "grt/integral/iostats.h"
#include "grt/common/const.h"
#include "grt/common/model.h"


/**
 * 用于判断拟合情况，当|F(k)| < s * max{|F(k)|}，调整该区间的F(k)积分和，
 * 防止其在小幅值处采样过多。
 */
#define REF_AMP_SCALE 1e-6 

/**
 * 自适应划分区间的最小dk
 */
#define SA_MIN_DK 1e-7


/** 辛普森积分的区间范围 */
typedef enum {
    SIMPSON_As2H_A = -1,
    SIMPSON_A_Ap2H = 1,
    SIMPSON_Ap2H_Ap4H = 2,
    SIMPSON_A_ApH = 3,
    SIMPSON_ApH_Ap2H = 4
} SIMPSON_INTV;




// 区间结构体
typedef struct { 
    real_t k3[3];
    cplxChnlGrid F3[3]; 
    cplxChnlGrid Fz3[3]; 
    cplxChnlGrid F3_raw[3];
    cplxChnlGrid Fz3_raw[3];
} KInterval;

// 区间栈结构体
typedef struct { 
    KInterval *data; 
    size_t size; 
    size_t capacity; 
} KIntervalStack;

/** 初始化区间栈 */
static void stack_init(KIntervalStack *stack, size_t init_capacity) {
    stack->data = GRT_SAFE_MALLOC(init_capacity * sizeof(KInterval));
    stack->size = 0;
    stack->capacity = init_capacity;
}

/** 从栈顶部压入元素 */
static void stack_push(KIntervalStack *stack, KInterval item) {
    // 扩容
    if(stack->size >= stack->capacity){
        stack->capacity *= 2;
        stack->data = GRT_SAFE_REALLOC(stack->data, stack->capacity * sizeof(KInterval));
    }
    stack->data[stack->size++] = item;
}

/** 从栈顶部拿走元素 */
static KInterval stack_pop(KIntervalStack *stack) {
    if(stack->size == 0) {
        GRTRaiseError("Pop from empty stack\n");
    }
    return stack->data[--stack->size];
}




/** 
 * 三点辛普森积分计算，基于区间内三个点的函数值，定义二次曲线，再根据stats定义积分区间,
 * itype 的 0-2 对应 q/w/v 普通项，3-4 对应近场振幅及其径向权重
 * 假设区间内三点为[a, a+h, a+2h]，stats的取值对应的积分区间为
 *    + stats = -1,  [a-2h, a]
 *    + stats = 1,   [a, a+2h]
 *    + stats = 2,   [a+2h, a+4h]
 *    + stats = 3,   [a, a+h]
 *    + stats = 4,   [a+h, a+2h]
 * 
 */
static cplx_t simpson(const KInterval *item_pt, int im, int itype, bool isuiz, SIMPSON_INTV stats) {
    cplx_t Fint = 0.0;
    real_t klen = item_pt->k3[2] -  item_pt->k3[0];
    const cplxChnlGrid *F3 = (isuiz)? item_pt->Fz3 : item_pt->F3;
    
    // 普通项使用 F*sqrt(k)，额外两个通道分别检查近场振幅 H/sqrt(k) 及其径向权重 H*sqrt(k)
    cplx_t F[3];
    for(int i=0; i<3; ++i){
        real_t sk = sqrt(item_pt->k3[i]);
        if(itype < GRT_CHANNEL_NUM){
            F[i] = F3[i][im][itype] * sk;
        } else {
            cplx_t H = F3[i][im][0] + F3[i][im][2];
            F[i] = (itype == GRT_CHANNEL_NUM)? H / sk : H * sk;
        }
    }

    if(stats == SIMPSON_A_Ap2H){
        Fint = klen * (F[0] + 4.0*F[1] + F[2]) / 6.0;
    }
    else if(stats == SIMPSON_As2H_A){
        Fint = klen * (19.0*F[0] - 20.0*F[1] + 7.0*F[2]) / 6.0;
    }
    else if(stats == SIMPSON_Ap2H_Ap4H){
        Fint = klen * (7.0*F[0] - 20.0*F[1] + 19.0*F[2]) / 6.0;
    }
    else if(stats == SIMPSON_A_ApH){
        Fint = klen * (5.0*F[0] + 8.0*F[1] - F[2]) / 24.0;
    }
    else if(stats == SIMPSON_ApH_Ap2H){
        Fint = klen * ( - F[0] + 8.0*F[1] + 5.0*F[2]) / 24.0;
    }
    else{
        GRTRaiseError("wrong simpson stats (%d).\n", stats);
    }

    return Fint;
}

/** 更新普通核函数及近场组合 H=qm+vm 的幅值参考 */
static void update_maxabs(const cplxChnlGrid F, real_t maxabsF[GRT_GTYPES_MAX], real_t maxabsNear[GRT_SRC_M_NUM]){
    GRT_LOOP_ChnlGrid(im, c){
        real_t amp = fabs(F[im][c]);
        if(amp > maxabsF[GRT_SRC_M_GTYPES[im]])  maxabsF[GRT_SRC_M_GTYPES[im]] = amp;
    }
    for(int im=0; im<GRT_SRC_M_NUM; ++im){
        if(GRT_SRC_M_ORDERS[im] == 0)  continue;
        real_t amp = fabs(F[im][0] + F[im][2]);
        if(amp > maxabsNear[im])  maxabsNear[im] = amp;
    }
}


/** 检查普通项和近场项的拟合情况，返回 true 表示通过 */
static bool check_fit(
    const KInterval *ptKitv, const KInterval *ptKitvL, const KInterval *ptKitvR, real_t kref,
    bool isuiz, bool calc_uir, const real_t maxabsQWV[GRT_GTYPES_MAX], const real_t maxabsNear[GRT_SRC_M_NUM], real_t tol)
{
    const cplxChnlGrid *F3L = (isuiz)? ptKitvL->Fz3 : ptKitvL->F3;
    const cplxChnlGrid *F3R = (isuiz)? ptKitvR->Fz3 : ptKitvR->F3;
    const KInterval *intervals[3] = {ptKitv, ptKitvL, ptKitvR};

    for(int im=0; im<GRT_SRC_M_NUM; ++im){
        for(int itype=0; itype<GRT_CHANNEL_NUM+1+(int)calc_uir; ++itype){
            bool isnear = itype >= GRT_CHANNEL_NUM;
            if(GRT_SRC_M_ORDERS[im] == 0 && (isnear || GRT_QWV_CODES[itype] == 'v'))  continue;
            if(isnear && maxabsNear[im] == 0.0)  continue;

            // 高波数段普通项仅检查 v，近场组合仍需独立检查
            if(!isnear && ptKitv->k3[0] > kref && GRT_QWV_CODES[itype] != 'v')  continue;

            real_t ref_amp = REF_AMP_SCALE * (isnear? maxabsNear[im] : maxabsQWV[GRT_SRC_M_GTYPES[im]]);
            bool islowamp = true;
            for(int d=0; d<3; ++d){
                cplx_t FL = (isnear)? F3L[d][im][0] + F3L[d][im][2] : F3L[d][im][itype];
                islowamp = islowamp && (fabs(FL) < ref_amp);
                if(d>0){
                    cplx_t FR = (isnear)? F3R[d][im][0] + F3R[d][im][2] : F3R[d][im][itype];
                    islowamp = islowamp && (fabs(FR) < ref_amp);
                }
            }

            cplx_t S11 = simpson(ptKitv,  im, itype, isuiz, SIMPSON_A_ApH);
            cplx_t S12 = simpson(ptKitv,  im, itype, isuiz, SIMPSON_ApH_Ap2H);
            cplx_t S21 = simpson(ptKitvL, im, itype, isuiz, SIMPSON_A_Ap2H);
            cplx_t S22 = simpson(ptKitvR, im, itype, isuiz, SIMPSON_A_Ap2H);
            cplx_t coarse[3] = {S11 + S12, S11, S12};
            cplx_t fine[3]   = {S21 + S22, S21, S22};

            // 分别比较整个区间、左半区间和右半区间的拟合积分
            for(int i=0; i<3; ++i){
                real_t S_ref = fabs(coarse[i] + fine[i]);
                if(islowamp){
                    real_t k1 = intervals[i]->k3[0], k2 = intervals[i]->k3[2];
                    real_t kcoef = (itype == GRT_CHANNEL_NUM)? 2.0*(sqrt(k2) - sqrt(k1))
                                                             : RTWOTHIRD*(k2*sqrt(k2) - k1*sqrt(k1));
                    S_ref = ref_amp * kcoef;
                }
                // 乘法比较允许两个积分同时为零，避免相对误差中的 0/0
                if(fabs(coarse[i] - fine[i]) > tol * S_ref)  return false;
            }
        }
    }
    return true;
}

/**
 * 根据该区间内采样的三个点，拟合普通项和近场项的慢变振幅，计算区间内的积分
 * 将Bessel函数近似为 sqrt(2/(\pi kr)) cos(kr - (2m+1)/4 \pi)，
 * 普通项拟合 sqrt(k)*F(k,w)，近场项拟合 (qm+vm)/sqrt(k)，这样积分时可以避免计算超越函数
 * 
 */
static void interv_integ(const KInterval *ptKitv, size_t nr, real_t *rs, K_INTEG *K)
{
    // 震中距rs循环
    for(size_t ir = 0; ir < nr; ++ir){

        // 防御：r=0 应在参数入口已拒绝（SAFIM 不适于零震中距）
        if(GRT_IS_ZERO(rs[ir]))  continue;


        memset(K->SUM, 0, sizeof(cplxIntegGrid));

        // 该分段内的积分
        grt_int_Pk_sa_filon(ptKitv->k3, rs[ir], ptKitv->F3, false, K->SUM);

        GRT_LOOP_IntegGrid(im, v){
            int modr = GRT_SRC_M_ORDERS[im];
            if((modr==0 && v!=0 && v!=2))  continue;

            K->sumJ[ir][im][v] += K->SUM[im][v];
        }

        if(K->calc_upar){
            //----------------------------- ui_z --------------------------------------
            grt_int_Pk_sa_filon(ptKitv->k3, rs[ir], ptKitv->Fz3, false, K->SUM);

            GRT_LOOP_IntegGrid(im, v){
                int modr = GRT_SRC_M_ORDERS[im];
                if((modr==0 && v!=0 && v!=2))  continue;

                K->sumJz[ir][im][v] += K->SUM[im][v];
            }

            //----------------------------- ui_r --------------------------------------
            grt_int_Pk_sa_filon(ptKitv->k3, rs[ir], ptKitv->F3, true, K->SUM);
            
            GRT_LOOP_IntegGrid(im, v){
                int modr = GRT_SRC_M_ORDERS[im];
                if((modr==0 && v!=0 && v!=2))  continue;

                K->sumJr[ir][im][v] += K->SUM[im][v];
            }

        }

    }
}



real_t grt_sa_filon_integ(
    MODEL1D_STATE *mstat, real_t k0, real_t dk0, real_t tol, real_t kmax, real_t kref,
    size_t nr, real_t *rs, K_INTEG *K, FILE *fstats, GRT_KernelFunc kerfunc)
{   
    
    real_t kmin = k0 + dk0;
    if(kmin >= kmax)  return k0;

    // 从0开始，存储第二部分Filon积分的结果
    K_INTEG *K2 = grt_init_K_INTEG(K->calc_upar, nr);

    // 区间栈
    KIntervalStack stack;
    stack_init(&stack, 64);

    // 初始化一个总区间
    KInterval Kitv = { 
        .k3 = {
            kmin, 
            (kmin+kmax)*0.5, 
            kmax
        }, 
        .F3 = {{{0}}}, 
        .Fz3 = {{{0}}},
        .F3_raw = {{{0}}},
        .Fz3_raw = {{{0}}}
    };
    for(int i=0; i<3; ++i) {
        kerfunc(mstat, Kitv.k3[i], Kitv.F3[i], K->calc_upar, Kitv.Fz3[i]);
        if(mstat->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

        memcpy(Kitv.F3_raw[i], Kitv.F3[i], sizeof(cplxChnlGrid));
        if(K->calc_upar){
            memcpy(Kitv.Fz3_raw[i], Kitv.Fz3[i], sizeof(cplxChnlGrid));
        }

        if(K->applyDCM){
            grt_int_dcm_revision(kmax, Kitv.k3[i], Kitv.F3[i], K->QWV_kmax, K->calc_upar, Kitv.Fz3[i], K->QWVz_kmax);
        }
    }
    stack_push(&stack, Kitv);

    // 以下采样过程中的遇到的QWV最大绝对值，注意其会动态变化
    // 其中VF,HF使用一套，EX,DD,DS,SS使用一套，分别对应格林函数Gij和格林函数空间导数Gij,k
    // 仅用于对积分幅度进行参考
    real_t maxabsQWV[GRT_GTYPES_MAX]={0};
    real_t maxabsQWV_uiz[GRT_GTYPES_MAX]={0};
    real_t maxabsNear[GRT_SRC_M_NUM]={0};
    real_t maxabsNear_uiz[GRT_SRC_M_NUM]={0};

    // 记录第一个值
    if(fstats!=NULL)  grt_write_stats(fstats, Kitv.k3[0], (K->calc_upar)? Kitv.Fz3_raw[0] : Kitv.F3_raw[0]);

    // 自适应采样
    while(stack.size > 0) {
        Kitv = stack_pop(&stack);

        // 左右两个区间
        KInterval Kitv_left = { 
            .k3 = {
                Kitv.k3[0], 
                (Kitv.k3[0]+Kitv.k3[1])*0.5, 
                Kitv.k3[1]
            }, 
            .F3 = {{{0}}}, 
            .Fz3 = {{{0}}},
            .F3_raw = {{{0}}},
            .Fz3_raw = {{{0}}}
        };
        KInterval Kitv_right = { 
            .k3 = {
                Kitv.k3[1], 
                (Kitv.k3[1]+Kitv.k3[2])*0.5, 
                Kitv.k3[2]
            }, 
            .F3 = {{{0}}}, 
            .Fz3 = {{{0}}},
            .F3_raw = {{{0}}},
            .Fz3_raw = {{{0}}}
        };
        memcpy(Kitv_left.F3[0], Kitv.F3[0], sizeof(cplxChnlGrid));
        memcpy(Kitv_left.F3[2], Kitv.F3[1], sizeof(cplxChnlGrid));
        memcpy(Kitv_right.F3[0], Kitv.F3[1], sizeof(cplxChnlGrid));
        memcpy(Kitv_right.F3[2], Kitv.F3[2], sizeof(cplxChnlGrid));
        memcpy(Kitv_left.F3_raw[0], Kitv.F3_raw[0], sizeof(cplxChnlGrid));
        memcpy(Kitv_left.F3_raw[2], Kitv.F3_raw[1], sizeof(cplxChnlGrid));
        memcpy(Kitv_right.F3_raw[0], Kitv.F3_raw[1], sizeof(cplxChnlGrid));
        memcpy(Kitv_right.F3_raw[2], Kitv.F3_raw[2], sizeof(cplxChnlGrid));
        if(K->calc_upar){
            memcpy(Kitv_left.Fz3[0], Kitv.Fz3[0], sizeof(cplxChnlGrid));
            memcpy(Kitv_left.Fz3[2], Kitv.Fz3[1], sizeof(cplxChnlGrid));
            memcpy(Kitv_right.Fz3[0], Kitv.Fz3[1], sizeof(cplxChnlGrid));
            memcpy(Kitv_right.Fz3[2], Kitv.Fz3[2], sizeof(cplxChnlGrid));
            memcpy(Kitv_left.Fz3_raw[0], Kitv.Fz3_raw[0], sizeof(cplxChnlGrid));
            memcpy(Kitv_left.Fz3_raw[2], Kitv.Fz3_raw[1], sizeof(cplxChnlGrid));
            memcpy(Kitv_right.Fz3_raw[0], Kitv.Fz3_raw[1], sizeof(cplxChnlGrid));
            memcpy(Kitv_right.Fz3_raw[2], Kitv.Fz3_raw[2], sizeof(cplxChnlGrid));
        }
        
        kerfunc(mstat, Kitv_left.k3[1], Kitv_left.F3[1], K->calc_upar, Kitv_left.Fz3[1]);
        if(mstat->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
        memcpy(Kitv_left.F3_raw[1], Kitv_left.F3[1], sizeof(cplxChnlGrid));
        if(K->calc_upar){
            memcpy(Kitv_left.Fz3_raw[1], Kitv_left.Fz3[1], sizeof(cplxChnlGrid));
        }

        kerfunc(mstat, Kitv_right.k3[1], Kitv_right.F3[1], K->calc_upar, Kitv_right.Fz3[1]);
        if(mstat->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
        memcpy(Kitv_right.F3_raw[1], Kitv_right.F3[1], sizeof(cplxChnlGrid));
        if(K->calc_upar){
            memcpy(Kitv_right.Fz3_raw[1], Kitv_right.Fz3[1], sizeof(cplxChnlGrid));
        }

        if(K->applyDCM){
            grt_int_dcm_revision(kmax, Kitv_left.k3[1], Kitv_left.F3[1], K->QWV_kmax, K->calc_upar, Kitv_left.Fz3[1], K->QWVz_kmax);
            grt_int_dcm_revision(kmax, Kitv_right.k3[1], Kitv_right.F3[1], K->QWV_kmax, K->calc_upar, Kitv_right.Fz3[1], K->QWVz_kmax);
        }

        // 增加新值，并比较QWV最大绝对值
        for(int i=0; i<3; ++i){
            update_maxabs(Kitv_left.F3[i], maxabsQWV, maxabsNear);
            if(K->calc_upar)  update_maxabs(Kitv_left.Fz3[i], maxabsQWV_uiz, maxabsNear_uiz);
            if(i>0){
                update_maxabs(Kitv_right.F3[i], maxabsQWV, maxabsNear);
                if(K->calc_upar)  update_maxabs(Kitv_right.Fz3[i], maxabsQWV_uiz, maxabsNear_uiz);
            }  
        }

        // 检查区间是否符合要求
        bool goodtol = check_fit(&Kitv, &Kitv_left, &Kitv_right, kref, false, K->calc_upar, maxabsQWV, maxabsNear, tol);
        if(K->calc_upar){
            goodtol = goodtol && check_fit(&Kitv, &Kitv_left, &Kitv_right, kref, true, false, maxabsQWV_uiz, maxabsNear_uiz, tol);
        }
        
        // 区间不符合要求，且采样间隔还足够继续细分
        if(! goodtol && Kitv_left.k3[2] - Kitv_left.k3[0] > SA_MIN_DK) {
            // 推入右子区间（后进先出）
            stack_push(&stack, Kitv_right);
            // 推入左子区间（先处理）
            stack_push(&stack, Kitv_left);
        } else {
            // 由于栈的特性，最终记录的k值采样是按顺序的
            // 记录后四个采样值
            if(fstats!=NULL){
                for(int i=1; i<3; ++i){
                    grt_write_stats(fstats, Kitv_left.k3[i], (K->calc_upar)? Kitv_left.Fz3_raw[i] : Kitv_left.F3_raw[i]);
                }
                for(int i=1; i<3; ++i){
                    grt_write_stats(fstats, Kitv_right.k3[i], (K->calc_upar)? Kitv_right.Fz3_raw[i] : Kitv_right.F3_raw[i]);
                }
            }
            // 计算积分
            interv_integ(&Kitv, nr, rs, K2);
        }
    } // END sampling

    // 乘上总系数 sqrt(2.0/(PI*r)) / dk0,  除dks0是在该函数外还会再乘dk0, 并将结果加到原数组中
    for(size_t ir = 0; ir < nr; ++ir){
        // 防御：r=0 应在参数入口已拒绝（总系数含 1/√r）
        if(GRT_IS_ZERO(rs[ir]))  continue;
        real_t tmp = sqrt(2.0/(PI*rs[ir])) / dk0;

        GRT_LOOP_IntegGrid(im, v){
            K->sumJ[ir][im][v] += K2->sumJ[ir][im][v] * tmp;
            if(K->calc_upar){
                K->sumJz[ir][im][v] += K2->sumJz[ir][im][v] * tmp;
                K->sumJr[ir][im][v] += K2->sumJr[ir][im][v] * tmp;
            }
        }
    }


    BEFORE_RETURN:
    GRT_SAFE_FREE_PTR(stack.data);
    grt_free_K_INTEG(K2);

    return Kitv.k3[2]; // 最后k值
}
