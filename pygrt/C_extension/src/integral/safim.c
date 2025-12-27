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
} KInterval;

// 区间栈结构体
typedef struct { 
    KInterval *data; 
    size_t size; 
    size_t capacity; 
} KIntervalStack;

/** 初始化区间栈 */
static void stack_init(KIntervalStack *stack, size_t init_capacity) {
    stack->data = (KInterval*)malloc(init_capacity * sizeof(KInterval));
    stack->size = 0;
    stack->capacity = init_capacity;
}

/** 从栈顶部压入元素 */
static void stack_push(KIntervalStack *stack, KInterval item) {
    // 扩容
    if(stack->size >= stack->capacity){
        stack->capacity *= 2;
        stack->data = (KInterval*)realloc(stack->data, stack->capacity * sizeof(KInterval));
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
 * 假设区间内三点为[a, a+h, a+2h]，stats的取值对应的积分区间为
 *    + stats = -1,  [a-2h, a]
 *    + stats = 1,   [a, a+2h]
 *    + stats = 2,   [a+2h, a+4h]
 *    + stats = 3,   [a, a+h]
 *    + stats = 4,   [a+h, a+2h]
 * 
 */
static cplx_t simpson(const KInterval *item_pt, int im, int iqwv, bool isuiz, SIMPSON_INTV stats) {
    cplx_t Fint = 0.0;
    real_t klen = item_pt->k3[2] -  item_pt->k3[0];
    const cplxChnlGrid *F3 = (isuiz)? item_pt->Fz3 : item_pt->F3;
    
    // 使用F(k)*sqrt(k)来衡量积分值，这可以平衡后续计算F(k)*Jm(kr)k积分时的系数
    real_t sk[3];
    for(int i=0; i<3; ++i){
        sk[i] = sqrt(item_pt->k3[i]);
    }

    if(stats == SIMPSON_A_Ap2H){
        Fint = klen * (F3[0][im][iqwv]*sk[0] + 4.0*F3[1][im][iqwv]*sk[1] + F3[2][im][iqwv]*sk[2]) / 6.0;
    }
    else if(stats == SIMPSON_As2H_A){
        Fint = klen * (19.0*F3[0][im][iqwv]*sk[0] - 20.0*F3[1][im][iqwv]*sk[1] + 7.0*F3[2][im][iqwv]*sk[2]) / 6.0;
    }
    else if(stats == SIMPSON_Ap2H_Ap4H){
        Fint = klen * (7.0*F3[0][im][iqwv]*sk[0] - 20.0*F3[1][im][iqwv]*sk[1] + 19.0*F3[2][im][iqwv]*sk[2]) / 6.0;
    }
    else if(stats == SIMPSON_A_ApH){
        Fint = klen * (5.0*F3[0][im][iqwv]*sk[0] + 8.0*F3[1][im][iqwv]*sk[1] - F3[2][im][iqwv]*sk[2]) / 24.0;
    }
    else if(stats == SIMPSON_ApH_Ap2H){
        Fint = klen * ( - F3[0][im][iqwv]*sk[0] + 8.0*F3[1][im][iqwv]*sk[1] + 5.0*F3[2][im][iqwv]*sk[2]) / 24.0;
    }
    else{
        GRTRaiseError("wrong simpson stats (%d).\n", stats);
    }

    return Fint;
}

/** 比较QWV的最大绝对值 */
static void get_maxabsQWV(const cplxChnlGrid F, real_t maxabsF[GRT_GTYPES_MAX]){
    real_t tmp;
    GRT_LOOP_ChnlGrid(im, c){
        tmp = fabs(F[im][c]);
        if(tmp > maxabsF[GRT_SRC_M_GTYPES[im]]){
            maxabsF[GRT_SRC_M_GTYPES[im]] = tmp;
        }
    }
}


/** 检查区间是否符合要求，返回True表示通过 */
static bool check_fit(
    const KInterval *ptKitv, const KInterval *ptKitvL, const KInterval *ptKitvR, real_t kref,
    bool isuiz, real_t maxabsQWV[GRT_GTYPES_MAX], real_t tol)
{
    // 计算积分差异
    cplx_t S11, S12, S21, S22;

    // 核函数
    const cplxChnlGrid *F3L = (isuiz)? ptKitvL->Fz3 : ptKitvL->F3;
    const cplxChnlGrid *F3R = (isuiz)? ptKitvR->Fz3 : ptKitvR->F3;

    // 取近似积分 \int_k1^k2 k^0.5 dk
    real_t kcoef13 = RTWOTHIRD*( ptKitv->k3[2]*sqrt(ptKitv->k3[2]) - ptKitv->k3[0]*sqrt(ptKitv->k3[0]) );
    real_t kcoef12 = RTWOTHIRD*( ptKitvL->k3[2]*sqrt(ptKitvL->k3[2]) - ptKitvL->k3[0]*sqrt(ptKitvL->k3[0]) );
    real_t kcoef23 = RTWOTHIRD*( ptKitvR->k3[2]*sqrt(ptKitvR->k3[2]) - ptKitvR->k3[0]*sqrt(ptKitvR->k3[0]) );

    real_t S_dif, S_ref;
    bool badtol = false;

    GRT_LOOP_ChnlGrid(im, c){
        int igtyp = GRT_SRC_M_GTYPES[im];
        if(GRT_SRC_M_ORDERS[im]==0 && GRT_QWV_CODES[c]=='v')  continue;
        // qw和v分开采样?
        // if(isqw && GRT_QWV_CODES[c]=='v')  continue;
        // if(!isqw && GRT_QWV_CODES[c]!='v') continue;

        // k值在后段，只根据数值稳定的v分量判断是否拟合好
        if(ptKitv->k3[0] > kref && GRT_QWV_CODES[c]!='v')  continue;

        S11 = simpson(ptKitv, im, c, isuiz, SIMPSON_A_ApH);
        S12 = simpson(ptKitv, im, c, isuiz, SIMPSON_ApH_Ap2H);
        S21 = simpson(ptKitvL, im, c, isuiz, SIMPSON_A_Ap2H);
        S22 = simpson(ptKitvR, im, c, isuiz, SIMPSON_A_Ap2H);

        bool islowamp = true;
        real_t ref_amp = REF_AMP_SCALE*maxabsQWV[igtyp];
        // 比较当前区间内5个核函数幅值，是否都低于参考值
        for(int d=0; d<3; ++d){
            islowamp = islowamp && (fabs(F3L[d][im][c]) < ref_amp);
            if(d>0){
                islowamp = islowamp && (fabs(F3R[d][im][c]) < ref_amp);
            }
        }
        // 如果核函数幅值确实较低，则限制S_ref，此时正负号不重要
        // 三个拟合规则(a,b,c)
        // (a)
        if(islowamp){
            S_ref = ref_amp * kcoef13;
        } else {
            S_ref = fabs(S11 + S12 + S21 + S22);
        }
        S_dif = fabs(S11 + S12 - S21 - S22);
        badtol = badtol || (S_dif/S_ref > tol);  // 有一个不合格就继续采样
        if(badtol)  goto BEFORE_RETURN;
        // (b)
        if(islowamp){
            S_ref = ref_amp * kcoef12;
        } else {
            S_ref = fabs(S11 + S21);
        }
        S_dif = fabs(S11 - S21);
        badtol = badtol || (S_dif/S_ref > tol);  // 有一个不合格就继续采样
        if(badtol)  goto BEFORE_RETURN;
        // (c)
        if(islowamp){
            S_ref = ref_amp * kcoef23;
        } else {
            S_ref = fabs(S12 + S22);
        }
        S_dif = fabs(S12 - S22);
        badtol = badtol || (S_dif/S_ref > tol);  // 有一个不合格就继续采样
        if(badtol)  goto BEFORE_RETURN;
    }

    BEFORE_RETURN:
    return (! badtol);
}

/**
 * 根据该区间内采样的三个点，拟合二次函数，计算 F(k,w)Jm(kr)k 在区间内的积分
 * 将Bessel函数近似为 sqrt(2/(\pi kr)) cos(kr - (2m+1)/4 \pi)，
 * 以下实际拟合的二次函数是 sqrt(k)*F(k,w), 这样积分时可以避免计算超越函数
 * 
 */
static void interv_integ(const KInterval *ptKitv, size_t nr, real_t *rs, K_INTEG *K)
{
    // 震中距rs循环
    for(size_t ir=0; ir<nr; ++ir){

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
    GRT_MODEL1D *mod1d, real_t k0, real_t dk0, real_t tol, real_t kmax, real_t kref,
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
        .Fz3 = {{{0}}} 
    };
    for(int i=0; i<3; ++i) {
        kerfunc(mod1d, Kitv.k3[i], Kitv.F3[i], K->calc_upar, Kitv.Fz3[i]);
        if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;
        
        if(K->applyDCM){
            GRT_LOOP_ChnlGrid(im, c){
                Kitv.F3[i][im][c] -= K->QWV_kmax[im][c];
                if(K->calc_upar) Kitv.Fz3[i][im][c] -= K->QWVz_kmax[im][c];
            }
        }
    }
    stack_push(&stack, Kitv);

    // 以下采样过程中的遇到的QWV最大绝对值，注意其会动态变化
    // 其中VF,HF使用一套，EX,DD,DS,SS使用一套，分别对应格林函数Gij和格林函数空间导数Gij,k
    // 仅用于对积分幅度进行参考
    real_t maxabsQWV[GRT_GTYPES_MAX]={0};
    real_t maxabsQWV_uiz[GRT_GTYPES_MAX]={0};

    // 记录第一个值
    if(fstats!=NULL)  grt_write_stats(fstats, Kitv.k3[0], (K->calc_upar)? Kitv.Fz3[0] : Kitv.F3[0]);

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
            .Fz3 = {{{0}}} 
        };
        KInterval Kitv_right = { 
            .k3 = {
                Kitv.k3[1], 
                (Kitv.k3[1]+Kitv.k3[2])*0.5, 
                Kitv.k3[2]
            }, 
            .F3 = {{{0}}}, 
            .Fz3 = {{{0}}} 
        };
        memcpy(Kitv_left.F3[0], Kitv.F3[0], sizeof(cplxChnlGrid));
        memcpy(Kitv_left.F3[2], Kitv.F3[1], sizeof(cplxChnlGrid));
        memcpy(Kitv_right.F3[0], Kitv.F3[1], sizeof(cplxChnlGrid));
        memcpy(Kitv_right.F3[2], Kitv.F3[2], sizeof(cplxChnlGrid));
        if(K->calc_upar){
            memcpy(Kitv_left.Fz3[0], Kitv.Fz3[0], sizeof(cplxChnlGrid));
            memcpy(Kitv_left.Fz3[2], Kitv.Fz3[1], sizeof(cplxChnlGrid));
            memcpy(Kitv_right.Fz3[0], Kitv.Fz3[1], sizeof(cplxChnlGrid));
            memcpy(Kitv_right.Fz3[2], Kitv.Fz3[2], sizeof(cplxChnlGrid));
        }
        
        kerfunc(mod1d, Kitv_left.k3[1], Kitv_left.F3[1], K->calc_upar, Kitv_left.Fz3[1]);
        if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

        kerfunc(mod1d, Kitv_right.k3[1], Kitv_right.F3[1], K->calc_upar, Kitv_right.Fz3[1]);
        if(mod1d->stats==GRT_INVERSE_FAILURE)  goto BEFORE_RETURN;

        if(K->applyDCM){
            GRT_LOOP_ChnlGrid(im, c){
                Kitv_left.F3[1][im][c] -= K->QWV_kmax[im][c];
                Kitv_right.F3[1][im][c] -= K->QWV_kmax[im][c];
                if(K->calc_upar){
                    Kitv_left.Fz3[1][im][c] -= K->QWVz_kmax[im][c];
                    Kitv_right.Fz3[1][im][c] -= K->QWVz_kmax[im][c];
                }
            }
        }

        // 增加新值，并比较QWV最大绝对值
        for(int i=0; i<3; ++i){
            get_maxabsQWV(Kitv_left.F3[i], maxabsQWV);
            if(K->calc_upar)  get_maxabsQWV(Kitv_left.Fz3[i], maxabsQWV_uiz);
            if(i>0){
                get_maxabsQWV(Kitv_right.F3[i], maxabsQWV);
                if(K->calc_upar)  get_maxabsQWV(Kitv_right.Fz3[i], maxabsQWV_uiz);
            }  
        }

        // 检查区间是否符合要求
        real_t goodtol = check_fit(&Kitv, &Kitv_left, &Kitv_right, kref, false, maxabsQWV, tol);
        if(K->calc_upar){
            goodtol = goodtol && check_fit(&Kitv, &Kitv_left, &Kitv_right, kref, true, maxabsQWV_uiz, tol);
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
                    grt_write_stats(fstats, Kitv_left.k3[i], (K->calc_upar)? Kitv_left.Fz3[i] : Kitv_left.F3[i]);
                }
                for(int i=1; i<3; ++i){
                    grt_write_stats(fstats, Kitv_right.k3[i], (K->calc_upar)? Kitv_right.Fz3[i] : Kitv_right.F3[i]);
                }
            }
            // 计算积分
            interv_integ(&Kitv, nr, rs, K2);
        }
    } // END sampling

    // 乘上总系数 sqrt(2.0/(PI*r)) / dk0,  除dks0是在该函数外还会再乘dk0, 并将结果加到原数组中
    for(size_t ir=0; ir<nr; ++ir){
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
