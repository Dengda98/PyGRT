/**
 * @file   saroot.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-05
 * 
 *     使用自适应采样寻找久期函数零点
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/common/search.h"

#include "grt/modal/secular.h"
#include "grt/modal/saroot.h"


/**
 * 自适应划分区间的最小dk
 */
#define SA_MIN_DK 1e-8


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
    real_t x3[3];
    cplx_t F3[3]; 
} Interval;

// 区间栈结构体
typedef struct { 
    Interval *data; 
    int size; 
    int capacity; 
} IntervalStack;

/** 初始化区间栈 */
static void stack_init(IntervalStack *stack, int init_capacity) {
    stack->data = (Interval*)malloc(init_capacity * sizeof(Interval));
    stack->size = 0;
    stack->capacity = init_capacity;
}

/** 从栈顶部压入元素 */
static void stack_push(IntervalStack *stack, Interval item) {
    // 扩容
    if(stack->size >= stack->capacity){
        stack->capacity *= 2;
        stack->data = (Interval*)realloc(stack->data, stack->capacity * sizeof(Interval));
    }
    stack->data[stack->size++] = item;
}

/** 从栈顶部拿走元素 */
static Interval stack_pop(IntervalStack *stack) {
    if(stack->size == 0) {
        fprintf(stderr, "Pop from empty stack\n");
        exit(EXIT_FAILURE);
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
static cplx_t simpson(const Interval *item_pt, SIMPSON_INTV stats) {
    cplx_t Fint = 0.0;
    real_t len = item_pt->x3[2] -  item_pt->x3[0];
    const cplx_t *F3 = item_pt->F3;

    if(stats == SIMPSON_A_Ap2H){
        Fint = len * (F3[0] + 4.0*F3[1] + F3[2]) / 6.0;
    }
    else if(stats == SIMPSON_As2H_A){
        Fint = len * (19.0*F3[0] - 20.0*F3[1] + 7.0*F3[2]) / 6.0;
    }
    else if(stats == SIMPSON_Ap2H_Ap4H){
        Fint = len * (7.0*F3[0] - 20.0*F3[1] + 19.0*F3[2]) / 6.0;
    }
    else if(stats == SIMPSON_A_ApH){
        Fint = len * (5.0*F3[0] + 8.0*F3[1] - F3[2]) / 24.0;
    }
    else if(stats == SIMPSON_ApH_Ap2H){
        Fint = len * ( - F3[0] + 8.0*F3[1] + 5.0*F3[2]) / 24.0;
    }
    else{
        fprintf(stderr, "wrong simpson stats (%d).\n", stats);
        exit(EXIT_FAILURE);
    }

    return Fint;
}


/** 检查区间是否符合要求，返回True表示通过 */
static bool check_fit(
    const Interval *ptCitv, const Interval *ptCitvL, const Interval *ptCitvR, real_t tol)
{
    // 计算积分差异
    cplx_t S11, S12, S21, S22;

    real_t S_dif, S_ref;
    bool badtol = false;

    S11 = simpson(ptCitv, SIMPSON_A_ApH);
    S12 = simpson(ptCitv, SIMPSON_ApH_Ap2H);
    S21 = simpson(ptCitvL, SIMPSON_A_Ap2H);
    S22 = simpson(ptCitvR, SIMPSON_A_Ap2H);

    // 三个拟合规则(a,b,c)
    // (a)
    S_ref = fabs(S11 + S12 + S21 + S22);
    S_dif = fabs(S11 + S12 - S21 - S22);
    badtol = badtol || (S_dif/S_ref > tol);  // 有一个不合格就继续采样
    if(badtol)  goto BEFORE_RETURN;
    // (b)
    S_ref = fabs(S11 + S21);
    S_dif = fabs(S11 - S21);
    badtol = badtol || (S_dif/S_ref > tol);  // 有一个不合格就继续采样
    if(badtol)  goto BEFORE_RETURN;
    // (c)
    S_ref = fabs(S12 + S22);
    S_dif = fabs(S12 - S22);
    badtol = badtol || (S_dif/S_ref > tol);  // 有一个不合格就继续采样
    if(badtol)  goto BEFORE_RETURN;


    BEFORE_RETURN:
    return (! badtol);
}



void grt_sa_secular_roots(GRT_MODEL1D *mod1d, EIGENV_METHOD *eigmet, EIGENV *eigv, const int secRaylType, const size_t iref)
{
    cplx_t root_sec = -1.0;

    // 区间栈
    IntervalStack stack;
    stack_init(&stack, 64);
    
    // 初始化一个总区间
    Interval Citv;
    Interval last_Citv_right = {.x3 = {-1.0}};
    
    // 在[cmin, cmax]之间初始化少量采样一些区间
    {
        real_t cmin=eigmet->cmin + eigmet->vgap;
        real_t cmax=eigmet->cmax - eigmet->vgap;
        size_t N = 61; // 要求是奇数
        real_t dc = (cmax - cmin)/(N - 1);
        for(size_t i = 0; i < N-2; i+=2){
            for(int p=0; p<3; ++p) {
                Citv.x3[2-p] = cmin + dc*(N-1-i-p);
                grt_secular_function(mod1d, Citv.x3[2-p], secRaylType, iref, eigmet->wtype, &Citv.F3[2-p]);
            }
            stack_push(&stack, Citv);
        }
    }

    // 记录第一个值
    if(eigmet->print_sec){
        fprintf(stdout, "# secRaylType=%d, iref=%zu, f=%f, tol=%e, c1=%f, c2=%f\n", secRaylType, iref, creal(mod1d->omega)/PI2, eigmet->satol, eigmet->cmin, eigmet->cmax);
        fprintf(stdout, "%.16e %.16e %.16e \n", Citv.x3[0], GRT_CMPLX_SPLIT(Citv.F3[0]));
    }

    cplx_t secCheck[6];
    real_t secAbsCheck[6];
    real_t cCheck[6];
    real_t dc, lastdc=9e30; 

    // 自适应采样
    // 注意这里是k值降序采样
    while(stack.size > 0) {
        // 阶数过多，提前退出
        // if(nroot >= maxnroot)  break;

        Citv = stack_pop(&stack);
        dc = Citv.x3[2] - Citv.x3[0];

        // 频散数组满载时，不必要在更高速（低波数）的区间搜索
        if(!eigmet->print_sec && eigv->n==eigmet->nmode && eigv->c_roots != NULL && Citv.x3[0] >= eigv->c_roots[eigv->n-1]) break;

        // 左右两个区间
        Interval Citv_left = { 
            .x3 = {
                Citv.x3[0], 
                (Citv.x3[0]+Citv.x3[1])*0.5, 
                Citv.x3[1]
            }, 
            .F3 = {0}, 
        };
        Interval Citv_right = { 
            .x3 = {
                Citv.x3[1], 
                (Citv.x3[1]+Citv.x3[2])*0.5, 
                Citv.x3[2]
            }, 
            .F3 = {0}, 
        };
        Citv_left.F3[0] = Citv.F3[0];
        Citv_left.F3[2] = Citv.F3[1];
        Citv_right.F3[0] = Citv.F3[1];
        Citv_right.F3[2] = Citv.F3[2];
        grt_secular_function(mod1d, Citv_left.x3[1], secRaylType, iref, eigmet->wtype, &Citv_left.F3[1]);
        grt_secular_function(mod1d, Citv_right.x3[1], secRaylType, iref, eigmet->wtype, &Citv_right.F3[1]);

        // 检查区间是否符合要求
        real_t goodtol = check_fit(&Citv, &Citv_left, &Citv_right, eigmet->satol);

        // 当且仅当 dk 在范围内时才判断 goodtol，否则对于较大dk一律划分子区间，较小dk一律结束区间
        if(dc > eigv->dcmax || (dc > eigv->dcmin && (!goodtol || dc > 2.0*lastdc)) ){
            // 推入右子区间（后进先出）
            stack_push(&stack, Citv_right);

            // 推入左子区间（先处理）
            stack_push(&stack, Citv_left);
        } else {
            
            lastdc = dc;

            // 由于栈的特性，最终记录的k值采样是按顺序的
            // 记录后四个采样值 
            if(eigmet->print_sec){
                for(int i=1; i<3; ++i){
                    fprintf(stdout, "%.16e %.16e %.16e \n", Citv_left.x3[i], GRT_CMPLX_SPLIT(Citv_left.F3[i]));
                }
                for(int i=1; i<3; ++i){
                    fprintf(stdout, "%.16e %.16e %.16e \n", Citv_right.x3[i], GRT_CMPLX_SPLIT(Citv_right.F3[i]));
                }
                continue;
            }
            
            // 找零点
            secCheck[0] = last_Citv_right.F3[1];
            secAbsCheck[0] = fabs(secCheck[0]);
            cCheck[0] = last_Citv_right.x3[1];
            real_t root_c = -1.0;
            for(int i=0; i<3; ++i){
                secCheck[i+1] = Citv_left.F3[i];
                secAbsCheck[i+1] = fabs(secCheck[i+1]);
                cCheck[i+1] = Citv_left.x3[i];
            }
            for(int i=1; i<3; ++i){
                secCheck[i+3] = Citv_right.F3[i];
                secAbsCheck[i+3] = fabs(secCheck[i+3]);
                cCheck[i+3] = Citv_right.x3[i];
            }
            
            // 检查零点，久期函数绝对值是否形成波谷
            // 假设6点中仅含有一个波谷
            bool hasroot = false;
            for(int i=2; i<6; ++i){
                if(cCheck[i] < 0.0)  continue;  // 未初始化的右边界点
                if(cCheck[i] < root_c)  continue;  // 该零点已搜索过

                hasroot = secAbsCheck[i-2] >= secAbsCheck[i-1] && secAbsCheck[i-1] <= secAbsCheck[i];
                if(!hasroot)  continue;

                // real_t dk0 = 4.0*(kCheck[5] - kCheck[4]);
                // if(dk0 > dk)  dk0 = dk;
                // MYINT is = bisearch_root(
                //     mod1d, w, kCheck[i-1], dk0*0.5, secCheck[i-1], secAbsCheck[i-1], 
                //     secRaylType, iref, isRayl, &root_sec, &root_k
                // );
                grt_goldensection_search_root(
                    mod1d, cCheck[i-2], cCheck[i], 
                    secRaylType, iref, eigmet->wtype, eigmet->rtol, &root_sec, &root_c
                );

                if(root_c > 0.0 && !grt_check_vel_in_mod(mod1d, root_c, eigmet->vgap)){
                    if(eigv->c_roots == NULL){
                        // 需初始化动态内存
                        eigv->c_roots = (real_t *)realloc(eigv->c_roots, sizeof(real_t)*1);
                        eigv->c_roots_iref = (size_t *)realloc(eigv->c_roots_iref, sizeof(size_t)*1);
                        eigv->c_roots[0] = root_c;
                        eigv->c_roots_iref[0] = iref;
                        eigv->n = 1;
                    } else {
                        // 排除很接近的值
                        size_t idx = grt_findClosest_real_t(eigv->c_roots, eigv->n, root_c);
                        if(fabs(root_c - eigv->c_roots[idx]) > eigmet->cgap){
                            size_t capacity = eigmet->nmode;
                            // 扩容动态内存
                            if(capacity <= 0){
                                capacity = eigv->n+1;
                            }
                            eigv->c_roots = (real_t *)realloc(eigv->c_roots, sizeof(real_t)*capacity);
                            eigv->c_roots_iref = (size_t *)realloc(eigv->c_roots_iref, sizeof(size_t)*capacity);

                            ssize_t pos = grt_insertOrdered(eigv->c_roots, &eigv->n, capacity, &root_c, sizeof(real_t), true, grt_compare_real_t);

                            // 更新对应 c_roots_iref
                            if(pos >= 0){
                                memmove(eigv->c_roots_iref+(pos+1), eigv->c_roots_iref+pos, (eigv->n-1 - pos)*sizeof(size_t));
                                eigv->c_roots_iref[pos] = iref;
                            }
                        } 
                    }
                }
            }
            
            last_Citv_right = Citv_right;
            
        } // END if
    } // END while

    GRT_SAFE_FREE_PTR(stack.data);
}



// void grt_sa_secular_roots(
//     const GRT_MODEL1D *mod1d, const real_t *V_sort, const real_t w, 
//     const real_t c10, const real_t tol, const real_t c20, const real_t dcmin, const real_t dcmax,
//     const bool isRayl, const MYINT secRaylType, const MYINT iref, const bool print_sec,
//     const MYINT nmode, real_t **pt_k_roots, MYINT **pt_k_roots_iref, MYINT *Nroot)
// {
//     MYINT inv_stats=GRT_INVERSE_SUCCESS;
//     MYINT nroot = *Nroot;
//     cplx_t root_sec = -1.0;

//     // 区间栈
//     KIntervalStack stack;
//     stack_init(&stack, 64);

    
//     // 初始化一个总区间
//     KInterval Kitv;
//     KInterval last_Kitv_left = {.k3 = {-1.0}};
    
//     // 根据速度分布，初始化一些区间
//     real_t kmin, kmax;
//     MYINT uniq_n = 0;
//     kmin = w/c20 + GRT_MIN_K_GAP;
//     // 速度降序循环，这样初始化的区间k在栈中会降序排列，但每个区间内还是以升序排列
//     for(MYINT u=mod1d->n*2-1; u >= 0; --u){
//         // 跳过重复速度
//         if( uniq_n > 0 && V_sort[u] == V_sort[mod1d->n*2 - uniq_n])  continue;
//         uniq_n++;
        
//         // 只在区间内搜索
//         // 由于是速度降序循环，一定是先碰到非零速度，再到可能存在的 0 速度
//         // 而传入的 c10, c20 必须非 0 ，故会成功的触发 break
//         if(V_sort[u] >= c20)  continue;
//         if(V_sort[u] < c10)  break;
        
//         // 判断使 c10 始终在范围内
//         if(u==0 || V_sort[u-1] < c10){
//             kmax = w/c10 + GRT_MIN_K_GAP;
//         } else {
//             kmax = w/V_sort[u] + GRT_MIN_K_GAP;
//         }

//         KInterval Kitv1 = { 
//             .k3 = {
//                 kmin, 
//                 (kmin+kmax)*0.5, 
//                 kmax
//             },
//             .F3 = {0}, 
//         };
//         Kitv = Kitv1;
        
//         for(MYINT i=0; i<3; ++i) {
//             grt_secular_function(mod1d, w, Kitv.k3[i], secRaylType, iref, isRayl, &Kitv.F3[i], &inv_stats);
//         }
//         stack_push(&stack, Kitv);

//         kmin = kmax;
//     }

//     // 确定 dc/c 的范围
//     real_t rdcmin = dcmin / c20;
//     real_t rdcmax = dcmax / c20;

//     // 记录第一个值
//     if(print_sec){
//         fprintf(stdout, "# secRaylType=%d, iref=%d, f=%f, tol=%e, c1=%f, c2=%f\n", secRaylType, iref, w/PI2, tol, c10, c20);
//         fprintf(stdout, "%.16e %.16e %.16e \n", w/Kitv.k3[2], creal(Kitv.F3[2]), cimag(Kitv.F3[2]));
//     }

//     cplx_t secCheck[6];
//     real_t secAbsCheck[6];
//     real_t kCheck[6];
//     // 4096相当于初始化一个比较细的dk
//     real_t dk, lastdk = (w/c10 - w/c20)/4096; 

//     // 自适应采样
//     // 注意这里是k值降序采样
//     while(stack.size > 0) {
//         // 阶数过多，提前退出
//         // if(nroot >= maxnroot)  break;

//         Kitv = stack_pop(&stack);
//         dk = Kitv.k3[2] - Kitv.k3[0];

//         // 频散数组满载时，不必要在更高速（低波数）的区间搜索
//         if(!print_sec && nroot==nmode && *pt_k_roots != NULL && Kitv.k3[2] <= (*pt_k_roots)[nroot-1]) break;

//         // 左右两个区间
//         KInterval Kitv_left = { 
//             .k3 = {
//                 Kitv.k3[0], 
//                 (Kitv.k3[0]+Kitv.k3[1])*0.5, 
//                 Kitv.k3[1]
//             }, 
//             .F3 = {0}, 
//         };
//         KInterval Kitv_right = { 
//             .k3 = {
//                 Kitv.k3[1], 
//                 (Kitv.k3[1]+Kitv.k3[2])*0.5, 
//                 Kitv.k3[2]
//             }, 
//             .F3 = {0}, 
//         };
//         Kitv_left.F3[0] = Kitv.F3[0];
//         Kitv_left.F3[2] = Kitv.F3[1];
//         Kitv_right.F3[0] = Kitv.F3[1];
//         Kitv_right.F3[2] = Kitv.F3[2];
//         grt_secular_function(mod1d, w, Kitv_left.k3[1], secRaylType, iref, isRayl, &Kitv_left.F3[1], &inv_stats);
//         grt_secular_function(mod1d, w, Kitv_right.k3[1], secRaylType, iref, isRayl, &Kitv_right.F3[1], &inv_stats);

//         // 检查区间是否符合要求
//         real_t goodtol = check_fit(&Kitv, &Kitv_left, &Kitv_right, tol);

//         // 区间不符合要求或区间相比上一个区间过宽，且采样间隔还足够继续细分
//         // if( (!goodtol || dk > 2.0*lastdk) && dk > 2.0*GRT_MIN_K_GAP) {

//         // 当且仅当 dk 在范围内时才判断 goodtol，否则对于较大dk一律划分子区间，较小dk一律结束区间
//         if(dk/Kitv.k3[2] > rdcmax || (dk/Kitv.k3[2] > rdcmin && (!goodtol || dk > 2.0*lastdk)) ){
//             // 推入左子区间（后进先出）
//             stack_push(&stack, Kitv_left);

//             // 推入右子区间（先处理）
//             stack_push(&stack, Kitv_right);
//         } else {
            
//             lastdk = dk;

//             // 由于栈的特性，最终记录的k值采样是按顺序的
//             // 记录后四个采样值 
//             if(print_sec){
//                 for(MYINT i=1; i>=0; --i){
//                     fprintf(stdout, "%.16e %.16e %.16e \n", w/Kitv_right.k3[i], creal(Kitv_right.F3[i]), cimag(Kitv_right.F3[i]));
//                 }
//                 for(MYINT i=2; i>=0; --i){
//                     fprintf(stdout, "%.16e %.16e %.16e \n", w/Kitv_left.k3[i], creal(Kitv_left.F3[i]), cimag(Kitv_left.F3[i]));
//                 }
//                 continue;
//             }
            
//             // 找零点
//             real_t root_k = -1.0;
//             for(MYINT i=0; i<3; ++i){
//                 secCheck[i] = Kitv_left.F3[i];
//                 secAbsCheck[i] = fabs(secCheck[i]);
//                 kCheck[i] = Kitv_left.k3[i];
//             }
//             for(MYINT i=1; i<3; ++i){
//                 secCheck[i+2] = Kitv_right.F3[i];
//                 secAbsCheck[i+2] = fabs(secCheck[i+2]);
//                 kCheck[i+2] = Kitv_right.k3[i];
//             }
//             secCheck[5] = last_Kitv_left.F3[1];
//             secAbsCheck[5] = fabs(secCheck[5]);
//             kCheck[5] = last_Kitv_left.k3[1];
            
//             // 检查零点，久期函数绝对值是否形成波谷
//             // 假设6点中仅含有一个波谷
//             // fprintf(stderr, "# kCheck[6]\n");
//             // for(int i=0; i<6; ++i){fprintf(stderr, "# %.16e\n", kCheck[i]);}
//             bool hasroot = false;
//             for(MYINT i=2; i<6; ++i){
//                 if(kCheck[i] < 0.0)  continue;  // 未初始化的右边界点
//                 if(kCheck[i] < root_k)  continue;  // 该零点已搜索过

//                 hasroot = secAbsCheck[i-2] >= secAbsCheck[i-1] && secAbsCheck[i-1] <= secAbsCheck[i];
//                 if(!hasroot)  continue;

//                 // real_t dk0 = 4.0*(kCheck[5] - kCheck[4]);
//                 // if(dk0 > dk)  dk0 = dk;
//                 // MYINT is = bisearch_root(
//                 //     mod1d, w, kCheck[i-1], dk0*0.5, secCheck[i-1], secAbsCheck[i-1], 
//                 //     secRaylType, iref, isRayl, &root_sec, &root_k
//                 // );
//                 grt_goldensection_search_root(
//                     mod1d, w, kCheck[i-2], kCheck[i],
//                     secRaylType, iref, isRayl, &root_sec, &root_k
//                 );

//                 if(root_k > 0.0 && !grt_check_vel_in_mod(mod1d, w/root_k, GRT_MIN_MODEL_V_GAP)){
//                     if(*pt_k_roots == NULL){
//                         // 需初始化动态内存
//                         *pt_k_roots = (real_t *)realloc(*pt_k_roots, sizeof(real_t)*1);
//                         *pt_k_roots_iref = (MYINT *)realloc(*pt_k_roots_iref, sizeof(MYINT)*1);
//                         (*pt_k_roots)[0] = root_k;
//                         (*pt_k_roots_iref)[0] = iref;
//                         nroot = 1;
//                     } else {
//                         // 排除很接近的值
//                         MYINT idx = grt_findClosest_real_t(*pt_k_roots, nroot, root_k);
//                         if(fabs(root_k - (*pt_k_roots)[idx]) > GRT_MIN_K_GAP){
//                             MYINT capacity = nmode;
//                             // 扩容动态内存
//                             if(capacity <= 0){
//                                 capacity = nroot+1;
//                             }
//                             *pt_k_roots = (real_t *)realloc(*pt_k_roots, sizeof(real_t)*capacity);
//                             *pt_k_roots_iref = (MYINT *)realloc(*pt_k_roots_iref, sizeof(MYINT)*capacity);

//                             MYINT pos = grt_insertOrdered(*pt_k_roots, &nroot, capacity, &root_k, sizeof(real_t), false, grt_compare_real_t);

//                             // 更新对应 k_roots_iref
//                             if(pos >= 0){
//                                 memmove(*pt_k_roots_iref+(pos+1), *pt_k_roots_iref+pos, (nroot-1 - pos)*sizeof(MYINT));
//                                 (*pt_k_roots_iref)[pos] = iref;
//                             }
//                         } 
//                     }
//                 }
//             }
            
//             last_Kitv_left = Kitv_left;
            
//         } // END if
//     } // END while

//     GRT_SAFE_FREE_PTR(stack.data);

//     *Nroot = nroot;
// }


