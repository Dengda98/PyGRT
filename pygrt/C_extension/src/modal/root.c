/**
 * @file   root.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-05
 * 
 *     搜寻久期函数零点
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/common/search.h"
#include "grt/common/progressbar.h"

#include "grt/modal/root.h"
#include "grt/modal/secular.h"



/** 
 * 黄金分割法确定久期函数零点 
 * 
 * @param[in]      mod1d        模型结构体指针
 * @param[in]      w            圆频率
 * @param[in]      c1           左区间相速度
 * @param[in]      c2           右区间相速度
 * @param[in]      iref         久期函数层位
 * @param[in]      wtype        频散类型
 * @param[in]      rtol         久期函数零点附近的幅值阈值
 * @param[out]     root_sec     零点的久期函数值
 * @param[out]     root_c       零点的相速度
 * 
 */
static size_t grt_goldensection_search_root(
    GRT_MODEL1D *mod1d, const real_t c1, const real_t c2, 
    const size_t iref, const DISPER_TYPE wtype,
    const real_t rtol, const real_t cgap, cplx_t *root_sec, real_t *root_c)
{
    // 细化阶段使用高精度矩阵计算
    GRT_USE_HIGH_PRECISION = true;

    real_t c[2] = {c1, c2};
    real_t x[2];
    x[0] = c[1] - GOLDEN_RATIO*(c[1] - c[0]);
    x[1] = c[0] + GOLDEN_RATIO*(c[1] - c[0]);
    cplx_t sec[2];
    real_t secABS[2];

    for(int i=0; i<2; ++i){
        grt_secular_function(mod1d, x[i], iref, wtype, &sec[i]);
        secABS[i] = fabs(sec[i]);
    }
    size_t is=0;

    #if __DEBUG_SECULAR__ == 1
    fprintf(stderr, "# f=%f, root_sec2=%.9e%+.9ej|%.9e%+.9ej, c1=%.16e, c2=%.16e, dc=%.16e, root_c=%.9e, %d\n", 
            creal(mod1d->omega)/PI2, creal(sec[0]), cimag(sec[0]), creal(sec[1]), cimag(sec[1]), x[0], x[1], x[1]-x[0], 
            x[0], is);
    fflush(stderr);
    #endif
    
    do{
        if(secABS[1] > secABS[0]){
            c[1] = x[1];
            sec[1] = sec[0];
            secABS[1] = secABS[0];
            x[1] = x[0];
            x[0] = c[1] - GOLDEN_RATIO*(c[1] - c[0]);
            grt_secular_function(mod1d, x[0], iref, wtype, &sec[0]);
            secABS[0] = fabs(sec[0]);
        } else {
            c[0] = x[0];
            sec[0] = sec[1];
            secABS[0] = secABS[1];
            x[0] = x[1];
            x[1] = c[0] + GOLDEN_RATIO*(c[1] - c[0]);
            grt_secular_function(mod1d, x[1], iref, wtype, &sec[1]);
            secABS[1] = fabs(sec[1]);
        }

        if(secABS[0] < rtol && c[1] - c[0] < 0.5*cgap)  break;

        // 接近达到浮点数精度极限
        if(c[1] - c[0] < 1e-15)  break;

        ++is;
        #if __DEBUG_SECULAR__ == 1
        fprintf(stderr, "# f=%f, root_sec2=%.9e%+.9ej|%.9e%+.9ej, secABS=%.9e, c1=%.16e, c2=%.16e, dc=%.19e, root_c=%.9e, %d\n", 
            creal(mod1d->omega)/PI2, creal(sec[0]), cimag(sec[0]), creal(sec[1]), cimag(sec[1]), secABS[0], x[0], x[1], x[1]-x[0], 
            x[0], is);
        fflush(stderr);
        #endif
        
    } while(is < 800);

    #if __DEBUG_SECULAR__ == 1
    fprintf(stderr, "# f=%f, root_sec2=%.9e%+.9ej|%.9e%+.9ej, dc=%e, root_c=%.9e, %d\n", 
            creal(mod1d->omega)/PI2, creal(sec[0]), cimag(sec[0]), creal(sec[1]), cimag(sec[1]), x[1]-x[0], 
            x[0], is);
    fflush(stderr);
    #endif

    // 是否为零点
    if(secABS[0] < rtol){
        *root_sec = sec[0];
        *root_c = x[0];
        #if __DEBUG_SECULAR__ == 1
        fprintf(stderr, "#      f=%f, root_sec=%.9e%+.9ej, root_c=%.9e, %d\n", creal(mod1d->omega)/PI2, creal(sec[0]), cimag(sec[0]), x[0], is);
        fflush(stderr);
        #endif
    } else {
        *root_sec = -999.0;
        *root_c = -1.0;
    }

    // 恢复
    GRT_USE_HIGH_PRECISION = false;

    return is;
}


// ---------------------------------------------------------------------------------------------------
//                                   自适应采样寻找久期函数零点
// ---------------------------------------------------------------------------------------------------

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


/**
 * 使用自适应采样寻找久期函数零点
 * 
 * @param[in]      mod1d        模型结构体指针
 * @param[in]      V_sort       升序的模型速度
 * @param[in]      w            圆频率
 * @param[in]      c10          相速度搜索下界
 * @param[in]      tol          自适应收敛精度
 * @param[in]      c20          相速度搜索上界
 * @param[in]      isRayl       Rayleigh or Love
 * @param[in]      iref         久期函数层位
 * @param[in]      print_sec    是否打印久期函数
 * @param[in]      nmode        最大零点个数, 如果 nmode<=0 则找全部零点
 * @param[out]     pt_c_roots   存储零点的本征值
 * @param[out]     pt_c_roots_iref   存储所使用的久期函数层位
 * @param[out]     Nroot        最终零点数量
 * 
 */
static void grt_adaptive_step_secular_roots(
    GRT_MODEL1D *mod1d, EIGENV_INFO *eigmet, EIGENV *eigv, const size_t iref,
    const real_t *cpred, const size_t npred)
{
    cplx_t root_sec = -1.0;

    // 区间栈
    IntervalStack stack;
    stack_init(&stack, 64);
    
    // 初始化一个总区间
    Interval Citv;
    Interval last_Citv_right = {.x3 = {-1.0}};
    
    // 根据 cpred 在[cmin, cmax]之间初始化采样一些区间
    for(size_t i = 0; i < npred - 1; ++i){
        real_t c1 = cpred[npred - i - 2];
        real_t c2 = cpred[npred - i - 1];
        if(fabs(c2 - c1) < eigmet->cgap) continue;

        Citv.x3[0] = c1;
        Citv.x3[1] = 0.5 * (c1 + c2);
        Citv.x3[2] = c2;
        for(size_t p = 0; p < 3; ++p){
            grt_secular_function(mod1d, Citv.x3[p], iref, eigmet->wtype, &Citv.F3[p]);
        }
        stack_push(&stack, Citv);
    }


    // 记录第一个值
    if(eigmet->print_sec){
        fprintf(stdout, "# iref=%zu, f=%f, tol=%e, c1=%f, c2=%f\n", iref, creal(mod1d->omega)/PI2, eigmet->satol, eigmet->cmin, eigmet->cmax);
        fprintf(stdout, "%.16e %.16e %.16e \n", Citv.x3[0], GRT_CMPLX_SPLIT(Citv.F3[0]));
    }

    cplx_t secCheck[6];
    real_t secAbsCheck[6];
    real_t cCheck[6];
    real_t dc, refdc, lastdc=9e30; 

    // 自适应采样
    // 注意这里是k值降序采样
    while(stack.size > 0) {
        // 阶数过多，提前退出
        // if(nroot >= maxnroot)  break;

        Citv = stack_pop(&stack);
        dc = Citv.x3[2] - Citv.x3[0];
        
        // 根据初始采样给出参考步长
        {
            size_t idx = grt_findLessEqualClosest_real_t(cpred, npred, Citv.x3[0]);
            if(idx == npred-1) idx--;
            refdc = GRT_MAX((cpred[idx+1] - cpred[idx])*1e-2, eigmet->cgap);
            // refdc = eigmet->cgap;
        }

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
        grt_secular_function(mod1d, Citv_left.x3[1], iref, eigmet->wtype, &Citv_left.F3[1]);
        grt_secular_function(mod1d, Citv_right.x3[1], iref, eigmet->wtype, &Citv_right.F3[1]);

        // 检查区间是否符合要求
        real_t goodtol = check_fit(&Citv, &Citv_left, &Citv_right, eigmet->satol);

        // 当且仅当 dk 在范围内时才判断 goodtol，否则对于较大dk一律划分子区间，较小dk一律结束区间
        if(dc > refdc && (!goodtol || dc > 2.0*lastdc)){
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
            bool isTrough = false;
            for(int i=2; i<6; ++i){
                if(cCheck[i] < 0.0)  continue;  // 未初始化的右边界点
                if(cCheck[i] < root_c)  continue;  // 该零点已搜索过

                isTrough = secAbsCheck[i-2] >= secAbsCheck[i-1] && secAbsCheck[i-1] <= secAbsCheck[i];
                if(!isTrough)  continue;

                grt_goldensection_search_root(
                    mod1d, cCheck[i-2], cCheck[i], 
                    iref, eigmet->wtype, eigmet->rtol, eigmet->cgap, &root_sec, &root_c
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


// ---------------------------------------------------------------------------------------------------
//                                   固定间隔采样寻找久期函数零点
// ---------------------------------------------------------------------------------------------------

static void grt_fixed_step_secular_roots(
    GRT_MODEL1D *mod1d, EIGENV_INFO *eigmet, EIGENV *eigv, const size_t iref, 
    const real_t *cpred, const size_t npred)
{
    real_t cmin=cpred[0];
    real_t cmax=cpred[npred-1];
    real_t dc = eigmet->uniform_dc;

    cplx_t root_sec = -1.0;
    real_t root_c = 0.0;

    real_t c3[3] = {0};
    cplx_t f3[3] = {0};
    for(int i = 0; i < 2; ++i){
        c3[i] = cmin + dc*i;
        grt_secular_function(mod1d, c3[i], iref, eigmet->wtype, &f3[i]);
    }
    c3[2] = c3[1] + dc;

    if(eigmet->print_sec){
        fprintf(stdout, "# iref=%zu, f=%f, tol=%e, c1=%f, c2=%f\n", iref, creal(mod1d->omega)/PI2, eigmet->satol, eigmet->cmin, eigmet->cmax);
        fprintf(stdout, "%.16e %.16e %.16e \n", c3[0], GRT_CMPLX_SPLIT(f3[0]));
        fprintf(stdout, "%.16e %.16e %.16e \n", c3[1], GRT_CMPLX_SPLIT(f3[1]));
    }
    
    while(c3[2] < cmax) {
        grt_secular_function(mod1d, c3[2], iref, eigmet->wtype, &f3[2]);

        if(eigmet->print_sec){
            fprintf(stdout, "%.16e %.16e %.16e \n", c3[2], GRT_CMPLX_SPLIT(f3[2]));
            goto UPDATE_C3_F3;
        }

        // 找零点
        bool isTrough = (fabs(f3[0]) >= fabs(f3[1])) && (fabs(f3[1]) <= fabs(f3[2]));
        if(! isTrough) goto UPDATE_C3_F3;

        grt_goldensection_search_root(
            mod1d, c3[0], c3[2], 
            iref, eigmet->wtype, eigmet->rtol, eigmet->cgap, &root_sec, &root_c
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
        
        UPDATE_C3_F3:
        for(int i = 0; i < 2; ++i){
            c3[i] = c3[i+1];
            f3[i] = f3[i+1];
        }
        c3[2] = c3[1] + dc;
    }
    
}






// 搜索零点的函数指针
typedef void (*SearchFunc)(
    GRT_MODEL1D *mod1d, EIGENV_INFO *eigmet, EIGENV *eigv, const size_t iref, 
    const real_t *cpred, const size_t npred);


static real_t grt_get_approx_nroots(GRT_MODEL1D *mod1d, const real_t freq, DISPER_TYPE wtype, const real_t cphase)
{
    real_t res = 0;
    grt_mod1d_xa_xb(mod1d, PI2*freq/cphase);

    for(size_t iy = 0; iy < mod1d->n - 1; ++iy){
        res += fabs(cimag(((mod1d->Vb[iy] > 0.0) ? mod1d->xb[iy] : 0.0))) * mod1d->Thk[iy];
        if(wtype == GRT_DISPERSION_RAYL){
            res += fabs(cimag(mod1d->Va[iy])) * mod1d->Thk[iy];
        }
        // printf("xa=" GRT_CMPLX_FMT ", xb=" GRT_CMPLX_FMT "\n", GRT_CMPLX_SPLIT(xa), GRT_CMPLX_SPLIT(xb));
    }
    res *= 2.0*freq;
    res /= cphase;
    return res;
}


// static void get_dc_min_max(GRT_MODEL1D *mod1d, const real_t freq, const real_t cmin, const real_t cmax, real_t *pt_dcmin, real_t *pt_dcmax)
// {
//     size_t approx_n = grt_get_approx_maxnroots(mod1d, freq, cmax);
//     *pt_dcmax = (cmax - cmin) / approx_n;
//     *pt_dcmin = *pt_dcmax * 1e-4;
//     // printf("# freq=%e, dcmin=%e, dcmax=%e, n=%zu, cmax=%f\n", creal(mod1d->omega)/PI2, *pt_dcmin, *pt_dcmax, approx_n, cmax);
// }

/** 将对单个频率的处理包装成函数，方便被其它函数并行调用 */
static void get_secular_roots_single_freq(GRT_MODEL1D *mod1d, EIGENV_INFO *eigmet, EIGENV *eigv, SearchFunc searchfunc)
{
    // 根据经验公式确定初始采样
    real_t *cpred = NULL;
    size_t npred = 0;
    DISPER_TYPE wtype = eigmet->wtype;
    real_t freq = creal(mod1d->omega)/PI2;
    real_t cmin=eigmet->cmin + eigmet->vgap;
    real_t cmax=eigmet->cmax - eigmet->vgap;
    real_t realN = grt_get_approx_nroots(mod1d, freq, wtype, cmax) + 1;
    real_t dc = (cmax - cmin) / realN;
    real_t c1, c2, rn1, rn2;
    c1 = cmin;
    rn1 = grt_get_approx_nroots(mod1d, freq, wtype, cmin);
    c2 = cmax;
    while(c1 < cmax){
        cpred = (real_t *)realloc(cpred, sizeof(real_t)*(npred+1));
        cpred[npred++] = c1;
        c2 = GRT_MIN(c1 + dc, cmax);
        rn2 = grt_get_approx_nroots(mod1d, freq, wtype, c2);
        while(rn2 > 0.0 && fabs(rn2 - rn1) > 0.5){
            c2 = c1 + (c2 - c1)*0.5;
            rn2 = grt_get_approx_nroots(mod1d, freq, wtype, c2);
        }
        c1 = c2; rn1 = rn2;
    }
    cpred = (real_t *)realloc(cpred, sizeof(real_t)*(npred+1));
    cpred[npred++] = c2;
    if(eigmet->wtype == GRT_DISPERSION_RAYL && mod1d->Vb[0] > 0.0){
        real_t target = 0.8*mod1d->Vb[0];
        cpred = (real_t *)realloc(cpred, sizeof(real_t)*(npred+1));
        grt_insertOrdered(cpred, &npred, npred+1, &target, sizeof(real_t), true, grt_compare_real_t);
    }

    int secRaylType = 0; // 使用哪种Rayleigh波的secular function
    size_t iref = 0; // 使用哪一层的secular function

    // 对于首层
    iref = 0;
    if(eigmet->wtype == GRT_DISPERSION_RAYL){
        // Rayleigh 波
        // 搜索“基阶”进行补充
        secRaylType = 0;
        searchfunc(mod1d, eigmet, eigv, secRaylType, iref, cpred, npred);
    } 
    else if(eigmet->wtype == GRT_DISPERSION_LOVE){
        // Love 波 跳过液体层
        if(! mod1d->isLiquid[iref]){
            secRaylType = 1;
            searchfunc(mod1d, eigmet, eigv, secRaylType, iref, cpred, npred);
        }
    }
    else {
        GRTRaiseError("Wrong execution.");
    }

    // 检查低速层
    for(iref=1; iref<mod1d->n-1; ++iref){
        // 跳过两侧固体的递增 Vp 和 递增 Vs
        if( ! mod1d->isLiquid[iref-1] && ! mod1d->isLiquid[iref] && mod1d->Vb[iref-1] <= mod1d->Vb[iref] && mod1d->Va[iref-1] <= mod1d->Va[iref]){
            continue;
        }
        // 跳过两侧液体的递增 Vp
        if( mod1d->isLiquid[iref-1] && mod1d->isLiquid[iref] && mod1d->Va[iref-1] <= mod1d->Va[iref]){
            continue;
        }
        // 对于 Love 波跳过液体层
        if(eigmet->wtype == GRT_DISPERSION_LOVE && mod1d->isLiquid[iref]){
            continue;
        }
        secRaylType = 1;
        searchfunc(mod1d, eigmet, eigv, secRaylType, iref, cpred, npred);
    }

    GRT_SAFE_FREE_PTR(cpred);
}


void grt_get_secular_roots(GRT_MODEL1D *mod1d, EIGENV_INFO *eigmet, const bool print_progressbar)
{
    real_t *freqs = eigmet->freqs;
    size_t nf = eigmet->nf;

    SearchFunc searchfunc = NULL;
    if(eigmet->uniform_dc > 0.0){
        searchfunc = grt_fixed_step_secular_roots;
    } else {
        searchfunc = grt_adaptive_step_secular_roots;
    }

    mod1d->omgref = PI2*freqs[nf-1];
    
    // 进度条变量 
    int progress=0;

#define __CALL_SECULAR_ROOT(mod1d, iw) \
    mod1d->omega = PI2*freqs[iw]; \
    get_secular_roots_single_freq(mod1d, eigmet, eigmet->eigv+iw, searchfunc);

    if(eigmet->print_sec){
        __CALL_SECULAR_ROOT(mod1d, 0);
        goto FINISH;
    }

#ifdef _OPENMP
    // 如果外部已在并行区域中，此处串行
    if (omp_in_parallel()) {
        for(size_t iw = 0; iw < nf; ++iw){
            __CALL_SECULAR_ROOT(mod1d, iw);
            // 记录进度条变量 
            progress++;
            if(print_progressbar) grt_printprogressBar("Computing Dispersion Curves: ", progress*100/nf);
        }
    }
    // 否则正常使用并行
    else{
        #pragma omp parallel for schedule(guided) default(shared) 
        for(size_t iw = 0; iw < nf; ++iw){
            GRT_MODEL1D *local_mod1d = grt_copy_mod1d(mod1d);
            __CALL_SECULAR_ROOT(local_mod1d, iw);

            // 记录进度条变量 
            #pragma omp critical
            {   
                progress++;
                if(print_progressbar) grt_printprogressBar("Computing Dispersion Curves: ", progress*100/nf);
            }
            grt_free_mod1d(local_mod1d);
        }
    }

#else
    // 串行
    for(size_t iw = 0; iw < nf; ++iw){
        __CALL_SECULAR_ROOT(mod1d, iw);
        // 记录进度条变量 
        progress++;
        if(print_progressbar) grt_printprogressBar("Computing Dispersion Curves: ", progress*100/nf);
    }
#endif

FINISH:
}
