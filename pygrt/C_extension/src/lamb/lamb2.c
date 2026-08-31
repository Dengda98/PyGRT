/**
 * @file   lamb2.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 *    使用广义闭合解求解第二类 Lamb 问题，参考：
 *
 *        张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社
 *
 *    本文件中的求解顺序对应第 7.4 节：先构造多项式，再做部分分式
 *    展开，最后组合基本积分。多项式系数没有逐项硬编码，而是直接
 *    按表 7.2.1、表 7.2.2、表 7.3.1 和表 7.3.2 的代数式生成。
 */

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "grt/common/checkerror.h"
#include "grt/lamb/elliptic.h"
#include "grt/lamb/lamb2.h"
#include "grt/lamb/lamb_util.h"


#define LAMB2_POLY_SIZE 16
#define LAMB2_TAIL_SIZE 6
#define LAMB2_SMALL 1e-12


/** 多项式，c[i] 是 x^i 的系数 */
typedef struct {
    cplx_t c[LAMB2_POLY_SIZE];
    int degree;
} LAMB2_POLY;

/** 一个分母为三个二次因子的部分分式展开 */
typedef struct {
    cplx_t pair[3][2];
    cplx_t tail[LAMB2_TAIL_SIZE];
    int ntail;
} LAMB2_PF;

/** 基本积分所属的积分项 */
typedef enum {
    LAMB2_P_TERM,
    LAMB2_S_TERM,
    LAMB2_SP_TERM,
} LAMB2_TERM;

/** 与材料、角度和 Rayleigh 根有关的常量 */
typedef struct {
    real_t nu;
    real_t k;
    real_t k2;
    real_t kp;
    real_t kp2;
    real_t theta;
    real_t phi;
    real_t st;
    real_t ct;
    real_t sf;
    real_t cf;
    real_t theta_c;
    cplx_t y[3];
    cplx_t yp[3];
} LAMB2_VARS;

/** 某一个时间点和某一类积分所需的变量 */
typedef struct {
    LAMB2_TERM term;
    real_t tbar;
    real_t kp2;
    real_t m;
    real_t n;
    real_t xi1;
    real_t xi2;
    real_t z1;
    real_t z2;
    real_t z1sq;
    real_t z2sq;
    real_t m_elliptic;
    real_t c_main;
    real_t c_residue;
    real_t c1;
    real_t c2;
} LAMB2_CTX;


static LAMB2_POLY poly_const(const cplx_t value)
{
    LAMB2_POLY result = {0};
    result.c[0] = value;
    result.degree = 0;
    return result;
}

static LAMB2_POLY poly_x(void)
{
    LAMB2_POLY result = {0};
    result.c[1] = 1.0;
    result.degree = 1;
    return result;
}

static void poly_trim(LAMB2_POLY *poly)
{
    while(poly->degree > 0 && cabs(poly->c[poly->degree]) < 1e-30){
        poly->c[poly->degree] = 0.0;
        --poly->degree;
    }
}

static LAMB2_POLY poly_add(const LAMB2_POLY a, const LAMB2_POLY b)
{
    LAMB2_POLY result = {0};
    result.degree = GRT_MAX(a.degree, b.degree);
    for(int i=0; i<=result.degree; ++i){
        result.c[i] = a.c[i] + b.c[i];
    }
    poly_trim(&result);
    return result;
}

static LAMB2_POLY poly_sub(const LAMB2_POLY a, const LAMB2_POLY b)
{
    LAMB2_POLY result = {0};
    result.degree = GRT_MAX(a.degree, b.degree);
    for(int i=0; i<=result.degree; ++i){
        result.c[i] = a.c[i] - b.c[i];
    }
    poly_trim(&result);
    return result;
}

static LAMB2_POLY poly_scale(const LAMB2_POLY a, const cplx_t scale)
{
    LAMB2_POLY result = {0};
    result.degree = a.degree;
    for(int i=0; i<=a.degree; ++i){
        result.c[i] = a.c[i] * scale;
    }
    poly_trim(&result);
    return result;
}

static LAMB2_POLY poly_mul(const LAMB2_POLY a, const LAMB2_POLY b)
{
    LAMB2_POLY result = {0};
    if(a.degree + b.degree >= LAMB2_POLY_SIZE){
        GRTRaiseError("The polynomial degree is too large in lamb2.\n");
    }
    result.degree = a.degree + b.degree;
    for(int i=0; i<=a.degree; ++i){
        for(int j=0; j<=b.degree; ++j){
            result.c[i+j] += a.c[i] * b.c[j];
        }
    }
    poly_trim(&result);
    return result;
}

static LAMB2_POLY poly_pow(const LAMB2_POLY a, const int n)
{
    LAMB2_POLY result = poly_const(1.0);
    for(int i=0; i<n; ++i){
        result = poly_mul(result, a);
    }
    return result;
}

static void poly_divide(
    const LAMB2_POLY numerator, const LAMB2_POLY denominator,
    LAMB2_POLY *quotient, LAMB2_POLY *remainder)
{
    *quotient = poly_const(0.0);
    *remainder = numerator;
    poly_trim(remainder);

    if(denominator.degree <= 0 && cabs(denominator.c[0]) == 0.0){
        GRTRaiseError("The polynomial denominator is zero in lamb2.\n");
    }

    while(remainder->degree >= denominator.degree &&
          !(remainder->degree == 0 && cabs(remainder->c[0]) < 1e-30)){
        int offset = remainder->degree - denominator.degree;
        cplx_t factor = remainder->c[remainder->degree] / denominator.c[denominator.degree];
        quotient->c[offset] += factor;
        quotient->degree = GRT_MAX(quotient->degree, offset);
        for(int j=0; j<=denominator.degree; ++j){
            remainder->c[j+offset] -= factor * denominator.c[j];
        }
        poly_trim(remainder);
    }
    poly_trim(quotient);
}

/**
 * 按式 (7.2.2) 和式 (7.3.1) 的分母分解式构造多项式
 *
 * R'(x) = 16 k'^2 (x^2+y1)(x^2+y2)(x^2+y3)
 */
static LAMB2_POLY make_rayleigh_denominator(const cplx_t roots[3], const real_t kp2)
{
    LAMB2_POLY x2 = poly_pow(poly_x(), 2);
    LAMB2_POLY denominator = poly_const(16.0 * kp2);
    for(int i=0; i<3; ++i){
        denominator = poly_mul(denominator, poly_add(x2, poly_const(roots[i])));
    }
    return denominator;
}

/**
 * 按式 (7.4.1.2) 和式 (7.2.5) 将任意分子多项式展开为部分分式
 *
 * pair[i][0] * x + pair[i][1] 对应 (x^2+y_i)^{-1}，tail 对应多项式商
 */
static void make_partial_fraction(
    const LAMB2_POLY *numerator, const cplx_t roots[3], const real_t kp2,
    LAMB2_PF *pf)
{
    LAMB2_POLY denominator = make_rayleigh_denominator(roots, kp2);
    LAMB2_POLY quotient, remainder;
    poly_divide(*numerator, denominator, &quotient, &remainder);

    memset(pf, 0, sizeof(*pf));
    pf->ntail = quotient.degree + 1;
    if(pf->ntail > LAMB2_TAIL_SIZE){
        GRTRaiseError("The partial-fraction quotient is too large in lamb2.\n");
    }
    for(int i=0; i<pf->ntail; ++i){
        pf->tail[i] = quotient.c[i];
    }

    for(int i=0; i<3; ++i){
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        cplx_t delta = 16.0 * kp2 * (roots[i] - roots[j]) * (roots[i] - roots[k]);
        cplx_t yi = roots[i];
        cplx_t yi2 = yi * yi;
        cplx_t odd = remainder.c[1] - remainder.c[3] * yi + remainder.c[5] * yi2;
        cplx_t even = remainder.c[0] - remainder.c[2] * yi + remainder.c[4] * yi2;
        pf->pair[i][0] = odd / delta;
        pf->pair[i][1] = even / delta;
    }
}

static real_t clamp_elliptic_parameter(real_t m)
{
    if(m <= 0.0){
        if(m < -1e-9){
            GRTRaiseError("The elliptic parameter is out of range in lamb2: %e.\n", m);
        }
        return DBL_EPSILON;
    }
    if(m >= 1.0){
        if(m > 1.0 + 1e-9){
            GRTRaiseError("The elliptic parameter is out of range in lamb2: %e.\n", m);
        }
        return 1.0 - 16.0 * DBL_EPSILON;
    }
    return m;
}

static bool is_real_complex(const cplx_t value)
{
    return fabs(cimag(value)) <= LAMB2_SMALL * (1.0 + fabs(creal(value)));
}

static int interval_indicator(const real_t x, const real_t left, const real_t right)
{
    return x > left && x < right;
}

static int F1_P(const cplx_t z, const real_t z2sq)
{
    return is_real_complex(z) && interval_indicator(creal(z), -z2sq, 0.0);
}

static int F1_S(const cplx_t z, const real_t z1, const real_t z2sq)
{
    real_t left = GRT_MAX(-z1, -z2sq);
    return is_real_complex(z) && interval_indicator(creal(z), left, 0.0);
}

static int F2(const cplx_t z, const real_t z2sq)
{
    real_t x = creal(z);
    real_t y = cimag(z);
    if(!interval_indicator(x, -z2sq, 0.0) || y <= 0.0){
        return 0;
    }
    real_t ymax = sqrt((z2sq + x) * (1.0 - x));
    return y < ymax;
}

static real_t positive_sqrt(const real_t value, const char *name)
{
    if(value < 0.0){
        if(value > -1e-10){
            return 0.0;
        }
        GRTRaiseError("The discriminant %s is negative in lamb2: %e.\n", name, value);
    }
    return sqrt(value);
}

/** 式 (7.2.11)、(7.3.7) 和 (7.3.18) 中的分式线性变换参数 */
static void make_context_P(const real_t tbar, const LAMB2_VARS *V, LAMB2_CTX *ctx)
{
    ctx->term = LAMB2_P_TERM;
    ctx->tbar = tbar;
    ctx->kp2 = V->kp2;
    ctx->m = tbar * V->ct;
    ctx->n = sqrt(tbar*tbar - V->k2) * V->st;

    real_t d = (tbar*tbar + V->k2*V->ct*V->ct - 1.0) / ctx->m;
    real_t discriminant = d*d + 4.0*V->kp2;
    real_t root_gap = positive_sqrt(discriminant, "P");
    ctx->xi1 = 0.5 * (d + root_gap);
    ctx->xi2 = 0.5 * (d - root_gap);
    ctx->z1sq = ctx->xi1 / (-ctx->xi2);
    ctx->z2sq = (ctx->xi1 - ctx->m) / (ctx->m - ctx->xi2);
    ctx->z1 = sqrt(ctx->z1sq);
    ctx->z2 = sqrt(ctx->z2sq);
    ctx->m_elliptic = clamp_elliptic_parameter(ctx->z2sq / ctx->z1sq);
    ctx->c_main = 1.0 / sqrt(ctx->xi1 * (ctx->m - ctx->xi2));
    ctx->c_residue = 1.0 / sqrt(ctx->xi2 * (ctx->xi2 - ctx->m));
    ctx->c1 = 0.0;
    ctx->c2 = 0.0;
}

static void make_context_S(const real_t tbar, const LAMB2_VARS *V, LAMB2_CTX *ctx)
{
    ctx->term = LAMB2_S_TERM;
    ctx->tbar = tbar;
    ctx->kp2 = V->kp2;
    ctx->m = tbar * V->ct;
    ctx->n = sqrt(tbar*tbar - 1.0) * V->st;

    real_t d = (tbar*tbar - V->k2 + V->ct*V->ct) / ctx->m;
    real_t discriminant = d*d - 4.0*V->kp2;
    real_t root_gap = positive_sqrt(discriminant, "S");
    ctx->xi1 = 0.5 * (d + root_gap);
    ctx->xi2 = 0.5 * (d - root_gap);
    ctx->z1sq = ctx->xi1 / ctx->xi2;
    ctx->z2sq = (ctx->xi1 - ctx->m) / (ctx->m - ctx->xi2);
    ctx->z1 = sqrt(ctx->z1sq);
    ctx->z2 = sqrt(ctx->z2sq);
    ctx->m_elliptic = clamp_elliptic_parameter(ctx->z2sq / (ctx->z1sq + ctx->z2sq));
    ctx->c_main = 1.0 / sqrt(ctx->m * (ctx->xi1 - ctx->xi2));
    ctx->c_residue = 1.0 / sqrt(ctx->xi2 * (ctx->m - ctx->xi2));
    ctx->c1 = 0.0;
    ctx->c2 = 0.0;
}

static void make_context_SP(const real_t tbar, const LAMB2_VARS *V, LAMB2_CTX *ctx)
{
    ctx->term = LAMB2_SP_TERM;
    ctx->tbar = tbar;
    ctx->kp2 = V->kp2;
    ctx->m = tbar * V->ct;
    ctx->n = sqrt(1.0 - tbar*tbar) * V->st;

    real_t d = (tbar*tbar - V->k2 + V->ct*V->ct) / ctx->m;
    real_t discriminant = d*d - 4.0*V->kp2;
    real_t root_gap = positive_sqrt(discriminant, "S-P");
    ctx->xi1 = 0.5 * (d + root_gap);
    ctx->xi2 = 0.5 * (d - root_gap);
    ctx->z1sq = ctx->xi1 / ctx->xi2;
    ctx->z2sq = (ctx->xi1 - ctx->m) / (ctx->xi2 - ctx->m);
    ctx->z1 = sqrt(ctx->z1sq);
    ctx->z2 = sqrt(ctx->z2sq);
    ctx->m_elliptic = clamp_elliptic_parameter(
        ctx->m * (ctx->xi1 - ctx->xi2) / (ctx->xi2 * (ctx->xi1 - ctx->m)));
    ctx->c_main = 1.0 / sqrt(ctx->xi2 * (ctx->xi1 - ctx->m));
    ctx->c_residue = 0.0;
    ctx->c1 = ctx->z2sq - 1.0;
    ctx->c2 = ctx->z2sq - ctx->z1sq;
}

static LAMB2_POLY f1_S(const LAMB2_POLY x2, const LAMB2_POLY gamma3,
    const LAMB2_POLY gp, const LAMB2_POLY etaa2, const LAMB2_POLY z)
{
    LAMB2_POLY first = poly_mul(poly_sub(x2, z), gamma3);
    LAMB2_POLY second = poly_scale(poly_mul(poly_mul(x2, gp), etaa2), 16.0);
    second = poly_mul(second, z);
    return poly_add(first, second);
}

static LAMB2_POLY f2_S(const LAMB2_POLY x2, const LAMB2_POLY gamma,
    const LAMB2_POLY gp, const LAMB2_POLY z)
{
    return poly_add(poly_mul(gp, poly_sub(x2, z)), poly_mul(gamma, z));
}

/**
 * 按表 7.2.1 和表 7.2.2 构造 P 波项的 M^(1)、M^(2) 及其导数
 *
 * 式 (7.2.1) 给出 q、p、eta_beta 和 gamma 的 x 表达式，表中的各项
 * 因而都可以直接用多项式运算生成。这里的导数是式 (7.1.6a)-(7.1.6c)
 * 在 x 变量下的多项式形式。
 */
static void build_P_polynomials(
    const real_t tbar, const LAMB2_VARS *V,
    LAMB2_POLY M[2][3][3], LAMB2_POLY dM[3][2][3][3],
    LAMB2_POLY dM_receiver[3][2][3][3])
{
    LAMB2_POLY x = poly_x();
    LAMB2_POLY x2 = poly_pow(x, 2);
    LAMB2_POLY x3 = poly_pow(x, 3);
    LAMB2_POLY g = poly_sub(x2, poly_const(V->k2));
    LAMB2_POLY etab2 = poly_add(x2, poly_const(V->kp2));
    LAMB2_POLY gamma = poly_add(poly_scale(x2, 2.0), poly_const(1.0 - 2.0*V->k2));
    LAMB2_POLY gamma2 = poly_pow(gamma, 2);
    LAMB2_POLY gamma3 = poly_pow(gamma, 3);
    LAMB2_POLY q = poly_add(poly_scale(x, V->ct / V->st), poly_const(-tbar / V->st));
    LAMB2_POLY Q = poly_add(poly_sub(x2, poly_scale(x, 2.0*tbar*V->ct)),
        poly_const(tbar*tbar - V->k2*V->st*V->st));
    LAMB2_POLY p2 = poly_scale(Q, 1.0 / (V->st*V->st));
    LAMB2_POLY q2 = poly_pow(q, 2);
    LAMB2_POLY eps = poly_sub(poly_scale(q2, V->cf*V->cf), poly_scale(p2, V->sf*V->sf));
    LAMB2_POLY epbar = poly_sub(poly_scale(q2, V->sf*V->sf), poly_scale(p2, V->cf*V->cf));
    LAMB2_POLY zeta = poly_scale(poly_add(q2, p2), V->sf*V->cf);
    LAMB2_POLY kappa = poly_sub(poly_scale(q2, V->cf*V->cf), poly_scale(p2, 3.0*V->sf*V->sf));
    LAMB2_POLY kappabar = poly_sub(poly_scale(q2, V->sf*V->sf), poly_scale(p2, 3.0*V->cf*V->cf));
    LAMB2_POLY kappap = poly_add(kappa, poly_scale(p2, 2.0));
    LAMB2_POLY kappabarp = poly_add(kappabar, poly_scale(p2, 2.0));
    LAMB2_POLY q2_minus_p2 = poly_sub(q2, p2);
    real_t receiver_coef = 1.0 - 2.0*V->k2;

    memset(M, 0, sizeof(LAMB2_POLY) * 2 * 3 * 3);
    memset(dM, 0, sizeof(LAMB2_POLY) * 3 * 2 * 3 * 3);
    memset(dM_receiver, 0, sizeof(LAMB2_POLY) * 3 * 2 * 3 * 3);

    /* 表 7.2.1，xi=1 */
    M[0][0][0] = poly_scale(poly_mul(poly_mul(poly_mul(x2, g), etab2), eps), 8.0);
    M[0][0][1] = poly_scale(poly_mul(poly_mul(poly_mul(x2, g), etab2), zeta), 8.0);
    M[0][0][2] = poly_scale(poly_mul(poly_mul(poly_mul(x3, g), etab2), q), 8.0*V->cf);
    M[0][1][0] = M[0][0][1];
    M[0][1][1] = poly_scale(poly_mul(poly_mul(poly_mul(x2, g), etab2), epbar), 8.0);
    M[0][1][2] = poly_scale(poly_mul(poly_mul(poly_mul(x3, g), etab2), q), 8.0*V->sf);
    M[0][2][0] = poly_scale(poly_mul(poly_mul(x, gamma3), q), V->cf);
    M[0][2][1] = poly_scale(poly_mul(poly_mul(x, gamma3), q), V->sf);
    M[0][2][2] = poly_mul(x2, gamma3);

    /* 表 7.2.1，xi=2 */
    M[1][0][0] = poly_scale(poly_mul(poly_mul(poly_mul(x, etab2), gamma2), eps), 2.0);
    M[1][0][1] = poly_scale(poly_mul(poly_mul(poly_mul(x, etab2), gamma2), zeta), 2.0);
    M[1][0][2] = poly_scale(poly_mul(poly_mul(poly_mul(x2, etab2), gamma2), q), 2.0*V->cf);
    M[1][1][0] = M[1][0][1];
    M[1][1][1] = poly_scale(poly_mul(poly_mul(poly_mul(x, etab2), gamma2), epbar), 2.0);
    M[1][1][2] = poly_scale(poly_mul(poly_mul(poly_mul(x2, etab2), gamma2), q), 2.0*V->sf);
    M[1][2][0] = poly_scale(poly_mul(poly_mul(poly_mul(x2, g), etab2), poly_mul(gamma, q)), 4.0*V->cf);
    M[1][2][1] = poly_scale(poly_mul(poly_mul(poly_mul(x2, g), etab2), poly_mul(gamma, q)), 4.0*V->sf);
    M[1][2][2] = poly_scale(poly_mul(poly_mul(x3, g), poly_mul(etab2, gamma)), 4.0);

    for(int xi=0; xi<2; ++xi){
        LAMB2_POLY x2qgetab2 = poly_mul(poly_mul(x2, q), poly_mul(g, etab2));
        LAMB2_POLY x3getab2 = poly_mul(poly_mul(x3, g), etab2);
        LAMB2_POLY xqetab2gamma2 = poly_mul(poly_mul(x, q), poly_mul(etab2, gamma2));
        LAMB2_POLY x2getab2gamma = poly_mul(poly_mul(x2, g), poly_mul(etab2, gamma));
        LAMB2_POLY x3qgetab2gamma = poly_mul(poly_mul(x3, q), poly_mul(g, poly_mul(etab2, gamma)));
        real_t factor = (xi == 0) ? 8.0 : 2.0;

        /* 表 7.2.2，k'=1 */
        if(xi == 0){
            dM[0][xi][0][0] = poly_scale(poly_mul(x2qgetab2, kappa), -factor*V->cf);
            dM[0][xi][0][1] = poly_scale(poly_mul(x2qgetab2, kappap), -factor*V->sf);
            dM[0][xi][1][1] = poly_scale(poly_mul(x2qgetab2, kappabarp), -factor*V->cf);
            dM[0][xi][0][2] = poly_scale(poly_mul(x3getab2, eps), -factor);
            dM[0][xi][1][2] = poly_scale(poly_mul(x3getab2, zeta), -factor);
            dM[0][xi][2][0] = poly_scale(poly_mul(x, poly_mul(gamma3, eps)), -1.0);
            dM[0][xi][2][1] = poly_scale(poly_mul(x, poly_mul(gamma3, zeta)), -1.0);
            dM[0][xi][2][2] = poly_scale(poly_mul(poly_mul(x2, q), gamma3), -V->cf);
        } else {
            dM[0][xi][0][0] = poly_scale(poly_mul(xqetab2gamma2, kappa), -factor*V->cf);
            dM[0][xi][0][1] = poly_scale(poly_mul(xqetab2gamma2, kappap), -factor*V->sf);
            dM[0][xi][1][1] = poly_scale(poly_mul(xqetab2gamma2, kappabarp), -factor*V->cf);
            dM[0][xi][0][2] = poly_scale(poly_mul(poly_mul(x2, poly_mul(etab2, gamma2)), eps), -factor);
            dM[0][xi][1][2] = poly_scale(poly_mul(poly_mul(x2, poly_mul(etab2, gamma2)), zeta), -factor);
            dM[0][xi][2][0] = poly_scale(poly_mul(x2getab2gamma, eps), -4.0);
            dM[0][xi][2][1] = poly_scale(poly_mul(x2getab2gamma, zeta), -4.0);
            dM[0][xi][2][2] = poly_scale(poly_mul(x3qgetab2gamma, poly_const(V->cf)), -4.0);
        }
        dM[0][xi][1][0] = dM[0][xi][0][1];

        /* 表 7.2.2，k'=2 */
        dM[1][xi][0][0] = dM[0][xi][0][1];
        dM[1][xi][0][1] = dM[0][xi][1][1];
        dM[1][xi][0][2] = dM[0][xi][1][2];
        dM[1][xi][1][0] = dM[1][xi][0][1];
        if(xi == 0){
            dM[1][xi][1][1] = poly_scale(poly_mul(x2qgetab2, kappabar), -factor*V->sf);
            dM[1][xi][1][2] = poly_scale(poly_mul(x3getab2, epbar), -factor);
            dM[1][xi][2][1] = poly_scale(poly_mul(x, poly_mul(gamma3, epbar)), -1.0);
            dM[1][xi][2][2] = poly_scale(poly_mul(poly_mul(x2, q), gamma3), -V->sf);
        } else {
            dM[1][xi][1][1] = poly_scale(poly_mul(xqetab2gamma2, kappabar), -factor*V->sf);
            dM[1][xi][1][2] = poly_scale(poly_mul(poly_mul(x2, poly_mul(etab2, gamma2)), epbar), -factor);
            dM[1][xi][2][1] = poly_scale(poly_mul(x2getab2gamma, epbar), -4.0);
            dM[1][xi][2][2] = poly_scale(poly_mul(x3qgetab2gamma, poly_const(V->sf)), -4.0);
        }
        dM[1][xi][2][0] = dM[0][xi][2][1];

        /* 式 (7.1.6c)，k'=3 */
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                dM[2][xi][i][j] = poly_scale(poly_mul(x, M[xi][i][j]), -1.0);
            }
        }

        /*
         * 接收点导数的前两个方向由水平平移不变性得到，
         * ∂_a G = -∂_{a'} G，a=1,2。其第三个方向使用
         * Johnson (1974) 接收点表达式，并按式 (7.2.3a) 同样有理化
         * 具体关系和式 (12)-(13) 见 lamb2_derivatives.md
         */
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                dM_receiver[0][xi][i][j] = poly_scale(dM[0][xi][i][j], -1.0);
                dM_receiver[1][xi][i][j] = poly_scale(dM[1][xi][i][j], -1.0);
            }
            dM_receiver[2][xi][0][i] = dM[0][xi][2][i];
        }
        dM_receiver[2][xi][1][0] = dM_receiver[2][xi][0][1];
        dM_receiver[2][xi][1][1] = dM[1][xi][2][1];
        dM_receiver[2][xi][1][2] = dM[1][xi][2][2];

        LAMB2_POLY receiver_q;
        LAMB2_POLY receiver_33;
        /* 本文档式 (12) 的 P 波接收点竖向分子 */
        if(xi == 0){
            receiver_q = poly_mul(poly_mul(poly_mul(x2, q), poly_mul(g, etab2)), q2_minus_p2);
            receiver_33 = poly_mul(poly_mul(x3, poly_mul(g, etab2)), q2_minus_p2);
            dM_receiver[2][xi][2][0] = poly_scale(receiver_q, -8.0*receiver_coef*V->cf);
            dM_receiver[2][xi][2][1] = poly_scale(receiver_q, -8.0*receiver_coef*V->sf);
            dM_receiver[2][xi][2][2] = poly_scale(receiver_33, -8.0*receiver_coef);
        } else {
            receiver_q = poly_mul(poly_mul(poly_mul(x, q), poly_mul(etab2, gamma2)), q2_minus_p2);
            receiver_33 = poly_mul(poly_mul(x2, poly_mul(etab2, gamma2)), q2_minus_p2);
            dM_receiver[2][xi][2][0] = poly_scale(receiver_q, -2.0*receiver_coef*V->cf);
            dM_receiver[2][xi][2][1] = poly_scale(receiver_q, -2.0*receiver_coef*V->sf);
            dM_receiver[2][xi][2][2] = poly_scale(receiver_33, -2.0*receiver_coef);
        }
    }
}

/**
 * 按表 7.3.1 和表 7.3.2 构造 S 波项的 N^(1)、N^(2) 及其导数
 *
 * 式 (7.3.1) 统一处理 S1、S2 和 S-P 三段积分。S1/S2 使用同一组
 * 多项式，二者的路径差异已经包含在第 7.3.5 节的基本积分中。
 */
static void build_S_polynomials(
    const real_t tbar, const LAMB2_VARS *V,
    LAMB2_POLY N[2][3][3], LAMB2_POLY dN[3][2][3][3],
    LAMB2_POLY dN_receiver[3][2][3][3])
{
    LAMB2_POLY x = poly_x();
    LAMB2_POLY x2 = poly_pow(x, 2);
    LAMB2_POLY x3 = poly_pow(x, 3);
    LAMB2_POLY gp = poly_sub(x2, poly_const(1.0));
    LAMB2_POLY etaa2 = poly_sub(x2, poly_const(V->kp2));
    LAMB2_POLY gamma = poly_sub(poly_scale(x2, 2.0), poly_const(1.0));
    LAMB2_POLY gamma2 = poly_pow(gamma, 2);
    LAMB2_POLY gamma3 = poly_pow(gamma, 3);
    LAMB2_POLY q = poly_add(poly_scale(x, V->ct / V->st), poly_const(-tbar / V->st));
    LAMB2_POLY Q = poly_add(poly_sub(x2, poly_scale(x, 2.0*tbar*V->ct)),
        poly_const(tbar*tbar - V->st*V->st));
    LAMB2_POLY p2 = poly_scale(Q, 1.0 / (V->st*V->st));
    LAMB2_POLY q2 = poly_pow(q, 2);
    LAMB2_POLY eps = poly_sub(poly_scale(q2, V->cf*V->cf), poly_scale(p2, V->sf*V->sf));
    LAMB2_POLY epbar = poly_sub(poly_scale(q2, V->sf*V->sf), poly_scale(p2, V->cf*V->cf));
    LAMB2_POLY zeta = poly_scale(poly_add(q2, p2), V->sf*V->cf);
    LAMB2_POLY q2_minus_p2 = poly_sub(q2, p2);
    LAMB2_POLY kappa = poly_sub(poly_scale(q2, V->cf*V->cf), poly_scale(p2, 3.0*V->sf*V->sf));
    LAMB2_POLY kappabar = poly_sub(poly_scale(q2, V->sf*V->sf), poly_scale(p2, 3.0*V->cf*V->cf));
    LAMB2_POLY kappap = poly_add(kappa, poly_scale(p2, 2.0));
    LAMB2_POLY kappabarp = poly_add(kappabar, poly_scale(p2, 2.0));
    LAMB2_POLY common1 = poly_sub(gamma3, poly_scale(poly_mul(poly_mul(x2, gp), etaa2), 16.0));
    real_t receiver_coef = 1.0 - 2.0*V->k2;

    memset(N, 0, sizeof(LAMB2_POLY) * 2 * 3 * 3);
    memset(dN, 0, sizeof(LAMB2_POLY) * 3 * 2 * 3 * 3);
    memset(dN_receiver, 0, sizeof(LAMB2_POLY) * 3 * 2 * 3 * 3);

    /* 表 7.3.1，xi=1 */
    N[0][0][0] = f1_S(x2, gamma3, gp, etaa2, epbar);
    N[0][0][1] = poly_mul(common1, zeta);
    N[0][0][2] = poly_scale(poly_mul(poly_mul(x, q), gamma3), -V->cf);
    N[0][1][0] = N[0][0][1];
    N[0][1][1] = f1_S(x2, gamma3, gp, etaa2, eps);
    N[0][1][2] = poly_scale(poly_mul(poly_mul(x, q), gamma3), -V->sf);
    N[0][2][0] = poly_scale(poly_mul(poly_mul(poly_mul(x3, gp), etaa2), q), -8.0*V->cf);
    N[0][2][1] = poly_scale(poly_mul(poly_mul(poly_mul(x3, gp), etaa2), q), -8.0*V->sf);
    N[0][2][2] = poly_scale(poly_mul(poly_mul(x2, poly_pow(gp, 2)), etaa2), -8.0);

    /* 表 7.3.1，xi=2 */
    N[1][0][0] = poly_scale(poly_mul(poly_mul(poly_mul(x, etaa2), gamma), f2_S(x2, gamma, gp, epbar)), 4.0);
    N[1][0][1] = poly_scale(poly_mul(poly_mul(poly_mul(x, etaa2), gamma),
        poly_mul(poly_sub(gamma, gp), zeta)), -4.0);
    N[1][0][2] = poly_scale(poly_mul(poly_mul(poly_mul(x2, gp), etaa2), poly_mul(gamma, q)), -4.0*V->cf);
    N[1][1][0] = N[1][0][1];
    N[1][1][1] = poly_scale(poly_mul(poly_mul(poly_mul(x, etaa2), gamma), f2_S(x2, gamma, gp, eps)), 4.0);
    N[1][1][2] = poly_scale(poly_mul(poly_mul(poly_mul(x2, gp), etaa2), poly_mul(gamma, q)), -4.0*V->sf);
    N[1][2][0] = poly_scale(poly_mul(poly_mul(poly_mul(x2, etaa2), gamma2), q), -2.0*V->cf);
    N[1][2][1] = poly_scale(poly_mul(poly_mul(poly_mul(x2, etaa2), gamma2), q), -2.0*V->sf);
    N[1][2][2] = poly_scale(poly_mul(poly_mul(x, gp), poly_mul(etaa2, gamma2)), -2.0);

    for(int xi=0; xi<2; ++xi){
        LAMB2_POLY xetaagamma = poly_mul(x, poly_mul(etaa2, gamma));

        if(xi == 0){
            /* 表 7.3.2，xi=1，k'=1 */
            dN[0][xi][0][0] = poly_scale(poly_mul(f1_S(x2, gamma3, gp, etaa2, kappabarp), poly_mul(q, poly_const(V->cf))), -1.0);
            dN[0][xi][0][1] = poly_scale(poly_mul(poly_mul(common1, kappap), poly_mul(q, poly_const(V->sf))), -1.0);
            dN[0][xi][0][2] = poly_mul(poly_mul(x, poly_mul(eps, gamma3)), poly_const(1.0));
            dN[0][xi][1][1] = poly_scale(poly_mul(f1_S(x2, gamma3, gp, etaa2, kappa), poly_mul(q, poly_const(V->cf))), -1.0);
            dN[0][xi][1][2] = poly_mul(poly_mul(x, poly_mul(zeta, gamma3)), poly_const(1.0));
            dN[0][xi][2][0] = poly_scale(poly_mul(poly_mul(poly_mul(x3, gp), etaa2), eps), 8.0);
            dN[0][xi][2][1] = poly_scale(poly_mul(poly_mul(poly_mul(x3, gp), etaa2), zeta), 8.0);
            dN[0][xi][2][2] = poly_scale(poly_mul(poly_mul(poly_mul(x2, poly_pow(gp, 2)), etaa2), poly_mul(q, poly_const(V->cf))), 8.0);

            /* 表 7.3.2，xi=1，k'=2 */
            dN[1][xi][0][0] = poly_scale(poly_mul(f1_S(x2, gamma3, gp, etaa2, kappabar), poly_mul(q, poly_const(V->sf))), -1.0);
            dN[1][xi][0][1] = poly_scale(poly_mul(poly_mul(common1, kappabarp), poly_mul(q, poly_const(V->cf))), -1.0);
            dN[1][xi][0][2] = dN[0][xi][1][2];
            dN[1][xi][1][1] = poly_scale(poly_mul(f1_S(x2, gamma3, gp, etaa2, kappap), poly_mul(q, poly_const(V->sf))), -1.0);
            dN[1][xi][1][2] = poly_mul(poly_mul(x, poly_mul(epbar, gamma3)), poly_const(1.0));
            dN[1][xi][2][0] = dN[0][xi][2][1];
            dN[1][xi][2][1] = poly_scale(poly_mul(poly_mul(poly_mul(x3, gp), etaa2), epbar), 8.0);
            dN[1][xi][2][2] = poly_scale(poly_mul(poly_mul(poly_mul(x2, poly_pow(gp, 2)), etaa2), poly_mul(q, poly_const(V->sf))), 8.0);
        } else {
            /* 表 7.3.2，xi=2，k'=1 */
            dN[0][xi][0][0] = poly_scale(poly_mul(poly_mul(xetaagamma, f2_S(x2, gamma, gp, kappabarp)), poly_mul(q, poly_const(V->cf))), -4.0);
            dN[0][xi][0][1] = poly_scale(poly_mul(poly_mul(xetaagamma, poly_mul(poly_sub(gp, gamma), kappap)), poly_mul(q, poly_const(V->sf))), -4.0);
            dN[0][xi][0][2] = poly_scale(poly_mul(poly_mul(poly_mul(x2, gp), etaa2), poly_mul(gamma, eps)), 4.0);
            dN[0][xi][1][1] = poly_scale(poly_mul(poly_mul(xetaagamma, f2_S(x2, gamma, gp, kappa)), poly_mul(q, poly_const(V->cf))), -4.0);
            dN[0][xi][1][2] = poly_scale(poly_mul(poly_mul(poly_mul(x2, gp), etaa2), poly_mul(gamma, zeta)), 4.0);
            dN[0][xi][2][0] = poly_scale(poly_mul(poly_mul(x2, poly_mul(etaa2, gamma2)), eps), 2.0);
            dN[0][xi][2][1] = poly_scale(poly_mul(poly_mul(x2, poly_mul(etaa2, gamma2)), zeta), 2.0);
            dN[0][xi][2][2] = poly_scale(poly_mul(poly_mul(poly_mul(x, gp), etaa2), poly_mul(gamma2, poly_mul(q, poly_const(V->cf)))), 2.0);

            /* 表 7.3.2，xi=2，k'=2 */
            dN[1][xi][0][0] = poly_scale(poly_mul(poly_mul(xetaagamma, f2_S(x2, gamma, gp, kappabar)), poly_mul(q, poly_const(V->sf))), -4.0);
            dN[1][xi][0][1] = poly_scale(
                poly_mul(poly_mul(xetaagamma, poly_mul(poly_sub(gp, gamma), kappabarp)),
                    poly_mul(q, poly_const(V->cf))), -4.0);
            dN[1][xi][0][2] = dN[0][xi][1][2];
            dN[1][xi][1][1] = poly_scale(poly_mul(poly_mul(xetaagamma, f2_S(x2, gamma, gp, kappap)), poly_mul(q, poly_const(V->sf))), -4.0);
            dN[1][xi][1][2] = poly_scale(poly_mul(poly_mul(poly_mul(x2, gp), etaa2), poly_mul(gamma, epbar)), 4.0);
            dN[1][xi][2][0] = dN[0][xi][2][1];
            dN[1][xi][2][1] = poly_scale(poly_mul(poly_mul(x2, poly_mul(etaa2, gamma2)), epbar), 2.0);
            dN[1][xi][2][2] = poly_scale(poly_mul(poly_mul(poly_mul(x, gp), etaa2), poly_mul(gamma2, poly_mul(q, poly_const(V->sf)))), 2.0);
        }

        dN[0][xi][1][0] = dN[0][xi][0][1];
        dN[1][xi][1][0] = dN[1][xi][0][1];

        /* 式 (7.3.1) 中的 N_(ij,3')^(xi)=-x N_ij^(xi) */
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                dN[2][xi][i][j] = poly_scale(poly_mul(x, N[xi][i][j]), -1.0);
            }
        }

        /* 接收点 1、2 方向沿用式 (4) 的水平平移关系，第三方向按式 (7)、(9) 构造 */
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                dN_receiver[0][xi][i][j] = poly_scale(dN[0][xi][i][j], -1.0);
                dN_receiver[1][xi][i][j] = poly_scale(dN[1][xi][i][j], -1.0);
            }
            dN_receiver[2][xi][0][i] = dN[0][xi][2][i];
        }
        dN_receiver[2][xi][1][0] = dN_receiver[2][xi][0][1];
        dN_receiver[2][xi][1][1] = dN[1][xi][2][1];
        dN_receiver[2][xi][1][2] = dN[1][xi][2][2];

        LAMB2_POLY receiver_q;
        LAMB2_POLY receiver_33;
        /* 本文档式 (16) 的 S 波和 S-P 波接收点竖向分子 */
        if(xi == 0){
            receiver_q = poly_mul(poly_mul(x2, q), gamma3);
            receiver_33 = poly_mul(poly_mul(x, gamma3), q2_minus_p2);
            dN_receiver[2][xi][2][0] = poly_scale(receiver_q, -receiver_coef*V->cf);
            dN_receiver[2][xi][2][1] = poly_scale(receiver_q, -receiver_coef*V->sf);
            dN_receiver[2][xi][2][2] = poly_scale(receiver_33, receiver_coef);
        } else {
            receiver_q = poly_mul(poly_mul(poly_mul(x3, q), poly_mul(gp, etaa2)), gamma);
            receiver_33 = poly_mul(poly_mul(x2, poly_mul(gp, etaa2)), poly_mul(gamma, q2_minus_p2));
            dN_receiver[2][xi][2][0] = poly_scale(receiver_q, -4.0*receiver_coef*V->cf);
            dN_receiver[2][xi][2][1] = poly_scale(receiver_q, -4.0*receiver_coef*V->sf);
            dN_receiver[2][xi][2][2] = poly_scale(receiver_33, 4.0*receiver_coef);
        }

    }
}

/** 式 (7.2.8)-(7.2.10b) 和式 (7.3.5) 的 U 基本积分 */
static cplx_t basic_U(const int number, const cplx_t c, const LAMB2_CTX *ctx)
{
    if(ctx->term == LAMB2_SP_TERM){
        /* 式 (7.3.4) 和 7.3.4.2：这些 U 项对最终虚部没有贡献 */
        return 0.0;
    }

    if(number == 1){
        real_t n = ctx->n;
        cplx_t m1 = sqrt(c) - I*ctx->m;
        cplx_t m2 = sqrt(c) + I*ctx->m;
        cplx_t A1 = sqrt((m1 - n) / (m1 + n));
        cplx_t A2 = sqrt((m2 + n) / (m2 - n));
        return A1/(m1 - n)*atan(A1) - A2/(m2 + n)*atan(A2);
    }
    if(number == 2){
        real_t n = ctx->n;
        cplx_t m1 = sqrt(c) - I*ctx->m;
        cplx_t m2 = sqrt(c) + I*ctx->m;
        cplx_t A1 = sqrt((m1 - n) / (m1 + n));
        cplx_t A2 = sqrt((m2 + n) / (m2 - n));
        return I/sqrt(c) * (A1/(m1 - n)*atan(A1) + A2/(m2 + n)*atan(A2));
    }

    real_t m = ctx->m;
    real_t n = ctx->n;
    if(number == 3){
        return I * HALFPI;
    }
    if(number == 4){
        return I * HALFPI * m;
    }
    if(number == 5){
        return I * PI/4.0 * (2.0*m*m - n*n);
    }
    if(number == 6){
        return I * PI/4.0 * m * (2.0*m*m - 3.0*n*n);
    }

    GRTRaiseError("Wrong U basic-integral number in lamb2: %d.\n", number);
}

static void calculate_H(const LAMB2_CTX *ctx, real_t H[5])
{
    real_t m = ctx->m_elliptic;
    real_t a;
    if(ctx->term == LAMB2_P_TERM){
        a = 1.0 / ctx->z2sq;
    } else if(ctx->term == LAMB2_S_TERM){
        a = -(ctx->z2sq + 1.0) / ctx->z2sq;
    } else {
        a = -ctx->c1 / ctx->c2;
    }

    real_t K = grt_ellipticK(m);
    real_t E = grt_ellipticE(m);
    real_t H0 = K;
    real_t H1 = grt_ellipticPi(-1.0/a, m) / a;
    real_t Hminus1 = ((a*m + 1.0)*K - E) / m;
    real_t gamma1 = 3.0*m*a*a + 2.0*a*(m + 1.0) + 1.0;
    real_t gamma2 = m + 1.0 + 3.0*m*a;
    real_t gamma3 = a*(a + 1.0)*(a*m + 1.0);
    real_t H2 = (gamma1*H1 - m*Hminus1) / (2.0*gamma3);
    real_t H3 = (3.0*gamma1*H2 - 2.0*gamma2*H1 + m*H0) / (4.0*gamma3);
    real_t H4 = (5.0*gamma1*H3 - 4.0*gamma2*H2 + 3.0*m*H1) / (6.0*gamma3);

    H[0] = H0;
    H[1] = H1;
    H[2] = H2;
    H[3] = H3;
    H[4] = H4;
}

static real_t tail_V_P(const int number, const LAMB2_CTX *ctx)
{
    /* 式 (7.2.29a)-(7.2.29e) */
    real_t H[5];
    calculate_H(ctx, H);
    real_t xi1 = ctx->xi1;
    real_t xi2 = ctx->xi2;
    real_t z2 = ctx->z2;
    real_t z2m2 = 1.0/(z2*z2);
    real_t z2m4 = z2m2*z2m2;
    real_t z2m6 = z2m4*z2m2;
    real_t z2m8 = z2m6*z2m2;
    real_t beta11 = xi1 - xi2;
    real_t beta12 = xi1 - 2.0*xi2;
    real_t beta13 = xi1 - 3.0*xi2;
    real_t beta15 = xi1 - 5.0*xi2;
    real_t beta35 = 3.0*xi1 - 5.0*xi2;
    real_t result;

    if(number == 3){
        result = ctx->c_main * H[0];
    } else if(number == 4){
        result = ctx->c_main * (xi2*H[0] + beta11*z2m2*H[1]);
    } else if(number == 5){
        result = ctx->c_main * (xi2*xi2*H[0] - beta11*beta13*z2m2*H[1]
            + 2.0*beta11*beta11*z2m4*H[2]);
    } else if(number == 6){
        result = ctx->c_main * (xi2*xi2*xi2*H[0]
            - 3.0*xi2*beta11*beta12*z2m2*H[1]
            - 3.0*beta11*beta11*beta13*z2m4*H[2]
            + 4.0*beta11*beta11*beta11*z2m6*H[3]);
    } else if(number == 7){
        result = ctx->c_main * (xi2*xi2*xi2*xi2*H[0]
            - 2.0*xi2*xi2*beta11*beta35*z2m2*H[1]
            + beta11*beta11*(beta15*beta15 + 4.0*ctx->kp2)*z2m4*H[2]
            - 8.0*beta11*beta11*beta11*beta13*z2m6*H[3]
            + 8.0*beta11*beta11*beta11*beta11*z2m8*H[4]);
    } else {
        GRTRaiseError("Wrong P-wave V basic-integral number in lamb2: %d.\n", number);
    }
    return result;
}

static real_t tail_V_S(const int number, const LAMB2_CTX *ctx)
{
    /* 式 (7.3.17a)-(7.3.17e) */
    real_t H[5];
    calculate_H(ctx, H);
    real_t xi1 = ctx->xi1;
    real_t xi2 = ctx->xi2;
    real_t z2 = ctx->z2;
    real_t z2m2 = 1.0/(z2*z2);
    real_t z2m4 = z2m2*z2m2;
    real_t z2m6 = z2m4*z2m2;
    real_t z2m8 = z2m6*z2m2;
    real_t beta11 = xi1 - xi2;
    real_t beta12 = xi1 - 2.0*xi2;
    real_t beta13 = xi1 - 3.0*xi2;
    real_t beta15 = xi1 - 5.0*xi2;
    real_t beta35 = 3.0*xi1 - 5.0*xi2;
    real_t result;

    if(number == 3){
        result = ctx->c_main * H[0];
    } else if(number == 4){
        result = ctx->c_main * (xi2*H[0] - beta11*z2m2*H[1]);
    } else if(number == 5){
        result = ctx->c_main * (xi2*xi2*H[0] + beta11*beta13*z2m2*H[1]
            + 2.0*beta11*beta11*z2m4*H[2]);
    } else if(number == 6){
        result = ctx->c_main * (xi2*xi2*xi2*H[0]
            + 3.0*xi2*beta11*beta12*z2m2*H[1]
            - 3.0*beta11*beta11*beta13*z2m4*H[2]
            - 4.0*beta11*beta11*beta11*z2m6*H[3]);
    } else if(number == 7){
        result = ctx->c_main * (xi2*xi2*xi2*xi2*H[0]
            + 2.0*xi2*xi2*beta11*beta35*z2m2*H[1]
            + beta11*beta11*(beta15*beta15 - 4.0*ctx->kp2)*z2m4*H[2]
            + 8.0*beta11*beta11*beta11*beta13*z2m6*H[3]
            + 8.0*beta11*beta11*beta11*beta11*z2m8*H[4]);
    } else {
        GRTRaiseError("Wrong S-wave V basic-integral number in lamb2: %d.\n", number);
    }
    return result;
}

static real_t tail_V_SP(const int number, const LAMB2_CTX *ctx)
{
    /* 式 (7.3.22a)-(7.3.22e)，M1-M4 按式 (E.7)、(E.9a)-(E.9c) */
    real_t H[5];
    calculate_H(ctx, H);
    real_t xi1 = ctx->xi1;
    real_t xi2 = ctx->xi2;
    real_t z2 = ctx->z2;
    real_t beta11 = xi1 - xi2;
    real_t beta12 = xi1 - 2.0*xi2;
    real_t beta13 = xi1 - 3.0*xi2;
    real_t beta14 = xi1 - 4.0*xi2;
    real_t beta15 = xi1 - 5.0*xi2;
    real_t beta17 = xi1 - 7.0*xi2;
    real_t beta35 = 3.0*xi1 - 5.0*xi2;
    real_t c1 = ctx->c1;
    real_t c2 = ctx->c2;
    real_t c1m1 = 1.0/c2;
    real_t c1m2 = c1m1*c1m1;
    real_t c1m3 = c1m2*c1m1;
    real_t c1m4 = c1m3*c1m1;
    real_t M1 = PI/(2.0*sqrt(c1)*sqrt(c1-c2));
    real_t M2 = PI*(2.0*c1-c2)/(4.0*pow(c1*(c1-c2), 1.5));
    real_t M3 = PI*(8.0*c1*c1 - 8.0*c1*c2 + 3.0*c2*c2)
        /(16.0*pow(c1*(c1-c2), 2.5));
    real_t M4 = PI*(2.0*c1-c2)*(8.0*c1*c1 - 8.0*c1*c2 + 5.0*c2*c2)
        /(32.0*pow(c1*(c1-c2), 3.5));
    real_t result;

    if(number == 3){
        result = -ctx->c_main * H[0];
    } else if(number == 4){
        result = -ctx->c_main * (xi2*H[0] + beta11*c1m1*H[1]
            + beta11*z2*M1);
    } else if(number == 5){
        result = -ctx->c_main * (xi2*xi2*H[0] - beta11*beta13*c1m1*H[1]
            + 2.0*beta11*beta11*c1m2*H[2]
            + 2.0*beta11*z2*(xi2*M1 - beta11*M2));
    } else if(number == 6){
        result = -ctx->c_main * (xi2*xi2*xi2*H[0]
            - 3.0*xi2*beta11*beta12*c1m1*H[1]
            - 3.0*beta11*beta11*beta13*c1m2*H[2]
            + 4.0*beta11*beta11*beta11*c1m3*H[3]
            + beta11*z2*(3.0*xi2*xi2*M1 + beta11*beta17*M2
                         + 4.0*beta11*beta11*M3));
    } else if(number == 7){
        result = -ctx->c_main * (xi2*xi2*xi2*xi2*H[0]
            - 2.0*xi2*xi2*beta11*beta35*c1m1*H[1]
            + beta11*beta11*(beta15*beta15 - 4.0*ctx->kp2)*c1m2*H[2]
            - 8.0*beta11*beta11*beta11*beta13*c1m3*H[3]
            + 8.0*beta11*beta11*beta11*beta11*c1m4*H[4]
            + 4.0*beta11*z2*(xi2*xi2*xi2*M1 + xi2*beta11*beta14*M2
                         - beta11*beta11*beta15*M3 - 2.0*beta11*beta11*beta11*M4));
    } else {
        GRTRaiseError("Wrong S-P-wave V basic-integral number in lamb2: %d.\n", number);
    }
    return result;
}

typedef struct {
    cplx_t h1p;
    cplx_t h2p;
    cplx_t h1m;
    cplx_t h2m;
    cplx_t D;
    cplx_t s[2];
    cplx_t z0[2];
} LAMB2_V_AUX;

static void make_V_aux(const cplx_t c, const LAMB2_CTX *ctx, LAMB2_V_AUX *aux)
{
    /* 式 (7.2.15)、式 (7.3.12) 和式 (7.3.20) 中的 z0、s 及 h 变量 */
    cplx_t root = sqrt(-c);
    cplx_t gap = ctx->xi1 - ctx->xi2;
    cplx_t h2psq;

    aux->h1p = c + ctx->xi1*ctx->xi1;
    aux->h2p = c + ctx->xi2*ctx->xi2;
    aux->h1m = c - ctx->xi1*ctx->xi1;
    aux->h2m = c - ctx->xi2*ctx->xi2;
    h2psq = aux->h2p * aux->h2p;
    aux->D = gap*gap - 2.0*(c + ctx->xi1*ctx->xi2);
    aux->s[0] = -(c + ctx->xi1*ctx->xi2 + gap*root)
        * (c + ctx->xi1*ctx->xi2 + gap*root) / h2psq;
    aux->s[1] = -(c + ctx->xi1*ctx->xi2 - gap*root)
        * (c + ctx->xi1*ctx->xi2 - gap*root) / h2psq;
    aux->z0[0] = (ctx->xi1 + root) / (ctx->xi2 + root);
    aux->z0[1] = (ctx->xi1 - root) / (ctx->xi2 - root);
}

static cplx_t residue_P(
    const int number, const int sign, const cplx_t c,
    const LAMB2_CTX *ctx, const LAMB2_V_AUX *aux)
{
    int other = 1 - sign;
    cplx_t z = aux->z0[sign];
    cplx_t denominator = (c + ctx->xi2*ctx->xi2) * (z - aux->z0[other]);
    denominator *= sqrt(z*z + ctx->z1sq) * sqrt(z*z + ctx->z2sq);
    cplx_t numerator = (number == 1) ? (z - 1.0)*(ctx->xi2*z - ctx->xi1)
                                     : (z - 1.0)*(z - 1.0);
    return ctx->c_residue * numerator / denominator;
}

static cplx_t residue_S(
    const int number, const int sign, const cplx_t c,
    const LAMB2_CTX *ctx, const LAMB2_V_AUX *aux)
{
    int other = 1 - sign;
    cplx_t z = aux->z0[sign];
    cplx_t denominator = (c + ctx->xi2*ctx->xi2) * (z - aux->z0[other]);
    denominator *= sqrt(ctx->z1sq - z*z) * sqrt(z*z + ctx->z2sq);
    cplx_t numerator = (number == 1) ? (z - 1.0)*(ctx->xi2*z - ctx->xi1)
                                     : (z - 1.0)*(z - 1.0);
    return ctx->c_residue * numerator / denominator;
}

static cplx_t V_P_pair(const int number, const cplx_t c, const LAMB2_CTX *ctx)
{
    /* 式 (7.2.24)、式 (7.2.26)-(7.2.28) */
    LAMB2_V_AUX aux;
    make_V_aux(c, ctx, &aux);

    cplx_t h2psq = aux.h2p * aux.h2p;
    cplx_t zeta0 = ctx->xi2 / aux.h2p;
    cplx_t eta0 = 1.0 / aux.h2p;
    cplx_t zeta[2];
    cplx_t eta[2];
    cplx_t zeta_bar[2];
    cplx_t eta_bar[2];
    cplx_t base = 0.0;
    cplx_t branch = 0.0;
    cplx_t alpha11, alpha12, alpha21, alpha22;
    real_t m = ctx->m_elliptic;

    for(int i=0; i<2; ++i){
        int other = 1 - i;
        cplx_t s = aux.s[i];
        cplx_t denominator = s * (aux.s[other] - s) * h2psq;
        zeta[i] = (ctx->xi1 + s*ctx->xi2)
            * (s*aux.h2p + aux.h1p) / denominator;
        zeta_bar[i] = (ctx->xi1 - ctx->xi2)
            * (s*aux.h2m + aux.h1m) / ((aux.s[other] - s)*h2psq);
        eta[i] = (aux.h2p*s*s - aux.D*s + aux.h1p) / denominator;
        eta_bar[i] = 2.0*(ctx->xi1 - ctx->xi2)*(s*ctx->xi2 + ctx->xi1)
            / ((aux.s[other] - s)*h2psq);
    }

    cplx_t n1 = ctx->z2sq / aux.s[0];
    cplx_t n2 = ctx->z2sq / aux.s[1];
    cplx_t elliptic = 0.0;
    if(number == 1){
        elliptic = zeta0*grt_ellipticK(m)
            + zeta[0]*grt_ellipticPi(n1, m)
            + zeta[1]*grt_ellipticPi(n2, m);
    } else {
        elliptic = eta0*grt_ellipticK(m)
            + eta[0]*grt_ellipticPi(n1, m)
            + eta[1]*grt_ellipticPi(n2, m);
    }

    alpha11 = sqrt(ctx->z1sq - aux.s[0]);
    alpha12 = sqrt(ctx->z1sq - aux.s[1]);
    alpha21 = sqrt(ctx->z2sq - aux.s[0]);
    alpha22 = sqrt(ctx->z2sq - aux.s[1]);
    if(number == 1){
        branch = ctx->z1*zeta_bar[0]/(alpha11*alpha21)
            * atan(I*ctx->z2*alpha11/(ctx->z1*alpha21))
            + ctx->z1*zeta_bar[1]/(alpha12*alpha22)
            * atan(I*ctx->z2*alpha12/(ctx->z1*alpha22));
    } else {
        branch = ctx->z1*eta_bar[0]/(alpha11*alpha21)
            * atan(I*ctx->z2*alpha11/(ctx->z1*alpha21))
            + ctx->z1*eta_bar[1]/(alpha12*alpha22)
            * atan(I*ctx->z2*alpha12/(ctx->z1*alpha22));
    }

    if(is_real_complex(c)){
        real_t value = ctx->c_main * creal(elliptic);
        for(int sign=0; sign<2; ++sign){
            real_t d = creal(residue_P(number, sign, c, ctx, &aux));
            if(creal(c) < 0.0){
                value -= PI * F1_P(aux.z0[sign], ctx->z2sq) * d;
            } else if(creal(c) > 0.0){
                value -= 2.0*PI * F2(aux.z0[sign], ctx->z2sq) * d;
            }
        }
        return I*value;
    }

    base = I*ctx->c_main*elliptic;
    if(number == 1){
        base -= I*ctx->c_main*branch;
    } else {
        base += I*ctx->c_main*branch;
    }
    base -= 2.0*PI*I * (
        F2(aux.z0[0], ctx->z2sq)*residue_P(number, 0, c, ctx, &aux)
        + F2(aux.z0[1], ctx->z2sq)*residue_P(number, 1, c, ctx, &aux));
    return base;
}

static cplx_t V_S_pair(const int number, const cplx_t c, const LAMB2_CTX *ctx)
{
    /* 式 (7.3.13)-(7.3.16)，这里的 V^S 已包含 S1-H(theta-theta_c)S2 */
    LAMB2_V_AUX aux;
    make_V_aux(c, ctx, &aux);

    cplx_t h2psq = aux.h2p * aux.h2p;
    cplx_t zeta0 = ctx->xi2 / aux.h2p;
    cplx_t eta0 = 1.0 / aux.h2p;
    cplx_t zeta[2];
    cplx_t eta[2];
    cplx_t zeta_bar[2];
    cplx_t eta_bar[2];
    cplx_t elliptic = 0.0;
    cplx_t branch = 0.0;
    real_t m = ctx->m_elliptic;

    for(int i=0; i<2; ++i){
        int other = 1 - i;
        cplx_t s = aux.s[i];
        cplx_t denominator = (aux.s[other] - s) * h2psq * (s - ctx->z2sq);
        zeta[i] = (ctx->xi1 + s*ctx->xi2)
            * (s*aux.h2p + aux.h1p) / denominator;
        zeta_bar[i] = (ctx->xi1 - ctx->xi2)
            * (s*aux.h2m + aux.h1m) / ((aux.s[other] - s)*h2psq);
        eta[i] = (aux.h2p*s*s - aux.D*s + aux.h1p) / denominator;
        eta_bar[i] = 2.0*(ctx->xi1 - ctx->xi2)*(s*ctx->xi2 + ctx->xi1)
            / ((aux.s[other] - s)*h2psq);
    }

    cplx_t n1 = ctx->z2sq / (ctx->z2sq - aux.s[0]);
    cplx_t n2 = ctx->z2sq / (ctx->z2sq - aux.s[1]);
    if(number == 1){
        elliptic = zeta0*grt_ellipticK(m)
            + zeta[0]*grt_ellipticPi(n1, m)
            + zeta[1]*grt_ellipticPi(n2, m);
    } else {
        elliptic = eta0*grt_ellipticK(m)
            + eta[0]*grt_ellipticPi(n1, m)
            + eta[1]*grt_ellipticPi(n2, m);
    }

    cplx_t alpha11 = sqrt(aux.s[0] + ctx->z1sq);
    cplx_t alpha12 = sqrt(aux.s[0] - ctx->z2sq);
    cplx_t alpha21 = sqrt(aux.s[1] + ctx->z1sq);
    cplx_t alpha22 = sqrt(aux.s[1] - ctx->z2sq);
    if(number == 1){
        branch = zeta_bar[0]/(alpha11*alpha12)
            * atan(ctx->z2*alpha11/(ctx->z1*alpha12))
            + zeta_bar[1]/(alpha21*alpha22)
            * atan(ctx->z2*alpha21/(ctx->z1*alpha22));
    } else {
        branch = eta_bar[0]/(alpha11*alpha12)
            * atan(ctx->z2*alpha11/(ctx->z1*alpha12))
            + eta_bar[1]/(alpha21*alpha22)
            * atan(ctx->z2*alpha21/(ctx->z1*alpha22));
    }

    if(is_real_complex(c)){
        real_t value = ctx->c_main * creal(elliptic);
        for(int sign=0; sign<2; ++sign){
            real_t d = creal(residue_S(number, sign, c, ctx, &aux));
            if(creal(c) < 0.0){
                value -= PI * F1_S(aux.z0[sign], ctx->z1, ctx->z2sq) * d;
            } else if(creal(c) > 0.0){
                value -= 2.0*PI * F2(aux.z0[sign], ctx->z2sq) * d;
            }
        }
        return I*value;
    }

    cplx_t result = I*ctx->c_main*elliptic;
    if(number == 1){
        result -= branch / sqrt(ctx->xi2*(ctx->m - ctx->xi2));
    } else {
        result += branch / sqrt(ctx->xi2*(ctx->m - ctx->xi2));
    }
    result -= 2.0*PI*I * (
        F2(aux.z0[0], ctx->z2sq)*residue_S(number, 0, c, ctx, &aux)
        + F2(aux.z0[1], ctx->z2sq)*residue_S(number, 1, c, ctx, &aux));
    return result;
}

static cplx_t V_SP_pair(const int number, const cplx_t c, const LAMB2_CTX *ctx)
{
    /* 式 (7.3.20) 和式 (7.3.21) */
    LAMB2_V_AUX aux;
    make_V_aux(c, ctx, &aux);
    cplx_t h2psq = aux.h2p * aux.h2p;
    cplx_t zeta0 = ctx->xi2 / aux.h2p;
    cplx_t eta0 = 1.0 / aux.h2p;
    cplx_t z0p2 = aux.z0[0] * aux.z0[0];
    cplx_t z0m2 = aux.z0[1] * aux.z0[1];
    cplx_t denominator = h2psq * (z0p2 - z0m2);
    cplx_t denominator_reverse = -denominator;
    cplx_t lambda_p, lambda_m, lambda_bar_p, lambda_bar_m;
    cplx_t lambda_prime_p, lambda_prime_m;
    cplx_t lambda_bar_prime_p, lambda_bar_prime_m;
    cplx_t mu1p, mu2p, mu1m, mu2m;
    cplx_t Pi_p, Pi_m;

    lambda_p = (ctx->xi2*z0p2 - ctx->xi1)
        * ((c + ctx->xi2*ctx->xi2)*z0p2 - (c + ctx->xi1*ctx->xi1)) / denominator;
    lambda_m = (ctx->xi2*z0m2 - ctx->xi1)
        * ((c + ctx->xi2*ctx->xi2)*z0m2 - (c + ctx->xi1*ctx->xi1)) / denominator_reverse;
    lambda_prime_p = (ctx->xi1 - ctx->xi2)
        * ((c - ctx->xi2*ctx->xi2)*z0p2 - (c - ctx->xi1*ctx->xi1)) / denominator;
    lambda_prime_m = (ctx->xi1 - ctx->xi2)
        * ((c - ctx->xi2*ctx->xi2)*z0m2 - (c - ctx->xi1*ctx->xi1)) / denominator_reverse;

    lambda_bar_p = ((c + ctx->xi2*ctx->xi2)*z0p2*z0p2
        + ((ctx->xi1 - ctx->xi2)*(ctx->xi1 - ctx->xi2)
           - 2.0*(c + ctx->xi1*ctx->xi2))*z0p2
        + (c + ctx->xi1*ctx->xi1)) / denominator;
    lambda_bar_m = ((c + ctx->xi2*ctx->xi2)*z0m2*z0m2
        + ((ctx->xi1 - ctx->xi2)*(ctx->xi1 - ctx->xi2)
           - 2.0*(c + ctx->xi1*ctx->xi2))*z0m2
        + (c + ctx->xi1*ctx->xi1)) / denominator_reverse;
    lambda_bar_prime_p = 2.0*(ctx->xi1 - ctx->xi2)
        * (ctx->xi2*z0p2 - ctx->xi1) / denominator;
    lambda_bar_prime_m = 2.0*(ctx->xi1 - ctx->xi2)
        * (ctx->xi2*z0m2 - ctx->xi1) / denominator_reverse;

    mu1p = sqrt(ctx->z1sq - z0p2);
    mu2p = sqrt(ctx->z2sq - z0p2);
    mu1m = sqrt(ctx->z1sq - z0m2);
    mu2m = sqrt(ctx->z2sq - z0m2);
    Pi_p = grt_ellipticPi(
        (ctx->z2sq - ctx->z1sq)/(ctx->z2sq - z0p2), ctx->m_elliptic);
    Pi_m = grt_ellipticPi(
        (ctx->z2sq - ctx->z1sq)/(ctx->z2sq - z0m2), ctx->m_elliptic);

    cplx_t bracket;
    if(number == 1){
        bracket = zeta0*grt_ellipticK(ctx->m_elliptic)
            + lambda_p*Pi_p/(mu2p*mu2p)
            + lambda_m*Pi_m/(mu2m*mu2m)
            + PI*lambda_prime_p*ctx->z2/(2.0*mu1p*mu2p)
            + PI*lambda_prime_m*ctx->z2/(2.0*mu1m*mu2m);
    } else {
        bracket = eta0*grt_ellipticK(ctx->m_elliptic)
            + lambda_bar_p*Pi_p/(mu2p*mu2p)
            + lambda_bar_m*Pi_m/(mu2m*mu2m)
            - PI*lambda_bar_prime_p*ctx->z2/(2.0*mu1p*mu2p)
            - PI*lambda_bar_prime_m*ctx->z2/(2.0*mu1m*mu2m);
    }
    return -I*ctx->c_main*bracket;
}

static cplx_t basic_V(
    const int number, const cplx_t c, const LAMB2_CTX *ctx)
{
    if(number <= 2){
        if(ctx->term == LAMB2_P_TERM){
            return V_P_pair(number, c, ctx);
        }
        if(ctx->term == LAMB2_S_TERM){
            return V_S_pair(number, c, ctx);
        }
        return V_SP_pair(number, c, ctx);
    }

    real_t value;
    if(ctx->term == LAMB2_P_TERM){
        value = tail_V_P(number, ctx);
    } else if(ctx->term == LAMB2_S_TERM){
        value = tail_V_S(number, ctx);
    } else {
        value = tail_V_SP(number, ctx);
    }
    return I*value;
}

static cplx_t evaluate_partial_fraction(
    const LAMB2_PF *pf, const cplx_t roots[3], const LAMB2_CTX *ctx,
    const bool use_V)
{
    cplx_t result = 0.0;
    for(int i=0; i<3; ++i){
        if(use_V){
            result += pf->pair[i][0]*basic_V(1, roots[i], ctx)
                + pf->pair[i][1]*basic_V(2, roots[i], ctx);
        } else {
            result += pf->pair[i][0]*basic_U(1, roots[i], ctx)
                + pf->pair[i][1]*basic_U(2, roots[i], ctx);
        }
    }
    for(int i=0; i<pf->ntail; ++i){
        int number = i + 3;
        if(use_V){
            result += pf->tail[i]*basic_V(number, 0.0, ctx);
        } else {
            result += pf->tail[i]*basic_U(number, 0.0, ctx);
        }
    }
    return result;
}

static real_t evaluate_numerator(
    const LAMB2_POLY *numerator, const cplx_t roots[3],
    const LAMB2_CTX *ctx, const bool use_V)
{
    LAMB2_PF pf;
    make_partial_fraction(numerator, roots, ctx->kp2, &pf);
    return cimag(evaluate_partial_fraction(&pf, roots, ctx, use_V));
}

static void evaluate_matrix(
    const LAMB2_POLY numerator[2][3][3], const cplx_t roots[3],
    const LAMB2_CTX *ctx, const bool use_V, real_t result[3][3])
{
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            result[i][j] = evaluate_numerator(&numerator[0][i][j], roots, ctx, false)
                + evaluate_numerator(&numerator[1][i][j], roots, ctx, use_V);
        }
    }
}

static void evaluate_derivatives(
    const LAMB2_POLY numerator[3][2][3][3], const cplx_t roots[3],
    const LAMB2_CTX *ctx, const bool use_V, real_t result[3][3][3])
{
    for(int k=0; k<3; ++k){
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                result[k][i][j] = evaluate_numerator(&numerator[k][0][i][j], roots, ctx, false)
                    + evaluate_numerator(&numerator[k][1][i][j], roots, ctx, use_V);
            }
        }
    }
}

/**
 * 计算一个时间点的 F、F_(,k') 和 F_(,k)
 *
 * 这里直接对应式 (7.4.1.4)。P 项和 S 项分别组合，S 项的 V 基本积分
 * 已按式 (7.3.9) 采用 S1-H(theta-theta_c)S2 的统一表达。接收点导数
 * 的积分分子按接收点边界关系构造
 */
static void evaluate_lamb2_time(
    const real_t tbar, const LAMB2_VARS *V,
    real_t F[3][3], real_t Fk_source[3][3][3],
    real_t Fk_receiver[3][3][3])
{
    memset(F, 0, sizeof(real_t) * 3 * 3);
    memset(Fk_source, 0, sizeof(real_t) * 3 * 3 * 3);
    memset(Fk_receiver, 0, sizeof(real_t) * 3 * 3 * 3);

    if(tbar > V->k){
        LAMB2_POLY M[2][3][3];
        LAMB2_POLY dM_source[3][2][3][3];
        LAMB2_POLY dM_receiver[3][2][3][3];
        LAMB2_CTX ctx;
        make_context_P(tbar, V, &ctx);
        build_P_polynomials(tbar, V, M, dM_source, dM_receiver);
        evaluate_matrix(M, V->y, &ctx, true, F);
        evaluate_derivatives(dM_source, V->y, &ctx, true, Fk_source);
        evaluate_derivatives(dM_receiver, V->y, &ctx, true, Fk_receiver);
    }

    if(tbar > 1.0){
        LAMB2_POLY N[2][3][3];
        LAMB2_POLY dN_source[3][2][3][3];
        LAMB2_POLY dN_receiver[3][2][3][3];
        LAMB2_CTX ctx;
        make_context_S(tbar, V, &ctx);
        build_S_polynomials(tbar, V, N, dN_source, dN_receiver);
        real_t value[3][3];
        real_t dvalue_source[3][3][3];
        real_t dvalue_receiver[3][3][3];
        evaluate_matrix(N, V->yp, &ctx, true, value);
        evaluate_derivatives(dN_source, V->yp, &ctx, true, dvalue_source);
        evaluate_derivatives(dN_receiver, V->yp, &ctx, true, dvalue_receiver);
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                F[i][j] += value[i][j];
                for(int k=0; k<3; ++k){
                    Fk_source[k][i][j] += dvalue_source[k][i][j];
                    Fk_receiver[k][i][j] += dvalue_receiver[k][i][j];
                }
            }
        }
    }

    real_t t_sp = cos(V->theta - V->theta_c);
    if(V->theta > V->theta_c && tbar > t_sp && tbar < 1.0){
        LAMB2_POLY N[2][3][3];
        LAMB2_POLY dN_source[3][2][3][3];
        LAMB2_POLY dN_receiver[3][2][3][3];
        LAMB2_CTX ctx;
        make_context_SP(tbar, V, &ctx);
        build_S_polynomials(tbar, V, N, dN_source, dN_receiver);
        real_t value[3][3];
        real_t dvalue_source[3][3][3];
        real_t dvalue_receiver[3][3][3];
        evaluate_matrix(N, V->yp, &ctx, true, value);
        evaluate_derivatives(dN_source, V->yp, &ctx, true, dvalue_source);
        evaluate_derivatives(dN_receiver, V->yp, &ctx, true, dvalue_receiver);
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                F[i][j] -= value[i][j];
                for(int k=0; k<3; ++k){
                    Fk_source[k][i][j] -= dvalue_source[k][i][j];
                    Fk_receiver[k][i][j] -= dvalue_receiver[k][i][j];
                }
            }
        }
    }
}

static real_t derivative_three_points(
    const real_t x0, const real_t x1, const real_t x2,
    const real_t f0, const real_t f1, const real_t f2, const real_t x)
{
    real_t result = f0 * (2.0*x - x1 - x2) / ((x0 - x1)*(x0 - x2));
    result += f1 * (2.0*x - x0 - x2) / ((x1 - x0)*(x1 - x2));
    result += f2 * (2.0*x - x0 - x1) / ((x2 - x0)*(x2 - x1));
    return result;
}

static void differentiate_Fk(
    const real_t *ts, const int nt, const real_t (*Fk)[3][3][3],
    real_t (*dG)[3][3][3])
{
    if(nt == 1){
        memset(dG, 0, sizeof(real_t) * 3 * 3 * 3);
        return;
    }

    for(int k=0; k<3; ++k){
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                if(nt == 2){
                    real_t value = (Fk[1][k][i][j] - Fk[0][k][i][j]) / (ts[1] - ts[0]);
                    dG[0][k][i][j] = value;
                    dG[1][k][i][j] = value;
                } else {
                    dG[0][k][i][j] = derivative_three_points(
                        ts[0], ts[1], ts[2], Fk[0][k][i][j], Fk[1][k][i][j], Fk[2][k][i][j], ts[0]);
                    for(int n=1; n<nt-1; ++n){
                        dG[n][k][i][j] = derivative_three_points(
                            ts[n-1], ts[n], ts[n+1], Fk[n-1][k][i][j], Fk[n][k][i][j], Fk[n+1][k][i][j], ts[n]);
                    }
                    dG[nt-1][k][i][j] = derivative_three_points(
                        ts[nt-3], ts[nt-2], ts[nt-1], Fk[nt-3][k][i][j],
                        Fk[nt-2][k][i][j], Fk[nt-1][k][i][j], ts[nt-1]);
                }
            }
        }
    }
}

static void print_lamb2_header(void)
{
    char *name = NULL;
    printf("#%13s", "tbar");
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            GRT_SAFE_ASPRINTF(&name, "G%d%d", i+1, j+1);
            printf("%14s", name);
        }
    }
    for(int k=0; k<3; ++k){
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                GRT_SAFE_ASPRINTF(&name, "G%d%d,%d'", i+1, j+1, k+1);
                printf("%14s", name);
            }
        }
    }
    for(int k=0; k<3; ++k){
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                GRT_SAFE_ASPRINTF(&name, "G%d%d,%d", i+1, j+1, k+1);
                printf("%14s", name);
            }
        }
    }
    GRT_SAFE_FREE_PTR(name);
    printf("\n");
}

static void print_lamb2_derivatives(const real_t dG[3][3][3])
{
    for(int k=0; k<3; ++k){
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                printf("%14.6e", dG[k][i][j]);
            }
        }
    }
}

static void solve_lamb2(
    const real_t nu, const real_t *ts, const int nt,
    const real_t theta, const real_t azimuth,
    real_t (*G)[3][3], real_t (*dG_source)[3][3][3],
    real_t (*dG_receiver)[3][3][3])
{
    if(nu <= 0.0 || nu >= 0.5){
        GRTRaiseError("possion ratio (%lf) is out of bound.", nu);
    }
    if(ts == NULL || nt <= 0){
        GRTRaiseError("The time series for lamb2 should not be empty.\n");
    }
    if(theta <= 0.0 || theta >= 90.0){
        GRTRaiseError("theta should be in (0, 90) degree for lamb2.\n");
    }
    if(azimuth < 0.0 || azimuth > 360.0){
        GRTRaiseError("azimuth should be in [0, 360] degree for lamb2.\n");
    }
    for(int i=0; i<nt; ++i){
        if(ts[i] < 0.0){
            GRTRaiseError("The time series for lamb2 should be nonnegative.\n");
        }
        if(i > 0 && ts[i] <= ts[i-1]){
            GRTRaiseError("The time series for lamb2 should be strictly increasing.\n");
        }
    }

    LAMB2_VARS V = {0};
    V.nu = nu;
    V.k2 = 0.5 * (1.0 - 2.0*nu) / (1.0 - nu);
    V.k = sqrt(V.k2);
    V.kp2 = 1.0 - V.k2;
    V.kp = sqrt(V.kp2);
    V.theta = theta * DEG1;
    V.phi = azimuth * DEG1;
    V.st = sin(V.theta);
    V.ct = cos(V.theta);
    V.sf = sin(V.phi);
    V.cf = cos(V.phi);
    V.theta_c = asin(V.k);
    grt_rayleigh1_roots(V.nu, V.y);
    for(int i=0; i<3; ++i){
        V.yp[i] = V.y[i] - V.kp2;
    }

    real_t (*F)[3][3] = calloc((size_t)nt, sizeof(*F));
    real_t (*Fk_source)[3][3][3] = calloc((size_t)nt, sizeof(*Fk_source));
    real_t (*Fk_receiver)[3][3][3] = calloc((size_t)nt, sizeof(*Fk_receiver));
    real_t (*dG_source_tmp)[3][3][3] = calloc((size_t)nt, sizeof(*dG_source_tmp));
    real_t (*dG_receiver_tmp)[3][3][3] = calloc((size_t)nt, sizeof(*dG_receiver_tmp));
    if(F == NULL || Fk_source == NULL || Fk_receiver == NULL
        || dG_source_tmp == NULL || dG_receiver_tmp == NULL){
        GRT_SAFE_FREE_PTR(F);
        GRT_SAFE_FREE_PTR(Fk_source);
        GRT_SAFE_FREE_PTR(Fk_receiver);
        GRT_SAFE_FREE_PTR(dG_source_tmp);
        GRT_SAFE_FREE_PTR(dG_receiver_tmp);
        GRTRaiseError("Cannot allocate lamb2 output arrays.\n");
    }

    for(int i=0; i<nt; ++i){
        evaluate_lamb2_time(ts[i], &V, F[i], Fk_source[i], Fk_receiver[i]);
    }
    differentiate_Fk(ts, nt, Fk_source, dG_source_tmp);
    differentiate_Fk(ts, nt, Fk_receiver, dG_receiver_tmp);

    bool isprint = (G == NULL && dG_source == NULL && dG_receiver == NULL);
    if(isprint){
        print_lamb2_header();
    }
    for(int n=0; n<nt; ++n){
        if(G != NULL){
            memcpy(G[n], F[n], sizeof(real_t) * 3 * 3);
        }
        if(dG_source != NULL){
            memcpy(dG_source[n], dG_source_tmp[n], sizeof(real_t) * 3 * 3 * 3);
        }
        if(dG_receiver != NULL){
            memcpy(dG_receiver[n], dG_receiver_tmp[n], sizeof(real_t) * 3 * 3 * 3);
        }
        if(isprint){
            printf("%14.6e", ts[n]);
            for(int i=0; i<3; ++i){
                for(int j=0; j<3; ++j){
                    printf("%14.6e", F[n][i][j]);
                }
            }
            print_lamb2_derivatives(dG_source_tmp[n]);
            print_lamb2_derivatives(dG_receiver_tmp[n]);
            printf("\n");
        }
    }

    GRT_SAFE_FREE_PTR(F);
    GRT_SAFE_FREE_PTR(Fk_source);
    GRT_SAFE_FREE_PTR(Fk_receiver);
    GRT_SAFE_FREE_PTR(dG_source_tmp);
    GRT_SAFE_FREE_PTR(dG_receiver_tmp);
}

void grt_solve_lamb2(
    const real_t nu, const real_t *ts, const int nt,
    const real_t theta, const real_t azimuth,
    real_t (*G)[3][3], real_t (*dG_source)[3][3][3],
    real_t (*dG_receiver)[3][3][3])
{
    solve_lamb2(nu, ts, nt, theta, azimuth, G, dG_source, dG_receiver);
}
