/**
 * @file   lamb3.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-09
 *
 *    使用广义闭合解求解第三类 Lamb 问题，参考：
 *
 *        张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社
 *
 *    反射项和转换项的多项式系数在时间循环前预计算，并按到时跳过
 *    时间序列不会用到的 P/S/PS/SP 项；运行时只组合时间幂次与基本积分。
 *    源点和接收点导数均按第 8.1 节的矩阵关系实现，第 7 章中的基本积分
 *    直接复用 lamb_basic.c 的实现
 */

#include <float.h>

#include "grt/lamb/elliptic.h"
#include "grt/lamb/lamb3.h"
#include "grt/lamb/lamb_basic.h"
#include "grt/lamb/lamb_util.h"
#include "grt/lamb/qromb.h"

// 预先计算好的多项式系数
#include "lamb3_coeffs.c_"

#define LAMB3_TAIL_SIZE 7
#define LAMB3_QROMB_EPS 1e-6
#define LAMB3_RATIO_EPS 1e-8
/* 反射射线角过小时，闭合公式中的变量变换开始失去有效数字 */
#define LAMB3_SMALL_R_WARNING_RATIO 1e-2

typedef LAMB_POLY LAMB3_POLY;

/** 一个有理分式的部分分式系数 */
typedef struct {
    cplx_t pair[3][2];            ///< pair[i][0/1]，i=0,1,2 为根索引，0/1 为 x/常数系数
    cplx_t pole[4];               ///< pole[i]，i=0,1,2,3 分别为 x^(-1)、x^(-2)、x^(-3)、x^(-4) 的主部系数
    cplx_t tail[LAMB3_TAIL_SIZE]; ///< tail[m]，m 为多项式商的 x 次数
    int npole;                    ///< x=0 处主部的有效项数
    int ntail;                    ///< 多项式商的有效项数
} LAMB3_PF;

/** 按时间幂次保存的部分分式系数 */
typedef struct {
    cplx_t pair[3][2][LAMB3_TIME_SIZE];            ///< pair[i][0/1][r]，i 为根索引，0/1 为 x/常数项，r 为 tbar 次数
    cplx_t pole[4][LAMB3_TIME_SIZE];               ///< pole[i][r]，i 为负幂阶数减一，r 为 tbar 次数
    cplx_t tail[LAMB3_TAIL_SIZE][LAMB3_TIME_SIZE]; ///< tail[m][r]，m 为 x 次数，r 为 tbar 次数
    int npole;                                     ///< x=0 处主部的有效项数
    int ntail;                                     ///< 多项式商的有效项数
    int time_degree;                               ///< tbar 的最高次数
} LAMB3_PF_COEFF;

/** 反射项的一组部分分式系数 */
typedef struct {
    LAMB3_PF_COEFF M[2][3][3];                 ///< M[xi][i][j]，xi=0/1 对应 U/V，i,j=0,1,2 分别为接收点和源点分量
    LAMB3_PF_COEFF dM[3][2][3][3];             ///< dM[k][xi][i][j]，k=0,1,2 为源点坐标导数方向，xi=0/1 对应 U/V，i,j 为矩阵分量
} LAMB3_REFLECTION_PF_SET;

/** 转换项的一组部分分式系数 */
typedef struct {
    LAMB3_PF_COEFF M[3][3];                 ///< M[i][j]，i,j=0,1,2 分别为接收点和源点分量
    LAMB3_PF_COEFF dM[3][3][3];             ///< dM[k][i][j]，k=0,1,2 为源点坐标导数方向，i,j 为矩阵分量
    LAMB3_PF_COEFF receiver_vertical[3][3]; ///< receiver_vertical[i][j]，i,j 为矩阵分量
} LAMB3_CONVERSION_PF_SET;

/** 一个几何排列对应的全部部分分式系数 */
typedef struct {
    LAMB3_REFLECTION_PF_SET P;  ///< P 波反射项的部分分式结果
    LAMB3_REFLECTION_PF_SET S;  ///< S 波反射项的部分分式结果
    LAMB3_CONVERSION_PF_SET PS; ///< P 入射、S 转换项的部分分式结果
    LAMB3_CONVERSION_PF_SET SP; ///< S 入射、P 转换项的部分分式结果
} LAMB3_PF_COEFFICIENTS;

/** 式 (8.3.12)-(8.3.17) 中的 PS 基本积分参数 */
typedef struct {
    real_t xi1;        ///< 四次式分解参数 xi1
    real_t xi2;        ///< 四次式分解参数 xi2
    real_t z1sq;       ///< z1^2
    real_t z2sq;       ///< z2^2
    real_t z1;         ///< z1=sqrt(z1^2)
    real_t z2;         ///< z2=sqrt(z2^2)
    real_t m;          ///< 椭圆积分参数 m=z2^2/z1^2
    real_t cps;        ///< PS 基本积分的公共系数
    real_t asym_delta; ///< xi1、xi2 接近时的渐近展开参数
} LAMB3_PS_CTX;

/** 一个时间点的 PS 基本积分 */
typedef struct {
    cplx_t pair[3][2]; ///< pair[i][0/1]，i=0,1,2 为根索引，0/1 为二次因子的两个基本积分
    real_t V[11];      ///< V[j+4] 保存阶数 j=-4,...,6 的基本积分
} LAMB3_PS_BASIS;

/** 一个时间点和一组反射积分对应的基本积分 */
typedef struct {
    cplx_t pair[3][2];            ///< pair[i][0/1]，i=0,1,2 为根索引，0/1 为二次因子的两个基本积分
    cplx_t tail[LAMB3_TAIL_SIZE]; ///< tail[m]，m 为 x 次数
} LAMB3_REFLECTION_BASIS;

/** 一次时间点评估得到的各个 F 项 */
typedef struct {
    real_t F[3][3];              ///< F[i][j] 是接收点分量 i 和源点分量 j 的矩阵
    real_t Fk_source[3][3][3];   ///< Fk_source[k'][i][j] 是源点坐标导数
    real_t Fk_receiver[3][3][3]; ///< Fk_receiver[k][i][j] 是接收点坐标导数
} LAMB3_F;

/* phi -> phi + pi 后各矩阵分量的符号，索引为 [接收点分量][源点分量] */
static const int LAMB3_PHI_PI_SIGN[3][3] = {
    {1, 1, -1},
    {1, 1, -1},
    {-1, -1, 1},
};

static LAMB3_POLY make_denominator(const cplx_t roots[3], const int pole_order, const bool plus, const cplx_t scale) {
    LAMB3_POLY denominator = grt_lamb_poly_const(scale);
    LAMB3_POLY x = grt_lamb_poly_x();
    LAMB3_POLY xpow = grt_lamb_poly_const(1.0);
    for (int i = 0; i < pole_order; ++i) {
        xpow = grt_lamb_poly_mul(xpow, x);
    }
    denominator = grt_lamb_poly_mul(denominator, xpow);
    for (int i = 0; i < 3; ++i) {
        denominator = grt_lamb_poly_mul(denominator, grt_lamb_poly_factor(roots[i], plus));
    }
    return denominator;
}

/**
 * 对 x^p 乘三个二次因子的分母作部分分式分解
 *
 * plus=true  对应 (x^2+y_i)
 * plus=false 对应 (x^2-y_i)
 * pole[i] 对应 x^(-(i+1))，tail[m] 对应 x^m
 */
static void make_partial_fraction(const LAMB3_POLY *numerator, const cplx_t roots[3], const int pole_order, const bool plus, const cplx_t scale,
                                  LAMB3_PF *pf) {
    LAMB3_POLY denominator = make_denominator(roots, pole_order, plus, scale);
    LAMB3_POLY quotient;
    LAMB3_POLY remainder;
    grt_lamb_poly_divide(*numerator, denominator, &quotient, &remainder);

    memset(pf, 0, sizeof(*pf));
    pf->npole = pole_order;
    pf->ntail = quotient.degree + 1;
    if (pf->npole > 4 || pf->ntail > LAMB3_TAIL_SIZE) {
        GRTRaiseError("The partial-fraction degree is too large in lamb3.\n");
    }
    for (int i = 0; i < pf->ntail; ++i) {
        pf->tail[i] = quotient.c[i];
    }

    /* x=0 处的主部系数 */
    LAMB3_POLY even_denominator = grt_lamb_poly_const(scale);
    for (int i = 0; i < 3; ++i) {
        even_denominator = grt_lamb_poly_mul(even_denominator, grt_lamb_poly_factor(roots[i], plus));
    }
    cplx_t series[LAMB3_POLY_SIZE] = {0};
    for (int n = 0; n < pole_order; ++n) {
        cplx_t value = n <= remainder.degree ? remainder.c[n] : 0.0;
        for (int j = 1; j <= n; ++j) {
            value -= even_denominator.c[j] * series[n - j];
        }
        series[n] = value / even_denominator.c[0];
        pf->pole[pole_order - n - 1] = series[n];
    }

    /* 每个二次因子的两个简单极点 */
    for (int i = 0; i < 3; ++i) {
        LAMB3_POLY quotient_factor = grt_lamb_poly_const(scale);
        LAMB3_POLY xpow = grt_lamb_poly_const(1.0);
        LAMB3_POLY x = grt_lamb_poly_x();
        for (int j = 0; j < pole_order; ++j) {
            xpow = grt_lamb_poly_mul(xpow, x);
        }
        quotient_factor = grt_lamb_poly_mul(quotient_factor, xpow);
        for (int j = 0; j < 3; ++j) {
            if (i != j) {
                quotient_factor = grt_lamb_poly_mul(quotient_factor, grt_lamb_poly_factor(roots[j], plus));
            }
        }
        cplx_t c = plus ? sqrt(-roots[i]) : sqrt(roots[i]);
        if (fabs(c) < 1e-14) {
            GRTRaiseError("A zero partial-fraction root is not supported in lamb3.\n");
        }
        cplx_t value_plus = grt_lamb_poly_eval(&remainder, c) / grt_lamb_poly_eval(&quotient_factor, c);
        cplx_t value_minus = grt_lamb_poly_eval(&remainder, -c) / grt_lamb_poly_eval(&quotient_factor, -c);
        pf->pair[i][0] = (value_plus - value_minus) / (2.0 * c);
        pf->pair[i][1] = 0.5 * (value_plus + value_minus);
    }
}

static void make_partial_fraction_coeff(const LAMB3_POLY_COEFF *numerator, const cplx_t roots[3], const int pole_order, const bool plus,
                                        const cplx_t scale, LAMB3_PF_COEFF *result) {
    memset(result, 0, sizeof(*result));
    result->npole = pole_order;
    result->time_degree = numerator->time_degree;
    for (int r = 0; r <= numerator->time_degree; ++r) {
        LAMB3_POLY polynomial = {0};
        polynomial.degree = numerator->degree;
        for (int m = 0; m <= numerator->degree; ++m) {
            polynomial.c[m] = numerator->c[r][m];
        }
        LAMB3_PF pf;
        make_partial_fraction(&polynomial, roots, pole_order, plus, scale, &pf);
        result->ntail = GRT_MAX(result->ntail, pf.ntail);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 2; ++j) {
                result->pair[i][j][r] = pf.pair[i][j];
            }
        }
        for (int i = 0; i < pf.npole; ++i) {
            result->pole[i][r] = pf.pole[i];
        }
        for (int i = 0; i < pf.ntail; ++i) {
            result->tail[i][r] = pf.tail[i];
        }
    }
}

static void make_partial_fraction_coeff_x(const LAMB3_POLY_COEFF *numerator, const cplx_t roots[3], const int pole_order, const bool plus,
                                          const cplx_t scale_x, const cplx_t scale_denominator, LAMB3_PF_COEFF *result) {
    LAMB3_POLY_COEFF shifted = {0};
    shifted.degree = numerator->degree + 1;
    shifted.time_degree = numerator->time_degree;
    for (int r = 0; r <= numerator->time_degree; ++r) {
        for (int m = 0; m <= numerator->degree; ++m) {
            shifted.c[r][m + 1] = scale_x * numerator->c[r][m];
        }
    }
    make_partial_fraction_coeff(&shifted, roots, pole_order, plus, scale_denominator, result);
}

static void make_reflection_pf_set(const LAMB3_COEFF_SET *coeffs, const cplx_t roots[3], const real_t kp2, LAMB3_REFLECTION_PF_SET *result) {
    memset(result, 0, sizeof(*result));
    for (int xi = 0; xi < 2; ++xi) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                make_partial_fraction_coeff(&coeffs->M[xi][i][j], roots, 0, true, 16.0 * kp2, &result->M[xi][i][j]);
                for (int k = 0; k < 2; ++k) {
                    make_partial_fraction_coeff(&coeffs->dM[k][xi][i][j], roots, 0, true, 16.0 * kp2, &result->dM[k][xi][i][j]);
                }
                make_partial_fraction_coeff_x(&coeffs->M[xi][i][j], roots, 0, true, -1.0, 16.0 * kp2, &result->dM[2][xi][i][j]);
            }
        }
    }
}

static void make_conversion_pf_set(const LAMB3_CONVERSION_COEFF_SET *coeffs, const cplx_t roots[3], LAMB3_CONVERSION_PF_SET *result) {
    memset(result, 0, sizeof(*result));
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const bool reuse_M = i == 1 && j == 0;
            if (!reuse_M) {
                make_partial_fraction_coeff(&coeffs->M[i][j], roots, 3, false, 1.0, &result->M[i][j]);
            }
            for (int k = 0; k < 2; ++k) {
                const bool reuse_derivative = (i == 1 && j == 0) || (k == 1 && (i == 0 || (i == 2 && j == 0)));
                if (!reuse_derivative) {
                    make_partial_fraction_coeff(&coeffs->dM[k][i][j], roots, 4, false, 1.0, &result->dM[k][i][j]);
                }
            }
            /* 源点竖直导数使用系数阶段已生成的 dM[2] */
            if (!reuse_M) {
                make_partial_fraction_coeff(&coeffs->dM[2][i][j], roots, 4, false, 1.0, &result->dM[2][i][j]);
            }
        }
    }

    /* 按式 (8.3.1) 复用矩阵对称性和两个水平导数方向的关系 */
    result->M[1][0] = result->M[0][1];
    result->dM[0][1][0] = result->dM[0][0][1];
    result->dM[1][0][0] = result->dM[0][0][1];
    result->dM[1][0][1] = result->dM[0][1][1];
    result->dM[1][0][2] = result->dM[0][1][2];
    result->dM[1][1][0] = result->dM[1][0][1];
    result->dM[1][2][0] = result->dM[0][2][1];
    result->dM[2][1][0] = result->dM[2][0][1];
}

static void make_conversion_receiver_pf_set(const LAMB3_CONVERSION_COEFF_SET *partner, const cplx_t roots[3], LAMB3_CONVERSION_PF_SET *result) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            /* M_,3^R = [M_,3'^partner(phi+pi)^T] */
            make_partial_fraction_coeff(&partner->dM[2][j][i], roots, 4, false, LAMB3_PHI_PI_SIGN[j][i], &result->receiver_vertical[i][j]);
        }
    }
}

static void make_coefficients(const LAMB3_VARS *V, const bool need_P, const bool need_S, const bool need_PS, const bool need_SP,
                              LAMB3_PF_COEFFICIENTS *coefficients) {
    /* 原始和部分分式系数体积较大，统一放在堆上，并按到时跳过未用到的项 */
    if (need_P || need_S) {
        LAMB3_COEFF_SET *reflection_raw = calloc(1, sizeof(*reflection_raw));
        if (need_P) {
            make_lamb3_P_coefficients(V, reflection_raw);
            make_reflection_pf_set(reflection_raw, V->rayleigh, V->kp2, &coefficients->P);
        }
        if (need_S) {
            make_lamb3_S_coefficients(V, reflection_raw);
            make_reflection_pf_set(reflection_raw, V->rayleigh_shifted, V->kp2, &coefficients->S);
        }
        GRT_SAFE_FREE_PTR(reflection_raw);
    }

    if (need_PS || need_SP) {
        LAMB3_CONVERSION_COEFF_SET *conversion_raw = calloc(1, sizeof(*conversion_raw));
        /* 接收点竖向导数对应交换深度后的 PS/SP 项 */
        LAMB3_VARS reciprocal = *V;
        reciprocal.source_depth = V->receiver_depth;
        reciprocal.receiver_depth = V->source_depth;
        if (need_PS) {
            make_lamb3_PS_coefficients(V, 1.0, 1, conversion_raw);
            make_conversion_pf_set(conversion_raw, V->rayleigh_ps, &coefficients->PS);
            make_lamb3_SP_coefficients(&reciprocal, conversion_raw);
            make_conversion_receiver_pf_set(conversion_raw, V->rayleigh_ps, &coefficients->PS);
        }
        if (need_SP) {
            make_lamb3_SP_coefficients(V, conversion_raw);
            make_conversion_pf_set(conversion_raw, V->rayleigh_ps, &coefficients->SP);
            make_lamb3_PS_coefficients(&reciprocal, 1.0, 1, conversion_raw);
            make_conversion_receiver_pf_set(conversion_raw, V->rayleigh_ps, &coefficients->SP);
        }
        GRT_SAFE_FREE_PTR(conversion_raw);
    }
}

/** 计算实系数四次方程的四个复根 */
static void quartic_roots(const real_t coefficient[5], cplx_t roots[4]) {
    real_t leading = coefficient[0];
    if (leading == 0.0) {
        GRTRaiseError("The leading coefficient of the PS quartic is zero.\n");
    }

    real_t a = coefficient[1] / leading;
    real_t b = coefficient[2] / leading;
    real_t c = coefficient[3] / leading;
    real_t d = coefficient[4] / leading;
    real_t p = b - 3.0 * a * a / 8.0;
    real_t q = a * a * a / 8.0 - a * b / 2.0 + c;
    real_t r = -3.0 * pow(a, 4) / 256.0 + a * a * b / 16.0 - a * c / 4.0 + d;
    real_t shift = a / 4.0;

    if (fabs(q) < 1e-14 * (1.0 + fabs(p) + fabs(r))) {
        cplx_t discriminant = p * p - 4.0 * r;
        cplx_t u = 0.5 * (-p + sqrt(discriminant));
        cplx_t v = 0.5 * (-p - sqrt(discriminant));
        roots[0] = sqrt(u) - shift;
        roots[1] = -sqrt(u) - shift;
        roots[2] = sqrt(v) - shift;
        roots[3] = -sqrt(v) - shift;
        return;
    }

    cplx_t resolvent[3];
    grt_lamb_cubic_roots(2.0 * p, p * p - 4.0 * r, -q * q, resolvent);
    cplx_t U = 0.0;
    real_t best = -DBL_MAX;
    for (int i = 0; i < 3; ++i) {
        if (fabs(cimag(resolvent[i])) < 1e-9 * (1.0 + fabs(creal(resolvent[i]))) && creal(resolvent[i]) > best) {
            best = creal(resolvent[i]);
            U = creal(resolvent[i]);
        }
    }
    if (best <= 0.0) {
        for (int i = 0; i < 3; ++i) {
            if (creal(resolvent[i]) > best) {
                best = creal(resolvent[i]);
                U = resolvent[i];
            }
        }
    }
    cplx_t u = sqrt(U);
    if (fabs(u) < 1e-14) {
        GRTRaiseError("The PS quartic resolvent is singular in lamb3.\n");
    }
    cplx_t v = 0.5 * (p + U - q / u);
    cplx_t w = 0.5 * (p + U + q / u);
    roots[0] = 0.5 * (-u + sqrt(u * u - 4.0 * v)) - shift;
    roots[1] = 0.5 * (-u - sqrt(u * u - 4.0 * v)) - shift;
    roots[2] = 0.5 * (u + sqrt(u * u - 4.0 * w)) - shift;
    roots[3] = 0.5 * (u - sqrt(u * u - 4.0 * w)) - shift;
}

static real_t snell_residual(const real_t horizontal_split, const real_t R, const real_t p_depth, const real_t s_depth, const real_t k) {
    real_t rp = hypot(horizontal_split, p_depth);
    real_t rs = hypot(R - horizontal_split, s_depth);
    return (R - horizontal_split) / rs - k * horizontal_split / rp;
}

/** 按式 (4.5.20)-(4.5.21) 求 PS 或 SP 的归一化到时 */
static real_t conversion_arrival(const real_t R, const real_t p_depth, const real_t s_depth, const real_t k, const real_t direct_distance) {
    real_t left = GRT_MAX(1e-12 * (1.0 + R), DBL_MIN);
    real_t right = R - left;
    real_t fleft = snell_residual(left, R, p_depth, s_depth, k);
    real_t fright = snell_residual(right, R, p_depth, s_depth, k);
    if (fleft * fright > 0.0) {
        GRTRaiseError("Cannot bracket the conversion-wave arrival in lamb3.\n");
    }
    for (int i = 0; i < 100; ++i) {
        real_t middle = 0.5 * (left + right);
        real_t fmiddle = snell_residual(middle, R, p_depth, s_depth, k);
        if (fabs(fmiddle) < 1e-14 || right - left < 1e-13 * (1.0 + R)) {
            left = right = middle;
            break;
        }
        if (fleft * fmiddle > 0.0) {
            left = middle;
            fleft = fmiddle;
        } else {
            right = middle;
            fright = fmiddle;
        }
    }
    real_t split = 0.5 * (left + right);
    real_t rp = hypot(split, p_depth);
    real_t rs = hypot(R - split, s_depth);
    return (k * rp + rs) / direct_distance;
}

static void make_vars(const real_t nu, const real_t R, const real_t source_depth, const real_t receiver_depth, const real_t azimuth, LAMB3_VARS *V) {
    memset(V, 0, sizeof(*V));
    V->nu = nu;
    V->k2 = 0.5 * (1.0 - 2.0 * nu) / (1.0 - nu);
    V->k = sqrt(V->k2);
    V->kp2 = 1.0 - V->k2;
    V->kp = sqrt(V->kp2);
    V->R = R;
    V->source_depth = source_depth;
    V->receiver_depth = receiver_depth;
    V->r = hypot(R, source_depth - receiver_depth);
    V->rp = hypot(R, source_depth + receiver_depth);
    V->theta = atan2(R, source_depth - receiver_depth);
    V->theta_ref = atan2(R, source_depth + receiver_depth);
    V->phi = azimuth * DEG1;
    V->st = sin(V->theta);
    V->ct = cos(V->theta);
    V->st_ref = sin(V->theta_ref);
    V->ct_ref = cos(V->theta_ref);
    V->sf = sin(V->phi);
    V->cf = cos(V->phi);
    V->theta_c = asin(V->k);
    V->varsigma = V->r / V->rp;
    V->r_over_R = V->r / V->R;
    V->kap1 = 1.0 / V->kp;
    V->kap2 = 1.0 / V->kp2;
    V->reflection_P = V->k * V->rp / V->r;
    V->reflection_S = V->rp / V->r;
    V->t_sps = cos(V->theta_ref - V->theta_c) * V->reflection_S;
    V->conversion_scale = V->kp2 * V->r / (16.0 * V->R);
    V->conversion_dscale = V->kp * V->kp2 * V->r / (32.0 * V->R);
    V->use_angle_ratio = fabs(V->cf) > LAMB3_RATIO_EPS;
    V->angle_ratio = V->use_angle_ratio ? V->sf / V->cf : 0.0;
    V->supercritical = V->theta_ref > V->theta_c;

    if (source_depth > 0.0 && receiver_depth > 0.0) {
        V->tps = conversion_arrival(R, source_depth, receiver_depth, V->k, V->r);
        V->tsp = conversion_arrival(R, receiver_depth, source_depth, V->k, V->r);
    } else {
        V->tps = INFINITY;
        V->tsp = INFINITY;
    }

    grt_rayleigh1_roots(nu, V->rayleigh);
    for (int i = 0; i < 3; ++i) {
        V->rayleigh_shifted[i] = V->rayleigh[i] - V->kp2;
    }
    grt_lamb_cubic_roots(8.0 * nu * nu - 8.0 * nu + 3.0, 8.0 * nu - 5.0, 1.0, V->rayleigh_ps);
}



/** 式 (8.4.2) 中 PP/SS1 的 V9 数积，积分变量为 x=m+i n cos(phi) */
typedef struct {
    real_t m;    ///< 路径起点 m=sbar*cos(theta_ref)
    real_t n;    ///< 虚轴半长 n
    real_t kp2;  ///< k'^2，PP 取 x^2+k'^2，SS1 取 x^2-k'^2
    int    sign; ///< +1 对应 PP，-1 对应 SS1
} LAMB3_V9_PATH_CTX;

/** 式 (8.4.2.2) 中 SS2 的 V9 数积，x 从 m 积到 k' */
typedef struct {
    real_t m;  ///< 路径起点 m
    real_t n;  ///< Q_S=(x-m)^2+n^2 中的 n
    real_t kp; ///< 路径终点 k'
} LAMB3_V9_SS2_CTX;

static cplx_t positive_branch(const cplx_t value) {
    cplx_t result = sqrt(value);
    if (creal(result) < 0.0) {
        result = -result;
    }
    return result;
}

/** PP 看 |xi2|^6，SS 看 |xi1|^6，与 8.4.2.1-8.4.2.2 的判据一致 */
static bool use_reflection_V9_quadrature(const LAMB_BASIC_CONTEXT *ctx) {
    if (ctx->term == LAMB_BASIC_P_TERM) {
        return pow(fabs(ctx->xi2), 6.0) > 1e8;
    }
    if (ctx->term == LAMB_BASIC_S_TERM) {
        return pow(fabs(ctx->xi1), 6.0) > 1e8;
    }
    return false;
}

/** PP/SS1：x=m+i n cos(phi)，d x / sqrt(Q) = i d phi */
static real_t reflection_V9_path_integrand(real_t phi, void *userdata) {
    const LAMB3_V9_PATH_CTX *C = userdata;
    cplx_t x = C->m + I * C->n * cos(phi);
    cplx_t denom = positive_branch(x * x + C->sign * C->kp2);
    if (cabs(denom) < 1e-30) {
        return 0.0;
    }
    cplx_t x2 = x * x;
    return creal(x2 * x2 * x2 / denom);
}

/** SS2：x=m+(k'-m) sin^2(phi)，去掉 x=k' 处 sqrt(k'^2-x^2) 的端点奇点 */
static real_t reflection_V9_SS2_integrand(real_t phi, void *userdata) {
    const LAMB3_V9_SS2_CTX *C = userdata;
    real_t sine = sin(phi);
    real_t x = C->m + (C->kp - C->m) * sine * sine;
    real_t QS = (x - C->m) * (x - C->m) + C->n * C->n;
    real_t x2 = x * x;
    return x2 * x2 * x2 / sqrt(QS) * 2.0 * sqrt(C->kp - C->m) * sine / sqrt(x + C->kp);
}

/**
 * 按 8.4.2 退回到引入 xi 之前的 V9 数积
 * PP 只积虚轴路径，SS 再按 (7.3.9) 组合 V9^{SS1}-H V9^{SS2}
 */
static cplx_t numerical_reflection_V9(const LAMB_BASIC_CONTEXT *ctx) {
    LAMB3_V9_PATH_CTX path = {
        .m = ctx->m,
        .n = ctx->n,
        .kp2 = ctx->kp2,
        .sign = ctx->term == LAMB_BASIC_P_TERM ? 1 : -1,
    };
    real_t value = grt_lamb_qromb(reflection_V9_path_integrand, 0.0, HALFPI, LAMB3_QROMB_EPS, &path);

    /* SS2 仅在 m<k' 时非零，此时 1<sbar<k'/cos(theta_ref) 且 theta_ref>theta_c */
    if (ctx->term == LAMB_BASIC_S_TERM) {
        real_t kp = sqrt(ctx->kp2);
        if (kp > ctx->m) {
            LAMB3_V9_SS2_CTX SS2 = {.m = ctx->m, .n = ctx->n, .kp = kp};
            value += grt_lamb_qromb(reflection_V9_SS2_integrand, 0.0, HALFPI, LAMB3_QROMB_EPS, &SS2);
        }
    }
    return I * value;
}

static cplx_t basic_U(const int number, const LAMB_BASIC_CONTEXT *ctx) {
    if (number <= 6) {
        return grt_lamb_basic_U(number, 0.0, ctx);
    }
    if (ctx->term == LAMB_BASIC_SP_TERM) {
        return 0.0;
    }

    int exponent = number - 3;
    real_t m = ctx->m;
    real_t n = ctx->n;
    real_t value = 0.0;
    for (int even_power = 0; even_power <= exponent; even_power += 2) {
        int half_power = even_power / 2;
        real_t moment = PI / 2.0 * tgamma((real_t)half_power + 0.5) / (sqrt(PI) * tgamma((real_t)half_power + 1.0));
        value += tgamma((real_t)exponent + 1.0) / (tgamma((real_t)even_power + 1.0) * tgamma((real_t)(exponent - even_power) + 1.0)) *
                 pow(m, exponent - even_power) * pow(n, even_power) * (half_power % 2 == 0 ? 1.0 : -1.0) * moment;
    }
    return I * value;
}

static void calculate_H7(const LAMB_BASIC_CONTEXT *ctx, real_t H[7]) {
    grt_lamb_calculate_H(ctx, grt_ellipticK(ctx->m_elliptic), H);
    real_t a;
    if (ctx->term == LAMB_BASIC_P_TERM) {
        a = 1.0 / ctx->z2sq;
    } else if (ctx->term == LAMB_BASIC_S_TERM) {
        a = -(ctx->z2sq + 1.0) / ctx->z2sq;
    } else {
        a = -ctx->c1 / ctx->c2;
    }
    real_t m = ctx->m_elliptic;
    real_t gamma1 = 3.0 * m * a * a + 2.0 * a * (m + 1.0) + 1.0;
    real_t gamma2 = m + 1.0 + 3.0 * m * a;
    real_t gamma3 = a * (a + 1.0) * (a * m + 1.0);
    H[5] = (7.0 * gamma1 * H[4] - 6.0 * gamma2 * H[3] + 5.0 * m * H[2]) / (8.0 * gamma3);
    H[6] = (9.0 * gamma1 * H[5] - 8.0 * gamma2 * H[4] + 7.0 * m * H[3]) / (10.0 * gamma3);
}

/** 计算反射项的 V8、V9 基本积分 */
static real_t reflection_V_high(const int number, const LAMB_BASIC_CONTEXT *ctx) {
    real_t H[7];
    calculate_H7(ctx, H);
    real_t xi1 = ctx->xi1;
    real_t xi2 = ctx->xi2;
    real_t z2 = ctx->z2;
    real_t z2m2 = 1.0 / (z2 * z2);
    real_t z2m4 = z2m2 * z2m2;
    real_t z2m6 = z2m4 * z2m2;
    real_t z2m8 = z2m6 * z2m2;
    real_t z2m10 = z2m8 * z2m2;
    real_t z2m12 = z2m10 * z2m2;
    real_t beta11 = xi1 - xi2;
    real_t beta13 = xi1 - 3.0 * xi2;
    real_t beta23 = 2.0 * xi1 - 3.0 * xi2;
    real_t beta57 = 5.0 * xi1 - 7.0 * xi2;
    real_t value;

    if (ctx->term == LAMB_BASIC_P_TERM || ctx->term == LAMB_BASIC_S_TERM) {
        /* P、S 两项只在奇次组合的符号上不同 */
        const real_t reflection_sign = ctx->term == LAMB_BASIC_P_TERM ? -1.0 : 1.0;
        real_t A18 = xi1 * xi1 - 8.0 * xi1 * xi2 + 11.0 * xi2 * xi2;
        real_t A110 = xi1 * xi1 - 10.0 * xi1 * xi2 + 17.0 * xi2 * xi2;
        real_t A167 = xi1 * xi1 - 6.0 * xi1 * xi2 + 7.0 * xi2 * xi2;
        real_t A326 = 3.0 * xi1 * xi1 - 26.0 * xi1 * xi2 + 43.0 * xi2 * xi2;
        real_t B133 = xi1 * xi1 * xi1 - 33.0 * xi1 * xi1 * xi2 + 183.0 * xi1 * xi2 * xi2 - 231.0 * xi2 * xi2 * xi2;
        if (number == 8) {
            value = pow(xi2, 5) * H[0] + 5.0 * reflection_sign * pow(xi2, 3) * beta11 * beta23 * z2m2 * H[1] +
                    5.0 * xi2 * beta11 * beta11 * A18 * z2m4 * H[2] - 5.0 * reflection_sign * pow(beta11, 3) * A110 * z2m6 * H[3] -
                    20.0 * pow(beta11, 4) * beta13 * z2m8 * H[4] - 16.0 * reflection_sign * pow(beta11, 5) * z2m10 * H[5];
        } else if (number == 9) {
            value = pow(xi2, 6) * H[0] + 3.0 * reflection_sign * pow(xi2, 4) * beta11 * beta57 * z2m2 * H[1] +
                    15.0 * xi2 * xi2 * beta11 * beta11 * A167 * z2m4 * H[2] + reflection_sign * pow(beta11, 3) * B133 * z2m6 * H[3] +
                    6.0 * pow(beta11, 4) * A326 * z2m8 * H[4] + 48.0 * reflection_sign * pow(beta11, 5) * beta13 * z2m10 * H[5] +
                    32.0 * pow(beta11, 6) * z2m12 * H[6];
        } else {
            GRTRaiseError("Wrong %s-wave high-order V number in lamb3: %d.\n", ctx->term == LAMB_BASIC_P_TERM ? "P" : "S", number);
        }
    } else {
        real_t c1 = ctx->c1;
        real_t c2 = ctx->c2;
        real_t c2m1 = 1.0 / c2;
        real_t c2m2 = c2m1 * c2m1;
        real_t c2m3 = c2m2 * c2m1;
        real_t c2m4 = c2m3 * c2m1;
        real_t c2m5 = c2m4 * c2m1;
        real_t c2m6 = c2m5 * c2m1;
        real_t M1 = PI / (2.0 * sqrt(c1) * sqrt(c1 - c2));
        real_t M2 = PI * (2.0 * c1 - c2) / (4.0 * pow(c1 * (c1 - c2), 1.5));
        real_t M3 = PI * (8.0 * c1 * c1 - 8.0 * c1 * c2 + 3.0 * c2 * c2) / (16.0 * pow(c1 * (c1 - c2), 2.5));
        real_t M4 = PI * (2.0 * c1 - c2) * (8.0 * c1 * c1 - 8.0 * c1 * c2 + 5.0 * c2 * c2) / (32.0 * pow(c1 * (c1 - c2), 3.5));
        real_t M5 = PI * (128.0 * pow(c1, 4) - 256.0 * pow(c1, 3) * c2 + 288.0 * c1 * c1 * c2 * c2 - 160.0 * c1 * pow(c2, 3) + 35.0 * pow(c2, 4)) /
                    (256.0 * pow(c1 * (c1 - c2), 4.5));
        real_t M6 = PI * (2.0 * c1 - c2) *
                    (128.0 * pow(c1, 4) - 256.0 * pow(c1, 3) * c2 + 352.0 * c1 * c1 * c2 * c2 - 224.0 * c1 * pow(c2, 3) + 63.0 * pow(c2, 4)) /
                    (512.0 * pow(c1 * (c1 - c2), 5.5));
        real_t A18 = xi1 * xi1 - 8.0 * xi1 * xi2 + 11.0 * xi2 * xi2;
        real_t A110 = xi1 * xi1 - 10.0 * xi1 * xi2 + 17.0 * xi2 * xi2;
        real_t A167 = xi1 * xi1 - 6.0 * xi1 * xi2 + 7.0 * xi2 * xi2;
        real_t A326 = 3.0 * xi1 * xi1 - 26.0 * xi1 * xi2 + 43.0 * xi2 * xi2;
        real_t A122 = xi1 * xi1 - 22.0 * xi1 * xi2 + 61.0 * xi2 * xi2;
        real_t A336 = 3.0 * xi1 * xi1 - 36.0 * xi1 * xi2 + 73.0 * xi2 * xi2;
        real_t B133 = xi1 * xi1 * xi1 - 33.0 * xi1 * xi1 * xi2 + 183.0 * xi1 * xi2 * xi2 - 231.0 * xi2 * xi2 * xi2;
        real_t beta25 = 2.0 * xi1 - 5.0 * xi2;
        real_t beta313 = 3.0 * xi1 - 13.0 * xi2;
        real_t beta111 = xi1 - 11.0 * xi2;
        real_t beta14 = xi1 - 4.0 * xi2;
        if (number == 8) {
            value = -(pow(xi2, 5) * H[0] - 5.0 * pow(xi2, 3) * beta11 * beta23 * c2m1 * H[1] + 5.0 * xi2 * beta11 * beta11 * A18 * c2m2 * H[2] +
                      5.0 * pow(beta11, 3) * A110 * c2m3 * H[3] - 20.0 * pow(beta11, 4) * beta13 * c2m4 * H[4] + 16.0 * pow(beta11, 5) * c2m5 * H[5] +
                      beta11 * z2 *
                          (5.0 * pow(xi2, 4) * M1 + 10.0 * xi2 * xi2 * beta11 * beta13 * M2 + beta11 * beta11 * A122 * M3 +
                           4.0 * pow(beta11, 3) * beta313 * M4 + 16.0 * pow(beta11, 4) * M5));
        } else if (number == 9) {
            value = -(pow(xi2, 6) * H[0] - 3.0 * pow(xi2, 4) * beta11 * beta57 * c2m1 * H[1] +
                      15.0 * xi2 * xi2 * beta11 * beta11 * A167 * c2m2 * H[2] - pow(beta11, 3) * B133 * c2m3 * H[3] +
                      6.0 * pow(beta11, 4) * A326 * c2m4 * H[4] - 48.0 * pow(beta11, 5) * beta13 * c2m5 * H[5] + 32.0 * pow(beta11, 6) * c2m6 * H[6] +
                      2.0 * beta11 * z2 *
                          (3.0 * pow(xi2, 5) * M1 + 5.0 * pow(xi2, 3) * beta25 * beta11 * M2 + xi2 * beta11 * beta11 * A336 * M3 -
                           3.0 * pow(beta11, 3) * beta13 * beta111 * M4 - 16.0 * pow(beta11, 4) * beta14 * M5 - 16.0 * pow(beta11, 5) * M6));
        } else {
            GRTRaiseError("Wrong S-P-wave high-order V number in lamb3: %d.\n", number);
        }
    }
    return ctx->c_main * value;
}

static cplx_t basic_V(const int number, const LAMB_BASIC_CONTEXT *ctx) {
    /* 判据成立时不要走 V9 的椭圆组合，直接积引入 xi 之前的式子 */
    if (number == 9 && use_reflection_V9_quadrature(ctx)) {
        return numerical_reflection_V9(ctx);
    }
    if (number >= 3 && number <= 7) {
        real_t K = grt_ellipticK(ctx->m_elliptic);
        real_t H[5];
        grt_lamb_calculate_H(ctx, K, H);
        real_t value;
        if (ctx->term == LAMB_BASIC_P_TERM) {
            value = grt_lamb_tail_V_P(number, ctx, H);
        } else if (ctx->term == LAMB_BASIC_S_TERM) {
            value = grt_lamb_tail_V_S(number, ctx, H);
        } else {
            value = grt_lamb_tail_V_SP(number, ctx, H);
        }
        return I * value;
    }
    if (number == 8 || number == 9) {
        return I * reflection_V_high(number, ctx);
    }
    GRTRaiseError("Wrong V basic-integral number in lamb3: %d.\n", number);
}

static void make_reflection_basis(const int ntail, const cplx_t roots[3], const LAMB_BASIC_CONTEXT *ctx, const bool use_U,
                                  LAMB3_REFLECTION_BASIS *basis) {
    memset(basis, 0, sizeof(*basis));
    real_t K = use_U ? 0.0 : grt_ellipticK(ctx->m_elliptic);
    for (int i = 0; i < 3; ++i) {
        if (use_U) {
            grt_lamb_make_U_pair(roots[i], ctx, basis->pair[i]);
        } else if (ctx->term == LAMB_BASIC_P_TERM) {
            grt_lamb_make_V_P_pair(roots[i], ctx, K, basis->pair[i]);
        } else if (ctx->term == LAMB_BASIC_S_TERM) {
            grt_lamb_make_V_S_pair(roots[i], ctx, K, basis->pair[i]);
        } else {
            grt_lamb_make_V_SP_pair(roots[i], ctx, K, basis->pair[i]);
        }
    }
    for (int i = 0; i < ntail; ++i) {
        basis->tail[i] = use_U ? basic_U(i + 3, ctx) : basic_V(i + 3, ctx);
    }
}

static cplx_t evaluate_reflection_pf(const LAMB3_PF_COEFF *coeff, const LAMB3_REFLECTION_BASIS *basis, const real_t t) {
    cplx_t value = 0.0;
    for (int i = 0; i < 3; ++i) {
        value += grt_lamb_eval_time_coeff(coeff->pair[i][0], coeff->time_degree, t) * basis->pair[i][0];
        value += grt_lamb_eval_time_coeff(coeff->pair[i][1], coeff->time_degree, t) * basis->pair[i][1];
    }
    for (int i = 0; i < coeff->ntail; ++i) {
        value += grt_lamb_eval_time_coeff(coeff->tail[i], coeff->time_degree, t) * basis->tail[i];
    }
    return value;
}

static real_t evaluate_reflection_component(const LAMB3_PF_COEFF *coeff_U, const LAMB3_PF_COEFF *coeff_V, const LAMB3_REFLECTION_BASIS *basis_U,
                                            const LAMB3_REFLECTION_BASIS *basis_V, const bool use_U, const real_t t) {
    cplx_t value = 0.0;
    if (use_U) {
        value += evaluate_reflection_pf(coeff_U, basis_U, t);
    }
    value += evaluate_reflection_pf(coeff_V, basis_V, t);
    return cimag(value);
}

/**
 * 评估一组反射项
 * F[i][j] 的索引分别表示接收点分量和源点分量
 * Fk_source[k'][i][j] 的索引依次表示源点坐标方向、接收点分量和源点分量
 * Fk_receiver[k][i][j] 的索引依次表示接收点坐标方向、接收点分量和源点分量
 */
static void evaluate_reflection_matrix(const LAMB3_PF_COEFF numerator[2][3][3], const LAMB3_REFLECTION_BASIS *basis_U,
                                       const LAMB3_REFLECTION_BASIS *basis_V, const bool use_U, const real_t t, const LAMB3_VARS *V,
                                       real_t result[3][3]) {
    result[0][0] = evaluate_reflection_component(&numerator[0][0][0], &numerator[1][0][0], basis_U, basis_V, use_U, t);
    result[0][1] = evaluate_reflection_component(&numerator[0][0][1], &numerator[1][0][1], basis_U, basis_V, use_U, t);
    result[0][2] = evaluate_reflection_component(&numerator[0][0][2], &numerator[1][0][2], basis_U, basis_V, use_U, t);
    result[1][0] = result[0][1];
    result[1][1] = evaluate_reflection_component(&numerator[0][1][1], &numerator[1][1][1], basis_U, basis_V, use_U, t);
    result[2][0] = -result[0][2];
    result[2][2] = evaluate_reflection_component(&numerator[0][2][2], &numerator[1][2][2], basis_U, basis_V, use_U, t);
    if (V->use_angle_ratio) {
        result[1][2] = result[0][2] * V->angle_ratio;
        result[2][1] = -result[1][2];
    } else {
        result[1][2] = evaluate_reflection_component(&numerator[0][1][2], &numerator[1][1][2], basis_U, basis_V, use_U, t);
        result[2][1] = evaluate_reflection_component(&numerator[0][2][1], &numerator[1][2][1], basis_U, basis_V, use_U, t);
    }
}

/** 计算反射项相对源点 x1' 的空间导数矩阵 */
static void evaluate_reflection_source_x1(const LAMB3_PF_COEFF numerator[2][3][3], const LAMB3_REFLECTION_BASIS *basis_U,
                                          const LAMB3_REFLECTION_BASIS *basis_V, const bool use_U, const real_t t, real_t result[3][3]) {
    result[0][0] = evaluate_reflection_component(&numerator[0][0][0], &numerator[1][0][0], basis_U, basis_V, use_U, t);
    result[0][1] = evaluate_reflection_component(&numerator[0][0][1], &numerator[1][0][1], basis_U, basis_V, use_U, t);
    result[0][2] = evaluate_reflection_component(&numerator[0][0][2], &numerator[1][0][2], basis_U, basis_V, use_U, t);
    result[1][0] = result[0][1];
    result[1][1] = evaluate_reflection_component(&numerator[0][1][1], &numerator[1][1][1], basis_U, basis_V, use_U, t);
    result[1][2] = evaluate_reflection_component(&numerator[0][1][2], &numerator[1][1][2], basis_U, basis_V, use_U, t);
    result[2][0] = -result[0][2];
    result[2][1] = -result[1][2];
    result[2][2] = evaluate_reflection_component(&numerator[0][2][2], &numerator[1][2][2], basis_U, basis_V, use_U, t);
}

/** 计算反射项相对源点 x2' 的空间导数矩阵，P 项可复用 x1' 结果 */
static void evaluate_reflection_source_x2(const LAMB3_PF_COEFF numerator[2][3][3], const real_t source_x1[3][3],
                                          const LAMB3_REFLECTION_BASIS *basis_U, const LAMB3_REFLECTION_BASIS *basis_V,
                                          const bool use_U, const bool is_S, const real_t t, const LAMB3_VARS *V,
                                          real_t result[3][3]) {
    if (is_S) {
        result[0][0] = evaluate_reflection_component(&numerator[0][0][0], &numerator[1][0][0], basis_U, basis_V, use_U, t);
        result[0][1] = evaluate_reflection_component(&numerator[0][0][1], &numerator[1][0][1], basis_U, basis_V, use_U, t);
    } else {
        result[0][0] = source_x1[0][1];
        result[0][1] = source_x1[1][1];
    }
    result[0][2] = source_x1[1][2];
    result[1][0] = result[0][1];
    result[1][1] = evaluate_reflection_component(&numerator[0][1][1], &numerator[1][1][1], basis_U, basis_V, use_U, t);
    result[1][2] = evaluate_reflection_component(&numerator[0][1][2], &numerator[1][1][2], basis_U, basis_V, use_U, t);
    result[2][0] = -result[0][2];
    result[2][1] = -result[1][2];
    if (V->use_angle_ratio) {
        result[2][2] = source_x1[2][2] * V->angle_ratio;
    } else {
        result[2][2] = evaluate_reflection_component(&numerator[0][2][2], &numerator[1][2][2], basis_U, basis_V, use_U, t);
    }
}

/** 评估一组反射项并利用第 8.1、8.2 节的矩阵关系减少重复积分 */
static void evaluate_reflection_set(const LAMB3_REFLECTION_PF_SET *coeffs, const LAMB3_REFLECTION_BASIS *basis_U,
                                    const LAMB3_REFLECTION_BASIS *basis_V, const bool use_U, const bool is_S, const real_t t,
                                    const LAMB3_VARS *V, real_t F[3][3], real_t Fk_source[3][3][3],
                                    real_t Fk_receiver[3][3][3]) {
    evaluate_reflection_matrix(coeffs->M, basis_U, basis_V, use_U, t, V, F);
    evaluate_reflection_source_x1(coeffs->dM[0], basis_U, basis_V, use_U, t, Fk_source[0]);
    evaluate_reflection_source_x2(coeffs->dM[1], Fk_source[0], basis_U, basis_V, use_U, is_S, t, V, Fk_source[1]);
    evaluate_reflection_matrix(coeffs->dM[2], basis_U, basis_V, use_U, t, V, Fk_source[2]);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Fk_receiver[0][i][j] = -Fk_source[0][i][j];
            Fk_receiver[1][i][j] = -Fk_source[1][i][j];
            Fk_receiver[2][i][j] = LAMB3_PHI_PI_SIGN[j][i] * Fk_source[2][j][i];
        }
    }
}

static void select_quartic_pairs(const cplx_t roots[4], real_t *q1, real_t *r1, real_t *q2, real_t *r2) {
    bool used[4] = {false, false, false, false};
    real_t q[2] = {0.0, 0.0};
    real_t r[2] = {0.0, 0.0};
    int npair = 0;
    for (int i = 0; i < 4 && npair < 2; ++i) {
        if (used[i]) {
            continue;
        }
        int partner = -1;
        real_t distance = DBL_MAX;
        for (int j = i + 1; j < 4; ++j) {
            if (!used[j]) {
                real_t candidate = fabs(roots[j] - conj(roots[i]));
                if (candidate < distance) {
                    distance = candidate;
                    partner = j;
                }
            }
        }
        if (partner < 0) {
            GRTRaiseError("Cannot pair the roots of the PS quartic in lamb3.\n");
        }
        used[i] = true;
        used[partner] = true;
        q[npair] = -0.5 * (creal(roots[i]) + creal(roots[partner]));
        r[npair] = 0.5 * (fabs(roots[i]) * fabs(roots[i]) + fabs(roots[partner]) * fabs(roots[partner]));
        ++npair;
    }
    if (r[0] >= r[1]) {
        *q1 = q[0];
        *r1 = r[0];
        *q2 = q[1];
        *r2 = r[1];
    } else {
        *q1 = q[1];
        *r1 = r[1];
        *q2 = q[0];
        *r2 = r[0];
    }
}

/** 根据式 (8.3.1)-(8.3.5) 计算 PS 项的四次根式参数 */
static void make_PS_context(const real_t tbar, const real_t p_depth, const real_t s_depth, const LAMB3_VARS *V, LAMB3_PS_CTX *ctx) {
    real_t z = s_depth / V->R;
    real_t zp = p_depth / V->R;
    real_t zplus = z + zp;
    real_t zminus = z - zp;
    real_t ttilde = tbar * V->r_over_R;
    real_t kap1 = V->kap1;
    real_t kap2 = V->kap2;
    real_t coefficient[5] = {
        zplus * zplus + 1.0,           -4.0 * ttilde * kap1 * zplus, 2.0 * (zplus * zminus + 1.0 + 2.0 * (ttilde * ttilde - 1.0) * kap2),
        -4.0 * ttilde * kap1 * zminus, zminus * zminus + 1.0,
    };
    cplx_t roots[4];
    quartic_roots(coefficient, roots);

    real_t q1, r1, q2, r2;
    select_quartic_pairs(roots, &q1, &r1, &q2, &r2);
    real_t denominator = q2 - q1;
    if (fabs(denominator) < 1e-14) {
        GRTRaiseError("The PS quartic has coincident quadratic factors in lamb3.\n");
    }
    real_t delta = (r1 - r2) * (r1 - r2) - 4.0 * (q2 - q1) * (q1 * r2 - q2 * r1);
    real_t root_delta = grt_lamb_positive_sqrt(delta, "PS quartic");
    ctx->xi1 = (r1 - r2 + root_delta) / (2.0 * denominator);
    ctx->xi2 = (r1 - r2 - root_delta) / (2.0 * denominator);
    if (ctx->xi1 <= ctx->xi2) {
        GRT_SWAP(real_t, ctx->xi1, ctx->xi2);
    }

    real_t b1 = (ctx->xi2 + q1) / (ctx->xi2 - ctx->xi1);
    real_t b2 = (ctx->xi2 + q2) / (ctx->xi2 - ctx->xi1);
    real_t c1 = (ctx->xi1 + q1) / (ctx->xi1 - ctx->xi2);
    real_t c2 = (ctx->xi1 + q2) / (ctx->xi1 - ctx->xi2);
    if (b1 <= 0.0 || b2 <= 0.0 || c1 <= 0.0 || c2 <= 0.0) {
        GRTRaiseError("The PS quartic factorization is not positive in lamb3.\n");
    }
    ctx->z1sq = c2 / b2;
    ctx->z2sq = c1 / b1;
    ctx->z1 = sqrt(ctx->z1sq);
    ctx->z2 = sqrt(ctx->z2sq);
    ctx->m = grt_lamb_clamp_elliptic_parameter(ctx->z2sq / ctx->z1sq, 1e-8, "lamb3");
    ctx->cps = -1.0 / ((ctx->xi1 - ctx->xi2) * sqrt(b1 * c2 * (zplus * zplus + 1.0)));
    ctx->asym_delta = (q2 - q1) * (q1 * r2 - q2 * r1) / ((r1 - r2) * (r1 - r2));
}

static cplx_t PS_pair(const int number, const cplx_t c, const LAMB3_PS_CTX *P) {
    LAMB_BASIC_CONTEXT ctx = {0};
    ctx.term = LAMB_BASIC_P_TERM;
    ctx.xi1 = P->xi1;
    ctx.xi2 = P->xi2;
    ctx.z1sq = P->z1sq;
    ctx.z2sq = P->z2sq;
    ctx.z1 = P->z1;
    ctx.z2 = P->z2;
    ctx.m_elliptic = P->m;

    LAMB_BASIC_V_AUX aux;
    grt_lamb_make_V_aux(-c, &ctx, &aux);
    cplx_t h2psq = aux.h2p * aux.h2p;
    cplx_t zeta0 = ctx.xi2 / aux.h2p;
    cplx_t eta0 = 1.0 / aux.h2p;
    cplx_t zeta[2];
    cplx_t eta[2];
    cplx_t zeta_bar[2];
    cplx_t eta_bar[2];
    for (int i = 0; i < 2; ++i) {
        int other = 1 - i;
        cplx_t s = aux.s[i];
        cplx_t denominator = s * (aux.s[other] - s) * h2psq;
        zeta[i] = (ctx.xi1 + s * ctx.xi2) * (s * aux.h2p + aux.h1p) / denominator;
        zeta_bar[i] = (ctx.xi1 - ctx.xi2) * (s * aux.h2m + aux.h1m) / ((aux.s[other] - s) * h2psq);
        eta[i] = (aux.h2p * s * s - aux.D * s + aux.h1p) / denominator;
        eta_bar[i] = 2.0 * (ctx.xi1 - ctx.xi2) * (s * ctx.xi2 + ctx.xi1) / ((aux.s[other] - s) * h2psq);
    }

    cplx_t n1 = ctx.z2sq / aux.s[0];
    cplx_t n2 = ctx.z2sq / aux.s[1];
    const bool use_zeta = number == 1;
    const cplx_t *main_coeff = use_zeta ? zeta : eta;
    const cplx_t *branch_coeff = use_zeta ? zeta_bar : eta_bar;
    const cplx_t main_coeff0 = use_zeta ? zeta0 : eta0;
    cplx_t elliptic = main_coeff0 * grt_ellipticK(ctx.m_elliptic) + main_coeff[0] * grt_ellipticPi(n1, ctx.m_elliptic) +
                      main_coeff[1] * grt_ellipticPi(n2, ctx.m_elliptic);

    cplx_t alpha11 = sqrt(ctx.z1sq - aux.s[0]);
    cplx_t alpha12 = sqrt(ctx.z1sq - aux.s[1]);
    cplx_t alpha21 = sqrt(ctx.z2sq - aux.s[0]);
    cplx_t alpha22 = sqrt(ctx.z2sq - aux.s[1]);
    cplx_t branch = ctx.z1 * branch_coeff[0] / (alpha11 * alpha21) * atan(I * ctx.z2 * alpha11 / (ctx.z1 * alpha21)) +
                    ctx.z1 * branch_coeff[1] / (alpha12 * alpha22) * atan(I * ctx.z2 * alpha12 / (ctx.z1 * alpha22));

    if (grt_lamb_is_real(c)) {
        return I * P->cps * creal(elliptic);
    }
    return I * P->cps * (elliptic + (use_zeta ? -branch : branch));
}

static void calculate_PS_H(const LAMB3_PS_CTX *P, const int alpha, real_t H[7]) {
    real_t a = alpha == 1 ? P->xi1 * P->xi1 / (P->xi2 * P->xi2 * P->z2sq) : 1.0 / P->z2sq;
    real_t m = P->m;
    real_t K = grt_ellipticK(m);
    real_t E = grt_ellipticE(m);
    real_t Hminus1 = ((a * m + 1.0) * K - E) / m;
    H[0] = K;
    H[1] = grt_ellipticPi(-1.0 / a, m) / a;
    real_t gamma1 = 3.0 * m * a * a + 2.0 * a * (m + 1.0) + 1.0;
    real_t gamma2 = m + 1.0 + 3.0 * m * a;
    real_t gamma3 = a * (a + 1.0) * (a * m + 1.0);
    H[2] = (gamma1 * H[1] - m * Hminus1) / (2.0 * gamma3);
    H[3] = (3.0 * gamma1 * H[2] - 2.0 * gamma2 * H[1] + m * H[0]) / (4.0 * gamma3);
    H[4] = (5.0 * gamma1 * H[3] - 4.0 * gamma2 * H[2] + 3.0 * m * H[1]) / (6.0 * gamma3);
    H[5] = (7.0 * gamma1 * H[4] - 6.0 * gamma2 * H[3] + 5.0 * m * H[2]) / (8.0 * gamma3);
    H[6] = (9.0 * gamma1 * H[5] - 8.0 * gamma2 * H[4] + 7.0 * m * H[3]) / (10.0 * gamma3);
}

static real_t PS_Vj_with_H(const int j, const LAMB3_PS_CTX *P, const real_t H1[7], const real_t H2[7]) {
    real_t xi1 = P->xi1;
    real_t xi2 = P->xi2;
    real_t z2 = P->z2;
    real_t z2m2 = 1.0 / (z2 * z2);
    real_t z2m4 = z2m2 * z2m2;
    real_t z2m6 = z2m4 * z2m2;
    real_t z2m8 = z2m6 * z2m2;
    real_t z2m10 = z2m8 * z2m2;
    real_t z2m12 = z2m10 * z2m2;
    real_t xi2z2m2 = 1.0 / (xi2 * xi2 * z2 * z2);
    real_t xi2z2m4 = xi2z2m2 * xi2z2m2;
    real_t xi2z2m6 = xi2z2m4 * xi2z2m2;
    real_t xi2z2m8 = xi2z2m6 * xi2z2m2;
    real_t beta11 = xi1 - xi2;
    real_t beta12 = xi1 - 2.0 * xi2;
    real_t beta13 = xi1 - 3.0 * xi2;
    real_t beta21 = 2.0 * xi1 - xi2;
    real_t beta23 = 2.0 * xi1 - 3.0 * xi2;
    real_t beta31 = 3.0 * xi1 - xi2;
    real_t beta35 = 3.0 * xi1 - 5.0 * xi2;
    real_t beta53 = 5.0 * xi1 - 3.0 * xi2;
    real_t beta57 = 5.0 * xi1 - 7.0 * xi2;
    real_t A18 = xi1 * xi1 - 8.0 * xi1 * xi2 + 11.0 * xi2 * xi2;
    real_t A110 = xi1 * xi1 - 10.0 * xi1 * xi2 + 17.0 * xi2 * xi2;
    real_t A167 = xi1 * xi1 - 6.0 * xi1 * xi2 + 7.0 * xi2 * xi2;
    real_t A114 = xi1 * xi1 - 14.0 * xi1 * xi2 + 25.0 * xi2 * xi2;
    real_t B133 = xi1 * xi1 * xi1 - 33.0 * xi1 * xi1 * xi2 + 183.0 * xi1 * xi2 * xi2 - 231.0 * xi2 * xi2 * xi2;
    real_t A326 = 3.0 * xi1 * xi1 - 26.0 * xi1 * xi2 + 43.0 * xi2 * xi2;
    real_t value;

    if ((j == -4 || j == -3 || j == -2) && fabs(P->asym_delta) < 1e-3) {
        real_t K = H2[0];
        real_t E = grt_ellipticE(P->m);
        real_t I0 = K;
        real_t I1 = (K - E) / P->m;
        real_t I2 = (2.0 * (P->m + 1.0) * I1 - I0) / (3.0 * P->m);
        real_t I3 = (4.0 * (P->m + 1.0) * I2 - 3.0 * I1) / (5.0 * P->m);
        real_t delta = P->asym_delta;
        real_t z22 = P->z2 * P->z2;
        real_t z24 = z22 * z22;
        real_t z26 = z24 * z22;
        if (j == -4) {
            value = pow(P->xi1, -4.0) * (I0 - 2.0 * (1.0 - delta) * (3.0 - 5.0 * delta) * z22 * I1 +
                                         (1.0 - 6.0 * delta) * (1.0 - 10.0 * delta) * z24 * I2 - 10.0 * delta * delta * z26 * I3);
        } else if (j == -3) {
            value = pow(P->xi1, -3.0) * (I0 - 3.0 * (1.0 - delta) * (1.0 - 2.0 * delta) * z22 * I1 - 3.0 * (1.0 - 6.0 * delta) * delta * z24 * I2);
        } else {
            value = pow(P->xi1, -2.0) * (I0 - (1.0 - delta) * (1.0 - 3.0 * delta) * z22 * I1 + 3.0 * delta * delta * z24 * I2);
        }
        return P->cps * value;
    }

    if (j == -4) {
        value = pow(xi2, -4) * (H1[0] - 2.0 * beta11 * beta53 * xi2z2m2 * H1[1] +
                                beta11 * beta11 * (25.0 * xi1 * xi1 - 14.0 * xi1 * xi2 + xi2 * xi2) * xi2z2m4 * H1[2] -
                                8.0 * xi1 * xi1 * pow(beta11, 3) * beta31 * xi2z2m6 * H1[3] + 8.0 * pow(xi1, 4) * pow(beta11, 4) * xi2z2m8 * H1[4]);
    } else if (j == -3) {
        value = pow(xi2, -3) * (H1[0] - 3.0 * beta11 * beta21 * xi2z2m2 * H1[1] + 3.0 * xi1 * beta11 * beta11 * beta31 * xi2z2m4 * H1[2] -
                                4.0 * pow(xi1, 3) * pow(beta11, 3) * xi2z2m6 * H1[3]);
    } else if (j == -2) {
        value = pow(xi2, -2) * (H1[0] - beta11 * beta31 * xi2z2m2 * H1[1] + 2.0 * xi1 * xi1 * beta11 * beta11 * xi2z2m4 * H1[2]);
    } else if (j == -1) {
        value = pow(xi2, -1) * (H1[0] - xi1 * beta11 * xi2z2m2 * H1[1]);
    } else if (j == 0) {
        value = H2[0];
    } else if (j == 1) {
        value = xi2 * H2[0] + beta11 * z2m2 * H2[1];
    } else if (j == 2) {
        value = xi2 * xi2 * H2[0] - beta11 * beta13 * z2m2 * H2[1] + 2.0 * beta11 * beta11 * z2m4 * H2[2];
    } else if (j == 3) {
        value = pow(xi2, 3) * H2[0] - 3.0 * xi2 * beta11 * beta12 * z2m2 * H2[1] - 3.0 * beta11 * beta11 * beta13 * z2m4 * H2[2] +
                4.0 * pow(beta11, 3) * z2m6 * H2[3];
    } else if (j == 4) {
        value = pow(xi2, 4) * H2[0] - 2.0 * xi2 * xi2 * beta11 * beta35 * z2m2 * H2[1] + beta11 * beta11 * A114 * z2m4 * H2[2] -
                8.0 * pow(beta11, 3) * beta13 * z2m6 * H2[3] + 8.0 * pow(beta11, 4) * z2m8 * H2[4];
    } else if (j == 5) {
        value = pow(xi2, 5) * H2[0] - 5.0 * pow(xi2, 3) * beta11 * beta23 * z2m2 * H2[1] + 5.0 * xi2 * beta11 * beta11 * A18 * z2m4 * H2[2] +
                5.0 * pow(beta11, 3) * A110 * z2m6 * H2[3] - 20.0 * pow(beta11, 4) * beta13 * z2m8 * H2[4] + 16.0 * pow(beta11, 5) * z2m10 * H2[5];
    } else if (j == 6) {
        value = pow(xi2, 6) * H2[0] - 3.0 * pow(xi2, 4) * beta11 * beta57 * z2m2 * H2[1] + 15.0 * xi2 * xi2 * beta11 * beta11 * A167 * z2m4 * H2[2] -
                pow(beta11, 3) * B133 * z2m6 * H2[3] + 6.0 * pow(beta11, 4) * A326 * z2m8 * H2[4] - 48.0 * pow(beta11, 5) * beta13 * z2m10 * H2[5] +
                32.0 * pow(beta11, 6) * z2m12 * H2[6];
    } else {
        GRTRaiseError("Wrong PS basic-integral number in lamb3: %d.\n", j);
    }
    return P->cps * value;
}

static void make_PS_basis(const cplx_t roots[3], const LAMB3_PS_CTX *P, LAMB3_PS_BASIS *basis) {
    memset(basis, 0, sizeof(*basis));
    for (int i = 0; i < 3; ++i) {
        basis->pair[i][0] = PS_pair(1, roots[i], P);
        basis->pair[i][1] = PS_pair(2, roots[i], P);
    }
    real_t H1[7];
    real_t H2[7];
    calculate_PS_H(P, 1, H1);
    calculate_PS_H(P, 2, H2);
    /* basis->V[j+4] 对应式 (8.3.12)-(8.3.17) 中的 V_j，j=-4,...,6 */
    for (int j = -4; j <= 6; ++j) {
        basis->V[j + 4] = PS_Vj_with_H(j, P, H1, H2);
    }
}

static cplx_t evaluate_conversion_pf(const LAMB3_PF_COEFF *coeff, const LAMB3_PS_BASIS *basis, const real_t t) {
    cplx_t value = 0.0;
    for (int i = 0; i < 3; ++i) {
        value += grt_lamb_eval_time_coeff(coeff->pair[i][0], coeff->time_degree, t) * basis->pair[i][0];
        value += grt_lamb_eval_time_coeff(coeff->pair[i][1], coeff->time_degree, t) * basis->pair[i][1];
    }
    for (int i = 0; i < coeff->npole; ++i) {
        /* pole[i] 对应 V_(3-i)，即阶数 -(i+1) 的主部积分 */
        value += I * grt_lamb_eval_time_coeff(coeff->pole[i], coeff->time_degree, t) * basis->V[3 - i];
    }
    for (int i = 0; i < coeff->ntail; ++i) {
        value += I * grt_lamb_eval_time_coeff(coeff->tail[i], coeff->time_degree, t) * basis->V[i + 4];
    }
    return value;
}

/** 评估转换项基本矩阵并利用 M10=M01 关系减少一次积分 */
static void evaluate_conversion_matrix(const LAMB3_PF_COEFF numerator[3][3], const LAMB3_PS_BASIS *basis, const real_t t, const real_t scale,
                                       real_t result[3][3]) {
    result[0][0] = scale * cimag(evaluate_conversion_pf(&numerator[0][0], basis, t));
    result[0][1] = scale * cimag(evaluate_conversion_pf(&numerator[0][1], basis, t));
    result[0][2] = scale * cimag(evaluate_conversion_pf(&numerator[0][2], basis, t));
    result[1][0] = result[0][1];
    result[1][1] = scale * cimag(evaluate_conversion_pf(&numerator[1][1], basis, t));
    result[1][2] = scale * cimag(evaluate_conversion_pf(&numerator[1][2], basis, t));
    result[2][0] = scale * cimag(evaluate_conversion_pf(&numerator[2][0], basis, t));
    result[2][1] = scale * cimag(evaluate_conversion_pf(&numerator[2][1], basis, t));
    result[2][2] = scale * cimag(evaluate_conversion_pf(&numerator[2][2], basis, t));
}

/** 按式 (8.3.1) 评估转换项的三个源点导数矩阵 */
static void evaluate_conversion_source_derivatives(const LAMB3_PF_COEFF numerator[3][3][3], const LAMB3_PS_BASIS *basis, const real_t t,
                                                   const real_t scale, real_t result[3][3][3]) {
    result[0][0][0] = scale * cimag(evaluate_conversion_pf(&numerator[0][0][0], basis, t));
    result[0][0][1] = scale * cimag(evaluate_conversion_pf(&numerator[0][0][1], basis, t));
    result[0][0][2] = scale * cimag(evaluate_conversion_pf(&numerator[0][0][2], basis, t));
    result[0][1][0] = result[0][0][1];
    result[0][1][1] = scale * cimag(evaluate_conversion_pf(&numerator[0][1][1], basis, t));
    result[0][1][2] = scale * cimag(evaluate_conversion_pf(&numerator[0][1][2], basis, t));
    result[0][2][0] = scale * cimag(evaluate_conversion_pf(&numerator[0][2][0], basis, t));
    result[0][2][1] = scale * cimag(evaluate_conversion_pf(&numerator[0][2][1], basis, t));
    result[0][2][2] = scale * cimag(evaluate_conversion_pf(&numerator[0][2][2], basis, t));

    result[1][0][0] = result[0][0][1];
    result[1][0][1] = result[0][1][1];
    result[1][0][2] = result[0][1][2];
    result[1][1][0] = result[1][0][1];
    result[1][1][1] = scale * cimag(evaluate_conversion_pf(&numerator[1][1][1], basis, t));
    result[1][1][2] = scale * cimag(evaluate_conversion_pf(&numerator[1][1][2], basis, t));
    result[1][2][0] = result[0][2][1];
    result[1][2][1] = scale * cimag(evaluate_conversion_pf(&numerator[1][2][1], basis, t));
    result[1][2][2] = scale * cimag(evaluate_conversion_pf(&numerator[1][2][2], basis, t));

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[2][i][j] = scale * cimag(evaluate_conversion_pf(&numerator[2][i][j], basis, t));
        }
    }
}

/** 计算一个 P 或 S 反射项，并累加其矩阵和两个坐标导数 */
static void evaluate_reflection_wave(const real_t sbar, const real_t sbar2, const LAMB_BASIC_VARS *basic_vars, const LAMB3_VARS *V,
                                     const LAMB3_REFLECTION_PF_SET *coeffs, const cplx_t roots[3], const LAMB_BASIC_TERM term,
                                     real_t F[3][3], real_t Fk_source[3][3][3], real_t Fk_receiver[3][3][3]) {
    LAMB_BASIC_CONTEXT ctx = {0};
    if (term == LAMB_BASIC_P_TERM) {
        grt_lamb_make_context_P(sbar, sbar2, basic_vars, &ctx);
    } else {
        grt_lamb_make_context_S(sbar, sbar2, basic_vars, &ctx);
    }
    LAMB3_REFLECTION_BASIS basis_U;
    LAMB3_REFLECTION_BASIS basis_V;
    make_reflection_basis(LAMB3_TAIL_SIZE, roots, &ctx, true, &basis_U);
    make_reflection_basis(LAMB3_TAIL_SIZE, roots, &ctx, false, &basis_V);
    real_t value[3][3];
    real_t source[3][3][3];
    real_t receiver[3][3][3];
    evaluate_reflection_set(coeffs, &basis_U, &basis_V, true, term == LAMB_BASIC_S_TERM, sbar, V, value, source, receiver);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F[i][j] += value[i][j];
            for (int k = 0; k < 3; ++k) {
                Fk_source[k][i][j] += source[k][i][j];
                Fk_receiver[k][i][j] += receiver[k][i][j];
            }
        }
    }
}

/**
 * 计算一个 PS 或 SP 转换项，并写入其矩阵和两个坐标导数
 * F[i][j] 的索引分别表示接收点分量和源点分量
 * Fk_source[k'][i][j] 的索引依次表示源点坐标方向、接收点分量和源点分量
 * Fk_receiver[k][i][j] 的索引依次表示接收点坐标方向、接收点分量和源点分量
 */
static void evaluate_conversion_term(const real_t tbar, const LAMB3_VARS *V, const LAMB3_CONVERSION_PF_SET *coeffs,
                                     const bool swap_depth, real_t F[3][3], real_t Fk_source[3][3][3], real_t Fk_receiver[3][3][3]) {
    const real_t p_depth = swap_depth ? V->receiver_depth : V->source_depth;
    const real_t s_depth = swap_depth ? V->source_depth : V->receiver_depth;
    LAMB3_PS_CTX P;
    make_PS_context(tbar, p_depth, s_depth, V, &P);
    LAMB3_PS_BASIS basis;
    make_PS_basis(V->rayleigh_ps, &P, &basis);
    const real_t scale = V->conversion_scale;
    const real_t derivative_scale = V->conversion_dscale;
    evaluate_conversion_matrix(coeffs->M, &basis, tbar, scale, F);
    evaluate_conversion_source_derivatives(coeffs->dM, &basis, tbar, derivative_scale, Fk_source);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Fk_receiver[0][i][j] = -Fk_source[0][i][j];
            Fk_receiver[1][i][j] = -Fk_source[1][i][j];
            Fk_receiver[2][i][j] = derivative_scale * cimag(evaluate_conversion_pf(&coeffs->receiver_vertical[i][j], &basis, tbar));
        }
    }
}

/**
 * 计算反射项和 sPs 项
 * F[i][j] 的索引分别表示接收点分量和源点分量
 * Fk_source[k'][i][j] 的索引依次表示源点坐标方向、接收点分量和源点分量
 * Fk_receiver[k][i][j] 的索引依次表示接收点坐标方向、接收点分量和源点分量
 * Fsps[i][j] 的索引分别表示接收点分量和源点分量
 * Fk_sps[k'][i][j] 的索引依次表示源点坐标方向、接收点分量和源点分量
 * Fk_sps_receiver[k][i][j] 的索引依次表示接收点坐标方向、接收点分量和源点分量
 */
static void reflection_terms(const real_t tbar, const LAMB3_VARS *V, const LAMB3_REFLECTION_PF_SET *P_coeffs, const LAMB3_REFLECTION_PF_SET *S_coeffs,
                             real_t F[3][3], real_t Fk_source[3][3][3], real_t Fk_receiver[3][3][3], real_t Fsps[3][3], real_t Fk_sps[3][3][3],
                             real_t Fk_sps_receiver[3][3][3]) {
    memset(F, 0, sizeof(real_t) * 3 * 3);
    memset(Fk_source, 0, sizeof(real_t) * 3 * 3 * 3);
    memset(Fk_receiver, 0, sizeof(real_t) * 3 * 3 * 3);
    memset(Fsps, 0, sizeof(real_t) * 3 * 3);
    memset(Fk_sps, 0, sizeof(real_t) * 3 * 3 * 3);
    memset(Fk_sps_receiver, 0, sizeof(real_t) * 3 * 3 * 3);

    real_t sbar = V->varsigma * tbar;
    real_t sbar2 = sbar * sbar;
    LAMB_BASIC_VARS W = {V->k, V->k2, V->kp2, V->st_ref, V->ct_ref};

    if (sbar > V->k) {
        evaluate_reflection_wave(sbar, sbar2, &W, V, P_coeffs, V->rayleigh, LAMB_BASIC_P_TERM, F, Fk_source, Fk_receiver);
    }

    if (sbar > 1.0) {
        evaluate_reflection_wave(sbar, sbar2, &W, V, S_coeffs, V->rayleigh_shifted, LAMB_BASIC_S_TERM, F, Fk_source, Fk_receiver);
    }

    if (sbar < 1.0) {
        /* sPs 变量替换只在反射射线角超过临界角时有效 */
        if (V->supercritical && tbar > V->t_sps && tbar < V->reflection_S) {
            LAMB_BASIC_CONTEXT ctx = {0};
            grt_lamb_make_context_SP(sbar, sbar2, &W, &ctx);
            LAMB3_REFLECTION_BASIS basis_V;
            make_reflection_basis(LAMB3_TAIL_SIZE, V->rayleigh_shifted, &ctx, false, &basis_V);
            evaluate_reflection_set(S_coeffs, NULL, &basis_V, false, true, sbar, V, Fsps, Fk_sps, Fk_sps_receiver);
        }
    }
}

/** 计算直达 P 波的源点竖直导数矩阵 */
static void add_direct_P_vertical(const real_t scale, const real_t st, const real_t ct, const real_t sf, const real_t cf, const real_t st2,
                                  const real_t ct2, const real_t a11, const real_t a53, const real_t A, const real_t B, const real_t C,
                                  real_t result[3][3]) {
    result[0][0] += -scale * A * ct;
    result[0][1] += -scale * a53 * st2 * ct * sf * cf;
    result[0][2] += scale * C * st * cf;
    result[1][0] += -scale * a53 * st2 * ct * sf * cf;
    result[1][1] += -scale * B * ct;
    result[1][2] += scale * C * st * sf;
    result[2][0] += scale * C * st * cf;
    result[2][1] += scale * C * st * sf;
    result[2][2] += -scale * (a53 * ct2 - 3.0 * a11) * ct;
}

static void set_direct_P_vertical(const real_t tbar, const LAMB3_VARS *V, real_t result[3][3]) {
    real_t st = V->st;
    real_t ct = V->ct;
    real_t sf = V->sf;
    real_t cf = V->cf;
    real_t st2 = st * st;
    real_t ct2 = ct * ct;
    real_t a11 = tbar * tbar - V->k2;
    real_t a53 = 5.0 * tbar * tbar - 3.0 * V->k2;
    real_t A = a53 * st2 * cf * cf - a11;
    real_t B = a53 * st2 * sf * sf - a11;
    real_t C = a53 * ct2 - a11;

    add_direct_P_vertical(tbar, st, ct, sf, cf, st2, ct2, a11, a53, A, B, C, result);
}

/** 计算直达 S 波的源点竖直导数矩阵 */
static void add_direct_S_vertical(const real_t scale, const real_t st, const real_t ct, const real_t sf, const real_t cf, const real_t st2,
                                  const real_t ct2, const real_t b53, real_t result[3][3]) {
    real_t vertical_horizontal = 5.0 * scale * scale * st2 - 4.0 * scale * scale - 3.0 * st2 + 2.0;
    real_t vertical_11 = 5.0 * scale * scale * st2 * cf * cf - scale * scale - 3.0 * st2 * cf * cf - 1.0;
    real_t vertical_22 = 5.0 * scale * scale * st2 * sf * sf - scale * scale - 3.0 * st2 * sf * sf - 1.0;
    real_t vertical_33 = 5.0 * scale * scale * ct2 - 3.0 * scale * scale - 3.0 * ct2 + 1.0;

    result[0][0] += scale * vertical_11 * ct;
    result[0][1] += scale * b53 * st2 * ct * sf * cf;
    result[0][2] += scale * vertical_horizontal * st * cf;
    result[1][0] += scale * b53 * st2 * ct * sf * cf;
    result[1][1] += scale * vertical_22 * ct;
    result[1][2] += scale * vertical_horizontal * st * sf;
    result[2][0] += scale * vertical_horizontal * st * cf;
    result[2][1] += scale * vertical_horizontal * st * sf;
    result[2][2] += scale * vertical_33 * ct;
}

static void set_direct_S_vertical(const real_t tbar, const LAMB3_VARS *V, real_t result[3][3]) {
    real_t st = V->st;
    real_t ct = V->ct;
    real_t sf = V->sf;
    real_t cf = V->cf;
    real_t st2 = st * st;
    real_t ct2 = ct * ct;
    real_t b53 = 5.0 * tbar * tbar - 3.0;

    add_direct_S_vertical(tbar, st, ct, sf, cf, st2, ct2, b53, result);
}

static void set_direct_P(const real_t tbar, const LAMB3_VARS *V, real_t F[3][3], real_t Fk_source[3][3][3]) {
    real_t st = V->st;
    real_t ct = V->ct;
    real_t sf = V->sf;
    real_t cf = V->cf;
    real_t st2 = st * st;
    real_t ct2 = ct * ct;
    real_t a11 = tbar * tbar - V->k2;
    real_t a31 = 3.0 * tbar * tbar - V->k2;
    real_t a53 = 5.0 * tbar * tbar - 3.0 * V->k2;

    F[0][0] += a31 * st2 * cf * cf - a11;
    F[0][1] += a31 * st2 * sf * cf;
    F[0][2] += -a31 * st * ct * cf;
    F[1][0] += a31 * st2 * sf * cf;
    F[1][1] += a31 * st2 * sf * sf - a11;
    F[1][2] += -a31 * st * ct * sf;
    F[2][0] += -a31 * st * ct * cf;
    F[2][1] += -a31 * st * ct * sf;
    F[2][2] += a31 * ct2 - a11;

    real_t A = a53 * st2 * cf * cf - a11;
    real_t B = a53 * st2 * sf * sf - a11;
    real_t C = a53 * ct2 - a11;
    real_t D = a53 * st2 * cf * cf - 3.0 * a11;
    real_t E = a53 * st2 * sf * sf - 3.0 * a11;
    real_t scale = tbar;
    Fk_source[0][0][0] += scale * D * st * cf;
    Fk_source[0][0][1] += scale * A * st * sf;
    Fk_source[0][0][2] += -scale * A * ct;
    Fk_source[0][1][0] += scale * A * st * sf;
    Fk_source[0][1][1] += scale * B * st * cf;
    Fk_source[0][1][2] += -scale * a53 * st2 * ct * sf * cf;
    Fk_source[0][2][0] += -scale * A * ct;
    Fk_source[0][2][1] += -scale * a53 * st2 * ct * sf * cf;
    Fk_source[0][2][2] += scale * C * st * cf;
    Fk_source[1][0][0] += scale * A * st * sf;
    Fk_source[1][0][1] += scale * B * st * cf;
    Fk_source[1][0][2] += -scale * a53 * st2 * ct * sf * cf;
    Fk_source[1][1][0] += scale * B * st * cf;
    Fk_source[1][1][1] += scale * E * st * sf;
    Fk_source[1][1][2] += -scale * B * ct;
    Fk_source[1][2][0] += -scale * a53 * st2 * ct * sf * cf;
    Fk_source[1][2][1] += -scale * B * ct;
    Fk_source[1][2][2] += scale * C * st * sf;
    add_direct_P_vertical(tbar, st, ct, sf, cf, st2, ct2, a11, a53, A, B, C, Fk_source[2]);
}

static void set_direct_S(const real_t tbar, const LAMB3_VARS *V, real_t F[3][3], real_t Fk_source[3][3][3]) {
    real_t st = V->st;
    real_t ct = V->ct;
    real_t sf = V->sf;
    real_t cf = V->cf;
    real_t st2 = st * st;
    real_t ct2 = ct * ct;
    real_t b11 = tbar * tbar - 1.0;
    real_t b31 = 3.0 * tbar * tbar - 1.0;
    real_t b53 = 5.0 * tbar * tbar - 3.0;
    real_t bbar11 = tbar * tbar + 1.0;

    F[0][0] += -b31 * st2 * cf * cf + bbar11;
    F[0][1] += -b31 * st2 * sf * cf;
    F[0][2] += b31 * st * ct * cf;
    F[1][0] += -b31 * st2 * sf * cf;
    F[1][1] += -b31 * st2 * sf * sf + bbar11;
    F[1][2] += b31 * st * ct * sf;
    F[2][0] += b31 * st * ct * cf;
    F[2][1] += b31 * st * ct * sf;
    F[2][2] += b31 * st2 - 2.0 * b11;

    real_t A = b53 * st2 * cf * cf - b31;
    real_t B = b53 * st2 * cf * cf - b11;
    real_t C = b53 * st2 * sf * sf - bbar11;
    real_t D = b53 * st2 * sf * sf - b11;
    real_t E = b53 * ct2 - bbar11;
    real_t Abar = b53 * st2 * cf * cf - bbar11;
    real_t D31 = b53 * st2 * sf * sf - b31;
    real_t scale = tbar;
    Fk_source[0][0][0] += -scale * A * st * cf;
    Fk_source[0][0][1] += -scale * B * st * sf;
    Fk_source[0][0][2] += scale * B * ct;
    Fk_source[0][1][0] += -scale * B * st * sf;
    Fk_source[0][1][1] += -scale * C * st * cf;
    Fk_source[0][1][2] += scale * b53 * st2 * ct * sf * cf;
    Fk_source[0][2][0] += scale * B * ct;
    Fk_source[0][2][1] += scale * b53 * st2 * ct * sf * cf;
    Fk_source[0][2][2] += -scale * E * st * cf;
    Fk_source[1][0][0] += -scale * Abar * st * sf;
    Fk_source[1][0][1] += -scale * D * st * cf;
    Fk_source[1][0][2] += scale * b53 * st2 * ct * sf * cf;
    Fk_source[1][1][0] += -scale * D * st * cf;
    Fk_source[1][1][1] += -scale * D31 * st * sf;
    Fk_source[1][1][2] += scale * D * ct;
    Fk_source[1][2][0] += scale * b53 * st2 * ct * sf * cf;
    Fk_source[1][2][1] += scale * D * ct;
    Fk_source[1][2][2] += -scale * E * st * sf;
    add_direct_S_vertical(tbar, st, ct, sf, cf, st2, ct2, b53, Fk_source[2]);
}

static void set_direct_receiver(const real_t tbar, const LAMB3_VARS *V, const real_t dF_source[3][3][3], real_t dF_receiver[3][3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            dF_receiver[0][i][j] = -dF_source[0][i][j];
            dF_receiver[1][i][j] = -dF_source[1][i][j];
        }
    }

    LAMB3_VARS W = *V;
    W.ct = -V->ct;
    W.sf = -V->sf;
    W.cf = -V->cf;
    real_t vertical[3][3] = {0};
    if (tbar > V->k) {
        set_direct_P_vertical(tbar, &W, vertical);
    }
    if (tbar > 1.0) {
        set_direct_S_vertical(tbar, &W, vertical);
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            dF_receiver[2][i][j] = vertical[j][i];
        }
    }
}

static void conversion_terms(const real_t tbar, const LAMB3_VARS *V, const LAMB3_CONVERSION_PF_SET *PS_coeffs,
                             const LAMB3_CONVERSION_PF_SET *SP_coeffs, real_t Fps[3][3], real_t Fk_ps[3][3][3], real_t Fk_ps_receiver[3][3][3],
                             real_t Fsp[3][3], real_t Fk_sp[3][3][3], real_t Fk_sp_receiver[3][3][3]) {
    memset(Fps, 0, sizeof(real_t) * 3 * 3);
    memset(Fk_ps, 0, sizeof(real_t) * 3 * 3 * 3);
    memset(Fk_ps_receiver, 0, sizeof(real_t) * 3 * 3 * 3);
    memset(Fsp, 0, sizeof(real_t) * 3 * 3);
    memset(Fk_sp, 0, sizeof(real_t) * 3 * 3 * 3);
    memset(Fk_sp_receiver, 0, sizeof(real_t) * 3 * 3 * 3);
    if (tbar <= V->tps && tbar <= V->tsp) {
        return;
    }

    if (tbar > V->tps) {
        evaluate_conversion_term(tbar, V, PS_coeffs, false, Fps, Fk_ps, Fk_ps_receiver);
    }

    if (tbar > V->tsp) {
        evaluate_conversion_term(tbar, V, SP_coeffs, true, Fsp, Fk_sp, Fk_sp_receiver);
    }
}

static real_t shift_lamb3_boundary(const real_t tbar, const LAMB3_VARS *V, const real_t tbar_eps) {
    // 与 lamb1 一致，精确命中波前时使用下一阶段的右侧值
    if (tbar == V->k || tbar == 1.0 || tbar == V->reflection_P || tbar == V->reflection_S || tbar == V->tps || tbar == V->tsp) {
        return tbar + tbar_eps;
    }
    if (V->supercritical && tbar == V->t_sps) {
        return tbar + tbar_eps;
    }
    return tbar;
}

static void evaluate_time(const real_t tbar, const LAMB3_VARS *V, const LAMB3_REFLECTION_PF_SET *P_reflection,
                          const LAMB3_REFLECTION_PF_SET *S_reflection, const LAMB3_CONVERSION_PF_SET *PS_conversion,
                          const LAMB3_CONVERSION_PF_SET *SP_conversion, LAMB3_F *out) {
    memset(out, 0, sizeof(*out));
    real_t direct[3][3] = {0};
    real_t direct_source[3][3][3] = {0};
    real_t direct_receiver[3][3][3] = {0};
    if (tbar > V->k) {
        set_direct_P(tbar, V, direct, direct_source);
    }
    if (tbar > 1.0) {
        set_direct_S(tbar, V, direct, direct_source);
    }
    set_direct_receiver(tbar, V, direct_source, direct_receiver);

    real_t reflection[3][3];
    real_t reflection_source[3][3][3];
    real_t reflection_receiver[3][3][3];
    real_t sps[3][3];
    real_t sps_source[3][3][3];
    real_t sps_receiver[3][3][3];
    reflection_terms(tbar, V, P_reflection, S_reflection, reflection, reflection_source, reflection_receiver, sps, sps_source, sps_receiver);

    real_t ps[3][3];
    real_t ps_source[3][3][3];
    real_t ps_receiver[3][3][3];
    real_t sp[3][3];
    real_t sp_source[3][3][3];
    real_t sp_receiver[3][3][3];
    conversion_terms(tbar, V, PS_conversion, SP_conversion, ps, ps_source, ps_receiver, sp, sp_source, sp_receiver);

    real_t varsigma = V->varsigma;
    real_t reflection_indicator = V->supercritical ? 1.0 : 0.0;
    /* 将第三类公式中的整体 2.0 归入与 lamb1/lamb2 一致的输出归一化 */
    const real_t output_scale = 1.0 / 2.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            out->F[i][j] += output_scale *
                            (QUARTERPI * direct[i][j] +
                             varsigma * (reflection[i][j] - reflection_indicator * sps[i][j]) - 4.0 * (ps[i][j] + sp[i][j]));
            for (int k = 0; k < 3; ++k) {
                out->Fk_source[k][i][j] += output_scale *
                                           (QUARTERPI * direct_source[k][i][j] +
                                            varsigma * (reflection_source[k][i][j] - reflection_indicator * sps_source[k][i][j]) -
                                            4.0 * (ps_source[k][i][j] + sp_source[k][i][j]));
                out->Fk_receiver[k][i][j] += output_scale *
                                             (QUARTERPI * direct_receiver[k][i][j] +
                                              varsigma * (reflection_receiver[k][i][j] - reflection_indicator * sps_receiver[k][i][j]) -
                                              4.0 * (ps_receiver[k][i][j] + sp_receiver[k][i][j]));
            }
        }
    }
}

void grt_solve_lamb3(
    const real_t nu, const real_t *ts, const int nt, const real_t R, const real_t source_depth, const real_t receiver_depth,
    const real_t azimuth, real_t (*G)[3][3], real_t (*dG_source)[3][3][3], real_t (*dG_receiver)[3][3][3])
{
    if (nu <= 0.0 || nu >= 0.5) {
        GRTRaiseError("poisson ratio (%lf) is out of bound in lamb3.\n", nu);
    }
    if (nu <= LAMB_NU_WARNING_MARGIN || nu >= 0.5 - LAMB_NU_WARNING_MARGIN) {
        GRTRaiseWarning("Poisson ratio (%lf) is close to the boundary of (0, 0.5); calculation is very likely to fail.", nu);
    }
    if (ts == NULL || nt <= 0) {
        GRTRaiseError("The time series for lamb3 should not be empty.\n");
    }
    if (R <= 0.0) {
        GRTRaiseError("The horizontal distance R should be positive in lamb3.\n");
    }
    if (source_depth <= 0.0 || receiver_depth <= 0.0) {
        GRTRaiseError("Source and receiver depths should be strictly positive in lamb3.\n");
    }
    if (azimuth < 0.0 || azimuth > 360.0) {
        GRTRaiseError("azimuth should be in [0, 360] degree for lamb3.\n");
    }
    for (int i = 0; i < nt; ++i) {
        if (ts[i] < 0.0) {
            GRTRaiseError("The time series for lamb3 should be nonnegative.\n");
        }
        if (i > 0 && ts[i] <= ts[i - 1]) {
            GRTRaiseError("The time series for lamb3 should be strictly increasing.\n");
        }
    }

    LAMB3_VARS V;
    make_vars(nu, R, source_depth, receiver_depth, azimuth, &V);
    real_t horizontal_distance_ratio = R / V.rp;
    real_t source_depth_ratio = source_depth / V.r;
    real_t receiver_depth_ratio = receiver_depth / V.r;
    if (horizontal_distance_ratio <= LAMB3_SMALL_R_WARNING_RATIO) {
        GRTRaiseWarning(
            "The horizontal distance ratio R/r'=%e is small in lamb3; calculation is very likely to fail.", horizontal_distance_ratio);
    }
    if (source_depth_ratio <= LAMB_SURFACE_DEPTH_WARNING_RATIO || receiver_depth_ratio <= LAMB_SURFACE_DEPTH_WARNING_RATIO) {
        GRTRaiseWarning(
            "Source or receiver depth ratio to direct distance is close to zero (source/r=%e, receiver/r=%e); calculation is very likely to fail.",
            source_depth_ratio, receiver_depth_ratio);
    }
    const real_t tbar_eps = nt > 1 ? GRT_MIN(1e-8, (ts[1] - ts[0]) * 1e-5) : 1e-8;
    /* 末点若正好落在波前上会被右移，用略大的 tEnd 判断以免漏构造系数 */
    const real_t tEnd = ts[nt - 1] + tbar_eps;
    const bool need_P = tEnd >= V.reflection_P;
    const bool need_S = tEnd >= V.reflection_S || (V.supercritical && ts[0] < V.reflection_S && tEnd >= V.t_sps);
    const bool need_PS = tEnd >= V.tps;
    const bool need_SP = tEnd >= V.tsp;
    /* 大型部分分式系数工作区放在堆上，避免占用线程栈 */
    LAMB3_PF_COEFFICIENTS *coefficients = calloc(1, sizeof(*coefficients));
    make_coefficients(&V, need_P, need_S, need_PS, need_SP, coefficients);

    bool isprint = G == NULL && dG_source == NULL && dG_receiver == NULL;
    real_t(*F)[3][3] = G != NULL ? G : calloc((size_t)nt, sizeof(*F));
    real_t(*Fk_source)[3][3][3] = calloc((size_t)nt, sizeof(*Fk_source));
    real_t(*Fk_receiver)[3][3][3] = calloc((size_t)nt, sizeof(*Fk_receiver));
    real_t(*dG_source_tmp)[3][3][3] = dG_source != NULL ? dG_source : calloc((size_t)nt, sizeof(*dG_source_tmp));
    real_t(*dG_receiver_tmp)[3][3][3] = dG_receiver != NULL ? dG_receiver : calloc((size_t)nt, sizeof(*dG_receiver_tmp));

    for (int i = 0; i < nt; ++i) {
        LAMB3_F value;
        real_t tbar = shift_lamb3_boundary(ts[i], &V, tbar_eps);
        evaluate_time(tbar, &V, &coefficients->P, &coefficients->S, &coefficients->PS, &coefficients->SP, &value);
        memcpy(F[i], value.F, sizeof(value.F));
        memcpy(Fk_source[i], value.Fk_source, sizeof(value.Fk_source));
        memcpy(Fk_receiver[i], value.Fk_receiver, sizeof(value.Fk_receiver));
    }
    grt_lamb_differentiate_Fk(ts, nt, Fk_source, dG_source_tmp);
    grt_lamb_differentiate_Fk(ts, nt, Fk_receiver, dG_receiver_tmp);

    GRT_SAFE_FREE_PTR(coefficients);

    if (isprint) {
        grt_lamb_print_green_series(stdout, ts, nt, F);
    }

    if (G == NULL) {
        GRT_SAFE_FREE_PTR(F);
    }
    GRT_SAFE_FREE_PTR(Fk_source);
    GRT_SAFE_FREE_PTR(Fk_receiver);
    if (dG_source == NULL) {
        GRT_SAFE_FREE_PTR(dG_source_tmp);
    }
    if (dG_receiver == NULL) {
        GRT_SAFE_FREE_PTR(dG_receiver_tmp);
    }
}
