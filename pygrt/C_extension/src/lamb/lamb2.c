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
 *    展开，最后组合基本积分。多项式的 x 系数由 SymPy 预先推导，
 *    部分分式结果在时间循环前预计算，运行时只组合时间幂次和基本积分。
 */

#include "grt/lamb/elliptic.h"
#include "grt/lamb/lamb2.h"
#include "grt/lamb/lamb_basic.h"
#include "grt/lamb/lamb_util.h"

// 预先计算好的多项式系数
#include "lamb2_coeffs.c_"

#define LAMB2_TAIL_SIZE 6
#define LAMB2_U_TAIL_SIZE 4
#define LAMB2_V_TAIL_SIZE 5
#define LAMB2_RATIO_EPS 1e-8
/* 水平距离相对于直线距离过小时，闭合解的角向展开开始失去有效数字 */
#define LAMB2_SMALL_R_WARNING_RATIO 1e-3

typedef LAMB_POLY LAMB2_POLY;

/** 一个分母为三个二次因子的部分分式展开 */
typedef struct {
    cplx_t pair[3][2];            ///< pair[i][0/1]，i=0,1,2 为根索引，0/1 为 x/常数系数
    cplx_t tail[LAMB2_TAIL_SIZE]; ///< tail[m]，m 为多项式商的 x 次数
    int ntail;                    ///< 多项式商的有效项数
} LAMB2_PF;

/** 按时间幂次保存部分分式展开结果 */
typedef struct {
    cplx_t pair[3][2][LAMB2_TIME_SIZE];            ///< pair[i][0/1][r]，i 为根索引，0/1 为 x/常数项，r 为 tbar 次数
    cplx_t tail[LAMB2_TAIL_SIZE][LAMB2_TIME_SIZE]; ///< tail[m][r]，m 为 x 次数，r 为 tbar 次数
    int ntail;                                     ///< 多项式商的有效项数
    int time_degree;                               ///< tbar 的最高次数
} LAMB2_PF_COEFF;

/** 一类积分项的全部部分分式展开结果 */
typedef struct {
    LAMB2_PF_COEFF M[2][3][3];              ///< M[xi][i][j]，xi=0/1 对应 U/V，i,j=0,1,2 分别为接收点和源点分量
    LAMB2_PF_COEFF dM[3][2][3][3];          ///< dM[k][xi][i][j]，k=0,1,2 为源点坐标导数方向，xi=0/1 对应 U/V，i,j 为矩阵分量
    LAMB2_PF_COEFF receiver_vertical[2][3]; ///< receiver_vertical[xi][j]，xi=0/1 对应 U/V，固定为接收点分量 i=2 的 [2][j] 分量
} LAMB2_PF_SET;

/** P 波项和 S 波项所需的全部部分分式展开结果 */
typedef struct {
    LAMB2_PF_SET P; ///< P 波反射项的部分分式结果
    LAMB2_PF_SET S; ///< S 波反射项的部分分式结果
} LAMB2_PF_COEFFICIENTS;

/** 一个时间点上可复用的基本积分值 */
typedef struct {
    cplx_t pair[3][2];            ///< pair[i][0/1] 是第 i 个二次因子的两个基本积分
    cplx_t tail[LAMB2_TAIL_SIZE]; ///< tail[m] 是 x^m 对应的基本积分
} LAMB2_BASIS;

/** 评估一个时间点的部分分式系数所需的公共参数 */
typedef struct {
    const LAMB2_BASIS *basis_U; ///< U 基本积分，SP 项中为 NULL
    const LAMB2_BASIS *basis_V; ///< V 基本积分
    real_t tbar;                ///< 当前无量纲时间
} LAMB2_EVAL_CTX;

/**
 * 按式 (7.2.2) 和式 (7.3.1) 的分母分解式构造多项式
 *
 * R'(x) = 16 k'^2 (x^2+y1)(x^2+y2)(x^2+y3)
 */
static LAMB2_POLY make_rayleigh_denominator(const cplx_t roots[3], const real_t kp2) {
    LAMB2_POLY denominator = {0};
    denominator.c[0] = 16.0 * kp2;
    denominator.degree = 0;
    for (int i = 0; i < 3; ++i) {
        LAMB2_POLY next = {0};
        next.degree = denominator.degree + 2;
        for (int j = 0; j <= denominator.degree; ++j) {
            next.c[j] += denominator.c[j] * roots[i];
            next.c[j + 2] += denominator.c[j];
        }
        denominator = next;
    }
    return denominator;
}

/**
 * 按式 (7.4.1.2) 和式 (7.2.5) 将任意分子多项式展开为部分分式
 *
 * pair[i][0] * x + pair[i][1] 对应 (x^2+y_i)^{-1}，tail 对应多项式商
 */
static void make_partial_fraction(const LAMB2_POLY *numerator, const LAMB2_POLY *denominator, const cplx_t roots[3], const real_t kp2, LAMB2_PF *pf) {
    LAMB2_POLY quotient, remainder;
    grt_lamb_poly_divide(*numerator, *denominator, &quotient, &remainder);

    memset(pf, 0, sizeof(*pf));
    pf->ntail = quotient.degree + 1;
    if (pf->ntail > LAMB2_TAIL_SIZE) {
        GRTRaiseError("The partial-fraction quotient is too large in lamb2.\n");
    }
    for (int i = 0; i < pf->ntail; ++i) {
        pf->tail[i] = quotient.c[i];
    }

    for (int i = 0; i < 3; ++i) {
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

/** 将一个时间幂次的多项式系数展开为部分分式 */
static void make_partial_fraction_coeff(const LAMB2_POLY_COEFF *numerator, const LAMB2_POLY *denominator, const cplx_t roots[3], const real_t kp2,
                                        LAMB2_PF_COEFF *result) {
    memset(result, 0, sizeof(*result));
    result->time_degree = numerator->time_degree;
    for (int r = 0; r <= numerator->time_degree; ++r) {
        LAMB2_POLY polynomial = {0};
        polynomial.degree = numerator->degree;
        for (int m = 0; m <= numerator->degree; ++m) {
            polynomial.c[m] = numerator->c[r][m];
        }

        LAMB2_PF pf;
        make_partial_fraction(&polynomial, denominator, roots, kp2, &pf);
        if (pf.ntail > result->ntail) {
            result->ntail = pf.ntail;
        }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 2; ++j) {
                result->pair[i][j][r] = pf.pair[i][j];
            }
        }
        for (int i = 0; i < pf.ntail; ++i) {
            result->tail[i][r] = pf.tail[i];
        }
    }
}

/** 将 -x 乘到多项式后展开为部分分式 */
static void make_partial_fraction_coeff_neg_x(const LAMB2_POLY_COEFF *numerator, const LAMB2_POLY *denominator, const cplx_t roots[3],
                                              const real_t kp2, LAMB2_PF_COEFF *result) {
    LAMB2_POLY_COEFF shifted = {0};
    shifted.degree = numerator->degree + 1;
    shifted.time_degree = numerator->time_degree;
    for (int r = 0; r <= numerator->time_degree; ++r) {
        for (int m = 0; m <= numerator->degree; ++m) {
            shifted.c[r][m + 1] = -numerator->c[r][m];
        }
    }
    make_partial_fraction_coeff(&shifted, denominator, roots, kp2, result);
}

/** 按已知的矩阵比例关系复用一组部分分式系数 */
static void scale_partial_fraction_coeff(const LAMB2_PF_COEFF *source, const real_t scale, LAMB2_PF_COEFF *target) {
    *target = *source;
    for (int r = 0; r <= target->time_degree; ++r) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 2; ++j) {
                target->pair[i][j][r] *= scale;
            }
        }
        for (int i = 0; i < target->ntail; ++i) {
            target->tail[i][r] *= scale;
        }
    }
}

/** 为一类积分项预计算全部部分分式展开结果 */
static void make_partial_fraction_set(const LAMB2_COEFF_SET *coeffs, const cplx_t roots[3], const real_t kp2, const bool is_S,
                                      const bool use_angle_ratio, const real_t angle_ratio, LAMB2_PF_SET *result) {
    LAMB2_POLY denominator = make_rayleigh_denominator(roots, kp2);
    memset(result, 0, sizeof(*result));

    for (int xi = 0; xi < 2; ++xi) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                const bool reuse_M = (i == 1 && j == 0) ||
                                     (use_angle_ratio && ((i == 1 && j == 2) || (i == 2 && j == 1)));
                if (!reuse_M) {
                    make_partial_fraction_coeff(&coeffs->M[xi][i][j], &denominator, roots, kp2, &result->M[xi][i][j]);
                    make_partial_fraction_coeff_neg_x(&coeffs->M[xi][i][j], &denominator, roots, kp2, &result->dM[2][xi][i][j]);
                }
                for (int k = 0; k < 2; ++k) {
                    const bool reuse_derivative = (i == 1 && j == 0) ||
                                                  (k == 1 && ((i == 0 && j == 2) || (i == 2 && j == 0))) ||
                                                  (k == 1 && !is_S && i == 0 && j <= 1) ||
                                                  (k == 1 && use_angle_ratio && i == 2 && j == 2);
                    if (!reuse_derivative) {
                        make_partial_fraction_coeff(&coeffs->dM[xi][k][i][j], &denominator, roots, kp2, &result->dM[k][xi][i][j]);
                    }
                }
            }
        }

        for (int j = 0; j < 3; ++j) {
            make_partial_fraction_coeff(&coeffs->receiver_vertical[xi][j], &denominator, roots, kp2, &result->receiver_vertical[xi][j]);
        }

        result->M[xi][1][0] = result->M[xi][0][1];
        result->dM[2][xi][1][0] = result->dM[2][xi][0][1];
        for (int k = 0; k < 2; ++k) {
            result->dM[k][xi][1][0] = result->dM[k][xi][0][1];
        }
        result->dM[1][xi][0][2] = result->dM[0][xi][1][2];
        result->dM[1][xi][2][0] = result->dM[0][xi][2][1];
        if (!is_S) {
            result->dM[1][xi][0][0] = result->dM[0][xi][0][1];
            result->dM[1][xi][0][1] = result->dM[0][xi][1][1];
            result->dM[1][xi][1][0] = result->dM[1][xi][0][1];
        }
        if (use_angle_ratio) {
            scale_partial_fraction_coeff(&result->M[xi][0][2], angle_ratio, &result->M[xi][1][2]);
            scale_partial_fraction_coeff(&result->M[xi][2][0], angle_ratio, &result->M[xi][2][1]);
            scale_partial_fraction_coeff(&result->dM[0][xi][2][2], angle_ratio, &result->dM[1][xi][2][2]);
            scale_partial_fraction_coeff(&result->dM[2][xi][0][2], angle_ratio, &result->dM[2][xi][1][2]);
            scale_partial_fraction_coeff(&result->dM[2][xi][2][0], angle_ratio, &result->dM[2][xi][2][1]);
        }
    }

    for (int xi = 0; xi < 2; ++xi) {
        int max_ntail = (xi == 0) ? LAMB2_U_TAIL_SIZE : LAMB2_V_TAIL_SIZE;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (result->M[xi][i][j].ntail > max_ntail) {
                    GRTRaiseError("The M partial-fraction tail is too large in lamb2.\n");
                }
                for (int k = 0; k < 3; ++k) {
                    if (result->dM[k][xi][i][j].ntail > max_ntail) {
                        GRTRaiseError("The derivative partial-fraction tail is too large in lamb2.\n");
                    }
                }
            }
        }
        for (int j = 0; j < 3; ++j) {
            if (result->receiver_vertical[xi][j].ntail > max_ntail) {
                GRTRaiseError("The receiver partial-fraction tail is too large in lamb2.\n");
            }
        }
    }
}

/** 缓存当前时间点和一组路径对应的基本积分 */
static void make_basis(const int ntail, const cplx_t roots[3], const LAMB_BASIC_CONTEXT *ctx, const bool use_V, LAMB2_BASIS *basis) {
    real_t K = use_V ? grt_ellipticK(ctx->m_elliptic) : 0.0;
    for (int i = 0; i < 3; ++i) {
        if (use_V) {
            if (ctx->term == LAMB_BASIC_P_TERM) {
                grt_lamb_make_V_P_pair(roots[i], ctx, K, basis->pair[i]);
            } else if (ctx->term == LAMB_BASIC_S_TERM) {
                grt_lamb_make_V_S_pair(roots[i], ctx, K, basis->pair[i]);
            } else {
                grt_lamb_make_V_SP_pair(roots[i], ctx, K, basis->pair[i]);
            }
        } else {
            grt_lamb_make_U_pair(roots[i], ctx, basis->pair[i]);
        }
    }
    if (use_V) {
        real_t H[5];
        grt_lamb_calculate_H(ctx, K, H);
        for (int i = 0; i < ntail; ++i) {
            if (ctx->term == LAMB_BASIC_P_TERM) {
                basis->tail[i] = I * grt_lamb_tail_V_P(i + 3, ctx, H);
            } else if (ctx->term == LAMB_BASIC_S_TERM) {
                basis->tail[i] = I * grt_lamb_tail_V_S(i + 3, ctx, H);
            } else {
                basis->tail[i] = I * grt_lamb_tail_V_SP(i + 3, ctx, H);
            }
        }
    } else {
        for (int i = 0; i < ntail; ++i) {
            basis->tail[i] = grt_lamb_basic_U(i + 3, 0.0, ctx);
        }
    }
}

static inline real_t evaluate_partial_fraction_coeff(const LAMB2_PF_COEFF *coeff, const LAMB2_BASIS *basis, const LAMB2_EVAL_CTX *ctx) {
    cplx_t value = 0.0;
    for (int i = 0; i < 3; ++i) {
        value += grt_lamb_eval_time_coeff(coeff->pair[i][0], coeff->time_degree, ctx->tbar) * basis->pair[i][0];
        value += grt_lamb_eval_time_coeff(coeff->pair[i][1], coeff->time_degree, ctx->tbar) * basis->pair[i][1];
    }
    for (int i = 0; i < coeff->ntail; ++i) {
        value += grt_lamb_eval_time_coeff(coeff->tail[i], coeff->time_degree, ctx->tbar) * basis->tail[i];
    }
    return cimag(value);
}

static inline real_t evaluate_component(const LAMB2_PF_COEFF *numerator_U, const LAMB2_PF_COEFF *numerator_V, const LAMB2_EVAL_CTX *ctx) {
    if (ctx->basis_U != NULL) {
        real_t result = evaluate_partial_fraction_coeff(numerator_U, ctx->basis_U, ctx);
        result += evaluate_partial_fraction_coeff(numerator_V, ctx->basis_V, ctx);
        return result;
    }
    return evaluate_partial_fraction_coeff(numerator_V, ctx->basis_V, ctx);
}

/**
 * numerator[xi][i][j] 的索引依次为 U/V、接收点分量和源点分量
 * result[i][j] 保存对应的位移矩阵分量
 */
static void evaluate_base_matrix(const LAMB2_PF_COEFF numerator[2][3][3], const LAMB2_VARS *V, const LAMB2_EVAL_CTX *ctx, real_t result[3][3]) {
    result[0][0] = evaluate_component(&numerator[0][0][0], &numerator[1][0][0], ctx);
    result[0][1] = evaluate_component(&numerator[0][0][1], &numerator[1][0][1], ctx);
    result[0][2] = evaluate_component(&numerator[0][0][2], &numerator[1][0][2], ctx);
    result[1][1] = evaluate_component(&numerator[0][1][1], &numerator[1][1][1], ctx);
    result[2][0] = evaluate_component(&numerator[0][2][0], &numerator[1][2][0], ctx);
    result[2][2] = evaluate_component(&numerator[0][2][2], &numerator[1][2][2], ctx);
    result[1][0] = result[0][1];

    if (V->use_angle_ratio) {
        result[1][2] = result[0][2] * V->angle_ratio;
        result[2][1] = result[2][0] * V->angle_ratio;
    } else {
        result[1][2] = evaluate_component(&numerator[0][1][2], &numerator[1][1][2], ctx);
        result[2][1] = evaluate_component(&numerator[0][2][1], &numerator[1][2][1], ctx);
    }
}

static inline real_t evaluate_source_component(const LAMB2_PF_SET *coeffs, const int direction, const int i, const int j, const LAMB2_EVAL_CTX *ctx) {
    return evaluate_component(&coeffs->dM[direction][0][i][j], &coeffs->dM[direction][1][i][j], ctx);
}

/**
 * 计算源点导数并利用表内的矩阵关系减少重复分量
 * result[k'][i][j] 的索引依次为源点坐标方向、接收点分量和源点分量
 */
static void evaluate_source_derivatives(const LAMB2_PF_SET *coeffs, const LAMB2_VARS *V, const bool is_S, const LAMB2_EVAL_CTX *ctx, real_t result[3][3][3]) {
    real_t d11 = evaluate_source_component(coeffs, 0, 0, 0, ctx);
    real_t d12 = evaluate_source_component(coeffs, 0, 0, 1, ctx);
    real_t d13 = evaluate_source_component(coeffs, 0, 0, 2, ctx);
    real_t d22 = evaluate_source_component(coeffs, 0, 1, 1, ctx);
    real_t d23 = evaluate_source_component(coeffs, 0, 1, 2, ctx);
    real_t d31 = evaluate_source_component(coeffs, 0, 2, 0, ctx);
    real_t d32 = evaluate_source_component(coeffs, 0, 2, 1, ctx);
    real_t d33 = evaluate_source_component(coeffs, 0, 2, 2, ctx);

    result[0][0][0] = d11;
    result[0][0][1] = d12;
    result[0][0][2] = d13;
    result[0][1][0] = d12;
    result[0][1][1] = d22;
    result[0][1][2] = d23;
    result[0][2][0] = d31;
    result[0][2][1] = d32;
    result[0][2][2] = d33;

    /* S 项的 [1][0][0]、[1][0][1] 没有 P 项中的对称关系 */
    const real_t d_second_00 = is_S ? evaluate_source_component(coeffs, 1, 0, 0, ctx) : d12;
    const real_t d_second_01 = is_S ? evaluate_source_component(coeffs, 1, 0, 1, ctx) : d22;
    result[1][0][0] = d_second_00;
    result[1][0][1] = d_second_01;
    result[1][0][2] = d23;
    result[1][1][0] = d_second_01;
    result[1][1][1] = evaluate_source_component(coeffs, 1, 1, 1, ctx);
    result[1][1][2] = evaluate_source_component(coeffs, 1, 1, 2, ctx);
    result[1][2][0] = d32;
    result[1][2][1] = evaluate_source_component(coeffs, 1, 2, 1, ctx);
    const bool use_ratio = is_S && V->use_angle_ratio;
    result[1][2][2] = use_ratio ? result[0][2][2] * V->angle_ratio : evaluate_source_component(coeffs, 1, 2, 2, ctx);

    evaluate_base_matrix(coeffs->dM[2], V, ctx, result[2]);
}

/**
 * 利用平移关系和边界关系恢复接收点导数
 * source[k'][i][j] 的索引依次为源点坐标方向、接收点分量和源点分量
 * result[k][i][j] 的索引依次为接收点坐标方向、接收点分量和源点分量
 */
static void evaluate_receiver_derivatives(const LAMB2_PF_SET *coeffs, const LAMB2_EVAL_CTX *ctx, const real_t source[3][3][3], real_t result[3][3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[0][i][j] = -source[0][i][j];
            result[1][i][j] = -source[1][i][j];
        }
        result[2][0][i] = source[0][2][i];
        result[2][1][i] = source[1][2][i];
    }

    for (int j = 0; j < 3; ++j) {
        result[2][2][j] = evaluate_component(&coeffs->receiver_vertical[0][j], &coeffs->receiver_vertical[1][j], ctx);
    }
}

/**
 * 计算一个 P、S 或 S-P 项，并按 sign 累加到当前时间点的输出
 * F[i][j] 的索引分别表示接收点分量和源点分量
 * Fk_source[k'][i][j] 的索引依次表示源点坐标方向、接收点分量和源点分量
 * Fk_receiver[k][i][j] 的索引依次表示接收点坐标方向、接收点分量和源点分量
 */
static void evaluate_lamb2_term(const real_t tbar, const real_t tbar2, const LAMB_BASIC_VARS *basic_vars,
                                const LAMB2_PF_SET *coeffs, const cplx_t roots[3], const LAMB_BASIC_TERM term, const bool use_U,
                                const real_t sign, const LAMB2_VARS *V, real_t F[3][3],
                                real_t Fk_source[3][3][3], real_t Fk_receiver[3][3][3]) {
    LAMB_BASIC_CONTEXT ctx = {0};
    if (term == LAMB_BASIC_P_TERM) {
        grt_lamb_make_context_P(tbar, tbar2, basic_vars, &ctx);
    } else if (term == LAMB_BASIC_S_TERM) {
        grt_lamb_make_context_S(tbar, tbar2, basic_vars, &ctx);
    } else {
        grt_lamb_make_context_SP(tbar, tbar2, basic_vars, &ctx);
    }

    LAMB2_BASIS basis_U = {0};
    LAMB2_BASIS basis_V;
    const LAMB2_BASIS *basis_U_ptr = NULL;
    if (use_U) {
        make_basis(LAMB2_U_TAIL_SIZE, roots, &ctx, false, &basis_U);
        basis_U_ptr = &basis_U;
    }
    make_basis(LAMB2_V_TAIL_SIZE, roots, &ctx, true, &basis_V);
    LAMB2_EVAL_CTX eval = {basis_U_ptr, &basis_V, tbar};

    real_t value[3][3];
    real_t dvalue_source[3][3][3];
    real_t dvalue_receiver[3][3][3];
    evaluate_base_matrix(coeffs->M, V, &eval, value);
    evaluate_source_derivatives(coeffs, V, term == LAMB_BASIC_S_TERM, &eval, dvalue_source);
    evaluate_receiver_derivatives(coeffs, &eval, dvalue_source, dvalue_receiver);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F[i][j] += sign * value[i][j];
            for (int k = 0; k < 3; ++k) {
                Fk_source[k][i][j] += sign * dvalue_source[k][i][j];
                Fk_receiver[k][i][j] += sign * dvalue_receiver[k][i][j];
            }
        }
    }
}

/**
 * 计算一个时间点的 F、F_(,k') 和 F_(,k)
 *
 * 这里直接对应式 (7.4.1.4)。P 项、S 项和 S-P 项共用同一套评估流程，
 * 仅由到时条件、根组、基本积分类型和累加符号区分。接收点导数的积分
 * 分子按接收点边界关系构造
 */
static void evaluate_lamb2_time(const real_t tbar, const LAMB2_VARS *V, const LAMB2_PF_COEFFICIENTS *coefficients, real_t F[3][3],
                                real_t Fk_source[3][3][3], real_t Fk_receiver[3][3][3]) {
    memset(F, 0, sizeof(real_t) * 3 * 3);
    memset(Fk_source, 0, sizeof(real_t) * 3 * 3 * 3);
    memset(Fk_receiver, 0, sizeof(real_t) * 3 * 3 * 3);

    const real_t tbar2 = tbar * tbar;
    const LAMB_BASIC_VARS basic_vars = {V->k, V->k2, V->kp2, V->st, V->ct};

    if (tbar > V->k) {
        evaluate_lamb2_term(tbar, tbar2, &basic_vars, &coefficients->P, V->y, LAMB_BASIC_P_TERM, true, 1.0, V, F,
                            Fk_source, Fk_receiver);
    }
    if (tbar > 1.0) {
        evaluate_lamb2_term(tbar, tbar2, &basic_vars, &coefficients->S, V->yp, LAMB_BASIC_S_TERM, true, 1.0, V, F,
                            Fk_source, Fk_receiver);
    }
    if (V->supercritical && tbar > V->t_sp && tbar < 1.0) {
        evaluate_lamb2_term(tbar, tbar2, &basic_vars, &coefficients->S, V->yp, LAMB_BASIC_SP_TERM, false, -1.0, V, F,
                            Fk_source, Fk_receiver);
    }
}

/**
 * 将地下源、地表接收的结果转为地表源、地下接收
 *
 * 对偶问题已按方位角 phi+pi 计算，这里只需转置 G_ij，
 * 并将源点导数与接收点导数对调后再转置 i,j
 */
static void apply_lamb2_surface_source_reciprocity(
    const int nt, real_t (*G)[3][3], real_t (*dG_source)[3][3][3], real_t (*dG_receiver)[3][3][3])
{
    for (int n = 0; n < nt; ++n) {
        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 3; ++j) {
                const real_t tmp = G[n][i][j];
                G[n][i][j] = G[n][j][i];
                G[n][j][i] = tmp;
            }
        }
        for (int k = 0; k < 3; ++k) {
            real_t source[3][3], receiver[3][3];
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    source[i][j] = dG_source[n][k][i][j];
                    receiver[i][j] = dG_receiver[n][k][i][j];
                }
            }
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    dG_source[n][k][i][j] = receiver[j][i];
                    dG_receiver[n][k][i][j] = source[j][i];
                }
            }
        }
    }
}


static real_t shift_lamb2_boundary(const real_t tbar, const LAMB2_VARS *V, const real_t tbar_eps) {
    // 与 lamb1 一致，精确命中波前时使用下一阶段的右侧值
    if (tbar == V->k || tbar == 1.0 || (V->supercritical && tbar == V->t_sp)) {
        return tbar + tbar_eps;
    }
    return tbar;
}

void grt_solve_lamb2(
    const real_t nu, const real_t *ts, const int nt, const real_t R, const real_t source_depth, const real_t receiver_depth,
    const real_t azimuth, real_t (*G)[3][3], real_t (*dG_source)[3][3][3], real_t (*dG_receiver)[3][3][3])
{
    if (nu <= 0.0 || nu >= 0.5) {
        GRTRaiseError("poisson ratio (%lf) is out of bound.", nu);
    }
    if (nu <= LAMB_NU_WARNING_MARGIN || nu >= 0.5 - LAMB_NU_WARNING_MARGIN) {
        GRTRaiseWarning("Poisson ratio (%lf) is close to the boundary of (0, 0.5); calculation is very likely to fail.", nu);
    }
    if (ts == NULL || nt <= 0) {
        GRTRaiseError("The time series for lamb2 should not be empty.\n");
    }
    if (R <= 0.0) {
        GRTRaiseError("The horizontal distance R should be positive in lamb2.\n");
    }
    if (source_depth < 0.0 || receiver_depth < 0.0) {
        GRTRaiseError("The source and receiver depths should be nonnegative in lamb2.\n");
    }
    const bool buried_source = source_depth > 0.0 && receiver_depth == 0.0;
    const bool surface_source = source_depth == 0.0 && receiver_depth > 0.0;
    if (!buried_source && !surface_source) {
        GRTRaiseError("lamb2 requires exactly one of source and receiver depths to be strictly positive, and the other to be zero.\n");
    }
    const real_t buried_depth = buried_source ? source_depth : receiver_depth;
    real_t solve_azimuth = azimuth;
    if (surface_source) {
        solve_azimuth = azimuth + 180.0;
        if (solve_azimuth >= 360.0) {
            solve_azimuth -= 360.0;
        }
    }
    real_t direct_distance = hypot(R, buried_depth);
    real_t horizontal_distance_ratio = R / direct_distance;
    real_t buried_depth_ratio = buried_depth / direct_distance;
    if (horizontal_distance_ratio <= LAMB2_SMALL_R_WARNING_RATIO) {
        GRTRaiseWarning(
            "The horizontal distance ratio R/r=%e is small in lamb2; calculation is very likely to fail.", horizontal_distance_ratio);
    }
    if (buried_depth_ratio <= LAMB_SURFACE_DEPTH_WARNING_RATIO) {
        GRTRaiseWarning(
            "The lamb2 underground depth ratio h/r=%e is close to the free surface; calculation is very likely to fail.", buried_depth_ratio);
    }
    if (azimuth < 0.0 || azimuth > 360.0) {
        GRTRaiseError("azimuth should be in [0, 360] degree for lamb2.\n");
    }
    for (int i = 0; i < nt; ++i) {
        if (ts[i] < 0.0) {
            GRTRaiseError("The time series for lamb2 should be nonnegative.\n");
        }
        if (i > 0 && ts[i] <= ts[i - 1]) {
            GRTRaiseError("The time series for lamb2 should be strictly increasing.\n");
        }
    }
    LAMB2_VARS V = {0};
    V.nu = nu;
    V.k2 = 0.5 * (1.0 - 2.0 * nu) / (1.0 - nu);
    V.k = sqrt(V.k2);
    V.kp2 = 1.0 - V.k2;
    V.kp = sqrt(V.kp2);
    V.theta = atan2(R, buried_depth);
    V.phi = solve_azimuth * DEG1;
    V.st = sin(V.theta);
    V.ct = cos(V.theta);
    V.sf = sin(V.phi);
    V.cf = cos(V.phi);
    V.theta_c = asin(V.k);
    V.t_sp = cos(V.theta - V.theta_c);
    V.use_angle_ratio = fabs(V.cf) > LAMB2_RATIO_EPS;
    V.angle_ratio = V.use_angle_ratio ? V.sf / V.cf : 0.0;
    V.supercritical = V.theta > V.theta_c;
    const real_t tbar_eps = nt > 1 ? GRT_MIN(1e-8, (ts[1] - ts[0]) * 1e-5) : 1e-8;
    grt_rayleigh1_roots(V.nu, V.y);
    for (int i = 0; i < 3; ++i) {
        V.yp[i] = V.y[i] - V.kp2;
    }

    /* 末点若正好落在波前上会被右移，用略大的 tEnd 判断以免漏构造系数 */
    const real_t tEnd = ts[nt - 1] + tbar_eps;
    bool need_P = tEnd >= V.k;
    bool need_S = tEnd >= 1.0 || (V.supercritical && ts[0] < 1.0 && tEnd >= V.t_sp);
    /* 大型多项式系数工作区放在堆上，避免占用线程栈 */
    LAMB2_COEFF_SET *coefficients = calloc(1, sizeof(*coefficients));
    LAMB2_PF_COEFFICIENTS *pf_coefficients = calloc(1, sizeof(*pf_coefficients));
    if (need_P) {
        make_lamb2_P_coefficients(&V, coefficients);
        make_partial_fraction_set(coefficients, V.y, V.kp2, false, V.use_angle_ratio, V.angle_ratio, &pf_coefficients->P);
    }
    if (need_S) {
        make_lamb2_S_coefficients(&V, coefficients);
        make_partial_fraction_set(coefficients, V.yp, V.kp2, true, V.use_angle_ratio, V.angle_ratio, &pf_coefficients->S);
    }
    GRT_SAFE_FREE_PTR(coefficients);

    bool isprint = (G == NULL && dG_source == NULL && dG_receiver == NULL);
    real_t(*F)[3][3] = G != NULL ? G : calloc((size_t)nt, sizeof(*F));
    real_t(*Fk_source)[3][3][3] = calloc((size_t)nt, sizeof(*Fk_source));
    real_t(*Fk_receiver)[3][3][3] = calloc((size_t)nt, sizeof(*Fk_receiver));
    real_t(*dG_source_tmp)[3][3][3] = dG_source != NULL ? dG_source : calloc((size_t)nt, sizeof(*dG_source_tmp));
    real_t(*dG_receiver_tmp)[3][3][3] = dG_receiver != NULL ? dG_receiver : calloc((size_t)nt, sizeof(*dG_receiver_tmp));

    for (int i = 0; i < nt; ++i) {
        real_t tbar = shift_lamb2_boundary(ts[i], &V, tbar_eps);
        evaluate_lamb2_time(tbar, &V, pf_coefficients, F[i], Fk_source[i], Fk_receiver[i]);
    }
    grt_lamb_differentiate_Fk(ts, nt, Fk_source, dG_source_tmp);
    grt_lamb_differentiate_Fk(ts, nt, Fk_receiver, dG_receiver_tmp);
    if (surface_source) {
        apply_lamb2_surface_source_reciprocity(nt, F, dG_source_tmp, dG_receiver_tmp);
    }

    GRT_SAFE_FREE_PTR(pf_coefficients);

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
