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
 *    展开，最后组合基本积分；多项式系数由未展开谱域矩阵统一构造，
 *    部分分式结果在时间循环前预计算，运行时只组合时间幂次和基本积分
 */

#include <string.h>

#include "grt/lamb/elliptic.h"
#include "grt/lamb/lamb2.h"
#include "grt/lamb/lamb_basic.h"
#include "grt/lamb/lamb_poly.h"
#include "grt/lamb/lamb_util.h"

// 多项式系数构造
#include "lamb2_coeffs.c_"

#define LAMB2_TAIL_SIZE 6
#define LAMB2_RATIO_EPS 1e-8
/* 水平距离相对于直线距离过小时，闭合解的角向展开开始失去有效数字 */
#define LAMB2_SMALL_R_WARNING_RATIO 1e-3

enum {
    LAMB2_PHASE_P  = 1u << 0,  ///< 直达 P 波
    LAMB2_PHASE_S  = 1u << 1,  ///< 直达 S 波
    LAMB2_PHASE_SP = 1u << 2,  ///< 地下源到地表接收点的 SP 波
    LAMB2_PHASE_PS = 1u << 3,  ///< 地表源到地下接收点的 PS 波
};

static const GRT_LAMB_PHASE_OPTION LAMB2_PHASE_OPTIONS[] = {
    {"P", LAMB2_PHASE_P},
    {"S", LAMB2_PHASE_S},
    {"SP", LAMB2_PHASE_SP},
    {"PS", LAMB2_PHASE_PS},
};

/** 一个分母为三个二次因子的部分分式展开 */
typedef struct {
    cplx_t pair[3][2];            ///< pair[i][0/1]，i=0,1,2 为根索引，0/1 为 x/常数系数
    cplx_t tail[LAMB2_TAIL_SIZE]; ///< tail[m]，m 为多项式商的 x 次数
    int ntail;                    ///< 多项式商的有效项数
} LAMB2_PF;

/** 按时间幂次保存部分分式展开结果 */
typedef struct {
    cplx_t pair[3][2][LAMB2_PF_TBAR_SIZE];            ///< pair[i][0/1][r]，i 为根索引，0/1 为 x/常数项，r 为 tbar 次数
    cplx_t tail[LAMB2_TAIL_SIZE][LAMB2_PF_TBAR_SIZE]; ///< tail[m][r]，m 为 x 次数，r 为 tbar 次数
    int ntail;                                     ///< 多项式商的有效项数
    int time_degree;                               ///< tbar 的最高次数
} LAMB2_PF_COEFF;

/** 一类积分项的全部部分分式展开结果 */
typedef struct {
    LAMB2_PF_COEFF M[2][3][3];              ///< M[xi][i][j]，xi=0/1 对应 U/V，i,j=0,1,2 分别为接收点和源点分量
    LAMB2_PF_COEFF dMs[2][3][3][3];         ///< dMs[xi][kp][i][j]，源点坐标导数
    LAMB2_PF_COEFF receiver_vertical[2][3]; ///< receiver_vertical[xi][j] 是接收点竖向导数的 [2][j] 分量
    LAMB2_PF_COEFF mixed[2][3][3][3][3];   ///< mixed[xi][k][kp][i][j]，混合导数对应的二次时间积分分子
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
static LAMB_X_POLY make_rayleigh_denominator(const cplx_t roots[3], const real_t kp2) {
    LAMB_X_POLY denominator = {0};
    denominator.c[0] = 16.0 * kp2;
    denominator.degree = 0;
    for (int i = 0; i < 3; ++i) {
        LAMB_X_POLY next = {0};
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
static void make_partial_fraction(const LAMB_X_POLY *numerator, const LAMB_X_POLY *denominator,
                                  const cplx_t roots[3], const real_t kp2, LAMB2_PF *pf) {
    LAMB_X_POLY quotient, remainder;
    grt_lamb_x_poly_divide(*numerator, *denominator, &quotient, &remainder);

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
static void make_partial_fraction_coeff(const LAMB_TBAR_X_POLY *numerator, const LAMB_X_POLY *denominator, const cplx_t roots[3], const real_t kp2,
                                        LAMB2_PF_COEFF *result) {
    if (numerator->tbar_degree >= LAMB2_PF_TBAR_SIZE) {
        GRTRaiseError("The tbar degree is too large in lamb2 partial fractions.\n");
    }
    memset(result, 0, sizeof(*result));
    result->time_degree = numerator->tbar_degree;
    for (int r = 0; r <= numerator->tbar_degree; ++r) {
        LAMB_X_POLY polynomial = {0};
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

/** 为一类积分项预计算全部基本和混合导数的部分分式展开结果 */
static void make_partial_fraction_set(const LAMB2_COEFF_SET *coeffs, const LAMB_X_POLY *denominator, const cplx_t roots[3], const real_t kp2,
                                      const bool use_angle_ratio, const real_t angle_ratio, const bool need_mixed, LAMB2_PF_SET *result) {
    memset(result, 0, sizeof(*result));

    for (int xi = 0; xi < 2; ++xi) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                make_partial_fraction_coeff(&coeffs->M[xi][i][j], denominator, roots, kp2, &result->M[xi][i][j]);
                for (int k = 0; k < 3; ++k) {
                    make_partial_fraction_coeff(
                        &coeffs->dMs[xi][k][i][j], denominator, roots, kp2, &result->dMs[xi][k][i][j]);
                    if (need_mixed) {
                        for (int kp = 0; kp < 3; ++kp) {
                            make_partial_fraction_coeff(&coeffs->d2M[xi][k][kp][i][j], denominator, roots, kp2,
                                                        &result->mixed[xi][k][kp][i][j]);
                        }
                    }
                }
            }
        }

        for (int j = 0; j < 3; ++j) {
            make_partial_fraction_coeff(&coeffs->receiver_vertical[xi][j], denominator, roots, kp2,
                                        &result->receiver_vertical[xi][j]);
        }

        if (use_angle_ratio) {
            scale_partial_fraction_coeff(&result->M[xi][0][2], angle_ratio, &result->M[xi][1][2]);
            scale_partial_fraction_coeff(&result->M[xi][2][0], angle_ratio, &result->M[xi][2][1]);
            scale_partial_fraction_coeff(&result->dMs[xi][0][2][2], angle_ratio, &result->dMs[xi][1][2][2]);
            scale_partial_fraction_coeff(&result->dMs[xi][2][0][2], angle_ratio, &result->dMs[xi][2][1][2]);
            scale_partial_fraction_coeff(&result->dMs[xi][2][2][0], angle_ratio, &result->dMs[xi][2][2][1]);
        }
    }

    for (int xi = 0; xi < 2; ++xi) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (result->M[xi][i][j].ntail > LAMB2_TAIL_SIZE) {
                    GRTRaiseError("The M partial-fraction tail is too large in lamb2.\n");
                }
                for (int k = 0; k < 3; ++k) {
                    if (result->dMs[xi][k][i][j].ntail > LAMB2_TAIL_SIZE) {
                        GRTRaiseError("The derivative partial-fraction tail is too large in lamb2.\n");
                    }
                }
                if (need_mixed) {
                    for (int receiver_direction = 0; receiver_direction < 3; ++receiver_direction) {
                        for (int source_direction = 0; source_direction < 3; ++source_direction) {
                            if (result->mixed[xi][receiver_direction][source_direction][i][j].ntail > LAMB2_TAIL_SIZE) {
                                GRTRaiseError("The mixed partial-fraction tail is too large in lamb2.\n");
                            }
                        }
                    }
                }
            }
        }
        for (int j = 0; j < 3; ++j) {
            if (result->receiver_vertical[xi][j].ntail > LAMB2_TAIL_SIZE) {
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
            if (i + 3 <= 7) {
                if (ctx->term == LAMB_BASIC_P_TERM) {
                    basis->tail[i] = I * grt_lamb_tail_V_P(i + 3, ctx, H);
                } else if (ctx->term == LAMB_BASIC_S_TERM) {
                    basis->tail[i] = I * grt_lamb_tail_V_S(i + 3, ctx, H);
                } else {
                    basis->tail[i] = I * grt_lamb_tail_V_SP(i + 3, ctx, H);
                }
            } else {
                basis->tail[i] = I * grt_lamb_tail_V8(ctx);
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

/** 评估三个坐标方向的导数分子 */
static void evaluate_derivatives(const LAMB2_PF_COEFF numerator[2][3][3][3], const LAMB2_EVAL_CTX *ctx, real_t result[3][3][3]) {
    for (int direction = 0; direction < 3; ++direction) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result[direction][i][j] = evaluate_component(&numerator[0][direction][i][j], &numerator[1][direction][i][j], ctx);
            }
        }
    }
}


/** 利用平移关系和自由表面关系恢复接收点导数 */
static void evaluate_receiver_derivatives(const LAMB2_PF_SET *coeffs, const LAMB2_EVAL_CTX *ctx,
                                          const real_t source[3][3][3], real_t result[3][3][3]) {
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
 * 计算一个 P、S 或 SP 项，并按 sign 累加到当前时间点的输出
 * F[i][j] 的索引分别表示接收点分量和源点分量
 * Fk_source[k'][i][j] 的索引依次表示源点坐标方向、接收点分量和源点分量
 * Fk_receiver[k][i][j] 的索引依次表示接收点坐标方向、接收点分量和源点分量
 * Fkk[k][k'][i][j] 的索引依次表示接收点方向、源点方向、接收点分量和源点分量
 * Fkk 是混合空间导数对应的二次时间积分项，最终需要连续求两次时间导数
 */
static void evaluate_lamb2_term(const real_t tbar, const real_t tbar2, const LAMB_BASIC_VARS *basic_vars,
                                const LAMB2_PF_SET *coeffs, const cplx_t roots[3], const LAMB_BASIC_TERM term, const bool use_U,
                                const real_t sign, const LAMB2_VARS *V, const bool need_mixed, real_t F[3][3],
                                real_t Fk_source[3][3][3], real_t Fk_receiver[3][3][3],
                                real_t Fkk[3][3][3][3]) {
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
        make_basis(LAMB2_TAIL_SIZE, roots, &ctx, false, &basis_U);
        basis_U_ptr = &basis_U;
    }
    make_basis(LAMB2_TAIL_SIZE, roots, &ctx, true, &basis_V);
    LAMB2_EVAL_CTX eval = {basis_U_ptr, &basis_V, tbar};

    real_t value[3][3];
    real_t dvalue_source[3][3][3];
    real_t dvalue_receiver[3][3][3];
    evaluate_base_matrix(coeffs->M, V, &eval, value);
    evaluate_derivatives(coeffs->dMs, &eval, dvalue_source);
    evaluate_receiver_derivatives(coeffs, &eval, dvalue_source, dvalue_receiver);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F[i][j] += sign * value[i][j];
            for (int k = 0; k < 3; ++k) {
                Fk_source[k][i][j] += sign * dvalue_source[k][i][j];
                Fk_receiver[k][i][j] += sign * dvalue_receiver[k][i][j];
                if (need_mixed) {
                    for (int kp = 0; kp < 3; ++kp) {
                        Fkk[k][kp][i][j] += sign * evaluate_component(
                            &coeffs->mixed[0][k][kp][i][j], &coeffs->mixed[1][k][kp][i][j], &eval);
                    }
                }
            }
        }
    }
}

/**
 * 计算一个时间点的 F、F_(,k')、F_(,k) 和 F_(,kk')
 *
 * 这里直接对应式 (7.4.1.4)，P 项、S 项和 SP 项共用同一套评估流程，
 * 仅由到时条件、根组、基本积分类型和累加符号区分；接收点导数的积分
 * 分子按接收点边界关系构造
 */
static void evaluate_lamb2_time(const real_t tbar, const LAMB2_VARS *V, const LAMB2_PF_COEFFICIENTS *coefficients,
                                const unsigned int phase_mask, real_t F[3][3],
                                real_t Fk_source[3][3][3], real_t Fk_receiver[3][3][3], const bool need_mixed,
                                real_t Fkk[3][3][3][3]) {
    memset(F, 0, sizeof(real_t) * 3 * 3);
    memset(Fk_source, 0, sizeof(real_t) * 3 * 3 * 3);
    memset(Fk_receiver, 0, sizeof(real_t) * 3 * 3 * 3);
    if (need_mixed) {
        memset(Fkk, 0, sizeof(real_t) * 3 * 3 * 3 * 3);
    }

    const real_t tbar2 = tbar * tbar;
    const LAMB_BASIC_VARS basic_vars = {V->k, V->k2, V->kp2, V->st, V->ct};

    if ((phase_mask & LAMB2_PHASE_P) != 0u && tbar > V->k) {
        evaluate_lamb2_term(tbar, tbar2, &basic_vars, &coefficients->P, V->rayleigh_roots, LAMB_BASIC_P_TERM, true, 1.0, V, need_mixed, F,
                            Fk_source, Fk_receiver, Fkk);
    }
    if ((phase_mask & LAMB2_PHASE_S) != 0u && tbar > 1.0) {
        evaluate_lamb2_term(tbar, tbar2, &basic_vars, &coefficients->S, V->shifted_rayleigh_roots, LAMB_BASIC_S_TERM, true, 1.0, V, need_mixed, F,
                            Fk_source, Fk_receiver, Fkk);
    }
    if ((phase_mask & (LAMB2_PHASE_SP | LAMB2_PHASE_PS)) != 0u && V->supercritical && tbar > V->tSP && tbar < 1.0) {
        evaluate_lamb2_term(tbar, tbar2, &basic_vars, &coefficients->S, V->shifted_rayleigh_roots, LAMB_BASIC_SP_TERM, false, -1.0, V, need_mixed, F,
                            Fk_source, Fk_receiver, Fkk);
    }
}

/**
 * 将地下源、地表接收的结果转为地表源、地下接收
 *
 * 对偶问题已按方位角 phi+pi 计算，这里只需转置 G_ij，
 * 并将源点导数与接收点导数对调后再转置 i,j，混合导数同时交换 k,k'
 */
static void apply_lamb2_surface_source_reciprocity(
    const int nt, real_t (*G)[3][3], real_t (*dG_source)[3][3][3], real_t (*dG_receiver)[3][3][3],
    real_t (*dG_mixed)[3][3][3][3])
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
        if (dG_mixed != NULL) {
            real_t mixed[3][3][3][3];
            for (int k = 0; k < 3; ++k) {
                for (int kp = 0; kp < 3; ++kp) {
                    for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            mixed[k][kp][i][j] = dG_mixed[n][k][kp][i][j];
                        }
                    }
                }
            }
            for (int k = 0; k < 3; ++k) {
                for (int kp = 0; kp < 3; ++kp) {
                    for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            dG_mixed[n][k][kp][i][j] = mixed[kp][k][j][i];
                        }
                    }
                }
            }
        }
    }
}


static real_t shift_lamb2_boundary(const real_t tbar, const LAMB2_VARS *V, const real_t tbar_eps) {
    // 与 lamb1 一致，精确命中波前时使用下一阶段的右侧值
    if (tbar == V->k || tbar == 1.0 || (V->supercritical && tbar == V->tSP)) {
        return tbar + tbar_eps;
    }
    return tbar;
}


void grt_compute_lamb2_travt(
    const real_t nu, const real_t R, const real_t depsrc, const real_t deprcv,
    real_t *tP, real_t *tSP)
{
    if (nu <= 0.0 || nu >= 0.5) {
        GRTRaiseError("poisson ratio (%lf) is out of bound.", nu);
    }
    if (R <= 0.0) {
        GRTRaiseError("The horizontal distance R should be positive in lamb2.\n");
    }
    if (depsrc < 0.0 || deprcv < 0.0 || (depsrc > 0.0) == (deprcv > 0.0)) {
        GRTRaiseError("lamb2 requires exactly one of source and receiver depths to be strictly positive.\n");
    }

    const real_t k = sqrt(0.5 * (1.0 - 2.0 * nu) / (1.0 - nu));
    const real_t buried_depth = depsrc > 0.0 ? depsrc : deprcv;
    const real_t theta = atan2(R, buried_depth);
    const real_t theta_c = asin(k);
    if (tP != NULL) {
        *tP = k;
    }
    if (tSP != NULL) {
        *tSP = theta > theta_c ? cos(theta - theta_c) : -1.0;
    }
}


void grt_solve_lamb2(
    const real_t nu, const real_t *ts, const int nt, const real_t R, const real_t depsrc, const real_t deprcv,
    const real_t azimuth, const char *phase_list, real_t (*G)[3][3], real_t (*dG_source)[3][3][3], real_t (*dG_receiver)[3][3][3],
    real_t (*dG_mixed)[3][3][3][3])
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
    if (depsrc < 0.0 || deprcv < 0.0) {
        GRTRaiseError("The source and receiver depths should be nonnegative in lamb2.\n");
    }
    const bool buried_source = depsrc > 0.0 && deprcv == 0.0;
    const bool surface_source = depsrc == 0.0 && deprcv > 0.0;
    if (!buried_source && !surface_source) {
        GRTRaiseError("lamb2 requires exactly one of source and receiver depths to be strictly positive, and the other to be zero.\n");
    }
    unsigned int phase_mask = grt_lamb_parse_phase_list(
        phase_list, LAMB2_PHASE_OPTIONS, sizeof(LAMB2_PHASE_OPTIONS) / sizeof(LAMB2_PHASE_OPTIONS[0]),
        "P, S, SP, PS");
    const unsigned int unavailable_phase = buried_source ? LAMB2_PHASE_PS : LAMB2_PHASE_SP;
    const bool has_unavailable_phase = (phase_mask & unavailable_phase) != 0u;
    if (phase_list != NULL && has_unavailable_phase) {
        GRTRaiseWarning(
            "Phase %s is unavailable for this lamb2 source-receiver geometry and is ignored; "
            "use %s for the converted phase.",
            buried_source ? "PS" : "SP", buried_source ? "SP" : "PS");
    }
    phase_mask &= ~unavailable_phase;
    if (phase_list != NULL && has_unavailable_phase && phase_mask == 0u) {
        GRTRaiseWarning(
            "No available Lamb phase remains for this source-receiver geometry; the output is all zeros.");
    }
    const real_t buried_depth = buried_source ? depsrc : deprcv;
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
    /* 允许时间轴从发震时刻之前开始，因果解在发震前保持为零 */
    for (int i = 0; i < nt; ++i) {
        if (i > 0 && ts[i] <= ts[i - 1]) {
            GRTRaiseError("The time series for lamb2 should be strictly increasing.\n");
        }
    }
    real_t tP;
    real_t tSP;
    grt_compute_lamb2_travt(nu, R, depsrc, deprcv, &tP, &tSP);
    LAMB2_VARS V = {0};
    V.nu = nu;
    V.k = tP;
    V.k2 = V.k * V.k;
    V.kp2 = 1.0 - V.k2;
    V.kp = sqrt(V.kp2);
    V.theta = atan2(R, buried_depth);
    V.phi = solve_azimuth * DEG1;
    V.st = sin(V.theta);
    V.ct = cos(V.theta);
    V.sf = sin(V.phi);
    V.cf = cos(V.phi);
    V.theta_c = asin(V.k);
    V.tSP = tSP;
    V.use_angle_ratio = fabs(V.cf) > LAMB2_RATIO_EPS;
    V.angle_ratio = V.use_angle_ratio ? V.sf / V.cf : 0.0;
    V.supercritical = V.tSP >= 0.0;
    const real_t tbar_eps = nt > 1 ? GRT_MIN(1e-8, (ts[1] - ts[0]) * 1e-5) : 1e-8;
    grt_rayleigh1_roots(V.nu, V.rayleigh_roots);
    for (int i = 0; i < 3; ++i) {
        V.shifted_rayleigh_roots[i] = V.rayleigh_roots[i] - V.kp2;
    }

    /* 末点若正好落在波前上会被右移，用略大的 tEnd 判断以免漏构造系数 */
    const real_t tEnd = ts[nt - 1] + tbar_eps;
    const bool need_P = (phase_mask & LAMB2_PHASE_P) != 0u && tEnd >= V.k;
    const bool need_S = (phase_mask & LAMB2_PHASE_S) != 0u ||
                        ((phase_mask & (LAMB2_PHASE_SP | LAMB2_PHASE_PS)) != 0u &&
                         (tEnd >= 1.0 || (V.supercritical && ts[0] < 1.0 && tEnd >= V.tSP)));
    const bool need_mixed = dG_mixed != NULL;
    /* 大型多项式系数工作区放在堆上，避免占用线程栈 */
    LAMB2_COEFF_SET *coefficients = GRT_SAFE_CALLOC(1, sizeof(*coefficients));
    LAMB2_PF_COEFFICIENTS *pf_coefficients = GRT_SAFE_CALLOC(1, sizeof(*pf_coefficients));
    if (need_P) {
        const LAMB_X_POLY denominator = make_rayleigh_denominator(V.rayleigh_roots, V.kp2);
        make_lamb2_P_coefficients(&V, need_mixed, coefficients);
        make_partial_fraction_set(coefficients, &denominator, V.rayleigh_roots, V.kp2, V.use_angle_ratio, V.angle_ratio,
                                  need_mixed, &pf_coefficients->P);
    }
    if (need_S) {
        const LAMB_X_POLY denominator = make_rayleigh_denominator(V.shifted_rayleigh_roots, V.kp2);
        make_lamb2_S_coefficients(&V, need_mixed, coefficients);
        make_partial_fraction_set(coefficients, &denominator, V.shifted_rayleigh_roots, V.kp2, V.use_angle_ratio, V.angle_ratio,
                                  need_mixed, &pf_coefficients->S);
    }
    GRT_SAFE_FREE_PTR(coefficients);

    const bool isprint = G == NULL && dG_source == NULL && dG_receiver == NULL && dG_mixed == NULL;
    real_t(*F)[3][3] = G != NULL ? G : GRT_SAFE_CALLOC((size_t)nt, sizeof(*F));
    real_t(*Fk_source)[3][3][3] = GRT_SAFE_CALLOC((size_t)nt, sizeof(*Fk_source));
    real_t(*Fk_receiver)[3][3][3] = GRT_SAFE_CALLOC((size_t)nt, sizeof(*Fk_receiver));
    real_t(*Fkk)[3][3][3][3] = need_mixed ? GRT_SAFE_CALLOC((size_t)nt, sizeof(*Fkk)) : NULL;
    real_t(*dG_source_tmp)[3][3][3] = dG_source != NULL ? dG_source : GRT_SAFE_CALLOC((size_t)nt, sizeof(*dG_source_tmp));
    real_t(*dG_receiver_tmp)[3][3][3] = dG_receiver != NULL ? dG_receiver : GRT_SAFE_CALLOC((size_t)nt, sizeof(*dG_receiver_tmp));

    for (int i = 0; i < nt; ++i) {
        real_t tbar = shift_lamb2_boundary(ts[i], &V, tbar_eps);
        evaluate_lamb2_time(tbar, &V, pf_coefficients, phase_mask, F[i], Fk_source[i], Fk_receiver[i], need_mixed,
                            need_mixed ? Fkk[i] : NULL);
    }
    grt_lamb_differentiate_Fk(ts, nt, Fk_source, dG_source_tmp);
    grt_lamb_differentiate_Fk(ts, nt, Fk_receiver, dG_receiver_tmp);
    if (need_mixed) {
        grt_lamb_differentiate_Fkk(ts, nt, Fkk, dG_mixed);
    }
    if (surface_source) {
        apply_lamb2_surface_source_reciprocity(
            nt, F, dG_source_tmp, dG_receiver_tmp, dG_mixed);
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
    GRT_SAFE_FREE_PTR(Fkk);
    if (dG_source == NULL) {
        GRT_SAFE_FREE_PTR(dG_source_tmp);
    }
    if (dG_receiver == NULL) {
        GRT_SAFE_FREE_PTR(dG_receiver_tmp);
    }
}
