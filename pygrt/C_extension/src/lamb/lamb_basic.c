/**
 * @file   lamb_basic.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-09
 *
 *    第二类和第三类 Lamb 问题共用的基本积分
 */

#include <math.h>

#include "grt/common/checkerror.h"
#include "grt/lamb/elliptic.h"
#include "grt/lamb/lamb_basic.h"
#include "grt/lamb/lamb_util.h"

static int interval_indicator(const real_t x, const real_t left, const real_t right) {
    return x > left && x < right;
}

static int F1_P(const cplx_t z, const real_t z2sq) {
    return grt_lamb_is_real(z) && interval_indicator(creal(z), -z2sq, 0.0);
}

static int F1_S(const cplx_t z, const real_t z1, const real_t z2sq) {
    real_t left = GRT_MAX(-z1, -z2sq);
    return grt_lamb_is_real(z) && interval_indicator(creal(z), left, 0.0);
}

static int F2(const cplx_t z, const real_t z2sq) {
    real_t x = creal(z);
    real_t y = cimag(z);
    if (!interval_indicator(x, -z2sq, 0.0) || y <= 0.0) {
        return 0;
    }
    real_t ymax = sqrt((z2sq + x) * (1.0 - x));
    return y < ymax;
}

/** 式 (7.2.11)、(7.3.7) 和 (7.3.18) 中的分式线性变换参数 */
void grt_lamb_make_context_P(const real_t tbar, const real_t tbar2, const LAMB_BASIC_VARS *V, LAMB_BASIC_CONTEXT *ctx) {
    ctx->term = LAMB_BASIC_P_TERM;
    ctx->tbar = tbar;
    ctx->kp2 = V->kp2;
    ctx->m = tbar * V->ct;
    ctx->n = sqrt(tbar2 - V->k2) * V->st;

    real_t d = (tbar2 + V->k2 * V->ct * V->ct - 1.0) / ctx->m;
    real_t discriminant = d * d + 4.0 * V->kp2;
    real_t root_gap = grt_lamb_positive_sqrt(discriminant, "P");
    ctx->xi1 = 0.5 * (d + root_gap);
    ctx->xi2 = 0.5 * (d - root_gap);
    ctx->z1sq = ctx->xi1 / (-ctx->xi2);
    ctx->z2sq = (ctx->xi1 - ctx->m) / (ctx->m - ctx->xi2);
    ctx->z1 = sqrt(ctx->z1sq);
    ctx->z2 = sqrt(ctx->z2sq);
    ctx->m_elliptic = grt_lamb_clamp_elliptic_parameter(ctx->z2sq / ctx->z1sq, 1e-9, "Lamb basic integrals");
    ctx->c_main = 1.0 / sqrt(ctx->xi1 * (ctx->m - ctx->xi2));
    ctx->c_residue = 1.0 / sqrt(ctx->xi2 * (ctx->xi2 - ctx->m));
    ctx->c1 = 0.0;
    ctx->c2 = 0.0;
}

void grt_lamb_make_context_S(const real_t tbar, const real_t tbar2, const LAMB_BASIC_VARS *V, LAMB_BASIC_CONTEXT *ctx) {
    ctx->term = LAMB_BASIC_S_TERM;
    ctx->tbar = tbar;
    ctx->kp2 = V->kp2;
    ctx->m = tbar * V->ct;
    ctx->n = sqrt(tbar2 - 1.0) * V->st;

    real_t d = (tbar2 - V->k2 + V->ct * V->ct) / ctx->m;
    real_t discriminant = d * d - 4.0 * V->kp2;
    real_t root_gap = grt_lamb_positive_sqrt(discriminant, "S");
    ctx->xi1 = 0.5 * (d + root_gap);
    ctx->xi2 = 0.5 * (d - root_gap);
    ctx->z1sq = ctx->xi1 / ctx->xi2;
    ctx->z2sq = (ctx->xi1 - ctx->m) / (ctx->m - ctx->xi2);
    ctx->z1 = sqrt(ctx->z1sq);
    ctx->z2 = sqrt(ctx->z2sq);
    ctx->m_elliptic = grt_lamb_clamp_elliptic_parameter(ctx->z2sq / (ctx->z1sq + ctx->z2sq), 1e-9, "Lamb basic integrals");
    ctx->c_main = 1.0 / sqrt(ctx->m * (ctx->xi1 - ctx->xi2));
    ctx->c_residue = 1.0 / sqrt(ctx->xi2 * (ctx->m - ctx->xi2));
    ctx->c1 = 0.0;
    ctx->c2 = 0.0;
}

void grt_lamb_make_context_SP(const real_t tbar, const real_t tbar2, const LAMB_BASIC_VARS *V, LAMB_BASIC_CONTEXT *ctx) {
    ctx->term = LAMB_BASIC_SP_TERM;
    ctx->tbar = tbar;
    ctx->kp2 = V->kp2;
    ctx->m = tbar * V->ct;
    ctx->n = sqrt(1.0 - tbar2) * V->st;

    real_t d = (tbar2 - V->k2 + V->ct * V->ct) / ctx->m;
    real_t discriminant = d * d - 4.0 * V->kp2;
    real_t root_gap = grt_lamb_positive_sqrt(discriminant, "S-P");
    ctx->xi1 = 0.5 * (d + root_gap);
    ctx->xi2 = 0.5 * (d - root_gap);
    ctx->z1sq = ctx->xi1 / ctx->xi2;
    ctx->z2sq = (ctx->xi1 - ctx->m) / (ctx->xi2 - ctx->m);
    ctx->z1 = sqrt(ctx->z1sq);
    ctx->z2 = sqrt(ctx->z2sq);
    ctx->m_elliptic = grt_lamb_clamp_elliptic_parameter(
        ctx->m * (ctx->xi1 - ctx->xi2) / (ctx->xi2 * (ctx->xi1 - ctx->m)), 1e-9, "Lamb basic integrals");
    ctx->c_main = 1.0 / sqrt(ctx->xi2 * (ctx->xi1 - ctx->m));
    ctx->c_residue = 0.0;
    ctx->c1 = ctx->z2sq - 1.0;
    ctx->c2 = ctx->z2sq - ctx->z1sq;
}

cplx_t grt_lamb_basic_U(const int number, const cplx_t c, const LAMB_BASIC_CONTEXT *ctx) {
    if (ctx->term == LAMB_BASIC_SP_TERM) {
        /* 式 (7.3.4) 和 7.3.4.2：这些 U 项对最终虚部没有贡献 */
        return 0.0;
    }

    if (number == 1) {
        real_t n = ctx->n;
        cplx_t m1 = sqrt(c) - I * ctx->m;
        cplx_t m2 = sqrt(c) + I * ctx->m;
        cplx_t A1 = sqrt((m1 - n) / (m1 + n));
        cplx_t A2 = sqrt((m2 + n) / (m2 - n));
        return A1 / (m1 - n) * atan(A1) - A2 / (m2 + n) * atan(A2);
    }
    if (number == 2) {
        real_t n = ctx->n;
        cplx_t m1 = sqrt(c) - I * ctx->m;
        cplx_t m2 = sqrt(c) + I * ctx->m;
        cplx_t A1 = sqrt((m1 - n) / (m1 + n));
        cplx_t A2 = sqrt((m2 + n) / (m2 - n));
        return I / sqrt(c) * (A1 / (m1 - n) * atan(A1) + A2 / (m2 + n) * atan(A2));
    }

    real_t m = ctx->m;
    real_t n = ctx->n;
    if (number == 3) {
        return I * HALFPI;
    }
    if (number == 4) {
        return I * HALFPI * m;
    }
    if (number == 5) {
        return I * PI / 4.0 * (2.0 * m * m - n * n);
    }
    if (number == 6) {
        return I * PI / 4.0 * m * (2.0 * m * m - 3.0 * n * n);
    }

    GRTRaiseError("Wrong U basic-integral number in the Lamb basic integrals: %d.\n", number);
}

/** 同时计算同一分母对应的两个 U 基本积分 */
void grt_lamb_make_U_pair(const cplx_t c, const LAMB_BASIC_CONTEXT *ctx, cplx_t result[2]) {
    if (ctx->term == LAMB_BASIC_SP_TERM) {
        result[0] = 0.0;
        result[1] = 0.0;
        return;
    }

    real_t n = ctx->n;
    cplx_t sqrt_c = sqrt(c);
    cplx_t m1 = sqrt_c - I * ctx->m;
    cplx_t m2 = sqrt_c + I * ctx->m;
    cplx_t A1 = sqrt((m1 - n) / (m1 + n));
    cplx_t A2 = sqrt((m2 + n) / (m2 - n));
    cplx_t first = A1 / (m1 - n) * atan(A1);
    cplx_t second = A2 / (m2 + n) * atan(A2);

    result[0] = first - second;
    result[1] = I / sqrt_c * (first + second);
}

void grt_lamb_calculate_H(const LAMB_BASIC_CONTEXT *ctx, const real_t K, real_t H[5]) {
    real_t m = ctx->m_elliptic;
    real_t a;
    if (ctx->term == LAMB_BASIC_P_TERM) {
        a = 1.0 / ctx->z2sq;
    } else if (ctx->term == LAMB_BASIC_S_TERM) {
        a = -(ctx->z2sq + 1.0) / ctx->z2sq;
    } else {
        a = -ctx->c1 / ctx->c2;
    }

    real_t E = grt_ellipticE(m);
    real_t H0 = K;
    real_t H1 = grt_ellipticPi(-1.0 / a, m) / a;
    real_t Hminus1 = ((a * m + 1.0) * K - E) / m;
    real_t gamma1 = 3.0 * m * a * a + 2.0 * a * (m + 1.0) + 1.0;
    real_t gamma2 = m + 1.0 + 3.0 * m * a;
    real_t gamma3 = a * (a + 1.0) * (a * m + 1.0);
    real_t H2 = (gamma1 * H1 - m * Hminus1) / (2.0 * gamma3);
    real_t H3 = (3.0 * gamma1 * H2 - 2.0 * gamma2 * H1 + m * H0) / (4.0 * gamma3);
    real_t H4 = (5.0 * gamma1 * H3 - 4.0 * gamma2 * H2 + 3.0 * m * H1) / (6.0 * gamma3);

    H[0] = H0;
    H[1] = H1;
    H[2] = H2;
    H[3] = H3;
    H[4] = H4;
}

real_t grt_lamb_tail_V_P(const int number, const LAMB_BASIC_CONTEXT *ctx, const real_t H[5]) {
    /* 式 (7.2.29a)-(7.2.29e) */
    real_t xi1 = ctx->xi1;
    real_t xi2 = ctx->xi2;
    real_t z2 = ctx->z2;
    real_t z2m2 = 1.0 / (z2 * z2);
    real_t z2m4 = z2m2 * z2m2;
    real_t z2m6 = z2m4 * z2m2;
    real_t z2m8 = z2m6 * z2m2;
    real_t beta11 = xi1 - xi2;
    real_t beta12 = xi1 - 2.0 * xi2;
    real_t beta13 = xi1 - 3.0 * xi2;
    real_t beta15 = xi1 - 5.0 * xi2;
    real_t beta35 = 3.0 * xi1 - 5.0 * xi2;
    real_t result;

    if (number == 3) {
        result = ctx->c_main * H[0];
    } else if (number == 4) {
        result = ctx->c_main * (xi2 * H[0] + beta11 * z2m2 * H[1]);
    } else if (number == 5) {
        result = ctx->c_main * (xi2 * xi2 * H[0] - beta11 * beta13 * z2m2 * H[1] + 2.0 * beta11 * beta11 * z2m4 * H[2]);
    } else if (number == 6) {
        result = ctx->c_main * (xi2 * xi2 * xi2 * H[0] - 3.0 * xi2 * beta11 * beta12 * z2m2 * H[1] - 3.0 * beta11 * beta11 * beta13 * z2m4 * H[2] +
                                4.0 * beta11 * beta11 * beta11 * z2m6 * H[3]);
    } else if (number == 7) {
        result = ctx->c_main * (xi2 * xi2 * xi2 * xi2 * H[0] - 2.0 * xi2 * xi2 * beta11 * beta35 * z2m2 * H[1] +
                                beta11 * beta11 * (beta15 * beta15 + 4.0 * ctx->kp2) * z2m4 * H[2] -
                                8.0 * beta11 * beta11 * beta11 * beta13 * z2m6 * H[3] + 8.0 * beta11 * beta11 * beta11 * beta11 * z2m8 * H[4]);
    } else {
        GRTRaiseError("Wrong P-wave V basic-integral number in the Lamb basic integrals: %d.\n", number);
    }
    return result;
}

real_t grt_lamb_tail_V_S(const int number, const LAMB_BASIC_CONTEXT *ctx, const real_t H[5]) {
    /* 式 (7.3.17a)-(7.3.17e) */
    real_t xi1 = ctx->xi1;
    real_t xi2 = ctx->xi2;
    real_t z2 = ctx->z2;
    real_t z2m2 = 1.0 / (z2 * z2);
    real_t z2m4 = z2m2 * z2m2;
    real_t z2m6 = z2m4 * z2m2;
    real_t z2m8 = z2m6 * z2m2;
    real_t beta11 = xi1 - xi2;
    real_t beta12 = xi1 - 2.0 * xi2;
    real_t beta13 = xi1 - 3.0 * xi2;
    real_t beta15 = xi1 - 5.0 * xi2;
    real_t beta35 = 3.0 * xi1 - 5.0 * xi2;
    real_t result;

    if (number == 3) {
        result = ctx->c_main * H[0];
    } else if (number == 4) {
        result = ctx->c_main * (xi2 * H[0] - beta11 * z2m2 * H[1]);
    } else if (number == 5) {
        result = ctx->c_main * (xi2 * xi2 * H[0] + beta11 * beta13 * z2m2 * H[1] + 2.0 * beta11 * beta11 * z2m4 * H[2]);
    } else if (number == 6) {
        result = ctx->c_main * (xi2 * xi2 * xi2 * H[0] + 3.0 * xi2 * beta11 * beta12 * z2m2 * H[1] - 3.0 * beta11 * beta11 * beta13 * z2m4 * H[2] -
                                4.0 * beta11 * beta11 * beta11 * z2m6 * H[3]);
    } else if (number == 7) {
        result = ctx->c_main * (xi2 * xi2 * xi2 * xi2 * H[0] + 2.0 * xi2 * xi2 * beta11 * beta35 * z2m2 * H[1] +
                                beta11 * beta11 * (beta15 * beta15 - 4.0 * ctx->kp2) * z2m4 * H[2] +
                                8.0 * beta11 * beta11 * beta11 * beta13 * z2m6 * H[3] + 8.0 * beta11 * beta11 * beta11 * beta11 * z2m8 * H[4]);
    } else {
        GRTRaiseError("Wrong S-wave V basic-integral number in the Lamb basic integrals: %d.\n", number);
    }
    return result;
}

real_t grt_lamb_tail_V_SP(const int number, const LAMB_BASIC_CONTEXT *ctx, const real_t H[5]) {
    /* 式 (7.3.22a)-(7.3.22e)，M1-M4 按式 (E.7)、(E.9a)-(E.9c) */
    real_t xi1 = ctx->xi1;
    real_t xi2 = ctx->xi2;
    real_t z2 = ctx->z2;
    real_t beta11 = xi1 - xi2;
    real_t beta12 = xi1 - 2.0 * xi2;
    real_t beta13 = xi1 - 3.0 * xi2;
    real_t beta14 = xi1 - 4.0 * xi2;
    real_t beta15 = xi1 - 5.0 * xi2;
    real_t beta17 = xi1 - 7.0 * xi2;
    real_t beta35 = 3.0 * xi1 - 5.0 * xi2;
    real_t c1 = ctx->c1;
    real_t c2 = ctx->c2;
    real_t c1m1 = 1.0 / c2;
    real_t c1m2 = c1m1 * c1m1;
    real_t c1m3 = c1m2 * c1m1;
    real_t c1m4 = c1m3 * c1m1;
    real_t c1c1mc2 = c1 * (c1 - c2);
    real_t sqrt_c1c1mc2 = sqrt(c1c1mc2);
    real_t c1c1mc2_2 = c1c1mc2 * c1c1mc2;
    real_t c1c1mc2_3 = c1c1mc2_2 * c1c1mc2;
    real_t M1 = PI / (2.0 * sqrt_c1c1mc2);
    real_t M2 = PI * (2.0 * c1 - c2) / (4.0 * c1c1mc2 * sqrt_c1c1mc2);
    real_t M3 = PI * (8.0 * c1 * c1 - 8.0 * c1 * c2 + 3.0 * c2 * c2) / (16.0 * c1c1mc2_2 * sqrt_c1c1mc2);
    real_t M4 = PI * (2.0 * c1 - c2) * (8.0 * c1 * c1 - 8.0 * c1 * c2 + 5.0 * c2 * c2) / (32.0 * c1c1mc2_3 * sqrt_c1c1mc2);
    real_t result;

    if (number == 3) {
        result = -ctx->c_main * H[0];
    } else if (number == 4) {
        result = -ctx->c_main * (xi2 * H[0] + beta11 * c1m1 * H[1] + beta11 * z2 * M1);
    } else if (number == 5) {
        result = -ctx->c_main * (xi2 * xi2 * H[0] - beta11 * beta13 * c1m1 * H[1] + 2.0 * beta11 * beta11 * c1m2 * H[2] +
                                 2.0 * beta11 * z2 * (xi2 * M1 - beta11 * M2));
    } else if (number == 6) {
        result = -ctx->c_main * (xi2 * xi2 * xi2 * H[0] - 3.0 * xi2 * beta11 * beta12 * c1m1 * H[1] - 3.0 * beta11 * beta11 * beta13 * c1m2 * H[2] +
                                 4.0 * beta11 * beta11 * beta11 * c1m3 * H[3] +
                                 beta11 * z2 * (3.0 * xi2 * xi2 * M1 + beta11 * beta17 * M2 + 4.0 * beta11 * beta11 * M3));
    } else if (number == 7) {
        result = -ctx->c_main *
                 (xi2 * xi2 * xi2 * xi2 * H[0] - 2.0 * xi2 * xi2 * beta11 * beta35 * c1m1 * H[1] +
                  beta11 * beta11 * (beta15 * beta15 - 4.0 * ctx->kp2) * c1m2 * H[2] - 8.0 * beta11 * beta11 * beta11 * beta13 * c1m3 * H[3] +
                  8.0 * beta11 * beta11 * beta11 * beta11 * c1m4 * H[4] +
                  4.0 * beta11 * z2 *
                      (xi2 * xi2 * xi2 * M1 + xi2 * beta11 * beta14 * M2 - beta11 * beta11 * beta15 * M3 - 2.0 * beta11 * beta11 * beta11 * M4));
    } else {
        GRTRaiseError("Wrong S-P-wave V basic-integral number in the Lamb basic integrals: %d.\n", number);
    }
    return result;
}

void grt_lamb_make_V_aux(const cplx_t c, const LAMB_BASIC_CONTEXT *ctx, LAMB_BASIC_V_AUX *aux) {
    /* 式 (7.2.15)、式 (7.3.12) 和式 (7.3.20) 中的 z0、s 及 h 变量 */
    cplx_t root = sqrt(-c);
    cplx_t gap = ctx->xi1 - ctx->xi2;
    cplx_t h2psq;

    aux->h1p = c + ctx->xi1 * ctx->xi1;
    aux->h2p = c + ctx->xi2 * ctx->xi2;
    aux->h1m = c - ctx->xi1 * ctx->xi1;
    aux->h2m = c - ctx->xi2 * ctx->xi2;
    h2psq = aux->h2p * aux->h2p;
    aux->D = gap * gap - 2.0 * (c + ctx->xi1 * ctx->xi2);
    aux->s[0] = -(c + ctx->xi1 * ctx->xi2 + gap * root) * (c + ctx->xi1 * ctx->xi2 + gap * root) / h2psq;
    aux->s[1] = -(c + ctx->xi1 * ctx->xi2 - gap * root) * (c + ctx->xi1 * ctx->xi2 - gap * root) / h2psq;
    aux->z0[0] = (ctx->xi1 + root) / (ctx->xi2 + root);
    aux->z0[1] = (ctx->xi1 - root) / (ctx->xi2 - root);
}

static void residue_P_pair(const int sign, const cplx_t c, const LAMB_BASIC_CONTEXT *ctx, const LAMB_BASIC_V_AUX *aux, cplx_t result[2]) {
    int other = 1 - sign;
    cplx_t z = aux->z0[sign];
    cplx_t denominator = (c + ctx->xi2 * ctx->xi2) * (z - aux->z0[other]);
    denominator *= sqrt(z * z + ctx->z1sq) * sqrt(z * z + ctx->z2sq);
    cplx_t scale = ctx->c_residue / denominator;
    result[0] = scale * (z - 1.0) * (ctx->xi2 * z - ctx->xi1);
    result[1] = scale * (z - 1.0) * (z - 1.0);
}

static void residue_S_pair(const int sign, const cplx_t c, const LAMB_BASIC_CONTEXT *ctx, const LAMB_BASIC_V_AUX *aux, cplx_t result[2]) {
    int other = 1 - sign;
    cplx_t z = aux->z0[sign];
    cplx_t denominator = (c + ctx->xi2 * ctx->xi2) * (z - aux->z0[other]);
    denominator *= sqrt(ctx->z1sq - z * z) * sqrt(z * z + ctx->z2sq);
    cplx_t scale = ctx->c_residue / denominator;
    result[0] = scale * (z - 1.0) * (ctx->xi2 * z - ctx->xi1);
    result[1] = scale * (z - 1.0) * (z - 1.0);
}

void grt_lamb_make_V_P_pair(const cplx_t c, const LAMB_BASIC_CONTEXT *ctx, const real_t K, cplx_t result[2]) {
    /* 式 (7.2.24)、式 (7.2.26)-(7.2.28) */
    LAMB_BASIC_V_AUX aux;
    grt_lamb_make_V_aux(c, ctx, &aux);

    cplx_t h2psq = aux.h2p * aux.h2p;
    cplx_t zeta0 = ctx->xi2 / aux.h2p;
    cplx_t eta0 = 1.0 / aux.h2p;
    cplx_t zeta[2];
    cplx_t eta[2];
    cplx_t zeta_bar[2];
    cplx_t eta_bar[2];
    cplx_t branch[2];
    cplx_t residue[2][2];
    cplx_t alpha11, alpha12, alpha21, alpha22;
    real_t m = ctx->m_elliptic;

    for (int i = 0; i < 2; ++i) {
        int other = 1 - i;
        cplx_t s = aux.s[i];
        cplx_t denominator = s * (aux.s[other] - s) * h2psq;
        zeta[i] = (ctx->xi1 + s * ctx->xi2) * (s * aux.h2p + aux.h1p) / denominator;
        zeta_bar[i] = (ctx->xi1 - ctx->xi2) * (s * aux.h2m + aux.h1m) / ((aux.s[other] - s) * h2psq);
        eta[i] = (aux.h2p * s * s - aux.D * s + aux.h1p) / denominator;
        eta_bar[i] = 2.0 * (ctx->xi1 - ctx->xi2) * (s * ctx->xi2 + ctx->xi1) / ((aux.s[other] - s) * h2psq);
    }

    cplx_t n1 = ctx->z2sq / aux.s[0];
    cplx_t n2 = ctx->z2sq / aux.s[1];
    cplx_t Pi1 = grt_ellipticPi(n1, m);
    cplx_t Pi2 = grt_ellipticPi(n2, m);
    cplx_t elliptic[2] = {
        zeta0 * K + zeta[0] * Pi1 + zeta[1] * Pi2,
        eta0 * K + eta[0] * Pi1 + eta[1] * Pi2,
    };

    alpha11 = sqrt(ctx->z1sq - aux.s[0]);
    alpha12 = sqrt(ctx->z1sq - aux.s[1]);
    alpha21 = sqrt(ctx->z2sq - aux.s[0]);
    alpha22 = sqrt(ctx->z2sq - aux.s[1]);
    cplx_t atan1 = atan(I * ctx->z2 * alpha11 / (ctx->z1 * alpha21));
    cplx_t atan2 = atan(I * ctx->z2 * alpha12 / (ctx->z1 * alpha22));
    branch[0] = ctx->z1 * zeta_bar[0] / (alpha11 * alpha21) * atan1 + ctx->z1 * zeta_bar[1] / (alpha12 * alpha22) * atan2;
    branch[1] = ctx->z1 * eta_bar[0] / (alpha11 * alpha21) * atan1 + ctx->z1 * eta_bar[1] / (alpha12 * alpha22) * atan2;
    for (int sign = 0; sign < 2; ++sign) {
        residue_P_pair(sign, c, ctx, &aux, residue[sign]);
    }

    if (grt_lamb_is_real(c)) {
        real_t value[2] = {
            ctx->c_main * creal(elliptic[0]),
            ctx->c_main * creal(elliptic[1]),
        };
        for (int sign = 0; sign < 2; ++sign) {
            real_t d[2] = {
                creal(residue[sign][0]),
                creal(residue[sign][1]),
            };
            if (creal(c) < 0.0) {
                real_t factor = PI * F1_P(aux.z0[sign], ctx->z2sq);
                value[0] -= factor * d[0];
                value[1] -= factor * d[1];
            } else if (creal(c) > 0.0) {
                real_t factor = 2.0 * PI * F2(aux.z0[sign], ctx->z2sq);
                value[0] -= factor * d[0];
                value[1] -= factor * d[1];
            }
        }
        result[0] = I * value[0];
        result[1] = I * value[1];
        return;
    }

    result[0] = I * ctx->c_main * (elliptic[0] - branch[0]);
    result[1] = I * ctx->c_main * (elliptic[1] + branch[1]);
    real_t factor0 = F2(aux.z0[0], ctx->z2sq);
    real_t factor1 = F2(aux.z0[1], ctx->z2sq);
    result[0] -= 2.0 * PI * I * (factor0 * residue[0][0] + factor1 * residue[1][0]);
    result[1] -= 2.0 * PI * I * (factor0 * residue[0][1] + factor1 * residue[1][1]);
}

void grt_lamb_make_V_S_pair(const cplx_t c, const LAMB_BASIC_CONTEXT *ctx, const real_t K, cplx_t result[2]) {
    /* 式 (7.3.13)-(7.3.16)，这里的 V^S 已包含 S1-H(theta-theta_c)S2 */
    LAMB_BASIC_V_AUX aux;
    grt_lamb_make_V_aux(c, ctx, &aux);

    cplx_t h2psq = aux.h2p * aux.h2p;
    cplx_t zeta0 = ctx->xi2 / aux.h2p;
    cplx_t eta0 = 1.0 / aux.h2p;
    cplx_t zeta[2];
    cplx_t eta[2];
    cplx_t zeta_bar[2];
    cplx_t eta_bar[2];
    real_t m = ctx->m_elliptic;

    for (int i = 0; i < 2; ++i) {
        int other = 1 - i;
        cplx_t s = aux.s[i];
        cplx_t denominator = (aux.s[other] - s) * h2psq * (s - ctx->z2sq);
        zeta[i] = (ctx->xi1 + s * ctx->xi2) * (s * aux.h2p + aux.h1p) / denominator;
        zeta_bar[i] = (ctx->xi1 - ctx->xi2) * (s * aux.h2m + aux.h1m) / ((aux.s[other] - s) * h2psq);
        eta[i] = (aux.h2p * s * s - aux.D * s + aux.h1p) / denominator;
        eta_bar[i] = 2.0 * (ctx->xi1 - ctx->xi2) * (s * ctx->xi2 + ctx->xi1) / ((aux.s[other] - s) * h2psq);
    }

    cplx_t n1 = ctx->z2sq / (ctx->z2sq - aux.s[0]);
    cplx_t n2 = ctx->z2sq / (ctx->z2sq - aux.s[1]);
    cplx_t Pi1 = grt_ellipticPi(n1, m);
    cplx_t Pi2 = grt_ellipticPi(n2, m);
    cplx_t elliptic[2] = {
        zeta0 * K + zeta[0] * Pi1 + zeta[1] * Pi2,
        eta0 * K + eta[0] * Pi1 + eta[1] * Pi2,
    };

    cplx_t alpha11 = sqrt(aux.s[0] + ctx->z1sq);
    cplx_t alpha12 = sqrt(aux.s[0] - ctx->z2sq);
    cplx_t alpha21 = sqrt(aux.s[1] + ctx->z1sq);
    cplx_t alpha22 = sqrt(aux.s[1] - ctx->z2sq);
    cplx_t atan1 = atan(ctx->z2 * alpha11 / (ctx->z1 * alpha12));
    cplx_t atan2 = atan(ctx->z2 * alpha21 / (ctx->z1 * alpha22));
    cplx_t branch[2] = {
        zeta_bar[0] / (alpha11 * alpha12) * atan1 + zeta_bar[1] / (alpha21 * alpha22) * atan2,
        eta_bar[0] / (alpha11 * alpha12) * atan1 + eta_bar[1] / (alpha21 * alpha22) * atan2,
    };
    cplx_t residue[2][2];
    for (int sign = 0; sign < 2; ++sign) {
        residue_S_pair(sign, c, ctx, &aux, residue[sign]);
    }

    if (grt_lamb_is_real(c)) {
        real_t value[2] = {
            ctx->c_main * creal(elliptic[0]),
            ctx->c_main * creal(elliptic[1]),
        };
        for (int sign = 0; sign < 2; ++sign) {
            real_t d[2] = {
                creal(residue[sign][0]),
                creal(residue[sign][1]),
            };
            if (creal(c) < 0.0) {
                real_t factor = PI * F1_S(aux.z0[sign], ctx->z1, ctx->z2sq);
                value[0] -= factor * d[0];
                value[1] -= factor * d[1];
            } else if (creal(c) > 0.0) {
                real_t factor = 2.0 * PI * F2(aux.z0[sign], ctx->z2sq);
                value[0] -= factor * d[0];
                value[1] -= factor * d[1];
            }
        }
        result[0] = I * value[0];
        result[1] = I * value[1];
        return;
    }

    real_t branch_factor = 1.0 / sqrt(ctx->xi2 * (ctx->m - ctx->xi2));
    result[0] = I * ctx->c_main * elliptic[0] - branch[0] * branch_factor;
    result[1] = I * ctx->c_main * elliptic[1] + branch[1] * branch_factor;
    real_t factor0 = F2(aux.z0[0], ctx->z2sq);
    real_t factor1 = F2(aux.z0[1], ctx->z2sq);
    result[0] -= 2.0 * PI * I * (factor0 * residue[0][0] + factor1 * residue[1][0]);
    result[1] -= 2.0 * PI * I * (factor0 * residue[0][1] + factor1 * residue[1][1]);
}

void grt_lamb_make_V_SP_pair(const cplx_t c, const LAMB_BASIC_CONTEXT *ctx, const real_t K, cplx_t result[2]) {
    /* 式 (7.3.20) 和式 (7.3.21) */
    LAMB_BASIC_V_AUX aux;
    grt_lamb_make_V_aux(c, ctx, &aux);
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

    lambda_p = (ctx->xi2 * z0p2 - ctx->xi1) * ((c + ctx->xi2 * ctx->xi2) * z0p2 - (c + ctx->xi1 * ctx->xi1)) / denominator;
    lambda_m = (ctx->xi2 * z0m2 - ctx->xi1) * ((c + ctx->xi2 * ctx->xi2) * z0m2 - (c + ctx->xi1 * ctx->xi1)) / denominator_reverse;
    lambda_prime_p = (ctx->xi1 - ctx->xi2) * ((c - ctx->xi2 * ctx->xi2) * z0p2 - (c - ctx->xi1 * ctx->xi1)) / denominator;
    lambda_prime_m = (ctx->xi1 - ctx->xi2) * ((c - ctx->xi2 * ctx->xi2) * z0m2 - (c - ctx->xi1 * ctx->xi1)) / denominator_reverse;

    lambda_bar_p = ((c + ctx->xi2 * ctx->xi2) * z0p2 * z0p2 +
                    ((ctx->xi1 - ctx->xi2) * (ctx->xi1 - ctx->xi2) - 2.0 * (c + ctx->xi1 * ctx->xi2)) * z0p2 + (c + ctx->xi1 * ctx->xi1)) /
                   denominator;
    lambda_bar_m = ((c + ctx->xi2 * ctx->xi2) * z0m2 * z0m2 +
                    ((ctx->xi1 - ctx->xi2) * (ctx->xi1 - ctx->xi2) - 2.0 * (c + ctx->xi1 * ctx->xi2)) * z0m2 + (c + ctx->xi1 * ctx->xi1)) /
                   denominator_reverse;
    lambda_bar_prime_p = 2.0 * (ctx->xi1 - ctx->xi2) * (ctx->xi2 * z0p2 - ctx->xi1) / denominator;
    lambda_bar_prime_m = 2.0 * (ctx->xi1 - ctx->xi2) * (ctx->xi2 * z0m2 - ctx->xi1) / denominator_reverse;

    mu1p = sqrt(ctx->z1sq - z0p2);
    mu2p = sqrt(ctx->z2sq - z0p2);
    mu1m = sqrt(ctx->z1sq - z0m2);
    mu2m = sqrt(ctx->z2sq - z0m2);
    Pi_p = grt_ellipticPi((ctx->z2sq - ctx->z1sq) / (ctx->z2sq - z0p2), ctx->m_elliptic);
    Pi_m = grt_ellipticPi((ctx->z2sq - ctx->z1sq) / (ctx->z2sq - z0m2), ctx->m_elliptic);

    cplx_t bracket1 = zeta0 * K + lambda_p * Pi_p / (mu2p * mu2p) + lambda_m * Pi_m / (mu2m * mu2m) +
                      PI * lambda_prime_p * ctx->z2 / (2.0 * mu1p * mu2p) + PI * lambda_prime_m * ctx->z2 / (2.0 * mu1m * mu2m);
    cplx_t bracket2 = eta0 * K + lambda_bar_p * Pi_p / (mu2p * mu2p) + lambda_bar_m * Pi_m / (mu2m * mu2m) -
                      PI * lambda_bar_prime_p * ctx->z2 / (2.0 * mu1p * mu2p) - PI * lambda_bar_prime_m * ctx->z2 / (2.0 * mu1m * mu2m);
    result[0] = -I * ctx->c_main * bracket1;
    result[1] = -I * ctx->c_main * bracket2;
}
