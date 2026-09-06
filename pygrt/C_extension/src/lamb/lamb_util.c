/**
 * @file   lamb_util.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    一些使用广义闭合解求解 Lamb 问题过程中可能用到的辅助函数
 */

#include <float.h>
#include <string.h>

#include "grt/lamb/lamb_util.h"


static bool is_derivative_suboption(const char *text)
{
    return text[0] == '+' && (text[1] == 's' || text[1] == 'r');
}


static char *copy_path(const char *start, const size_t length)
{
    char *path = calloc(length + 1, sizeof(*path));
    if (path == NULL) {
        GRTRaiseError("Cannot allocate a Lamb derivative output path.\n");
    }
    memcpy(path, start, length);
    return path;
}


void grt_lamb_parse_derivative_paths(const char *argument, char **source_path, char **receiver_path)
{
    if (argument == NULL || argument[0] == '\0') {
        GRTRaiseError("The -S argument should contain +s<path> and/or +r<path>.\n");
    }

    const char *cursor = argument;
    while (cursor[0] != '\0') {
        if (!is_derivative_suboption(cursor)) {
            GRTRaiseError("The -S argument should contain +s<path> and/or +r<path>.\n");
        }

        char **target = cursor[1] == 's' ? source_path : receiver_path;
        if (*target != NULL) {
            GRTRaiseError("The -S argument contains a duplicated derivative output suboption.\n");
        }

        const char *path_start = cursor + 2;
        const char *next = path_start;
        while (next[0] != '\0' && !is_derivative_suboption(next)) {
            ++next;
        }
        if (next == path_start) {
            GRTRaiseError("The -S argument contains an empty derivative output path.\n");
        }

        *target = copy_path(path_start, (size_t)(next - path_start));
        cursor = next;
    }
}


/** 创建常数多项式 */
LAMB_POLY grt_lamb_poly_const(const cplx_t value)
{
    LAMB_POLY result = {0};
    result.c[0] = value;
    result.degree = 0;
    return result;
}

LAMB_POLY grt_lamb_poly_x(void)
{
    LAMB_POLY result = {0};
    result.c[1] = 1.0;
    result.degree = 1;
    return result;
}

void grt_lamb_poly_trim(LAMB_POLY *poly)
{
    while (poly->degree > 0 && cabs(poly->c[poly->degree]) < LAMB_POLY_EPS) {
        poly->c[poly->degree] = 0.0;
        --poly->degree;
    }
}

LAMB_POLY grt_lamb_poly_mul(const LAMB_POLY a, const LAMB_POLY b)
{
    LAMB_POLY result = {0};
    if (a.degree + b.degree >= LAMB_POLY_SIZE) {
        GRTRaiseError("The polynomial degree is too large in the Lamb utilities.\n");
    }
    result.degree = a.degree + b.degree;
    for (int i = 0; i <= a.degree; ++i) {
        for (int j = 0; j <= b.degree; ++j) {
            result.c[i + j] += a.c[i] * b.c[j];
        }
    }
    grt_lamb_poly_trim(&result);
    return result;
}

LAMB_POLY grt_lamb_poly_factor(const cplx_t root, const bool plus)
{
    LAMB_POLY result = {0};
    result.c[0] = plus ? root : -root;
    result.c[2] = 1.0;
    result.degree = 2;
    return result;
}

cplx_t grt_lamb_poly_eval(const LAMB_POLY *poly, const cplx_t x)
{
    cplx_t result = 0.0;
    for (int i = poly->degree; i >= 0; --i) {
        result = result * x + poly->c[i];
    }
    return result;
}

void grt_lamb_poly_divide(
    const LAMB_POLY numerator, const LAMB_POLY denominator,
    LAMB_POLY *quotient, LAMB_POLY *remainder)
{
    *quotient = grt_lamb_poly_const(0.0);
    *remainder = numerator;
    grt_lamb_poly_trim(remainder);

    if (denominator.degree <= 0 && cabs(denominator.c[0]) == 0.0) {
        GRTRaiseError("The polynomial denominator is zero in the Lamb utilities.\n");
    }

    while (remainder->degree >= denominator.degree &&
           !(remainder->degree == 0 && cabs(remainder->c[0]) < LAMB_POLY_EPS)) {
        int offset = remainder->degree - denominator.degree;
        cplx_t factor = remainder->c[remainder->degree] / denominator.c[denominator.degree];
        quotient->c[offset] += factor;
        quotient->degree = GRT_MAX(quotient->degree, offset);
        for (int j = 0; j <= denominator.degree; ++j) {
            remainder->c[j + offset] -= factor * denominator.c[j];
        }
        grt_lamb_poly_trim(remainder);
    }
    grt_lamb_poly_trim(quotient);
}

bool grt_lamb_is_real(const cplx_t value)
{
    return fabs(cimag(value)) <= 1e-12 * (1.0 + fabs(creal(value)));
}

real_t grt_lamb_positive_sqrt(const real_t value, const char *name)
{
    if (value < 0.0) {
        if (value > -1e-10) {
            return 0.0;
        }
        GRTRaiseError("The discriminant %s is negative in the Lamb calculation: %e.\n", name, value);
    }
    return sqrt(value);
}

real_t grt_lamb_clamp_elliptic_parameter(const real_t value, const real_t tolerance, const char *name)
{
    if (value <= 0.0) {
        if (value < -tolerance) {
            GRTRaiseError("The elliptic parameter is out of range in %s: %e.\n", name, value);
        }
        return DBL_EPSILON;
    }
    if (value >= 1.0) {
        if (value > 1.0 + tolerance) {
            GRTRaiseError("The elliptic parameter is out of range in %s: %e.\n", name, value);
        }
        return 1.0 - 16.0 * DBL_EPSILON;
    }
    return value;
}

/**
 * 求解一元三次方程的根， \f$ x^3 + ax^2 + bx + c = 0 \f$
 *
 * @param[in]      a       系数 a
 * @param[in]      b       系数 b
 * @param[in]      c       系数 c
 * @param[out]     roots    三个复根
 */
void grt_lamb_cubic_roots(const real_t a, const real_t b, const real_t c, cplx_t roots[3])
{
    real_t Q = (a * a - 3.0 * b) / 9.0;
    real_t R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;
    real_t Q3 = Q * Q * Q;
    real_t R2 = R * R;

    roots[0] = roots[1] = roots[2] = 0.0;
    if (Q > 0.0 && Q3 > R2) {
        real_t ratio = R / sqrt(Q3);
        ratio = GRT_MAX(-1.0, GRT_MIN(1.0, ratio));
        real_t angle = acos(ratio);
        roots[0] = -2.0 * sqrt(Q) * cos(angle / 3.0) - a / 3.0;
        roots[1] = -2.0 * sqrt(Q) * cos((angle - 2.0 * PI) / 3.0) - a / 3.0;
        roots[2] = -2.0 * sqrt(Q) * cos((angle + 2.0 * PI) / 3.0) - a / 3.0;
    } else {
        real_t discriminant = R2 - Q3;
        real_t A = pow(fabs(R) + sqrt(GRT_MAX(0.0, discriminant)), 1.0 / 3.0);
        A = R > 0.0 ? -A : A;
        real_t B = A == 0.0 ? 0.0 : Q / A;
        roots[0] = -0.5 * (A + B) - a / 3.0 + I * sqrt(3.0) / 2.0 * (A - B);
        roots[1] = -0.5 * (A + B) - a / 3.0 - I * sqrt(3.0) / 2.0 * (A - B);
        roots[2] = A + B - a / 3.0;
    }
}

cplx_t grt_lamb_eval_time_coeff(const cplx_t *coefficient, const int degree, const real_t t)
{
    cplx_t result = 0.0;
    for (int i = degree; i >= 0; --i) {
        result = result * t + coefficient[i];
    }
    return result;
}

real_t grt_lamb_derivative_three_points(
    const real_t x0, const real_t x1, const real_t x2,
    const real_t f0, const real_t f1, const real_t f2, const real_t x)
{
    real_t result = f0 * (2.0 * x - x1 - x2) / ((x0 - x1) * (x0 - x2));
    result += f1 * (2.0 * x - x0 - x2) / ((x1 - x0) * (x1 - x2));
    result += f2 * (2.0 * x - x0 - x1) / ((x2 - x0) * (x2 - x1));
    return result;
}

void grt_lamb_differentiate_Fk(
    const real_t *ts, const int nt, const real_t (*Fk)[3][3][3],
    real_t (*dG)[3][3][3])
{
    if (nt == 1) {
        memset(dG, 0, sizeof(real_t) * 3 * 3 * 3);
        return;
    }

    for (int k = 0; k < 3; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (nt == 2) {
                    real_t value = (Fk[1][k][i][j] - Fk[0][k][i][j]) / (ts[1] - ts[0]);
                    dG[0][k][i][j] = value;
                    dG[1][k][i][j] = value;
                } else {
                    dG[0][k][i][j] = grt_lamb_derivative_three_points(
                        ts[0], ts[1], ts[2], Fk[0][k][i][j], Fk[1][k][i][j], Fk[2][k][i][j], ts[0]);
                    for (int n = 1; n < nt - 1; ++n) {
                        dG[n][k][i][j] = grt_lamb_derivative_three_points(
                            ts[n - 1], ts[n], ts[n + 1], Fk[n - 1][k][i][j], Fk[n][k][i][j], Fk[n + 1][k][i][j], ts[n]);
                    }
                    dG[nt - 1][k][i][j] = grt_lamb_derivative_three_points(
                        ts[nt - 3], ts[nt - 2], ts[nt - 1], Fk[nt - 3][k][i][j], Fk[nt - 2][k][i][j], Fk[nt - 1][k][i][j], ts[nt - 1]);
                }
            }
        }
    }
}

static void print_component_header(FILE *fp, const char *name)
{
    fprintf(fp, "%14s", name);
}


void grt_lamb_print_green_series(FILE *fp, const real_t *ts, const int nt, const real_t (*G)[3][3])
{
    fprintf(fp, "#%13s", "tbar");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            char name[16];
            snprintf(name, sizeof(name), "G%d%d", i + 1, j + 1);
            print_component_header(fp, name);
        }
    }
    fprintf(fp, "\n");

    for (int n = 0; n < nt; ++n) {
        fprintf(fp, "%14.6e", ts[n]);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                fprintf(fp, "%14.6e", G[n][i][j]);
            }
        }
        fprintf(fp, "\n");
    }
}


void grt_lamb_print_derivative_series(FILE *fp, const real_t *ts, const int nt, const real_t (*dG)[3][3][3], const bool source)
{
    fprintf(fp, "#%13s", "tbar");
    for (int k = 0; k < 3; ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                char name[16];
                snprintf(name, sizeof(name), "G%d%d,%d%s", i + 1, j + 1, k + 1, source ? "'" : "");
                print_component_header(fp, name);
            }
        }
    }
    fprintf(fp, "\n");

    for (int n = 0; n < nt; ++n) {
        fprintf(fp, "%14.6e", ts[n]);
        for (int k = 0; k < 3; ++k) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    fprintf(fp, "%14.6e", dG[n][k][i][j]);
                }
            }
        }
        fprintf(fp, "\n");
    }
}

void grt_rayleigh1_roots(real_t nu, cplx_t y3[3])
{
    real_t a, b, c;
    real_t nu2, nu3, nu4;
    nu2 = nu*nu;
    nu3 = nu2*nu;
    nu4 = nu3*nu;
    real_t snu = 1.0 - nu;
    real_t snu2 = snu*snu;
    real_t snu3 = snu2*snu;
    a = -0.5 * (2.0*nu2 + 1.0)/snu;
    b = 0.25 * (4.0*nu3 - 4.0*nu2 + 4.0*nu - 1.0)/snu2;
    c = -0.125*nu4/snu3;
    grt_lamb_cubic_roots(a, b, c, y3);
}

void grt_rayleigh2_roots(real_t m, cplx_t y3[3])
{
    real_t a, b, c;
    a = 0.5*(2.0*m - 3.0)/(1.0 - m);
    b = 0.5/(1.0 - m);
    c = - 0.0625/(1 - m);
    grt_lamb_cubic_roots(a, b, c, y3);
}

cplx_t grt_evalpoly2(const cplx_t *C, const int n, const cplx_t y, const int offset)
{
    cplx_t res = 0.0;
    cplx_t p = 1.0;
    for(int i=0; i<=n; ++i){
        res += C[2*i+offset] * p;
        p *= y;
    }
    return res;
}
