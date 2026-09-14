/**
 * @file   lamb_poly.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-09
 *
 *    Lamb 问题闭合解中两类多项式的运算
 */

#include <string.h>

#include "grt/lamb/lamb_poly.h"


LAMB_X_POLY grt_lamb_x_poly_const(const cplx_t value)
{
    LAMB_X_POLY result = {0};
    result.c[0] = value;
    result.degree = 0;
    return result;
}


LAMB_X_POLY grt_lamb_x_poly_x(void)
{
    LAMB_X_POLY result = {0};
    result.c[1] = 1.0;
    result.degree = 1;
    return result;
}


void grt_lamb_x_poly_trim(LAMB_X_POLY *poly)
{
    while (poly->degree > 0 && cabs(poly->c[poly->degree]) < LAMB_X_POLY_EPS) {
        poly->c[poly->degree] = 0.0;
        --poly->degree;
    }
}


LAMB_X_POLY grt_lamb_x_poly_mul(const LAMB_X_POLY a, const LAMB_X_POLY b)
{
    LAMB_X_POLY result = {0};
    if (a.degree + b.degree >= LAMB_POLY_X_SIZE) {
        GRTRaiseError("The polynomial degree is too large in the Lamb utilities.\n");
    }
    result.degree = a.degree + b.degree;
    for (int i = 0; i <= a.degree; ++i) {
        for (int j = 0; j <= b.degree; ++j) {
            result.c[i + j] += a.c[i] * b.c[j];
        }
    }
    grt_lamb_x_poly_trim(&result);
    return result;
}


LAMB_X_POLY grt_lamb_x_poly_factor(const cplx_t root, const bool plus)
{
    LAMB_X_POLY result = {0};
    result.c[0] = plus ? root : -root;
    result.c[2] = 1.0;
    result.degree = 2;
    return result;
}


cplx_t grt_lamb_x_poly_eval(const LAMB_X_POLY *poly, const cplx_t x)
{
    cplx_t result = 0.0;
    for (int i = poly->degree; i >= 0; --i) {
        result = result * x + poly->c[i];
    }
    return result;
}


void grt_lamb_x_poly_divide(
    const LAMB_X_POLY numerator, const LAMB_X_POLY denominator,
    LAMB_X_POLY *quotient, LAMB_X_POLY *remainder)
{
    *quotient = grt_lamb_x_poly_const(0.0);
    *remainder = numerator;
    grt_lamb_x_poly_trim(remainder);

    if (denominator.degree <= 0 && cabs(denominator.c[0]) == 0.0) {
        GRTRaiseError("The polynomial denominator is zero in the Lamb utilities.\n");
    }

    while (remainder->degree >= denominator.degree &&
           !(remainder->degree == 0 && cabs(remainder->c[0]) < LAMB_X_POLY_EPS)) {
        int offset = remainder->degree - denominator.degree;
        cplx_t factor = remainder->c[remainder->degree] / denominator.c[denominator.degree];
        quotient->c[offset] += factor;
        quotient->degree = GRT_MAX(quotient->degree, offset);
        for (int j = 0; j <= denominator.degree; ++j) {
            remainder->c[j + offset] -= factor * denominator.c[j];
        }
        grt_lamb_x_poly_trim(remainder);
    }
    grt_lamb_x_poly_trim(quotient);
}


void grt_lamb_tbar_x_poly_zero(LAMB_TBAR_X_POLY *poly)
{
    memset(poly, 0, sizeof(*poly));
}


void grt_lamb_tbar_x_poly_const(LAMB_TBAR_X_POLY *poly, const real_t value)
{
    grt_lamb_tbar_x_poly_zero(poly);
    poly->c[0][0] = value;
}


void grt_lamb_tbar_x_poly_x(LAMB_TBAR_X_POLY *poly)
{
    grt_lamb_tbar_x_poly_zero(poly);
    poly->degree = 1;
    poly->c[0][1] = 1.0;
}


void grt_lamb_tbar_x_poly_tbar(LAMB_TBAR_X_POLY *poly)
{
    grt_lamb_tbar_x_poly_zero(poly);
    poly->tbar_degree = 1;
    poly->c[1][0] = 1.0;
}


void grt_lamb_tbar_x_poly_trim(LAMB_TBAR_X_POLY *poly)
{
    int degree = 0;
    int tbar_degree = 0;

    for (int r = 0; r < LAMB_POLY_TBAR_SIZE; ++r) {
        for (int m = 0; m < LAMB_POLY_X_SIZE; ++m) {
            if (poly->c[r][m] != 0.0) {
                degree = GRT_MAX(degree, m);
                tbar_degree = GRT_MAX(tbar_degree, r);
            }
        }
    }
    poly->degree = degree;
    poly->tbar_degree = tbar_degree;
}


void grt_lamb_tbar_x_poly_add(const LAMB_TBAR_X_POLY *left, const LAMB_TBAR_X_POLY *right,
                             LAMB_TBAR_X_POLY *result)
{
    const int degree = GRT_MAX(left->degree, right->degree);
    const int tbar_degree = GRT_MAX(left->tbar_degree, right->tbar_degree);
    if (degree >= LAMB_POLY_X_SIZE || tbar_degree >= LAMB_POLY_TBAR_SIZE) {
        GRTRaiseError("The tbar-x polynomial sum is too large.\n");
    }

    LAMB_TBAR_X_POLY sum = {0};
    sum.degree = degree;
    sum.tbar_degree = tbar_degree;
    for (int r = 0; r <= tbar_degree; ++r) {
        for (int m = 0; m <= degree; ++m) {
            sum.c[r][m] = left->c[r][m] + right->c[r][m];
        }
    }
    grt_lamb_tbar_x_poly_trim(&sum);
    *result = sum;
}


void grt_lamb_tbar_x_poly_add_scaled(const LAMB_TBAR_X_POLY *base, const real_t scale,
                                    const LAMB_TBAR_X_POLY *term, LAMB_TBAR_X_POLY *result)
{
    LAMB_TBAR_X_POLY scaled = *term;
    for (int r = 0; r <= term->tbar_degree; ++r) {
        for (int m = 0; m <= term->degree; ++m) {
            scaled.c[r][m] *= scale;
        }
    }
    grt_lamb_tbar_x_poly_add(base, &scaled, result);
}


void grt_lamb_tbar_x_poly_scale(const LAMB_TBAR_X_POLY *source, const real_t scale,
                               LAMB_TBAR_X_POLY *result)
{
    *result = *source;
    for (int r = 0; r <= source->tbar_degree; ++r) {
        for (int m = 0; m <= source->degree; ++m) {
            result->c[r][m] *= scale;
        }
    }
    grt_lamb_tbar_x_poly_trim(result);
}


void grt_lamb_tbar_x_poly_mul(const LAMB_TBAR_X_POLY *left, const LAMB_TBAR_X_POLY *right,
                             LAMB_TBAR_X_POLY *result)
{
    const int degree = left->degree + right->degree;
    const int tbar_degree = left->tbar_degree + right->tbar_degree;
    if (degree >= LAMB_POLY_X_SIZE || tbar_degree >= LAMB_POLY_TBAR_SIZE) {
        GRTRaiseError("The tbar-x polynomial product is too large.\n");
    }

    LAMB_TBAR_X_POLY product = {0};
    product.degree = degree;
    product.tbar_degree = tbar_degree;
    for (int r1 = 0; r1 <= left->tbar_degree; ++r1) {
        for (int m1 = 0; m1 <= left->degree; ++m1) {
            for (int r2 = 0; r2 <= right->tbar_degree; ++r2) {
                for (int m2 = 0; m2 <= right->degree; ++m2) {
                    product.c[r1 + r2][m1 + m2] += left->c[r1][m1] * right->c[r2][m2];
                }
            }
        }
    }
    grt_lamb_tbar_x_poly_trim(&product);
    *result = product;
}


void grt_lamb_tbar_x_poly_product2_scaled(const LAMB_TBAR_X_POLY *first, const LAMB_TBAR_X_POLY *second,
                                         const real_t scale, LAMB_TBAR_X_POLY *result)
{
    grt_lamb_tbar_x_poly_mul(first, second, result);
    grt_lamb_tbar_x_poly_scale(result, scale, result);
}


void grt_lamb_tbar_x_poly_product3(const LAMB_TBAR_X_POLY *first, const LAMB_TBAR_X_POLY *second,
                                  const LAMB_TBAR_X_POLY *third, LAMB_TBAR_X_POLY *result)
{
    LAMB_TBAR_X_POLY product;
    grt_lamb_tbar_x_poly_mul(first, second, &product);
    grt_lamb_tbar_x_poly_mul(&product, third, result);
}


void grt_lamb_tbar_x_poly_product3_scaled(const LAMB_TBAR_X_POLY *first, const LAMB_TBAR_X_POLY *second,
                                         const LAMB_TBAR_X_POLY *third, const real_t scale,
                                         LAMB_TBAR_X_POLY *result)
{
    grt_lamb_tbar_x_poly_product3(first, second, third, result);
    grt_lamb_tbar_x_poly_scale(result, scale, result);
}


void grt_lamb_tbar_x_poly_product4_scaled(const LAMB_TBAR_X_POLY *first, const LAMB_TBAR_X_POLY *second,
                                         const LAMB_TBAR_X_POLY *third, const LAMB_TBAR_X_POLY *fourth,
                                         const real_t scale, LAMB_TBAR_X_POLY *result)
{
    LAMB_TBAR_X_POLY product;
    grt_lamb_tbar_x_poly_product3(first, second, third, &product);
    grt_lamb_tbar_x_poly_mul(&product, fourth, result);
    grt_lamb_tbar_x_poly_scale(result, scale, result);
}


void grt_lamb_tbar_x_poly_product5_scaled(const LAMB_TBAR_X_POLY *first, const LAMB_TBAR_X_POLY *second,
                                         const LAMB_TBAR_X_POLY *third, const LAMB_TBAR_X_POLY *fourth,
                                         const LAMB_TBAR_X_POLY *fifth, const real_t scale,
                                         LAMB_TBAR_X_POLY *result)
{
    LAMB_TBAR_X_POLY product;
    grt_lamb_tbar_x_poly_product4_scaled(first, second, third, fourth, 1.0, &product);
    grt_lamb_tbar_x_poly_mul(&product, fifth, result);
    grt_lamb_tbar_x_poly_scale(result, scale, result);
}


void grt_lamb_tbar_x_poly_product6_scaled(const LAMB_TBAR_X_POLY *first, const LAMB_TBAR_X_POLY *second,
                                         const LAMB_TBAR_X_POLY *third, const LAMB_TBAR_X_POLY *fourth,
                                         const LAMB_TBAR_X_POLY *fifth, const LAMB_TBAR_X_POLY *sixth,
                                         const real_t scale, LAMB_TBAR_X_POLY *result)
{
    LAMB_TBAR_X_POLY product;
    grt_lamb_tbar_x_poly_product5_scaled(first, second, third, fourth, fifth, 1.0, &product);
    grt_lamb_tbar_x_poly_mul(&product, sixth, result);
    grt_lamb_tbar_x_poly_scale(result, scale, result);
}
