/**
 * @file   lamb_util.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    一些使用广义闭合解求解 Lamb 问题过程中可能用到的辅助函数
 */

#include <float.h>
#include <ctype.h>
#include <string.h>

#include "grt/lamb/lamb_util.h"


unsigned int grt_lamb_parse_phase_list(
    const char *phase_list, const GRT_LAMB_PHASE_OPTION *options,
    const size_t option_count, const char *available_names)
{
    unsigned int all_mask = 0u;
    for (size_t i = 0; i < option_count; ++i) {
        all_mask |= options[i].bit;
    }
    if (phase_list == NULL) {
        return all_mask;
    }

    char *copy = strdup(phase_list);
    if (copy == NULL) {
        GRTRaiseError("Unable to allocate the Lamb phase selection.\n");
    }

    unsigned int mask = 0u;
    char *token = strtok(copy, ",");
    while (token != NULL) {
        while (isspace((unsigned char)*token)) {
            ++token;
        }
        char *end = token + strlen(token);
        while (end > token && isspace((unsigned char)end[-1])) {
            --end;
            *end = '\0';
        }

        if (token[0] == '\0') {
            GRTRaiseWarning(
                "An empty Lamb phase name is ignored. Available phases: %s.", available_names);
        } else {
            unsigned int phase = 0u;
            for (size_t i = 0; i < option_count; ++i) {
                if (strcasecmp(token, options[i].name) == 0) {
                    phase = options[i].bit;
                    break;
                }
            }
            if (phase == 0u) {
                GRTRaiseWarning(
                    "Unsupported Lamb phase %s is ignored. Available phases: %s.",
                    token, available_names);
            } else if ((mask & phase) != 0u) {
                GRTRaiseWarning("Duplicated Lamb phase %s is ignored; it is recorded only once.", token);
            } else {
                mask |= phase;
            }
        }

        token = strtok(NULL, ",");
    }

    GRT_SAFE_FREE_PTR(copy);
    if (mask == 0u) {
        GRTRaiseWarning(
            "No valid Lamb phase remains; the output is all zeros. Available phases: %s.", available_names);
    }
    return mask;
}


static bool is_derivative_suboption(const char *text, const bool allow_mixed)
{
    return text[0] == '+' && (text[1] == 's' || text[1] == 'r' || (allow_mixed && text[1] == 'm'));
}


static char *copy_path(const char *start, const size_t length)
{
    char *path = GRT_SAFE_CALLOC(length + 1, sizeof(*path));
    memcpy(path, start, length);
    return path;
}


static void parse_derivative_paths(
    const char *argument, char **source_path, char **receiver_path, char **mixed_path, const bool allow_mixed)
{
    const char *option_names = allow_mixed ? "+s<path>, +r<path> or +m<path>" : "+s<path> and/or +r<path>";
    if (argument == NULL || argument[0] == '\0') {
        GRTRaiseError("The -S argument should contain %s.\n", option_names);
    }

    const char *cursor = argument;
    while (cursor[0] != '\0') {
        if (!is_derivative_suboption(cursor, allow_mixed)) {
            GRTRaiseError("The -S argument should contain %s.\n", option_names);
        }

        char **target = cursor[1] == 's' ? source_path : cursor[1] == 'r' ? receiver_path : mixed_path;
        if (target == NULL) {
            GRTRaiseError("The derivative output path pointer is NULL.\n");
        }
        if (*target != NULL) {
            GRTRaiseError("The -S argument contains a duplicated derivative output suboption.\n");
        }

        const char *path_start = cursor + 2;
        const char *next = path_start;
        while (next[0] != '\0' && !is_derivative_suboption(next, allow_mixed)) {
            ++next;
        }
        if (next == path_start) {
            GRTRaiseError("The -S argument contains an empty derivative output path.\n");
        }

        *target = copy_path(path_start, (size_t)(next - path_start));
        cursor = next;
    }

    if (source_path != NULL && receiver_path != NULL && *source_path != NULL && *receiver_path != NULL &&
        strcmp(*source_path, *receiver_path) == 0) {
        GRTRaiseError("The source and receiver derivative paths must be different.\n");
    }
    if (allow_mixed && source_path != NULL && mixed_path != NULL && *source_path != NULL && *mixed_path != NULL &&
        strcmp(*source_path, *mixed_path) == 0) {
        GRTRaiseError("The source and mixed derivative paths must be different.\n");
    }
    if (allow_mixed && receiver_path != NULL && mixed_path != NULL && *receiver_path != NULL && *mixed_path != NULL &&
        strcmp(*receiver_path, *mixed_path) == 0) {
        GRTRaiseError("The receiver and mixed derivative paths must be different.\n");
    }
}


void grt_lamb_parse_derivative_paths(const char *argument, char **source_path, char **receiver_path)
{
    parse_derivative_paths(argument, source_path, receiver_path, NULL, false);
}


void grt_lamb_parse_derivative_paths_with_mixed(
    const char *argument, char **source_path, char **receiver_path, char **mixed_path)
{
    parse_derivative_paths(argument, source_path, receiver_path, mixed_path, true);
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


void grt_lamb_differentiate_Fkk(
    const real_t *ts, const int nt, const real_t (*Fkk)[3][3][3][3],
    real_t (*dG)[3][3][3][3])
{
    real_t (*input)[3][3][3] = GRT_SAFE_CALLOC((size_t)nt, sizeof(*input));
    real_t (*first)[3][3][3] = GRT_SAFE_CALLOC((size_t)nt, sizeof(*first));
    real_t (*second)[3][3][3] = GRT_SAFE_CALLOC((size_t)nt, sizeof(*second));

    for (int k = 0; k < 3; ++k) {
        for (int n = 0; n < nt; ++n) {
            memcpy(input[n], Fkk[n][k], sizeof(input[n]));
        }
        grt_lamb_differentiate_Fk(ts, nt, input, first);
        grt_lamb_differentiate_Fk(ts, nt, first, second);
        for (int n = 0; n < nt; ++n) {
            for (int kp = 0; kp < 3; ++kp) {
                memcpy(dG[n][k][kp], second[n][kp], sizeof(second[n][kp]));
            }
        }
    }

    GRT_SAFE_FREE_PTR(input);
    GRT_SAFE_FREE_PTR(first);
    GRT_SAFE_FREE_PTR(second);
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


void grt_lamb_print_mixed_derivative_series(FILE *fp, const real_t *ts, const int nt, const real_t (*dG)[3][3][3][3])
{
    fprintf(fp, "#%13s", "tbar");
    for (int k = 0; k < 3; ++k) {
        for (int kp = 0; kp < 3; ++kp) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    char name[20];
                    snprintf(name, sizeof(name), "G%d%d,%d,%d'", i + 1, j + 1, k + 1, kp + 1);
                    print_component_header(fp, name);
                }
            }
        }
    }
    fprintf(fp, "\n");

    for (int n = 0; n < nt; ++n) {
        fprintf(fp, "%14.6e", ts[n]);
        for (int k = 0; k < 3; ++k) {
            for (int kp = 0; kp < 3; ++kp) {
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        fprintf(fp, "%14.6e", dG[n][k][kp][i][j]);
                    }
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
