/**
 * @file   qromb.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-09
 *
 *    C implementation of the qromb algorithm described in Numerical Recipes
 *    in FORTRAN 77, section 4.3.
 */

#include "grt/lamb/qromb.h"


#define LAMB_QROMB_JMAX 20
#define LAMB_QROMB_K 5


static real_t polint(
    const real_t *xa, const real_t *ya, const int n, const real_t x, real_t *dy)
{
    real_t c[LAMB_QROMB_K];
    real_t d[LAMB_QROMB_K];
    int ns = 0;
    real_t dif = fabs(x - xa[0]);

    for(int i = 0; i < n; ++i){
        real_t dift = fabs(x - xa[i]);
        if(dift < dif){
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }

    real_t y = ya[ns];
    --ns;
    for(int m = 1; m < n; ++m){
        for(int i = 0; i < n - m; ++i){
            real_t ho = xa[i] - x;
            real_t hp = xa[i + m] - x;
            real_t w = c[i + 1] - d[i];
            real_t den = ho - hp;
            if(den == 0.0){
                GRTRaiseError("Repeated abscissas in Lamb Romberg interpolation.\n");
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        if(2 * ns < n - m){
            *dy = c[ns + 1];
        } else {
            *dy = d[ns];
            --ns;
        }
        y += *dy;
    }
    return y;
}


real_t grt_lamb_qromb(
    GRT_LAMB_INTEGRAND func, real_t a, real_t b, real_t eps, void *userdata)
{
    if(func == NULL){
        GRTRaiseError("The Lamb Romberg integrand should not be NULL.\n");
    }
    if(eps <= 0.0){
        GRTRaiseError("The Lamb Romberg tolerance should be positive.\n");
    }
    if(a == b){
        return 0.0;
    }

    real_t s[LAMB_QROMB_JMAX];
    real_t h[LAMB_QROMB_JMAX];
    h[0] = 1.0;

    size_t it = 1;
    for(int j = 0; j < LAMB_QROMB_JMAX; ++j){
        if(j == 0){
            s[j] = 0.5 * (b - a) * (func(a, userdata) + func(b, userdata));
        } else {
            real_t del = (b - a) / (real_t)it;
            real_t x = a + 0.5 * del;
            real_t sum = 0.0;
            for(size_t i = 0; i < it; ++i){
                sum += func(x, userdata);
                x += del;
            }
            s[j] = 0.5 * (s[j - 1] + (b - a) * sum / (real_t)it);
        }

        if(j >= LAMB_QROMB_K - 1){
            int first = j - (LAMB_QROMB_K - 1);
            real_t xh[LAMB_QROMB_K];
            real_t ys[LAMB_QROMB_K];
            for(int i = 0; i < LAMB_QROMB_K; ++i){
                xh[i] = h[first + i];
                ys[i] = s[first + i];
            }
            real_t ds = 0.0;
            real_t ss = polint(xh, ys, LAMB_QROMB_K, 0.0, &ds);
            if(fabs(ds) <= eps * (1.0 + fabs(ss))){
                return ss;
            }
        }

        if(j + 1 < LAMB_QROMB_JMAX){
            s[j + 1] = s[j];
            h[j + 1] = 0.25 * h[j];
        }
        if(j > 0){
            it *= 2;
        }
    }

    GRTRaiseError("Too many steps in Lamb Romberg integration.\n");
}
