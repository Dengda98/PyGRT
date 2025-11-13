/**
 * @file   lamb1.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    使用广义闭合解求解第一类 Lamb 问题，参考：
 * 
 *        张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社
 */

#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "grt/lamb/elliptic.h"
#include "grt/lamb/lamb_util.h"
#include "grt/common/checkerror.h"


typedef struct {
    // 常量，不随时间变化
    MYREAL nu;
    MYREAL k;      ///< k = beta/alpha;
    MYREAL kk;     ///< kk = k*k;
    MYREAL kp;     ///< kp = sqrt(1 - k^2);
    MYREAL kpkp;   ///< kpkp = kp*kp;
    MYREAL h1;     ///< h1 = 2 * kk - 1;
    MYREAL h2;     ///< h2 = kk * kpkp;
    MYREAL h3;     ///< h3 = kpkp * h1^2;
    MYREAL b3;     ///< b3 = 3 * kp^2 - 1;
    MYREAL b6;     ///< b6 = 6 * kp^2 - 1;
    MYREAL c;      ///< c = 6 * h2 - 1;
    MYREAL sf;     ///< sin(phi)
    MYREAL cf;     ///< cos(phi)

    // 一元三次方程的三个根，其中 y30 为正根
    MYCOMPLEX ys[3];
    MYCOMPLEX ysp[3];

    // Rayleigh 波项
    MYREAL kpa;
    MYREAL kpakpa;
    MYREAL RaylR;
    MYREAL RaylQ[3][3];
} VARS;



static void ckim_P(MYREAL tbar, VARS *V, MYCOMPLEX ckim[3][6][10])
{
    MYREAL kk = V->kk;
    MYREAL h1 = V->h1;
    MYREAL h2 = V->h2;
    MYREAL h3 = V->h3;
    MYREAL b3 = V->b3;
    MYREAL b6 = V->b6;
    MYREAL c = V->c;

    MYREAL tbtb = tbar*tbar;
    MYREAL A = 5.0 + 8.0*tbtb - 12.0*kk*tbtb - 24.0*h2;
    MYREAL B = h1 * (1.0 + 5.0*tbtb - 6.0*kk*tbtb - 8.0*h2);
    MYREAL d = 2.0*h1 - tbtb;

    MYREAL aa = tbtb - kk;

    // 表 6.6.1
    // (ξ, i, m), ξ=1,2,  i=1...5,  m=0...9

    ckim[1][1][6] = -8.0*tbtb; 
    ckim[1][1][4] = -8.0*tbtb*h1;
    ckim[1][1][2] =  8.0*tbtb*h2;

    ckim[1][2][8] = -8.0;
    ckim[1][2][6] = -8.0*(h1 - aa);
    ckim[1][2][4] =  8.0*(h1*aa + h2);
    ckim[1][2][2] = -8.0*h2*aa;

    ckim[1][3][7] =  8.0*I*tbar;
    ckim[1][3][5] =  8.0*I*tbar*h1;
    ckim[1][3][3] = -8.0*I*tbar*h2;

    ckim[1][4][7] =  8.0*I*tbar;
    ckim[1][4][5] =  12.0*I*tbar*h1;
    ckim[1][4][3] =  6.0*I*tbar*h1*h1;
    ckim[1][4][1] =  I*tbar*h1*h1*h1;
    
    ckim[1][5][8] =  8.0;
    ckim[1][5][6] =  12.0*h1;
    ckim[1][5][4] =  6.0*h1*h1;
    ckim[1][5][2] =  h1*h1*h1;


    ckim[2][1][7] = -8.0*I*tbtb; 
    ckim[2][1][5] =  8.0*I*tbtb*b3;
    ckim[2][1][3] =  2.0*I*tbtb*h1*b6;
    ckim[2][1][1] =  2.0*I*tbtb*h3;

    ckim[2][2][9] = -8.0*I; 
    ckim[2][2][7] = -8.0*I*d; 
    ckim[2][2][5] = -2.0*I*A;
    ckim[2][2][3] = -2.0*I*B;
    ckim[2][2][1] = -2.0*I*aa*h3;

    ckim[2][3][8] = -8.0*tbar; 
    ckim[2][3][6] =  8.0*tbar*b3; 
    ckim[2][3][4] =  2.0*tbar*h1*b6; 
    ckim[2][3][2] =  2.0*tbar*h3; 
    
    ckim[2][4][8] = -8.0*tbar; 
    ckim[2][4][6] = -12.0*tbar*h1; 
    ckim[2][4][4] =  4.0*tbar*c; 
    ckim[2][4][2] =  4.0*tbar*h1*h2; 

    ckim[2][5][9] =  8.0*I; 
    ckim[2][5][7] =  12.0*I*h1; 
    ckim[2][5][5] = -4.0*I*c;
    ckim[2][5][3] = -4.0*I*h1*h2;
}

static void ckim_S1(MYREAL tbar, VARS *V, MYCOMPLEX ckim[3][6][10])
{
    MYREAL kk = V->kk;
    MYREAL kpkp = V->kpkp;

    MYREAL tbtb = tbar*tbar;
    MYREAL bb = tbtb - 1.0;
    MYREAL A = 15.0 - 10.0*tbtb + 8.0*kk*tbtb - 16.0*kk;
    MYREAL B = 9.0 - 10.0*tbtb + 16.0*kk*bb;
    MYREAL C0 = 4.0 - 3.0*tbtb + 2.0*kk*tbtb - 3.0*kk;


    MYREAL h[6][10] = {0};
    for(int i=0; i<6; ++i){
        for(int j=0; j<10; ++j){
            h[i][j] = i + j*kpkp;
        }
    }

    // 表 6.6.2
    // (ξ, i, m), ξ=1,2,  i=0...5,  m=0...9

    ckim[1][0][8] =  8.0; 
    ckim[1][0][6] =  12.0; 
    ckim[1][0][4] =  6.0; 
    ckim[1][0][2] =  1.0; 

    ckim[1][1][8] = -8.0; 
    ckim[1][1][6] =  4.0*(2.0*bb - h[1][4]);
    ckim[1][1][4] = -2.0*A;
    ckim[1][1][2] = -B;
    ckim[1][1][0] = -bb;

    ckim[1][2][6] = -8.0*tbtb;
    ckim[1][2][4] = -4.0*tbtb*h[1][4];
    ckim[1][2][2] = -2.0*tbtb*(8.0*kpkp - 3.0);
    ckim[1][2][0] =  tbtb;

    ckim[1][3][7] = -8.0*I*tbar;
    ckim[1][3][5] = -12.0*I*tbar;
    ckim[1][3][3] = -6.0*I*tbar;
    ckim[1][3][1] = -I*tbar;

    ckim[1][4][7] = -8.0*I*tbar;
    ckim[1][4][5] = -8.0*I*tbar*h[1][1];
    ckim[1][4][3] = -8.0*I*tbar*kpkp;

    ckim[1][5][8] = -8.0;
    ckim[1][5][6] = -8.0*h[2][1];
    ckim[1][5][4] = -8.0*h[1][2];
    ckim[1][5][2] = -8.0*kpkp;


    ckim[2][0][9] =  8.0;
    ckim[2][0][7] =  4.0*h[3][2];
    ckim[2][0][5] =  4.0*h[1][3];
    ckim[2][0][3] =  4.0*kpkp;

    ckim[2][1][9] = -8.0;
    ckim[2][1][7] =  4.0*(2.0*tbtb - h[3][2]);
    ckim[2][1][5] = -4.0*C0;
    ckim[2][1][3] =  4.0*bb*kpkp;

    ckim[2][2][7] = -8.0*tbtb;
    ckim[2][2][5] = -4.0*tbtb*h[1][2];
    ckim[2][2][3] = -4.0*tbtb*kpkp;

    ckim[2][3][8] = -8.0*I*tbar;
    ckim[2][3][6] = -4.0*I*tbar*h[3][2];
    ckim[2][3][4] = -4.0*I*tbar*h[1][3];
    ckim[2][3][2] = -4.0*I*tbar*kpkp;

    ckim[2][4][8] = -8.0*I*tbar;
    ckim[2][4][6] = -8.0*I*tbar*h[1][1];
    ckim[2][4][4] = -2.0*I*tbar*h[1][4];
    ckim[2][4][2] = -2.0*I*tbar*kpkp;

    ckim[2][5][9] = -8.0;
    ckim[2][5][7] = -8.0*h[2][1];
    ckim[2][5][5] = -2.0*h[5][8];
    ckim[2][5][3] = -2.0*h[1][5];
    ckim[2][5][1] = -2.0*kpkp;
}

static void ckim_S2(MYREAL tbar, VARS *V, MYCOMPLEX ckim[3][6][10])
{
    ckim_S1(tbar, V, ckim);

    for(int k=1; k<=2; ++k){
        for(int i=0; i<6; ++i){
            ckim[k][i][1] *= -I;
            ckim[k][i][5] *= -I;
            ckim[k][i][9] *= -I;

            ckim[k][i][3] *= I;
            ckim[k][i][7] *= I;

            ckim[k][i][2] *= -1;
            ckim[k][i][6] *= -1;
        }
    }
}

static void Cmat(VARS *V, MYCOMPLEX ckim[3][6][10], MYCOMPLEX C[3][3][3][10])
{
    MYREAL sf = V->sf;
    MYREAL cf = V->cf;
    // 构建矩阵 C (3x3) 
    for(int k=1; k<=2; ++k){
        for(int m=0; m<=9; ++m){
            C[0][0][k][m] = ckim[k][0][m] + ckim[k][1][m]*cf*cf + ckim[k][2][m]*sf*sf;
            C[0][1][k][m] = (ckim[k][1][m] - ckim[k][2][m])*sf*cf;
            C[0][2][k][m] = ckim[k][3][m]*cf;

            C[1][0][k][m] = C[0][1][k][m];
            C[1][1][k][m] = ckim[k][0][m] + ckim[k][1][m]*sf*sf + ckim[k][2][m]*cf*cf;
            C[1][2][k][m] = ckim[k][3][m]*sf;

            C[2][0][k][m] = ckim[k][4][m]*cf;
            C[2][1][k][m] = ckim[k][4][m]*sf;
            C[2][2][k][m] = ckim[k][5][m];
        }
    }
}

static void uv(int sgn, MYREAL kpkp, MYCOMPLEX ys[3], MYCOMPLEX C[3][3][3][10], MYCOMPLEX u[3][3][10], MYCOMPLEX v[3][3][11])
{
    MYCOMPLEX y1, y2, y3;
    y1 = ys[0]*sgn;
    y2 = ys[1]*sgn;
    y3 = ys[2]*sgn;

    // 构建系数 u_ij,m 和 v_ij,n
    // i,j = 1,2,3   m=1...9   n=1...10
    MYCOMPLEX delta[4] = {0};
    delta[1] = -16.0 * sgn * kpkp *(y1 - y2)*(y1 - y3);
    delta[2] = -16.0 * sgn * kpkp *(y2 - y1)*(y2 - y3);
    delta[3] = -16.0 * sgn * kpkp *(y3 - y1)*(y3 - y2);

    for(int i1=0; i1<3; ++i1){
        for(int i2=0; i2<3; ++i2){
            u[i1][i2][1] = grt_evalpoly2(C[i1][i2][1], 3, y1, 1) / delta[1];
            u[i1][i2][2] = grt_evalpoly2(C[i1][i2][1], 4, y1, 0) / delta[1];
            u[i1][i2][3] = grt_evalpoly2(C[i1][i2][1], 3, y2, 1) / delta[2];
            u[i1][i2][4] = grt_evalpoly2(C[i1][i2][1], 4, y2, 0) / delta[2];
            u[i1][i2][5] = grt_evalpoly2(C[i1][i2][1], 3, y3, 1) / delta[3];
            u[i1][i2][6] = grt_evalpoly2(C[i1][i2][1], 4, y3, 0) / delta[3];
            u[i1][i2][7] = -1.0/(16.0*kpkp) * sgn * (C[i1][i2][1][6] + (y1+y2+y3)*C[i1][i2][1][8]);
            u[i1][i2][8] = -1.0/(16.0*kpkp) * sgn * C[i1][i2][1][7];
            u[i1][i2][9] = -1.0/(16.0*kpkp) * sgn * C[i1][i2][1][8];

            v[i1][i2][1] = grt_evalpoly2(C[i1][i2][2], 4, y1, 1) / delta[1];
            v[i1][i2][2] = grt_evalpoly2(C[i1][i2][2], 4, y1, 0) / delta[1];
            v[i1][i2][3] = grt_evalpoly2(C[i1][i2][2], 4, y2, 1) / delta[2];
            v[i1][i2][4] = grt_evalpoly2(C[i1][i2][2], 4, y2, 0) / delta[2];
            v[i1][i2][5] = grt_evalpoly2(C[i1][i2][2], 4, y3, 1) / delta[3];
            v[i1][i2][6] = grt_evalpoly2(C[i1][i2][2], 4, y3, 0) / delta[3];
            v[i1][i2][7] = -1.0/(16.0*kpkp) * sgn * (C[i1][i2][2][6] + (y1+y2+y3)*C[i1][i2][2][8]);
            v[i1][i2][8] = -1.0/(16.0*kpkp) * sgn * (C[i1][i2][2][7] + (y1+y2+y3)*C[i1][i2][2][9]);
            // !!!!!!!!!!!!!
            v[i1][i2][9]  = -1.0/(16.0*kpkp) * sgn * C[i1][i2][2][8];
            v[i1][i2][10] = -1.0/(16.0*kpkp) * sgn * C[i1][i2][2][9];
        }
    }
}

static MYCOMPLEX U_P(int n, MYREAL tbar, MYCOMPLEX c, VARS *V)
{
    MYREAL kk = V->kk;
    MYREAL tbtb = tbar*tbar;
    MYREAL aa = tbtb - kk;
    MYREAL a = sqrt(aa);  // 要求 tbar >= k

    MYCOMPLEX xi = sqrt(aa - c + 0*I);

    if(n==1){
        if(cimag(c) == 0.0 && creal(c) < aa){
            return 0.5/xi * log(fabs((creal(xi)+a)/(creal(xi)-a)));
        } else {
            return I/xi *  atan(a/(I*xi));
        }
    }
    else if(n == 2){
        if(cimag(c) == 0.0 && creal(c) < aa && creal(c) > 0.0){
            return 0.0;
        } else {
            return -0.5 * M_PI/c * sqrt(c / (c - aa));
        }
    }
    else if(n == 3){
        return M_PI * 0.5;
    }
    else if(n == 4){
        return a;
    }
    else if(n == 5){
        return M_PI*aa*0.25;
    }
    GRTRaiseError("Wrong execution in function %s.", __func__);
}

static MYCOMPLEX V_P(int n, MYREAL tbar, MYCOMPLEX c, VARS *V)
{
    MYREAL kp = V->kp;
    MYREAL kk = V->kk;
    MYREAL kpkp = V->kpkp;
    MYREAL tbtb = tbar*tbar;
    MYREAL aa = tbtb - kk;
    MYREAL a = sqrt(aa);  // 要求 tbar >= k
    MYREAL bb = tbtb - 1.0;

    MYCOMPLEX xi = sqrt(aa - c + 0*I);
    MYCOMPLEX et = sqrt(kpkp - c + 0*I);
    
    MYREAL m1 = aa/kpkp;
    MYREAL m1p = kpkp/aa;
    MYREAL m2 = bb/aa;
    MYCOMPLEX n1 = bb/(xi*xi);

    if(tbar > 1.0){
        if(n == 1){
            if(cimag(c) == 0.0 && ((creal(c)) > aa || creal(c) < kpkp)){
                return 0.5/(xi*et) * (log(fabs((c + a*kp + xi*et)*(xi + et)/((c + a*kp - xi*et)*(xi - et)))) - I*M_PI);
            } else {
                return -I/(xi*et) * (atan(I*kp*xi/(a*et)) + M_PI*0.5);
            }
        }
        else if(n == 2){
            return - 1.0/(a*c) * ( grt_ellipticPi(kpkp/c, m1p) + I*c/(xi*xi)*grt_ellipticPi(n1, m2) );
        }
        else if(n == 3){
            return 1.0/a*( grt_ellipticK(m1p) - I*grt_ellipticK(m2) );
        }
        else if(n == 4){
            return 0.5*(log((a + kp)/(a - kp)) - I*M_PI);
        }
        else if(n == 5){
            return -a * ( grt_ellipticE(m1p) - grt_ellipticK(m1p) + I*grt_ellipticE(m2) );
        }
        else if(n == 6){
            return 0.25*(kpkp + aa)*log((a + kp)/(a - kp)) - 0.5*kp*a - 0.25*I*M_PI*(kpkp + aa);
        }
    }
    else {
        if(n == 1){
            if(cimag(c) == 0.0 && ((creal(c)) < aa || creal(c) > kpkp)){
                return 0.5/(xi*et) * log(fabs((c + a*kp + xi*et)*(xi + et)/((c + a*kp - xi*et)*(xi - et))));
            } else {
                return -I/(xi*et) * atan(I*a*et/(kp*xi));
            }
        }
        else if(n == 2){
            return - 1.0/(kp*c) * grt_ellipticPi(aa/c, m1);
        }
        else if(n == 3){
            return 1.0/kp*grt_ellipticK(m1);
        }
        else if(n == 4){
            return 0.5*log((a + kp)/(kp - a));
        }
        else if(n == 5){
            return -kp * ( grt_ellipticE(m1) - grt_ellipticK(m1) );
        }
        else if(n == 6){
            return 0.25*(kpkp + aa)*log((a + kp)/(kp - a)) - 0.5*kp*a;
        }
    }

    GRTRaiseError("Wrong execution in function %s.", __func__);
}

static MYCOMPLEX U_S1(int n, MYREAL tbar, MYCOMPLEX c, VARS *V)
{
    MYREAL tbtb = tbar*tbar;
    MYREAL bb = tbtb - 1.0;
    MYREAL b = sqrt(bb);  // 要求 tbar >= 1.0

    MYCOMPLEX xip = sqrt(bb - c + 0*I);

    if(n==1){
        if(cimag(c) == 0.0 && creal(c) < bb){
            return 0.5/xip * log(fabs((creal(xip)+b)/(creal(xip)-b)));
        } else {
            return I/xip *  atan(b/(I*xip));
        }
    }
    else if(n == 2){
        if(cimag(c) == 0.0 && creal(c) < bb && creal(c) > 0.0){
            return 0.0;
        } else {
            return -0.5 * M_PI/c * sqrt(c / (c - bb));
        }
    }
    else if(n == 3){
        return M_PI * 0.5;
    }
    else if(n == 4){
        return b;
    }
    else if(n == 5){
        return M_PI*bb*0.25;
    }
    GRTRaiseError("Wrong execution in function %s.", __func__);
}

static MYCOMPLEX V_S1(int n, MYREAL tbar, MYCOMPLEX c, VARS *V)
{
    MYREAL kp = V->kp;
    MYREAL kk = V->kk;
    MYREAL kpkp = V->kpkp;
    MYREAL tbtb = tbar*tbar;
    MYREAL aa = tbtb - kk;
    MYREAL a = sqrt(aa);  // 要求 tbar >= k
    MYREAL bb = tbtb - 1.0;
    MYREAL b = sqrt(bb);  // 要求 tbar >= 1.0

    MYCOMPLEX xip = sqrt(bb - c + 0*I);
    MYCOMPLEX etp = sqrt(kpkp + c + 0*I);
    
    MYREAL m2 = bb/aa;
    MYCOMPLEX n2 = bb/(xip*xip);

    if(n==1){
        if(cimag(c) == 0.0 && creal(c) > - kpkp && creal(c) < bb){
            return 0.5/(xip*etp) * log(fabs((kp*xip + b*etp)/(kp*xip - b*etp)));
        } else {
            return I/(xip*etp) *  atan(b*etp/(I*kp*xip));
        }
    }
    else if(n == 2){
        return 1.0/(a*xip*xip) * grt_ellipticPi(n2, m2);
    }
    else if(n == 3){
        return 1.0/a * grt_ellipticK(m2);
    }
    else if(n == 4){
        return atan(b/kp);
    }
    else if(n == 5){
        return a*grt_ellipticE(m2) - kpkp/a*grt_ellipticK(m2);
    }
    else if(n == 6){
        return 0.5*(b*kp + (bb - kpkp)*atan(b/kp));
    }
    GRTRaiseError("Wrong execution in function %s.", __func__);
}

static MYCOMPLEX U_S2(int n, MYREAL tbar, MYCOMPLEX c, VARS *V)
{
    MYREAL kp = V->kp;
    MYREAL kk = V->kk;
    MYREAL tbtb = tbar*tbar;
    MYREAL aa = tbtb - kk;
    MYREAL a = sqrt(aa);  // 要求 tbar >= k
    MYREAL bb = tbtb - 1.0;
    MYREAL b = sqrt(bb);  // 要求 tbar >= 1.0
    MYREAL f = kp + a;
    MYCOMPLEX d = sqrt(c*(bb - c));

    MYCOMPLEX xip = sqrt(bb - c + 0*I);

    if(n == 1){
        return 1.0/(I*xip) * atan((a - b)*I*xip/(c + b*(a-b)));
    }
    else if(n == 2){
        return 1.0/d * atan(kp*xip*xip/(a*d));
    }
    else if(n == 3){
        return log(f/b);
    }
    else if(n == 4){
        return 0.5*(f-b)*(f-b)/f;
    }
    else if(n == 5){
        MYREAL f4 = f*f*f*f;
        MYREAL b4 = bb*bb;
        return 0.125*(f4 - b4)/(f*f) - 0.5*bb*log(f/b);
    }
    GRTRaiseError("Wrong execution in function %s.", __func__);
}

static MYCOMPLEX V_S2(int n, MYREAL tbar, MYCOMPLEX c, VARS *V)
{
    MYREAL kp = V->kp;
    MYREAL kk = V->kk;
    MYREAL kpkp = V->kpkp;
    MYREAL tbtb = tbar*tbar;
    MYREAL aa = tbtb - kk;
    MYREAL a = sqrt(aa);  // 要求 tbar >= k
    MYREAL bb = tbtb - 1.0;
    MYREAL b = sqrt(bb);  // 要求 tbar >= 1.0

    MYCOMPLEX xip = sqrt(bb - c + 0*I);
    MYCOMPLEX etp = sqrt(kpkp + c + 0*I);
    
    MYREAL m1p = kpkp/aa;
    MYCOMPLEX n3 = kpkp/(etp*etp);

    if(n==1){
        if(cimag(c) == 0.0 && creal(c) > - kpkp && creal(c) < bb){
            return 0.5/(xip*etp) * log(fabs((kp*xip + b*etp)/(kp*xip - b*etp)));
        } else {
            return - I/(xip*etp) *  atan(I*kp*xip/(b*etp));
        }
    }
    else if(n == 2){
        return 1.0/(a*etp*etp) * grt_ellipticPi(n3, m1p);
    }
    else if(n == 3){
        return 1.0/a * grt_ellipticK(m1p);
    }
    else if(n == 4){
        return atan(kp/b);
    }
    else if(n == 5){
        return a*grt_ellipticE(m1p) - bb/a*grt_ellipticK(m1p);
    }
    else if(n == 6){
        return 0.5*(b*kp + (kpkp - bb)*atan(kp/b));
    }
    GRTRaiseError("Wrong execution in function %s.", __func__);
}

static MYCOMPLEX U_SP(int n, MYREAL tbar, MYCOMPLEX c, VARS *V)
{
    MYREAL kk = V->kk;
    MYREAL kp = V->kp;
    MYREAL tbtb = tbar*tbar;
    MYREAL aa = tbtb - kk;
    MYREAL a = sqrt(aa);  // 要求 tbar >= k
    MYREAL bb = tbtb - 1.0;
    MYREAL bp = sqrt(- bb);  // 要求 tbar <= 1.0
    MYREAL f = kp + a;
    MYCOMPLEX d = sqrt(c*(bb - c));

    MYCOMPLEX xip = sqrt(bb - c + 0*I);

    if(n == 1){
        return - 1.0/(I*xip) * atan(a/(I*xip));
    }
    else if(n == 2){
        return - 1.0/d * atan((kp*f+bb)*c/(kp*f*d));
    }
    else if(n == 3){
        return log(f/bp);
    }
    else if(n == 4){
        return 0.5*(f*f + bb)/f;
    }
    else if(n == 5){
        MYREAL f4 = f*f*f*f;
        MYREAL b4 = bb*bb;
        return 0.125*(f4 - b4)/(f*f) - 0.5*bb*log(f/bp);
    }
    GRTRaiseError("Wrong execution in function %s.", __func__);
}

static MYCOMPLEX V_SP(int n, MYREAL tbar, MYCOMPLEX c, VARS *V)
{
    MYREAL kk = V->kk;
    MYREAL kp = V->kp;
    MYREAL kpkp = V->kpkp;
    MYREAL tbtb = tbar*tbar;
    MYREAL aa = tbtb - kk;
    MYREAL bb = tbtb - 1.0;

    MYCOMPLEX etp = sqrt(kpkp + c + 0*I);
    
    MYREAL m1 = aa/kpkp;
    MYCOMPLEX n4 = aa/(etp*etp);

    if(n==1){
        if(cimag(c) == 0.0 && creal(c) > - kpkp && creal(c) < bb){
            return 0.0;
        } else {
            // 这里未使用 xip 替代，以方便处理相位
            return M_PI_2 / (c - bb) * sqrt((c - bb)/(c + kpkp));  
        }
    }
    else if(n == 2){
        return 1.0/(kp*etp*etp) * grt_ellipticPi(n4, m1);
    }
    else if(n == 3){
        return 1.0/kp * grt_ellipticK(m1);
    }
    else if(n == 4){
        return M_PI_2;
    }
    else if(n == 5){
        return kp*grt_ellipticE(m1);
    }
    else if(n == 6){
        return M_PI_4 * (kpkp - bb);
    }
    GRTRaiseError("Wrong execution in function %s.", __func__);
}

static void build_raw(MYREAL tbar, VARS *V, MYCOMPLEX u[3][3], 
    MYCOMPLEX C[3][3][3][10], 
    MYCOMPLEX (*Ufunc)(int, MYREAL, MYCOMPLEX, VARS *), 
    MYCOMPLEX (*Vfunc)(int, MYREAL, MYCOMPLEX, VARS *), 
    int sgn, MYCOMPLEX ys[3])
{
    // 计算 u, v 系数
    MYCOMPLEX uc[3][3][10] = {0};
    MYCOMPLEX vc[3][3][11] = {0};
    uv(sgn, V->kpkp, ys, C, uc, vc);

    MYCOMPLEX y1, y2, y3;
    y1 = ys[0];
    y2 = ys[1];
    y3 = ys[2];

    // 计算 UP, VP 函数
    MYCOMPLEX UF[6][4] = {0};
    UF[1][1] = Ufunc(1, tbar, y1, V);
    UF[1][2] = Ufunc(1, tbar, y2, V);
    UF[1][3] = Ufunc(1, tbar, y3, V);
    UF[2][1] = Ufunc(2, tbar, y1, V);
    UF[2][2] = Ufunc(2, tbar, y2, V);
    UF[2][3] = Ufunc(2, tbar, y3, V);
    for(int i=3; i<=5; ++i){
        UF[i][1] = Ufunc(i, tbar, y3, V);
    }

    MYCOMPLEX VF[7][4] = {0};
    VF[1][1] = Vfunc(1, tbar, y1, V);
    VF[1][2] = Vfunc(1, tbar, y2, V);
    VF[1][3] = Vfunc(1, tbar, y3, V);
    VF[2][1] = Vfunc(2, tbar, y1, V);
    VF[2][2] = Vfunc(2, tbar, y2, V);
    VF[2][3] = Vfunc(2, tbar, y3, V);
    for(int i=3; i<=6; ++i){
        VF[i][1] = Vfunc(i, tbar, y3, V);
    }

    // 组合成 P 波项
    for(int i1=0; i1<3; ++i1){
        for(int i2=0; i2<3; ++i2){
            MYCOMPLEX tmp = 0.0;
            tmp +=  uc[i1][i2][1]*UF[1][1] + uc[i1][i2][2]*UF[2][1]
                  + uc[i1][i2][3]*UF[1][2] + uc[i1][i2][4]*UF[2][2]
                  + uc[i1][i2][5]*UF[1][3] + uc[i1][i2][6]*UF[2][3];
            for(int q=7; q <= 9; ++q){
                tmp += uc[i1][i2][q]*UF[q-4][1];
            }
            tmp +=  vc[i1][i2][1]*VF[1][1] + vc[i1][i2][2]*VF[2][1]
                  + vc[i1][i2][3]*VF[1][2] + vc[i1][i2][4]*VF[2][2]
                  + vc[i1][i2][5]*VF[1][3] + vc[i1][i2][6]*VF[2][3];
            for(int q=7; q <= 10; ++q){
                tmp += vc[i1][i2][q]*VF[q-4][1];
            }
            u[i1][i2] = tmp;
        }
    }
}

static void build_P(MYREAL tbar, VARS *V, MYREAL u[3][3])
{
    MYREAL k = V->k;
    if(tbar < k) return;

    MYCOMPLEX cu[3][3];
    MYCOMPLEX ckim[3][6][10] = {0};
    MYCOMPLEX C[3][3][3][10] = {0};
    ckim_P(tbar, V, ckim);
    Cmat(V, ckim, C);
    build_raw(tbar, V, cu, C, U_P, V_P, 1, V->ys);
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            u[i][j] = creal(cu[i][j]);
        }
    }
}

static void build_S1(MYREAL tbar, VARS *V, MYREAL u[3][3])
{
    if(tbar < 1.0) return;

    MYCOMPLEX cu[3][3];
    MYCOMPLEX ckim[3][6][10] = {0};
    MYCOMPLEX C[3][3][3][10] = {0};
    ckim_S1(tbar, V, ckim);
    Cmat(V, ckim, C);
    build_raw(tbar, V, cu, C, U_S1, V_S1, 1, V->ysp);
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            u[i][j] = creal(cu[i][j]);
        }
    }
}

static void build_S2_SP(MYREAL tbar, VARS *V, MYREAL us2[3][3], MYREAL usp[3][3])
{
    MYREAL k = V->k;
    if(tbar < k) return;
    
    MYCOMPLEX cu[3][3];
    MYCOMPLEX ckim[3][6][10] = {0};
    MYCOMPLEX C[3][3][3][10] = {0};
    ckim_S2(tbar, V, ckim);
    Cmat(V, ckim, C);

    if(tbar < 1.0){
        build_raw(tbar, V, cu, C, U_SP, V_SP, -1, V->ysp);
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                usp[i][j] = - cimag(cu[i][j]);
            }
        }
    }
    else {
        build_raw(tbar, V, cu, C, U_S2, V_S2, -1, V->ysp);
        for(int i=0; i<3; ++i){
            for(int j=0; j<3; ++j){
                us2[i][j] = - cimag(cu[i][j]);
            }
        }
    }
   
}

static void build_R(MYREAL tbar, VARS *V, MYREAL u[3][3])
{
    MYREAL kpakpa = V->kpakpa;
    MYREAL kpa = V->kpa;
    if(tbar < kpa)  return;

    MYREAL coef = M_PI_4 * tbar * V->RaylR / sqrt(tbar*tbar - kpakpa);

    for(int i1=0; i1<3; ++i1){
        for(int i2=0; i2<3; ++i2){
            u[i1][i2] = V->RaylQ[i1][i2] * coef;
        }
    }
}





void grt_solve_lamb1(
    const MYREAL nu, const MYREAL *ts, const int nt, const MYREAL azimuth, MYREAL (*u)[3][3])
{
    // 检查泊松比范围
    if(nu <= 0.0 || nu >= 0.5){
        GRTRaiseError("possion ratio (%lf) is out of bound.", nu);
    }

    // 根据情况判断是打印在屏幕还是记录到内存中
    bool isprint = (u == NULL);

    // 先打印标题
    if(isprint){
        char *stmp = NULL;
        printf("#");
        printf("%13s", "tbar");
        for(int i1=0; i1<3; ++i1){
        for(int i2=0; i2<3; ++i2){
            GRT_SAFE_ASPRINTF(&stmp, "G%d%d", i1+1, i2+1);
            printf("%14s", stmp);
        }}
        GRT_SAFE_FREE_PTR(stmp);
        printf("\n");
    }

    MYREAL phi = azimuth * DEG1;
    
    MYREAL tbar_eps = GRT_MIN(1e-8, (ts[1]-ts[0]) * 1e-5);

    // 初始化相关变量
    VARS V0 = {0};
    VARS *V = &V0;
    V0.kk = 0.5 * (1.0 - 2.0*nu)/(1.0 - nu);
    V0.k = sqrt(V0.kk);
    V0.nu = nu;
    V0.kpkp = 1.0 - V0.kk;
    V0.kp = sqrt(V0.kpkp);
    V0.h1 = 2.0 * V0.kk - 1.0;
    V0.h2 = V0.kk * V0.kpkp;
    V0.h3 = V0.kpkp * V0.h1 * V0.h1;
    V0.b3 = 3.0 * V0.kpkp - 1.0;
    V0.b6 = 6.0 * V0.kpkp - 1.0;
    V0.c = 6.0 * V0.h2 - 1.0;
    V0.sf = sin(phi);
    V0.cf = cos(phi);

    // 求一元三次方程的根
    {
        MYREAL a, b, c;
        MYREAL nu2, nu3, nu4;
        nu2 = V0.nu*V0.nu;
        nu3 = nu2*V0.nu;
        nu4 = nu3*V0.nu;
        MYREAL snu = 1.0 - V0.nu;
        MYREAL snu2 = snu*snu;
        MYREAL snu3 = snu2*snu;
        a = -0.5 * (2.0*nu2 + 1.0)/snu;
        b = 0.25 * (4.0*nu3 - 4.0*nu2 + 4.0*V0.nu - 1.0)/snu2;
        c = -0.125*nu4/snu3;
        grt_roots3(a, b, c, V0.ys);
        for(int i=0; i<3; ++i){
            V0.ysp[i] = V0.ys[i] - V0.kpkp;
        }
    }

    // 另一种形式的Rayleigh波函数
    {   
        MYCOMPLEX y3[3];
        MYREAL a, b, c;
        MYREAL m = V0.kk;
        a = 0.5*(2.0*m - 3.0)/(1 - m);
        b = 0.5/(1 - m);
        c = - 0.0625/(1 - m);
        grt_roots3(a, b, c, y3);
        
        V0.RaylQ[0][2] = V->cf;
        V0.RaylQ[1][2] = V->sf;
        V0.RaylQ[2][0] = - V0.RaylQ[0][2];
        V0.RaylQ[2][1] = - V0.RaylQ[1][2];

        V0.kpakpa = creal(y3[2]);
        V0.kpa = sqrt(V0.kpakpa);

        MYREAL u0 = sqrt(V0.kpakpa - V0.kk);
        MYREAL v0 = sqrt(V0.kpakpa - 1.0);

        MYREAL R1, R2;
        R1 = (1.0 - 2.0*V0.kpakpa)*u0*v0 + 2.0*u0*u0*v0*v0;
        R2 = 2.0*(1.0 - 2.0*V0.kpakpa)*u0*v0 + 2.0*u0*u0*v0*v0 + V0.kpakpa*(u0*u0 + v0*v0);
        V0.RaylR = R1 / R2;
    }

    for(int i=0; i < nt; ++i){
        MYREAL up[3][3] = {0};
        MYREAL us1[3][3] = {0};
        MYREAL us2[3][3] = {0};
        MYREAL usp[3][3] = {0};
        MYREAL uR[3][3] = {0};
        MYREAL tbar = ts[i];

        // 跳过一些震相到时处的奇点
        if(tbar == 1.0 || tbar == V0.k || tbar == V0.kpa)   tbar += tbar_eps;

        build_P(tbar,  V, up);
        build_S1(tbar, V, us1);
        build_S2_SP(tbar, V, us2, usp);
        build_R(tbar, V, uR);

        if(isprint){
            printf("%14.6e", tbar);
            for(int i1=0; i1<3; ++i1){
            for(int i2=0; i2<3; ++i2){
                printf("%14.6e", up[i1][i2] + us1[i1][i2] + us2[i1][i2] + usp[i1][i2] + uR[i1][i2]);
            }}
            printf("\n");
        } else {
            for(int i1=0; i1<3; ++i1){
            for(int i2=0; i2<3; ++i2){
                u[i][i1][i2] = up[i1][i2] + us1[i1][i2] + us2[i1][i2] + usp[i1][i2] + uR[i1][i2];
            }}
        }
    }
}