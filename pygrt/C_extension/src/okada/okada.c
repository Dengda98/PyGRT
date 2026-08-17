/**
 * @file   okada.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * Okada 静态位移解析解（https://www.bosai.go.jp/e/dc3d.html）
 *
 * 本文件按 Okada 的 Fortran 源码改写，内部数组仍保持原程序的 12 个量的顺序
 *
 */

#include <math.h>
#include <string.h>

#include "grt/okada/okada.h"

/** Okada 介质常数及倾角三角函数 */
typedef struct {
    real_t alp1;
    real_t alp2;
    real_t alp3;
    real_t alp4;
    real_t alp5;
    real_t sd;
    real_t cd;
    real_t sdsd;
    real_t cdcd;
    real_t sdcd;
    real_t s2d;
    real_t c2d;
} OKADA_MEDIUM;

/** 点源公式中重复使用的几何量 */
typedef struct {
    real_t x;
    real_t y;
    real_t d;
    real_t p;
    real_t q;
    real_t s;
    real_t t;
    real_t xy;
    real_t x2;
    real_t y2;
    real_t d2;
    real_t r;
    real_t r2;
    real_t r3;
    real_t r5;
    real_t qr;
    real_t qrx;
    real_t a3;
    real_t a5;
    real_t b3;
    real_t c3;
    real_t uy;
    real_t vy;
    real_t wy;
    real_t uz;
    real_t vz;
    real_t wz;
} OKADA_POINT;

/** 有限断层角点公式中重复使用的几何量 */
typedef struct {
    real_t xi;
    real_t et;
    real_t q;
    real_t xi2;
    real_t et2;
    real_t q2;
    real_t r;
    real_t r2;
    real_t r3;
    real_t r5;
    real_t y;
    real_t d;
    real_t tt;
    real_t alx;
    real_t ale;
    real_t x11;
    real_t y11;
    real_t x32;
    real_t y32;
    real_t ey;
    real_t ez;
    real_t fy;
    real_t fz;
    real_t gy;
    real_t gz;
    real_t hy;
    real_t hz;
} OKADA_FINITE;

static const real_t OKADA_EPS = 1e-6;
static const real_t OKADA_PI2 = 6.283185307179586476925286766559;

/** 将 12 个位移及偏导结果清零
 *
 * @param[out]  u   待清零的结果数组
 */
static void okada_zero_u(real_t u[12])
{
    memset(u, 0, sizeof(real_t) * 12);
}

/** 将一个公式贡献按给定系数累加到结果数组
 *
 * @param[in,out] u       累加结果数组
 * @param[in]     du      待累加的公式贡献
 * @param[in]     scale   累加系数
 */
static void okada_add_scaled(real_t u[12], const real_t du[12], real_t scale)
{
    for(int i = 0; i < 12; ++i) u[i] += scale * du[i];
}

/** 将内部 12 元数组整理为公开 API 的位移和偏导矩阵
 *
 * @param[in]   u           内部结果数组
 * @param[out]  displacement 位移结果
 * @param[out]  upar        位移偏导结果
 */
static void okada_copy_output(const real_t u[12], real_t displacement[3], real_t upar[3][3])
{
    displacement[0] = u[0];
    displacement[1] = u[1];
    displacement[2] = u[2];

    upar[0][0] = u[3];
    upar[0][1] = u[4];
    upar[0][2] = u[5];
    upar[1][0] = u[6];
    upar[1][1] = u[7];
    upar[1][2] = u[8];
    upar[2][0] = u[9];
    upar[2][1] = u[10];
    upar[2][2] = u[11];
}

/** 初始化介质常数和断层倾角相关三角函数
 *
 * 对应 Fortran 源码中的 DCCON0
 *
 * @param[out]  c       Okada 介质常数
 * @param[in]   alpha   均匀介质常数
 * @param[in]   dip     断层倾角，单位为度
 */
static void DCCON0(OKADA_MEDIUM *c, real_t alpha, real_t dip)
{
    c->alp1 = (1.0 - alpha) / 2.0;
    c->alp2 = alpha / 2.0;
    c->alp3 = (1.0 - alpha) / alpha;
    c->alp4 = 1.0 - alpha;
    c->alp5 = alpha;

    c->sd = sin(dip * PI / 180.0);
    c->cd = cos(dip * PI / 180.0);
    if(fabs(c->cd) < OKADA_EPS){
        c->cd = 0.0;
        c->sd = (c->sd > 0.0) ? 1.0 : ((c->sd < 0.0) ? -1.0 : 0.0);
    }
    c->sdsd = c->sd * c->sd;
    c->cdcd = c->cd * c->cd;
    c->sdcd = c->sd * c->cd;
    c->s2d = 2.0 * c->sdcd;
    c->c2d = c->cdcd - c->sdsd;
}

/** 初始化点源公式使用的几何量
 *
 * 对应 Fortran 源码中的 DCCON1
 *
 * @param[in]   c       Okada 介质常数
 * @param[out]  p       点源几何量
 * @param[in]   x       局部 X 坐标
 * @param[in]   y       局部 Y 坐标
 * @param[in]   d       震源到观测点的深度方向距离
 */
static void DCCON1(const OKADA_MEDIUM *c, OKADA_POINT *p, real_t x, real_t y, real_t d)
{
    if(fabs(x) < OKADA_EPS) x = 0.0;
    if(fabs(y) < OKADA_EPS) y = 0.0;
    if(fabs(d) < OKADA_EPS) d = 0.0;

    p->x = x;
    p->y = y;
    p->d = d;
    p->p = y * c->cd + d * c->sd;
    p->q = y * c->sd - d * c->cd;
    p->s = p->p * c->sd + p->q * c->cd;
    p->t = p->p * c->cd - p->q * c->sd;
    p->xy = x * y;
    p->x2 = x * x;
    p->y2 = y * y;
    p->d2 = d * d;
    p->r2 = p->x2 + p->y2 + p->d2;
    p->r = sqrt(p->r2);
    if(p->r == 0.0) return;

    p->r3 = p->r * p->r2;
    p->r5 = p->r3 * p->r2;
    p->a3 = 1.0 - 3.0 * p->x2 / p->r2;
    p->a5 = 1.0 - 5.0 * p->x2 / p->r2;
    p->b3 = 1.0 - 3.0 * p->y2 / p->r2;
    p->c3 = 1.0 - 3.0 * p->d2 / p->r2;
    p->qr = 3.0 * p->q / p->r5;
    p->qrx = 5.0 * p->qr * x / p->r2;
    p->uy = c->sd - 5.0 * y * p->q / p->r2;
    p->uz = c->cd + 5.0 * d * p->q / p->r2;
    p->vy = p->s - 5.0 * y * p->p * p->q / p->r2;
    p->vz = p->t + 5.0 * d * p->p * p->q / p->r2;
    p->wy = p->uy + c->sd;
    p->wz = p->uz + c->cd;
}

/** 计算点源真实源项的公式贡献
 *
 * 对应 Fortran 源码中的 UA0，四类 potency 分别对应走向滑动、倾向滑动、张裂和爆炸源
 *
 * @param[in]   c       Okada 介质常数
 * @param[in]   p       点源几何量
 * @param[in]   pot1    走向滑动 potency
 * @param[in]   pot2    倾向滑动 potency
 * @param[in]   pot3    张裂 potency
 * @param[in]   pot4    爆炸源 potency
 * @param[out]  u       12 个位移及偏导结果
 */
static void UA0(const OKADA_MEDIUM *c, const OKADA_POINT *p,
    real_t pot1, real_t pot2, real_t pot3, real_t pot4, real_t u[12])
{
    const real_t alp1 = c->alp1, alp2 = c->alp2;
    const real_t sd = c->sd, cd = c->cd, s2d = c->s2d, c2d = c->c2d;
    const real_t x = p->x, y = p->y, d = p->d, p0 = p->p;
    const real_t q = p->q, s = p->s, t = p->t, xy = p->xy;
    const real_t x2 = p->x2, r3 = p->r3, r5 = p->r5;
    const real_t qr = p->qr, qrx = p->qrx, a3 = p->a3, a5 = p->a5;
    const real_t b3 = p->b3, c3 = p->c3, uy = p->uy, vy = p->vy;
    const real_t wy = p->wy, uz = p->uz, vz = p->vz, wz = p->wz;
    real_t du[12];
    okada_zero_u(u);

    // 走向滑动贡献
    if(pot1 != 0.0){
        du[0] = alp1 * q / r3 + alp2 * x2 * qr;
        du[1] = alp1 * x / r3 * sd + alp2 * xy * qr;
        du[2] = -alp1 * x / r3 * cd + alp2 * x * d * qr;
        du[3] = x * qr * (-alp1 + alp2 * (1.0 + a5));
        du[4] = alp1 * a3 / r3 * sd + alp2 * y * qr * a5;
        du[5] = -alp1 * a3 / r3 * cd + alp2 * d * qr * a5;
        du[6] = alp1 * (sd / r3 - y * qr) + alp2 * 3.0 * x2 / r5 * uy;
        du[7] = 3.0 * x / r5 * (-alp1 * y * sd + alp2 * (y * uy + q));
        du[8] = 3.0 * x / r5 * (alp1 * y * cd + alp2 * d * uy);
        du[9] = alp1 * (cd / r3 + d * qr) + alp2 * 3.0 * x2 / r5 * uz;
        du[10] = 3.0 * x / r5 * (alp1 * d * sd + alp2 * y * uz);
        du[11] = 3.0 * x / r5 * (-alp1 * d * cd + alp2 * (d * uz - q));
        okada_add_scaled(u, du, pot1 / OKADA_PI2);
    }

    // 倾向滑动贡献
    if(pot2 != 0.0){
        du[0] = alp2 * x * p0 * qr;
        du[1] = alp1 * s / r3 + alp2 * y * p0 * qr;
        du[2] = -alp1 * t / r3 + alp2 * d * p0 * qr;
        du[3] = alp2 * p0 * qr * a5;
        du[4] = -alp1 * 3.0 * x * s / r5 - alp2 * y * p0 * qrx;
        du[5] = alp1 * 3.0 * x * t / r5 - alp2 * d * p0 * qrx;
        du[6] = alp2 * 3.0 * x / r5 * vy;
        du[7] = alp1 * (s2d / r3 - 3.0 * y * s / r5)
            + alp2 * (3.0 * y / r5 * vy + p0 * qr);
        du[8] = -alp1 * (c2d / r3 - 3.0 * y * t / r5) + alp2 * 3.0 * d / r5 * vy;
        du[9] = alp2 * 3.0 * x / r5 * vz;
        du[10] = alp1 * (c2d / r3 + 3.0 * d * s / r5) + alp2 * 3.0 * y / r5 * vz;
        du[11] = alp1 * (s2d / r3 - 3.0 * d * t / r5)
            + alp2 * (3.0 * d / r5 * vz - p0 * qr);
        okada_add_scaled(u, du, pot2 / OKADA_PI2);
    }

    // 张裂贡献
    if(pot3 != 0.0){
        du[0] = alp1 * x / r3 - alp2 * x * q * qr;
        du[1] = alp1 * t / r3 - alp2 * y * q * qr;
        du[2] = alp1 * s / r3 - alp2 * d * q * qr;
        du[3] = alp1 * a3 / r3 - alp2 * q * qr * a5;
        du[4] = -alp1 * 3.0 * x * t / r5 + alp2 * y * q * qrx;
        du[5] = -alp1 * 3.0 * x * s / r5 + alp2 * d * q * qrx;
        du[6] = -alp1 * 3.0 * xy / r5 - alp2 * x * qr * wy;
        du[7] = alp1 * (c2d / r3 - 3.0 * y * t / r5) - alp2 * (y * wy + q) * qr;
        du[8] = alp1 * (s2d / r3 - 3.0 * y * s / r5) - alp2 * d * qr * wy;
        du[9] = alp1 * 3.0 * x * d / r5 - alp2 * x * qr * wz;
        du[10] = -alp1 * (s2d / r3 - 3.0 * d * t / r5) - alp2 * y * qr * wz;
        du[11] = alp1 * (c2d / r3 + 3.0 * d * s / r5) - alp2 * (d * wz - q) * qr;
        okada_add_scaled(u, du, pot3 / OKADA_PI2);
    }

    // 爆炸源贡献
    if(pot4 != 0.0){
        du[0] = -alp1 * x / r3;
        du[1] = -alp1 * y / r3;
        du[2] = -alp1 * d / r3;
        du[3] = -alp1 * a3 / r3;
        du[4] = alp1 * 3.0 * xy / r5;
        du[5] = alp1 * 3.0 * x * d / r5;
        du[6] = du[4];
        du[7] = -alp1 * b3 / r3;
        du[8] = alp1 * 3.0 * y * d / r5;
        du[9] = -du[5];
        du[10] = -du[8];
        du[11] = alp1 * c3 / r3;
        okada_add_scaled(u, du, pot4 / OKADA_PI2);
    }
}

/** 计算点源镜像源的 UB0 公式贡献
 *
 * @param[in]   c       Okada 介质常数
 * @param[in]   p       点源几何量
 * @param[in]   z       观测点 Z 坐标
 * @param[in]   pot1    走向滑动 potency
 * @param[in]   pot2    倾向滑动 potency
 * @param[in]   pot3    张裂 potency
 * @param[in]   pot4    爆炸源 potency
 * @param[out]  u       12 个位移及偏导结果
 */
static void UB0(const OKADA_MEDIUM *c, const OKADA_POINT *p, real_t z,
    real_t pot1, real_t pot2, real_t pot3, real_t pot4, real_t u[12])
{
    const real_t alp3 = c->alp3;
    const real_t sd = c->sd, sdcd = c->sdcd, sdsd = c->sdsd;
    const real_t x = p->x, y = p->y, d = p->d, p0 = p->p;
    const real_t q = p->q, qr = p->qr, qrx = p->qrx, xy = p->xy;
    const real_t x2 = p->x2, y2 = p->y2, d2 = p->d2;
    const real_t r = p->r, r2 = p->r2, r3 = p->r3, r5 = p->r5;
    const real_t a3 = p->a3, a5 = p->a5, b3 = p->b3, c3 = p->c3;
    const real_t uy = p->uy, vy = p->vy, uz = p->uz, vz = p->vz;
    const real_t wy = p->wy, wz = p->wz;
    real_t du[12];
    const real_t cc = d + z;
    const real_t rd = r + d;
    const real_t d12 = 1.0 / (r * rd * rd);
    const real_t d32 = d12 * (2.0 * r + d) / r2;
    const real_t d33 = d12 * (3.0 * r + d) / (r2 * rd);
    const real_t d53 = d12 * (8.0 * r2 + 9.0 * r * d + 3.0 * d2) / (r2 * r2 * rd);
    const real_t d54 = d12 * (5.0 * r2 + 4.0 * r * d + d2) / r3 * d12;
    const real_t fi1 = y * (d12 - x2 * d33);
    const real_t fi2 = x * (d12 - y2 * d33);
    const real_t fi3 = x / r3 - fi2;
    const real_t fi4 = -xy * d32;
    const real_t fi5 = 1.0 / (r * rd) - x2 * d32;
    const real_t fj1 = -3.0 * xy * (d33 - x2 * d54);
    const real_t fj2 = 1.0 / r3 - 3.0 * d12 + 3.0 * x2 * y2 * d54;
    const real_t fj3 = a3 / r3 - fj2;
    const real_t fj4 = -3.0 * xy / r5 - fj1;
    const real_t fk1 = -y * (d32 - x2 * d53);
    const real_t fk2 = -x * (d32 - y2 * d53);
    const real_t fk3 = -3.0 * x * d / r5 - fk2;

    okada_zero_u(u);
    // 走向滑动贡献
    if(pot1 != 0.0){
        du[0] = -x2 * qr - alp3 * fi1 * sd;
        du[1] = -xy * qr - alp3 * fi2 * sd;
        du[2] = -cc * x * qr - alp3 * fi4 * sd;
        du[3] = -x * qr * (1.0 + a5) - alp3 * fj1 * sd;
        du[4] = -y * qr * a5 - alp3 * fj2 * sd;
        du[5] = -cc * qr * a5 - alp3 * fk1 * sd;
        du[6] = -3.0 * x2 / r5 * uy - alp3 * fj2 * sd;
        du[7] = -3.0 * xy / r5 * uy - x * qr - alp3 * fj4 * sd;
        du[8] = -3.0 * cc * x / r5 * uy - alp3 * fk2 * sd;
        du[9] = -3.0 * x2 / r5 * uz + alp3 * fk1 * sd;
        du[10] = -3.0 * xy / r5 * uz + alp3 * fk2 * sd;
        du[11] = 3.0 * x / r5 * (-cc * uz + alp3 * y * sd);
        okada_add_scaled(u, du, pot1 / OKADA_PI2);
    }

    // 倾向滑动贡献
    if(pot2 != 0.0){
        du[0] = -x * p0 * qr + alp3 * fi3 * sdcd;
        du[1] = -y * p0 * qr + alp3 * fi1 * sdcd;
        du[2] = -cc * p0 * qr + alp3 * fi5 * sdcd;
        du[3] = -p0 * qr * a5 + alp3 * fj3 * sdcd;
        du[4] = y * p0 * qrx + alp3 * fj1 * sdcd;
        du[5] = cc * p0 * qrx + alp3 * fk3 * sdcd;
        du[6] = -3.0 * x / r5 * vy + alp3 * fj1 * sdcd;
        du[7] = -3.0 * y / r5 * vy - p0 * qr + alp3 * fj2 * sdcd;
        du[8] = -3.0 * cc / r5 * vy + alp3 * fk1 * sdcd;
        du[9] = -3.0 * x / r5 * vz - alp3 * fk3 * sdcd;
        du[10] = -3.0 * y / r5 * vz - alp3 * fk1 * sdcd;
        du[11] = -3.0 * cc / r5 * vz + alp3 * a3 / r3 * sdcd;
        okada_add_scaled(u, du, pot2 / OKADA_PI2);
    }

    // 张裂贡献
    if(pot3 != 0.0){
        du[0] = x * q * qr - alp3 * fi3 * sdsd;
        du[1] = y * q * qr - alp3 * fi1 * sdsd;
        du[2] = cc * q * qr - alp3 * fi5 * sdsd;
        du[3] = q * qr * a5 - alp3 * fj3 * sdsd;
        du[4] = -y * q * qrx - alp3 * fj1 * sdsd;
        du[5] = -cc * q * qrx - alp3 * fk3 * sdsd;
        du[6] = x * qr * wy - alp3 * fj1 * sdsd;
        du[7] = qr * (y * wy + q) - alp3 * fj2 * sdsd;
        du[8] = cc * qr * wy - alp3 * fk1 * sdsd;
        du[9] = x * qr * wz + alp3 * fk3 * sdsd;
        du[10] = y * qr * wz + alp3 * fk1 * sdsd;
        du[11] = cc * qr * wz - alp3 * a3 / r3 * sdsd;
        okada_add_scaled(u, du, pot3 / OKADA_PI2);
    }

    // 爆炸源贡献
    if(pot4 != 0.0){
        du[0] = alp3 * x / r3;
        du[1] = alp3 * y / r3;
        du[2] = alp3 * d / r3;
        du[3] = alp3 * a3 / r3;
        du[4] = -alp3 * 3.0 * xy / r5;
        du[5] = -alp3 * 3.0 * x * d / r5;
        du[6] = du[4];
        du[7] = alp3 * b3 / r3;
        du[8] = -alp3 * 3.0 * y * d / r5;
        du[9] = -du[5];
        du[10] = -du[8];
        du[11] = -alp3 * c3 / r3;
        okada_add_scaled(u, du, pot4 / OKADA_PI2);
    }
}

/** 计算点源镜像源的 UC0 公式贡献
 *
 * @param[in]   c       Okada 介质常数
 * @param[in]   p       点源几何量
 * @param[in]   z       观测点 Z 坐标
 * @param[in]   pot1    走向滑动 potency
 * @param[in]   pot2    倾向滑动 potency
 * @param[in]   pot3    张裂 potency
 * @param[in]   pot4    爆炸源 potency
 * @param[out]  u       12 个位移及偏导结果
 */
static void UC0(const OKADA_MEDIUM *c, const OKADA_POINT *p, real_t z,
    real_t pot1, real_t pot2, real_t pot3, real_t pot4, real_t u[12])
{
    const real_t alp4 = c->alp4, alp5 = c->alp5;
    const real_t sd = c->sd, cd = c->cd, sdcd = c->sdcd, sdsd = c->sdsd;
    const real_t s2d = c->s2d, c2d = c->c2d;
    const real_t x = p->x, y = p->y, d = p->d, p0 = p->p, q = p->q;
    const real_t s = p->s, t = p->t, xy = p->xy, x2 = p->x2;
    const real_t y2 = p->y2, d2 = p->d2, r2 = p->r2, r3 = p->r3, r5 = p->r5;
    const real_t a3 = p->a3, a5 = p->a5, c3 = p->c3, qr = p->qr, qrx = p->qrx;
    real_t du[12];
    const real_t cc = d + z;
    const real_t q2 = q * q;
    const real_t r7 = r5 * r2;
    const real_t a7 = 1.0 - 7.0 * x2 / r2;
    const real_t b5 = 1.0 - 5.0 * y2 / r2;
    const real_t b7 = 1.0 - 7.0 * y2 / r2;
    const real_t c5 = 1.0 - 5.0 * d2 / r2;
    const real_t c7 = 1.0 - 7.0 * d2 / r2;
    const real_t d7 = 2.0 - 7.0 * q2 / r2;
    const real_t qr5 = 5.0 * q / r2;
    const real_t qr7 = 7.0 * q / r2;
    const real_t dr5 = 5.0 * d / r2;

    okada_zero_u(u);
    // 走向滑动贡献
    if(pot1 != 0.0){
        du[0] = -alp4 * a3 / r3 * cd + alp5 * cc * qr * a5;
        du[1] = 3.0 * x / r5 * (alp4 * y * cd + alp5 * cc * (sd - y * qr5));
        du[2] = 3.0 * x / r5 * (-alp4 * y * sd + alp5 * cc * (cd + d * qr5));
        du[3] = alp4 * 3.0 * x / r5 * (2.0 + a5) * cd - alp5 * cc * qrx * (2.0 + a7);
        du[4] = 3.0 / r5 * (alp4 * y * a5 * cd + alp5 * cc * (a5 * sd - y * qr5 * a7));
        du[5] = 3.0 / r5 * (-alp4 * y * a5 * sd + alp5 * cc * (a5 * cd + d * qr5 * a7));
        du[6] = du[4];
        du[7] = 3.0 * x / r5 * (alp4 * b5 * cd
            - alp5 * 5.0 * cc / r2 * (2.0 * y * sd + q * b7));
        du[8] = 3.0 * x / r5 * (-alp4 * b5 * sd
            + alp5 * 5.0 * cc / r2 * (d * b7 * sd - y * c7 * cd));
        du[9] = 3.0 / r5 * (-alp4 * d * a5 * cd + alp5 * cc * (a5 * cd + d * qr5 * a7));
        du[10] = 15.0 * x / r7 * (alp4 * y * d * cd + alp5 * cc * (d * b7 * sd - y * c7 * cd));
        du[11] = 15.0 * x / r7 * (-alp4 * y * d * sd + alp5 * cc * (2.0 * d * cd - q * c7));
        okada_add_scaled(u, du, pot1 / OKADA_PI2);
    }

    // 倾向滑动贡献
    if(pot2 != 0.0){
        du[0] = alp4 * 3.0 * x * t / r5 - alp5 * cc * p0 * qrx;
        du[1] = -alp4 / r3 * (c2d - 3.0 * y * t / r2)
            + alp5 * 3.0 * cc / r5 * (s - y * p0 * qr5);
        du[2] = -alp4 * a3 / r3 * sdcd + alp5 * 3.0 * cc / r5 * (t + d * p0 * qr5);
        du[3] = alp4 * 3.0 * t / r5 * a5 - alp5 * 5.0 * cc * p0 * qr / r2 * a7;
        du[4] = 3.0 * x / r5 * (alp4 * (c2d - 5.0 * y * t / r2)
            - alp5 * 5.0 * cc / r2 * (s - y * p0 * qr7));
        du[5] = 3.0 * x / r5 * (alp4 * (2.0 + a5) * sdcd
            - alp5 * 5.0 * cc / r2 * (t + d * p0 * qr7));
        du[6] = du[4];
        du[7] = 3.0 / r5 * (alp4 * (2.0 * y * c2d + t * b5)
            + alp5 * cc * (s2d - 10.0 * y * s / r2 - p0 * qr5 * b7));
        du[8] = 3.0 / r5 * (alp4 * y * a5 * sdcd
            - alp5 * cc * ((3.0 + a5) * c2d + y * p0 * dr5 * qr7));
        du[9] = 3.0 * x / r5 * (-alp4 * (s2d - t * dr5)
            - alp5 * 5.0 * cc / r2 * (t + d * p0 * qr7));
        du[10] = 3.0 / r5 * (-alp4 * (d * b5 * c2d + y * c5 * s2d)
            - alp5 * cc * ((3.0 + a5) * c2d + y * p0 * dr5 * qr7));
        du[11] = 3.0 / r5 * (-alp4 * d * a5 * sdcd
            - alp5 * cc * (s2d - 10.0 * d * t / r2 + p0 * qr5 * c7));
        okada_add_scaled(u, du, pot2 / OKADA_PI2);
    }

    // 张裂贡献
    if(pot3 != 0.0){
        du[0] = 3.0 * x / r5 * (-alp4 * s + alp5 * (cc * q * qr5 - z));
        du[1] = alp4 / r3 * (s2d - 3.0 * y * s / r2)
            + alp5 * 3.0 / r5 * (cc * (t - y + y * q * qr5) - y * z);
        du[2] = -alp4 / r3 * (1.0 - a3 * sdsd)
            - alp5 * 3.0 / r5 * (cc * (s - d + d * q * qr5) - d * z);
        du[3] = -alp4 * 3.0 * s / r5 * a5 + alp5 * (cc * qr * qr5 * a7 - 3.0 * z / r5 * a5);
        du[4] = 3.0 * x / r5 * (-alp4 * (s2d - 5.0 * y * s / r2)
            - alp5 * 5.0 / r2 * (cc * (t - y + y * q * qr7) - y * z));
        du[5] = 3.0 * x / r5 * (alp4 * (1.0 - (2.0 + a5) * sdsd)
            + alp5 * 5.0 / r2 * (cc * (s - d + d * q * qr7) - d * z));
        du[6] = du[4];
        du[7] = 3.0 / r5 * (-alp4 * (2.0 * y * s2d + s * b5)
            - alp5 * (cc * (2.0 * sdsd + 10.0 * y * (t - y) / r2 - q * qr5 * b7) + z * b5));
        du[8] = 3.0 / r5 * (alp4 * y * (1.0 - a5 * sdsd)
            + alp5 * (cc * (3.0 + a5) * s2d - y * dr5 * (cc * d7 + z)));
        du[9] = 3.0 * x / r5 * (-alp4 * (c2d + s * dr5)
            + alp5 * (5.0 * cc / r2 * (s - d + d * q * qr7) - 1.0 - z * dr5));
        du[10] = 3.0 / r5 * (alp4 * (d * b5 * s2d - y * c5 * c2d)
            + alp5 * (cc * ((3.0 + a5) * s2d - y * dr5 * d7) - y * (1.0 + z * dr5)));
        du[11] = 3.0 / r5 * (-alp4 * d * (1.0 - a5 * sdsd)
            - alp5 * (cc * (c2d + 10.0 * d * (s - d) / r2 - q * qr5 * c7) + z * (1.0 + c5)));
        okada_add_scaled(u, du, pot3 / OKADA_PI2);
    }

    // 爆炸源贡献
    if(pot4 != 0.0){
        du[0] = alp4 * 3.0 * x * d / r5;
        du[1] = alp4 * 3.0 * y * d / r5;
        du[2] = alp4 * c3 / r3;
        du[3] = alp4 * 3.0 * d / r5 * a5;
        du[4] = -alp4 * 15.0 * xy * d / r7;
        du[5] = -alp4 * 3.0 * x / r5 * c5;
        du[6] = du[4];
        du[7] = alp4 * 3.0 * d / r5 * b5;
        du[8] = -alp4 * 3.0 * y / r5 * c5;
        du[9] = du[5];
        du[10] = du[8];
        du[11] = alp4 * 3.0 * d / r5 * (2.0 + c5);
        okada_add_scaled(u, du, pot4 / OKADA_PI2);
    }
}

/** 组合点源真实源和镜像源贡献
 *
 * @param[in]   c       Okada 介质常数
 * @param[in]   x       局部 X 坐标
 * @param[in]   y       局部 Y 坐标
 * @param[in]   z       观测点 Z 坐标
 * @param[in]   depth   点源深度
 * @param[in]   pot1    走向滑动 potency
 * @param[in]   pot2    倾向滑动 potency
 * @param[in]   pot3    张裂 potency
 * @param[in]   pot4    爆炸源 potency
 * @param[out]  out     内部 12 元结果数组
 * @return              0 表示正常，1 表示奇异点
 */
static int dc3d0_internal(const OKADA_MEDIUM *c, real_t x, real_t y, real_t z, real_t depth,
    real_t pot1, real_t pot2, real_t pot3, real_t pot4, real_t out[12])
{
    OKADA_POINT p;
    real_t dua[12], dub[12], duc[12];
    real_t zz = z;
    okada_zero_u(out);

    // 真实源贡献，Fortran 中对应 DD=DEPTH+Z
    DCCON1(c, &p, x, y, depth + z);
    if(p.r == 0.0) return 1;
    UA0(c, &p, pot1, pot2, pot3, pot4, dua);
    for(int i = 0; i < 12; ++i){
        if(i < 9) out[i] -= dua[i];
        else out[i] += dua[i];
    }

    // 镜像源贡献，Fortran 中对应 DD=DEPTH-Z
    DCCON1(c, &p, x, y, depth - z);
    UA0(c, &p, pot1, pot2, pot3, pot4, dua);
    UB0(c, &p, zz, pot1, pot2, pot3, pot4, dub);
    UC0(c, &p, zz, pot1, pot2, pot3, pot4, duc);
    for(int i = 0; i < 12; ++i){
        real_t du = dua[i] + dub[i] + zz * duc[i];
        if(i >= 9) du += duc[i - 9];
        out[i] += du;
    }
    return 0;
}

/** 计算点源静态位移和位移偏导
 *
 * 对应 Fortran 源码中的 DC3D0，公开结果不在此处进行单位换算
 *
 * @param[in]   alpha   均匀介质常数
 * @param[in]   x       局部 X 坐标
 * @param[in]   y       局部 Y 坐标
 * @param[in]   z       观测点 Z 坐标，向上为正
 * @param[in]   depth   点源深度
 * @param[in]   dip     点源倾角，单位为度
 * @param[in]   pot1    走向滑动 potency
 * @param[in]   pot2    倾向滑动 potency
 * @param[in]   pot3    张裂 potency
 * @param[in]   pot4    爆炸源 potency
 * @param[out]  u       局部 X、Y、Z 位移
 * @param[out]  upar    位移偏导 [偏导坐标][位移分量]
 * @return              0 表示正常，1 表示奇异点，2 表示输入观测点 Z 为正
 */
int grt_okada_dc3d0(
    real_t alpha, real_t x, real_t y, real_t z, real_t depth, real_t dip,
    real_t pot1, real_t pot2, real_t pot3, real_t pot4,
    real_t u[3], real_t upar[3][3])
{
    OKADA_MEDIUM c;
    real_t out[12];
    int iret = 0;

    okada_zero_u(out);
    if(z > 0.0){
        iret = 2;
    } else {
        DCCON0(&c, alpha, dip);
        iret = dc3d0_internal(&c, x, y, z, depth, pot1, pot2, pot3, pot4, out);
    }
    okada_copy_output(out, u, upar);
    return iret;
}

/** 初始化有限断层角点公式使用的几何量
 *
 * 对应 Fortran 源码中的 DCCON2
 *
 * @param[in]   c       Okada 介质常数
 * @param[out]  f       有限断层几何量
 * @param[in]   xi      角点相对 X 坐标
 * @param[in]   et      角点相对倾向坐标
 * @param[in]   q       角点到断层面的法向距离
 * @param[in]   kxi     X 方向对数项分支标志
 * @param[in]   ket     倾向方向对数项分支标志
 */
static void DCCON2(const OKADA_MEDIUM *c, OKADA_FINITE *f,
    real_t xi, real_t et, real_t q, int kxi, int ket)
{
    if(fabs(xi) < OKADA_EPS) xi = 0.0;
    if(fabs(et) < OKADA_EPS) et = 0.0;
    if(fabs(q) < OKADA_EPS) q = 0.0;

    f->xi = xi;
    f->et = et;
    f->q = q;
    f->xi2 = xi * xi;
    f->et2 = et * et;
    f->q2 = q * q;
    f->r2 = f->xi2 + f->et2 + f->q2;
    f->r = sqrt(f->r2);
    if(f->r == 0.0) return;
    f->r3 = f->r * f->r2;
    f->r5 = f->r3 * f->r2;
    f->y = et * c->cd + q * c->sd;
    f->d = et * c->sd - q * c->cd;

    f->tt = (q == 0.0) ? 0.0 : atan(xi * et / (q * f->r));

    if(kxi){
        f->alx = -log(f->r - xi);
        f->x11 = 0.0;
        f->x32 = 0.0;
    } else {
        real_t rxi = f->r + xi;
        f->alx = log(rxi);
        f->x11 = 1.0 / (f->r * rxi);
        f->x32 = (f->r + rxi) * f->x11 * f->x11 / f->r;
    }

    if(ket){
        f->ale = -log(f->r - et);
        f->y11 = 0.0;
        f->y32 = 0.0;
    } else {
        real_t ret = f->r + et;
        f->ale = log(ret);
        f->y11 = 1.0 / (f->r * ret);
        f->y32 = (f->r + ret) * f->y11 * f->y11 / f->r;
    }

    f->ey = c->sd/f->r - f->y*q/f->r3;
    f->ez = c->cd/f->r + f->d*q/f->r3;
    f->fy = f->d/f->r3 + f->xi2*f->y32*c->sd;
    f->fz = f->y/f->r3 + f->xi2*f->y32*c->cd;
    f->gy = 2.0*f->x11*c->sd - f->y*q*f->x32;
    f->gz = 2.0*f->x11*c->cd + f->d*q*f->x32;
    f->hy = f->d*q*f->x32 + xi*q*f->y32*c->sd;
    f->hz = f->y*q*f->x32 + xi*q*f->y32*c->cd;
}

/** 计算有限断层真实源项的 UA 公式贡献
 *
 * @param[in]   c       Okada 介质常数
 * @param[in]   f       有限断层角点几何量
 * @param[in]   disl1   走向滑动位错
 * @param[in]   disl2   倾向滑动位错
 * @param[in]   disl3   张裂位错
 * @param[out]  u       12 个位移及偏导结果
 */
static void UA(const OKADA_MEDIUM *c, const OKADA_FINITE *f,
    real_t disl1, real_t disl2, real_t disl3, real_t u[12])
{
    const real_t alp1 = c->alp1, alp2 = c->alp2, sd = c->sd, cd = c->cd;
    const real_t xi = f->xi, et = f->et, q = f->q, q2 = f->q2;
    const real_t r = f->r, r3 = f->r3, y = f->y, d = f->d;
    const real_t tt = f->tt, alx = f->alx, ale = f->ale;
    const real_t x11 = f->x11, y11 = f->y11, y32 = f->y32;
    const real_t ey = f->ey, ez = f->ez, fy = f->fy, fz = f->fz;
    const real_t gy = f->gy, gz = f->gz, hy = f->hy, hz = f->hz;
    real_t du[12];
    const real_t xy = xi * y11;
    const real_t qx = q * x11;
    const real_t qy = q * y11;
    okada_zero_u(u);

    // 走向滑动贡献
    if(disl1 != 0.0){
        du[0] = 0.5 * tt + alp2 * xi * qy;
        du[1] = alp2 * q / r;
        du[2] = alp1 * ale - alp2 * q * qy;
        du[3] = -alp1 * qy - alp2 * xi * xi * q * y32;
        du[4] = -alp2 * xi * q / r3;
        du[5] = alp1 * xy + alp2 * xi * q2 * y32;
        du[6] = alp1 * xy * sd + alp2 * xi * fy + 0.5 * d * x11;
        du[7] = alp2 * ey;
        du[8] = alp1 * (cd / r + qy * sd) - alp2 * q * fy;
        du[9] = alp1 * xy * cd + alp2 * xi * fz + 0.5 * y * x11;
        du[10] = alp2 * ez;
        du[11] = -alp1 * (sd / r - qy * cd) - alp2 * q * fz;
        okada_add_scaled(u, du, disl1 / OKADA_PI2);
    }

    // 倾向滑动贡献
    if(disl2 != 0.0){
        du[0] = alp2 * q / r;
        du[1] = 0.5 * tt + alp2 * et * qx;
        du[2] = alp1 * alx - alp2 * q * qx;
        du[3] = -alp2 * xi * q / r3;
        du[4] = -0.5 * qy - alp2 * et * q / r3;
        du[5] = alp1 / r + alp2 * q2 / r3;
        du[6] = alp2 * ey;
        du[7] = alp1 * d * x11 + 0.5 * xy * sd + alp2 * et * gy;
        du[8] = alp1 * y * x11 - alp2 * q * gy;
        du[9] = alp2 * ez;
        du[10] = alp1 * y * x11 + 0.5 * xy * cd + alp2 * et * gz;
        du[11] = -alp1 * d * x11 - alp2 * q * gz;
        okada_add_scaled(u, du, disl2 / OKADA_PI2);
    }

    // 张裂贡献
    if(disl3 != 0.0){
        du[0] = -alp1 * ale - alp2 * q * qy;
        du[1] = -alp1 * alx - alp2 * q * qx;
        du[2] = 0.5 * tt - alp2 * (et * qx + xi * qy);
        du[3] = -alp1 * xy + alp2 * xi * q2 * y32;
        du[4] = -alp1 / r + alp2 * q2 / r3;
        du[5] = -alp1 * qy - alp2 * q * q2 * y32;
        du[6] = -alp1 * (cd / r + qy * sd) - alp2 * q * fy;
        du[7] = -alp1 * y * x11 - alp2 * q * gy;
        du[8] = alp1 * (d * x11 + xy * sd) + alp2 * q * hy;
        du[9] = alp1 * (sd / r - qy * cd) - alp2 * q * fz;
        du[10] = alp1 * d * x11 - alp2 * q * gz;
        du[11] = alp1 * (y * x11 + xy * cd) + alp2 * q * hz;
        okada_add_scaled(u, du, disl3 / OKADA_PI2);
    }
}

/** 计算有限断层镜像源的 UB 公式贡献
 *
 * @param[in]   c       Okada 介质常数
 * @param[in]   f       有限断层角点几何量
 * @param[in]   disl1   走向滑动位错
 * @param[in]   disl2   倾向滑动位错
 * @param[in]   disl3   张裂位错
 * @param[out]  u       12 个位移及偏导结果
 */
static void UB(const OKADA_MEDIUM *c, const OKADA_FINITE *f,
    real_t disl1, real_t disl2, real_t disl3, real_t u[12])
{
    const real_t alp3 = c->alp3;
    const real_t sd = c->sd, cd = c->cd, sdcd = c->sdcd, sdsd = c->sdsd, cdcd = c->cdcd;
    const real_t xi = f->xi, et = f->et, q = f->q, xi2 = f->xi2, q2 = f->q2;
    const real_t r = f->r, y = f->y, d = f->d, tt = f->tt;
    const real_t ale = f->ale, x11 = f->x11, y11 = f->y11;
    const real_t y32 = f->y32, ey = f->ey, ez = f->ez;
    const real_t fy = f->fy, fz = f->fz, gy = f->gy, gz = f->gz;
    const real_t hy = f->hy, hz = f->hz, r3 = f->r3;
    real_t du[12];
    const real_t rd = r + d;
    const real_t d11 = 1.0 / (r * rd);
    const real_t aj2 = xi * y / rd * d11;
    const real_t aj5 = -(d + y * y / rd) * d11;
    real_t ai3, ai4, ak1, ak3, aj3, aj6;
    real_t ai1, ai2, ak2, ak4, aj1, aj4;
    const real_t xy = xi * y11;
    const real_t qx = q * x11;
    const real_t qy = q * y11;

    // 根据倾角选择一般情况或垂直断层的辅助量公式
    if(cd != 0.0){
        if(xi == 0.0){
            ai4 = 0.0;
        } else {
            real_t x0 = sqrt(xi2 + q2);
            ai4 = 1.0 / cdcd * (
                xi / rd * sdcd
                + 2.0 * atan((et * (x0 + q * cd) + x0 * (r + x0) * sd)
                    / (xi * (r + x0) * cd)));
        }
        ai3 = (y * cd / rd - ale + sd * log(rd)) / cdcd;
        ak1 = xi * (d11 - y11 * sd) / cd;
        ak3 = (q * y11 - y * d11) / cd;
        aj3 = (ak1 - aj2 * sd) / cd;
        aj6 = (ak3 - aj5 * sd) / cd;
    } else {
        real_t rd2 = rd*rd;
        ai3 = (et / rd + y * q / rd2 - ale) / 2.0;
        ai4 = xi * y / rd2 / 2.0;
        ak1 = xi * q / rd * d11;
        ak3 = sd / rd * (xi2 * d11 - 1.0);
        aj3 = -xi / rd2 * (q2 * d11 - 0.5);
        aj6 = -y / rd2 * (xi2 * d11 - 0.5);
    }

    ai1 = -xi / rd * cd - ai4 * sd;
    ai2 = log(rd) + ai3 * sd;
    ak2 = 1.0 / r + ak3 * sd;
    ak4 = xy * cd - ak1 * sd;
    aj1 = aj5 * cd - aj6 * sd;
    aj4 = -xy - aj2 * cd + aj3 * sd;

    okada_zero_u(u);
    // 走向滑动贡献
    if(disl1 != 0.0){
        du[0] = -xi * qy - tt - alp3 * ai1 * sd;
        du[1] = -q / r + alp3 * y / rd * sd;
        du[2] = q * qy - alp3 * ai2 * sd;
        du[3] = xi2 * q * y32 - alp3 * aj1 * sd;
        du[4] = xi * q / r3 - alp3 * aj2 * sd;
        du[5] = -xi * q2 * y32 - alp3 * aj3 * sd;
        du[6] = -xi * fy - d * x11 + alp3 * (xy + aj4) * sd;
        du[7] = -ey + alp3 * (1.0 / r + aj5) * sd;
        du[8] = q * fy - alp3 * (qy - aj6) * sd;
        du[9] = -xi * fz - y * x11 + alp3 * ak1 * sd;
        du[10] = -ez + alp3 * y * d11 * sd;
        du[11] = q * fz + alp3 * ak2 * sd;
        okada_add_scaled(u, du, disl1 / OKADA_PI2);
    }

    // 倾向滑动贡献
    if(disl2 != 0.0){
        du[0] = -q / r + alp3 * ai3 * sdcd;
        du[1] = -et * qx - tt - alp3 * xi / rd * sdcd;
        du[2] = q * qx + alp3 * ai4 * sdcd;
        du[3] = xi * q / r3 + alp3 * aj4 * sdcd;
        du[4] = et * q / r3 + qy + alp3 * aj5 * sdcd;
        du[5] = -q2 / r3 + alp3 * aj6 * sdcd;
        du[6] = -ey + alp3 * aj1 * sdcd;
        du[7] = -et * gy - xy * sd + alp3 * aj2 * sdcd;
        du[8] = q * gy + alp3 * aj3 * sdcd;
        du[9] = -ez - alp3 * ak3 * sdcd;
        du[10] = -et * gz - xy * cd - alp3 * xi * d11 * sdcd;
        du[11] = q * gz - alp3 * ak4 * sdcd;
        okada_add_scaled(u, du, disl2 / OKADA_PI2);
    }

    // 张裂贡献
    if(disl3 != 0.0){
        du[0] = q * qy - alp3 * ai3 * sdsd;
        du[1] = q * qx + alp3 * xi / rd * sdsd;
        du[2] = et * qx + xi * qy - tt - alp3 * ai4 * sdsd;
        du[3] = -xi * q2 * y32 - alp3 * aj4 * sdsd;
        du[4] = -q2 / r3 - alp3 * aj5 * sdsd;
        du[5] = q * q2 * y32 - alp3 * aj6 * sdsd;
        du[6] = q * fy - alp3 * aj1 * sdsd;
        du[7] = q * gy - alp3 * aj2 * sdsd;
        du[8] = -q * hy - alp3 * aj3 * sdsd;
        du[9] = q * fz + alp3 * ak3 * sdsd;
        du[10] = q * gz + alp3 * xi * d11 * sdsd;
        du[11] = -q * hz + alp3 * ak4 * sdsd;
        okada_add_scaled(u, du, disl3 / OKADA_PI2);
    }
}

/** 计算有限断层镜像源的 UC 公式贡献
 *
 * @param[in]   c       Okada 介质常数
 * @param[in]   f       有限断层角点几何量
 * @param[in]   z       观测点 Z 坐标
 * @param[in]   disl1   走向滑动位错
 * @param[in]   disl2   倾向滑动位错
 * @param[in]   disl3   张裂位错
 * @param[out]  u       12 个位移及偏导结果
 */
static void UC(const OKADA_MEDIUM *c, const OKADA_FINITE *f, real_t z,
    real_t disl1, real_t disl2, real_t disl3, real_t u[12])
{
    const real_t alp4 = c->alp4, alp5 = c->alp5;
    const real_t sd = c->sd, cd = c->cd, sdcd = c->sdcd, cdcd = c->cdcd, sdsd = c->sdsd;
    const real_t xi = f->xi, et = f->et, q = f->q, xi2 = f->xi2, et2 = f->et2;
    const real_t r = f->r, r2 = f->r2, r3 = f->r3, r5 = f->r5;
    const real_t y = f->y, d = f->d, x11 = f->x11, y11 = f->y11;
    const real_t x32 = f->x32, y32 = f->y32;
    real_t du[12];
    const real_t cc = d + z;
    const real_t x53 = (8.0 * r2 + 9.0 * r * xi + 3.0 * xi2) * x11 * x11 * x11 / r2;
    const real_t y53 = (8.0 * r2 + 9.0 * r * et + 3.0 * et2) * y11 * y11 * y11 / r2;
    const real_t h = q * cd - z;
    const real_t z32 = sd / r3 - h * y32;
    const real_t z53 = 3.0 * sd / r5 - h * y53;
    const real_t y0 = y11 - xi2 * y32;
    const real_t z0 = z32 - xi2 * z53;
    const real_t ppy = cd / r3 + q * y32 * sd;
    const real_t ppz = sd / r3 - q * y32 * cd;
    const real_t qq = z * y32 + z32 + z0;
    const real_t qqy = 3.0 * cc * d / r5 - qq * sd;
    const real_t qqz = 3.0 * cc * y / r5 - qq * cd + q * y32;
    const real_t xy = xi * y11;
    const real_t qy = q * y11;
    const real_t qr = 3.0 * q / r5;
    const real_t cdr = (cc + d) / r3;
    const real_t yy0 = y / r3 - y0 * cd;

    okada_zero_u(u);
    // 走向滑动贡献
    if(disl1 != 0.0){
        du[0] = alp4 * xy * cd - alp5 * xi * q * z32;
        du[1] = alp4 * (cd / r + 2.0 * qy * sd) - alp5 * cc * q / r3;
        du[2] = alp4 * qy * cd - alp5 * (cc * et / r3 - z * y11 + xi2 * z32);
        du[3] = alp4 * y0 * cd - alp5 * q * z0;
        du[4] = -alp4 * xi * (cd / r3 + 2.0 * q * y32 * sd) + alp5 * cc * xi * qr;
        du[5] = -alp4 * xi * q * y32 * cd + alp5 * xi * (3.0 * cc * et / r5 - qq);
        du[6] = -alp4 * xi * ppy * cd - alp5 * xi * qqy;
        du[7] = alp4 * 2.0 * (d / r3 - y0 * sd) * sd - y / r3 * cd
            - alp5 * (cdr * sd - et / r3 - cc * y * qr);
        du[8] = -alp4 * q / r3 + yy0 * sd
            + alp5 * (cdr * cd + cc * d * qr - (y0 * cd + q * z0) * sd);
        du[9] = alp4 * xi * ppz * cd - alp5 * xi * qqz;
        du[10] = alp4 * 2.0 * (y / r3 - y0 * cd) * sd + d / r3 * cd
            - alp5 * (cdr * cd + cc * d * qr);
        du[11] = yy0 * cd - alp5 * (cdr * sd - cc * y * qr - y0 * sdsd + q * z0 * cd);
        okada_add_scaled(u, du, disl1 / OKADA_PI2);
    }

    // 倾向滑动贡献
    if(disl2 != 0.0){
        du[0] = alp4 * cd / r - qy * sd - alp5 * cc * q / r3;
        du[1] = alp4 * y * x11 - alp5 * cc * et * q * x32;
        du[2] = -d * x11 - xy * sd - alp5 * cc * (x11 - q * q * x32);
        du[3] = -alp4 * xi / r3 * cd + alp5 * cc * xi * qr + xi * q * y32 * sd;
        du[4] = -alp4 * y / r3 + alp5 * cc * et * qr;
        du[5] = d / r3 - y0 * sd + alp5 * cc / r3 * (1.0 - 3.0 * q * q / r2);
        du[6] = -alp4 * et / r3 + y0 * sdsd - alp5 * (cdr * sd - cc * y * qr);
        du[7] = alp4 * (x11 - y * y * x32)
            - alp5 * cc * ((d + 2.0 * q * cd) * x32 - y * et * q * x53);
        du[8] = xi * ppy * sd + y * d * x32
            + alp5 * cc * ((y + 2.0 * q * sd) * x32 - y * q * q * x53);
        du[9] = -q / r3 + y0 * sdcd - alp5 * (cdr * cd + cc * d * qr);
        du[10] = alp4 * y * d * x32
            - alp5 * cc * ((y - 2.0 * q * sd) * x32 + d * et * q * x53);
        du[11] = -xi * ppz * sd + x11 - d * d * x32
            - alp5 * cc * ((d - 2.0 * q * cd) * x32 - d * q * q * x53);
        okada_add_scaled(u, du, disl2 / OKADA_PI2);
    }

    // 张裂贡献
    if(disl3 != 0.0){
        du[0] = -alp4 * (sd / r + qy * cd) - alp5 * (z * y11 - q * q * z32);
        du[1] = alp4 * 2.0 * xy * sd + d * x11 - alp5 * cc * (x11 - q * q * x32);
        du[2] = alp4 * (y * x11 + xy * cd) + alp5 * q * (cc * et * x32 + xi * z32);
        du[3] = alp4 * xi / r3 * sd + xi * q * y32 * cd
            + alp5 * xi * (3.0 * cc * et / r5 - 2.0 * z32 - z0);
        du[4] = alp4 * 2.0 * y0 * sd - d / r3 + alp5 * cc / r3 * (1.0 - 3.0 * q * q / r2);
        du[5] = -alp4 * yy0 - alp5 * (cc * et * qr - q * z0);
        du[6] = alp4 * (q / r3 + y0 * sdcd)
            + alp5 * (z / r3 * cd + cc * d * qr - q * z0 * sd);
        du[7] = -alp4 * 2.0 * xi * ppy * sd - y * d * x32
            + alp5 * cc * ((y + 2.0 * q * sd) * x32 - y * q * q * x53);
        du[8] = -alp4 * (xi * ppy * cd - x11 + y * y * x32)
            + alp5 * (cc * ((d + 2.0 * q * cd) * x32 - y * et * q * x53) + xi * qqy);
        du[9] = -et / r3 + y0 * cdcd
            - alp5 * (z / r3 * sd - cc * y * qr - y0 * sdsd + q * z0 * cd);
        du[10] = alp4 * 2.0 * xi * ppz * sd - x11 + d * d * x32
            - alp5 * cc * ((d - 2.0 * q * cd) * x32 - d * q * q * x53);
        du[11] = alp4 * (xi * ppz * cd + y * d * x32)
            + alp5 * (cc * ((y - 2.0 * q * sd) * x32 + d * et * q * x53) + xi * qqz);
        okada_add_scaled(u, du, disl3 / OKADA_PI2);
    }
}

/** 判断有限断层角点是否落在奇异位置
 *
 * @param[in]   q       观测点到断层面的法向距离
 * @param[in]   xi1     X 方向第一个角点坐标
 * @param[in]   xi2     X 方向第二个角点坐标
 * @param[in]   et1     倾向方向第一个角点坐标
 * @param[in]   et2     倾向方向第二个角点坐标
 * @return              角点位于奇异位置时返回 true
 */
static bool okada_finite_singular(real_t q, real_t xi1, real_t xi2, real_t et1, real_t et2)
{
    return q == 0.0 && (
        (xi1*xi2 <= 0.0 && et1*et2 == 0.0)
        || (et1*et2 <= 0.0 && xi1*xi2 == 0.0));
}

/** 累加有限断层真实源的四个角点贡献
 *
 * @param[in]       c       Okada 介质常数
 * @param[in]       xi      两个 X 方向角点坐标
 * @param[in]       et      两个倾向方向角点坐标
 * @param[in]       q       观测点到断层面的法向距离
 * @param[in]       kxi     X 方向对数项分支标志
 * @param[in]       ket     倾向方向对数项分支标志
 * @param[in]       disl1   走向滑动位错
 * @param[in]       disl2   倾向滑动位错
 * @param[in]       disl3   张裂位错
 * @param[in,out]   u       累加结果数组
 */
static void okada_add_finite_real(const OKADA_MEDIUM *c,
    real_t xi[2], real_t et[2], real_t q, int kxi[2], int ket[2],
    real_t disl1, real_t disl2, real_t disl3, real_t u[12])
{
    const real_t sd = c->sd, cd = c->cd;
    for(int k = 0; k < 2; ++k){
        for(int j = 0; j < 2; ++j){
            OKADA_FINITE f;
            real_t dua[12], du[12];
            DCCON2(c, &f, xi[j], et[k], q, kxi[k], ket[j]);
            UA(c, &f, disl1, disl2, disl3, dua);

            for(int i = 0; i < 12; i += 3){
                du[i] = -dua[i];
                du[i+1] = -dua[i+1] * cd + dua[i+2] * sd;
                du[i+2] = -dua[i+1] * sd - dua[i+2] * cd;
                if(i == 9){
                    du[i] = -du[i];
                    du[i+1] = -du[i+1];
                    du[i+2] = -du[i+2];
                }
            }

            real_t sign = (j + k == 1) ? -1.0 : 1.0;
            okada_add_scaled(u, du, sign);
        }
    }
}

/** 累加有限断层镜像源的四个角点贡献
 *
 * @param[in]       c       Okada 介质常数
 * @param[in]       xi      两个 X 方向角点坐标
 * @param[in]       et      两个倾向方向角点坐标
 * @param[in]       q       观测点到断层面的法向距离
 * @param[in]       z       观测点 Z 坐标
 * @param[in]       kxi     X 方向对数项分支标志
 * @param[in]       ket     倾向方向对数项分支标志
 * @param[in]       disl1   走向滑动位错
 * @param[in]       disl2   倾向滑动位错
 * @param[in]       disl3   张裂位错
 * @param[in,out]   u       累加结果数组
 */
static void okada_add_finite_image(const OKADA_MEDIUM *c,
    real_t xi[2], real_t et[2], real_t q, real_t z, int kxi[2], int ket[2],
    real_t disl1, real_t disl2, real_t disl3, real_t u[12])
{
    const real_t sd = c->sd, cd = c->cd;
    for(int k = 0; k < 2; ++k){
        for(int j = 0; j < 2; ++j){
            OKADA_FINITE f;
            real_t dua[12], dub[12], duc[12], du[12];
            DCCON2(c, &f, xi[j], et[k], q, kxi[k], ket[j]);
            UA(c, &f, disl1, disl2, disl3, dua);
            UB(c, &f, disl1, disl2, disl3, dub);
            UC(c, &f, z, disl1, disl2, disl3, duc);

            for(int i = 0; i < 12; i += 3){
                du[i] = dua[i] + dub[i] + z*duc[i];
                du[i+1] = (dua[i+1] + dub[i+1] + z * duc[i+1]) * cd
                    - (dua[i+2] + dub[i+2] + z * duc[i+2]) * sd;
                du[i+2] = (dua[i+1] + dub[i+1] - z * duc[i+1]) * sd
                    + (dua[i+2] + dub[i+2] - z * duc[i+2]) * cd;
                if(i == 9){
                    du[9] += duc[0];
                    du[10] += duc[1] * cd - duc[2] * sd;
                    du[11] += -duc[1] * sd - duc[2] * cd;
                }
            }

            real_t sign = (j + k == 1) ? -1.0 : 1.0;
            okada_add_scaled(u, du, sign);
        }
    }
}

/** 计算矩形有限断层的静态位移和位移偏导
 *
 * 对应 Fortran 源码中的 DC3D，公开结果不在此处进行单位换算
 *
 * @param[in]   alpha   均匀介质常数
 * @param[in]   x       局部 X 坐标
 * @param[in]   y       局部 Y 坐标
 * @param[in]   z       观测点 Z 坐标，向上为正
 * @param[in]   depth   断层参考点深度
 * @param[in]   dip     断层倾角，单位为度
 * @param[in]   al1     沿走向坐标范围下界
 * @param[in]   al2     沿走向坐标范围上界
 * @param[in]   aw1     沿上倾方向坐标范围下界
 * @param[in]   aw2     沿上倾方向坐标范围上界
 * @param[in]   disl1   走向滑动位错，正值为左旋
 * @param[in]   disl2   倾向滑动位错，正值为逆冲
 * @param[in]   disl3   张裂位错，正值为张开
 * @param[out]  u       局部 X、Y、Z 位移
 * @param[out]  upar    位移偏导 [偏导坐标][位移分量]
 * @return              0 表示正常，1 表示奇异点，2 表示输入观测点 Z 为正
 */
int grt_okada_dc3d(
    real_t alpha, real_t x, real_t y, real_t z, real_t depth, real_t dip,
    real_t al1, real_t al2, real_t aw1, real_t aw2,
    real_t disl1, real_t disl2, real_t disl3,
    real_t u[3], real_t upar[3][3])
{
    OKADA_MEDIUM c;
    real_t out[12];
    real_t xi[2], et[2];
    int kxi[2] = {0, 0};
    int ket[2] = {0, 0};
    int iret = 0;

    okada_zero_u(out);
    if(z > 0.0){
        okada_copy_output(out, u, upar);
        return 2;
    }

    DCCON0(&c, alpha, dip);
    const real_t sd = c.sd;
    const real_t cd = c.cd;
    xi[0] = x - al1;
    xi[1] = x - al2;
    if(fabs(xi[0]) < OKADA_EPS) xi[0] = 0.0;
    if(fabs(xi[1]) < OKADA_EPS) xi[1] = 0.0;

    // 计算真实断层面对应的角点坐标和分支标志
    real_t d = depth + z;
    real_t p = y * cd + d * sd;
    real_t q = y * sd - d * cd;
    et[0] = p-aw1;
    et[1] = p-aw2;
    if(fabs(q) < OKADA_EPS) q = 0.0;
    if(fabs(et[0]) < OKADA_EPS) et[0] = 0.0;
    if(fabs(et[1]) < OKADA_EPS) et[1] = 0.0;

    if(okada_finite_singular(q, xi[0], xi[1], et[0], et[1])){
        okada_copy_output(out, u, upar);
        return 1;
    }

    real_t r12 = sqrt(xi[0]*xi[0]+et[1]*et[1]+q*q);
    real_t r21 = sqrt(xi[1]*xi[1]+et[0]*et[0]+q*q);
    real_t r22 = sqrt(xi[1]*xi[1]+et[1]*et[1]+q*q);
    if(xi[0] < 0.0 && r21+xi[1] < OKADA_EPS) kxi[0] = 1;
    if(xi[0] < 0.0 && r22+xi[1] < OKADA_EPS) kxi[1] = 1;
    if(et[0] < 0.0 && r12+et[1] < OKADA_EPS) ket[0] = 1;
    if(et[0] < 0.0 && r22+et[1] < OKADA_EPS) ket[1] = 1;

    okada_add_finite_real(&c, xi, et, q, kxi, ket, disl1, disl2, disl3, out);

    // 重新计算镜像断层面对应的角点坐标和分支标志
    d = depth - z;
    p = y * cd + d * sd;
    q = y * sd - d * cd;
    et[0] = p-aw1;
    et[1] = p-aw2;
    if(fabs(q) < OKADA_EPS) q = 0.0;
    if(fabs(et[0]) < OKADA_EPS) et[0] = 0.0;
    if(fabs(et[1]) < OKADA_EPS) et[1] = 0.0;

    if(okada_finite_singular(q, xi[0], xi[1], et[0], et[1])){
        okada_copy_output(out, u, upar);
        return 1;
    }

    kxi[0] = kxi[1] = ket[0] = ket[1] = 0;
    r12 = sqrt(xi[0]*xi[0]+et[1]*et[1]+q*q);
    r21 = sqrt(xi[1]*xi[1]+et[0]*et[0]+q*q);
    r22 = sqrt(xi[1]*xi[1]+et[1]*et[1]+q*q);
    if(xi[0] < 0.0 && r21+xi[1] < OKADA_EPS) kxi[0] = 1;
    if(xi[0] < 0.0 && r22+xi[1] < OKADA_EPS) kxi[1] = 1;
    if(et[0] < 0.0 && r12+et[1] < OKADA_EPS) ket[0] = 1;
    if(et[0] < 0.0 && r22+et[1] < OKADA_EPS) ket[1] = 1;

    okada_add_finite_image(&c, xi, et, q, z, kxi, ket, disl1, disl2, disl3, out);
    okada_copy_output(out, u, upar);
    return iret;
}
