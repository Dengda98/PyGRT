/**
 * @file   lamb_basic.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-09
 *
 *    第二类和第三类 Lamb 问题共用的基本积分
 */

#pragma once

#include "grt/common/const.h"

/** P、S 和 S-P 基本积分对应的积分项 */
typedef enum {
    LAMB_BASIC_P_TERM,
    LAMB_BASIC_S_TERM,
    LAMB_BASIC_SP_TERM,
} LAMB_BASIC_TERM;

/** 两类 Lamb 求解器共用的材料和角度参数 */
typedef struct {
    real_t k;   ///< 波速比 k=beta/alpha
    real_t k2;  ///< k^2
    real_t kp2; ///< kp^2=1-k^2
    real_t st;  ///< 当前积分路径射线角的正弦 sin(theta)
    real_t ct;  ///< 当前积分路径射线角的余弦 cos(theta)
} LAMB_BASIC_VARS;

/** V 基本积分两个二次因子分解所需的辅助量 */
typedef struct {
    cplx_t h1p;   ///< h1p=c+xi1^2
    cplx_t h2p;   ///< h2p=c+xi2^2
    cplx_t h1m;   ///< h1m=c-xi1^2
    cplx_t h2m;   ///< h2m=c-xi2^2
    cplx_t D;     ///< D=(xi1-xi2)^2-2*(c+xi1*xi2)
    cplx_t s[2];  ///< s[0/1] 对应 root=sqrt(-c) 的正号/负号分支
    cplx_t z0[2]; ///< z0[0/1] 对应 root=sqrt(-c) 的正号/负号分支
} LAMB_BASIC_V_AUX;

/** 某一个时间点和某一类积分所需的变量 */
typedef struct {
    LAMB_BASIC_TERM term; ///< 当前基本积分类型，取 P、S 或 S-P
    real_t tbar;          ///< 当前无量纲时间
    real_t kp2;           ///< kp^2=1-k^2
    real_t m;             ///< 第一分式线性变换参数 m=tbar*cos(theta)
    real_t n;             ///< 第二分式线性变换参数，随 P、S 或 S-P 项而变化
    real_t xi1;           ///< 二次因子分解中的较大根 xi1
    real_t xi2;           ///< 二次因子分解中的较小根 xi2
    real_t z1;            ///< z1=sqrt(z1sq)
    real_t z2;            ///< z2=sqrt(z2sq)
    real_t z1sq;          ///< z1^2
    real_t z2sq;          ///< z2^2
    real_t m_elliptic;    ///< 椭圆积分参数
    real_t c_main;        ///< 主积分项的公共系数
    real_t c_residue;     ///< 留数项的公共系数，S-P 项中为零
    real_t c1;            ///< S-P 项的辅助量 c1=z2^2-1
    real_t c2;            ///< S-P 项的辅助量 c2=z2^2-z1^2
} LAMB_BASIC_CONTEXT;

/** 构造 P 波反射项在一个无量纲时间点上的基本积分上下文 */
void grt_lamb_make_context_P(const real_t tbar, const real_t tbar2, const LAMB_BASIC_VARS *V, LAMB_BASIC_CONTEXT *ctx);

/** 构造 S 波反射项在一个无量纲时间点上的基本积分上下文 */
void grt_lamb_make_context_S(const real_t tbar, const real_t tbar2, const LAMB_BASIC_VARS *V, LAMB_BASIC_CONTEXT *ctx);

/** 构造 S-P 转换项在一个无量纲时间点上的基本积分上下文 */
void grt_lamb_make_context_SP(const real_t tbar, const real_t tbar2, const LAMB_BASIC_VARS *V, LAMB_BASIC_CONTEXT *ctx);

/** 计算指定阶数的 U 基本积分 */
cplx_t grt_lamb_basic_U(const int number, const cplx_t c, const LAMB_BASIC_CONTEXT *ctx);

/** 同时计算同一分母对应的两个 U 基本积分 */
void grt_lamb_make_U_pair(const cplx_t c, const LAMB_BASIC_CONTEXT *ctx, cplx_t result[2]);

/** 根据椭圆积分参数递推计算 H0 至 H4 基本积分组合 */
void grt_lamb_calculate_H(const LAMB_BASIC_CONTEXT *ctx, const real_t K, real_t H[5]);

/** 计算 P 波反射项的 V 基本积分尾项 */
real_t grt_lamb_tail_V_P(const int number, const LAMB_BASIC_CONTEXT *ctx, const real_t H[5]);

/** 计算 S 波反射项的 V 基本积分尾项 */
real_t grt_lamb_tail_V_S(const int number, const LAMB_BASIC_CONTEXT *ctx, const real_t H[5]);

/** 计算 S-P 转换项的 V 基本积分尾项 */
real_t grt_lamb_tail_V_SP(const int number, const LAMB_BASIC_CONTEXT *ctx, const real_t H[5]);

/** 计算 V 基本积分所需的 h、s 和 z0 辅助量 */
void grt_lamb_make_V_aux(const cplx_t c, const LAMB_BASIC_CONTEXT *ctx, LAMB_BASIC_V_AUX *aux);

/** 计算 P 波反射项对应的两个 V 基本积分 */
void grt_lamb_make_V_P_pair(const cplx_t c, const LAMB_BASIC_CONTEXT *ctx, const real_t K, cplx_t result[2]);

/** 计算 S 波反射项对应的两个 V 基本积分 */
void grt_lamb_make_V_S_pair(const cplx_t c, const LAMB_BASIC_CONTEXT *ctx, const real_t K, cplx_t result[2]);

/** 计算 S-P 转换项对应的两个 V 基本积分 */
void grt_lamb_make_V_SP_pair(const cplx_t c, const LAMB_BASIC_CONTEXT *ctx, const real_t K, cplx_t result[2]);
