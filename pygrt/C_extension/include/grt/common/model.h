/**
 * @file   model.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * `GRT_MODEL1D` 结构体相关操作函数
 */

#pragma once

#include <complex.h>
#include <stdbool.h>
#include "grt/common/const.h"
#include "grt/common/RT_matrix.h"


/** 1D 模型结构体，包括多个水平层，以及复数形式的弹性参数 */
typedef struct {
    size_t n;  ///< 层数，注意包括了震源和接收点的虚拟层，(n>=3)
    real_t depsrc; ///< 震源深度 km
    real_t deprcv; ///< 接收点深度 km
    size_t isrc; ///< 震源所在虚拟层位, isrc>=1
    size_t ircv; ///< 接收点所在虚拟层位, ircv>=1, ircv != isrc
    bool ircvup; ///< 接收点位于浅层, ircv < isrc
    bool io_depth; ///< 读取的模型首列为每层顶界面深度

    cplx_t omega;   ///< 圆频率
    real_t k;   ///< 波数
    cplx_t c_phase;   ///< 当前相速度

    real_t *Thk; ///< Thk[n], 最后一层厚度不使用(当作正无穷), km
    real_t *Dep; ///< Dep[n], 每一层顶界面深度，第一层必须为 0.0
    real_t *Va;  ///< Va[n]   P波速度  km/s
    real_t *Vb;  ///< Vb[n]   S波速度  km/s
    real_t *Rho; ///< Rho[n]  密度  g/cm^3
    real_t *Qa; ///< Qa[n]     P波Q值
    real_t *Qb; ///< Qb[n]     S波Q值
    real_t *Qainv; ///<   1/Q_p
    real_t *Qbinv; ///<   1/Q_s

    cplx_t *mu;       ///< mu[n] \f$ V_b^2 * \rho \f$
    cplx_t *lambda;   ///< lambda[n] \f$ V_a^2 * \rho - 2*\mu \f$
    cplx_t *delta;    ///< delta[n] \f$ (\lambda+\mu)/(\lambda+3*\mu) \f$
    cplx_t *atna;
    cplx_t *atnb;
    cplx_t *xa;
    cplx_t *xb;
    cplx_t *caca;
    cplx_t *cbcb;

    /* 状态变量，非 0 为异常值 */
    int stats;

    /* 根据震源和台站划分的广义 R/T 系数 */
    RT_MATRIX M_AL;
    RT_MATRIX M_BL;
    RT_MATRIX M_RS;
    RT_MATRIX M_FA;
    RT_MATRIX M_FB;

    /* 自由表面的反射系数矩阵，仅 RU, RUL 有用 */
    RT_MATRIX M_top;

    /* 接收点处的接收矩阵 (转为位移u和位移导数uiz的(B_m, C_m, P_m)系分量) */
    cplx_t R_EV[2][2];
    cplx_t R_EVL;
    cplx_t uiz_R_EV[2][2];
    cplx_t uiz_R_EVL;

    /* 震源处的震源系数 \f$ P_m, SV_m, SH_m  */
    cplxChnlGrid src_coefD;
    cplxChnlGrid src_coefU;

} GRT_MODEL1D;


/**
 * 打印 GRT_MODEL1D 模型参数信息，主要用于调试程序 
 * 
 * @param[in]    mod1d    `GRT_MODEL1D` 结构体指针
 * 
 */
void grt_print_mod1d(const GRT_MODEL1D *mod1d);

/**
 * 释放 `GRT_MODEL1D` 结构体指针 
 * 
 * @param[out]     mod1d      `GRT_MODEL1D` 结构体指针
 */
void grt_free_mod1d(GRT_MODEL1D *mod1d);

/**
 * 初始化 GRT_MODEL1D 模型内存空间 
 * 
 * @param[in]    n        模型层数 
 * 
 * @return    `GRT_MODEL1D` 结构体指针
 * 
 */
GRT_MODEL1D * grt_init_mod1d(size_t n);

/**
 * 复制 `GRT_MODEL1D` 结构体
 * 
 * @param[in]     mod1d1    `GRT_MODEL1D` 源结构体指针
 * @return        复制好的 `GRT_MODEL1D` 结构体指针
 * 
 */
GRT_MODEL1D * grt_copy_mod1d(const GRT_MODEL1D *mod1d1);

/**
 * 根据不同的 omega， 计算衰减系数，更新弹性模量
 * 
 * @param[in,out]     mod1d     `MODEL1D` 结构体指针
 * @param[in]         omega     复数频率
 */
void grt_attenuate_mod1d(GRT_MODEL1D *mod1d, cplx_t omega);

/**
 * 根据记录好的圆频率和波数，计算相速度和每层的 xa, xb, caca, cbcb
 * 
 * @param[in,out]      mod1d    模型结构体指针
 * @param[in]          k        波数
 */
void grt_mod1d_xa_xb(GRT_MODEL1D *mod1d, const real_t k);


/**
 * 扩容 `GRT_MODEL1D` 结构体
 * 
 * @param[in,out]     mod1d     `MODEL1D` 结构体指针
 * @param[in]         n         新层数
 */
void grt_realloc_mod1d(GRT_MODEL1D *mod1d, size_t n);

/**
 * 从文件中读取模型文件
 * 
 * @param[in]    modelpath      模型文件路径
 * @param[in]    depsrc         震源深度
 * @param[in]    deprcv         接收深度
 * @param[in]    allowLiquid    是否允许液体层
 * 
 * @return    `GRT_MODEL1D` 结构体指针
 * 
 */
GRT_MODEL1D * grt_read_mod1d_from_file(const char *modelpath, real_t depsrc, real_t deprcv, bool allowLiquid);


/**
 * 从模型文件中判断各个量的大致精度（字符串长度），以确定浮点数输出位数
 * 
 * @param[in]    modelpath      模型文件路径
 * @param[out]   diglen         每一列的最大字符串长度
 * 
 */
void grt_get_model_diglen_from_file(const char *modelpath, size_t diglen[6]);

/**
 * 浮点数比较，检查模型中是否存在该速度（不论Vp,Vs）
 * 
 * @param[in]   mod1d    模型
 * @param[in]   vel      输入速度
 * @param[in]   tol      浮点数比较精度
 * 
 * @return    是否存在
 */
bool grt_check_vel_in_mod(const GRT_MODEL1D *mod1d, const real_t vel, const real_t tol);

/**
 * 计算最大最小速度（非零值）
 * 
 * @param    mod1d   (in)`GRT_MODEL1D` 结构体指针
 * @param    vmin    (out)最小速度
 * @param    vmax    (out)最大速度
 * 
 */
void grt_get_mod1d_vmin_vmax(const GRT_MODEL1D *mod1d, real_t *vmin, real_t *vmax);