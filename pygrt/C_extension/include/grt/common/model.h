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


/** 1D 模型结构体，包括多个水平层，以及复数形式的弹性参数 */
typedef struct {
    MYINT n;  ///< 层数，注意包括了震源和接收点的虚拟层，(n>=3)
    MYREAL depsrc; ///< 震源深度 km
    MYREAL deprcv; ///< 接收点深度 km
    MYINT isrc; ///< 震源所在虚拟层位, isrc>=1
    MYINT ircv; ///< 接收点所在虚拟层位, ircv>=1, ircv != isrc
    bool ircvup; ///< 接收点位于浅层, ircv < isrc
    bool io_depth; ///< 读取的模型首列为每层顶界面深度

    MYREAL *Thk; ///< Thk[n], 最后一层厚度不使用(当作正无穷), km
    MYREAL *Dep; ///< Dep[n], 每一层顶界面深度，第一层必须为 0.0
    MYREAL *Va;  ///< Va[n]   P波速度  km/s
    MYREAL *Vb;  ///< Vb[n]   S波速度  km/s
    MYREAL *Rho; ///< Rho[n]  密度  g/cm^3
    MYREAL *Qa; ///< Qa[n]     P波Q值
    MYREAL *Qb; ///< Qb[n]     S波Q值
    MYREAL *Qainv; ///<   1/Q_p
    MYREAL *Qbinv; ///<   1/Q_s

    MYCOMPLEX *mu;       ///< mu[n] \f$ V_b^2 * \rho \f$
    MYCOMPLEX *lambda;   ///< lambda[n] \f$ V_a^2 * \rho - 2*\mu \f$
    MYCOMPLEX *delta;    ///< delta[n] \f$ (\lambda+\mu)/(\lambda+3*\mu) \f$
    MYCOMPLEX *atna;
    MYCOMPLEX *atnb;

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
GRT_MODEL1D * grt_init_mod1d(MYINT n);

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
void grt_attenuate_mod1d(GRT_MODEL1D *mod1d, MYCOMPLEX omega);

/**
 * 根据相速度和层位，计算 iy 层的 (c/vp)^2, (c/vs)^2 以及归一化垂直波数
 * 
 * @param[in]      mod1d1       `GRT_MODEL1D` 结构体指针
 * @param[in]      iy           层位索引
 * @param[in]      c_phase      相速度
 * @param[out]     pt_caca      \f$ (\frac{c}{V_P})^2 \f$
 * @param[out]     pt_xa        \f$ \sqrt{1 - (\frac{c}{V_P})^2} \f$
 * @param[out]     pt_cbcb      \f$ (\frac{c}{V_S})^2 \f$
 * @param[out]     pt_xb        \f$ \sqrt{1 - (\frac{c}{V_S})^2} \f$
 */
void grt_get_mod1d_xa_xb(
    const GRT_MODEL1D *mod1d, const MYINT iy, const MYCOMPLEX c_phase, 
    MYCOMPLEX *pt_caca, MYCOMPLEX *pt_xa, MYCOMPLEX *pt_cbcb, MYCOMPLEX *pt_xb);


/**
 * 扩容 `GRT_MODEL1D` 结构体
 * 
 * @param[in,out]     mod1d     `MODEL1D` 结构体指针
 * @param[in]         n         新层数
 */
void grt_realloc_mod1d(GRT_MODEL1D *mod1d, MYINT n);

/**
 * 从文件中读取模型文件
 * 
 * @param[in]    command        命令名称
 * @param[in]    modelpath      模型文件路径
 * @param[in]    depsrc         震源深度
 * @param[in]    deprcv         接收深度
 * @param[in]    allowLiquid    是否允许液体层
 * 
 * @return    `GRT_MODEL1D` 结构体指针
 * 
 */
GRT_MODEL1D * grt_read_mod1d_from_file(const char *command, const char *modelpath, double depsrc, double deprcv, bool allowLiquid);


/**
 * 从模型文件中判断各个量的大致精度（字符串长度），以确定浮点数输出位数
 * 
 * @param[in]    command        命令名称
 * @param[in]    modelpath      模型文件路径
 * @param[out]   diglen         每一列的最大字符串长度
 * 
 */
void grt_get_model_diglen_from_file(const char *command, const char *modelpath, MYINT diglen[6]);

/**
 * 浮点数比较，检查模型中是否存在该速度（不论Vp,Vs）
 * 
 * @param[in]   mod1d    模型
 * @param[in]   vel      输入速度
 * @param[in]   tol      浮点数比较精度
 * 
 * @return    是否存在
 */
bool grt_check_vel_in_mod(const GRT_MODEL1D *mod1d, const MYREAL vel, const MYREAL tol);

/**
 * 计算最大最小速度（非零值）
 * 
 * @param    mod1d   (in)`GRT_MODEL1D` 结构体指针
 * @param    vmin    (out)最小速度
 * @param    vmax    (out)最大速度
 * 
 */
void grt_get_mod1d_vmin_vmax(const GRT_MODEL1D *mod1d, double *vmin, double *vmax);