/**
 * @file   model.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * `MODEL1D` 结构体 以及 `PYMODEL1D` 结构体 的相关操作函数                  
 */

#pragma once

#include <complex.h>
#include <stdbool.h>
#include "const.h"

/**
 * 单个水平层的结构体
 */
typedef struct _LAYER {
    MYREAL thk;  ///< 层厚，最后一层厚度不使用(当作正无穷), km
    MYREAL Va;  ///<   P波速度  km/s
    MYREAL Vb;  ///<   S波速度  km/s
    MYREAL Rho; ///<   密度  g/cm^3
    MYREAL Qainv; ///<   1/Q_p
    MYREAL Qbinv; ///<   1/Q_s
    MYCOMPLEX mu;   ///< \f$ V_b^2 * \rho \f$
    MYCOMPLEX kaka;  ///<  \f$ (\omega/V_a)^2 \f$
    MYCOMPLEX kbkb;  ///<  \f$ (\omega/V_b)^2 \f$
} LAYER;

/**
 * 1D模型结构体，包括多个水平层
 */
typedef struct _MODEL1D {
    LAYER * lays; ///< 不同层的物性参数
    MYINT n;  ///< 层数，注意包括了震源和接收点的虚拟层，(n>=3)
    MYINT isrc; ///< 震源所在虚拟层位, isrc>=1
    MYINT ircv; ///< 接收点所在虚拟层位, ircv>=1, ircv != isrc
    MYINT imin; ///< 较浅的虚拟层位，min(isrc, ircv)
    MYINT imax; ///< 较深的虚拟层位，max(isrc, ircv)
    bool ircvup; ///< 接收点是否浅于震源层, ircv < isrc
    
} MODEL1D;


/**
 * 用于从Python传递数据的中间结构体
 */
typedef struct _PYMODEL1D {
    MYINT n;  ///< 层数，注意包括了震源和接收点的虚拟层，(n>=3)
    MYREAL depsrc; ///< 震源深度 km
    MYREAL deprcv; ///< 接收点深度 km
    MYINT isrc; ///< 震源所在虚拟层位, isrc>=1
    MYINT ircv; ///< 接收点所在虚拟层位, ircv>=1, ircv != isrc
    bool ircvup; ///< 接收点位于浅层, ircv < isrc

    MYREAL *Thk; ///< Thk[n], 最后一层厚度不使用(当作正无穷), km
    MYREAL *Va;  ///< Va[n]   P波速度  km/s
    MYREAL *Vb;  ///< Vb[n]   S波速度  km/s
    MYREAL *Rho; ///< Rho[n]  密度  g/cm^3
    MYREAL *Qa; ///< Qa[n]     P波Q值
    MYREAL *Qb; ///< Qb[n]     S波Q值
} PYMODEL1D;

/**
 * 格林函数的频谱，仅作为将结果返回至Python使用 
 */
typedef struct _GRN {
    MYINT nf; ///< 频点个数
    MYREAL *Re; ///< Re[nf] 实部
    MYREAL *Im; ///< Im[nf] 虚部
    // float *Re;
    // float *Im;
} GRN;



/**
 * 打印模型参数信息，主要用于调试程序 
 * 
 * @param    mod1d  (in)`MODEL1D` 结构体指针
 * 
 */
void print_mod1d(const MODEL1D *mod1d);

/**
 * 打印PYMODEL1D模型参数信息，主要用于调试程序 
 * 
 * @param    pymod  (in)`PYMODEL1D` 结构体指针
 * 
 */
void print_pymod(const PYMODEL1D *pymod);

/**
 * 释放 `PYMODEL1D` 结构体指针 
 * 
 * @param     pymod   (out)`PYMODEL1D` 结构体指针
 */
void free_pymod(PYMODEL1D *pymod);

/**
 * 初始化模型内存空间 
 * 
 * @param    n       (in)模型层数 
 * 
 * @return    `MODEL1D` 结构体指针
 * 
 */
MODEL1D * init_mod1d(MYINT n);

/**
 * 将 `PYMODEL1D` 结构体信息转到 `MODEL1D` 结构体中
 * 
 * @param    pymod1d    (in)`PYMODEL1D` 结构体指针
 * @param    mod1d      (out) `MODEL1D` 结构体指针
 * 
 */
void get_mod1d(const PYMODEL1D *pymod1d, MODEL1D *mod1d);

/**
 * 复制 `MODEL1D` 结构体，要求都以初始化
 * 
 * @param    mod1d1    (in)`MODEL1D` 源结构体指针
 * @param    mod1d2    (out)`MODEL1D` 目标结构体指针
 * 
 */
void copy_mod1d(const MODEL1D *mod1d1, MODEL1D *mod1d2);

/**
 * 释放 `MODEL1D` 结构体指针 
 * 
 * @param     mod1d   (out)`MODEL1D` 结构体指针
 */
void free_mod1d(MODEL1D *mod1d);


/**
 * 根据不同的omega， 更新模型的 \f$ k_a^2, k_b^2 \f$ 参数
 * 
 * @param     mod1d     (inout)`MODEL1D` 结构体指针
 * @param     omega     (in)复数频率
 */
void update_mod1d_omega(MODEL1D *mod1d, MYCOMPLEX omega);

/**
 * 初始化PYMODEL1D模型内存空间 
 * 
 * @param    n       (in)模型层数 
 * 
 * @return    `PYMODEL1D` 结构体指针
 * 
 */
PYMODEL1D * init_pymod(MYINT n);

/**
 * 从文件中读取模型文件，以PYMODEL1D结构体形式
 * 
 * @param    command       (in)命令名称
 * @param    modelpath     (in)模型文件路径
 * @param    depsrc        (in)震源深度
 * @param    deprcv        (in)接收深度
 * 
 * @return    `PYMODEL1D` 结构体指针
 * 
 */
PYMODEL1D * read_pymod_from_file(const char *command, const char *modelpath, double depsrc, double deprcv);