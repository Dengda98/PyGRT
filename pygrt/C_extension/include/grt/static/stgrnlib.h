/**
 * @file   stgrnlib.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * 静态格林函数库 STGRNLIB
 *
 * nc 维度为 depsrc×deprcv×north×east；
 * C 侧深度/水平坐标数组为 depsrcs、deprcvs、norths、easts；
 * 物理上格林函数还依赖震中距 r = hypot(north, east)；
 * 用 -R 建库时通常为 nnorth=1, norths=[0], easts=R（震中距序列）
 *
 */

#pragma once

#include <stdbool.h>
#include "grt/common/const.h"
#include "grt/integral/integ_process.h"

/**
 * 静态格林函数库
 *
 * 数组布局：
 *   - depsrcs[ndepsrc], deprcvs[ndeprcv] 升序 (km)
 *   - norths[nnorth], easts[neast] 为 north / east 坐标 (km)
 *   - u[is][ir][ipt][im][c]，ipt = ieast + inorth*neast，nr = nnorth*neast
 *   - 分量符号约定与 nc 文件 / Python dict 一致（Z 已取反等）
 */
typedef struct {
    size_t ndepsrc;             ///< 震源深度数量
    real_t *depsrcs;            ///< depsrcs[ndepsrc] (km)

    size_t ndeprcv;             ///< 接收深度数量
    real_t *deprcvs;            ///< deprcvs[ndeprcv] (km)

    size_t nnorth;              ///< north 方向点数
    size_t neast;               ///< east 方向点数（-R 建库时即震中距点数）
    real_t *norths;             ///< north 坐标 norths[nnorth]
    real_t *easts;              ///< east 坐标 easts[neast]（-R 建库时即震中距）

    bool calc_upar;             ///< 是否含位移空间导数

    real_t *src_va;             ///< 各震源深度处 P 波速 src_va[ndepsrc]
    real_t *src_vb;             ///< 各震源深度处 S 波速
    real_t *src_rho;            ///< 各震源深度处密度

    real_t *rcv_va;             ///< 各接收深度处 P 波速 rcv_va[ndeprcv]
    real_t *rcv_vb;             ///< 各接收深度处 S 波速
    real_t *rcv_rho;            ///< 各接收深度处密度

    realChnlGrid ***u;          ///< u[is][ir]，每个为 realChnlGrid[nr]
    realChnlGrid ***uiz;        ///< 可选，calc_upar=false 时为 NULL
    realChnlGrid ***uir;        ///< 可选，calc_upar=false 时为 NULL
} STGRNLIB;


/**
 * 按已填好的 ndepsrc/ndeprcv/nnorth/neast/calc_upar 申请 u/uiz/uir
 *
 * @param[in,out]  lib   STGRNLIB 结构体
 */
void grt_stgrnlib_allocate_u(STGRNLIB *lib);

/**
 * 按维度申请空壳 STGRNLIB：拷贝坐标轴，介质数组清零，并分配 u/uiz/uir
 *
 * 调用方随后填入介质与格林函数；最终用 grt_stgrnlib_free 释放
 *
 * @param[in]   ndepsrc    震源深度点数
 * @param[in]   depsrcs    震源深度数组 (km)
 * @param[in]   ndeprcv    接收深度点数
 * @param[in]   deprcvs    接收深度数组 (km)
 * @param[in]   nnorth     north 方向点数
 * @param[in]   norths     north 坐标数组 (km)
 * @param[in]   neast      east 方向点数
 * @param[in]   easts      east 坐标数组 (km)
 * @param[in]   calc_upar  是否分配位移偏导
 * @return      新分配的 STGRNLIB*
 */
STGRNLIB *grt_stgrnlib_alloc(
    size_t ndepsrc, const real_t *depsrcs,
    size_t ndeprcv, const real_t *deprcvs,
    size_t nnorth,  const real_t *norths,
    size_t neast,   const real_t *easts,
    bool calc_upar);

/**
 * 仅释放 u/uiz/uir
 *
 * @param[in,out]  lib   STGRNLIB 结构体
 */
void grt_stgrnlib_free_u(STGRNLIB *lib);

/**
 * 释放 STGRNLIB 内部所有堆内存，并 free(lib) 本身
 * （适用于 load_nc / create 返回的指针）
 *
 * @param[in,out]  lib   STGRNLIB 结构体，可为 NULL
 */
void grt_stgrnlib_free(STGRNLIB *lib);

/**
 * 从四维 nc 文件加载完整 STGRNLIB
 *
 * @param[in]   path   nc 文件路径
 * @return      新分配的 STGRNLIB*，调用方负责 grt_stgrnlib_free
 */
STGRNLIB *grt_stgrnlib_load_nc(const char *path);

/**
 * 将 STGRNLIB 写为四维 nc 文件
 *
 * @param[in]   lib    已填充的库
 * @param[in]   path   输出路径
 */
void grt_stgrnlib_save_nc(const STGRNLIB *lib, const char *path);

/**
 * 由内存数组创建 STGRNLIB（会拷贝一份数据）
 *
 * u / uiz / uir 布局为连续数组：
 *   [ndepsrc][ndeprcv][nr][SRC_M_NUM][CHANNEL_NUM]，nr = nnorth*neast
 * uiz/uir 在 calc_upar=false 时可传 NULL
 *
 * @param[in]   ndepsrc    震源深度点数
 * @param[in]   depsrcs    震源深度数组 (km)
 * @param[in]   ndeprcv    接收深度点数
 * @param[in]   deprcvs    接收深度数组 (km)
 * @param[in]   nnorth     north 方向点数
 * @param[in]   norths     north 坐标数组 (km)
 * @param[in]   neast      east 方向点数
 * @param[in]   easts      east 坐标数组 (km)
 * @param[in]   calc_upar  是否含位移偏导
 * @param[in]   src_va     各震源深度 P 波速
 * @param[in]   src_vb     各震源深度 S 波速
 * @param[in]   src_rho    各震源深度密度
 * @param[in]   rcv_va     各接收深度 P 波速
 * @param[in]   rcv_vb     各接收深度 S 波速
 * @param[in]   rcv_rho    各接收深度密度
 * @param[in]   u          位移格林函数连续数组
 * @param[in]   uiz        位移对 z 偏导连续数组，可为 NULL
 * @param[in]   uir        位移对 r 偏导连续数组，可为 NULL
 * @return      新分配的 STGRNLIB*，调用方负责 grt_stgrnlib_free
 */
STGRNLIB *grt_stgrnlib_create(
    size_t ndepsrc, const real_t *depsrcs,
    size_t ndeprcv, const real_t *deprcvs,
    size_t nnorth,  const real_t *norths,
    size_t neast,   const real_t *easts,
    bool calc_upar,
    const real_t *src_va,  const real_t *src_vb,  const real_t *src_rho,
    const real_t *rcv_va,  const real_t *rcv_vb,  const real_t *rcv_rho,
    const real_t *u,       const real_t *uiz,     const real_t *uir);

/**
 * 由 STGRNLIB 推算默认子断层尺寸 \f$\min(dr, dz)\f$
 *
 * - dr：震中距采样最小间距，r = hypot(north, east)
 * - dz：震源深度最小间隔（depsrcs）
 * 某一维仅 1 点时只用另一维；两者皆不可用则报错
 *
 * @param[in]   lib    STGRNLIB 结构体
 * @return      默认子断层边长 (km)
 */
real_t grt_stgrnlib_default_subfault_size(const STGRNLIB *lib);

/**
 * 循环计算多震源/接收深度静态格林函数并写入单个四维 nc
 *
 * Kproc 须已由 grt_prepare_static_grn 填好深度无关字段；
 * 循环内按各 (depsrc, deprcv) 的 hs 更新局部拷贝的 k0，避免积分过程改写模板
 *
 * @param[in]   modelpath     一维模型文件
 * @param[in]   ndepsrc       震源深度点数
 * @param[in]   depsrcs       震源深度数组 (km)
 * @param[in]   ndeprcv       接收深度点数
 * @param[in]   deprcvs       接收深度数组 (km)
 * @param[in]   nnorth        north 方向点数
 * @param[in]   norths        north 坐标数组 (km)
 * @param[in]   neast         east 方向点数
 * @param[in]   easts         east 坐标数组 (km)
 * @param[in]   k0            用户侧波数系数（同 greenfn -K+k）；循环内按 hs 缩放
 * @param[in]   Kproc         深度无关的波数积分参数模板
 * @param[in]   topbound      顶界面边界条件
 * @param[in]   botbound      底界面边界条件
 * @param[in]   calc_upar     是否计算位移偏导
 * @param[in]   outpath       输出 nc 路径
 * @param[in]   statsstr      波数积分中间文件目录；NULL 表示不输出；
 *                            多深度时应传 NULL
 */
void grt_compute_stgrnlib_to_nc(
    const char *modelpath,
    size_t ndepsrc, const real_t *depsrcs,
    size_t ndeprcv, const real_t *deprcvs,
    size_t nnorth,  const real_t *norths,
    size_t neast,   const real_t *easts,
    real_t k0,
    K_INTEG_PROCESS *Kproc,
    int topbound, int botbound,
    bool calc_upar,
    const char *outpath,
    const char *statsstr);
