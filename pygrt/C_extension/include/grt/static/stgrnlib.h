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

/**
 * 静态格林函数库
 *
 * 数组布局：
 *   - depsrcs[ndepsrc], deprcvs[ndeprcv] 严格升序 (km)
 *   - norths[nnorth], easts[neast] 为 north / east 坐标 (km)，各自严格升序
 *   - rs[nr] 为网格序震中距，ipt = ieast + inorth*neast，nr = nnorth*neast
 *   - sort_rs / sort_rs_idx 为 rs 的升序排列及回指网格 ipt 的索引
 *   - isUniform / dr：网格序 rs 是否已是等距升序及其步长（同 syn 查找加速）
 *   - u[is][ir][ipt][im][c]
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

    size_t nr;                  ///< = nnorth * neast
    real_t *rs;                 ///< 网格序震中距 rs[nr]
    real_t *sort_rs;            ///< 升序震中距 sort_rs[nr]
    size_t *sort_rs_idx;        ///< sort_rs[i] 对应网格 ipt = sort_rs_idx[i]
    bool isUniform;             ///< 网格序 rs 是否已是等距升序
    real_t dr;                  ///< 等距步长；非等距时无意义

    bool calc_upar;             ///< 是否含位移空间导数

    real_t *src_va;             ///< 各震源深度处 P 波速 src_va[ndepsrc]
    real_t *src_vb;             ///< 各震源深度处 S 波速
    real_t *src_rho;            ///< 各震源深度处密度

    real_t *rcv_va;             ///< 各接收深度处 P 波速 rcv_va[ndeprcv]
    real_t *rcv_vb;             ///< 各接收深度处 S 波速
    real_t *rcv_rho;            ///< 各接收深度处密度

    size_t nlayer;              ///< 建库模型层数
    real_t (*modarr)[GRT_MODARR_NCOL]; ///< 模型矩阵 [nlayer][6]：Thk/Va/Vb/Rho/Qa/Qb

    realChnlGrid ***u;          ///< u[is][ir]，每个为 realChnlGrid[nr]
    realChnlGrid ***uiz;        ///< 可选，calc_upar=false 时为 NULL
    realChnlGrid ***uir;        ///< 可选，calc_upar=false 时为 NULL
} STGRNLIB;


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
 * 释放 STGRNLIB 内部所有堆内存，并 free(lib) 本身
 * （适用于 load_nc / alloc 返回的指针）
 *
 * @param[in,out]  lib   STGRNLIB 结构体，可为 NULL
 */
void grt_stgrnlib_free(STGRNLIB *lib);

/**
 * 设置建库模型矩阵（会拷贝一份）
 *
 * @param[in,out]  lib      STGRNLIB
 * @param[in]      nlayer   层数，须 > 0
 * @param[in]      modarr   模型矩阵，每行 Thk/Va/Vb/Rho/Qa/Qb
 */
void grt_stgrnlib_set_modarr(
    STGRNLIB *lib, size_t nlayer, const real_t (*modarr)[GRT_MODARR_NCOL]);

/**
 * 由震中距与震源深度采样推断默认子断层尺寸 min(dr, dz)
 *
 * @param[in]   lib   STGRNLIB
 * @return      默认 dL=dW (km)
 */
real_t grt_stgrnlib_default_subfault_size(const STGRNLIB *lib);

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
