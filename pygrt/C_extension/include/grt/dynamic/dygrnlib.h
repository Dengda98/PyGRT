/**
 * @file   dygrnlib.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-10-03
 *
 *    动态格林函数库 DYGRNLIB：SAC 库扫描、读取和查询
 */

#pragma once

#include "grt/common/sacio.h"
#include "grt/common/search.h"
#include "grt/common/finite_fault.h"

/** 单个库节点，保存路径、几何、波形起点和震源介质 */
typedef struct {
    char *path;  ///< SAC 子目录
    real_t zs;   ///< 震源深度，km
    real_t zr;   ///< 接收深度，km
    real_t r;    ///< 震中距，km
    real_t b;    ///< 波形起点，s
    real_t vp;   ///< 震源点 P 波速度，km/s
    real_t vs;   ///< 震源点 S 波速度，km/s
    real_t rho;  ///< 震源点密度，g/cm^3
} DYGRNLIB_NODE;

/**
 * 单次角点合成使用的波形，由调用方独立拥有，使用后立即释放
 * gf[类型][源型][分量][样本]，类型 0/1/2 对应位移、深度偏导、距离偏导
 * 仅加载当前子源需要的基本源型，未加载的分量为 NULL，符号与库内 SAC 文件一致
 */
typedef struct {
    real_t *gf[3][GRT_SRC_M_NUM][3];  ///< 当前角点所需的格林函数波形
} DYGRNLIB_WAVEFORMS;

/**
 * 动态格林函数库，管理节点、公共时间采样、采样坐标轴及建库模型
 * 三个采样轴的所有组合必须完整，各轴可非等距
 * 节点按源深、接收深度、震中距升序排列，索引为 (is * ndeprcv + iz) * nr + ir
 * SAC 波形按需读取，不在库中缓存
 * 合成时直接使用库内时间采样，输出设置由调用方管理
 */
typedef struct {
    int nt;                             ///< 库内统一的时间样本数
    real_t dt;                          ///< 库内统一的采样间隔，s
    real_t wI;                          ///< 库内统一的虚频率，rad/s
    real_t stf_decay;                   ///< 每个采样点的阻尼因子 exp(-wI*dt)

    SACHEAD refhead;                    ///< 首个扫描节点的 SAC 参考头段
    DYGRNLIB_NODE *nodes;               ///< 完整三维网格的节点，震中距索引变化最快
    size_t nnode;                       ///< 节点数量，等于 ndepsrc * ndeprcv * nr

    real_t *depsrcs;                    ///< 升序震源深度轴，km
    real_t *deprcvs;                    ///< 升序接收深度轴，km
    real_t *rs;                         ///< 升序震中距轴，km
    size_t ndepsrc;                     ///< 震源深度采样数
    size_t ndeprcv;                     ///< 接收深度采样数
    size_t nr;                          ///< 震中距采样数

    size_t nlayer;                      ///< 模型层数
    real_t (*modarr)[GRT_MODARR_NCOL];  ///< 原始建库模型，不插虚拟层

    bool calc_upar;                     ///< 是否加载偏导
    bool direct_node_input;             ///< 是否直接以单个节点目录作为输入，与库的节点数量无关
} DYGRNLIB;

/**
 * 申请动态格林函数库，时间采样、虚频率、参考头段、坐标、节点和模型由调用方填入
 * @param[in]  calc_upar  是否包含位移空间导数
 * @return     新分配的动态库，调用方负责 grt_dygrnlib_free
 */
DYGRNLIB *grt_dygrnlib_alloc(bool calc_upar);

/**
 * 扫描 SAC 库的节点头段、采样轴和原始模型，波形按需加载
 * 输入路径自身包含基本格林函数时标记为单节点输入，否则扫描库根目录，单节点根目录仍为根目录输入
 * 每个节点仅读取一个 SAC 头段，并核对库子目录名中的源深、接收深度和震中距
 * 库根目录下所有子目录必须使用同一个模型名，且根目录中必须包含同名模型文件
 * 所有节点的 nt/dt 和虚频率必须一致，每个源深与接收深度组合必须包含相同的震中距采样
 * 时间采样或虚频率不一致、缺失组合或重复节点时报错
 * @param[in]  root       库根目录或单个格林函数节点目录
 * @param[in]  calc_upar  是否加载位移空间导数
 * @return     新分配的动态库，调用方负责 grt_dygrnlib_free
 */
DYGRNLIB *grt_dygrnlib_load(const char *root, bool calc_upar);

/**
 * 释放库拥有的节点、坐标轴和模型
 * @param[in,out]  lib  动态库，可为 NULL
 */
void grt_dygrnlib_free(DYGRNLIB *lib);

/**
 * 从三个采样轴取最小正间隔作为默认子断层尺寸
 * @param[in]  lib  动态格林函数库
 * @return     默认 dL=dW，km
 */
real_t grt_dygrnlib_default_subfault_size(const DYGRNLIB *lib);

/**
 * 查询源深、接收深度和震中距对应的角点及权重
 * @param[in]   lib      动态格林函数库
 * @param[in]   zs       震源深度，km
 * @param[in]   zr       接收深度，km
 * @param[in]   r        震中距，km
 * @param[in]   mode     精确、最近邻或线性插值查询
 * @param[out]  indices  非零权重角点在 nodes 数组中的索引，最多八个
 * @param[out]  weights  角点权重
 * @return         角点数量，超出范围或精确查询未匹配时返回 -1
 */
int grt_dygrnlib_find_corners(
    const DYGRNLIB *lib, real_t zs, real_t zr, real_t r, GRT_SAMPLE_MODE mode,
    size_t indices[8], real_t weights[8]);

/**
 * 按当前子源的基本源型读取一个角点的波形，不修改库或节点
 * 每次调用独立分配波形，可供各线程并行读取同一节点
 * @param[in]  lib     动态库，决定是否读取空间导数
 * @param[in]  index   节点在 nodes 数组中的索引
 * @param[in]  fault   当前震源断层，包含所需基本源型
 * @return     独立拥有的波形，调用方负责 grt_dygrnlib_free_waveforms
 */
DYGRNLIB_WAVEFORMS *grt_dygrnlib_read_node(const DYGRNLIB *lib, size_t index, const FINITE_FAULT *fault);

/**
 * 释放单次读取的全部波形，与库和节点的生命周期无关
 * @param[in,out]  waveforms  当前角点波形，可为 NULL
 */
void grt_dygrnlib_free_waveforms(DYGRNLIB_WAVEFORMS *waveforms);
