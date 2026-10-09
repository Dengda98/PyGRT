/**
 * @file   syn_output.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-10-04
 *
 *    动态合成结果的 SAC 头段和文件保存
 */

#pragma once

#include "grt/common/sacio.h"
#include "grt/common/rcv_point.h"

/** 动态合成的保存设置，路径和接收点均借用调用方数据 */
typedef struct {
    const char *root;               ///< 合成结果根目录
    bool rcv_subdirs;               ///< 是否按接收点建立子目录，-Q/-U 使用此方式
    const RCV_POINT *rcvs;           ///< 接收点集合
    const FINITE_FAULT *rcv_faults;  ///< 接收断层数组首地址，普通接收点为 NULL
    const char *channels;           ///< 分量名，ZRT 或 ZNE
    bool calc_upar;                 ///< 是否记录位移偏导
} DY_SYN_OUTPUT;

/**
 * 取格林函数点数和完整震源过程点数的最大值，所有接收点共用此长度
 * 震源过程点数包含采样对齐后的破裂延迟
 * @param[in]  grn_nt  格林函数点数，Lamb 模块使用指定的最小输出点数
 * @param[in]  dt      采样间隔，s
 * @param[in]  nfault  震源断层数量
 * @param[in]  faults  已初始化时间函数的震源断层
 * @return     统一输出点数
 */
int grt_syn_output_npts(int grn_nt, real_t dt, size_t nfault, const FINITE_FAULT *faults);

/**
 * 返回按名称前缀的数字索引排序的直属子目录；没有子目录时返回传入目录本身
 * 不读取 SAC 文件，调用方释放每个路径及数组
 * @param[in]   path   合成目录
 * @param[out]  count  接收点数量
 */
char **grt_syn_output_receiver_directories(const char *path, size_t *count);

/**
 * 填写源台几何、接收断层索引及介质头段，保留采样和震相信息
 * 有限源的震源深度传入 SAC_FLOAT_UNDEF，源介质保持未定义值
 * 无原始模型时保留调用方提供的介质头段
 * @param[in]      output  合成输出设置
 * @param[in]      ipt     接收点索引
 * @param[in,out]  hd      已初始化的 SAC 头段
 * @param[in]      depsrc  点源深度，km，有限源为 SAC_FLOAT_UNDEF
 * @param[in]      nlayer  原始模型层数
 * @param[in]      modarr  原始模型矩阵，无模型时为 NULL
 */
void grt_syn_output_set_receiver_header(const DY_SYN_OUTPUT *output, size_t ipt, SACHEAD *hd, real_t depsrc,
                                       size_t nlayer, const real_t (*modarr)[GRT_MODARR_NCOL]);

/**
 * 保存单个接收点位移和可选空间导数，并应用时间积分和微分
 * 极坐标单点直接保存在根目录，-Q/-U 使用索引和坐标命名子目录
 * @param[in]      output     合成输出设置
 * @param[in]      ir         接收点索引
 * @param[in,out]  trace      已填充头段的 SAC 模板，数据作为写出缓冲
 * @param[in]      data       位移在前、九个空间导数在后的连续双精度数组
 * @param[in]      dt         工作采样间隔，s，不从 SAC 头段回读
 * @param[in]      int_times  时间积分次数
 * @param[in]      dif_times  时间微分次数
 */
void grt_syn_output_save_receiver(const DY_SYN_OUTPUT *output, size_t ir, SACTRACE *trace, const real_t *data,
                                  real_t dt, int int_times, int dif_times);

/**
 * 在根目录保存源时间函数或有限源总标量矩率，长度覆盖完整破裂过程
 * @param[in]  root       合成结果根目录
 * @param[in]  dt         采样间隔，s
 * @param[in]  stf_decay  每个采样点的阻尼因子，写出时恢复物理时间函数
 * @param[in]  nfault     震源断层数量
 * @param[in]  faults     已初始化时间函数和动态标量矩的震源断层
 * @param[in]  finite     是否输出有限源总标量矩率
 */
void grt_syn_output_save_signal(const char *root, real_t dt, real_t stf_decay,
                                size_t nfault, const FINITE_FAULT *faults, bool finite);
