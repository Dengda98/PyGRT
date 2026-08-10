/**
 * @file   static_grn.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 * 以下代码实现的是 广义反射透射系数矩阵 计算静态格林函数，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *              透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 * 
 */



#pragma once


#include "grt/common/model.h"
#include "grt/integral/integ_process.h"


/**
 * 静态积分前准备：按震中距与用户参数填充深度无关的 K_INTEG_PROCESS 字段
 *
 * 不写入 Kproc->k0；调用方按 hs=max(|depsrc-deprcv|, GRT_MIN_DEPTH_GAP_SRC_RCV)
 * 自行设置 Kproc->k0 = k0_user * PI / hs
 *
 * @param[in]      nr            震中距数量
 * @param[in]      rs            震中距数组
 * @param[in]      Length        波数积分窗口长度因子（同 greenfn -L）；0 表示用默认值
 * @param[in]      filonLength   FIM 长度因子；<=0 表示关闭
 * @param[in]      safilonTol    SAFIM 容差；<=0 表示关闭
 * @param[in]      filonCut      FIM/SAFIM 分段偏移（同 greenfn -L+o）
 * @param[in]      keps          自动收敛判据（同 greenfn -K+e）
 * @param[in]      use_kmax_ref  是否以参考 kmax 为积分上限（同 greenfn -K+f）
 * @param[in]      convmet       波数积分收敛方法
 * @param[out]     Kproc         待填充的 K_INTEG_PROCESS
 */
void grt_prepare_static_grn(
    size_t nr, real_t *rs,
    real_t Length,
    real_t filonLength, real_t safilonTol, real_t filonCut,
    real_t keps, bool use_kmax_ref,
    int convmet,
    K_INTEG_PROCESS *Kproc);

/**
 * 积分计算Z, R, T三个分量静态格林函数的核心函数
 * 
 * @param[in,out]      mod1d            `MODEL1D` 结构体指针 
 * @param[in]      nr               震中距数量
 * @param[in]      rs               震中距数组 
 * @param[in,out]   Kproc            波数积分相关参数的结构体指针
 * @param[in]       calc_upar         是否计算位移u的空间导数
 * @param[out]      grn               浮点数数组，不同震源不同阶数的静态格林函数的Z、R、T分量
 * @param[out]      grn_uiz           浮点数数组，不同震源不同阶数的ui_z的Z、R、T分量
 * @param[out]      grn_uir           浮点数数组，不同震源不同阶数的ui_r的Z、R、T分量
 * 
 * @param[in]       statsstr           积分过程输出目录
 * 
 */
void grt_integ_static_grn(
    MODEL1D *mod1d, size_t nr, real_t *rs, K_INTEG_PROCESS *Kproc,
    bool calc_upar, 
    realChnlGrid grn[nr],
    realChnlGrid grn_uiz[nr],
    realChnlGrid grn_uir[nr],
    const char *statsstr);