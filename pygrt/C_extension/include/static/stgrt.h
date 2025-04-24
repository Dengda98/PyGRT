/**
 * @file   stgrt.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-03
 * 
 * 以下代码实现的是 广义反射透射系数矩阵+离散波数法 计算静态格林函数，参考：
 * 
 *         1. 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  
 *         2. 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、
 *              透射系数矩阵和离散波数方法[J]. 地球物理学报(3): 270-280.
 * 
 */



#pragma once


#include "common/model.h"


/**
 * 积分计算Z, R, T三个分量静态格林函数的核心函数
 * 
 * @param      pymod1d      (in)`PYMODEL1D` 结构体指针 
 * @param      nr           (in)震中距数量
 * @param      rs           (in)震中距数组 
 * @param      keps         (in)波数积分的收敛条件，要求在某震中距下所有格林函数都收敛，为负数代表不提前判断收敛，按照波数积分上限进行积分 
 * @param      k0           (in)波数积分的上限
 * @param      Length       (in)波数k积分间隔 \f$ dk=2\pi/(fabs(L)*r_{max}) \f$ 
 * @param      filonLength  (in)Filon积分间隔
 * @param      filonCut     (in)波数积分和Filon积分的分割点
 * @param      EXPgrn[nr][2]      (out)浮点数数组，爆炸源的Z、R分量频谱结果
 * @param      VFgrn[nr][2]       (out)浮点数数组，垂直力源的Z、R分量频谱结果
 * @param      HFgrn[nr][3]       (out)浮点数数组，水平力源的Z、R、T分量频谱结果
 * @param      DDgrn[nr][2]       (out)浮点数数组，45度倾滑的Z、R分量频谱结果
 * @param      DSgrn[nr][3]       (out)浮点数数组，90度倾滑的Z、R、T分量频谱结果
 * @param      SSgrn[nr][3]       (out)浮点数数组，90度走滑的Z、R、T分量频谱结果
 *
 * @param      calc_upar              (in)是否计算位移u的空间导数
 * @param      EXPgrn_uiz[nr][2]      (out)浮点数数组，爆炸源产生的ui_z(位移u对坐标z的偏导)的Z、R分量频谱结果，下同
 * @param      VFgrn_uiz[nr][2]       (out)浮点数数组，垂直力源的Z、R分量频谱结果
 * @param      HFgrn_uiz[nr][3]       (out)浮点数数组，水平力源的Z、R、T分量频谱结果
 * @param      DDgrn_uiz[nr][2]       (out)浮点数数组，45度倾滑的Z、R分量频谱结果
 * @param      DSgrn_uiz[nr][3]       (out)浮点数数组，90度倾滑的Z、R、T分量频谱结果
 * @param      SSgrn_uiz[nr][3]       (out)浮点数数组，90度走滑的Z、R、T分量频谱结果
 * @param      EXPgrn_uir[nr][2]      (out)浮点数数组，爆炸源产生的ui_r(位移u对坐标r的偏导)的Z、R分量频谱结果，下同
 * @param      VFgrn_uir[nr][2]       (out)浮点数数组，垂直力源的Z、R分量频谱结果
 * @param      HFgrn_uir[nr][3]       (out)浮点数数组，水平力源的Z、R、T分量频谱结果
 * @param      DDgrn_uir[nr][2]       (out)浮点数数组，45度倾滑的Z、R分量频谱结果
 * @param      DSgrn_uir[nr][3]       (out)浮点数数组，90度倾滑的Z、R、T分量频谱结果
 * @param      SSgrn_uir[nr][3]       (out)浮点数数组，90度走滑的Z、R、T分量频谱结果
 * 
 */
void integ_static_grn(
    PYMODEL1D *pymod1d, MYINT nr, MYREAL *rs, MYREAL vmin_ref, MYREAL keps, MYREAL k0, MYREAL Length,
    MYREAL filonLength, MYREAL filonCut, 

    // 返回值，维度2代表Z、R分量，维度3代表Z、R、T分量
    MYREAL EXPgrn[nr][2], // EXZ, EXR 的实部和虚部
    MYREAL VFgrn[nr][2],  // VFZ, VFR 的实部和虚部
    MYREAL HFgrn[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    MYREAL DDgrn[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    MYREAL DSgrn[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    MYREAL SSgrn[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]

    bool calc_upar,
    MYREAL EXPgrn_uiz[nr][2], // EXZ, EXR 的实部和虚部
    MYREAL VFgrn_uiz[nr][2],  // VFZ, VFR 的实部和虚部
    MYREAL HFgrn_uiz[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    MYREAL DDgrn_uiz[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    MYREAL DSgrn_uiz[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    MYREAL SSgrn_uiz[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]
    MYREAL EXPgrn_uir[nr][2], // EXZ, EXR 的实部和虚部
    MYREAL VFgrn_uir[nr][2],  // VFZ, VFR 的实部和虚部
    MYREAL HFgrn_uir[nr][3],  // HFZ, HFR, HFT 的实部和虚部
    MYREAL DDgrn_uir[nr][2],  // DDZ, DDR 的实部和虚部      [DD: 45-dip slip]
    MYREAL DSgrn_uir[nr][3],  // DSZ, DSR, DST 的实部和虚部 [DS: 90-dip slip]
    MYREAL SSgrn_uir[nr][3],  // SSZ, SSR, SST 的实部和虚部 [SS: strike slip]

    const char *statsstr // 积分结果输出
);