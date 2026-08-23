/**
 * @file   okada.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * Okada 静态位移解析解（https://www.bosai.go.jp/e/dc3d.html）
 *
 */

#pragma once

#include "grt/common/const.h"

/**
 * 计算均匀半空间中埋藏点源产生的静态位移及位移偏导
 *
 * @param[in]   alpha   均匀介质常数
 * @param[in]   x       Okada 坐标系中的观测点 X 坐标
 * @param[in]   y       Okada 坐标系中的观测点 Y 坐标
 * @param[in]   z       观测点 Z 坐标，向上为正
 * @param[in]   depth   点源深度
 * @param[in]   dip     震源倾角，单位为度
 * @param[in]   pot1    走向滑动 potency
 * @param[in]   pot2    倾向滑动 potency
 * @param[in]   pot3    张裂 potency
 * @param[in]   pot4    爆炸源 potency
 * @param[out]  u       位移，顺序为局部 X、Y、Z 分量
 * @param[out]  upar    位移偏导，第一维为偏导坐标，第二维为位移分量
 * @return              0 表示正常，1 表示奇异点，2 表示输入观测点 Z 为正
 */
int grt_okada_dc3d0(
    real_t alpha, real_t x, real_t y, real_t z, real_t depth, real_t dip,
    real_t pot1, real_t pot2, real_t pot3, real_t pot4,
    real_t u[3], real_t upar[3][3]);

/**
 * 计算均匀半空间中矩形有限断层产生的静态位移及位移偏导
 *
 * Okada 坐标系为右手坐标系，X 沿走向，Y 为上倾方向的水平投影，Z 向上
 * 观测点必须满足 z <= 0，所有长度必须使用相同单位
 *
 * @param[in]   alpha   均匀介质常数
 * @param[in]   x       Okada 坐标系中的观测点 X 坐标
 * @param[in]   y       Okada 坐标系中的观测点 Y 坐标
 * @param[in]   z       观测点 Z 坐标，向上为正
 * @param[in]   depth   断层参考点深度
 * @param[in]   dip     断层倾角，单位为度
 * @param[in]   al1     断层沿走向坐标范围下界
 * @param[in]   al2     断层沿走向坐标范围上界
 * @param[in]   aw1     断层沿上倾方向坐标范围下界
 * @param[in]   aw2     断层沿上倾方向坐标范围上界
 * @param[in]   disl1   走向滑动位错，正值表示左旋
 * @param[in]   disl2   倾向滑动位错，正值表示逆冲
 * @param[in]   disl3   张裂位错，正值表示张开
 * @param[out]  u       位移，顺序为局部 X、Y、Z 分量
 * @param[out]  upar    位移偏导，第一维为偏导坐标，第二维为位移分量
 * @return              0 表示正常，1 表示奇异点，2 表示输入观测点 Z 为正
 */
int grt_okada_dc3d(
    real_t alpha, real_t x, real_t y, real_t z, real_t depth, real_t dip,
    real_t al1, real_t al2, real_t aw1, real_t aw2,
    real_t disl1, real_t disl2, real_t disl3,
    real_t u[3], real_t upar[3][3]);
