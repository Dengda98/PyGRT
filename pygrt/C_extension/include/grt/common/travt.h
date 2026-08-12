/**
 * @file   travt.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-08
 * 
 *    计算一维均匀半无限层状介质的初至走时
 * 
 */

#pragma once 

#include "grt/common/const.h"

/**
 * 已知每层的厚度和速度，且震源和场点位于（虚拟）界面上,
 * 且不共享层位，即使深度相同，中间也考虑一个厚度为0的层。
 * 故当abs(isrc-ircv)==1时，说明两点位于同一物理层
 * 
 * @param[in]    Thk           每层厚度 
 * @param[in]    Vel0          每层速度
 * @param[in]    nlay          层数
 * @param[in]    isrc          震源所在层位
 * @param[in]    ircv          场点所在层位
 * @param[in]    dist          震中距
 */
real_t grt_compute_travt1d(
    const real_t *Thk, const real_t *Vel0, const size_t nlay, 
    const size_t isrc, const size_t ircv, const real_t dist);


/**
 * 从模型文件计算多个震中距的初至 P/S 走时
 *
 * 返回长度为 2*nr 的数组，按 [Tp0, Ts0, Tp1, Ts1, ...] 排列
 * 调用方需用 grt_free1d 释放；失败返回 NULL
 *
 * @param[in]    modelpath     一维分层模型文件路径
 * @param[in]    depsrc        震源深度 (km)
 * @param[in]    deprcv        台站深度 (km)
 * @param[in]    rs            震中距数组 (km)
 * @param[in]    nr            震中距个数
 */
real_t *grt_compute_travt1d_from_file(
    const char *modelpath,
    const real_t depsrc,
    const real_t deprcv,
    const real_t *rs,
    const size_t nr);
