/**
 * @file   modgrn.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 根据本征值计算面波解
 *   
 *                 
 */

#pragma once

#include "grt/common/model.h"
#include "grt/dynamic/grnspec.h"
#include "grt/modal/eigenfn.h"


/**
 * 根据频散结果，利用模态叠加法计算不同震源激发的面波格林函数
 * 
 * @param[in,out]      mod1d       模型结构体指针
 * @param[in]          wtype       频散类型
 * @param[in,out]      eigfnmet    本征函数数据结构体指针
 * @param[out]         grn         格林函数频谱
 */
void grt_modsum_grn_spec(MODEL1D *mod1d, const DISPER_TYPE wtype, EIGENFN_INFO *eigfnmet, GRNSPEC *grn);