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



void grt_modsum_grn_spec(MODEL1D *mod1d, const DISPER_TYPE wtype, EIGENFN_INFO *eigfnmet, GRNSPEC *grn);