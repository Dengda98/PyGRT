/**
 * @file   kernel.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-06
 * 
 *    动态或静态下计算核函数的函数指针
 * 
 */


#pragma once 


#include "grt/common/model.h"

/**
 * 计算核函数的函数指针，动态与静态的接口一致
 */
typedef void (*GRT_KernelFunc) (
    const GRT_MODEL1D *mod1d, MYCOMPLEX omega, MYREAL k, MYCOMPLEX QWV[GRT_SRC_M_NUM][GRT_QWV_NUM],
    bool calc_uiz, MYCOMPLEX QWV_uiz[GRT_SRC_M_NUM][GRT_QWV_NUM], MYINT *stats);