/**
 * @file   kernel.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-06
 * 
 *    动态或静态下计算核函数的函数指针
 * 
 */


#pragma once 


#include "common/model.h"

/**
 * 计算核函数的函数指针，动态与静态的接口一致
 */
typedef void (*KernelFunc) (
    const MODEL1D *mod1d, MYCOMPLEX omega, MYREAL k, MYCOMPLEX QWV[GRT_SRC_M_COUNTS][GRT_SRC_QWV_COUNTS],
    bool calc_uiz, MYCOMPLEX QWV_uiz[GRT_SRC_M_COUNTS][GRT_SRC_QWV_COUNTS]);