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
    const MODEL1D *mod1d, MYCOMPLEX omega, MYREAL k, MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM],
    bool calc_uiz, MYCOMPLEX QWV_uiz[SRC_M_NUM][QWV_NUM], MYINT *stats);