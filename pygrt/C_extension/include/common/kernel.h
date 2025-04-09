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


typedef void (*KernelFunc) (
    const MODEL1D *mod1d, MYCOMPLEX omega, MYREAL k,
    MYCOMPLEX EXP_qwv[3][3], MYCOMPLEX VF_qwv[3][3], MYCOMPLEX HF_qwv[3][3], MYCOMPLEX DC_qwv[3][3],
    bool calc_uiz,
    MYCOMPLEX EXP_uiz_qwv[3][3], MYCOMPLEX VF_uiz_qwv[3][3], MYCOMPLEX HF_uiz_qwv[3][3], MYCOMPLEX DC_uiz_qwv[3][3]);