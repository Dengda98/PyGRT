/**
 * @file   iostats.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 对积分过程中的核函数和积分值进行记录                  
 */

#pragma once 

#include <stdbool.h>
#include <stdio.h>

#include "const.h"



/**
 * 对离散波数法以及Filon积分的中间结果进行记录，其中涉及两种数组形状：
 *    + [3][3]. 存储的是核函数，第一个维度3代表阶数m=0,1,2，第二个维度3代表三类系数qm,wm,vm  
 *    + [3][4]. 存储的是该dk区间内的积分值，维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型
 * 
 * @param    f0     (out)二进制文件指针 
 * @param    k       (in)波数 
 * @param    EXP_qwv[3][3]    (in)爆炸源核函数
 * @param    VF_qwv[3][3]     (in)垂直力源核函数
 * @param    HF_qwv[3][3]     (in)水平力源核函数
 * @param    DC_qwv[3][3]     (in)双力偶源核函数
 * @param    EXP_J[3][4]      (in)爆炸源，该dk区间内的积分值，下同
 * @param    VF_J[3][4]       (in)垂直力源
 * @param    HF_J[3][4]       (in)水平力源
 * @param    DC_J[3][4]       (in)双力偶源
 * 
 * 
 * @note     文件记录的值均为波数积分的中间结果，与最终的结果还差一系列的系数，
 *           记录其值主要用于参考其变化趋势。
 */
void write_stats(
    FILE *f0, MYREAL k, 
    const MYCOMPLEX EXP_qwv[3][3], const MYCOMPLEX VF_qwv[3][3], 
    const MYCOMPLEX HF_qwv[3][3],  const MYCOMPLEX DC_qwv[3][3], 
    const MYCOMPLEX EXP_J[3][4], const MYCOMPLEX VF_J[3][4], 
    const MYCOMPLEX HF_J[3][4],  const MYCOMPLEX DC_J[3][4]);



/**
 * 对峰谷平均法的中间结果进行记录，其中[3][4]的数组形状代表存储的
 * 是最终的积分值，维度3代表阶数m=0,1,2，维度4代表4种类型的F(k,w)Jm(kr)k的类型
 * 
 * @param    f0     (out)二进制文件指针 
 * @param    k       (in)波数 
 * @param    maxNpt  (in)波峰+波谷的数量(本质是计算中提前预设的量，见ptam.c文件中PTA_method函数)
 * @param    EXPpt[3][4][maxNpt]    (in)爆炸源，最终收敛积分值使用的波峰波谷值，下同
 * @param    VFpt[3][4][maxNpt]     (in)垂直力源
 * @param    HFpt[3][4][maxNpt]     (in)水平力源
 * @param    DCpt[3][4][maxNpt]     (in)双力偶源
 * @param    kEXPpt[3][4][maxNpt]    (in)爆炸源，最终收敛积分值使用的波峰波谷值对应的波数k值，下同
 * @param    kVFpt[3][4][maxNpt]     (in)垂直力源
 * @param    kHFpt[3][4][maxNpt]     (in)水平力源
 * @param    kDCpt[3][4][maxNpt]     (in)双力偶源
 * 
 * @note     文件记录的积分值与最终的结果还差一系列的系数，
 *           记录其值主要用于参考其变化趋势。
 * 
 */
void write_stats_ptam(
    FILE *f0, MYREAL k, MYINT maxNpt, 
    const MYCOMPLEX EXPpt[3][4][maxNpt], const MYCOMPLEX VFpt[3][4][maxNpt],
    const MYCOMPLEX HFpt[3][4][maxNpt],  const MYCOMPLEX DCpt[3][4][maxNpt],
    const MYREAL kEXPpt[3][4][maxNpt], const MYREAL kVFpt[3][4][maxNpt],
    const MYREAL kHFpt[3][4][maxNpt],  const MYREAL kDCpt[3][4][maxNpt]);