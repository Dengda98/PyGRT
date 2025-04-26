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

#include "common/const.h"



/**
 * 将积分过程中计算的核函数写入文件
 * 
 * @param[out]    f0      二进制文件指针 
 * @param[in]     k       波数 
 * @param[in]     QWV     不同震源，不同阶数的核函数 \f$ q_m, w_m, v_m \f$
 * 
 * 
 * @note     文件记录的值均为波数积分的中间结果，与最终的结果还差一系列的系数，
 *           记录其值主要用于参考其变化趋势。
 */
void write_stats(
    FILE *f0, MYREAL k, const MYCOMPLEX QWV[SRC_M_NUM][QWV_NUM]);



/**
 * 记录峰谷平均法的峰谷位置
 * 
 * @param[out]    f0         二进制文件指针 
 * @param[in]     k          波数 
 * @param[in]     Kpt        最终收敛积分值使用的波峰波谷位置
 * @param[in]     Fpt        最终收敛积分值使用的波峰波谷幅值
 * 
 * @note     文件记录的积分值与最终的结果还差一系列的系数，
 *           记录其值主要用于参考其变化趋势。
 * 
 */
void write_stats_ptam(
    FILE *f0, MYREAL k, 
    MYREAL Kpt[SRC_M_NUM][INTEG_NUM][PTAM_MAX_PT],
    MYCOMPLEX Fpt[SRC_M_NUM][INTEG_NUM][PTAM_MAX_PT]);