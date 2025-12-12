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

#include "grt/common/const.h"



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
void grt_write_stats(FILE *f0, real_t k, const cplxChnlGrid QWV);


/**
 * 从二进制核函数文件读出一个数据块，写入到文本文件中
 * 
 * @param[in,out]     bf0    二进制文件指针，如果为NULL则打印标题
 * @param[out]        af0    文本文件指针
 * @param[in]         col0_name    第一列名称，一般为 k 或 c
 * 
 * @return   0表示读取成功，-1表示读取结果/失败
 */
int grt_extract_stats(FILE *bf0, FILE *af0, const char *col0_name);



/**
 * 记录峰谷平均法的峰谷位置
 * 
 * @param[out]    f0         二进制文件指针 
 * @param[in]     Kpt        最终收敛积分值使用的波峰波谷位置
 * @param[in]     Fpt        最终收敛积分值使用的波峰波谷幅值
 * 
 * @note     文件记录的积分值与最终的结果还差一系列的系数，
 *           记录其值主要用于参考其变化趋势。
 * 
 */
void grt_write_stats_ptam(
    FILE *f0, 
    realIntegGrid Kpt[GRT_PTAM_PT_MAX],
    cplxIntegGrid Fpt[GRT_PTAM_PT_MAX]);


/**
 * 从二进制峰谷位置文件读出一个数据块，写入到文本文件中
 * 
 * @param[in,out]     bf0    二进制文件指针，如果为NULL则打印标题
 * @param[out]        af0    文本文件指针
 * 
 * @return   0表示读取成功，-1表示读取结果/失败
 */
int grt_extract_stats_ptam(FILE *bf0, FILE *af0);

