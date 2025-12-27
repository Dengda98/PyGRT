/**
 * @file   sacio2.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-31
 * 
 *    在已有的sacio基础上进行部分函数的封装
 * 
 */

#pragma once 

#include <stdbool.h>
#include "grt/common/sacio.h"

/** 将 SAC 头段变量和数据体打包成一个结构体 */
typedef struct {
    SACHEAD hd;
    float *data;
} SACTRACE;

/**
 * 读取SAC文件
 * 
 * @param[in]       path          SAC文件路径
 * @param[in]       headonly      是否只读取头段变量
 * 
 * @return     SACTRACE 指针
 */
SACTRACE * grt_read_SACTRACE(const char *path, const bool headonly);

/**
 * 复制 SACTRACE
 * 
 * @param[in]    sac          源 SAC
 * @param[in]    zero_value   是否数据置零
 * @return    复制 SAC
 */
SACTRACE * grt_copy_SACTRACE(SACTRACE *sac, bool zero_value);

/**
 * 新建 SACTRACE
 * 
 * @param[in]     dt      时间间隔
 * @param[in]     nt      点数
 * @param[in]     b0      开始时刻
 * 
 * @return     SACTRACE 指针
 */
SACTRACE * grt_new_SACTRACE(float dt, int nt, float b0);

/** 将 SACTRACE 保存到本地 */
int grt_write_SACTRACE(const char *path, SACTRACE *sac);

/** 释放 SACTRACE 指针 */
void grt_free_SACTRACE(SACTRACE *sac);