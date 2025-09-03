/**
 * @file   sacio2.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-31
 * 
 *    在已有的sacio基础上进行部分函数的封装
 * 
 */

#pragma once 

#include "grt/common/sacio.h"

/**
 * 读取SAC头段变量
 * 
 * @param[in]       command        当前程序命令名称
 * @param[in]       name           SAC文件路径
 * @param[out]      hd             SAC头段变量结构体
 */
void grt_read_SAC_HEAD(const char *command, const char *name, SACHEAD *hd);


/**
 * 读取SAC文件
 * 
 * @param[in]       command       当前程序命令名称
 * @param[in]       name          SAC文件路径
 * @param[out]      hd            SAC头段变量结构体
 * @param[in,out]   arrout        预分配内存，不需要则设为NULL
 * 
 * @return     浮点数指针
 */
float * grt_read_SAC(const char *command, const char *name, SACHEAD *hd, float *arrout);