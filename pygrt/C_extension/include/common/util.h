/**
 * @file   util.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 其它辅助函数
 * 
 */

#pragma once 

#include "common/const.h"

/**
 * 指定分隔符，从一串字符串中分割出子字符串数组
 * 
 * @param[in]     string     原字符串
 * @param[in]     delim      分隔符
 * @param[out]    size       分割后的子字符串数组长度
 * 
 * @return   split    子字符串数组
 */
char ** string_split(const char *string, const char *delim, int *size);


/**
 * 从路径字符串中找到用/或\\分隔的最后一项
 * 
 * @param    path     路径字符串指针
 * 
 * @return   指向最后一项字符串的指针
 */
const char* get_basename(const char* path);