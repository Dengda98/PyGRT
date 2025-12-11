/**
 * @file   util.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 其它辅助函数
 * 
 */

#pragma once 

#include <stdio.h>
#include <stdbool.h>

#include "grt/common/const.h"
#include "grt/common/model.h"

/**
 * 指定分隔符，从一串字符串中分割出子字符串数组
 * 
 * @param[in]     string     原字符串
 * @param[in]     delim      分隔符
 * @param[out]    size       分割后的子字符串数组长度
 * 
 * @return   子字符串数组
 */
char ** grt_string_split(const char *string, const char *delim, size_t *size);

/**
 * 从文本文件中，将每行内容读入字符串数组
 * 
 * @param[in,out]     fp       文件指针
 * @param[out]        size     读入的字符串数组长度
 * 
 * @return   字符串数组
 * 
 */
char ** grt_string_from_file(FILE *fp, size_t *size);

/**
 * 判断字符串是否由特定的若个字符组成（充分条件）
 * 
 * @param[in]    str      待检查的字符串
 * @param[in]    alws     允许的字符集合
 * 
 * @return  是否符合
 */
bool grt_string_composed_of(const char *str, const char *alws);

/**
 * 指定分隔符，获得字符串的分割出的子字符串数。
 * 相当于是 grt_string_split 函数的简化版本
 * 
 * @param[in]     string     原字符串
 * @param[in]     delim      分隔符
 * 
 * @return   子字符串数
 */
int grt_string_ncols(const char *string, const char* delim);

/**
 * 从路径字符串中找到用/或\\分隔的最后一项
 * 
 * @param[in]    path     路径字符串指针
 * 
 * @return   指向最后一项字符串的指针
 */
const char* grt_get_basename(const char* path);


/**
 * 去除字符串首尾空白
 * 
 * @param[in,out]     str    字符串
 */
void grt_trim_whitespace(char* str);


/**
 * 检查是否为注释行或空行
 * 
 * @param[in]     line    读入一行的字符串
 */
bool grt_is_comment_or_empty(const char* line);


/**
 * 由于 Windows MSYS2 环境没有 getline 函数（即使定义了 _GNU_SOURCE）
 * 所以这里需要使用自定义的 getline 函数，参数与 POSIX 定义相同
 */
ssize_t grt_getline(char **lineptr, size_t *n, FILE *stream);