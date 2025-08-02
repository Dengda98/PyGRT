/**
 * @file   util.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 其它辅助函数
 * 
 */

#pragma once 



/**
 * 指定分隔符，从一串字符串中分割出子字符串数组
 */
char ** string_split(const char *string, char *delim, int *size);
