/**
 * @file   progressbar.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 以下代码实现进度条的输出
 * 
 */

#include "const.h"


#define _PROGRESSBAR_WIDTH_ 45  ///< 定义进度条的长度


/**
 * 根据百分比打印进度条  
 * 
 * @param    prefix     (in)进度条前缀字符串
 * @param    percentage (in)百分比(整数)
 */
void printprogressBar(const char *prefix, MYINT percentage);