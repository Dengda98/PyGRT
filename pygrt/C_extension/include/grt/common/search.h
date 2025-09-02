/**
 * @file   search.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 *                   
 */

#pragma once

#include <stdbool.h>
#include "common/const.h"

/**
 * 该函数对输入的整数数组进行线性搜索，找到目标值时返回其索引。
 * 如果目标值在数组中未找到，则返回 -1。
 *
 * @param[in] array   要搜索的整数数组。
 * @param[in] size    数组的大小（元素个数）。
 * @param[in] target  要查找的目标值。
 * 
 * @return idx    目标值的索引，如果未找到则返回 -1。
 *
 * @note 如果数组中存在多个目标值，该函数返回第一个匹配的索引。
 * 
 */
MYINT findElement_MYINT(const MYINT array[], MYINT size, MYINT target);

/**
 * 搜索浮点数数组中最接近目标值且小于目标值的索引。
 * 如果目标值在数组中未找到，则返回 -1。
 *
 * @param[in] array   要搜索的浮点数数组。
 * @param[in] size    数组的大小（元素个数）。
 * @param[in] target  要查找的目标值。
 * 
 * @return idx    目标值的索引，如果未找到则返回 -1。
 *
 * @note 如果数组中存在多个目标值，该函数返回第一个匹配的索引。
 * 
 */
MYINT findLessEqualClosest_MYREAL(const MYREAL array[], MYINT size, MYREAL target);

/**
 * 搜索浮点数数组中最接近目标值的索引。
 *
 * @param[in] array   要搜索的浮点数数组。
 * @param[in] size    数组的大小（元素个数）。
 * @param[in] target  要查找的目标值。
 * 
 * @return idx    目标值的索引
 *
 * @note 如果数组中存在多个目标值，该函数返回第一个匹配的索引。
 * 
 */
MYINT findClosest_MYREAL(const MYREAL array[], MYINT size, MYREAL target);

/**
 * 搜索浮点数数组的最大或最小值，返回其索引。
 *
 * @param[in] array   要搜索的浮点数数组。
 * @param[in] size    数组的大小（元素个数）。
 * @param[in] isMax   是否要找最大值，否则找最小值。
 * 
 * @return idx    目标值的索引。
 *
 * @note 如果数组中存在相同最值，该函数返回第一个匹配的索引。
 * 
 */
MYINT findMinMax_MYREAL(const MYREAL array[], MYINT size, bool isMax);