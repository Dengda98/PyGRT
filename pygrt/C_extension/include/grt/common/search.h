/**
 * @file   search.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 *                   
 */

#pragma once

#include <stdbool.h>
#include "grt/common/const.h"

// 定义 X 宏，为多个类型定义查找函数
#define __FOR_EACH_REAL \
    X(real_t)  X(float)  X(double)

#define __FOR_EACH_INT \
    X(size_t)

/**
 * 该函数对输入数组进行线性搜索，找到目标值时返回其索引。
 * 如果目标值在数组中未找到，则返回 -1。
 *
 * @param[in] array   要搜索的数组。
 * @param[in] size    数组的大小（元素个数）。
 * @param[in] target  要查找的目标值。
 * 
 * @return idx    目标值的索引，如果未找到则返回 -1。
 *
 * @note 如果数组中存在多个目标值，该函数返回第一个匹配的索引。
 * 
 */
#define X(T) \
ssize_t grt_findElement_##T(const T *array, size_t size, T target);

__FOR_EACH_REAL
__FOR_EACH_INT
#undef X



/**
 * 搜索数组中最接近目标值且小于目标值的索引。
 * 如果目标值在数组中未找到，则返回 -1。
 *
 * @param[in] array   要搜索的数组。
 * @param[in] size    数组的大小（元素个数）。
 * @param[in] target  要查找的目标值。
 * 
 * @return idx    目标值的索引，如果未找到则返回 -1。
 *
 * @note 如果数组中存在多个目标值，该函数返回第一个匹配的索引。
 * 
 */
#define X(T) \
ssize_t grt_findLessEqualClosest_##T(const T *array, size_t size, T target);

__FOR_EACH_REAL
#undef X



/**
 * 搜索数组中最接近目标值的索引。
 *
 * @param[in] array   要搜索的数组。
 * @param[in] size    数组的大小（元素个数）。
 * @param[in] target  要查找的目标值。
 * 
 * @return idx    目标值的索引
 *
 * @note 如果数组中存在多个目标值，该函数返回第一个匹配的索引。
 * 
 */
#define X(T) \
size_t grt_findClosest_##T(const T *array, size_t size, T target);

__FOR_EACH_REAL
#undef X


/**
 * 搜索数组的最值，返回其索引。
 *
 * @param[in] array   要搜索的数组。
 * @param[in] size    数组的大小（元素个数）。
 * 
 * @return idx    目标值的索引。
 *
 * @note 如果数组中存在相同最值，该函数返回第一个匹配的索引。
 * 
 */
#define X(T) \
size_t grt_findMin_##T(const T *array, size_t size);\
size_t grt_findMax_##T(const T *array, size_t size);\

__FOR_EACH_REAL
__FOR_EACH_INT
#undef X


/**
 * 比较函数
 * 
 * @param[in]   a    元素 a 地址 
 * @param[in]   b    元素 b 地址 
 * 
 * @return flag   比较结果，(1) a > b, (0) a == b, (-1) a < b
 */
#define X(T) \
int grt_compare_##T(const void *a, const void *b);\
int grt_argcompare_##T(const void *a, const void *b, void *arr);\

__FOR_EACH_REAL
__FOR_EACH_INT
#undef X
#undef __FOR_EACH_REAL
#undef __FOR_EACH_INT



/**
 * 在有序数组中插入元素，元素类型和数组类型需匹配
 * 
 * @param[in,out]   arr          有序数组地址
 * @param[in,out]   size         数组大小地址
 * @param[in]       capacity     数组最大容量
 * @param[in]       target       元素地址
 * @param[in]       elementSize  元素和数组内元素的字节长度
 * @param[in]       ascending    升序(true) 或 降序(false)
 * @param[in]       compare      比较函数
 * 
 * @return pos   插入位置的索引
 */
ssize_t grt_insertOrdered(
    void *arr, size_t *size, size_t capacity, const void *target, size_t elementSize, bool ascending,
    int (*compare)(const void *, const void *));

/** 元素比较函数，语义与 qsort 的比较函数相同 */
typedef int (*grt_compare_fn)(const void *a, const void *b);

/**
 * 计算任意元素类型数组的稳定升序排序索引（argsort）
 *
 * @param[in]  base          输入数组的首元素地址，不会被修改
 * @param[in]  n             数组元素个数
 * @param[in]  element_size  每个元素的字节长度，语义与 qsort 的 size 参数相同
 * @param[in]  compare       元素比较函数，传入 base 中两个元素的地址
 * @param[out] indices       输出索引数组，调用方需提前分配至少 n 个 size_t 空间
 *                           结果满足 base[indices[i]] 按 compare 升序排列
 *
 * @return 0 表示成功；-1 表示参数无效或字节数溢出；-2 表示内存分配失败
 *
 * @note 对 compare 返回 0 的元素，结果保持其在原数组中的相对顺序
 */
int grt_argsort(
    const void *base, size_t n, size_t element_size,
    grt_compare_fn compare, size_t *indices);

/**
 * 在升序数组上定位 1D 线性插值的左右邻点与权重
 *
 * n==1 时要求 |q-x[0]| 在容差内，此时 i0=i1=0、w=0；
 * q 越界（超出端点容差）返回 false
 *
 * @param[in]   x    升序坐标数组
 * @param[in]   n    长度
 * @param[in]   q    查询值
 * @param[out]  i0   左邻索引
 * @param[out]  i1   右邻索引
 * @param[out]  w    权重，使值 = (1-w)*x[i0] + w*x[i1]
 * @return      是否在范围内
 */
bool grt_locateLinearInterp(
    const real_t *x, size_t n, real_t q,
    size_t *i0, size_t *i1, real_t *w);
