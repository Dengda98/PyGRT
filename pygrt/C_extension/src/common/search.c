/**
 * @file   search.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 *                   
 */

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "grt/common/search.h"
#include "grt/common/const.h"

// 定义 X 宏，为多个类型定义查找函数
#define __FOR_EACH_REAL \
    X(real_t)  X(float)  X(double)

#define __FOR_EACH_INT \
    X(size_t)


#define X(T) \
ssize_t grt_findElement_##T(const T *array, size_t size, T target) {\
    for (size_t i = 0; i < size; ++i) {\
        if (array[i] == target) {\
            /** 找到目标元素，返回索引 */\
            return i;\
        }\
    }\
    /** 未找到目标元素，返回-1 */\
    return -1; \
}

__FOR_EACH_REAL
__FOR_EACH_INT
#undef X



#define X(T) \
ssize_t grt_findLessEqualClosest_##T(const T *array, size_t size, T target) {\
    ssize_t ires=-1;\
    T mindist=-1, dist=0;\
    for (size_t i = 0; i < size; ++i) {\
        dist = target-array[i];\
        if(dist >= 0 && (mindist < 0 || dist < mindist)){\
            ires = i;\
            mindist = dist;\
        }\
    }\
    return ires;\
}

__FOR_EACH_REAL
#undef X



#define X(T) \
size_t grt_findClosest_##T(const T *array, size_t size, T target) {\
    size_t ires=0;\
    T mindist=-1, dist=0;\
    for (size_t i = 0; i < size; ++i) {\
        dist = fabs(target-array[i]);\
        if(mindist < 0 || dist < mindist){\
            ires = i;\
            mindist = dist;\
        }\
    }\
    return ires;\
}

__FOR_EACH_REAL
#undef X


#define X(T) \
size_t grt_findMin_##T(const T *array, size_t size) {\
    T rval = array[0];\
    size_t idx=0;\
    for(size_t ir=0; ir<size; ++ir){\
        if(array[ir] < rval){\
            rval = array[ir];\
            idx = ir;\
        }\
    }\
    return idx;\
}\
size_t grt_findMax_##T(const T *array, size_t size) {\
    T rval = array[0];\
    size_t idx=0;\
    for(size_t ir=0; ir<size; ++ir){\
        if(array[ir] > rval){\
            rval = array[ir];\
            idx = ir;\
        }\
    }\
    return idx;\
}

__FOR_EACH_REAL
__FOR_EACH_INT
#undef X


#define X(T) \
int grt_compare_##T(const void *a, const void *b) {\
    T vala = *(T *)a;\
    T valb = *(T *)b;\
    if(vala > valb){\
        return 1;\
    } else if (vala < valb){\
        return -1;\
    } else {\
        return 0;\
    }\
}\
int grt_argcompare_##T(const void *a, const void *b, void *arr) {\
    size_t i1 = *(size_t *)a;\
    size_t i2 = *(size_t *)b;\
    T vala = ((T *)arr)[i1];\
    T valb = ((T *)arr)[i2];\
    if(vala > valb){\
        return 1;\
    } else if (vala < valb){\
        return -1;\
    } else {\
        return 0;\
    }\
}


__FOR_EACH_REAL
__FOR_EACH_INT
#undef X
#undef __FOR_EACH_REAL
#undef __FOR_EACH_INT


ssize_t grt_insertOrdered(
    void *arr, size_t *size, size_t capacity, const void *target, size_t elementSize, bool ascending,
    int (*compare)(const void *, const void *))
{    
    int sgn = (ascending)? 1 : -1;

    // 数组满载情况下，只可能插入更小(升序)或更大(降序)的数值
    if(*size == capacity && sgn*compare(target, arr+(*size-1)*elementSize) >= 0) return -1;

    // 找到插入位置
    size_t pos=*size;
    for(size_t i=0; i<*size; ++i){
        if(sgn*compare(target, arr+i*elementSize) < 0){
            pos = i;
            break;
        }
    }

    // 截断式插入，防止越界
    size_t lastpos = *size;
    if(lastpos >= capacity){
        lastpos = capacity-1;
    } else {
        ++(*size);
    }
    pos = GRT_MIN(pos, lastpos);

    // 移动插入位置后的元素
    memmove(arr + (pos + 1) * elementSize,
            arr + pos * elementSize,
            (lastpos - pos) * elementSize);

    // 插入新元素
    memcpy(arr + pos * elementSize, target, elementSize);

    return pos;
}



// 索引-元素地址配对结构体，用于 qsort 排序
// compare 保存在每个配对项中，使标准 qsort 的比较函数无需依赖全局状态
typedef struct {
    const unsigned char *element;
    size_t index;
    grt_compare_fn compare;
} GRTArgSortPair;

// 标准 qsort 比较函数：先按元素值升序，再按原始索引排序以保持稳定性
static int compare_argsort_pair(const void *a, const void *b)
{
    const GRTArgSortPair *pa = a;
    const GRTArgSortPair *pb = b;
    int result = pa->compare(pa->element, pb->element);

    if (result != 0) return result;
    return (pa->index > pb->index) - (pa->index < pb->index);
}

/** 计算任意元素类型数组的稳定升序排序索引（argsort）*/
int grt_argsort(
    const void *base, size_t n, size_t element_size,
    grt_compare_fn compare, size_t *indices)
{
    if (n == 0) return 0;
    if (base == NULL || element_size == 0 || compare == NULL || indices == NULL) return -1;
    if (n > (size_t)-1 / element_size ||
        n > (size_t)-1 / sizeof(GRTArgSortPair)) return -1;

    GRTArgSortPair *pairs = malloc(n * sizeof(*pairs));
    if (pairs == NULL) return -2;

    const unsigned char *bytes = base;
    for (size_t i = 0; i < n; i++) {
        pairs[i].element = bytes + i * element_size;
        pairs[i].index = i;
        pairs[i].compare = compare;
    }

    qsort(pairs, n, sizeof(*pairs), compare_argsort_pair);

    for (size_t i = 0; i < n; i++) {
        indices[i] = pairs[i].index;
    }

    free(pairs);
    return 0;
}


bool grt_locateLinearInterp(
    const real_t *x, size_t n, real_t q,
    size_t *i0, size_t *i1, real_t *w)
{
    const real_t atol = 1e-8;

    if(x == NULL || n == 0 || i0 == NULL || i1 == NULL || w == NULL){
        return false;
    }
    if(n == 1){
        if(fabs(q - x[0]) > atol) return false;
        *i0 = *i1 = 0;
        *w = 0.0;
        return true;
    }
    if(q < x[0] - atol || q > x[n - 1] + atol){
        return false;
    }
    if(q <= x[0]){
        *i0 = *i1 = 0;
        *w = 0.0;
        return true;
    }
    if(q >= x[n - 1]){
        *i0 = *i1 = n - 1;
        *w = 0.0;
        return true;
    }
    size_t i = 0;
    for(; i + 1 < n; ++i){
        if(q <= x[i + 1] + atol) break;
    }
    *i0 = i;
    *i1 = i + 1;
    real_t dx = x[*i1] - x[*i0];
    *w = (fabs(dx) < atol) ? 0.0 : (q - x[*i0]) / dx;
    return true;
}
