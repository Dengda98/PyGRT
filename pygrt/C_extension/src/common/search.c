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