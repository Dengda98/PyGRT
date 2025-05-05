/**
 * @file   search.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 *                   
 */

#include <stdlib.h>
#include <stdbool.h>

#include "common/search.h"
#include "common/const.h"

static bool _gt_(MYREAL a1, MYREAL a2) { return a1 > a2; }
static bool _lt_(MYREAL a1, MYREAL a2) { return a1 < a2; }


MYINT findElement_MYINT(const MYINT array[], MYINT size, MYINT target) {
    for (MYINT i = 0; i < size; ++i) {
        if (array[i] == target) {
            return i;  // 找到目标元素，返回索引
        }
    }
    return -1;  // 未找到目标元素，返回-1
}

MYINT findLessEqualClosest_MYREAL(const MYREAL array[], MYINT size, MYREAL target) {
    MYINT ires=-1;
    MYREAL mindist=-1.0, dist=0.0;
    for (MYINT i = 0; i < size; ++i) {
        dist = target-array[i];
        if(dist >= 0.0 && (mindist < 0.0 || dist < mindist)){
            ires = i;
            mindist = dist;
        }
    }
    return ires;
}

MYINT findClosest_MYREAL(const MYREAL array[], MYINT size, MYREAL target) {
    MYINT ires=0;
    MYREAL mindist=-1.0, dist=0.0;
    for (MYINT i = 0; i < size; ++i) {
        dist = fabs(target-array[i]);
        if(mindist < 0.0 || dist < mindist){
            ires = i;
            mindist = dist;
        }
    }
    return ires;
}

MYINT findMinMax_MYREAL(const MYREAL array[], MYINT size, bool isMax) {
    MYREAL rmax = array[0];
    MYINT idx=0;
    bool (*_func)(MYREAL, MYREAL) = (isMax)? _gt_ : _lt_;
    for(MYINT ir=0; ir<size; ++ir){
        if(_func(array[ir], rmax)){
            rmax = array[ir];
            idx = ir;
        }
    }
    return idx;
}

