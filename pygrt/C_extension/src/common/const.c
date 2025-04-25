/**
 * @file   const.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-25
 * 
 * 将全局变量放在该文件中
 */

#include "common/const.h"


/** 分别对应爆炸源(0阶)，垂直力源(0阶)，水平力源(1阶)，剪切源(0,1,2阶) */ 
const MYINT GRT_SRC_M_ORDERS[GRT_SRC_M_COUNTS] = {0, 0, 1, 0, 1, 2};

/** 不同震源，不同阶数的名称简写，用于命名 */
const char *GRT_SRC_M_NAME_ABBR[GRT_SRC_M_COUNTS] = {"EX", "VF", "HF", "DD", "DS", "SS"};

/** ZRT三分量代号 */
const char GRT_ZRTchs[] = {'Z', 'R', 'T'};

/** ZNE三分量代号 */
const char GRT_ZNEchs[] = {'Z', 'N', 'E'};
