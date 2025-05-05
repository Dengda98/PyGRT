/**
 * @file   const.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-25
 * 
 * 将全局变量放在该文件中
 */

#include "common/const.h"


/** 分别对应爆炸源(0阶)，垂直力源(0阶)，水平力源(1阶)，剪切源(0,1,2阶) */ 
const MYINT SRC_M_ORDERS[SRC_M_NUM] = {0, 0, 1, 0, 1, 2};

/** 不同震源类型使用的格林函数类型，0为Gij，1为格林函数导数Gij,k */
const MYINT SRC_M_GTYPES[SRC_M_NUM] = {1, 0, 0, 1, 1, 1};

/** 不同震源，不同阶数的名称简写，用于命名 */
const char *SRC_M_NAME_ABBR[SRC_M_NUM] = {"EX", "VF", "HF", "DD", "DS", "SS"};

/** q, w, v 名称代号 */
const char qwvchs[] = {'q', 'w', 'v'};

/** ZRT三分量代号 */
const char ZRTchs[] = {'Z', 'R', 'T'};

/** ZNE三分量代号 */
const char ZNEchs[] = {'Z', 'N', 'E'};
