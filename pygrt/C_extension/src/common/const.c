/**
 * @file   const.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-25
 * 
 * 将全局变量放在该文件中
 */

#include "grt/common/const.h"

/** 当前模块名，根据所调用模块进行切换 */
const char *GRT_MODULE_NAME = "grt";



/** 分别对应爆炸源(0阶)，垂直力源(0阶)，水平力源(1阶)，剪切源(0,1,2阶) */ 
const int GRT_SRC_M_ORDERS[GRT_SRC_M_NUM] = {0, 0, 1, 0, 1, 2};

/** 不同震源类型使用的格林函数类型，0为Gij，1为格林函数导数Gij,k */
const int GRT_SRC_M_GTYPES[GRT_SRC_M_NUM] = {1, 0, 0, 1, 1, 1};

/** 不同震源，不同阶数的名称简写，用于命名 */
const char *GRT_SRC_M_NAME_ABBR[GRT_SRC_M_NUM] = {"EX", "VF", "HF", "DD", "DS", "SS"};

/** q, w, v 名称代号 */
const char GRT_QWV_CODES[] = {'q', 'w', 'v'};

/** ZRT三分量代号 */
const char GRT_ZRT_CODES[] = {'Z', 'R', 'T'};

/** ZNE三分量代号 */
const char GRT_ZNE_CODES[] = {'Z', 'N', 'E'};



void grt_set_num_threads(int num_threads){
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif
}