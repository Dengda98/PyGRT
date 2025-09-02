/**
 * @file   util.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 其它辅助函数
 * 
 */

#pragma once 

#include <stdbool.h>

#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/common/myfftw.h"

/**
 * 指定分隔符，从一串字符串中分割出子字符串数组
 * 
 * @param[in]     string     原字符串
 * @param[in]     delim      分隔符
 * @param[out]    size       分割后的子字符串数组长度
 * 
 * @return   split    子字符串数组
 */
char ** string_split(const char *string, const char *delim, int *size);


/**
 * 从路径字符串中找到用/或\\分隔的最后一项
 * 
 * @param    path     路径字符串指针
 * 
 * @return   指向最后一项字符串的指针
 */
const char* get_basename(const char* path);



/**
 * 处理单个震中距对应的数据逆变换和SAC保存
 * 
 * @param[in]         command       模块名
 * @param[in]         pymod         模型结构体指针
 * @param[in]         s_output_dir  保存目录（调用前已创建）
 * @param[in]         s_modelname   模型名称
 * @param[in]         s_depsrc      震源深度字符串
 * @param[in]         s_deprcv      接收深度字符串
 * @param[in]         wI            虚频率
 * @param[in,out]     pt_fh         FFTW结构体
 * @param[in]         nr            震中距数量
 * @param[in]         s_dists       输入的震中距字符串数组
 * @param[in]         dists         震中距数组
 * @param[out]        travtPS       保存不同震中距的初至P、S
 * @param[in]         depsrc        震源深度
 * @param[in]         deprcv        接收深度
 * @param[in]         delayT0       延迟时间
 * @param[in]         delayV0       参考速度
 * @param[in]         calc_upar     是否计算位移偏导
 * @param[in]         doEX          是否保存爆炸源结果
 * @param[in]         doVF          是否保存垂直力源结果
 * @param[in]         doHF          是否保存水平力源结果
 * @param[in]         doDC          是否保存剪切力源结果
 * @param[in]         doDC          是否保存剪切力源结果
 * @param[in]         chalst        要保存的分量字符串
 * @param[in]         grn           格林函数频谱结果
 * @param[in]         grn_uiz       格林函数对z偏导频谱结果
 * @param[in]         grn_uir       格林函数对r偏导频谱结果
 * 
 */
void GF_freq2time_write_to_file(
    const char *command, const PYMODEL1D *pymod, 
    const char *s_output_dir, const char *s_modelname, const char *s_depsrc, const char *s_deprcv,    
    const MYREAL wI, FFTW_HOLDER *pt_fh,
    const MYINT nr, char *s_dists[nr], const MYREAL dists[nr], MYREAL travtPS[nr][2],
    const MYREAL depsrc, const MYREAL deprcv,
    const MYREAL delayT0, const MYREAL delayV0, const bool calc_upar,
    const bool doEX, const bool doVF, const bool doHF, const bool doDC, 
    const char *chalst,
    MYCOMPLEX *grn[nr][SRC_M_NUM][CHANNEL_NUM], 
    MYCOMPLEX *grn_uiz[nr][SRC_M_NUM][CHANNEL_NUM], 
    MYCOMPLEX *grn_uir[nr][SRC_M_NUM][CHANNEL_NUM]);