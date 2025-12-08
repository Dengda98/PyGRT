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
 * @return   子字符串数组
 */
char ** grt_string_split(const char *string, const char *delim, size_t *size);

/**
 * 从文本文件中，将每行内容读入字符串数组
 * 
 * @param[in,out]     fp       文件指针
 * @param[out]        size     读入的字符串数组长度
 * 
 * @return   字符串数组
 * 
 */
char ** grt_string_from_file(FILE *fp, size_t *size);

/**
 * 判断字符串是否由特定的若个字符组成（充分条件）
 * 
 * @param[in]    str      待检查的字符串
 * @param[in]    alws     允许的字符集合
 * 
 * @return  是否符合
 */
bool grt_string_composed_of(const char *str, const char *alws);

/**
 * 指定分隔符，获得字符串的分割出的子字符串数。
 * 相当于是 grt_string_split 函数的简化版本
 * 
 * @param[in]     string     原字符串
 * @param[in]     delim      分隔符
 * 
 * @return   子字符串数
 */
int grt_string_ncols(const char *string, const char* delim);

/**
 * 从路径字符串中找到用/或\\分隔的最后一项
 * 
 * @param[in]    path     路径字符串指针
 * 
 * @return   指向最后一项字符串的指针
 */
const char* grt_get_basename(const char* path);


/**
 * 去除字符串首尾空白
 * 
 * @param[in,out]     str    字符串
 */
void grt_trim_whitespace(char* str);


/**
 * 检查是否为注释行或空行
 * 
 * @param[in]     line    读入一行的字符串
 */
bool grt_is_comment_or_empty(const char* line);


/**
 * 由于 Windows MSYS2 环境没有 getline 函数（即使定义了 _GNU_SOURCE）
 * 所以这里需要使用自定义的 getline 函数，参数与 POSIX 定义相同
 */
ssize_t grt_getline(char **lineptr, size_t *n, FILE *stream);


/**
 * 处理单个震中距对应的数据逆变换和SAC保存
 * 
 * @param[in]         command       模块名
 * @param[in]         mod1d         模型结构体指针
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
void grt_GF_freq2time_write_to_file(
    const char *command, const GRT_MODEL1D *mod1d, 
    const char *s_output_dir, const char *s_modelname, const char *s_depsrc, const char *s_deprcv,    
    const real_t wI, GRT_FFTW_HOLDER *pt_fh,
    const size_t nr, char *s_dists[nr], const real_t dists[nr], real_t travtPS[nr][2],
    const real_t depsrc, const real_t deprcv,
    const real_t delayT0, const real_t delayV0, const bool calc_upar,
    const bool doEX, const bool doVF, const bool doHF, const bool doDC, 
    const char *chalst,
    pt_cplxChnlGrid grn[nr], 
    pt_cplxChnlGrid grn_uiz[nr], 
    pt_cplxChnlGrid grn_uir[nr]);