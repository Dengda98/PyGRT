/**
 * @file   signals.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12
 * 
 *    时间函数生成与信号处理
 *    采样参数、时间函数及工作数组均使用 real_t
 */


#pragma once

#include <stdbool.h>

#include "grt/common/const.h"


#define GRT_SIG_IMPULSE           'i'  ///< 脉冲信号代号
#define GRT_SIG_PARABOLA          'p'  ///< 抛物波代号
#define GRT_SIG_TRAPEZOID         't'  ///< 梯形波代号
#define GRT_SIG_ASYMMETRIC_COSINE 'c'  ///< 非对称余弦波代号
#define GRT_SIG_RICKER            'r'  ///< 雷克子波信号
#define GRT_SIG_CUSTOM            '0'  ///< 自定义时间函数代码


/**
 * 检查时间函数的类型设置和参数设置是否符合要求
 * 参数可附加 +d<delay> 指定延迟，单位为 s
 * 
 * @param[in]      tftype     时间函数类型
 * @param[in]      tfparams   时间函数参数，可包含 +d<delay> 后缀
 * 
 * @return     检查是否通过
 */
bool grt_check_tftype_tfparams(const char tftype, const char *tfparams);

/**
 * 根据类型和参数生成时间函数
 * 时间函数按矩形法进行面积归一化，即 dt 乘样本和为 1（雷克子波保留单位峰值）
 * 自定义时间函数按 dt 乘样本和归一化，非单位面积时警告
 * 积分相对幅值绝对积分接近零时，无法归一化并报错
 * 可在参数末尾追加 +d<delay> 指定延迟，单位为 s
 * 
 * @param[out]      TFnt       返回的点数，包含延迟对应的前导零
 * @param[in]       dt         时间间隔，s
 * @param[in]       tftype     时间函数类型
 * @param[in]       tfparams   时间函数参数，可包含 +d<delay> 后缀
 * 
 * @return     时间函数指针
 */
real_t * grt_get_time_function(int *TFnt, real_t dt, const char tftype, const char *tfparams);


/**
 * 解析完整时间函数选项，返回包含整数采样延迟的时间函数
 *
 * 选项格式为 -Dtftype[/tfparams][+d<delay>]
 *
 * @param[in]  option         完整 -D 时间函数选项，NULL 为脉冲
 * @param[in]  dt             采样间隔，s
 * @param[out] nt             包含延迟的样本数
 */
real_t *grt_time_function_from_option(const char *option, real_t dt, int *nt);

/**
 * 将时间量化到整数采样网格
 * @param[in] time  原始时间，s
 * @param[in] dt    最终输出采样间隔，s
 */
real_t grt_sample_aligned_time(real_t time, real_t dt);

/**
 * 在时域计算离散卷积
 * 循环卷积以 ny 为周期，要求 ny 不小于 nx 和 nh，较短输入补零
 * 线性卷积只返回前 ny 个样本，超出完整卷积长度的部分补零
 *
 * @param[in]    x            输入信号数组
 * @param[in]    nx           输入信号点数
 * @param[in]    h            卷积核数组
 * @param[in]    nh           卷积核点数
 * @param[out]   y            输出数组
 * @param[in]    ny           输出数组点数
 * @param[in]    iscircular   是否使用循环卷积
 *
 * 该函数只进行离散样本求和，不包含连续卷积所需的 dt 因子
 */
void grt_oaconvolve(const real_t *x, int nx, const real_t *h, int nh, real_t *y, int ny, bool iscircular);


/**
 * 计算某序列整个梯形积分值
 * 
 * @param[in]     x     信号数组 
 * @param[in]     nx    数组长度
 * @param[in]     dt    时间间隔
 * 
 * @return    积分结果
 */
real_t grt_trap_area(const real_t *x, int nx, real_t dt);


/**
 * 使用梯形法对时间序列积分
 * 
 * @param[in,out]     x     信号数组 
 * @param[in]         nx    数组长度
 * @param[in]         dt    时间间隔
 */
void grt_trap_integral(real_t *x, int nx, real_t dt);

/**
 * 对时间序列做中心一阶差分
 * 
 * @param[in,out]     x     信号数组 
 * @param[in]         nx    数组长度
 * @param[in]         dt    时间间隔
 */
void grt_differential(real_t *x, int nx, real_t dt);




/**
 * 生成抛物线波
 * 截止时刻向上对齐到采样网格，至少需要两个采样间隔，按矩形法进行面积归一化
 * 
 * @param[in]        dt        采样间隔
 * @param[in,out]    Tlen      信号时长，返回实际采样时长
 * @param[out]       Nt        返回的点数
 * 
 * @return   real_t 指针
 */
real_t * grt_get_parabola_wave(real_t dt, real_t *Tlen, int *Nt);



/**
 * 生成梯形波、三角波或矩形波
 * T1、T2、T3 分别为上升段、平台段、下降段时长，均须非负且总时长大于零
 * T2=0 时退化为三角波，T1=T3=0 时退化为矩形波
 * 三段时长分别向上对齐到采样网格，接近网格点时只消除浮点误差
 * 非零时段至少保留一个采样间隔，零时长段保持为零
 * 所有样本按分段函数直接取值，再按矩形法进行面积归一化以匹配离散卷积
 * 
 * @verbatim
 *   ^
 *   |
 *   |
 * 1-|       --------...--------
 *   |      /                   \ 
 *   |     /                     \ 
 *   |   ...                     ...
 *   |   /                         \
 *   |  /                           \
 *   | /                             \
 *   |------+------------------+------+---------------->
 *  O       T1              T1+T2  T1+T2+T3             T
 * 
 * @endverbatim
 * 
 * 
 * @param[in]        dt        采样间隔
 * @param[in,out]    T1        上升段时长，返回实际采样时长，s
 * @param[in,out]    T2        平台段时长，返回实际采样时长，s
 * @param[in,out]    T3        下降段时长，返回实际采样时长，s
 * @param[out]       Nt        返回的点数
 * 
 * @return   real_t 指针
 */
real_t * grt_get_trap_wave(real_t dt, real_t *T1, real_t *T2, real_t *T3, int *Nt);


/**
 * 生成非对称余弦滑移速率时间函数
 * T1、T2 分别为上升段和下降段时长，均须大于零
 * 两段时长分别向上对齐到采样网格，接近网格点时消除浮点误差
 * 上升段和下降段至少各保留一个采样间隔，按矩形法进行面积归一化
 *
 * @param[in]        dt        采样间隔，s
 * @param[in,out]    T1        上升段时长，返回实际采样时长，s
 * @param[in,out]    T2        下降段时长，返回实际采样时长，s
 * @param[out]       Nt        返回的点数
 *
 * @return   real_t 指针
 */
real_t *grt_get_asymmetric_cosine_wave(real_t dt, real_t *T1, real_t *T2, int *Nt);



/**
 * 生成雷克子波
 * 峰值时刻为 1/f0，保留解析波形的单位峰值
 * 
 * \f[ f(t)=(1-2 \pi^2 f_0^2 (t-t_0)^2 ) e^{ - \pi^2 f_0^2 (t-t_0)^2} \f]
 * 
 * @param[in]     dt        采样间隔
 * @param[in]     f0        主频
 * @param[out]    Nt        返回的点数
 * 
 * @return   real_t 指针
 */
real_t * grt_get_ricker_wave(real_t dt, real_t f0, int *Nt);


/**
 * 从文件中读入自定义时间函数，每个非注释行只能包含一列振幅值
 * 序列和由调用方根据采样间隔检查
 * 
 * @param[out]    Nt        返回的点数
 * @param[in]     tfparams  文件路径
 * 
 * @return   real_t 指针
 */
real_t * grt_get_custom_wave(int *Nt, const char *tfparams);

/**
 * 释放 C 侧 malloc 的一维数组
 * 
 * @param[out]     pt    指针
 */
void grt_free1d(void *pt);
