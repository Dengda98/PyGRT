/**
 * @file   signals.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12
 * 
 *                   
 */


#pragma once

#include <stdbool.h>


#define GRT_SIG_PARABOLA 'p'   ///< 抛物波代号
#define GRT_SIG_TRAPEZOID 't'  ///< 梯形波代号
#define GRT_SIG_RICKER   'r'   ///< 雷克子波信号
#define GRT_SIG_CUSTOM   '0'   ///< 自定义时间函数代码


/**
 * 检查时间函数的类型设置和参数设置是否符合要求
 * 
 * @param[in]      tftype     单个字符，指代时间函数类型
 * @param[in]      tfparams   时间函数参数
 * 
 * @return     检查是否通过
 */
bool grt_check_tftype_tfparams(const char tftype, const char *tfparams);

/**
 * 获得时间函数，要求提前运行check_tftype_tfparams函数以检查参数
 * 
 * @param[out]      TFnt       返回的点数
 * @param[in]       dt         时间间隔
 * @param[in]       tftype     单个字符，指代时间函数类型
 * @param[in]       tfparams   时间函数参数
 * 
 * @return     时间函数指针
 */
float * grt_get_time_function(int *TFnt, float dt, const char tftype, const char *tfparams);


/**
 * 时域线性卷积，要求提前运行check_tftype_tfparams函数以检查参数
 * 卷积结果会原地写入数组。
 * 
 * @param[in,out]  arr         待卷积的信号
 * @param[in]      nt          信号点数
 * @param[in]      dt          信号点时间间隔
 * @param[in]      tftype      单个字符，指代时间函数类型
 * @param[in]      tfparams    时间函数参数
 * @param[out]     TFarr       指向时间函数的指针的指针
 * @param[out]     TFnt        返回的时间函数点数
 */
void grt_linear_convolve_time_function(float *arr, int nt, float dt, const char tftype, const char *tfparams, float **TFarr, int *TFnt);


/**
 * 时间序列卷积函数，只卷积x的长度
 * 
 * @param[in]    x            长信号数组
 * @param[in]    nx           长信号点数
 * @param[in]    h            短信号数组
 * @param[in]    nh           短信号点数
 * @param[out]   y            输出数组
 * @param[in]    ny           输出数组点数
 * @param[in]    iscircular   是否使用循环卷积
 */
void grt_oaconvolve(float *x, int nx, float *h, int nh, float *y, int ny, bool iscircular);


/**
 * 计算某序列整个梯形积分值
 * 
 * @param[in]     x     信号数组 
 * @param[in]     nx    数组长度
 * @param[in]     dt    时间间隔
 * 
 * @return    积分结果
 */
float grt_trap_area(const float *x, int nx, float dt);


/**
 * 使用梯形法对时间序列积分
 * 
 * @param[in,out]     x     信号数组 
 * @param[in]         nx    数组长度
 * @param[in]         dt    时间间隔
 */
void grt_trap_integral(float *x, int nx, float dt);

/**
 * 对时间序列做中心一阶差分
 * 
 * @param[in,out]     x     信号数组 
 * @param[in]         nx    数组长度
 * @param[in]         dt    时间间隔
 */
void grt_differential(float *x, int nx, float dt);




/**
 * 生成抛物线波
 * 
 * @param[in]        dt        采样间隔
 * @param[in,out]    tlen      信号时长
 * @param[out]       Nt        返回的点数
 * 
 * @return   float指针
 */
float * grt_get_parabola_wave(float dt, float *Tlen, int *Nt);



/**
 * 生成梯形波或三角波
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
 *  O       T1                 T2     T3                T
 * 
 * @endverbatim
 * 
 * 
 * @param[in]        dt        采样间隔
 * @param[in,out]    T1        上坡截止时刻
 * @param[in,out]    T2        平台截止时刻
 * @param[in,out]    T3        下坡截止时刻
 * @param[out]       Nt        返回的点数
 * 
 * @return   float指针
 */
float * grt_get_trap_wave(float dt, float *T1, float *T2, float *T3, int *Nt);



/**
 * 生成雷克子波
 * 
 * \f[ f(t)=(1-2 \pi^2 f_0^2 (t-t_0)^2 ) e^{ - \pi^2 f_0^2 (t-t_0)^2} \f]
 * 
 * @param[in]     dt        采样间隔
 * @param[in]     f0        主频
 * @param[out]    Nt        返回的点数
 * 
 * @return   float指针
 */
float * grt_get_ricker_wave(float dt, float f0, int *Nt);


/**
 * 从文件中读入自定义时间函数
 * 
 * @param[out]    Nt        返回的点数
 * @param[in]     tfparams  文件路径
 * 
 * @return   float指针
 */
float * grt_get_custom_wave(int *Nt, const char *tfparams);

/**
 * 专用于在Python端释放C中申请的内存
 * 
 * @param[out]     pt    指针
 */
void grt_free1d(void *pt);