/**
 * @file   signals.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12
 * 
 *    时间函数生成与信号处理
 */


#pragma once

#include <stdbool.h>


#define GRT_SIG_IMPULSE  'i'   ///< 脉冲信号代号
#define GRT_SIG_PARABOLA 'p'   ///< 抛物波代号
#define GRT_SIG_TRAPEZOID 't'  ///< 梯形波代号
#define GRT_SIG_RICKER   'r'   ///< 雷克子波信号
#define GRT_SIG_CUSTOM   '0'   ///< 自定义时间函数代码


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
float * grt_get_time_function(int *TFnt, float dt, const char tftype, const char *tfparams);


/**
 * 解析完整时间函数选项，返回包含整数采样延迟的时间函数
 *
 * 选项格式为 -Dtftype[/tfparams][+d<delay>]
 *
 * @param[in]  option         完整 -D 时间函数选项，NULL 为脉冲
 * @param[in]  dt             采样间隔，s
 * @param[out] nt             包含延迟的样本数
 */
float *grt_time_function_from_option(const char *option, float dt, int *nt);

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
 *
 * 该函数只进行离散样本求和，不包含连续卷积所需的 dt 因子
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
 * 截止时刻向上对齐到采样网格，至少需要两个采样间隔，按矩形法进行面积归一化
 * 
 * @param[in]        dt        采样间隔
 * @param[in,out]    Tlen      信号时长，返回实际采样时长
 * @param[out]       Nt        返回的点数
 * 
 * @return   float指针
 */
float * grt_get_parabola_wave(float dt, float *Tlen, int *Nt);



/**
 * 生成梯形波、三角波或矩形波
 * T1=T2 时平台时长为零，退化为三角波
 * T1=0 时上坡时长为零，T2=T3 时下坡时长为零，两者同时满足时退化为矩形波
 * 截止时刻向上对齐到采样网格，接近网格点或相等的时刻按浮点容差处理
 * 非零上坡和下坡至少保留一个采样间隔，矩形波也至少保留一个采样间隔
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
 *  O       T1                 T2     T3                T
 * 
 * @endverbatim
 * 
 * 
 * @param[in]        dt        采样间隔
 * @param[in,out]    T1        上坡截止时刻，返回实际采样时刻
 * @param[in,out]    T2        平台截止时刻，返回实际采样时刻
 * @param[in,out]    T3        下坡截止时刻，返回实际采样时刻
 * @param[out]       Nt        返回的点数
 * 
 * @return   float指针
 */
float * grt_get_trap_wave(float dt, float *T1, float *T2, float *T3, int *Nt);



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
 * @return   float指针
 */
float * grt_get_ricker_wave(float dt, float f0, int *Nt);


/**
 * 从文件中读入自定义时间函数，每个非注释行只能包含一列振幅值
 * 序列和由调用方根据采样间隔检查
 * 
 * @param[out]    Nt        返回的点数
 * @param[in]     tfparams  文件路径
 * 
 * @return   float指针
 */
float * grt_get_custom_wave(int *Nt, const char *tfparams);

/**
 * 释放 C 侧 malloc 的一维数组
 * 
 * @param[out]     pt    指针
 */
void grt_free1d(void *pt);
