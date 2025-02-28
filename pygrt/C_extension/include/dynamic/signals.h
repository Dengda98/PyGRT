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
 * @param      tftype     单个字符，指代时间函数类型
 * @param      tfparams   时间函数参数
 * 
 * @return     检查是否通过
 */
bool check_tftype_tfparams(const char tftype, const char *tfparams);

/**
 * 获得时间函数，要求提前运行check_tftype_tfparams函数以检查参数
 * 
 * @param      TFnt    返回的点数
 * @param      dt      时间间隔
 * @param      tftype     单个字符，指代时间函数类型
 * @param      tfparams   时间函数参数
 * 
 * @return     时间函数指针
 */
float * get_time_function(int *TFnt, float dt, const char tftype, const char *tfparams);


/**
 * 时域线性卷积，要求提前运行check_tftype_tfparams函数以检查参数
 * 卷积结果会原地写入数组。
 * 
 * @param    arr    待卷积的信号
 * @param    nt     信号点数
 * @param    dt     信号点时间间隔
 * @param      tftype     单个字符，指代时间函数类型
 * @param      tfparams   时间函数参数
 * @param      TFarr      指向时间函数的指针的指针
 * @param      TFnt       返回的时间函数点数
 */
void linear_convolve_time_function(float *arr, int nt, float dt, const char tftype, const char *tfparams, float **TFarr, int *TFnt);


/**
 * 时间序列卷积函数，只卷积x的长度
 * 
 * @param    x      长信号数组
 * @param    nx     长信号点数
 * @param    h      短信号数组
 * @param    nh     短信号点数
 * @param    y      输出数组
 * @param    ny     输出数组点数
 * @param    iscircular     是否使用循环卷积
 */
void oaconvolve(float *x, int nx, float *h, int nh, float *y, int ny, bool iscircular);


/**
 * 计算某序列整个梯形积分值
 * 
 * @param     x     信号数组 
 * @param     nx    数组长度
 * @param     dt    时间间隔
 * 
 * @return    积分结果
 */
float trap_area(const float *x, int nx, float dt);


/**
 * 使用梯形法对时间序列积分
 * 
 * @param     x     信号数组 
 * @param     nx    数组长度
 * @param     dt    时间间隔
 */
void trap_integral(float *x, int nx, float dt);

/**
 * 对时间序列做中心一阶差分
 * 
 * @param     x     信号数组 
 * @param     nx    数组长度
 * @param     dt    时间间隔
 */
void differential(float *x, int nx, float dt);




/**
 * 生成抛物线波
 * 
 * @param    dt       (in)采样间隔
 * @param    tlen     (inout)信号时长
 * @param    Nt       (out)返回的点数
 * 
 * @return   float指针
 */
float * get_parabola_wave(float dt, float *Tlen, int *Nt);



/**
 * 生成梯形波或三角波
 * 
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
 * 
 * @param    dt       (in)采样间隔
 * @param    T1       (inout)上坡截止时刻
 * @param    T2       (inout)平台截止时刻
 * @param    T3       (inout)下坡截止时刻
 * @param    Nt       (out)返回的点数
 * 
 * @return   float指针
 */
float * get_trap_wave(float dt, float *T1, float *T2, float *T3, int *Nt);



/**
 * 生成雷克子波
 * 
 * \f[ f(t)=(1-2 \pi^2 \f_0^2 (t-t_0)^2 ) \exp{ - \pi^2 \f_0^2 (t-t_0)^2} \f]
 * 
 * @param    dt       (in)采样间隔
 * @param    f0       (in)主频
 * @param    Nt       (out)返回的点数
 * 
 * @return   float指针
 */
float * get_ricker_wave(float dt, float f0, int *Nt);


/**
 * 从文件中读入自定义时间函数
 * 
 * @param    Nt       (out)返回的点数
 * @param    tfparams  文件路径
 * 
 * @return   float指针
 */
float * get_custom_wave(int *Nt, const char *tfparams);