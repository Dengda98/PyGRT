/**
 * @file   sacio.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-31
 * 
 *    在已有的sacio基础上进行部分函数的封装
 * 
 */

#pragma once 

#include <stdbool.h>
#include "grt/common/const.h"
#include "grt/common/sacio.h_"

/**
 * PyGRT 自定义 SAC 头段字段的语义访问宏
 *
 * hd 为 SACHEAD 指针，GET 宏也接受 const 指针，SET 宏作为独立语句使用
 * 未定义值沿用 SAC_FLOAT_UNDEF 或 SAC_INT_UNDEF，字段布局和文件格式不变
 * 所有模块使用统一字段含义，每个实参只求值一次
 */
/** 公共虚频参数 */
#define GRT_SACHEAD_GET_IMAG_FREQ(hd)                 ((hd)->user0)           ///< 虚频系数 omega_I，rad/s
#define GRT_SACHEAD_SET_IMAG_FREQ(hd, value)          do { GRT_SACHEAD_GET_IMAG_FREQ(hd) = (value); } while(0)  ///< 设置虚频系数 omega_I，rad/s

/** 格林函数和合成波形的源点、接收点介质 */
#define GRT_SACHEAD_GET_RCV_VP(hd)                    ((hd)->user1)           ///< 接收点 P 波速度，km/s
#define GRT_SACHEAD_GET_RCV_VS(hd)                    ((hd)->user2)           ///< 接收点 S 波速度，km/s
#define GRT_SACHEAD_GET_RCV_RHO(hd)                   ((hd)->user3)           ///< 接收点密度，g/cm^3
#define GRT_SACHEAD_GET_RCV_QP_INV(hd)                ((hd)->user4)           ///< 接收点 P 波 Q 值的倒数
#define GRT_SACHEAD_GET_RCV_QS_INV(hd)                ((hd)->user5)           ///< 接收点 S 波 Q 值的倒数
#define GRT_SACHEAD_GET_SRC_VP(hd)                    ((hd)->user6)           ///< 震源点 P 波速度，km/s
#define GRT_SACHEAD_GET_SRC_VS(hd)                    ((hd)->user7)           ///< 震源点 S 波速度，km/s
#define GRT_SACHEAD_GET_SRC_RHO(hd)                   ((hd)->user8)           ///< 震源点密度，g/cm^3

#define GRT_SACHEAD_SET_RCV_VP(hd, value)                     do { GRT_SACHEAD_GET_RCV_VP(hd) = (value); } while(0)  ///< 设置接收点 P 波速度，km/s
#define GRT_SACHEAD_SET_RCV_VS(hd, value)                     do { GRT_SACHEAD_GET_RCV_VS(hd) = (value); } while(0)  ///< 设置接收点 S 波速度，km/s
#define GRT_SACHEAD_SET_RCV_RHO(hd, value)                    do { GRT_SACHEAD_GET_RCV_RHO(hd) = (value); } while(0)  ///< 设置接收点密度，g/cm^3
#define GRT_SACHEAD_SET_RCV_QP_INV(hd, value)                 do { GRT_SACHEAD_GET_RCV_QP_INV(hd) = (value); } while(0)  ///< 设置接收点 P 波 Q 值的倒数
#define GRT_SACHEAD_SET_RCV_QS_INV(hd, value)                 do { GRT_SACHEAD_GET_RCV_QS_INV(hd) = (value); } while(0)  ///< 设置接收点 S 波 Q 值的倒数
#define GRT_SACHEAD_SET_SRC_VP(hd, value)                     do { GRT_SACHEAD_GET_SRC_VP(hd) = (value); } while(0)  ///< 设置震源点 P 波速度，km/s
#define GRT_SACHEAD_SET_SRC_VS(hd, value)                     do { GRT_SACHEAD_GET_SRC_VS(hd) = (value); } while(0)  ///< 设置震源点 S 波速度，km/s
#define GRT_SACHEAD_SET_SRC_RHO(hd, value)                    do { GRT_SACHEAD_GET_SRC_RHO(hd) = (value); } while(0)  ///< 设置震源点密度，g/cm^3

/** 合成波形的接收点坐标、机制和有限接收断层索引 */
#define GRT_SACHEAD_GET_RCV_NORTH(hd)                 ((hd)->unused1)         ///< 接收点北向坐标，km
#define GRT_SACHEAD_GET_RCV_EAST(hd)                  ((hd)->unused2)         ///< 接收点东向坐标，km
#define GRT_SACHEAD_GET_RCV_STRIKE(hd)                ((hd)->unused3)         ///< 接收机制走向，degree
#define GRT_SACHEAD_GET_RCV_DIP(hd)                   ((hd)->unused4)         ///< 接收机制倾角，degree
#define GRT_SACHEAD_GET_RCV_RAKE(hd)                  ((hd)->unused5)         ///< 接收机制滑动角，degree
#define GRT_SACHEAD_GET_RCV_FAULT_INDEX(hd)           ((hd)->unused11)        ///< 所属接收断层索引，从 0 开始

#define GRT_SACHEAD_SET_RCV_NORTH(hd, value)                  do { GRT_SACHEAD_GET_RCV_NORTH(hd) = (value); } while(0)  ///< 设置接收点北向坐标，km
#define GRT_SACHEAD_SET_RCV_EAST(hd, value)                   do { GRT_SACHEAD_GET_RCV_EAST(hd) = (value); } while(0)  ///< 设置接收点东向坐标，km
#define GRT_SACHEAD_SET_RCV_STRIKE(hd, value)                 do { GRT_SACHEAD_GET_RCV_STRIKE(hd) = (value); } while(0)  ///< 设置接收机制走向，degree
#define GRT_SACHEAD_SET_RCV_DIP(hd, value)                    do { GRT_SACHEAD_GET_RCV_DIP(hd) = (value); } while(0)  ///< 设置接收机制倾角，degree
#define GRT_SACHEAD_SET_RCV_RAKE(hd, value)                   do { GRT_SACHEAD_GET_RCV_RAKE(hd) = (value); } while(0)  ///< 设置接收机制滑动角，degree
#define GRT_SACHEAD_SET_RCV_FAULT_INDEX(hd, value)            do { GRT_SACHEAD_GET_RCV_FAULT_INDEX(hd) = (value); } while(0)  ///< 设置所属接收断层索引，从 0 开始

/** 接收函数专用参数 */
#define GRT_SACHEAD_GET_RCVFN_RAYP(hd)                ((hd)->resp0)           ///< 接收函数水平射线参数，s/km
#define GRT_SACHEAD_GET_RCVFN_GAUSS_ALPHA(hd)         ((hd)->resp1)           ///< 接收函数高斯滤波参数，Hz
#define GRT_SACHEAD_GET_RCVFN_INCIDENT(hd)            ((hd)->resp2)           ///< 接收函数入射波类型，GRT_RCVFN_INCIDENT
#define GRT_SACHEAD_GET_RCVFN_OUTPUT(hd)              ((hd)->resp3)           ///< 接收函数输出类型，GRT_RCVFN_OUTPUT

#define GRT_SACHEAD_SET_RCVFN_RAYP(hd, value)                 do { GRT_SACHEAD_GET_RCVFN_RAYP(hd) = (value); } while(0)  ///< 设置接收函数水平射线参数，s/km
#define GRT_SACHEAD_SET_RCVFN_GAUSS_ALPHA(hd, value)          do { GRT_SACHEAD_GET_RCVFN_GAUSS_ALPHA(hd) = (value); } while(0)  ///< 设置接收函数高斯滤波参数，Hz
#define GRT_SACHEAD_SET_RCVFN_INCIDENT(hd, value)             do { GRT_SACHEAD_GET_RCVFN_INCIDENT(hd) = (value); } while(0)  ///< 设置接收函数入射波类型，GRT_RCVFN_INCIDENT
#define GRT_SACHEAD_SET_RCVFN_OUTPUT(hd, value)               do { GRT_SACHEAD_GET_RCVFN_OUTPUT(hd) = (value); } while(0)  ///< 设置接收函数输出类型，GRT_RCVFN_OUTPUT

/** 内存中的 SAC 记录，头段保留文件格式，波形使用 real_t，仅读写时转换为 float */
typedef struct {
    SACHEAD hd;     ///< SAC 格式头段
    real_t *data;   ///< 工作波形，只读取头段时为 NULL
} SACTRACE;

/**
 * 读取 SAC 文件，将 float 波形转换为 real_t
 * 
 * @param[in]       path          SAC文件路径
 * @param[in]       headonly      是否只读取头段变量
 * 
 * @return     SACTRACE 指针
 */
SACTRACE * grt_read_SACTRACE(const char *path, const bool headonly);

/**
 * 复制 SACTRACE
 * 
 * @param[in]    sac          源 SAC
 * @param[in]    zero_value   是否数据置零
 * @return    复制 SAC
 */
SACTRACE * grt_copy_SACTRACE(SACTRACE *sac, bool zero_value);

/**
 * 新建 SACTRACE
 * 
 * @param[in]     dt      时间间隔
 * @param[in]     nt      点数
 * @param[in]     b0      开始时刻
 * 
 * @return     SACTRACE 指针
 */
SACTRACE * grt_new_SACTRACE(real_t dt, int nt, real_t b0);

/**
 * 设置 SAC 起点，并根据采样点数和采样间隔更新终点
 *
 * @param[in,out]  hd     SAC 头段
 * @param[in]      begin  新的开始时刻，s
 */
void grt_sachead_set_begin(SACHEAD *hd, real_t begin);

/** 将 SACTRACE 保存到本地，临时转换为 float，不修改工作波形，写出失败直接报错 */
void grt_write_SACTRACE(const char *path, SACTRACE *sac);

/** 释放 SACTRACE 指针 */
void grt_free_SACTRACE(SACTRACE *sac);
