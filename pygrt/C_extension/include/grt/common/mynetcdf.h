/**
 * @file   mynetcdf.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 使用 netcdf 库的一些自定义宏
 * 
 */

#pragma once

#include <netcdf.h>

#include "grt/common/const.h"

#define GRT_NC_CHECK(call) ({\
    int status = (call); \
    if (status != NC_NOERR) { \
        GRTRaiseError("NetCDF error at %s:%d: %s\n", \
                __FILE__, __LINE__, nc_strerror(status)); \
    } \
})


#define GRT_NC_MYINT  NC_INT
#define GRT_NC_FUNC_MYINT(func)  func##_int

#ifdef GRT_USE_FLOAT 
    #define GRT_NC_MYREAL  NC_FLOAT
    #define GRT_NC_FUNC_MYREAL(func)  func##_float
#else 
    #define GRT_NC_MYREAL  NC_DOUBLE
    #define GRT_NC_FUNC_MYREAL(func)  func##_double
#endif

