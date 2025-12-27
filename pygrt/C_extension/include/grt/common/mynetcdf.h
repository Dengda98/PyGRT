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

#define NC_CHECK(call) ({\
    int status = (call); \
    if (status != NC_NOERR) { \
        GRTRaiseError("NetCDF error: %s\n", nc_strerror(status)); \
    } \
})

