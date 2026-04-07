/**
 * @file   grt.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 集合所有相关的头文件
 * 
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
#include <ctype.h>
#include <stdbool.h>

// ls grt/*/*.h | tr "" "\n" | sort | awk '{print "#include \""$1"\""}'
#include "grt/common/RT_matrix.h"
#include "grt/common/attenuation.h"
#include "grt/common/bessel.h"
#include "grt/common/checkerror.h"
#include "grt/common/colorstr.h"
#include "grt/common/const.h"
#include "grt/common/coord.h"
#include "grt/common/logo.h"
#include "grt/common/matrix.h"
#include "grt/common/model.h"
#include "grt/common/myfftw.h"
#include "grt/common/mynetcdf.h"
#include "grt/common/progressbar.h"
#include "grt/common/radiation.h"
#include "grt/common/sacio.h"
#include "grt/common/search.h"
#include "grt/common/travt.h"
#include "grt/common/util.h"
#include "grt/common/version.h"


#include "grt/dynamic/grn.h"
#include "grt/dynamic/grnspec.h"
#include "grt/dynamic/layer.h"
#include "grt/dynamic/signals.h"
#include "grt/dynamic/source.h"
#include "grt/dynamic/syn.h"


#include "grt/static/static_grn.h"
#include "grt/static/static_layer.h"
#include "grt/static/static_source.h"


#include "grt/integral/dcm.h"
#include "grt/integral/dwm.h"
#include "grt/integral/fim.h"
#include "grt/integral/integ_method.h"
#include "grt/integral/iostats.h"
#include "grt/integral/k_integ.h"
#include "grt/integral/kernel.h"
#include "grt/integral/ptam.h"
#include "grt/integral/quadratic.h"
#include "grt/integral/safim.h"


#include "grt/lamb/elliptic.h"
#include "grt/lamb/lamb1.h"
#include "grt/lamb/lamb_util.h"
