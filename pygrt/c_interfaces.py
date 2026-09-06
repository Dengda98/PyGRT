"""
    :file:     c_interfaces.py
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
    :date:     2024-07-24

    该文件包括 C 库的调用接口

"""

import os
from ctypes import POINTER, c_char_p, c_double, c_float, c_int, c_size_t, c_void_p, cdll


FPOINTER = POINTER(c_float)
IPOINTER = POINTER(c_int)
REAL = c_double
PREAL = POINTER(REAL)


libgrt = cdll.LoadLibrary(
    os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "C_extension/lib/libgrt.so",
    )
)
"""libgrt 库"""


C_grt_free = libgrt.grt_free1d
"""释放在 C 中申请的内存"""
C_grt_free.restype = None
C_grt_free.argtypes = [c_void_p]


C_grt_get_trap_wave = libgrt.grt_get_trap_wave
"""梯形波"""
C_grt_get_trap_wave.restype = FPOINTER
C_grt_get_trap_wave.argtypes = [
    c_float,
    FPOINTER,
    FPOINTER,
    FPOINTER,
    IPOINTER,
]


C_grt_get_parabola_wave = libgrt.grt_get_parabola_wave
"""抛物波"""
C_grt_get_parabola_wave.restype = FPOINTER
C_grt_get_parabola_wave.argtypes = [c_float, FPOINTER, IPOINTER]


C_grt_get_ricker_wave = libgrt.grt_get_ricker_wave
"""雷克子波"""
C_grt_get_ricker_wave.restype = FPOINTER
C_grt_get_ricker_wave.argtypes = [c_float, c_float, IPOINTER]


C_grt_solve_lamb1 = libgrt.grt_solve_lamb1
"""使用广义闭合解求解第一类 Lamb 问题"""
C_grt_solve_lamb1.restype = None
C_grt_solve_lamb1.argtypes = [
    REAL,
    PREAL,
    c_int,
    REAL,
    PREAL,
]


C_grt_solve_lamb2 = libgrt.grt_solve_lamb2
"""使用广义闭合解求解第二类 Lamb 问题"""
C_grt_solve_lamb2.restype = None
C_grt_solve_lamb2.argtypes = [
    REAL,
    PREAL,
    c_int,
    REAL,
    REAL,
    REAL,
    REAL,
    PREAL,
    PREAL,
    PREAL,
]


C_grt_compute_travt1d_from_file = libgrt.grt_compute_travt1d_from_file
"""从模型文件计算多个震中距的初至 P/S 走时"""
C_grt_compute_travt1d_from_file.restype = PREAL
C_grt_compute_travt1d_from_file.argtypes = [
    c_char_p,
    REAL,
    REAL,
    PREAL,
    c_size_t,
]
