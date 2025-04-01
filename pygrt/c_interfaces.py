"""
    :file:     c_interfaces.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包括 C库的调用接口  

"""


import os
from ctypes import *

from .c_structures import * 

FPOINTER = POINTER(c_float)
IPOINTER = POINTER(c_int)


c_PGRN = POINTER(c_GRN)

libgrt = cdll.LoadLibrary(
    os.path.join(
        os.path.abspath(os.path.dirname(__file__)), 
        "C_extension/lib/libgrt.so"))
"""libgrt库"""


C_integ_grn_spec = libgrt.integ_grn_spec
"""C库中计算格林函数的主函数 integ_grn_spec, 详见C API同名函数"""
C_integ_grn_spec.argtypes = [
    POINTER(c_PyModel1D), c_int, c_int, c_int, PREAL,       
    c_int, PREAL, REAL,
    REAL, REAL, REAL, c_bool, REAL, REAL,
    c_bool,

    POINTER(c_PGRN*2),
    POINTER(c_PGRN*2),
    POINTER(c_PGRN*3),
    POINTER(c_PGRN*2),
    POINTER(c_PGRN*3),
    POINTER(c_PGRN*3),

    c_bool,
    # uiz
    POINTER(c_PGRN*2),
    POINTER(c_PGRN*2),
    POINTER(c_PGRN*3),
    POINTER(c_PGRN*2),
    POINTER(c_PGRN*3),
    POINTER(c_PGRN*3),
    # uir
    POINTER(c_PGRN*2),
    POINTER(c_PGRN*2),
    POINTER(c_PGRN*3),
    POINTER(c_PGRN*2),
    POINTER(c_PGRN*3),
    POINTER(c_PGRN*3),

    c_char_p,
    c_int, 
    POINTER(c_int)
]


C_set_num_threads = libgrt.set_num_threads
"""设置多线程数"""
C_set_num_threads.restype = None 
C_set_num_threads.argtypes = [c_int]


def set_num_threads(n):
    r'''
        定义计算使用的多线程数

        :param       n:    线程数
    '''
    C_set_num_threads(n)


C_compute_travt1d = libgrt.compute_travt1d
"""计算1D层状半空间的初至波走时"""
C_compute_travt1d.restype = REAL 
C_compute_travt1d.argtypes = [
    PREAL, PREAL, c_int, 
    c_int, c_int, REAL
]



# -------------------------------------------------------------------
#                      C函数定义的时间函数
# -------------------------------------------------------------------
C_free = libgrt.free1d
"""释放在C中申请的内存"""
C_free.restype = None
C_free.argtypes = [c_void_p]

C_get_trap_wave = libgrt.get_trap_wave
"""梯形波"""
C_get_trap_wave.restype = FPOINTER
C_get_trap_wave.argtypes = [c_float, FPOINTER, FPOINTER, FPOINTER, IPOINTER]

C_get_parabola_wave = libgrt.get_parabola_wave
"""抛物波"""
C_get_parabola_wave.restype = FPOINTER
C_get_parabola_wave.argtypes = [c_float, FPOINTER, IPOINTER]

C_get_ricker_wave = libgrt.get_ricker_wave
"""雷克子波"""
C_get_ricker_wave.restype = FPOINTER
C_get_ricker_wave.argtypes = [c_float, c_float, IPOINTER]