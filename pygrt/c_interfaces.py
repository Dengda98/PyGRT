"""
    :file:     c_interfaces.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包括 C库的调用接口  

"""


import os
from ctypes import *

from .c_structures import * 


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