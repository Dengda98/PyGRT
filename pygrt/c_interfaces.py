"""
    :file:     c_interfaces.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包括 C库的调用接口  

"""


import os
from ctypes import \
    c_double, c_float, c_int, c_bool, c_char_p, c_void_p,\
    POINTER, cdll

from .c_structures import * 

FPOINTER = POINTER(c_float)
IPOINTER = POINTER(c_int)
DPOINTER = POINTER(c_double)


libgrt = cdll.LoadLibrary(
    os.path.join(
        os.path.abspath(os.path.dirname(__file__)), 
        "C_extension/lib/libgrt.so"))
"""libgrt库"""


C_integ_grn_spec = libgrt.integ_grn_spec
"""C库中计算格林函数的主函数 integ_grn_spec, 详见C API同名函数"""
C_integ_grn_spec.argtypes = [
    POINTER(c_PyModel1D), c_int, c_int, PREAL,       
    c_int, PREAL, REAL,
    REAL, REAL, REAL, REAL, REAL, REAL, REAL, REAL,
    c_bool,

    POINTER((PCPLX*CHANNEL_NUM)*SRC_M_NUM),

    c_bool,
    # uiz
    POINTER((PCPLX*CHANNEL_NUM)*SRC_M_NUM),
    # uir
    POINTER((PCPLX*CHANNEL_NUM)*SRC_M_NUM),

    c_char_p,
    c_int, 
    POINTER(c_int)
]


C_integ_static_grn = libgrt.integ_static_grn
"""计算静态格林函数"""
C_integ_static_grn.restype = None
C_integ_static_grn.argtypes = [
    POINTER(c_PyModel1D), c_int, PREAL, REAL, REAL, REAL, REAL, 
    REAL, REAL, REAL, 
    POINTER((REAL*CHANNEL_NUM)*SRC_M_NUM),
    c_bool,
    POINTER((REAL*CHANNEL_NUM)*SRC_M_NUM),
    POINTER((REAL*CHANNEL_NUM)*SRC_M_NUM),
    c_char_p
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


C_read_pymod_from_file = libgrt.read_pymod_from_file
"""读取模型文件并进行预处理"""
C_read_pymod_from_file.restype = POINTER(c_PyModel1D)
C_read_pymod_from_file.argtypes = [c_char_p, c_char_p, c_double, c_double, c_bool]

C_free_pymod = libgrt.free_pymod
"""释放C程序中申请的PYMODEL1D结构体内存"""
C_free_pymod.restype = None
C_free_pymod.argtypes = [POINTER(c_PyModel1D)]

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


# -------------------------------------------------------------------
#                      C函数定义的旋转函数
# -------------------------------------------------------------------
C_rot_zxy2zrt_vec = libgrt.rot_zxy2zrt_vec
"""直角坐标zxy到柱坐标zrt的矢量旋转"""
C_rot_zxy2zrt_vec.restype = None
C_rot_zxy2zrt_vec.argtypes = [c_double, DPOINTER]  # double, double[3]

C_rot_zxy2zrt_symtensor2odr = libgrt.rot_zxy2zrt_symtensor2odr
"""直角坐标zxy到柱坐标zrt的二阶对称张量旋转"""
C_rot_zxy2zrt_symtensor2odr.restype = None
C_rot_zxy2zrt_symtensor2odr.argtypes = [c_double, DPOINTER]  # double, double[6]

C_rot_zrt2zxy_upar = libgrt.rot_zrt2zxy_upar
"""柱坐标下的位移偏导 ∂u(z,r,t)/∂(z,r,t) 转到 直角坐标 ∂u(z,x,y)/∂(z,x,y)"""
C_rot_zrt2zxy_upar.restype = None
C_rot_zrt2zxy_upar.argtypes = [c_double, DPOINTER, DPOINTER, c_double]  # double, double[3], double[3][3], double


# -------------------------------------------------------------------
#                      C函数定义的衰减函数
# -------------------------------------------------------------------
C_py_attenuation_law = libgrt.py_attenuation_law
"""品质因子Q 对 波速的影响"""
C_py_attenuation_law.restype = None
C_py_attenuation_law.argtypes = [REAL, DPOINTER, DPOINTER]  # double, double[2], double[2]

