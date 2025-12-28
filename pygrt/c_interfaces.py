"""
    :file:     c_interfaces.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包括 C库的调用接口  

"""


import os
from ctypes import \
    c_double, c_float, c_int, c_size_t, c_bool, c_char_p, c_void_p,\
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


C_grt_integ_grn_spec = libgrt.grt_integ_grn_spec
"""C库中计算格林函数的主函数 integ_grn_spec, 详见C API同名函数"""
C_grt_integ_grn_spec.argtypes = [
    POINTER(c_GRT_MODEL1D), c_size_t, c_size_t, PREAL,       
    c_size_t, PREAL, REAL, c_bool, POINTER(c_K_INTEG_METHOD),
    c_bool, c_bool,
    POINTER((PCPLX*CHANNEL_NUM)*SRC_M_NUM),
    POINTER((PCPLX*CHANNEL_NUM)*SRC_M_NUM),
    POINTER((PCPLX*CHANNEL_NUM)*SRC_M_NUM),

    c_char_p,
    c_size_t, 
    POINTER(c_size_t)
]


C_grt_integ_static_grn = libgrt.grt_integ_static_grn
"""计算静态格林函数"""
C_grt_integ_static_grn.restype = None
C_grt_integ_static_grn.argtypes = [
    POINTER(c_GRT_MODEL1D), c_size_t, PREAL, POINTER(c_K_INTEG_METHOD),
    c_bool,
    POINTER((REAL*CHANNEL_NUM)*SRC_M_NUM),
    POINTER((REAL*CHANNEL_NUM)*SRC_M_NUM),
    POINTER((REAL*CHANNEL_NUM)*SRC_M_NUM),
    c_char_p
]


C_grt_set_num_threads = libgrt.grt_set_num_threads
"""设置多线程数"""
C_grt_set_num_threads.restype = None 
C_grt_set_num_threads.argtypes = [c_int]


def set_num_threads(n):
    r'''
        定义计算使用的多线程数

        :param       n:    线程数
    '''
    C_grt_set_num_threads(n)


C_grt_compute_travt1d = libgrt.grt_compute_travt1d
"""计算1D层状半空间的初至波走时"""
C_grt_compute_travt1d.restype = REAL 
C_grt_compute_travt1d.argtypes = [
    PREAL, PREAL, c_int, 
    c_int, c_int, REAL
]


C_grt_read_mod1d_from_file = libgrt.grt_read_mod1d_from_file
"""读取模型文件并进行预处理"""
C_grt_read_mod1d_from_file.restype = POINTER(c_GRT_MODEL1D)
C_grt_read_mod1d_from_file.argtypes = [c_char_p, c_double, c_double, c_bool]

C_grt_free_mod1d = libgrt.grt_free_mod1d
"""释放C程序中申请的 GRT_MODEL1D 结构体内存"""
C_grt_free_mod1d.restype = None
C_grt_free_mod1d.argtypes = [POINTER(c_GRT_MODEL1D)]

# -------------------------------------------------------------------
#                      C函数定义的时间函数
# -------------------------------------------------------------------
C_grt_free = libgrt.grt_free1d
"""释放在C中申请的内存"""
C_grt_free.restype = None
C_grt_free.argtypes = [c_void_p]

C_grt_get_trap_wave = libgrt.grt_get_trap_wave
"""梯形波"""
C_grt_get_trap_wave.restype = FPOINTER
C_grt_get_trap_wave.argtypes = [c_float, FPOINTER, FPOINTER, FPOINTER, IPOINTER]

C_grt_get_parabola_wave = libgrt.grt_get_parabola_wave
"""抛物波"""
C_grt_get_parabola_wave.restype = FPOINTER
C_grt_get_parabola_wave.argtypes = [c_float, FPOINTER, IPOINTER]

C_grt_get_ricker_wave = libgrt.grt_get_ricker_wave
"""雷克子波"""
C_grt_get_ricker_wave.restype = FPOINTER
C_grt_get_ricker_wave.argtypes = [c_float, c_float, IPOINTER]


# -------------------------------------------------------------------
#                      C函数定义的旋转函数
# -------------------------------------------------------------------
C_grt_rot_zxy2zrt_vec = libgrt.grt_rot_zxy2zrt_vec
"""直角坐标zxy到柱坐标zrt的矢量旋转"""
C_grt_rot_zxy2zrt_vec.restype = None
C_grt_rot_zxy2zrt_vec.argtypes = [c_double, DPOINTER]  # double, double[3]

C_grt_rot_zxy2zrt_symtensor2odr = libgrt.grt_rot_zxy2zrt_symtensor2odr
"""直角坐标zxy到柱坐标zrt的二阶对称张量旋转"""
C_grt_rot_zxy2zrt_symtensor2odr.restype = None
C_grt_rot_zxy2zrt_symtensor2odr.argtypes = [c_double, DPOINTER]  # double, double[6]

C_grt_rot_zrt2zxy_upar = libgrt.grt_rot_zrt2zxy_upar
"""柱坐标下的位移偏导 ∂u(z,r,t)/∂(z,r,t) 转到 直角坐标 ∂u(z,x,y)/∂(z,x,y)"""
C_grt_rot_zrt2zxy_upar.restype = None
C_grt_rot_zrt2zxy_upar.argtypes = [c_double, DPOINTER, DPOINTER, c_double]  # double, double[3], double[3][3], double


# -------------------------------------------------------------------
#                      C函数定义的衰减函数
# -------------------------------------------------------------------
C_grt_py_attenuation_law = libgrt.grt_py_attenuation_law
"""品质因子Q 对 波速的影响"""
C_grt_py_attenuation_law.restype = None
C_grt_py_attenuation_law.argtypes = [REAL, DPOINTER, DPOINTER]  # double, double[2], double[2]





# -------------------------------------------------------------------
#                      使用 C 函数求解 Lamb 问题
# -------------------------------------------------------------------
C_grt_solve_lamb1 = libgrt.grt_solve_lamb1
"""使用广义闭合解求解第一类 Lamb 问题"""
C_grt_solve_lamb1.restype = None
C_grt_solve_lamb1.argtypes = [
    REAL, PREAL, c_int, REAL, PREAL
]