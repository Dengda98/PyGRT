"""
    :file:     signals.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包括一些常见的时间信号，最高幅值均为1    


"""


import numpy as np  
import numpy.ctypeslib as npct
from ctypes import byref, cast

from .c_interfaces import *

__all__ = [
    "gen_triangle_wave",
    "gen_parabola_wave",
    "gen_trap_wave",
    "gen_ricker_wave",
]

def gen_triangle_wave(vlen, dt):
    '''
        generate triangle-shape wave  

        :param    vlen:    signal length (s)  
        :param    dt:      time interval (s)   

        :return: 
            - **wave** -    amplitude sequence
    '''
    return gen_trap_wave(vlen/2.0, vlen/2.0, vlen, dt)


def gen_parabola_wave(vlen, dt):
    '''
        generate parabola-shape wave    

        :param    vlen:    signal length (s)  
        :param    dt:      time interval (s)   
        
        :return: 
            - **wave** -    amplitude sequence
    '''
    ct1 = c_float(vlen)
    cnt = c_int(0)

    carr = C_grt_get_parabola_wave(dt, byref(ct1), byref(cnt))
    arr = npct.as_array(carr, shape=(cnt.value,)).copy()

    C_grt_free(carr)

    return arr

def gen_trap_wave(t1, t2, t3, dt):
    '''
        generate trapezoid-shape wave  

        :param    t1:      ramp-up cutoff time (s)  
        :param    t2:      plateau cutoff time (s)  
        :param    t3:      ramp-down cutoff time (s)  
        :param    dt:      time interval (s)   

        :return: 
            - **wave** -    amplitude sequence
    '''
    ct1 = c_float(t1)
    ct2 = c_float(t2)
    ct3 = c_float(t3)
    cnt = c_int(0)

    carr = C_grt_get_trap_wave(dt, byref(ct1), byref(ct2), byref(ct3), byref(cnt))
    arr = npct.as_array(carr, shape=(cnt.value,)).copy()

    C_grt_free(carr)

    return arr


def gen_ricker_wave(f0:float, dt:float):
    '''
        generate Ricker wavelet

        :param    f0:      center frequency (Hz)
        :param    dt:      time interval (s)

        :return:
            - **wave** -    amplitude sequence
    '''
    cnt = c_int(0)

    carr = C_grt_get_ricker_wave(dt, f0, byref(cnt))
    if cast(carr, c_void_p).value is None:
        raise ValueError("NULL pointer")
    arr = npct.as_array(carr, shape=(cnt.value,)).copy()

    C_grt_free(carr)

    return arr