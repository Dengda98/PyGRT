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
        生成三角信号  

        :param    vlen:    信号时长(s)  
        :param    dt:      采样间隔(s)   

        :return: 
            - **wave** -    波形幅值序列
    '''
    return gen_trap_wave(vlen/2.0, vlen/2.0, vlen, dt)


def gen_parabola_wave(vlen, dt):
    '''
        生成抛物线信号  

        :param    vlen:    信号时长(s)  
        :param    dt:      采样间隔(s)   
        
        :return: 
            - **wave** -    波形幅值序列
    '''
    ct1 = c_float(vlen)
    cnt = c_int(0)

    carr = C_get_parabola_wave(dt, byref(ct1), byref(cnt))
    arr = npct.as_array(carr, shape=(cnt.value,)).copy()

    C_free(carr)

    return arr

def gen_trap_wave(t1, t2, t3, dt):
    '''
        生成梯形信号  

        :param    t1:      上坡截止时刻(s)  
        :param    t2:      平台截止时刻(s)  
        :param    t3:      下坡截止时刻(s)  
        :param    dt:      采样间隔(s)   

        :return: 
            - **wave** -    波形幅值序列
    '''
    ct1 = c_float(t1)
    ct2 = c_float(t2)
    ct3 = c_float(t3)
    cnt = c_int(0)

    carr = C_get_trap_wave(dt, byref(ct1), byref(ct2), byref(ct3), byref(cnt))
    arr = npct.as_array(carr, shape=(cnt.value,)).copy()

    C_free(carr)

    return arr


def gen_ricker_wave(f0:float, dt:float):
    ''' 
        生成Ricker子波   

        :param    f0:      中心频率(Hz)  
        :param    dt:      采样间隔(s)   

        :return: 
            - **wave** -    波形幅值序列
    '''
    cnt = c_int(0)

    carr = C_get_ricker_wave(dt, f0, byref(cnt))
    if cast(carr, c_void_p).value is None:
        raise ValueError("NULL pointer")
    arr = npct.as_array(carr, shape=(cnt.value,)).copy()

    C_free(carr)

    return arr