"""
    :file:     signals.py
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
    :date:     2024-07-24

    震源时间函数使用面积归一化，另提供按解析峰值为1定标的 Ricker 卷积子波

"""

import numpy as np
import numpy.ctypeslib as npct
from ctypes import byref, cast

from .c_interfaces import *


__all__ = [
    "gen_triangle_wave",
    "gen_parabola_wave",
    "gen_trap_wave",
    "gen_asymmetric_cosine_wave",
    "gen_ricker_wave",
]


def gen_triangle_wave(vlen, dt):
    '''
        generate triangle-shape wave

        :param    vlen:    signal length (s)
        :param    dt:      time interval (s)

        :return:
            - **wave** -    amplitude sequence (float64)
    '''
    return gen_trap_wave(vlen/2.0, 0.0, vlen/2.0, dt)


def gen_parabola_wave(vlen, dt):
    '''
        generate parabola-shape wave

        :param    vlen:    signal length (s)
        :param    dt:      time interval (s)

        :return:
            - **wave** -    amplitude sequence (float64)
    '''
    ct1 = REAL(vlen)
    cnt = c_int(0)

    carr = C_grt_get_parabola_wave(dt, byref(ct1), byref(cnt))
    arr = npct.as_array(carr, shape=(cnt.value,)).copy()

    C_grt_free(carr)

    return arr


def gen_trap_wave(t1, t2, t3, dt):
    '''
        Generate an area-normalized trapezoid-shape wave.
        Durations are nonnegative and their sum must be positive.
        Each nonzero duration is rounded up to the sampling grid,
        with at least one sampling interval; zero durations stay zero.

        :param    t1:      Rise duration (s).
        :param    t2:      Plateau duration (s).
        :param    t3:      Fall duration (s).
        :param    dt:      time interval (s)

        :return:
            - **wave** -    amplitude sequence (float64)
    '''
    ct1 = REAL(t1)
    ct2 = REAL(t2)
    ct3 = REAL(t3)
    cnt = c_int(0)

    carr = C_grt_get_trap_wave(dt, byref(ct1), byref(ct2), byref(ct3), byref(cnt))
    arr = npct.as_array(carr, shape=(cnt.value,)).copy()

    C_grt_free(carr)

    return arr


def gen_asymmetric_cosine_wave(t1, t2, dt):
    '''
        Generate an area-normalized asymmetric cosine slip-rate function.
        Rise and fall durations are individually rounded up to the sampling grid,
        with at least one sampling interval for each branch.

        :param    t1:      Rise duration (s), positive.
        :param    t2:      Fall duration (s), positive.
        :param    dt:      Sampling interval (s), positive.

        :return:
            - **wave** -    Amplitude sequence (float64), with ``sum(wave)*dt = 1``.
    '''
    ct1 = REAL(t1)
    ct2 = REAL(t2)
    cnt = c_int(0)

    carr = C_grt_get_asymmetric_cosine_wave(dt, byref(ct1), byref(ct2), byref(cnt))
    arr = npct.as_array(carr, shape=(cnt.value,)).copy()

    C_grt_free(carr)

    return arr


def gen_ricker_wave(f0:float, dt:float):
    '''
        Generate a signed Ricker convolution wavelet with analytic peak amplitude 1.
        This wavelet is not an area-normalized unit-slip source process.

        :param    f0:      center frequency (Hz)
        :param    dt:      time interval (s)

        :return:
            - **wave** -    amplitude sequence (float64)
    '''
    cnt = c_int(0)

    carr = C_grt_get_ricker_wave(dt, f0, byref(cnt))
    if cast(carr, c_void_p).value is None:
        raise ValueError("NULL pointer")
    arr = npct.as_array(carr, shape=(cnt.value,)).copy()

    C_grt_free(carr)

    return arr
