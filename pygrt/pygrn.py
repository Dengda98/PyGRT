"""
    :file:     pygrn.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包括 Python端使用的格林函数 :class:`pygrt.pygrn.PyGreenFunction`

"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import numpy.ctypeslib as npct
from obspy import read, Stream, Trace, UTCDateTime
from obspy.io.sac import SACTrace
from scipy.fft import irfft
from typing import List, Dict



from ctypes import *
from .c_interfaces import *
from .c_structures import *


__all__ = [
    "PyGreenFunction",
]


class PyGreenFunction:
    def __init__(
            self, 
            name:str,
            nt:int, 
            dt:float, 
            freqs:np.ndarray,
            wI:float,
            dist:float,
            depsrc:float, 
            deprcv:float):
        ''' 
            Python端使用的格林函数类

            :param    name:          格林函数名称，震源类型(EX,VF,HF,DD,DS,SS)+三分量(Z,R,T)
            :param    nt:            时间点数  
            :param    dt:            采样间隔(s)  
            :param    freqs:         频率数组(Hz)
            :param    wI:          定义虚频率，omega = w - j*wI, wI = wI  
            :param    dist:          震中距(km)
            :param    depsrc:        震源深度(km)
            :param    deprcv:        台站深度(km)

        '''
        
        # 频率点
        self.freqs = freqs  # 未copy，共享内存  
        self.freqs.flags.writeable = False  # 不允许修改内部值  

        self.name = name
        self.nt = nt
        self.dt = dt 
        self.wI = wI 
        self.dist = dist 
        self.depsrc = depsrc
        self.deprcv = deprcv
        
        nf = len(self.freqs)
        
        # 频谱numpy数据 
        self.cmplx_grn = np.zeros((nf,), dtype=NPCT_CMPLX_TYPE)

        # 虚频率 
        self.wI = wI

        # 提前建立Trace时间序列  
        self.SACTrace = SACTrace(npts=nt, delta=dt, iztype='io') 
        sac = self.SACTrace
        sac.evdp = depsrc
        sac.stel = (-1)*deprcv
        sac.dist = dist
        sac.user0 = wI  # 记录虚频率
        sac.kstnm = 'SYN'
        sac.kcmpnm = name
    

    def plot_response(self):
        '''
            绘制频率响应图，包括幅度响应和相位响应
        '''
        amp = np.abs(self.cmplx_grn)
        phi = np.angle(self.cmplx_grn)

        freqs = self.freqs 

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw=dict(hspace=0.5)) 
        ax1.plot(freqs, amp, 'k', lw=0.6)
        ax1.set_xlabel("Frequency (Hz)")
        ax1.set_ylabel("Amplitude") 
        ax1.grid()


        ax2.plot(freqs, phi, 'k', lw=0.6)
        ax2.set_xlabel("Frequency (Hz)")
        ax2.set_ylabel("Phase") 
        ax2.set_yticks([-np.pi, 0, np.pi], ['$-\pi$', '$0$', '$\pi$'])
        ax2.grid()

        return fig, (ax1, ax2)

        
    def freq2time(self, T0:float, travtP:float, travtS:float, mult:float=1.0):
        '''
            将格林函数从频域转为时域，以 :class:`obspy.Trace` 的形式返回  

            :param    T0:      时域信号的起始时刻相对发震时刻的偏移量(s)，例如T0=5表示发震后5s开始记录波形 

            :return:
                - **tr**:      :class:`obspy.Trace` 类型的格林函数时间序列  
        '''

        self.cmplx_grn[:] *= mult

        freqs = self.freqs

        df = freqs[-1] - freqs[-2]
        sac = self.SACTrace
        nt = sac.npts
        dt = sac.delta
        wI = sac.user0

        T = nt*dt
        if not np.isclose(T*df, 1.0):
            raise ValueError(f"{sac.kcmpnm} length of window not match the freq interval.") 
        
        omegas = 2*np.pi*freqs

        cmlx_grn = self.cmplx_grn * np.exp(1j*omegas*T0)  # 时移

        # 实序列的傅里叶变换 
        data = irfft(cmlx_grn, nt, norm='backward') * (1/dt)  # *(1/dt)和连续傅里叶变换幅值保持一致
        # 抵消虚频率的影响
        data *= np.exp((np.arange(0,nt)*dt + T0)*wI)

        # 保存sac头段变量
        sac.o = 0.0
        sac.b = T0
        # 记录走时
        sac.kt0 = 'P'
        sac.t0 = travtP
        sac.kt1 = 'S'
        sac.t1 = travtS
        # 记录时域数据
        tr = sac.to_obspy_trace()
        tr.data = data
        

        return tr


    
    