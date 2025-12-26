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
            upsampling_n:int,
            freqs:np.ndarray,
            wI:float,
            dist:float,
            depsrc:float, 
            deprcv:float):
        ''' 
            :param    name:          source-type (EX,VF,HF,DD,DS,SS) + component (Z,R,T)
            :param    nt:            number of time points
            :param    dt:            time interval (s)  
            :param    upsampling_n:  upsampling factor 
            :param    freqs:         frequency array (Hz)
            :param    wI:            imaginary angular frequency `wI`，omega = w - j*wI
            :param    dist:          epicentral distance (km)
            :param    depsrc:        source depth (km)
            :param    deprcv:        receiver depth (km)
        '''
        
        # 频率点
        self.freqs = freqs  # 未copy，共享内存  
        self.freqs.flags.writeable = False  # 不允许修改内部值  

        self.name = name
        self.nt = nt
        self.dt = dt 
        self.upsampling_n = upsampling_n 
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
        self.SACTrace = SACTrace(npts=nt*upsampling_n, delta=dt/upsampling_n, iztype='io') 
        sac = self.SACTrace
        sac.evdp = depsrc
        sac.stel = (-1)*deprcv
        sac.dist = dist
        sac.user0 = wI  # 记录虚频率
        sac.kstnm = 'SYN'
        sac.kcmpnm = name
    

    def plot_response(self):
        '''
            plot the frequency response, including amplitude response and phase response
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
            Convert the Green's function from the frequency domain to the time domain 
            and return it in the form of :class:`obspy.Trace`

            :param    T0:      The offset (secs) of the starttime w.r.t the event origin,
                               for example, T0=5 denotes recording began 5 seconds after the origin.

            :return:
                - **tr**:      :class:`obspy.Trace` Green's function
        '''

        self.cmplx_grn[:] *= mult

        freqs = self.freqs

        df = freqs[-1] - freqs[-2]
        sac = self.SACTrace
        nt = sac.npts   # 可能考虑升采样的点数
        dt = sac.delta  # 可能考虑升采样的采样间隔
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


    
    