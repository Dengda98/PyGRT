"""
    :file:     signals.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包括一些常见的时间信号，最高幅值均为1    


"""


import numpy as np  

__all__ = [
    "gen_triangle_wave",
    "gen_square_wave",
    "gen_parabola_wave",
    "gen_dirac_wave",
    "gen_step_wave",
    "gen_ramp_wave",
    "gen_smooth_ramp_wave",
    "gen_ricker_wave",
]

def gen_triangle_wave(nlen, vlen, dt):
    '''
        生成三角信号  

        :param    nlen:    总时长(s)  
        :param    vlen:    信号时长(s)  
        :param    dt:      采样间隔(s)   

        :return: 
            - **t** -       时间序列 
            - **wave** -    波形幅值序列
    '''
    t = np.arange(0, nlen+dt, dt)
    triangle_wave = np.zeros_like(t)
    
    vnt = int(vlen/dt)+1
    phase = np.linspace(0, 1, vnt)
    triangle_wave[0:vnt] = (1 - 2 * np.abs(phase - 0.5))
    
    return t, triangle_wave

def gen_square_wave(nlen, vlen, dt):
    '''
        生成矩形信号  
        
        :param    nlen:    总时长(s)  
        :param    vlen:    信号时长(s)  
        :param    dt:      采样间隔(s)   

        :return: 
            - **t** -       时间序列 
            - **wave** -    波形幅值序列
    '''
    t = np.arange(0, nlen+dt, dt)
    square_wave = np.zeros_like(t)
    
    square_wave[0:int(vlen/dt)+1] = 1
    
    return t, square_wave

def gen_parabola_wave(nlen, vlen, dt):
    '''
        生成抛物线信号  

        :param    nlen:    总时长(s)  
        :param    vlen:    信号时长(s)  
        :param    dt:      采样间隔(s)   
        
        :return: 
            - **t** -       时间序列 
            - **wave** -    波形幅值序列
    '''
    t = np.arange(0, nlen+dt, dt)
    parabola_wave = np.zeros_like(t)
    
    vnt = int(vlen/dt)+1
    phase = np.linspace(0, 1, vnt)
    parabola_wave[0:vnt] = 1 - (phase - 0.5) ** 2 * 4
    
    return t, parabola_wave

def gen_dirac_wave(nlen, dt):
    '''
        生成脉冲信号      

        :param    nlen:    总时长(s)  
        :param    dt:      采样间隔(s)   

        :return: 
            - **t** -       时间序列 
            - **wave** -    波形幅值序列
    '''
    t = np.arange(0, nlen+dt, dt)
    wave = np.zeros_like(t) 
    wave[0] = 1
    return t, wave

def gen_step_wave(nlen, dt):
    '''
        生成阶跃信号(全1信号)   

        :param    nlen:    总时长(s)  
        :param    dt:      采样间隔(s)   

        :return: 
            - **t** -       时间序列 
            - **wave** -    波形幅值序列
    '''
    t = np.arange(0, nlen+dt, dt)
    step_wave = np.ones_like(t) 
    return t, step_wave

def gen_ramp_wave(nlen, t0, dt):
    '''
        生成斜坡信号  

        :param    nlen:    总时长(s)  
        :param    t0:      上坡时长(s)  
        :param    dt:      采样间隔(s)   

        :return: 
            - **t** -       时间序列 
            - **wave** -    波形幅值序列
    '''
    t = np.arange(0, nlen+dt, dt)
    ramp_wave = np.ones_like(t) 
    
    vnt = int(t0/dt)+1
    phase = np.linspace(0, 1, vnt)
    ramp_wave[:vnt] = phase
    return t, ramp_wave

def gen_smooth_ramp_wave(nlen, t0, dt):
    '''
        生成光滑斜坡信号  

        :param    nlen:    总时长(s)  
        :param    t0:      上坡(1-1/e)时长(s)  
        :param    dt:      采样间隔(s)  

        :return: 
            - **t** -       时间序列 
            - **wave** -    波形幅值序列
    '''
    t = np.arange(0, nlen+dt, dt)
    wave = (1 - np.exp(-t/t0))
    return t, wave

def gen_trap_wave(nlen, t1, t2, t3, dt):
    '''
        生成梯形信号  

        :param    nlen:    总时长(s)  
        :param    t1:      上坡截止时刻(s)  
        :param    t2:      平台截止时刻(s)  
        :param    t3:      下坡截止时刻(s)  
        :param    dt:      采样间隔(s)   

        :return: 
            - **t** -       时间序列 
            - **wave** -    波形幅值序列
    '''
    t = np.arange(0, nlen+dt, dt)
    wave = np.piecewise(
        t,
        [(t > 0 )* (t < t1), (t >= t1) * (t < t2), (t >= t2) * (t < t3), t >= t3],
        [lambda t: t / t1,
         1,
         lambda t: (t3 - t)/(t3 - t2), 
         0]
    )
    return t, wave

def gen_ricker_wave(nlen, f0:float, dt:float):
    ''' 
        生成Ricker子波   

        :param    nlen:    总时长(s)   
        :param    f0:      中心频率(Hz)  
        :param    dt:      采样间隔(s)   

        :return: 
            - **t** -       时间序列 
            - **wave** -    波形幅值序列
    '''

    t = np.arange(0, nlen+dt, dt)
    t0 = 1.0/f0
    a = np.pi**2 * f0**2 * (t-t0)**2
    wave = (1 - 2*a ) * np.exp(-a)

    return t, wave
