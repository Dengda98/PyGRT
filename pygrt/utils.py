"""
    :file:     utils.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包含一些数据处理操作上的补充:   
        1、双力偶源、单力源、爆炸源、矩张量源 通过格林函数合成理论地震图的函数\n
        2、Stream类型的时域卷积、微分、积分 (基于numpy和scipy)    \n
        3、Stream类型写到本地sac文件，自定义名称    \n
        4、读取波数积分和峰谷平均法过程文件  \n


"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from obspy import Stream, Trace 
from obspy.core import AttribDict
from copy import deepcopy
from scipy.signal import oaconvolve
from scipy.fft import rfft, irfft
import math 
import os

from .c_structures import *


__all__ = [
    "gen_syn_from_gf_DC",
    "gen_syn_from_gf_SF",
    "gen_syn_from_gf_EXP",
    "gen_syn_from_gf_MT",

    "stream_convolve",
    "stream_integral",
    "stream_diff",
    "stream_write_sac",

    "read_statsfile",
    "read_statsfile_ptam",
    "plot_statsdata",
    "plot_statsdata_ptam",
]


def gen_syn_from_gf_DC(st:Stream, M0:float, strike:float, dip:float, rake:float, az:float):
    '''
        双力偶源，角度单位均为度   

        :param    st:       计算好的时域格林函数, :class:`obspy.Stream` 类型
        :param    M0:       标量地震矩, 单位dyne*cm
        :param    strike:   走向，以北顺时针为正，0<=strike<=360
        :param    dip:      倾角，以水平面往下为正，0<=dip<=90
        :param    rake:     滑动角，在断层面相对于走向方向逆时针为正，-180<=rake<=180
        :param    az:       台站方位角，以北顺时针为正，0<=az<=360

        :return:
            - **stream** -  三分量ZRT地震图, :class:`obspy.Stream` 类型
    '''
    # 检查通道 
    all_sets = {'DDZ', 'DDR', 'DSZ', 'DSR', 'DST', 'SSZ', 'SSR', 'SST'}
    chs = [tr.stats.channel for tr in st]
    if not all_sets.issubset(chs):
        raise ValueError("Failed, No green functions of Double Couple source.")


    az_deg = az
    baz_deg = 180 + az
    if baz_deg > 360:
        baz_deg -= 360

    strike = math.radians(strike)
    dip    = math.radians(dip)
    rake   = math.radians(rake)
    az     = math.radians(az)
    theta = az - strike
    srak = math.sin(rake);   crak = math.cos(rake)
    sdip = math.sin(dip);    cdip = math.cos(dip)
    sdip2 = 2*sdip*cdip;     cdip2 = 2*cdip*cdip - 1
    sthe = math.sin(theta);  cthe = math.cos(theta)
    sthe2 = 2*sthe*cthe;     cthe2 = 2*cthe*cthe - 1

    M0 *= 1e-20

    # 公式(4.8.35)
    A0 = M0 * (0.5*sdip2*srak)
    A1 = M0 * (cdip*crak*cthe - cdip2*srak*sthe)
    A2 = M0 * (0.5*sdip2*srak*cthe2 + sdip*crak*sthe2)
    A4 = M0 * (- cdip2*srak*cthe - cdip*crak*sthe)
    A5 = M0 * (sdip*crak*cthe2 - 0.5*sdip2*srak*sthe2)

    trZ = st.select(channel='DDZ')[0].copy()
    trZ.stats.channel = 'DCZ'
    __check_trace_attr_sac(trZ, az=az_deg, baz=baz_deg)
    trZ.data = A0*st.select(channel='DDZ')[0].data + \
               A1*st.select(channel='DSZ')[0].data + \
               A2*st.select(channel='SSZ')[0].data
    
    trR = st.select(channel='DDR')[0].copy()
    trR.stats.channel = 'DCR'
    __check_trace_attr_sac(trR, az=az_deg, baz=baz_deg)
    trR.data = A0*st.select(channel='DDR')[0].data + \
               A1*st.select(channel='DSR')[0].data + \
               A2*st.select(channel='SSR')[0].data
    
    trT = st.select(channel='DST')[0].copy()
    trT.stats.channel = 'DCT'
    __check_trace_attr_sac(trT, az=az_deg, baz=baz_deg)
    trT.data = A4*st.select(channel='DST')[0].data + \
               A5*st.select(channel='SST')[0].data
    
    return Stream([trR, trT, trZ])


def gen_syn_from_gf_SF(st:Stream, S:float, fN:float, fE:float, fZ:float, az:float):
    '''
        单力源，力分量单位均为dyne   

        :param    st:    计算好的时域格林函数, :class:`obspy.Stream` 类型
        :param     S:    力的放大系数
        :param    fN:    北向力，向北为正  
        :param    fE:    东向力，向东为正  
        :param    fZ:    垂向力，向下为正  
        :param    az:    台站方位角，以北顺时针为正，0<=az<=360   

        :return:
            - **stream** - 三分量ZRT地震图, :class:`obspy.Stream` 类型
    '''
    # 检查通道 
    all_sets = {'VFZ', 'VFR', 'HFZ', 'HFR', 'HFT'}
    chs = [tr.stats.channel for tr in st]
    if not all_sets.issubset(chs):
        raise ValueError("Failed, No green functions of Single Force source.")
    
    az_deg = az
    baz_deg = 180 + az
    if baz_deg > 360:
        baz_deg -= 360

    az  = math.radians(az)
    saz = math.sin(az); caz = math.cos(az)

    S *= 1e-15
    fN *= S
    fE *= S
    fZ *= S

    # 公式(4.6.20)
    trZ = st.select(channel='VFZ')[0].copy()
    trZ.stats.channel = 'SFZ'
    __check_trace_attr_sac(trZ, az=az_deg, baz=baz_deg)
    trZ.data = (fN*caz + fE*saz)*st.select(channel='HFZ')[0].data + fZ*st.select(channel='VFZ')[0].data

    trR = st.select(channel='VFR')[0].copy()
    trR.stats.channel = 'SFR'
    __check_trace_attr_sac(trR, az=az_deg, baz=baz_deg)
    trR.data = (fN*caz + fE*saz)*st.select(channel='HFR')[0].data + fZ*st.select(channel='VFR')[0].data

    trT = st.select(channel='HFT')[0].copy()
    trT.stats.channel = 'SFT'
    __check_trace_attr_sac(trT, az=az_deg, baz=baz_deg)
    trT.data = (- fN*saz + fE*caz)*st.select(channel='HFT')[0].data

    return Stream([trR, trT, trZ])


def gen_syn_from_gf_EXP(st:Stream, M0:float, az:float):
    '''
        爆炸源

        :param    st:          计算好的时域格林函数, :class:`obspy.Stream` 类型
        :param    M0:          标量地震矩, 单位dyne*cm  
        :param    az:          台站方位角，以北顺时针为正，0<=az<=360 [不用于计算]  

        :return:
            - **stream** -       三分量ZRT地震图, :class:`obspy.Stream` 类型
    '''
    # 检查通道 
    all_sets = {'EXZ', 'EXR'}
    chs = [tr.stats.channel for tr in st]
    if not all_sets.issubset(chs):
        raise ValueError("Failed, No green functions of Explosion source.")

    az_deg = az
    baz_deg = 180 + az
    if baz_deg > 360:
        baz_deg -= 360

    M0 *= 1e-20

    trZ = st.select(channel='EXZ')[0].copy()
    trZ.stats.channel = 'EXZ'
    __check_trace_attr_sac(trZ, az=az_deg, baz=baz_deg)
    trZ.data = M0 * st.select(channel='EXZ')[0].data

    trR = st.select(channel='EXR')[0].copy()
    trR.stats.channel = 'EXR'
    __check_trace_attr_sac(trR, az=az_deg, baz=baz_deg)
    trR.data = M0 * st.select(channel='EXR')[0].data

    return Stream([trR, trZ])


def gen_syn_from_gf_MT(st:Stream, M0:float, MT:np.ndarray, az:float):
    ''' 
        矩张量源，单位为dyne*cm

        :param    st:          计算好的时域格林函数, :class:`obspy.Stream` 类型
        :param    M0:          标量地震矩
        :param    MT:          矩张量 (M11, M12, M13, M22, M23, M33),下标1,2,3分别代表北向，东向，垂直向上  
        :param    az:          台站方位角，以北顺时针为正，0<=az<=360   

        :return:
            - **stream** -       三分量ZRT地震图, :class:`obspy.Stream` 类型
    '''
    # 检查通道 
    all_sets = {'EXZ', 'EXR', 'DDZ', 'DDR', 'DSZ', 'DSR', 'DST', 'SSZ', 'SSR', 'SST'}
    chs = [tr.stats.channel for tr in st]
    if not all_sets.issubset(chs):
        raise ValueError("Failed, No green functions of Explosion or Double Couple source.")

    az_deg = az
    baz_deg = 180 + az
    if baz_deg > 360:
        baz_deg -= 360

    MT = np.array(MT)

    # 不使用*=，否则会原地修改MT
    MT = MT*M0

    # 去除各向同性分量，作为爆炸源
    Mexp = (MT[0] + MT[3] + MT[5])/3
    trR, trZ = gen_syn_from_gf_EXP(st, Mexp, az) # Mexp在函数内部会乘1e-20
    trT = trR.copy()
    trT.data[:] = 0.0

    Mexp *= 1e-20
    MT = MT*1e-20
    M11, M12, M13, M22, M23, M33 = MT
    M11 -= Mexp
    M22 -= Mexp
    M33 -= Mexp

    az  = math.radians(az)
    saz = math.sin(az); caz  = math.cos(az)
    saz2 = 2*saz*caz;   caz2 = 2*caz*caz - 1

    # 公式(4.9.7)，但去除了各向同性的量
    A0 = ((2*M33 - M11 - M22)/6 )
    A1 = (- (M13*caz + M23*saz))
    A2 = (0.5*(M11 - M22)*caz2+ M12*saz2)
    A4 = (M13*saz - M23*caz)
    A5 = (-0.5*(M11 - M22)*saz2 + M12*caz2)

    trZ.stats.channel = 'MTZ'
    __check_trace_attr_sac(trZ, az=az_deg, baz=baz_deg)
    trZ.data += A0*st.select(channel='DDZ')[0].data + \
                A1*st.select(channel='DSZ')[0].data + \
                A2*st.select(channel='SSZ')[0].data
    
    trR.stats.channel = 'MTR'
    __check_trace_attr_sac(trR, az=az_deg, baz=baz_deg)
    trR.data += A0*st.select(channel='DDR')[0].data + \
                A1*st.select(channel='DSR')[0].data + \
                A2*st.select(channel='SSR')[0].data
    
    trT.stats.channel = 'MTT'
    __check_trace_attr_sac(trT, az=az_deg, baz=baz_deg)
    trT.data += A4*st.select(channel='DST')[0].data + \
                A5*st.select(channel='SST')[0].data  
    
    return Stream([trR, trT, trZ]) 
    

def __check_trace_attr_sac(tr:Trace, **kwargs):
    '''
        临时函数，检查trace中是否有sac字典，并将kwargs内容填入  
    '''
    if hasattr(tr.stats, 'sac'):
        for k, v in kwargs.items():
            tr.stats.sac[k] = v 
    else: 
        tr.stats.sac = AttribDict(**kwargs)


# def stream_convolve(st0:Stream, timearr:np.ndarray, inplace=True):
#     '''
#         频域实现线性卷
#     '''
#     st = st0 if inplace else deepcopy(st0)
    
#     sacAttr = st[0].stats.sac
#     try:
#         wI = sacAttr['user0']  # 虚频率
#     except:
#         wI = 0.0
#     nt = sacAttr['npts']
#     dt = sacAttr['delta']

#     nt2 = len(timearr)
#     N = nt+nt2-1
#     nf = N//2 + 1

#     wI_exp1 = np.exp(-wI*np.arange(0,nt)*dt)
#     wI_exp2 = np.exp( wI*np.arange(0,nt)*dt)

#     fft_tf = np.ones((nf, ), dtype='c16') 
#     # if scale is None:
#     #     scale = 1.0/np.trapz(timearr, dx=dt)

#     timearr0 = timearr.copy()
#     timearr0.resize((N,))  # 填充0
#     timearr0[:nt] *= wI_exp1
#     # FFT 
#     fft_tf[:] = rfft(timearr0, N)
#     fft_tf[:] *= dt
    
#     # 对每一道做相同处理
#     for tr in st:
#         data = tr.data
#         # 虚频率 
#         data[:] *= wI_exp1

#         # FFT
#         fft_d = rfft(data, N)

#         # 卷积+系数
#         fft_d[:] *= fft_tf

#         # IFFT
#         data[:] = irfft(fft_d, N)[:nt]

#         # 虚频率 
#         data[:] *= wI_exp2

#     return st


def stream_convolve(st0:Stream, signal:np.ndarray, inplace=True):
    '''
        对stream中每一道信号做线性卷积  

        :param    st0:        记录多个Trace的 :class:`obspy.Stream` 类型
        :param    signal:     卷积信号
        :param    inplace:    是否做原地更改  

        :return:
            - **stream** -    处理后的结果, :class:`obspy.Stream` 类型
    '''
    st = st0 if inplace else deepcopy(st0)
    
    for tr in st:
        data = tr.data 
        data[:] = oaconvolve(data, signal, mode='full')[:data.shape[0]]

    return st


def stream_integral(st0:Stream, inplace=True):
    '''
        对stream中每一道信号做梯形积分
        
        :param    st0:        记录多个Trace的 :class:`obspy.Stream` 类型
        :param    inplace:    是否做原地更改  

        :return:
            - **stream** -    处理后的结果, :class:`obspy.Stream` 类型
    '''
    st = st0 if inplace else deepcopy(st0)
    for tr in st:
        dt = tr.stats.delta
        data = tr.data 
        lastx = data[0]
        data[0] = 0.0
        
        for i in range(1, len(data)):
            tmp = data[i]
            data[i] = 0.5*(data[i] + lastx)*dt + data[i-1]
            lastx = tmp

    return st


def stream_diff(st0:Stream, inplace=True):
    '''
        对stream中每一道信号做中心差分

        :param    st0:        记录多个Trace的 :class:`obspy.Stream` 类型
        :param    inplace:    是否做原地更改  

        :return:
            - **stream** -    处理后的结果, :class:`obspy.Stream` 类型
    '''
    st = st0 if inplace else deepcopy(st0)
    
    for tr in st:
        data = tr.data 
        data[:] = np.gradient(data, tr.stats.delta)

    return st


def stream_write_sac(st:Stream, path:str):
    '''
        将一系列Trace以SAC形式保存到本地，以发震时刻作为参考0时刻

        :param    st:         记录多个Trace的 :class:`obspy.Stream` 类型
        :param    path:       保存文件名，不要加后缀  

    '''
    # 新建对应文件夹
    parent_dir = os.path.dirname(path)
    if(len(parent_dir)>0):
        os.makedirs(parent_dir, exist_ok=True)

    # 每一道的保存路径为path+{channel}
    for tr in st:
        tr.write(".".join([path, tr.stats.channel]) + '.sac', format='SAC')







def read_statsfile(statsfile:str):
    '''
        读取单个频率下波数积分(或Filon积分)的记录文件  

        :param    statsfile:       文件路径  

        :return:
            - **data** -     `numpy.ndarray <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>`_ 自定义类型数组 
    '''
    data = np.fromfile(statsfile, 
        dtype=[
            ('k', NPCT_REAL_TYPE), 

            # 核函数 F(k, w)
            ('EXP_q0', NPCT_CMPLX_TYPE),
            ('EXP_w0', NPCT_CMPLX_TYPE),

            ('VF_q0', NPCT_CMPLX_TYPE),
            ('VF_w0', NPCT_CMPLX_TYPE),

            ('HF_q1', NPCT_CMPLX_TYPE),
            ('HF_w1', NPCT_CMPLX_TYPE),
            ('HF_v1', NPCT_CMPLX_TYPE),

            ('DC_q0', NPCT_CMPLX_TYPE),
            ('DC_w0', NPCT_CMPLX_TYPE),

            ('DC_q1', NPCT_CMPLX_TYPE),
            ('DC_w1', NPCT_CMPLX_TYPE),
            ('DC_v1', NPCT_CMPLX_TYPE),

            ('DC_q2', NPCT_CMPLX_TYPE),
            ('DC_w2', NPCT_CMPLX_TYPE),
            ('DC_v2', NPCT_CMPLX_TYPE),

            # ===============================
            # 被积函数 F(k, w)*Jm(kr)*k
            # 下划线后的两位数字，第一位m代表阶数，但不完全对应Bessel函数的阶数，
            # 因为存在Bessel函数的导数的情况，需要使用递推公式，故存在不对应情况；
            # 第二位p代表自定义的被积函数类型，基于式(5.6.22)，分别为  
            #       if m==0:
            #           if p==0:  - q0 * J1(kr) * k 
            #           if p==2:    w0 * J0(kr) * k
            #       else if m==1 or m==2:
            #           if p==0:    qm * Jm-1(kr) * k 
            #           if p==1:  - (qm + vm) * Jm(kr) * k / (kr)
            #           if p==2:    wm * Jm(kr) * k 
            #           if p==3:  - vm * Jm-1(kr) * k
            #               
            ('EXP_00', NPCT_CMPLX_TYPE),
            ('EXP_02', NPCT_CMPLX_TYPE),

            ('VF_00', NPCT_CMPLX_TYPE),
            ('VF_02', NPCT_CMPLX_TYPE),

            ('HF_10', NPCT_CMPLX_TYPE),
            ('HF_11', NPCT_CMPLX_TYPE),
            ('HF_12', NPCT_CMPLX_TYPE),
            ('HF_13', NPCT_CMPLX_TYPE),
            ('DC_00', NPCT_CMPLX_TYPE),
            ('DC_02', NPCT_CMPLX_TYPE),
            ('DC_10', NPCT_CMPLX_TYPE),
            ('DC_11', NPCT_CMPLX_TYPE),
            ('DC_12', NPCT_CMPLX_TYPE),
            ('DC_13', NPCT_CMPLX_TYPE),
            ('DC_20', NPCT_CMPLX_TYPE),
            ('DC_21', NPCT_CMPLX_TYPE),
            ('DC_22', NPCT_CMPLX_TYPE),
            ('DC_23', NPCT_CMPLX_TYPE),
        ]
    )

    return data


def read_statsfile_ptam(statsfile:str):
    '''
        读取单个频率下峰谷平均法的记录文件  

        :param    statsfile:       文件路径  

        :return:
            - **data** -     `numpy.ndarray <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>`_ 自定义类型数组 
    '''
    data = np.fromfile(statsfile, 
        dtype=[
            # 各格林函数数值积分的值(k上限位于不同的波峰波谷)
            # 名称中的两位数字的解释和`read_statsfile`函数中的解释相同，
            # 开头的sum表示这是波峰波谷位置处的数值积分的值(不含dk)，
            # 末尾的k表示对应积分值的波峰波谷位置的k值
            # 
            ('sum_EXP_00_k', NPCT_REAL_TYPE),
            ('sum_EXP_00',  NPCT_CMPLX_TYPE),
            ('sum_VF_00_k',  NPCT_REAL_TYPE),
            ('sum_VF_00',   NPCT_CMPLX_TYPE),
            ('sum_DC_00_k',  NPCT_REAL_TYPE),
            ('sum_DC_00',   NPCT_CMPLX_TYPE),

            ('sum_EXP_02_k', NPCT_REAL_TYPE),
            ('sum_EXP_02',  NPCT_CMPLX_TYPE),
            ('sum_VF_02_k',  NPCT_REAL_TYPE),
            ('sum_VF_02',   NPCT_CMPLX_TYPE),
            ('sum_DC_02_k',  NPCT_REAL_TYPE),
            ('sum_DC_02',   NPCT_CMPLX_TYPE),
            

            ('sum_HF_10_k',  NPCT_REAL_TYPE),
            ('sum_HF_10',   NPCT_CMPLX_TYPE),
            ('sum_DC_10_k',  NPCT_REAL_TYPE),
            ('sum_DC_10',   NPCT_CMPLX_TYPE),
            ('sum_HF_11_k',  NPCT_REAL_TYPE),
            ('sum_HF_11',   NPCT_CMPLX_TYPE),
            ('sum_DC_11_k',  NPCT_REAL_TYPE),
            ('sum_DC_11',   NPCT_CMPLX_TYPE),
            ('sum_HF_12_k',  NPCT_REAL_TYPE),
            ('sum_HF_12',   NPCT_CMPLX_TYPE),
            ('sum_DC_12_k',  NPCT_REAL_TYPE),
            ('sum_DC_12',   NPCT_CMPLX_TYPE),
            ('sum_HF_13_k',  NPCT_REAL_TYPE),
            ('sum_HF_13',   NPCT_CMPLX_TYPE),
            ('sum_DC_13_k',  NPCT_REAL_TYPE),
            ('sum_DC_13',   NPCT_CMPLX_TYPE),

            ('sum_DC_20_k',  NPCT_REAL_TYPE),
            ('sum_DC_20',   NPCT_CMPLX_TYPE),
            ('sum_DC_21_k',  NPCT_REAL_TYPE),
            ('sum_DC_21',   NPCT_CMPLX_TYPE),
            ('sum_DC_22_k',  NPCT_REAL_TYPE),
            ('sum_DC_22',   NPCT_CMPLX_TYPE),
            ('sum_DC_23_k',  NPCT_REAL_TYPE),
            ('sum_DC_23',   NPCT_CMPLX_TYPE),
        ]
    )

    return data


def plot_statsdata(statsdata:np.ndarray, srctype:str, mtype:str, qwvtype:str, ptype:str, RorI:bool=True):
    r'''
        根据 :func:`read_statsfile <pygrt.utils.read_statsfile>` 函数函数读取的数据，
        绘制核函数 :math:`F(k,\omega)`、被积函数 :math:`F(k,\omega)J_m(kr)k`，以及简单计算累积积分 :math:`\sum F(k,\omega)J_m(kr)k` 并绘制。

        .. note:: 并不是每个震源类型对应的每一阶每种积分类型都存在，详见 :ref:`grn_types`。

        :param    statsdata:         :func:`read_statsfile <pygrt.utils.read_statsfile>` 函数返回值 
        :param    srctype:           震源类型的缩写，包括EXP、VF、HF、DC  
        :param    qwvtype:           
        :param    mtype:             阶数(0,1,2)，不完全对应公式中Bessel函数的阶数，因为存在Bessel函数的导数，需要使用递推公式
        :param    ptype:             积分类型(0,1,2,3) 
        :param    RorI:              绘制实部还是虚部，默认实部

        :return:  
                - **fig** -                        matplotlib.Figure对象   
                - **(ax1, ax2, ax3)**  -           3个matplotlib.Axes对象的元组
    '''

    mtype = str(mtype)
    ptype = str(ptype)

    karr = statsdata['k'] 
    Fname = f'{srctype}_{qwvtype}{mtype}'
    Farr = statsdata[Fname]   # 核函数 F(k, w)
    FJname = f'{srctype}_{mtype}{ptype}'
    FJarr = statsdata[FJname]  # 被积函数 F(k, w) * Jm(kr) * k

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 9), gridspec_kw=dict(hspace=0.7))
    if RorI:
        ax1.plot(karr, np.real(Farr), 'k', lw=0.8, label='Real') 
    else:
        ax1.plot(karr, np.imag(Farr), 'k', lw=0.8, label='Imag') 

    ax1.set_xlabel('k /$km^{-1}$')
    ax1.set_title(f'F(k, w)   [{Fname}]')
    ax1.grid()
    ax1.legend(loc='lower left')


    if RorI:
        ax2.plot(karr, np.real(FJarr), 'k', lw=0.8, label='Real') 
    else:
        ax2.plot(karr, np.imag(FJarr), 'k', lw=0.8, label='Imag') 
    ax2.set_title(f'F(k, w)*Jm(kr)*k   [{FJname}]')
    ax2.set_xlabel('k /$km^{-1}$')
    ax2.grid()
    ax2.legend(loc='lower left')

    # 数值积分，不乘系数dk 
    # 以下的特殊处理是使用峰谷平均法时计算的中间积分值也会记录到该文件中，
    # 记录波峰波谷值的由另一个文件负责
    # 而常规积分阶段和峰谷平均法阶段使用的dk不一致，故先乘再除
    Parr = np.cumsum(FJarr[:-1] * np.diff(karr)) / (karr[1] - karr[0]) 

    if RorI:
        ax3.plot(karr[:-1], np.real(Parr), 'k', lw=0.8, label='Real') 
    else:
        ax3.plot(karr[:-1], np.imag(Parr), 'k', lw=0.8, label='Imag') 
    ax3.set_title('$\sum_k$ F(k, w)*Jm(kr)*k')
    ax3.set_xlabel("k /$km^{-1}$")
    ax3.grid()
    ax3.legend(loc='lower left')

    return fig, (ax1, ax2, ax3)


def plot_statsdata_ptam(statsdata:np.ndarray, statsdata_ptam:np.ndarray, srctype:str, mtype:str, ptype:str, RorI:bool=True):
    r'''
        根据 :func:`read_statsfile <pygrt.utils.read_statsfile>` 函数以及 
        :func:`read_statsfile_ptam <pygrt.utils.read_statsfile_ptam>` 函数读取的数据，
        简单计算并绘制累积积分以及峰谷平均法使用的波峰波谷的位置。

        .. note:: 并不是每个震源类型对应的每一阶每种积分类型都存在，详见 :ref:`grn_types`。

        :param    statsdata:         :func:`read_statsfile <pygrt.utils.read_statsfile>` 函数返回值 
        :param    statsdata_ptam:    :func:`read_statsfile <pygrt.utils.read_statsfile_ptam>` 函数返回值 
        :param    srctype:           震源类型的缩写，包括EXP、VF、HF、DC  
        :param    mtype:             阶数(0,1,2)，不完全对应公式中Bessel函数的阶数，因为存在Bessel函数的导数，需要使用递推公式
        :param    ptype:             积分类型(0,1,2,3) 
        :param    RorI:              绘制实部还是虚部，默认实部

        :return:  
                - **fig** -           matplotlib.Figure对象   
                - **ax**  -           matplotlib.Axes对象  
    '''

    fig, ax = plt.subplots(1,1, figsize=(8, 5))

    mtype = str(mtype)
    ptype = str(ptype)

    karr = statsdata['k'] 
    FJname = f'{srctype}_{mtype}{ptype}'
    FJarr = statsdata[FJname]  # 被积函数 F(k, w) * Jm(kr) * k
    # 波峰波谷位置，用红十字标记
    ptKarr = statsdata_ptam[f'sum_{srctype}_{mtype}{ptype}_k']
    ptFJarr = statsdata_ptam[f'sum_{srctype}_{mtype}{ptype}']

    # 数值积分，不乘系数dk 
    # 以下的特殊处理是使用峰谷平均法时计算的中间积分值也会记录到该文件中，
    # 记录波峰波谷值的由另一个文件负责
    # 而常规积分阶段和峰谷平均法阶段使用的dk不一致，故先乘再除
    Parr = np.cumsum(FJarr[:-1] * np.diff(karr)) / (karr[1] - karr[0]) 


    if RorI:
        ax.plot(karr[:-1], np.real(Parr), 'k', lw=0.8, label='Real') 
        ax.plot(ptKarr, np.real(ptFJarr), 'r+', markersize=6)
    else:
        ax.plot(karr[:-1], np.imag(Parr), 'k', lw=0.8, label='Imag') 
        ax.plot(ptKarr, np.imag(ptFJarr), 'r+', markersize=6)

    ax.set_title('$\sum_k$ F(k, w)*Jm(kr)*k')
    ax.set_xlabel("k /$km^{-1}$")
    ax.grid()
    ax.legend(loc='lower left')


    return fig, ax