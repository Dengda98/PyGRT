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
from ctypes import *
from typing import List
from copy import deepcopy

from numpy.typing import ArrayLike

from .c_interfaces import *


__all__ = [
    "gen_syn_from_gf_DC",
    "gen_syn_from_gf_SF",
    "gen_syn_from_gf_EXP",
    "gen_syn_from_gf_MT",

    "compute_strain",
    "compute_stress",

    "stream_convolve",
    "stream_integral",
    "stream_diff",
    "stream_write_sac",

    "read_statsfile",
    "read_statsfile_ptam",
    "plot_statsdata",
    "plot_statsdata_ptam",
]


#=================================================================================================================
#
#                                           根据辐射因子合成地震图
#
#=================================================================================================================

def _gen_syn_from_gf(st:Stream, calc_upar:bool, compute_type:str, M0:float, az:float, ZNE=False, **kwargs):
    r"""
        一个发起函数，根据不同震源参数，从格林函数中合成理论地震图

        :param    st:              计算好的时域格林函数, :class:`obspy.Stream` 类型
        :param    calc_upar:       是否计算位移u的空间导数
        :param    compute_type:    计算类型，应为以下之一: 
        :param    M0:              标量地震矩, 单位dyne*cm
        :param    az:              方位角(度)
        :param    ZNE:             是否以ZNE分量输出?
            'COMPUTE_EXP'(爆炸源), 'COMPUTE_SF'(单力源),
            'COMPUTE_DC'(双力偶源), 'COMPUTE_MT'(矩张量源)
    """
    chs = ['Z', 'R', 'T']
    sacin_prefixes = ["", "z", "r", ""]   # 输入通道名
    sacout_prefixes = ["", "z", "r", "t"]   # 输出通道名
    srcName = ["EX", "VF", "HF", "DD", "DS", "SS"]
    allchs = [tr.stats.channel for tr in st]

    s_compute_type = compute_type.split("_")[1][:2]

    baz = 180 + az
    if baz > 360:
        baz -= 360

    azrad = np.deg2rad(az)

    calcUTypes = 4 if calc_upar else 1

    stall = Stream()

    dist = st[0].stats.sac['dist']
    upar_scale:float = 1.0
    for ityp in range(calcUTypes):
        if ityp > 0:
            upar_scale = 1e-5
        if ityp == 3:
            upar_scale /= dist

        srcCoef = _set_source_coef(ityp==3, upar_scale, compute_type, M0, azrad, **kwargs)

        inpref = sacin_prefixes[ityp]
        outpref = sacout_prefixes[ityp]

        for c in range(3):
            ch = chs[c]
            tr:Trace = st[0].copy()
            tr.data[:] = 0.0
            tr.stats.channel = kcmpnm = f'{outpref}{s_compute_type}{ch}'
            __check_trace_attr_sac(tr, az=az, baz=baz, kcmpnm=kcmpnm)
            for k in range(6):
                coef = srcCoef[c, k]
                if coef==0.0:
                    continue

                # 读入数据
                channel = f'{inpref}{srcName[k]}{ch}'
                if channel not in allchs:
                    raise ValueError(f"Failed, channel=\"{channel}\" not exists.")
                    
                tr0 = st.select(channel=channel)[0].copy()
                
                tr.data += coef*tr0.data

            stall.append(tr)

    if ZNE:
        stall = _data_zrt2zne(stall)
            
    return stall


def _gen_syn_from_static_gf(grn:dict, calc_upar:bool, compute_type:str, M0:float, ZNE=False, **kwargs):
    r"""
        一个发起函数，根据不同震源参数，从静态格林函数中合成理论静态场

        :param    grn:             计算好的静态格林函数, 字典类型
        :param    calc_upar:       是否计算位移u的空间导数
        :param    compute_type:    计算类型，应为以下之一: 
        :param    M0:              标量地震矩, 单位dyne*cm
        :param    ZNE:             是否以ZNE分量输出?
            'COMPUTE_EXP'(爆炸源), 'COMPUTE_SF'(单力源),
            'COMPUTE_DC'(双力偶源), 'COMPUTE_MT'(矩张量源)
    """
    chs = ['Z', 'R', 'T']
    sacin_prefixes = ["", "z", "r", ""]   # 输入通道名
    srcName = ["EX", "VF", "HF", "DD", "DS", "SS"]
    allchs = list(grn.keys())

    s_compute_type = compute_type.split("_")[1][:2]

    calcUTypes = 4 if calc_upar else 1

    xarr:np.ndarray = grn['_xarr']
    yarr:np.ndarray = grn['_yarr']

    # 结果字典
    resDct = {}

    # 基本数据拷贝
    for k in grn.keys():
        if k[0] != '_':
            continue 
        resDct[k] = deepcopy(grn[k])

    XX = np.zeros((calcUTypes, 3, len(xarr), len(yarr)), dtype='f8')
    dblsyn = (c_double*3)()
    dblupar = (c_double*9)()

    for iy in range(len(yarr)):
        for ix in range(len(xarr)):
            # 震中距
            dist = max(np.sqrt(xarr[ix]**2 + yarr[iy]**2), 1e-5)

            # 方位角
            azrad = np.arctan2(yarr[iy], xarr[ix])

            upar_scale:float = 1.0
            for ityp in range(calcUTypes):
                if ityp > 0:
                    upar_scale = 1e-5
                if ityp == 3:
                    upar_scale /= dist

                srcCoef = _set_source_coef(ityp==3, upar_scale, compute_type, M0, azrad, **kwargs)

                inpref = sacin_prefixes[ityp]

                for c in range(3):
                    ch = chs[c]
                    for k in range(6):
                        coef = srcCoef[c, k]
                        if coef==0.0:
                            continue

                        # 读入数据
                        channel = f'{inpref}{srcName[k]}{ch}'
                        if channel not in allchs:
                            raise ValueError(f"Failed, channel=\"{channel}\" not exists.")
                            
                        XX[ityp, c, ix, iy] += coef*grn[channel][ix, iy]

            if ZNE:
                for i in range(3):
                    dblsyn[i] = XX[0, i, ix, iy]
                    if calc_upar:
                        for k in range(3):
                            dblupar[k + i*3] = XX[i+1, k, ix, iy]
                if calc_upar:
                    C_rot_zrt2zxy_upar(azrad, dblsyn, dblupar, dist*1e5)
                else:
                    C_rot_zxy2zrt_vec(-azrad, dblsyn)

                for i in range(3):
                    XX[0, i, ix, iy] = dblsyn[i]
                    if calc_upar:
                        for k in range(3):
                            XX[i+1, k, ix, iy] = dblupar[k + i*3]


    # 将XX数组分到字典中
    if ZNE:
        chs = ['Z', 'N', 'E']

    for ityp in range(calcUTypes):
        c1 = '' if ityp==0 else chs[ityp-1].lower()
        for c in range(3):
            resDct[f"{c1}{s_compute_type}{chs[c]}"] = XX[ityp, c]
                
    return resDct


def _data_zrt2zne(stall:Stream):
    r"""
        将位移分量和位移空间导数分量转为ZNE坐标系

        :param     stall:     柱坐标系(zrt)下合成地震图

        :return:
            - **stream** - :class:`obspy.Stream` 类型
    """

    chs = ['Z', 'R', 'T']

    midname = stall[0].stats.channel[-3:-1]

    synLst:List[Trace] = []  # 顺序要求Z, R, T
    uparLst:List[Trace] = [] # 顺序要求zXXZ, zXXR, zXXT, rXXZ, rXXR, rXXT, tXXZ, tXXR, tXXT
    stsyn_upar = Stream()
    for ch in chs:
        st = stall.select(channel=f"{midname}{ch}")
        if len(st) == 1:
            synLst.append(st[0])
        
        for ch2 in chs:
            st = stall.select(channel=f"{ch.lower()}{midname}{ch2}")
            if len(st) == 1:
                uparLst.append(st[0])

    if len(synLst) != 3:
        raise ValueError(f"WRONG! synLst should have 3 components.")
    if len(stsyn_upar) != 0 and len(stsyn_upar) != 9:
        raise ValueError(f"WRONG! stsyn_upar should have 0 or 9 components.")

    
    # 是否有空间导数
    doupar = (len(uparLst) == 9)
    
    nt = stall[0].stats.npts
    azrad = np.deg2rad(stall[0].stats.sac['az'])
    dist = stall[0].stats.sac['dist']

    dblsyn = (c_double * 3)()
    dbleupar = (c_double * 9)()

    # 对每一个时间点
    for n in range(nt):
        # 复制数据
        for i1 in range(3):
            dblsyn[i1] = synLst[i1].data[n]
            if doupar:
                for i2 in range(3):
                    dbleupar[i2 + i1*3] = uparLst[i2 + i1*3].data[n]
        
        if doupar:
            C_rot_zrt2zxy_upar(azrad, dblsyn, dbleupar, dist*1e5)
        else:
            C_rot_zxy2zrt_vec(-azrad, dblsyn)

        # 将结果写入原数组
        for i1 in range(3):
            synLst[i1].data[n] = dblsyn[i1]
            if doupar:
                for i2 in range(3):
                    uparLst[i2 + i1*3].data[n] = dbleupar[i2 + i1*3]

    znechs = ['Z', 'N', 'E']
    # 修改通道名
    for i1 in range(3):
        ch1 = znechs[i1]
        tr = synLst[i1]
        tr.stats.channel = tr.stats.sac['kcmpnm'] = f'{midname}{ch1}'
        if doupar:
            for i2 in range(3):
                ch2 = znechs[i2]
                tr = uparLst[i2 + i1*3]
                tr.stats.channel = tr.stats.sac['kcmpnm'] = f'{ch1.lower()}{midname}{ch2}'

    stres = Stream()
    stres.extend(synLst)
    if doupar:
        stres.extend(uparLst)

    return stres


def _set_source_coef(
    par_theta:bool, coef:float, compute_type:str, M0:float, azrad:float,
    fZ=None, fN=None, fE=None, 
    strike=None, dip=None, rake=None, 
    MT=None, **kwargs):
    r"""
        设置不同震源的方向因子矩阵

        :param    par_theta:       是否求对theta的偏导
        :param    coef:            比例系数
        :param    compute_type:    计算类型，应为以下之一: 
            'COMPUTE_EXP'(爆炸源), 'COMPUTE_SF'(单力源),
            'COMPUTE_DC'(双力偶源), 'COMPUTE_MT'(矩张量源)
        :param    M0:              地震矩
        :param    azrad:           方位角(弧度)

        - 其他参数根据计算类型可选:
            - 单力源需要: fZ, fN, fE, 
            - 双力偶源需要: strike, dip, rake
            - 矩张量源需要: MT=(Mxx, Mxy, Mxz, Myy, Myz, Mzz)
    """
    
    caz = np.cos(azrad)
    saz = np.sin(azrad)

    src_coef = np.zeros((3, 6), dtype='f8')
    
    # 计算乘法因子
    if compute_type == 'COMPUTE_SF':
        mult = 1e-15 * M0 * coef 
    else:
        mult = 1e-20 * M0 * coef 

    # 根据不同计算类型处理
    if compute_type == 'COMPUTE_EXP':
        # 爆炸源情况
        src_coef[0, 0] = src_coef[1, 0] = 0.0 if par_theta else mult  # Z/R分量
        src_coef[2, 0] = 0.0  # T分量
    
    elif compute_type == 'COMPUTE_SF':
        # 单力源情况
        # 计算各向异性系数
        A0 = fZ * mult
        A1 = (fN * caz + fE * saz) * mult
        A4 = (-fN * saz + fE * caz) * mult

        # 设置震源系数矩阵 (公式4.6.20)
        src_coef[0, 1] = src_coef[1, 1] = 0.0 if par_theta else A0  # VF, Z/R
        src_coef[0, 2] = src_coef[1, 2] = A4 if par_theta else A1  # HF, Z/R
        src_coef[2, 1] = 0.0  # VF, T
        src_coef[2, 2] = -A1 if par_theta else A4  # HF, T
    
    elif compute_type == 'COMPUTE_DC':
        # 双力偶源情况 (公式4.8.35)
        # 计算各种角度值(转为弧度)
        stkrad = np.deg2rad(strike)  # 走向角
        diprad = np.deg2rad(dip)    # 倾角
        rakrad = np.deg2rad(rake)   # 滑动角
        therad = azrad - stkrad  # 方位角与走向角差
        
        # 计算各种三角函数值
        srak = np.sin(rakrad);      crak = np.cos(rakrad)
        sdip = np.sin(diprad);      cdip = np.cos(diprad)
        sdip2 = 2.0 * sdip * cdip;  cdip2 = 2.0 * cdip**2 - 1.0
        sthe = np.sin(therad);      cthe = np.cos(therad)
        sthe2 = 2.0 * sthe * cthe;  cthe2 = 2.0 * cthe**2 - 1.0

        # 计算各向异性系数
        A0 = mult * (0.5 * sdip2 * srak)
        A1 = mult * (cdip * crak * cthe - cdip2 * srak * sthe)
        A2 = mult * (0.5 * sdip2 * srak * cthe2 + sdip * crak * sthe2)
        A4 = mult * (-cdip2 * srak * cthe - cdip * crak * sthe)
        A5 = mult * (sdip * crak * cthe2 - 0.5 * sdip2 * srak * sthe2)

        # 设置震源系数矩阵
        src_coef[0, 3] = src_coef[1, 3] = 0.0 if par_theta else A0  # DD, Z/R
        src_coef[0, 4] = src_coef[1, 4] = A4 if par_theta else A1  # DS, Z/R
        src_coef[0, 5] = src_coef[1, 5] = 2.0 * A5 if par_theta else A2  # SS, Z/R
        src_coef[2, 3] = 0.0  # DD, T
        src_coef[2, 4] = -A1 if par_theta else A4  # DS, T
        src_coef[2, 5] = -2.0 * A2 if par_theta else A5  # DS, T
    
    elif compute_type == 'COMPUTE_MT':
        # 矩张量源情况 (公式4.9.7，修改了各向同性项)
        # 初始化矩张量分量
        M11, M12, M13, M22, M23, M33 = MT
        
        # 计算各向同性部分并减去
        Mexp = (M11 + M22 + M33) / 3.0
        M11 -= Mexp
        M22 -= Mexp
        M33 -= Mexp
        
        # 计算方位角的2倍角三角函数
        caz2 = np.cos(2 * azrad)
        saz2 = np.sin(2 * azrad)
        
        # 计算各向异性系数
        A0 = mult * ((2.0 * M33 - M11 - M22) / 6.0)
        A1 = mult * (- (M13 * caz + M23 * saz))
        A2 = mult * (0.5 * (M11 - M22) * caz2 + M12 * saz2)
        A4 = mult * (M13 * saz - M23 * caz)
        A5 = mult * (-0.5 * (M11 - M22) * saz2 + M12 * caz2)

        # 设置震源系数矩阵
        src_coef[0, 0] = src_coef[1, 0] = 0.0 if par_theta else mult * Mexp  # EX, Z/R
        src_coef[0, 3] = src_coef[1, 3] = 0.0 if par_theta else A0  # DD, Z/R
        src_coef[0, 4] = src_coef[1, 4] = A4 if par_theta else A1  # DS, Z/R
        src_coef[0, 5] = src_coef[1, 5] = 2.0 * A5 if par_theta else A2  # SS, Z/R
        src_coef[2, 0] = 0.0  # EX, T
        src_coef[2, 3] = 0.0  # DD, T
        src_coef[2, 4] = -A1 if par_theta else A4  # DS, T
        src_coef[2, 5] = -2.0 * A2 if par_theta else A5  # DS, T


    return src_coef


def gen_syn_from_gf_DC(st:Stream|dict, M0:float, strike:float, dip:float, rake:float, az:float=-999, ZNE=False, calc_upar:bool=False):
    '''
        双力偶源，角度单位均为度   

        :param    st:       计算好的时域格林函数, :class:`obspy.Stream` 类型，或者静态格林函数（字典类型）
        :param    M0:       标量地震矩, 单位dyne*cm
        :param    strike:   走向，以北顺时针为正，0<=strike<=360
        :param    dip:      倾角，以水平面往下为正，0<=dip<=90
        :param    rake:     滑动角，在断层面相对于走向方向逆时针为正，-180<=rake<=180
        :param    az:       台站方位角，以北顺时针为正，0<=az<=360（静态情况不需要）
        :param    ZNE:             是否以ZNE分量输出?
        :param    calc_upar:     是否计算位移u的空间导数

        :return:
            - **stream** -  三分量地震图, :class:`obspy.Stream` 类型
    '''
    if isinstance(st, Stream):
        if az > 360 or az < -360:
            raise ValueError(f"WRONG azimuth ({az})")
        return _gen_syn_from_gf(st, calc_upar, "COMPUTE_DC", M0, az, ZNE, strike=strike, dip=dip, rake=rake)
    elif isinstance(st, dict):
        return _gen_syn_from_static_gf(st, calc_upar, "COMPUTE_DC", M0, ZNE, strike=strike, dip=dip, rake=rake)
    else:
        raise NotImplementedError


def gen_syn_from_gf_SF(st:Stream|dict, S:float, fN:float, fE:float, fZ:float, az:float=-999, ZNE=False, calc_upar:bool=False):
    '''
        单力源，力分量单位均为dyne   

        :param    st:    计算好的时域格林函数, :class:`obspy.Stream` 类型，或者静态格林函数（字典类型）
        :param     S:    力的放大系数
        :param    fN:    北向力，向北为正  
        :param    fE:    东向力，向东为正  
        :param    fZ:    垂向力，向下为正  
        :param    az:    台站方位角，以北顺时针为正，0<=az<=360 （静态情况不需要）  
        :param    ZNE:             是否以ZNE分量输出?
        :param    calc_upar:     是否计算位移u的空间导数

        :return:
            - **stream** - 三分量地震图, :class:`obspy.Stream` 类型
    '''
    if isinstance(st, Stream):
        if az > 360 or az < -360:
            raise ValueError(f"WRONG azimuth ({az})")
        return _gen_syn_from_gf(st, calc_upar, "COMPUTE_SF", S, az, ZNE, fN=fN, fE=fE, fZ=fZ)
    elif isinstance(st, dict):
        return _gen_syn_from_static_gf(st, calc_upar, "COMPUTE_SF", S, ZNE, fN=fN, fE=fE, fZ=fZ)
    else:
        raise NotImplementedError


def gen_syn_from_gf_EXP(st:Stream|dict, M0:float, az:float=-999, ZNE=False, calc_upar:bool=False):
    '''
        爆炸源

        :param    st:          计算好的时域格林函数, :class:`obspy.Stream` 类型，或者静态格林函数（字典类型）
        :param    M0:          标量地震矩, 单位dyne*cm  
        :param    az:          台站方位角，以北顺时针为正，0<=az<=360 [不用于计算] （静态情况不需要） 
        :param    ZNE:             是否以ZNE分量输出?
        :param    calc_upar:     是否计算位移u的空间导数

        :return:
            - **stream** -       三分量地震图, :class:`obspy.Stream` 类型
    '''
    if isinstance(st, Stream):
        if az > 360 or az < -360:
            raise ValueError(f"WRONG azimuth ({az})")
        return _gen_syn_from_gf(st, calc_upar, "COMPUTE_EXP", M0, az, ZNE)
    elif isinstance(st, dict):
        return _gen_syn_from_static_gf(st, calc_upar, "COMPUTE_EXP", M0, ZNE)
    else:
        raise NotImplementedError
    

def gen_syn_from_gf_MT(st:Stream|dict, M0:float, MT:ArrayLike, az:float=-999, ZNE=False, calc_upar:bool=False):
    ''' 
        矩张量源，单位为dyne*cm

        :param    st:          计算好的时域格林函数, :class:`obspy.Stream` 类型，或者静态格林函数（字典类型）
        :param    M0:          标量地震矩
        :param    MT:          矩张量 (M11, M12, M13, M22, M23, M33),下标1,2,3分别代表北向，东向，垂直向上  
        :param    az:          台站方位角，以北顺时针为正，0<=az<=360 （静态情况不需要）  
        :param    ZNE:             是否以ZNE分量输出?
        :param    calc_upar:     是否计算位移u的空间导数

        :return:
            - **stream** -       三分量地震图, :class:`obspy.Stream` 类型
    '''
    if isinstance(st, Stream):
        if az > 360 or az < -360:
            raise ValueError(f"WRONG azimuth ({az})")
        return _gen_syn_from_gf(st, calc_upar, "COMPUTE_MT", M0, az, ZNE, MT=MT)
    elif isinstance(st, dict):
        return _gen_syn_from_static_gf(st, calc_upar, "COMPUTE_MT", M0, ZNE, MT=MT)
    else:
        raise NotImplementedError


#=================================================================================================================
#
#                                           根据几何方程和本构方程合成应力和应变
#
#=================================================================================================================


def _compute_strain(st_syn:Stream):
    r"""
        Compute dynamic strain from synthetic spatial derivatives.

        :param     st_syn:      synthetic spatial derivatives.

        :return:
            - **stream** -  dynamic strain (unit: dyne/cm^2 = 0.1 Pa), in :class:`obspy.Stream` class.
    """

    midname = st_syn[0].stats.channel[-3:-1]

    # 检查是否每个trace是否具有相同midname
    for tr in st_syn:
        if tr.stats.channel[-3:-1] != midname:
            raise ValueError("WRONG INPUT! inconsistent component names.")
        
    zrtchs = ['Z', 'R', 'T']
    znechs = ['Z', 'N', 'E']
    chs = zrtchs

    # 判断是否有标志性的trace
    if len(st_syn.select(channel=f"n{midname}N")) > 0:
        chs = znechs

    dist = st_syn[0].stats.sac['dist']

    # ----------------------------------------------------------------------------------
    # 循环6个分量
    stres = Stream()
    for i1 in range(3):
        c1 = chs[i1]
        for i2 in range(i1, 3):
            c2 = chs[i2]

            channel = f"{c1.lower()}{midname}{c2}"
            st = st_syn.select(channel=channel)
            if len(st) == 0:
                raise NameError(f"{channel} not exists.")
            tr = st[0].copy()

            channel = f"{c2.lower()}{midname}{c1}"
            st = st_syn.select(channel=channel)
            if len(st) == 0:
                raise NameError(f"{channel} not exists.")
            tr.data = (tr.data + st[0].data) * 0.5

            # 特殊情况加上协变导数
            if c1=='R' and c2=='T':
                channel = f"{midname}T"
                st = st_syn.select(channel=channel)
                if len(st) == 0:
                    raise NameError(f"{channel} not exists.")
                tr.data -= 0.5*st[0].data / dist * 1e-5
            
            elif c1=='T' and c2=='T':
                channel = f"{midname}R"
                st = st_syn.select(channel=channel)
                if len(st) == 0:
                    raise NameError(f"{channel} not exists.")
                tr.data += st[0].data / dist * 1e-5

            # 修改通道名
            tr.stats.channel = tr.stats.sac['kcmpnm'] = f"{c1}{c2}"

            stres.append(tr)

    return stres


def _compute_static_strain(syn:dict):
    r"""
        Compute static strain from synthetic spatial derivatives.

        :param     syn:      synthetic spatial derivatives.

        :return:
            - **res** -  static strain, in dict class.
    """

    midname = ""
    # 检查是否每个分量是否具有相同midname
    for k in syn.keys():
        if k[0] == '_':
            continue
        if len(midname)==0:
            midname = k[-3:-1]

        if k[-3:-1] != midname:
            raise ValueError("WRONG INPUT! inconsistent component names.")
        
    zrtchs = ['Z', 'R', 'T']
    znechs = ['Z', 'N', 'E']
    chs = zrtchs

    # 判断是否有标志性的分量名
    if f"n{midname}N" in syn.keys():
        chs = znechs

    xarr:np.ndarray = syn['_xarr']
    yarr:np.ndarray = syn['_yarr']

    # 结果字典
    resDct = {}

    # 基本数据拷贝
    for k in syn.keys():
        if k[0] != '_':
            continue 
        resDct[k] = deepcopy(syn[k])

    # 6个分量建立数组
    for i1 in range(3):
        c1 = chs[i1]
        for i2 in range(i1, 3):
            c2 = chs[i2]
            channel = f"{c1}{c2}"
            resDct[channel] = np.zeros((len(xarr), len(yarr)), dtype='f8')


    for iy in range(len(yarr)):
        for ix in range(len(xarr)):
            # 震中距
            dist = max(np.sqrt(xarr[ix]**2 + yarr[iy]**2), 1e-5)

            # ----------------------------------------------------------------------------------
            # 循环6个分量
            for i1 in range(3):
                c1 = chs[i1]
                for i2 in range(i1, 3):
                    c2 = chs[i2]

                    channel = f"{c1.lower()}{midname}{c2}"
                    v12 = syn[channel][ix, iy]

                    channel = f"{c2.lower()}{midname}{c1}"
                    v21 = syn[channel][ix, iy]

                    val = 0.5*(v12 + v21)

                    # 特殊情况加上协变导数
                    if c1=='R' and c2=='T':
                        channel = f"{midname}T"
                        v0 = syn[channel][ix, iy]
                        val -= 0.5*v0 / dist * 1e-5

                    elif c1=='T' and c2=='T':
                        channel = f"{midname}R"
                        v0 = syn[channel][ix, iy]
                        val += v0 / dist * 1e-5

                    channel = f"{c1}{c2}"
                    resDct[channel][ix, iy] = val

    return resDct


def compute_strain(st:Stream|dict):
    if isinstance(st, Stream):
        return _compute_strain(st)
    elif isinstance(st, dict):
        return _compute_static_strain(st)
    else:
        raise NotImplementedError


def _compute_stress(st_syn:Stream):
    r"""
        Compute dynamic stress from synthetic spatial derivatives.

        :param     st_syn:      synthetic spatial derivatives.

        :return:
            - **stream** -  dynamic stress, in :class:`obspy.Stream` class.
    """

    # 由于有Q值的存在，lambda和mu变成了复数，需在频域进行

    midname = st_syn[0].stats.channel[-3:-1]

    # 检查是否每个trace是否具有相同midname
    for tr in st_syn:
        if tr.stats.channel[-3:-1] != midname:
            raise ValueError("WRONG INPUT! inconsistent component names.")
        
    zrtchs = ['Z', 'R', 'T']
    znechs = ['Z', 'N', 'E']
    chs = zrtchs
    rot2ZNE:bool = False

    # 判断是否有标志性的trace
    if len(st_syn.select(channel=f"n{midname}N")) > 0:
        chs = znechs
        rot2ZNE = True

    nt = st_syn[0].stats.npts
    dt = st_syn[0].stats.delta
    dist = st_syn[0].stats.sac['dist']
    df = 1.0/(nt*dt)
    nf = nt//2 + 1
    va = st_syn[0].stats.sac['user1']
    vb = st_syn[0].stats.sac['user2']
    rho = st_syn[0].stats.sac['user3']
    Qainv = st_syn[0].stats.sac['user4']
    Qbinv = st_syn[0].stats.sac['user5']

    # 计算不同频率下的拉梅系数
    mus = np.zeros((nf,), dtype='c16')
    lams = np.zeros((nf,), dtype='c16')
    omega = (REAL*2)(0.0, 0.0)
    atte = (REAL*2)(0.0, 0.0)
    for i in range(nf):
        freq = 0.01 if i==0 else df*i 
        w = 2.0*np.pi*freq 
        omega[0] = w
        C_py_attenuation_law(Qbinv, omega, atte)
        attb = atte[0] + atte[1]*1j
        mus[i] = vb*vb*attb*attb*rho*1e10
        C_py_attenuation_law(Qainv, omega, atte)
        atta = atte[0] + atte[1]*1j
        lams[i] = va*va*atta*atta*rho*1e10 - 2.0*mus[i]
    
    del omega, atte

    # ----------------------------------------------------------------------------------
    # 先计算体积应变u_kk = u_11 + u22 + u33 和 lamda的乘积
    lam_ukk = np.zeros((nf,), dtype='c16')
    for i in range(3):
        c = chs[i]
        channel = f"{c.lower()}{midname}{c}"
        st = st_syn.select(channel=channel)
        if len(st) == 0:
            raise NameError(f"{channel} not exists.")
        lam_ukk[:] += rfft(st[0].data, nt)

    # 加上协变导数
    if not rot2ZNE:
        channel = f"{midname}R"
        st = st_syn.select(channel=channel)
        if len(st) == 0:
            raise NameError(f"{channel} not exists.")
        lam_ukk[:] += rfft(st[0].data, nt) / dist * 1e-5

    lam_ukk[:] *= lams 

    # ----------------------------------------------------------------------------------
    # 循环6个分量
    stres = Stream()
    for i1 in range(3):
        c1 = chs[i1]
        for i2 in range(i1, 3):
            c2 = chs[i2]

            channel = f"{c1.lower()}{midname}{c2}"
            st = st_syn.select(channel=channel)
            if len(st) == 0:
                raise NameError(f"{channel} not exists.")
            tr = st[0].copy()
            fftarr = np.zeros((nf,), dtype='c16')

            channel = f"{c2.lower()}{midname}{c1}"
            st = st_syn.select(channel=channel)
            if len(st) == 0:
                raise NameError(f"{channel} not exists.")
            fftarr[:] = rfft(tr.data + st[0].data, nt) * mus

            # 对于对角线分量，需加上lambda * u_kk
            if c1==c2:
                fftarr[:] += lam_ukk

            # 特殊情况加上协变导数
            if c1=='R' and c2=='T':
                channel = f"{midname}T"
                st = st_syn.select(channel=channel)
                if len(st) == 0:
                    raise NameError(f"{channel} not exists.")
                fftarr[:] -= mus*rfft(st[0].data, nt) / dist * 1e-5
            
            elif c1=='T' and c2=='T':
                channel = f"{midname}R"
                st = st_syn.select(channel=channel)
                if len(st) == 0:
                    raise NameError(f"{channel} not exists.")
                fftarr[:] += 2.0*mus*rfft(st[0].data, nt) / dist * 1e-5

            # 修改通道名
            tr.stats.channel = tr.stats.sac['kcmpnm'] = f"{c1}{c2}"

            tr.data = irfft(fftarr, nt)

            stres.append(tr)

    return stres


def _compute_static_stress(syn:dict):
    r"""
        Compute static stress from synthetic spatial derivatives.

        :param     syn:      synthetic spatial derivatives.

        :return:
            - **res** -  static stress (unit: dyne/cm^2 = 0.1 Pa), in dict class.
    """

    midname = ""
    # 检查是否每个分量是否具有相同midname
    for k in syn.keys():
        if k[0] == '_':
            continue
        if len(midname)==0:
            midname = k[-3:-1]

        if k[-3:-1] != midname:
            raise ValueError("WRONG INPUT! inconsistent component names.")
        
    zrtchs = ['Z', 'R', 'T']
    znechs = ['Z', 'N', 'E']
    chs = zrtchs
    rot2ZNE:bool = False

    # 判断是否有标志性的分量名
    if f"n{midname}N" in syn.keys():
        chs = znechs
        rot2ZNE = True

    xarr:np.ndarray = syn['_xarr']
    yarr:np.ndarray = syn['_yarr']
    va = syn['_rcv_va']
    vb = syn['_rcv_vb']
    rho = syn['_rcv_rho']
    mu = vb*vb*rho*1e10
    lam = va*va*rho*1e10 - 2.0*mu

    # 结果字典
    resDct = {}

    # 基本数据拷贝
    for k in syn.keys():
        if k[0] != '_':
            continue 
        resDct[k] = deepcopy(syn[k])

    # 6个分量建立数组
    for i1 in range(3):
        c1 = chs[i1]
        for i2 in range(i1, 3):
            c2 = chs[i2]
            channel = f"{c1}{c2}"
            resDct[channel] = np.zeros((len(xarr), len(yarr)), dtype='f8')


    for iy in range(len(yarr)):
        for ix in range(len(xarr)):
            # 震中距
            dist = max(np.sqrt(xarr[ix]**2 + yarr[iy]**2), 1e-5)

            # ----------------------------------------------------------------------------------
            # 先计算体积应变u_kk = u_11 + u22 + u33 和 lamda的乘积
            lam_ukk = 0.0
            for i in range(3):
                c = chs[i]
                channel = f"{c.lower()}{midname}{c}"
                lam_ukk += syn[channel][ix, iy]
            
            # 加上协变导数
            if not rot2ZNE:
                channel = f"{midname}R"
                lam_ukk += syn[channel][ix, iy] / dist * 1e-5
            
            lam_ukk *= lam

            # ----------------------------------------------------------------------------------
            # 循环6个分量
            for i1 in range(3):
                c1 = chs[i1]
                for i2 in range(i1, 3):
                    c2 = chs[i2]

                    channel = f"{c1.lower()}{midname}{c2}"
                    v12 = syn[channel][ix, iy]

                    channel = f"{c2.lower()}{midname}{c1}"
                    v21 = syn[channel][ix, iy]

                    val = mu*(v12 + v21)

                    # 对于对角线分量，需加上lambda * u_kk
                    if c1==c2:
                        val += lam_ukk

                    # 特殊情况加上协变导数
                    if c1=='R' and c2=='T':
                        channel = f"{midname}T"
                        v0 = syn[channel][ix, iy]
                        val -= mu*v0 / dist * 1e-5

                    elif c1=='T' and c2=='T':
                        channel = f"{midname}R"
                        v0 = syn[channel][ix, iy]
                        val += 2.0*mu*v0 / dist * 1e-5

                    channel = f"{c1}{c2}"
                    resDct[channel][ix, iy] = val

    return resDct


def compute_stress(st:Stream|dict):
    if isinstance(st, Stream):
        return _compute_stress(st)
    elif isinstance(st, dict):
        return _compute_static_stress(st)
    else:
        raise NotImplementedError


def __check_trace_attr_sac(tr:Trace, **kwargs):
    '''
        临时函数，检查trace中是否有sac字典，并将kwargs内容填入  
    '''
    if hasattr(tr.stats, 'sac'):
        for k, v in kwargs.items():
            tr.stats.sac[k] = v 
    else: 
        tr.stats.sac = AttribDict(**kwargs)


#=================================================================================================================
#
#                                           卷积、微分、积分、保存SAC
#
#=================================================================================================================



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


def stream_convolve(st0:Stream, signal0:np.ndarray, inplace=True):
    '''
        对stream中每一道信号做线性卷积  

        :param    st0:        记录多个Trace的 :class:`obspy.Stream` 类型
        :param    signal0:    卷积信号
        :param    inplace:    是否做原地更改  

        :return:
            - **stream** -    处理后的结果, :class:`obspy.Stream` 类型
    '''
    st = st0 if inplace else deepcopy(st0)
    signal = deepcopy(signal0)
    
    for tr in st:
        data = tr.data 
        dt = tr.stats.delta
        
        fac = None
        user_wI = hasattr(tr.stats, "sac") and "user0" in tr.stats.sac
        # 使用虚频率先压制
        if user_wI:
            npts = tr.stats.npts
            wI = tr.stats.sac['user0']
            fac = np.exp(np.arange(0, npts)*dt*wI)
            signal = deepcopy(signal0)

            signal[:] /= fac[:len(signal)]
            data[:] /= fac

        data1 = np.pad(data, (len(signal)-1, 0), mode='wrap') # 强制循环卷
        data[:] = oaconvolve(data1, signal, mode='valid')[:data.shape[0]] * dt  # dt是连续卷积的系数

        if user_wI:
            data[:] *= fac

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
        :param    path:       保存文件名，不要加后缀，
                              例如path="GRN/out"表明sac文件将保存在GRN文件中，
                              文件名分别为outR.sac, outT.sac, outZ.sac，
                              若通道名前缀有"r","z"或"t"（表示位移空间导数），则该前缀
                              也会加到SAC文件名中，例如GRN/routR.sac

    '''
    # 新建对应文件夹
    parent_dir = os.path.dirname(path)
    if(len(parent_dir)>0):
        os.makedirs(parent_dir, exist_ok=True)

    filesub = os.path.basename(path)
    # 每一道的保存路径为path+{channel}
    for tr in st:
        prefix = ""
        if tr.stats.channel[0] in ['r', 't', 'z']:
            prefix = tr.stats.channel[0]
        filepath = os.path.join(parent_dir, f"{prefix}{filesub}{tr.stats.channel[-1]}.sac")
        tr.write(filepath, format='SAC')





#=================================================================================================================
#
#                                           积分过程文件读取及绘制
#
#=================================================================================================================


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