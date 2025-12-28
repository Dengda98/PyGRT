"""
    :file:     utils.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包含一些数据处理操作上的补充:   
        1、剪切源、单力源、爆炸源、矩张量源 通过格林函数合成理论地震图的函数\n
        2、Stream类型的时域卷积、微分、积分 (基于numpy和scipy)    \n
        3、Stream类型写到本地sac文件，自定义名称    \n
        4、读取波数积分和峰谷平均法过程文件  \n
        5、其它辅助函数  \n

"""


import numpy as np
import numpy.ctypeslib as npct
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from obspy import Stream, Trace 
from obspy.core import AttribDict
from copy import deepcopy
from scipy.signal import oaconvolve
from scipy.fft import rfft, irfft
from scipy.special import jv
from scipy.interpolate import interpn
import math 
import os
import glob
from typing import List, Union
from copy import deepcopy

from numpy.typing import ArrayLike

from .c_interfaces import *


__all__ = [
    "gen_syn_from_gf_DC",
    "gen_syn_from_gf_SF",
    "gen_syn_from_gf_EX",
    "gen_syn_from_gf_MT",

    "compute_strain",
    "compute_rotation",
    "compute_stress",

    "stream_convolve",
    "stream_integral",
    "stream_diff",
    "stream_write_sac",

    "read_kernels_freqs",
    "read_statsfile",
    "read_statsfile_ptam",
    "plot_statsdata",
    "plot_statsdata_ptam",

    "solve_lamb1"
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
                                    'COMPUTE_EX'(爆炸源), 'COMPUTE_SF'(单力源),
                                    'COMPUTE_DC'(剪切源), 'COMPUTE_MT'(矩张量源)
        :param    M0:              标量地震矩, 单位dyne*cm
        :param    az:              方位角(度)
        :param    ZNE:             是否以ZNE分量输出?
            
    """
    chs = ZRTchs
    sacin_prefixes = ["", "z", "r", ""]   # 输入通道名
    sacout_prefixes = ["", "z", "r", "t"]   # 输出通道名
    srcName = ["EX", "VF", "HF", "DD", "DS", "SS"]
    allchs = [tr.stats.channel for tr in st]

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

        srcRadi = _set_source_radi(ityp==3, upar_scale, compute_type, M0, azrad, **kwargs)

        inpref = sacin_prefixes[ityp]
        outpref = sacout_prefixes[ityp]

        for c in range(CHANNEL_NUM):
            ch = chs[c]
            tr:Trace = st[0].copy()
            tr.data[:] = 0.0
            tr.stats.channel = kcmpnm = f'{outpref}{ch}'
            __check_trace_attr_sac(tr, az=az, baz=baz, kcmpnm=kcmpnm)
            for k in range(SRC_M_NUM):
                coef = srcRadi[k, c]
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
                                    'COMPUTE_EX'(爆炸源), 'COMPUTE_SF'(单力源),
                                    'COMPUTE_DC'(剪切源), 'COMPUTE_MT'(矩张量源)
        :param    M0:              标量地震矩, 单位dyne*cm
        :param    ZNE:             是否以ZNE分量输出?
            
    """
    chs = ZRTchs
    sacin_prefixes = ["", "z", "r", ""]   # 输入通道名
    srcName = ["EX", "VF", "HF", "DD", "DS", "SS"]
    allchs = list(grn.keys())

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

                srcRadi = _set_source_radi(ityp==3, upar_scale, compute_type, M0, azrad, **kwargs)

                inpref = sacin_prefixes[ityp]

                for c in range(CHANNEL_NUM):
                    ch = chs[c]
                    for k in range(SRC_M_NUM):
                        coef = srcRadi[k, c]
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
                    C_grt_rot_zrt2zxy_upar(azrad, dblsyn, dblupar, dist*1e5)
                else:
                    C_grt_rot_zxy2zrt_vec(-azrad, dblsyn)

                for i in range(3):
                    XX[0, i, ix, iy] = dblsyn[i]
                    if calc_upar:
                        for k in range(3):
                            XX[i+1, k, ix, iy] = dblupar[k + i*3]


    # 将XX数组分到字典中
    if ZNE:
        chs = ZNEchs

    for ityp in range(calcUTypes):
        c1 = '' if ityp==0 else chs[ityp-1].lower()
        for c in range(3):
            resDct[f"{c1}{chs[c]}"] = XX[ityp, c].copy()
                
    return resDct


def _data_zrt2zne(stall:Stream):
    r"""
        将位移分量和位移空间导数分量转为ZNE坐标系

        :param     stall:     柱坐标系(zrt)下合成地震图

        :return:
            - **stream** - :class:`obspy.Stream` 类型
    """

    chs = ZRTchs

    synLst:List[Trace] = []  # 顺序要求Z, R, T
    uparLst:List[Trace] = [] # 顺序要求zXXZ, zXXR, zXXT, rXXZ, rXXR, rXXT, tXXZ, tXXR, tXXT
    stsyn_upar = Stream()
    for ch in chs:
        st = stall.select(channel=f"{ch}")
        if len(st) == 1:
            synLst.append(st[0])
        
        for ch2 in chs:
            st = stall.select(channel=f"{ch.lower()}{ch2}")
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
            C_grt_rot_zrt2zxy_upar(azrad, dblsyn, dbleupar, dist*1e5)
        else:
            C_grt_rot_zxy2zrt_vec(-azrad, dblsyn)

        # 将结果写入原数组
        for i1 in range(3):
            synLst[i1].data[n] = dblsyn[i1]
            if doupar:
                for i2 in range(3):
                    uparLst[i2 + i1*3].data[n] = dbleupar[i2 + i1*3]

    # 修改通道名
    for i1 in range(3):
        ch1 = ZNEchs[i1]
        tr = synLst[i1]
        tr.stats.channel = tr.stats.sac['kcmpnm'] = f'{ch1}'
        if doupar:
            for i2 in range(3):
                ch2 = ZNEchs[i2]
                tr = uparLst[i2 + i1*3]
                tr.stats.channel = tr.stats.sac['kcmpnm'] = f'{ch1.lower()}{ch2}'

    stres = Stream()
    stres.extend(synLst)
    if doupar:
        stres.extend(uparLst)

    return stres


def _set_source_radi(
    par_theta:bool, coef:float, compute_type:str, M0:float, azrad:float,
    fZ=None, fN=None, fE=None, 
    strike=None, dip=None, rake=None, 
    MT=None, **kwargs):
    r"""
        设置不同震源的方向因子矩阵

        :param    par_theta:       是否求对theta的偏导
        :param    coef:            比例系数
        :param    compute_type:    计算类型，应为以下之一: 
                                    'COMPUTE_EX'(爆炸源), 'COMPUTE_SF'(单力源),
                                    'COMPUTE_DC'(剪切源), 'COMPUTE_MT'(矩张量源)
        :param    M0:              地震矩
        :param    azrad:           方位角(弧度)

        - 其他参数根据计算类型可选:
            - 单力源需要: fZ, fN, fE, 
            - 剪切源需要: strike, dip, rake
            - 矩张量源需要: MT=(Mxx, Mxy, Mxz, Myy, Myz, Mzz)
    """
    
    caz = np.cos(azrad)
    saz = np.sin(azrad)

    src_coef = np.zeros((SRC_M_NUM, CHANNEL_NUM), dtype='f8')
    
    # 计算乘法因子
    if compute_type == 'COMPUTE_SF':
        mult = 1e-15 * M0 * coef 
    else:
        mult = 1e-20 * M0 * coef 

    # 根据不同计算类型处理
    if compute_type == 'COMPUTE_EX':
        # 爆炸源情况
        src_coef[0, 0] = src_coef[0, 1] = 0.0 if par_theta else mult  # Z/R分量
        src_coef[0, 2] = 0.0  # T分量
    
    elif compute_type == 'COMPUTE_SF':
        # 单力源情况
        # 计算各向异性系数
        A0 = fZ * mult
        A1 = (fN * caz + fE * saz) * mult
        A4 = (-fN * saz + fE * caz) * mult

        # 设置震源系数矩阵 (公式4.6.20)
        src_coef[1, 0] = src_coef[1, 1] = 0.0 if par_theta else A0  # VF, Z/R
        src_coef[2, 0] = src_coef[2, 1] = A4 if par_theta else A1  # HF, Z/R
        src_coef[1, 2] = 0.0  # VF, T
        src_coef[2, 2] = -A1 if par_theta else A4  # HF, T
    
    elif compute_type == 'COMPUTE_DC':
        # 剪切源情况 (公式4.8.35)
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
        src_coef[3, 0] = src_coef[3, 1] = 0.0 if par_theta else A0  # DD, Z/R
        src_coef[4, 0] = src_coef[4, 1] = A4 if par_theta else A1  # DS, Z/R
        src_coef[5, 0] = src_coef[5, 1] = 2.0 * A5 if par_theta else A2  # SS, Z/R
        src_coef[3, 2] = 0.0  # DD, T
        src_coef[4, 2] = -A1 if par_theta else A4  # DS, T
        src_coef[5, 2] = -2.0 * A2 if par_theta else A5  # DS, T
    
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
        src_coef[0, 0] = src_coef[0, 1] = 0.0 if par_theta else mult * Mexp  # EX, Z/R
        src_coef[3, 0] = src_coef[3, 1] = 0.0 if par_theta else A0  # DD, Z/R
        src_coef[4, 0] = src_coef[4, 1] = A4 if par_theta else A1  # DS, Z/R
        src_coef[5, 0] = src_coef[5, 1] = 2.0 * A5 if par_theta else A2  # SS, Z/R
        src_coef[0, 2] = 0.0  # EX, T
        src_coef[3, 2] = 0.0  # DD, T
        src_coef[4, 2] = -A1 if par_theta else A4  # DS, T
        src_coef[5, 2] = -2.0 * A2 if par_theta else A5  # DS, T


    return src_coef


def gen_syn_from_gf_DC(st:Union[Stream,dict], M0:float, strike:float, dip:float, rake:float, az:float=-999, ZNE=False, calc_upar:bool=False):
    '''
        Shear source, the unit of angles is all degrees(°)

        :param    st:       Green's functions in a :class:`obspy.Stream` (dynamic-case) or a dict (static-case)
        :param    M0:       scalar seismic moment (dyne*cm)
        :param    strike:   0 <= strike <= 360 (north=0, clockwise as positive)
        :param    dip:      0 <= dip <= 90
        :param    rake:     -180 <= rake <= 180 (on the fault plane, counterclockwise as positive)
        :param    az:       azimuth, 0 <= az <= 360 (not used for static case)
        :param    ZNE:          whether output in 'ZNE'-coord, default is 'ZRT'
        :param    calc_upar:    whether calculate the spatial derivatives of displacements.

        :return:
            - **stream** -  :class:`obspy.Stream`
    '''
    if isinstance(st, Stream):
        if az > 360 or az < -360:
            raise ValueError(f"WRONG azimuth ({az})")
        return _gen_syn_from_gf(st, calc_upar, "COMPUTE_DC", M0, az, ZNE, strike=strike, dip=dip, rake=rake)
    elif isinstance(st, dict):
        return _gen_syn_from_static_gf(st, calc_upar, "COMPUTE_DC", M0, ZNE, strike=strike, dip=dip, rake=rake)
    else:
        raise NotImplementedError


def gen_syn_from_gf_SF(st:Union[Stream,dict], S:float, fN:float, fE:float, fZ:float, az:float=-999, ZNE=False, calc_upar:bool=False):
    '''
        Single-force source (dyne)  

        :param    st:    Green's functions in a :class:`obspy.Stream` (dynamic-case) or a dict (static-case)
        :param     S:    scaling factor (dyne)
        :param    fN:    coefficient of Northward force   
        :param    fE:    coefficient of Eastward force
        :param    fZ:    coefficient of Vertical(Downward) force 
        :param    az:    azimuth, 0 <= az <= 360 (not used for static case)
        :param    ZNE:          whether output in 'ZNE'-coord, default is 'ZRT'
        :param    calc_upar:    whether calculate the spatial derivatives of displacements.

        :return:
            - **stream** - :class:`obspy.Stream`
    '''
    if isinstance(st, Stream):
        if az > 360 or az < -360:
            raise ValueError(f"WRONG azimuth ({az})")
        return _gen_syn_from_gf(st, calc_upar, "COMPUTE_SF", S, az, ZNE, fN=fN, fE=fE, fZ=fZ)
    elif isinstance(st, dict):
        return _gen_syn_from_static_gf(st, calc_upar, "COMPUTE_SF", S, ZNE, fN=fN, fE=fE, fZ=fZ)
    else:
        raise NotImplementedError


def gen_syn_from_gf_EX(st:Union[Stream,dict], M0:float, az:float=-999, ZNE=False, calc_upar:bool=False):
    '''
        Explosion

        :param    st:          Green's functions in a :class:`obspy.Stream` (dynamic-case) or a dict (static-case)
        :param    M0:          scalar seismic moment (dyne*cm)
        :param    az:          azimuth, 0 <= az <= 360 (not used for static case)
        :param    ZNE:          whether output in 'ZNE'-coord, default is 'ZRT'
        :param    calc_upar:    whether calculate the spatial derivatives of displacements.

        :return:
            - **stream** -       :class:`obspy.Stream`
    '''
    if isinstance(st, Stream):
        if az > 360 or az < -360:
            raise ValueError(f"WRONG azimuth ({az})")
        return _gen_syn_from_gf(st, calc_upar, "COMPUTE_EX", M0, az, ZNE)
    elif isinstance(st, dict):
        return _gen_syn_from_static_gf(st, calc_upar, "COMPUTE_EX", M0, ZNE)
    else:
        raise NotImplementedError
    

def gen_syn_from_gf_MT(st:Union[Stream,dict], M0:float, MT:ArrayLike, az:float=-999, ZNE=False, calc_upar:bool=False):
    ''' 
        Moment tensor

        :param    st:          Green's functions in a :class:`obspy.Stream` (dynamic-case) or a dict (static-case)
        :param    M0:          scalar seismic moment (dyne*cm)
        :param    MT:          coefficient of Moment tensor (M11, M12, M13, M22, M23, M33), subscripts 1,2,3 denote Northward,Eastward,Downward
        :param    az:          azimuth, 0 <= az <= 360 (not used for static case)
        :param    ZNE:          whether output in 'ZNE'-coord, default is 'ZRT'
        :param    calc_upar:    whether calculate the spatial derivatives of displacements.

        :return:
            - **stream** -     :class:`obspy.Stream`
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
#                                           根据几何方程和本构方程合成应力、应变、旋转张量
#
#=================================================================================================================


def _compute_strain_rotation(st_syn:Stream, Type:str):
    r"""
        Compute dynamic strain/rotation tensor from synthetic spatial derivatives.

        :param     st_syn:      synthetic spatial derivatives.
        :param     Type:        "strain" or "rotation"

        :return:
            - **stream** -  dynamic strain/rotation tensor, in :class:`obspy.Stream` class.
    """

    if Type == 'strain':
        sgn = 1
        i1_end = 3
        i2_offset = 0
    elif Type == 'rotation':
        sgn = -1
        i1_end = 2
        i2_offset = 1
    else:
        raise ValueError(f"{Type} not supported.")
        
    chs = ZRTchs

    # 判断是否有标志性的trace
    if len(st_syn.select(channel=f"nN")) > 0:
        chs = ZNEchs

    dist = st_syn[0].stats.sac['dist']

    # ----------------------------------------------------------------------------------
    # 循环6/3个分量
    stres = Stream()
    for i1 in range(i1_end):
        c1 = chs[i1]
        for i2 in range(i1+i2_offset, 3):
            c2 = chs[i2]

            channel = f"{c2.lower()}{c1}"
            st = st_syn.select(channel=channel)
            if len(st) == 0:
                raise NameError(f"{channel} not exists.")
            tr = st[0].copy()

            channel = f"{c1.lower()}{c2}"
            st = st_syn.select(channel=channel)
            if len(st) == 0:
                raise NameError(f"{channel} not exists.")
            tr.data = (tr.data + sgn*st[0].data) * 0.5

            # 特殊情况加上协变导数
            if c1=='R' and c2=='T':
                channel = f"T"
                st = st_syn.select(channel=channel)
                if len(st) == 0:
                    raise NameError(f"{channel} not exists.")
                tr.data -= 0.5*st[0].data / dist * 1e-5
            
            elif c1=='T' and c2=='T':
                channel = f"R"
                st = st_syn.select(channel=channel)
                if len(st) == 0:
                    raise NameError(f"{channel} not exists.")
                tr.data += st[0].data / dist * 1e-5

            # 修改通道名
            tr.stats.channel = tr.stats.sac['kcmpnm'] = f"{c1}{c2}"

            stres.append(tr)

    return stres


def _compute_static_strain_rotation(syn:dict, Type:str):
    r"""
        Compute static strain/rotation tensor from synthetic spatial derivatives.

        :param     syn:      synthetic spatial derivatives.
        :param     Type:        "strain" or "rotation"

        :return:
            - **res** -  static strain/rotation tensor, in dict class.
    """

    if Type == 'strain':
        sgn = 1
        i1_end = 3
        i2_offset = 0
    elif Type == 'rotation':
        sgn = -1
        i1_end = 2
        i2_offset = 1
    else:
        raise ValueError(f"{Type} not supported.")

    chs = ZRTchs

    # 判断是否有标志性的分量名
    if f"nN" in syn.keys():
        chs = ZNEchs

    xarr:np.ndarray = syn['_xarr']
    yarr:np.ndarray = syn['_yarr']

    # 结果字典
    resDct = {}

    # 基本数据拷贝
    for k in syn.keys():
        if k[0] != '_':
            continue 
        resDct[k] = deepcopy(syn[k])

    # 6/3个分量建立数组
    for i1 in range(i1_end):
        c1 = chs[i1]
        for i2 in range(i1+i2_offset, 3):
            c2 = chs[i2]
            channel = f"{c1}{c2}"
            resDct[channel] = np.zeros((len(xarr), len(yarr)), dtype='f8')


    for iy in range(len(yarr)):
        for ix in range(len(xarr)):
            # 震中距
            dist = max(np.sqrt(xarr[ix]**2 + yarr[iy]**2), 1e-5)

            # ----------------------------------------------------------------------------------
            # 循环6/3个分量
            for i1 in range(i1_end):
                c1 = chs[i1]
                for i2 in range(i1+i2_offset, 3):
                    c2 = chs[i2]

                    channel = f"{c2.lower()}{c1}"
                    v12 = syn[channel][ix, iy]

                    channel = f"{c1.lower()}{c2}"
                    v21 = syn[channel][ix, iy]

                    val = 0.5*(v12 + sgn*v21)

                    # 特殊情况加上协变导数
                    if c1=='R' and c2=='T':
                        channel = f"T"
                        v0 = syn[channel][ix, iy]
                        val -= 0.5*v0 / dist * 1e-5

                    elif c1=='T' and c2=='T':
                        channel = f"R"
                        v0 = syn[channel][ix, iy]
                        val += v0 / dist * 1e-5

                    channel = f"{c1}{c2}"
                    resDct[channel][ix, iy] = val

    return resDct


def compute_strain(st:Union[Stream,dict]):
    r"""
        Compute dynamic/static strain tensor from synthetic spatial derivatives.

        :param     st:      synthetic spatial derivatives
                            :class:`obspy.Stream` class for dynamic case, dict class for static case.

        :return:
            - **stres** -  dynamic/static strain tensor, in :class:`obspy.Stream` class or dict class.
    """
    if isinstance(st, Stream):
        return _compute_strain_rotation(st, "strain")
    elif isinstance(st, dict):
        return _compute_static_strain_rotation(st, "strain")
    else:
        raise NotImplementedError
    
def compute_rotation(st:Union[Stream,dict]):
    r"""
        Compute dynamic/static rotation tensor from synthetic spatial derivatives.

        :param     st:      synthetic spatial derivatives
                            :class:`obspy.Stream` class for dynamic case, dict class for static case.

        :return:
            - **stres** -  dynamic/static rotation tensor, in :class:`obspy.Stream` class or dict class.
    """
    if isinstance(st, Stream):
        return _compute_strain_rotation(st, "rotation")
    elif isinstance(st, dict):
        return _compute_static_strain_rotation(st, "rotation")
    else:
        raise NotImplementedError


def _compute_stress(st_syn:Stream):
    r"""
        Compute dynamic stress tensor from synthetic spatial derivatives.

        :param     st_syn:      synthetic spatial derivatives.

        :return:
            - **stream** -  dynamic stress tensor (unit: dyne/cm^2 = 0.1 Pa), in :class:`obspy.Stream` class.
    """

    # 由于有Q值的存在，lambda和mu变成了复数，需在频域进行

    chs = ZRTchs
    rot2ZNE:bool = False

    # 判断是否有标志性的trace
    if len(st_syn.select(channel=f"nN")) > 0:
        chs = ZNEchs
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
        C_grt_py_attenuation_law(Qbinv, omega, atte)
        attb = atte[0] + atte[1]*1j
        mus[i] = vb*vb*attb*attb*rho*1e10
        C_grt_py_attenuation_law(Qainv, omega, atte)
        atta = atte[0] + atte[1]*1j
        lams[i] = va*va*atta*atta*rho*1e10 - 2.0*mus[i]
    
    del omega, atte

    # ----------------------------------------------------------------------------------
    # 先计算体积应变u_kk = u_11 + u22 + u33 和 lamda的乘积
    lam_ukk = np.zeros((nf,), dtype='c16')
    for i in range(3):
        c = chs[i]
        channel = f"{c.lower()}{c}"
        st = st_syn.select(channel=channel)
        if len(st) == 0:
            raise NameError(f"{channel} not exists.")
        lam_ukk[:] += rfft(st[0].data, nt)

    # 加上协变导数
    if not rot2ZNE:
        channel = f"R"
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

            channel = f"{c2.lower()}{c1}"
            st = st_syn.select(channel=channel)
            if len(st) == 0:
                raise NameError(f"{channel} not exists.")
            tr = st[0].copy()
            fftarr = np.zeros((nf,), dtype='c16')

            channel = f"{c1.lower()}{c2}"
            st = st_syn.select(channel=channel)
            if len(st) == 0:
                raise NameError(f"{channel} not exists.")
            fftarr[:] = rfft(tr.data + st[0].data, nt) * mus

            # 对于对角线分量，需加上lambda * u_kk
            if c1==c2:
                fftarr[:] += lam_ukk

            # 特殊情况加上协变导数
            if c1=='R' and c2=='T':
                channel = f"T"
                st = st_syn.select(channel=channel)
                if len(st) == 0:
                    raise NameError(f"{channel} not exists.")
                fftarr[:] -= mus*rfft(st[0].data, nt) / dist * 1e-5
            
            elif c1=='T' and c2=='T':
                channel = f"R"
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
        Compute static stress tensor from synthetic spatial derivatives.

        :param     syn:      synthetic spatial derivatives.

        :return:
            - **res** -  static stress tensor (unit: dyne/cm^2 = 0.1 Pa), in dict class.
    """
 
    chs = ZRTchs
    rot2ZNE:bool = False

    # 判断是否有标志性的分量名
    if f"nN" in syn.keys():
        chs = ZNEchs
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
                channel = f"{c.lower()}{c}"
                lam_ukk += syn[channel][ix, iy]
            
            # 加上协变导数
            if not rot2ZNE:
                channel = f"R"
                lam_ukk += syn[channel][ix, iy] / dist * 1e-5
            
            lam_ukk *= lam

            # ----------------------------------------------------------------------------------
            # 循环6个分量
            for i1 in range(3):
                c1 = chs[i1]
                for i2 in range(i1, 3):
                    c2 = chs[i2]

                    channel = f"{c2.lower()}{c1}"
                    v12 = syn[channel][ix, iy]

                    channel = f"{c1.lower()}{c2}"
                    v21 = syn[channel][ix, iy]

                    val = mu*(v12 + v21)

                    # 对于对角线分量，需加上lambda * u_kk
                    if c1==c2:
                        val += lam_ukk

                    # 特殊情况加上协变导数
                    if c1=='R' and c2=='T':
                        channel = f"T"
                        v0 = syn[channel][ix, iy]
                        val -= mu*v0 / dist * 1e-5

                    elif c1=='T' and c2=='T':
                        channel = f"R"
                        v0 = syn[channel][ix, iy]
                        val += 2.0*mu*v0 / dist * 1e-5

                    channel = f"{c1}{c2}"
                    resDct[channel][ix, iy] = val

    return resDct


def compute_stress(st:Union[Stream,dict]):
    r"""
        Compute dynamic/static stress tensor from synthetic spatial derivatives.

        :param     st:      synthetic spatial derivatives
                            :class:`obspy.Stream` class for dynamic case, dict class for static case.

        :return:
            - **stres** -  dynamic/static stress tensor (unit: dyne/cm^2 = 0.1 Pa), in :class:`obspy.Stream` class or dict class.
    """
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
        convolve each trace with a signal

        :param    st0:        :class:`obspy.Stream`
        :param    signal0:    convolution signal
        :param    inplace:    whether change in-place  

        :return:
            - **stream** -    convolution result, :class:`obspy.Stream`
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
        Perform integration on each trace
        
        :param    st0:        :class:`obspy.Stream`
        :param    inplace:    whether change in-place  

        :return:
            - **stream** -    integration result, :class:`obspy.Stream`
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
        Perform central difference on each trace

        :param    st0:        :class:`obspy.Stream`
        :param    inplace:    whether change in-place  

        :return:
            - **stream** -    difference result, :class:`obspy.Stream`
    '''
    st = st0 if inplace else deepcopy(st0)
    
    for tr in st:
        data = tr.data 
        data[:] = np.gradient(data, tr.stats.delta)

    return st


def stream_write_sac(st:Stream, dir:str):
    '''
        save each trace to "dir/{channel}.sac"

        :param    st:         :class:`obspy.Stream`
        :param    dir:        saving directory

    '''
    # 新建对应文件夹
    os.makedirs(dir, exist_ok=True)

    # 每一道的保存路径为 dir/{channel}.sac
    for tr in st:
        filepath = os.path.join(dir, f"{tr.stats.channel}.sac")
        tr.write(filepath, format='SAC')





#=================================================================================================================
#
#                                           积分过程文件读取及绘制
#
#=================================================================================================================


def read_statsfile(statsfile:str):
    '''
        read a statsfile  

        :param    statsfile:       File path (Wildcards can be used to simplify input)

        :return:
            - **data** -     `numpy.ndarray <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>`_ custom type array 
    '''
    Lst = glob.glob(statsfile)
    if len(Lst) != 1:
        raise OSError(f"{statsfile} should only match one file, but {len(Lst)} matched.")
    statsfile = Lst[0]
    print(f"raed in {statsfile}.")

    basename = os.path.basename(statsfile)

    # 确定自定义数据类型  EX_q, EX_w, VF_q, ...
    dtype = [('k' if basename[0] == 'K' else 'c', NPCT_REAL_TYPE)]
    for im in range(SRC_M_NUM):
        modr = SRC_M_ORDERS[im]
        for c in range(QWV_NUM):
            if modr==0 and qwvchs[c] == 'v':
                continue 

            dtype.append((f"{SRC_M_NAME_ABBR[im]}_{qwvchs[c]}", NPCT_CMPLX_TYPE))


    data = np.fromfile(statsfile, dtype=dtype)

    return data


def read_kernels_freqs(statsdir:str, vels:Union[np.ndarray,None]=None, ktypes:Union[List[str],None]=None):
    r"""
        read all statsfiles in statsdir (except that of 0 frequency).
        If record wavenumber, interpolate to the phase velocity.

        :param        statsdir:     directory path
        :param        vels:         When a positive-order vels (km/s) is specified, files starting with `K_` are read 
                                    and linear interpolation from wavenumber to phase velocity is performed.
                                    Otherwise read the files starting with `C_`
        :param        ktype:        Specify the return of a series of kernel function names, 
                                    such as `EX_q`, `DS_w`, etc. By default, all are returned

        :return:
            - **kerDct**  -   kernel functions in a dict
    """

    dointerp = vels is not None

    if (dointerp) and not np.all(np.diff(vels) > 0):
        raise ValueError("vels must be in ascending order.")
    
    K_statspaths = glob.glob(os.path.join(statsdir, "K_*"))
    if len(K_statspaths) == 0 and dointerp:
        raise ValueError("You want to interpolate from k to c, but found 0 statsfiles recording k.")
    
    C_statspaths = glob.glob(os.path.join(statsdir, "C_*"))
    if len(C_statspaths) == 0 and not dointerp:
        raise ValueError("Found 0 statsfiles directly recording c.")
    
    statspaths = K_statspaths if dointerp else C_statspaths

    KLst = np.array(statspaths)
    freqs = np.array([float(s.split("_")[-1]) for s in KLst])
    # 根据freqs排序
    _idx = np.argsort(freqs)
    freqs[:] = freqs[_idx]
    KLst[:] = KLst[_idx]
    del _idx 

    # 去除零频
    if freqs[0] == 0.0:
        freqs = freqs[1:]
        KLst = KLst[1:]

    kerDct = {}
    kerDct['_vels'] = vels.copy() if dointerp else []
    kerDct['_freqs'] = freqs.copy()

    for i in range(len(freqs)):
        Kpath = KLst[i]
        freq = freqs[i]
        w = 2*np.pi*freq

        data = read_statsfile(Kpath)
        
        if dointerp:
            v = w/data['k']

            # 检查v范围
            v1 = np.min(v)
            v2 = np.max(v)
            if v1 > vels[0] or v2 < vels[-1]:
                raise ValueError(f"In freq={freq:.5e}, minV={v1:.5e}, maxV={v2:.5e}, insufficient wavenumber samples"
                                " to interpolate on vels.")
        else:
            if len(kerDct['_vels']) == 0:
                kerDct['_vels'] = data['c'].copy()

        for key in data.dtype.names:
            if key == 'k' or key == 'c':
                continue 
            if (ktypes is not None) and (key not in ktypes):
                continue 

            if key not in kerDct.keys():
                kerDct[key] = []

            if dointerp:
                # 如果越界会报错
                F = interpn((v,), data[key], vels)
                kerDct[key].append(F)
            else:
                kerDct[key].append(data[key])

    # 将每个核函数结果拼成2D数组
    for key in kerDct.keys():
        if key[0] == '_':
            continue
        kerDct[key] = np.vstack(kerDct[key])

    return kerDct


def read_statsfile_ptam(statsfile:str):
    '''
        read a statsfile from PTAM process  

        :param    statsfile:       File path (Wildcards can be used to simplify input)

        :return:
            - **data1** -     `numpy.ndarray <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>`_ custom type array, during DCM or (SA)FIM
            - **data2** -     `numpy.ndarray <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>`_ custom type array, during PTAM
            - **ptam_data** -   `numpy.ndarray <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>`_ custom type array, record the peak/trough from PTAM
            - **dist** -      epicentral distance from the filename (km)
    '''
    Lst = glob.glob(statsfile)
    if len(Lst) != 1:
        raise OSError(f"{statsfile} should only match one file, but {len(Lst)} matched.")
    statsfile = Lst[0]

    # 获得震中距
    dist = float(os.path.dirname(statsfile).split("_")[-1])

    # 从文件路径命名中，获得对应的K文件路径
    PTAMname = os.path.basename(statsfile)
    if "_" in PTAMname:  # 动态解
        splits = PTAMname.split("_")
        splits[-3] = "K"
        K_basename= "_".join(splits)
    else:
        K_basename = "K" # 静态解
        
    data1 = read_statsfile(os.path.join(os.path.dirname(os.path.dirname(statsfile)), K_basename))
    data2 = read_statsfile(os.path.join(os.path.dirname(statsfile), K_basename))

    # 确定自定义数据类型  sum_EX_0_k, sum_EX_0, sum_VF_0_k, ...
    # 各格林函数数值积分的值(k上限位于不同的波峰波谷)
    # 开头的sum表示这是波峰波谷位置处的数值积分的值(不含dk)，
    # 末尾的k表示对应积分值的波峰波谷位置的k值
    dtype = []
    for im in range(SRC_M_NUM):
        modr = SRC_M_ORDERS[im]
        for v in range(INTEG_NUM):
            if modr==0 and v!=0 and v!=2:
                continue 

            dtype.append((f"sum_{SRC_M_NAME_ABBR[im]}_{v}_k", NPCT_REAL_TYPE))
            dtype.append((f"sum_{SRC_M_NAME_ABBR[im]}_{v}", NPCT_CMPLX_TYPE))


    ptam_data = np.fromfile(statsfile, dtype=dtype)

    return data1, data2, ptam_data, dist



def _get_stats_Fname(statsdata:np.ndarray, karr:np.ndarray, dist:float, srctype:str, ptype:str):
    # 根据ptype获得对应的核函数
    krarr = karr*dist

    # 从数组中找到震源名称的索引
    try:
        _idx = SRC_M_NAME_ABBR.index(srctype)
        mtype = str(SRC_M_ORDERS[_idx])
    except:
        raise ValueError(f"{srctype} is an invalid name.")

    if mtype=='0':
        if ptype=='0':
            Fname = rf"$F(k,\omega)=q^{{({srctype})}}(k, \omega)$"
            Farr = statsdata[f'{srctype}_q']
            FJname = rf"$ - F(k,\omega)J_1(kr)k$"
            FJarr =  - jv(1, krarr) * Farr * karr
        elif ptype=='2':
            Fname = rf"$F(k,\omega)=w^{{({srctype})}}(k, \omega)$"
            FJname = rf"$F(k,\omega)J_0(kr)k$"
            Farr = statsdata[f'{srctype}_w']
            FJarr = jv(0, krarr) * Farr * karr
        else:
            raise ValueError(f"source {srctype}, m={mtype}, p={ptype} is not supported.")
        
    elif mtype in ['1', '2']:
        m = int(mtype)
        if ptype=='0':
            Fname = rf"$F(k,\omega)=q^{{({srctype})}}(k, \omega)$"
            Farr = statsdata[f'{srctype}_q']
            FJname = rf"$F(k,\omega)J_{m-1}(kr)k$"
            FJarr = jv(m-1, krarr) * Farr * karr
        elif ptype=='1':
            Fname = rf"$F(k,\omega)=q^{{({srctype})}}(k, \omega) + v^{{({srctype})}}(k, \omega)$"
            Farr = (statsdata[f'{srctype}_q'] + statsdata[f'{srctype}_v'])
            FJname = rf"$ - F(k,\omega) \dfrac{{{m}}}{{kr}} J_{m}(kr)k$"
            FJarr =  - jv(m, krarr) * Farr * m/dist
        elif ptype=='2':
            Fname = rf"$F(k,\omega)=w^{{({srctype})}}(k, \omega)$"
            Farr = statsdata[f'{srctype}_w']
            FJname = rf"$F(k,\omega)J_{m}(kr)k$"
            FJarr = jv(m, krarr) * Farr * karr
        elif ptype=='3':
            Fname = rf"$F(k,\omega)=v^{{({srctype})}}(k, \omega)$"
            Farr = statsdata[f'{srctype}_v']
            FJname = rf"$ - F(k,\omega)J_{m-1}(kr)k$"
            FJarr =  - jv(m-1, krarr) * Farr * karr
        else:
            raise ValueError(f"source {srctype}, m={mtype}, p={ptype} is not supported.")
        
    else:
        raise ValueError(f"source {srctype}, m={mtype}, p={ptype} is not supported.")
    
    return Fname, Farr, FJname, FJarr


def plot_statsdata(statsdata:np.ndarray, dist:float, srctype:str, ptype:str, RorI:Union[bool,int]=True,
                   fig:Union[Figure,None]=None, axs:Union[Axes,None]=None):
    r'''
        Based on the data read by the :func:`read_statsfile <pygrt.utils.read_statsfile>` function,
        plot the kernel function :math:`F(k,\omega)`, the integrand :math:`F(k,\omega)J_m(kr)k`, 
        and calculate the cumulative integral :math:`\sum F(k,\omega)J_m(kr)k` .

        .. note:: Not every source type corresponds to every order and every integration type, see :ref:`grn_types` for details.

        :param    statsdata:         return value of :func:`read_statsfile <pygrt.utils.read_statsfile>` function
        :param    dist:              epicentral distance (km)
        :param    srctype:           abbreviation of source type, including EX, VF, HF, DD, DS, SS
        :param    ptype:             integration type (0,1,2,3)
        :param    RorI:              whether to plot real or imaginary part, default is real part, pass 2 to plot both
        :param    fig:               user-defined matplotlib.Figure object, default is None
        :param    axs:               user-defined matplotlib.Axes object array (three elements), default is None

        :return:
                - **fig** -                        matplotlib.Figure object
                - **(ax1,ax2,ax3)** -              matplotlib.Axes object array
    '''

    ptype = str(ptype)

    karr = statsdata['k'] 
    dk = (karr[1] - karr[0])   # 假设均匀dk
    is_evendk = np.allclose(np.diff(karr), dk, atol=1e-10)  # 是否为均匀dk
    if not is_evendk:
        raise ValueError("Sorry, this function only supports even-distributed k.")
    
    if 0.5*np.pi/dk < dist:  # 对于bessel函数这种震荡函数，假设一个周期内至少取4个点
        print(f"WARNING! dist ({dist}) > PI/(2*dk) ({0.5*np.pi/dk:.5e}.)")

    Fname, Farr, FJname, FJarr = _get_stats_Fname(statsdata, karr, dist, srctype, ptype)
    
    if fig is None or axs is None:
        fig, axs = plt.subplots(3, 1, figsize=(8, 9), gridspec_kw=dict(hspace=0.7))
    
    # axs长度必须为三个
    if len(axs) != 3:
        raise ValueError("axs should have 3 elements.")

    ax1, ax2, ax3 = axs

    if isinstance(RorI, int) and RorI==2:
        ax1.plot(karr, np.real(Farr), lw=0.8, label='Real') 
        ax1.plot(karr, np.imag(Farr), lw=0.8, label='Imag') 
    else:
        if RorI:
            ax1.plot(karr, np.real(Farr), lw=0.8, label='Real') 
        else:
            ax1.plot(karr, np.imag(Farr), lw=0.8, label='Imag') 

    ax1.set_xlabel('k /$km^{-1}$')
    ax1.set_title(Fname)
    ax1.grid()
    ax1.legend(loc='lower left')

    if isinstance(RorI, int) and RorI==2:
        ax2.plot(karr, np.real(FJarr), lw=0.8, label='Real') 
        ax2.plot(karr, np.imag(FJarr), lw=0.8, label='Imag') 
    else:
        if RorI:
            ax2.plot(karr, np.real(FJarr), lw=0.8, label='Real') 
        else:
            ax2.plot(karr, np.imag(FJarr), lw=0.8, label='Imag') 
    ax2.set_title(FJname)
    ax2.set_xlabel('k /$km^{-1}$')
    ax2.grid()
    ax2.legend(loc='lower left')

    # 数值积分，不乘系数dk 
    Parr = np.cumsum(FJarr)

    if isinstance(RorI, int) and RorI==2:
        ax3.plot(karr, np.real(Parr), lw=0.8, label='Real') 
        ax3.plot(karr, np.imag(Parr), lw=0.8, label='Imag') 
    else:
        if RorI:
            ax3.plot(karr, np.real(Parr), lw=0.8, label='Real') 
        else:
            ax3.plot(karr, np.imag(Parr), lw=0.8, label='Imag') 
    ax3.set_title(f'$\sum_k$ {FJname}')
    ax3.set_xlabel("k /$km^{-1}$")
    ax3.grid()
    ax3.legend(loc='lower left')

    return fig, (ax1, ax2, ax3)


def plot_statsdata_ptam(statsdata1:np.ndarray, statsdata2:np.ndarray, statsdata_ptam:np.ndarray, 
                        dist:float, srctype:str, ptype:str, RorI:Union[bool,int]=True,
                        fig:Union[Figure,None]=None, axs:Union[Axes,None]=None):
    r'''
        Based on data read by the :func:`read_statsfile_ptam <pygrt.utils.read_statsfile_ptam>` function,
        simply calculate and plot the cumulative integral as well as the peak/trough positions used by PTAM.

        .. note:: Not every source type corresponds to every order and every integration type, see :ref:`grn_types` for details.

        :param    statsdata1:        integral process data during DWM or FIM
        :param    statsdata2:        integral process data during PTAM
        :param    statsdata_ptam:    peak/trough positions and amplitudes from PTAM
        :param    dist:              epicentral distance (km)
        :param    srctype:           abbreviation of source type, including EX, VF, HF, DD, DS, SS  
        :param    ptype:             integration type (0, 1, 2, 3)
        :param    RorI:              whether to plot real or imaginary part, default is real part, pass 2 to plot both
        :param    fig:               user-defined matplotlib.Figure object, default is None
        :param    axs:               user-defined matplotlib.Axes object array (three elements), default is None

        :return:  
                - **fig** -                        matplotlib.Figure object   
                - **(ax1, ax2, ax3)** -            matplotlib.Axes object array
    '''

    ptype = str(ptype)

    karr1 = statsdata1['k'] 
    dk1 = karr1[1] - karr1[0]
    Fname, Farr1, FJname, FJarr1 = _get_stats_Fname(statsdata1, karr1, dist, srctype, ptype)
    karr2 = statsdata2['k'] 
    dk2 = karr2[1] - karr2[0]
    Fname, Farr2, FJname, FJarr2 = _get_stats_Fname(statsdata2, karr2, dist, srctype, ptype)

    is_evendk = np.allclose(np.diff(karr1), dk1, atol=1e-10) and np.allclose(np.diff(karr2), dk2, atol=1e-10)  # 是否为均匀dk
    if not is_evendk:
        raise ValueError("Sorry, this function only supports even-distributed k.")

    # 将两个过程的结果拼起来
    Farr = np.hstack((Farr1, Farr2))
    karr = np.hstack((karr1, karr2))
    FJarr = np.hstack((FJarr1, FJarr2))

    if fig is None or axs is None:
        fig, axs = plt.subplots(3, 1, figsize=(8, 9), gridspec_kw=dict(hspace=0.7))
    
    # axs长度必须为三个
    if len(axs) != 3:
        raise ValueError("axs should have 3 elements.")

    ax1, ax2, ax3 = axs

    if isinstance(RorI, int) and RorI==2:
        ax1.plot(karr, np.real(Farr), lw=0.8, label='Real') 
        ax1.plot(karr, np.imag(Farr), lw=0.8, label='Imag') 
    else:
        if RorI:
            ax1.plot(karr, np.real(Farr), lw=0.8, label='Real') 
        else:
            ax1.plot(karr, np.imag(Farr), lw=0.8, label='Imag') 

    ax1.set_xlabel('k /$km^{-1}$')
    ax1.set_title(Fname)
    ax1.grid()
    ax1.legend(loc='lower left')

    if isinstance(RorI, int) and RorI==2:
        ax2.plot(karr, np.real(FJarr), lw=0.8, label='Real') 
        ax2.plot(karr, np.imag(FJarr), lw=0.8, label='Imag') 
    else:
        if RorI:
            ax2.plot(karr, np.real(FJarr), lw=0.8, label='Real') 
        else:
            ax2.plot(karr, np.imag(FJarr), lw=0.8, label='Imag') 
    ax2.set_title(FJname)
    ax2.set_xlabel('k /$km^{-1}$')
    ax2.grid()
    ax2.legend(loc='lower left')

    # 波峰波谷位置，用红十字标记
    ptKarr = statsdata_ptam[f'sum_{srctype}_{ptype}_k']
    ptFJarr = statsdata_ptam[f'sum_{srctype}_{ptype}']

    # 数值积分，不乘系数dk 
    Parr1 = np.cumsum(FJarr1) 
    Parr2 = np.cumsum(FJarr2)  
    Parr = np.hstack([Parr1, Parr2*dk2/dk1+Parr1[-1]])

    if isinstance(RorI, int) and RorI==2:
        ax3.plot(karr, np.real(Parr), lw=0.8, label='Real') 
        ax3.plot(ptKarr, np.real(ptFJarr), 'r+', markersize=6)
        ax3.plot(karr, np.imag(Parr), lw=0.8, label='Imag') 
        ax3.plot(ptKarr, np.imag(ptFJarr), 'r+', markersize=6)
    else:
        if RorI:
            ax3.plot(karr, np.real(Parr), lw=0.8, label='Real') 
            ax3.plot(ptKarr, np.real(ptFJarr), 'r+', markersize=6)
        else:
            ax3.plot(karr, np.imag(Parr), lw=0.8, label='Imag') 
            ax3.plot(ptKarr, np.imag(ptFJarr), 'r+', markersize=6)
    

    ax3.set_title(f'$\sum_k$ {FJname}')
    ax3.set_xlabel("k /$km^{-1}$")
    ax3.grid()
    ax3.legend(loc='lower left')


    return fig, (ax1, ax2, ax3)






def solve_lamb1(nu:float, ts:np.ndarray, azimuth:float):
    r"""
        solve the first-kind Lamb's problem using the generalized closed-form solution, see：

            张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社

        :param      nu:         possion ratio in (0, 0.5)
        :param      ts:         dimensionless time :math:`\bar{t}` ，:math:`\bar{t}=\dfrac{t}{r/\beta}`
        :param      azimuth:    azimuth in degree

        :return:    Normalized solution with shape (nt, 3, 3). To get the physical solution, 
                    divide by :math:`\pi^2 \mu r`
    """

    # 检查数据
    if np.any(ts < 0.0):
        raise ValueError("ts should be nonnegative.")
    if azimuth < 0.0 or azimuth > 360.0:
        raise ValueError("azimuth should be in [0, 360].")
    
    ts = np.array(ts).astype(NPCT_REAL_TYPE)
    
    # 定义结果数组
    nt = len(ts)
    u = np.zeros((nt, 3, 3), dtype=NPCT_REAL_TYPE)

    C_grt_solve_lamb1(nu, npct.as_ctypes(ts), nt, azimuth, npct.as_ctypes(u.ravel()))

    return u

