"""
    :file:     pymod.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包括 Python端使用的模型 :class:`pygrt.c_structures.c_PyModel1D`

"""


from __future__ import annotations
from multiprocessing import Value
import numpy as np
import numpy.ctypeslib as npct
from obspy import read, Stream, Trace, UTCDateTime
from scipy.fft import irfft, ifft
from obspy.core import AttribDict
from typing import List, Dict

from time import time
from copy import deepcopy

from ctypes import *
from .c_interfaces import *
from .c_structures import *
from .pygrn import PyGreenFunction


__all__ = [
    "PyModel1D",
]


class PyModel1D:
    def __init__(self, modarr0:np.ndarray, depsrc:float, deprcv:float):
        '''
            将震源和台站插入定义模型的数组，转为 :class:`PyModel1D <pygrt.pymod.PyModel1D>` 实例的形式  

            :param    modarr0:    模型数组，每行格式为[thickness(km), Vp(km/s), Vs(km/s), Rho(g/cm^3), Qp, Qs]  
            :param    depsrc:     震源深度(km)  
            :param    deprcv:     台站深度(km)  

        '''
        self.modarr:np.ndarray
        self.depsrc = depsrc 
        self.deprcv = deprcv 
        self.c_pymod1d:c_PyModel1D 

        if depsrc < 0:
            raise ValueError(f"depsrc ({depsrc}) < 0")
        if deprcv < 0:
            raise ValueError(f"deprcv ({deprcv}) < 0")


        modarr = modarr0.copy()

        # 最后一层，本质只为保证震源和台站能插入模型中
        modarr[-1, 0] = np.max([modarr[-1, 0], 9e10, deprcv, depsrc]) + 1.0 

        # 将震源和台站的虚拟层位插入模型
        dep = 0.0
        ircv = 0
        for i in range(modarr.shape[0]):
            dep += modarr[i,0]
            if dep > deprcv:
                ircv = i+1
                lay = np.copy(modarr[i,:])
                lay[0] = dep - deprcv
                modarr[i,0] -= lay[0]
                modarr = np.insert(modarr, ircv, lay, axis=0)
                break   

        dep = 0.0
        isrc = 0 
        for i in range(modarr.shape[0]):
            dep += modarr[i,0]
            if dep > depsrc:
                isrc = i+1
                lay = np.copy(modarr[i,:])
                lay[0] = dep - depsrc
                modarr[i,0] -= lay[0]
                modarr = np.insert(modarr, isrc, lay, axis=0)
                break 

        # 如果台站位于震源深度以下，需要调整层索引; 因为先插入的台站，再插入的震源 
        if depsrc < deprcv:
            ircv += 1

        # 调整层厚，拒绝0厚度层 
        # for i in range(modarr.shape[0]):
        #     if modarr[i, 0] == 0.0:
        #         modarr[i, 0] = 1e-5

            
        c_pymod1d = c_PyModel1D(
            modarr.shape[0],
            depsrc,
            deprcv,
            isrc, 
            ircv, 
            (ircv<isrc),

            npct.as_ctypes(modarr[:,0].astype(NPCT_REAL_TYPE)),
            npct.as_ctypes(modarr[:,1].astype(NPCT_REAL_TYPE)),
            npct.as_ctypes(modarr[:,2].astype(NPCT_REAL_TYPE)),
            npct.as_ctypes(modarr[:,3].astype(NPCT_REAL_TYPE)),
            npct.as_ctypes(modarr[:,4].astype(NPCT_REAL_TYPE)),
            npct.as_ctypes(modarr[:,5].astype(NPCT_REAL_TYPE)),
        )

        self.modarr = modarr
        self.c_pymod1d = c_pymod1d

        self.vmax = np.max(self.modarr[:, 1:3])
        self.vmin = np.min(self.modarr[:, 1:3])
        if self.vmin <= 0.0:
            raise ValueError("Zero Velocity ??")
        
    
    def compute_travt1d(self, dist:float):
        travtP = C_compute_travt1d(
            self.c_pymod1d.Thk,
            self.c_pymod1d.Va,
            self.c_pymod1d.n,
            self.c_pymod1d.isrc,
            self.c_pymod1d.ircv,
            dist
        )
        travtS = C_compute_travt1d(
            self.c_pymod1d.Thk,
            self.c_pymod1d.Vb,
            self.c_pymod1d.n,
            self.c_pymod1d.isrc,
            self.c_pymod1d.ircv,
            dist
        )

        return travtP, travtS


    def gen_gf_spectra(
        self, 
        distarr:np.ndarray|list[float]|float, 
        nt:int, 
        dt:float, 
        freqband:np.ndarray|list[float]=[-1,-1],
        zeta:float=0.8, 
        vmin_ref:float=0.0,
        keps:float=-1.0,  
        ampk:float=1.15,
        iwk0:bool=False,
        k0:float=5.0, 
        Length:float=0.0, 
        delayT0:float=0.0,
        delayV0:float=0.0,
        gf_source=['EXP', 'VF', 'HF', 'DC'],
        statsfile:str|None=None, 
        statsidxs:np.ndarray|list|None=None, 
        print_runtime:bool=True):
        
        r'''
            
            调用C库计算格林函数的主函数，以列表的形式返回，其中每个元素为对应震中距的格林函数频谱字典。
            

            :param    distarr:       多个震中距(km) 的数组, 或单个震中距的浮点数
            :param    nt:            时间点数，借助于 `SciPy`，nt不再要求是2的幂次
            :param    dt:            采样间隔(s)  
            :param    freqband:      频率范围(Hz)，以此确定待计算的离散频率点
            :param    zeta:          定义虚频率的系数 :math:`\zeta`， 虚频率 :math:`\tilde{\omega} = \omega - j*w_I, w_I = \zeta*\pi/T, T=nt*dt` , T为时窗长度。
                                     使用离散波数积分时为了避开附加源以及奇点的影响， :ref:`(Bouchon, 1981) <bouchon_1981>`  在频率上添加微小虚部，
                                     更多测试见 :ref:`(张海明, 2021) <zhang_book_2021>`
            :param    vmin_ref:      最小参考速度，默认vmin=max(minimum velocity, 0.1)，用于定义波数积分上限，小于0则在达到积分上限后使用峰谷平均法
                                    （默认当震源和场点深度差<=0.5km时自动使用峰谷平均法）
            :param    keps:          波数k积分收敛条件，见 :ref:`(Yao and Harkrider, 1983) <yao&harkrider_1983>`  :ref:`(初稿) <yao_init_manuscripts>`，
                                     为负数代表不提前判断收敛，按照波数积分上限进行积分
            :param    ampk:          影响波数k积分上限的系数，见下方
            :param    iwk0:          k0是否取随频率变化的线性关系，即 :math:` k_{0} = k_{0} * f/f_{max}` 
            :param    k0:            波数k积分的上限 :math:`\tilde{k_{max}}=\sqrt{(k_{0}*\pi/hs)^2 + (ampk*w/vmin_{ref})^2} ` , 波数k积分循环必须退出, hs=max(震源和台站深度差,1.0)
            :param    Length:        定义波数k积分的间隔 `dk=2\pi / (L*rmax)`, 选取要求见 :ref:`(Bouchon, 1981) <bouchon_1981>` 
                                     :ref:`(张海明, 2021) <zhang_book_2021>`，默认自动选择
            :param    gf_source:     待计算的震源类型
            :param    statsfile:     波数k积分（包括Filon积分和峰谷平均法）的过程记录文件，常用于debug或者观察积分过程中 :math:`F(k,\omega)` 和  :math:`F(k,\omega)J_m(kr)k` 的变化    
            :param    statsidxs:     仅输出特定频点的过程记录文件，建议给定频点，否则默认所有频率点的记录文件都输出，很占空间
            :param    print_runtime: 是否打印运行时间

            :return:
                - **gfspecList** -   列表，每个元素以字典的形式存有每个震中距的格林函数频谱(:class:`pygrt.pygrn.PyGreenFunction` ) 实例)
                
        '''

        depsrc = self.depsrc
        deprcv = self.deprcv

        calc_EXP:bool = 'EXP' in gf_source
        calc_VF:bool = 'VF' in gf_source
        calc_HF:bool = 'HF' in gf_source
        calc_DC:bool = 'DC' in gf_source

        if isinstance(distarr, list):
            distarr = np.array(distarr)
        elif isinstance(distarr, float) or isinstance(distarr, int):
            distarr = np.array([distarr*1.0]) 

        distarr = distarr.copy().astype(NPCT_REAL_TYPE)

        if np.any(distarr < 0):
            raise ValueError(f"distarr < 0")
        if nt < 0:
            raise ValueError(f"nt ({nt}) < 0")
        if dt < 0:
            raise ValueError(f"dt ({dt}) < 0")
        if zeta < 0:
            raise ValueError(f"zeta ({zeta}) < 0")
        if k0 < 0:
            raise ValueError(f"k0 ({k0}) < 0")

        nf = nt//2+1 
        df = 1/(nt*dt)
        fnyq = 1/(2*dt)
        # 确定频带范围 
        f1, f2 = freqband 
        if f1 >= f2 and f1 >= 0 and f2 >= 0:
            raise ValueError(f"freqband f1({f1}) >= f2({f2})")
        
        if f1 < 0:
            f1 = 0 
        if f2 < 0:
            f2 = fnyq+df
            
        f1 = max(0, f1) 
        f2 = min(f2, fnyq + df)
        nf1 = max(0, int(np.floor(f1/df)))
        nf2 = min(int(np.ceil(f2/df)), nf-1) 

        # 所有频点 
        freqs = (np.arange(0, nf)*df).astype(NPCT_REAL_TYPE) 

        # 虚频率 
        wI = zeta * np.pi/(nt*dt)

        # 避免绝对0震中距 
        nrs = len(distarr)
        for ir in range(nrs):
            if(distarr[ir] < 0.0):
                raise ValueError(f"r({distarr[ir]}) < 0")
            elif(distarr[ir] == 0.0):
                distarr[ir] = 1e-5 

        # 最大震中距
        rmax = np.max(distarr)
        
        # 转为C类型
        c_freqs = npct.as_ctypes(freqs)
        c_rs = npct.as_ctypes(np.array(distarr).astype(NPCT_REAL_TYPE) )

        # 参考最小速度
        if vmin_ref == 0.0:
            vmin_ref = max(self.vmin, 0.1)
            if abs(depsrc - deprcv) <= 0.5:
                vmin_ref = - abs(vmin_ref)  # 自动使用PTAM


        # 时窗长度
        winT = nt*dt 
        
        # 时窗最大截止时刻 
        tmax = delayT0 + winT
        if delayV0 > 0.0:
            tmax += rmax/delayV0

        # 设置波数积分间隔
        # 自动情况下给出保守值
        if Length == 0.0:
            Length = 15.0
            jus = (self.vmax*tmax)**2 - (depsrc - deprcv)**2
            if jus >= 0.0:
                Length = 1.0 + np.sqrt(jus)/rmax + 0.5  # 0.5作保守值
                if Length < 15.0:
                    Length = 15.0

            print(f"Length={Length:.2f}")

        # 初始化格林函数C结构体
        C_EXPgrn = ((c_PGRN*2)*nrs)() if calc_EXP else None
        C_VFgrn = ((c_PGRN*2)*nrs)() if calc_VF else None
        C_HFgrn = ((c_PGRN*3)*nrs)() if calc_HF else None
        C_DDgrn = ((c_PGRN*2)*nrs)() if calc_DC else None
        C_DSgrn = ((c_PGRN*3)*nrs)() if calc_DC else None
        C_SSgrn = ((c_PGRN*3)*nrs)() if calc_DC else None


        EXPgrn:List[List[PyGreenFunction]] = []
        VFgrn:List[List[PyGreenFunction]] = []
        DDgrn:List[List[PyGreenFunction]] = []
        HFgrn:List[List[PyGreenFunction]] = []
        DSgrn:List[List[PyGreenFunction]] = []
        SSgrn:List[List[PyGreenFunction]] = []

        
        for ir in range(nrs):
            dist = distarr[ir]
            EXPgrn.append([])
            VFgrn .append([])
            DDgrn .append([])
            HFgrn .append([])
            DSgrn .append([])
            SSgrn .append([])
            for i, comp in enumerate(['Z', 'R', 'T']):
                if i<2:
                    if calc_EXP:
                        grn = PyGreenFunction(f'EX{comp}', nt, dt, freqs, wI, dist, depsrc, deprcv)
                        EXPgrn[ir].append(grn)
                        C_EXPgrn[ir][i] = pointer(grn.c_grn)

                    if calc_VF:
                        grn = PyGreenFunction(f'VF{comp}', nt, dt, freqs, wI, dist, depsrc, deprcv)
                        VFgrn[ir].append(grn)
                        C_VFgrn[ir][i] = pointer(grn.c_grn)
                    
                    if calc_DC:
                        grn = PyGreenFunction(f'DD{comp}', nt, dt, freqs, wI, dist, depsrc, deprcv)
                        DDgrn[ir].append(grn)
                        C_DDgrn[ir][i] = pointer(grn.c_grn)
                
                if calc_HF:
                    grn = PyGreenFunction(f'HF{comp}', nt, dt, freqs, wI, dist, depsrc, deprcv)
                    HFgrn[ir].append(grn)
                    C_HFgrn[ir][i] = pointer(grn.c_grn)

                if calc_DC:
                    grn = PyGreenFunction(f'DS{comp}', nt, dt, freqs, wI, dist, depsrc, deprcv)
                    DSgrn[ir].append(grn)
                    C_DSgrn[ir][i] = pointer(grn.c_grn)
                    grn = PyGreenFunction(f'SS{comp}', nt, dt, freqs, wI, dist, depsrc, deprcv)
                    SSgrn[ir].append(grn)
                    C_SSgrn[ir][i] = pointer(grn.c_grn) 
                


        c_statsfile = None 
        if statsfile is not None:
            c_statsfile = c_char_p(statsfile.encode('utf-8'))

        nstatsidxs = 0 
        if statsidxs is not None:
            c_statsidxs = npct.as_ctypes(np.array(statsidxs).astype('i'))
            nstatsidxs = len(statsidxs)
        else:
            c_statsidxs = None


        # ===========================================
        # 打印参数设置 
        if print_runtime:
            print(f"vmin={self.vmin}")
            print(f"vmax={self.vmax}")
            print(f"vmin_ref={abs(vmin_ref)}", end="")
            if vmin_ref < 0.0:
                print(", using PTAM.")
            else:
                print("")
            print(f"Length={abs(Length)}", end="")
            if Length < 0.0:
                print(", using FIM.")
            else:
                print("")
            print(f"nt={nt}")
            print(f"dt={dt}")
            print(f"winT={winT}")
            print(f"zeta={zeta}")
            print(f"delayT0={delayT0}")
            print(f"delayV0={delayV0}")
            print(f"tmax={tmax}")
            print(f"k0={k0}")
            print(f"ampk={ampk}")
            print(f"iwk0={iwk0}")
            print(f"keps={keps}")
            print(f"maxfreq(Hz)={freqs[nf-1]}")
            print(f"f1(Hz)={freqs[nf1]}")
            print(f"f2(Hz)={freqs[nf2]}")
            print(f"distances(km)=", distarr)
            if nstatsidxs > 0:
                print(f"statsfile_index=", statsidxs)



        # 运行C库函数
        #/////////////////////////////////////////////////////////////////////////////////
        # 计算得到的格林函数的单位：
        #     单力源 HF[ZRT],VF[ZR]                  1e-15 cm/dyne
        #     爆炸源 EX[ZR]                          1e-20 cm/(dyne*cm)
        #     剪切源 DD[ZR],DS[ZRT],SS[ZRT]          1e-20 cm/(dyne*cm)
        #=================================================================================
        C_integ_grn_spec(
            self.c_pymod1d, nf1, nf2, nf, c_freqs, nrs, c_rs, wI, 
            vmin_ref, keps, ampk, iwk0, k0, Length, print_runtime,
            C_EXPgrn, C_VFgrn, C_HFgrn, C_DDgrn, C_DSgrn, C_SSgrn, 
            c_statsfile, nstatsidxs, c_statsidxs
        )
        #=================================================================================
        #/////////////////////////////////////////////////////////////////////////////////

    
        
        # 对应实际采集的地震信号，取向上为正(和理论推导使用的方向相反)
        dataLst = []
        for ir in range(nrs):
            stream = Stream()
            dist = distarr[ir]

            # 计算延迟
            delayT = delayT0 
            if delayV0 > 0.0:
                delayT += np.sqrt(dist**2 + (deprcv-depsrc)**2)/delayV0

            # 计算走时
            travtP, travtS = self.compute_travt1d(dist)

            for i, comp in enumerate(['Z', 'R', 'T']):
                sgn = -1 if comp=='Z' else 1
                if i<2:
                    if calc_EXP:
                        stream.append(EXPgrn[ir][i].freq2time(delayT, travtP, travtS, sgn ))
                    if calc_VF:
                        stream.append(VFgrn [ir][i].freq2time(delayT, travtP, travtS, sgn ))
                    if calc_DC:
                        stream.append(DDgrn [ir][i].freq2time(delayT, travtP, travtS, sgn ))
                
                if calc_HF:
                    stream.append(HFgrn [ir][i].freq2time(delayT, travtP, travtS, sgn ))

                if calc_DC:
                    stream.append(DSgrn [ir][i].freq2time(delayT, travtP, travtS, sgn ))
                    stream.append(SSgrn [ir][i].freq2time(delayT, travtP, travtS, sgn ))

            dataLst.append(stream)


        return dataLst  

    
