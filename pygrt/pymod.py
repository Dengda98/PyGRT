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
from typing import List, Dict, Union

from time import time
from copy import deepcopy

from ctypes import Array, pointer
from ctypes import _Pointer
from .c_interfaces import *
from .c_structures import *
from .pygrn import PyGreenFunction

PC_GRN2D = Array[Array[c_PGRN]]


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

        self.isrc = isrc
        self.ircv = ircv

        self.vmax = np.max(self.modarr[:, 1:3])
        self.vmin = np.min(self.modarr[:, 1:3])
        if self.vmin <= 0.0:
            raise ValueError("Zero Velocity ??")
        
    
    def compute_travt1d(self, dist:float):
        r"""
            调用C程序，计算初至P波和S波的走时

            :param       dist:    震中距

            :return:
              - **travtP**  -  初至P波走时(s)
              - **travtS**  -  初至S波走时(s)
        """
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


    def _init_grn(
        self,
        distarr:np.ndarray,
        calc_EXP:bool, calc_VF:bool, calc_HF:bool, calc_DC:bool,
        C_EXPgrn:PC_GRN2D, C_VFgrn:PC_GRN2D, C_HFgrn:PC_GRN2D, 
        C_DDgrn:PC_GRN2D, C_DSgrn:PC_GRN2D, C_SSgrn:PC_GRN2D, 
        nt:int, dt:float, freqs:np.ndarray, wI:float, prefix:str=''):

        '''
            建立各个震源对应的格林函数类
        '''

        depsrc = self.depsrc
        deprcv = self.deprcv

        EXPgrn:List[List[PyGreenFunction]] = []
        VFgrn:List[List[PyGreenFunction]] = []
        DDgrn:List[List[PyGreenFunction]] = []
        HFgrn:List[List[PyGreenFunction]] = []
        DSgrn:List[List[PyGreenFunction]] = []
        SSgrn:List[List[PyGreenFunction]] = []
        
        for ir in range(len(distarr)):
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
                        grn = PyGreenFunction(f'{prefix}EX{comp}', nt, dt, freqs, wI, dist, depsrc, deprcv)
                        EXPgrn[ir].append(grn)
                        C_EXPgrn[ir][i] = pointer(grn.c_grn)

                    if calc_VF:
                        grn = PyGreenFunction(f'{prefix}VF{comp}', nt, dt, freqs, wI, dist, depsrc, deprcv)
                        VFgrn[ir].append(grn)
                        C_VFgrn[ir][i] = pointer(grn.c_grn)
                    
                    if calc_DC:
                        grn = PyGreenFunction(f'{prefix}DD{comp}', nt, dt, freqs, wI, dist, depsrc, deprcv)
                        DDgrn[ir].append(grn)
                        C_DDgrn[ir][i] = pointer(grn.c_grn)
                
                if calc_HF:
                    grn = PyGreenFunction(f'{prefix}HF{comp}', nt, dt, freqs, wI, dist, depsrc, deprcv)
                    HFgrn[ir].append(grn)
                    C_HFgrn[ir][i] = pointer(grn.c_grn)

                if calc_DC:
                    grn = PyGreenFunction(f'{prefix}DS{comp}', nt, dt, freqs, wI, dist, depsrc, deprcv)
                    DSgrn[ir].append(grn)
                    C_DSgrn[ir][i] = pointer(grn.c_grn)
                    grn = PyGreenFunction(f'{prefix}SS{comp}', nt, dt, freqs, wI, dist, depsrc, deprcv)
                    SSgrn[ir].append(grn)
                    C_SSgrn[ir][i] = pointer(grn.c_grn) 

        return EXPgrn, VFgrn, HFgrn, DDgrn, DSgrn, SSgrn

    def gen_gf_spectra(self, *args, **kwargs):
        r"Bad function name, has already been removed. Use 'compute_grn' instead."
        raise NameError("Function 'gen_gf_spectra()' has been removed, use 'compute_grn' instead.")

    def compute_grn(
        self, 
        distarr:Union[np.ndarray,List[float],float], 
        nt:int, 
        dt:float, 
        freqband:Union[np.ndarray,List[float]]=[-1,-1],
        zeta:float=0.8, 
        vmin_ref:float=0.0,
        keps:float=-1.0,  
        ampk:float=1.15,
        k0:float=5.0, 
        Length:float=0.0, 
        filonLC:Union[np.ndarray,List[float]]=[0.0,0.0],
        delayT0:float=0.0,
        delayV0:float=0.0,
        calc_upar:bool=False,
        gf_source=['EXP', 'VF', 'HF', 'DC'],
        statsfile:Union[str,None]=None, 
        statsidxs:Union[np.ndarray,List[int],None]=None, 
        print_runtime:bool=True):
        
        r'''
            
            调用C库计算格林函数的主函数，以列表的形式返回，其中每个元素为对应震中距的格林函数 :class:`obspy.Stream` 类型。
            

            :param    distarr:       多个震中距(km) 的数组, 或单个震中距的浮点数
            :param    nt:            时间点数，借助于 `SciPy`，nt不再要求是2的幂次
            :param    dt:            采样间隔(s)  
            :param    freqband:      频率范围(Hz)，以此确定待计算的离散频率点
            :param    zeta:          定义虚频率的系数 :math:`\zeta` ， 虚频率 :math:`\tilde{\omega} = \omega - j*w_I, w_I = \zeta*\pi/T, T=nt*dt` , T为时窗长度。
                                     使用离散波数积分时为了避开附加源以及奇点的影响， :ref:`(Bouchon, 1981) <bouchon_1981>`  在频率上添加微小虚部，
                                     更多测试见 :ref:`(张海明, 2021) <zhang_book_2021>`
            :param    vmin_ref:      最小参考速度，默认vmin=max(minimum velocity, 0.1)，用于定义波数积分上限，小于0则在达到积分上限后使用峰谷平均法
                                    （默认当震源和场点深度差<=1km时自动使用峰谷平均法）
            :param    keps:          波数k积分收敛条件，见 :ref:`(Yao and Harkrider, 1983) <yao&harkrider_1983>`  :ref:`(初稿) <yao_init_manuscripts>`，
                                     为负数代表不提前判断收敛，按照波数积分上限进行积分
            :param    ampk:          影响波数k积分上限的系数，见下方
            :param    k0:            波数k积分的上限 :math:`\tilde{k_{max}}=\sqrt{(k_{0}*\pi/hs)^2 + (ampk*w/vmin_{ref})^2}` , 波数k积分循环必须退出, hs=max(震源和台站深度差,1.0)
            :param    Length:        定义波数k积分的间隔 `dk=2\pi / (L*rmax)`, 选取要求见 :ref:`(Bouchon, 1981) <bouchon_1981>` 
                                     :ref:`(张海明, 2021) <zhang_book_2021>`，默认自动选择
            :param    filonLC:       Filon积分的间隔 filonLength, 和波数积分和Filon积分的分割点filonCut, k*=<filonCut>/rmax
            :param    calc_upar:     是否计算位移u的空间导数
            :param    gf_source:     待计算的震源类型
            :param    statsfile:     波数k积分（包括Filon积分和峰谷平均法）的过程记录文件，常用于debug或者观察积分过程中 :math:`F(k,\omega)` 和  :math:`F(k,\omega)J_m(kr)k` 的变化    
            :param    statsidxs:     仅输出特定频点的过程记录文件，建议给定频点，否则默认所有频率点的记录文件都输出，很占空间
            :param    print_runtime: 是否打印运行时间

            :return:
                - **dataLst** -   列表，每个元素为 :class:`obspy.Stream` 类型 )
                
        '''

        depsrc = self.depsrc
        deprcv = self.deprcv

        calc_EXP:bool = 'EXP' in gf_source
        calc_VF:bool = 'VF' in gf_source
        calc_HF:bool = 'HF' in gf_source
        calc_DC:bool = 'DC' in gf_source

        if isinstance(distarr, float) or isinstance(distarr, int):
            distarr = np.array([distarr*1.0]) 

        distarr = np.array(distarr)
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
        
        if Length < 0.0:
            raise ValueError(f"Length ({Length}) < 0")
        if np.any(filonLC) < 0.0:
            raise ValueError(f"filonLC ({filonLC}) < 0") 
        
        filonLC = np.array(filonLC).astype(NPCT_REAL_TYPE)

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
            if abs(depsrc - deprcv) <= 1.0:
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

        # 位移u的空间导数
        C_EXPgrn_uiz = C_VFgrn_uiz = C_HFgrn_uiz = C_DDgrn_uiz = C_DSgrn_uiz = C_SSgrn_uiz = None
        C_EXPgrn_uir = C_VFgrn_uir = C_HFgrn_uir = C_DDgrn_uir = C_DSgrn_uir = C_SSgrn_uir = None
        if calc_upar:
            C_EXPgrn_uiz = ((c_PGRN*2)*nrs)() if calc_EXP else None
            C_VFgrn_uiz = ((c_PGRN*2)*nrs)() if calc_VF else None
            C_HFgrn_uiz = ((c_PGRN*3)*nrs)() if calc_HF else None
            C_DDgrn_uiz = ((c_PGRN*2)*nrs)() if calc_DC else None
            C_DSgrn_uiz = ((c_PGRN*3)*nrs)() if calc_DC else None
            C_SSgrn_uiz = ((c_PGRN*3)*nrs)() if calc_DC else None
            #
            C_EXPgrn_uir = ((c_PGRN*2)*nrs)() if calc_EXP else None
            C_VFgrn_uir = ((c_PGRN*2)*nrs)() if calc_VF else None
            C_HFgrn_uir = ((c_PGRN*3)*nrs)() if calc_HF else None
            C_DDgrn_uir = ((c_PGRN*2)*nrs)() if calc_DC else None
            C_DSgrn_uir = ((c_PGRN*3)*nrs)() if calc_DC else None
            C_SSgrn_uir = ((c_PGRN*3)*nrs)() if calc_DC else None


        EXPgrn, VFgrn, HFgrn, DDgrn, DSgrn, SSgrn = self._init_grn(
            distarr, calc_EXP, calc_VF, calc_HF, calc_DC, 
            C_EXPgrn, C_VFgrn, C_HFgrn, C_DDgrn, C_DSgrn, C_SSgrn, 
            nt, dt, freqs, wI)
        
        EXPgrn_uiz, VFgrn_uiz, HFgrn_uiz, DDgrn_uiz, DSgrn_uiz, SSgrn_uiz = ([] for _ in range(6))
        EXPgrn_uir, VFgrn_uir, HFgrn_uir, DDgrn_uir, DSgrn_uir, SSgrn_uir = ([] for _ in range(6))
        if calc_upar:
            EXPgrn_uiz, VFgrn_uiz, HFgrn_uiz, DDgrn_uiz, DSgrn_uiz, SSgrn_uiz = self._init_grn(
            distarr, calc_EXP, calc_VF, calc_HF, calc_DC, 
            C_EXPgrn_uiz, C_VFgrn_uiz, C_HFgrn_uiz, C_DDgrn_uiz, C_DSgrn_uiz, C_SSgrn_uiz, 
            nt, dt, freqs, wI, 'z')
            EXPgrn_uir, VFgrn_uir, HFgrn_uir, DDgrn_uir, DSgrn_uir, SSgrn_uir = self._init_grn(
            distarr, calc_EXP, calc_VF, calc_HF, calc_DC, 
            C_EXPgrn_uir, C_VFgrn_uir, C_HFgrn_uir, C_DDgrn_uir, C_DSgrn_uir, C_SSgrn_uir, 
            nt, dt, freqs, wI, 'r')


        c_statsfile = None 
        if statsfile is not None:
            os.makedirs(statsfile, exist_ok=True)
            c_statsfile = c_char_p(statsfile.encode('utf-8'))

        nstatsidxs = 0 
        if statsidxs is None:
            statsidxs = np.array([-1])

        statsidxs = np.array(statsidxs)
        c_statsidxs = npct.as_ctypes(np.array(statsidxs).astype('i'))
        nstatsidxs = len(statsidxs)


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
            if filonLC[0] > 0.0:
                print(f",{filonLC}, using FIM.")
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
            vmin_ref, keps, ampk, k0, Length, filonLC[0], filonLC[1], print_runtime,
            C_EXPgrn, C_VFgrn, C_HFgrn, C_DDgrn, C_DSgrn, C_SSgrn, 
            calc_upar, 
            C_EXPgrn_uiz, C_VFgrn_uiz, C_HFgrn_uiz, C_DDgrn_uiz, C_DSgrn_uiz, C_SSgrn_uiz, 
            C_EXPgrn_uir, C_VFgrn_uir, C_HFgrn_uir, C_DDgrn_uir, C_DSgrn_uir, C_SSgrn_uir, 
            c_statsfile, nstatsidxs, c_statsidxs
        )
        #=================================================================================
        #/////////////////////////////////////////////////////////////////////////////////

        # 震源和场点层的物性，写入sac头段变量
        rcv_va = self.modarr[self.ircv, 1]
        rcv_vb = self.modarr[self.ircv, 2]
        rcv_rho = self.modarr[self.ircv, 3]
        rcv_qainv = 1.0/self.modarr[self.ircv, 4]
        rcv_qbinv = 1.0/self.modarr[self.ircv, 5]
        src_va = self.modarr[self.isrc, 1]
        src_vb = self.modarr[self.isrc, 2]
        src_rho = self.modarr[self.isrc, 3]
        
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

                if calc_upar:
                    if i<2:
                        if calc_EXP:
                            stream.append(EXPgrn_uiz[ir][i].freq2time(delayT, travtP, travtS, sgn*(-1) ))
                            stream.append(EXPgrn_uir[ir][i].freq2time(delayT, travtP, travtS, sgn ))
                        if calc_VF:
                            stream.append(VFgrn_uiz [ir][i].freq2time(delayT, travtP, travtS, sgn*(-1) ))
                            stream.append(VFgrn_uir [ir][i].freq2time(delayT, travtP, travtS, sgn ))
                        if calc_DC:
                            stream.append(DDgrn_uiz [ir][i].freq2time(delayT, travtP, travtS, sgn*(-1) ))
                            stream.append(DDgrn_uir [ir][i].freq2time(delayT, travtP, travtS, sgn ))
                    
                    if calc_HF:
                        stream.append(HFgrn_uiz [ir][i].freq2time(delayT, travtP, travtS, sgn*(-1) ))
                        stream.append(HFgrn_uir [ir][i].freq2time(delayT, travtP, travtS, sgn ))

                    if calc_DC:
                        stream.append(DSgrn_uiz [ir][i].freq2time(delayT, travtP, travtS, sgn*(-1) ))
                        stream.append(DSgrn_uir [ir][i].freq2time(delayT, travtP, travtS, sgn ))
                        stream.append(SSgrn_uiz [ir][i].freq2time(delayT, travtP, travtS, sgn*(-1) ))
                        stream.append(SSgrn_uir [ir][i].freq2time(delayT, travtP, travtS, sgn ))

            # 在sac头段变量部分
            for tr in stream:
                SAC = tr.stats.sac
                SAC['user1'] = rcv_va
                SAC['user2'] = rcv_vb
                SAC['user3'] = rcv_rho
                SAC['user4'] = rcv_qainv
                SAC['user5'] = rcv_qbinv
                SAC['user6'] = src_va
                SAC['user7'] = src_vb
                SAC['user8'] = src_rho

            dataLst.append(stream)


        return dataLst  

    

    def compute_static_grn(
        self,
        xarr:Union[np.ndarray,List[float],float], 
        yarr:Union[np.ndarray,List[float],float], 
        vmin_ref:float=0.0,
        keps:float=-1.0,  
        k0:float=5.0, 
        Length:float=15.0, 
        filonLC:Union[np.ndarray,List[float]]=[0.0,0.0],
        calc_upar:bool=False,
        statsfile:Union[str,None]=None):

        r"""
            调用C库计算静态格林函数，以字典的形式返回

            :param       xarr:          北向坐标数组，或单个浮点数
            :param       yarr:          东向坐标数组，或单个浮点数
            :param       vmin_ref:      最小参考速度（具体数值不使用），小于0则在达到积分上限后使用峰谷平均法
                                       （默认当震源和场点深度差<=0.5km时自动使用峰谷平均法）
            :param       keps:          波数k积分收敛条件，见 :ref:`(Yao and Harkrider, 1983) <yao&harkrider_1983>`  :ref:`(初稿) <yao_init_manuscripts>`，
                                        为负数代表不提前判断收敛，按照波数积分上限进行积分
            :param       k0:            波数k积分的上限 :math:`\tilde{k_{max}}=(k_{0}*\pi/hs)^2` , 波数k积分循环必须退出, hs=max(震源和台站深度差,1.0)
            :param       Length:        定义波数k积分的间隔 `dk=2\pi / (L*rmax)`, 默认15；负数表示使用Filon积分
            :param       filonLC:       Filon积分的间隔 filonLength, 和波数积分和Filon积分的分割点filonCut, k*=<filonCut>/rmax
            :param       calc_upar:     是否计算位移u的空间导数
            :param       statsfile:     波数k积分（包括Filon积分和峰谷平均法）的过程记录文件，常用于debug或者观察积分过程中 :math:`F(k,\omega)` 和  :math:`F(k,\omega)J_m(kr)k` 的变化    

            :return:
                - **dataDct** -   字典形式的格林函数
        """

        if Length < 0.0:
            raise ValueError(f"Length ({Length}) < 0")
        if np.any(filonLC) < 0.0:
            raise ValueError(f"filonLC ({filonLC}) < 0") 
        

        depsrc = self.depsrc
        deprcv = self.deprcv

        if isinstance(xarr, float) or isinstance(xarr, int):
            xarr = np.array([xarr*1.0]) 
        xarr = np.array(xarr)

        if isinstance(yarr, float) or isinstance(yarr, int):
            yarr = np.array([yarr*1.0]) 
        yarr = np.array(yarr)

        nx = len(xarr)
        ny = len(yarr)
        nr = nx*ny
        rs = np.zeros((nr,), dtype=NPCT_REAL_TYPE)
        for iy in range(ny):
            for ix in range(nx):
                rs[ix + iy*nx] = max(np.sqrt(xarr[ix]**2 + yarr[iy]**2), 1e-5)
        c_rs = npct.as_ctypes(rs)

        # 参考最小速度
        if vmin_ref == 0.0:
            vmin_ref = max(self.vmin, 0.1)
            if abs(depsrc - deprcv) <= 1.0:
                vmin_ref = - abs(vmin_ref)  # 自动使用PTAM
        
        # 设置波数积分间隔
        if Length == 0.0:
            Length = 15.0

        # 积分状态文件
        c_statsfile = None 
        if statsfile is not None:
            os.makedirs(statsfile, exist_ok=True)
            c_statsfile = c_char_p(statsfile.encode('utf-8'))

        # 初始化格林函数
        EXPgrn = np.zeros((nr,2), dtype=NPCT_REAL_TYPE); C_EXPgrn = npct.as_ctypes(EXPgrn.reshape(-1))
        VFgrn = np.zeros((nr,2), dtype=NPCT_REAL_TYPE); C_VFgrn = npct.as_ctypes(VFgrn.reshape(-1))
        HFgrn = np.zeros((nr,3), dtype=NPCT_REAL_TYPE); C_HFgrn = npct.as_ctypes(HFgrn.reshape(-1))
        DDgrn = np.zeros((nr,2), dtype=NPCT_REAL_TYPE); C_DDgrn = npct.as_ctypes(DDgrn.reshape(-1))
        DSgrn = np.zeros((nr,3), dtype=NPCT_REAL_TYPE); C_DSgrn = npct.as_ctypes(DSgrn.reshape(-1))
        SSgrn = np.zeros((nr,3), dtype=NPCT_REAL_TYPE); C_SSgrn = npct.as_ctypes(SSgrn.reshape(-1))

        # 位移u的空间导数
        EXPgrn_uiz = np.zeros((nr,2), dtype=NPCT_REAL_TYPE); C_EXPgrn_uiz = npct.as_ctypes(EXPgrn_uiz.reshape(-1))
        VFgrn_uiz = np.zeros((nr,2), dtype=NPCT_REAL_TYPE); C_VFgrn_uiz = npct.as_ctypes(VFgrn_uiz.reshape(-1))
        HFgrn_uiz = np.zeros((nr,3), dtype=NPCT_REAL_TYPE); C_HFgrn_uiz = npct.as_ctypes(HFgrn_uiz.reshape(-1))
        DDgrn_uiz = np.zeros((nr,2), dtype=NPCT_REAL_TYPE); C_DDgrn_uiz = npct.as_ctypes(DDgrn_uiz.reshape(-1))
        DSgrn_uiz = np.zeros((nr,3), dtype=NPCT_REAL_TYPE); C_DSgrn_uiz = npct.as_ctypes(DSgrn_uiz.reshape(-1))
        SSgrn_uiz = np.zeros((nr,3), dtype=NPCT_REAL_TYPE); C_SSgrn_uiz = npct.as_ctypes(SSgrn_uiz.reshape(-1))

        EXPgrn_uir = np.zeros((nr,2), dtype=NPCT_REAL_TYPE); C_EXPgrn_uir = npct.as_ctypes(EXPgrn_uir.reshape(-1))
        VFgrn_uir = np.zeros((nr,2), dtype=NPCT_REAL_TYPE); C_VFgrn_uir = npct.as_ctypes(VFgrn_uir.reshape(-1))
        HFgrn_uir = np.zeros((nr,3), dtype=NPCT_REAL_TYPE); C_HFgrn_uir = npct.as_ctypes(HFgrn_uir.reshape(-1))
        DDgrn_uir = np.zeros((nr,2), dtype=NPCT_REAL_TYPE); C_DDgrn_uir = npct.as_ctypes(DDgrn_uir.reshape(-1))
        DSgrn_uir = np.zeros((nr,3), dtype=NPCT_REAL_TYPE); C_DSgrn_uir = npct.as_ctypes(DSgrn_uir.reshape(-1))
        SSgrn_uir = np.zeros((nr,3), dtype=NPCT_REAL_TYPE); C_SSgrn_uir = npct.as_ctypes(SSgrn_uir.reshape(-1))
        
        if not calc_upar:
            C_EXPgrn_uiz = C_VFgrn_uiz = C_HFgrn_uiz = C_DDgrn_uiz = C_DSgrn_uiz = C_SSgrn_uiz = None
            C_EXPgrn_uir = C_VFgrn_uir = C_HFgrn_uir = C_DDgrn_uir = C_DSgrn_uir = C_SSgrn_uir = None


        # 运行C库函数
        #/////////////////////////////////////////////////////////////////////////////////
        # 计算得到的格林函数的单位：
        #     单力源 HF[ZRT],VF[ZR]                  1e-15 cm/dyne
        #     爆炸源 EX[ZR]                          1e-20 cm/(dyne*cm)
        #     剪切源 DD[ZR],DS[ZRT],SS[ZRT]          1e-20 cm/(dyne*cm)
        #=================================================================================
        C_integ_static_grn(
            self.c_pymod1d, nr, c_rs, vmin_ref, keps, k0, Length, filonLC[0], filonLC[1],
            C_EXPgrn, C_VFgrn, C_HFgrn, C_DDgrn, C_DSgrn, C_SSgrn, 
            calc_upar, 
            C_EXPgrn_uiz, C_VFgrn_uiz, C_HFgrn_uiz, C_DDgrn_uiz, C_DSgrn_uiz, C_SSgrn_uiz, 
            C_EXPgrn_uir, C_VFgrn_uir, C_HFgrn_uir, C_DDgrn_uir, C_DSgrn_uir, C_SSgrn_uir, 
            c_statsfile
        )
        #=================================================================================
        #/////////////////////////////////////////////////////////////////////////////////

        # 震源和场点层的物性
        rcv_va = self.modarr[self.ircv, 1]
        rcv_vb = self.modarr[self.ircv, 2]
        rcv_rho = self.modarr[self.ircv, 3]
        src_va = self.modarr[self.isrc, 1]
        src_vb = self.modarr[self.isrc, 2]
        src_rho = self.modarr[self.isrc, 3]

        # 结果字典
        dataDct = {}
        dataDct['_xarr'] = xarr.copy()
        dataDct['_yarr'] = yarr.copy()
        dataDct['_src_va'] = src_va
        dataDct['_src_vb'] = src_vb
        dataDct['_src_rho'] = src_rho
        dataDct['_rcv_va'] = rcv_va
        dataDct['_rcv_vb'] = rcv_vb
        dataDct['_rcv_rho'] = rcv_rho

        # 整理结果，将每个格林函数以2d矩阵的形式存储，shape=(nx, ny)
        for i, ch in enumerate(['Z', 'R', 'T']):
            sgn = -1 if ch=='Z' else 1
            if i<2:
                dataDct[f'EX{ch}'] = sgn * EXPgrn[:,i].reshape((nx, ny), order='F')
                dataDct[f'VF{ch}'] = sgn * VFgrn[:,i].reshape((nx, ny), order='F')
                dataDct[f'DD{ch}'] = sgn * DDgrn[:,i].reshape((nx, ny), order='F')
            
            dataDct[f'HF{ch}'] = sgn * HFgrn[:,i].reshape((nx, ny), order='F')
            dataDct[f'DS{ch}'] = sgn * DSgrn[:,i].reshape((nx, ny), order='F')
            dataDct[f'SS{ch}'] = sgn * SSgrn[:,i].reshape((nx, ny), order='F')

            if calc_upar:
                if i<2:
                    dataDct[f'zEX{ch}'] = sgn * EXPgrn_uiz[:,i].reshape((nx, ny), order='F') * (-1)
                    dataDct[f'rEX{ch}'] = sgn * EXPgrn_uir[:,i].reshape((nx, ny), order='F')
                    dataDct[f'zVF{ch}'] = sgn * VFgrn_uiz[:,i].reshape((nx, ny), order='F') * (-1)
                    dataDct[f'rVF{ch}'] = sgn * VFgrn_uir[:,i].reshape((nx, ny), order='F')
                    dataDct[f'zDD{ch}'] = sgn * DDgrn_uiz[:,i].reshape((nx, ny), order='F') * (-1)
                    dataDct[f'rDD{ch}'] = sgn * DDgrn_uir[:,i].reshape((nx, ny), order='F')
                
                dataDct[f'zHF{ch}'] = sgn * HFgrn_uiz[:,i].reshape((nx, ny), order='F') * (-1)
                dataDct[f'rHF{ch}'] = sgn * HFgrn_uir[:,i].reshape((nx, ny), order='F')
                dataDct[f'zDS{ch}'] = sgn * DSgrn_uiz[:,i].reshape((nx, ny), order='F') * (-1)
                dataDct[f'rDS{ch}'] = sgn * DSgrn_uir[:,i].reshape((nx, ny), order='F')
                dataDct[f'zSS{ch}'] = sgn * SSgrn_uiz[:,i].reshape((nx, ny), order='F') * (-1)
                dataDct[f'rSS{ch}'] = sgn * SSgrn_uir[:,i].reshape((nx, ny), order='F')

        return dataDct