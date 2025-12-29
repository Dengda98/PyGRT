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
import tempfile

from time import time
from copy import deepcopy

from ctypes import Array, pointer
from ctypes import _Pointer
from .c_interfaces import *
from .c_structures import *
from .pygrn import PyGreenFunction

__all__ = [
    "PyModel1D",
]


class PyModel1D:
    def __init__(self, modarr0:np.ndarray, depsrc:float, deprcv:float, allowLiquid:bool=False):
        '''
            Create 1D model instance, and insert the imaginary layer of source and receiver.

            :param    modarr0:    model array, in the format of [thickness(km), Vp(km/s), Vs(km/s), Rho(g/cm^3), Qp, Qs]  
            :param    depsrc:     source depth (km)  
            :param    deprcv:     receiver depth (km)  
            :param    allowLiquid:    whether liquid layers are allowed

        '''
        self.depsrc:float = depsrc 
        self.deprcv:float = deprcv 
        self.c_mod1d:c_GRT_MODEL1D 
        self.hasLiquid:bool = allowLiquid  # 传入的模型是否有液体层

        # 将modarr写入临时数组
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmpfile:
            np.savetxt(tmpfile, modarr0, "%.15e")
            tmp_path = tmpfile.name  # 获取临时文件路径

        try:
            c_mod1d_ptr = C_grt_read_mod1d_from_file(tmp_path.encode("utf-8"), depsrc, deprcv, allowLiquid)
            self.c_mod1d = c_mod1d_ptr.contents  # 这部分内存在C中申请，需由C函数释放。占用不多，这里跳过
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
        
        self.isrc = self.c_mod1d.isrc
        self.ircv = self.c_mod1d.ircv

        va = npct.as_array(self.c_mod1d.Va, (self.c_mod1d.n,))
        vb = npct.as_array(self.c_mod1d.Vb, (self.c_mod1d.n,))
        if np.any(vb == 0.0):
            self.hasLiquid = True
        
        self.vmin = min(np.min(va), np.min(vb))
        # 最小非零速度
        nonzero_vb = vb[vb > 0]
        self.vmin = min(np.min(va), np.min(nonzero_vb)) if nonzero_vb.size else np.min(va)
        self.vmax = max(np.max(va), np.max(vb))

    
    def compute_travt1d(self, dist:float):
        r"""
            Call the C function to calculate the travel time of the first P-wave and S-wave

            :param       dist:    epicentral distance (km)

            :return:
              - **travtP**  -  first P-wave arrival (s)
              - **travtS**  -  first S-wave arrival (s)
        """
        travtP = C_grt_compute_travt1d(
            self.c_mod1d.Thk,
            self.c_mod1d.Va,
            self.c_mod1d.n,
            self.c_mod1d.isrc,
            self.c_mod1d.ircv,
            dist
        )
        travtS = C_grt_compute_travt1d(
            self.c_mod1d.Thk,
            self.c_mod1d.Vb,
            self.c_mod1d.n,
            self.c_mod1d.isrc,
            self.c_mod1d.ircv,
            dist
        )

        return travtP, travtS


    def _init_grn(
        self,
        distarr:np.ndarray,
        nt:int, dt:float, upsampling_n:int, freqs:np.ndarray, wI:float, prefix:str=''):

        '''
            建立各个震源对应的格林函数类
        '''

        depsrc = self.depsrc
        deprcv = self.deprcv
        nr = len(distarr)

        pygrnLst:List[List[List[PyGreenFunction]]] = []
        c_grnArr = (((PCPLX*CHANNEL_NUM)*SRC_M_NUM)*nr)()
        
        for ir in range(len(distarr)):
            dist = distarr[ir]
            pygrnLst.append([])
            for isrc in range(SRC_M_NUM):
                pygrnLst[ir].append([])
                for ic, comp in enumerate(ZRTchs):

                    pygrn = PyGreenFunction(f'{prefix}{SRC_M_NAME_ABBR[isrc]}{comp}', nt, dt, upsampling_n, freqs, wI, dist, depsrc, deprcv)
                    pygrnLst[ir][isrc].append(pygrn)
                    c_grnArr[ir][isrc][ic] = pygrn.cmplx_grn.ctypes.data_as(PCPLX)

        return pygrnLst, c_grnArr
    

    def gen_gf_spectra(self, *args, **kwargs):
        raise NameError("Function 'gen_gf_spectra()' has been removed, use 'compute_grn' instead.")

    def compute_grn(
        self, 
        distarr:Union[np.ndarray,List[float],float], 
        nt:int, 
        dt:float, 
        upsampling_n:int = 1,
        freqband:Union[np.ndarray,List[float]]=[-1,-1],
        zeta:float=0.8, 
        keepAllFreq:bool=False,
        vmin_ref:float=0.0,
        keps:float=-1.0,  
        ampk:float=1.15,
        k0:float=5.0, 
        Length:float=0.0, 
        filonLength:float=0.0,
        safilonTol:float=0.0,
        filonCut:float=0.0,
        converg_method:Union[str,None]=None,
        delayT0:float=0.0,
        delayV0:float=0.0,
        calc_upar:bool=False,
        gf_source=['EX', 'VF', 'HF', 'DC'],
        statsfile:Union[str,None]=None, 
        statsidxs:Union[np.ndarray,List[int],None]=None, 
        print_runtime:bool=True):
        
        r'''
            Call the C function to calculate the Green's functions at multiple distances and return them in a list, 
            where each element is in the form of :class: 'obspy.Stream' type.

            :param    distarr:       array of epicentral distances (km), or a single float
            :param    nt:            number of time points. with the help of `SciPy`, nt no longer needs to be a power of 2
            :param    dt:            time interval (s)  
            :param    upsampling_n:  upsampling factor 
            :param    freqband:      frequency range (Hz)
            :param    zeta:          zeta is used to define the imaginary angular frequency, 
                                     :math:`\tilde{\omega} = \omega - j*w_I, w_I = \zeta*\pi/T, T=nt*dt` .
                                     see Bouchon (1981) and 张海明 (2021) for more details and tests.
            :param    keepAllFreq:   calculate all frequency points, no matter how low the frequency is
            :param    vmin_ref:      minimum reference velocity (km/s). the default vmin=max(minimum velocity, 0.1), used to define the upper limit of k integral
            :param    keps:          automatic convergence condition, see Yao and Harkrider (1983) for more details.
                                     negative value denotes not use.
            :param    ampk:          The factor that affect the upper limit of the k integral, see below.
            :param    k0:            k0 used to define the upper limit :math:`\tilde{k_{max}}=\sqrt{(k_{0}*\pi/hs)^2 + (ampk*w/vmin_{ref})^2}` , hs=max(abs(depsrc-deprcv),1.0)
            :param    Length:        integration step `dk=2\pi / (L*rmax)`, see Bouchon (1981) and 张海明 (2021) for the criterion, default set automatically.
            :param    filonLength:   integration step of Fixed-Interval Filon's Integration Method
            :param    safilonTol:    precision of Self-Adaptive Filon's Integration Method
            :param    filonCut:      The splitting point of DWM and (SA)FIM, k*=<filonCut>/rmax, default is 0
            :param    converg_method:   The method of explicit convergence, you can set "DCM", "PTAM" or "none". Default use "DCM" when abs(depsrc-deprcv) <= 1.0 km
            :param    calc_upar:     whether calculate the spatial derivatives of displacements.
            :param    gf_source:     The source type to be calculated
            :param    statsfile:     directory path for saving the statsfile during k integral, used to debug or observe the variations of :math:`F(k,\omega)` and :math:`F(k,\omega)J_m(kr)k`    
            :param    statsidxs:     only output the statsfile at specific frequency indexes. It is recommended to specify the indexes; 
                                     otherwise, by default, statsfiles of all frequency will be output, which probably occupy a lot of disk space
            :param    print_runtime: whether print runtime and some other infomation.

            :return:
                - **dataLst** -   Green's Functions at multiple distances, in a list of :class:`obspy.Stream`
                
        '''

        depsrc = self.depsrc
        deprcv = self.deprcv

        calc_EX:bool = 'EX' in gf_source
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
        if vmin_ref < 0:
            raise ValueError(f"vmin_ref ({vmin_ref}) < 0")
        
        if Length < 0.0:
            raise ValueError(f"Length ({Length}) < 0")
        if filonLength < 0.0:
            raise ValueError(f"filonLength ({filonLength}) < 0") 
        if filonCut < 0.0:
            raise ValueError(f"filonCut ({filonCut}) < 0") 
        if safilonTol < 0.0:
            raise ValueError(f"filonCut ({safilonTol}) < 0") 
        
        # 只能设置一种filon积分方法
        if safilonTol > 0.0 and filonLength > 0.0:
            raise ValueError(f"You should only set one of filonLength and safilonTol.")
        
        # 只能设置规定的收敛方法
        if converg_method is not None and converg_method not in ['DCM', 'PTAM', 'none']:
            raise ValueError(f'Wrong converg_method ({converg_method})')
        
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
        nf1 = min(int(np.ceil(f1/df)), nf-1)
        nf2 = min(int(np.floor(f2/df)), nf-1)
        if nf2 < nf1:
            nf2 = nf1

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

        # 若不指定显式收敛方法，则根据情况自动使用PTAM
        if converg_method is None and abs(depsrc - deprcv) <= 1.0:
            converg_method = 'DCM'

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

        # 初始化格林函数
        pygrnLst, c_grnArr = self._init_grn(distarr, nt, dt, upsampling_n, freqs, wI, '')
        
        pygrnLst_uiz = []
        c_grnArr_uiz = None
        pygrnLst_uir = []
        c_grnArr_uir = None
        if calc_upar:
            pygrnLst_uiz, c_grnArr_uiz = self._init_grn(distarr, nt, dt, upsampling_n, freqs, wI, 'z')
            pygrnLst_uir, c_grnArr_uir = self._init_grn(distarr, nt, dt, upsampling_n, freqs, wI, 'r')


        c_statsfile = None 
        if statsfile is not None:
            os.makedirs(statsfile, exist_ok=True)
            c_statsfile = c_char_p(statsfile.encode('utf-8'))

            nstatsidxs = 0 
            if statsidxs is None:
                statsidxs = np.arange(nf)

            statsidxs = np.array(statsidxs)
            # 不能有负数
            if np.any(statsidxs < 0):
                raise ValueError("negative value in statsidxs is not supported.")
            
            c_statsidxs = npct.as_ctypes(np.array(statsidxs).astype(np.uint64))   # size_t
            nstatsidxs = len(statsidxs)
        else:
            c_statsfile = c_statsidxs = None
            nstatsidxs = 0


        # ===========================================
        # 打印参数设置 
        if print_runtime:
            print(f"k0={k0}")
            print(f"ampk={ampk}")
            print(f"keps={keps}")
            print(f"vmin={vmin_ref}")
            print(f"Length={Length}")
            print(f"kcut={filonCut}")
            print(f"filonLength={filonLength}")
            print(f"safilonTol={safilonTol}")
            print(f'converg_method={converg_method}')
            
            print(f"nt={nt}")
            print(f"dt={dt}")
            print(f"winT={winT}")
            print(f"zeta={zeta}")
            print(f"delayT0={delayT0}")
            print(f"delayV0={delayV0}")
            print(f"tmax={tmax}")
            
            print(f"maxfreq(Hz)={freqs[nf-1]}")
            print(f"f1(Hz)={freqs[nf1]}")
            print(f"f2(Hz)={freqs[nf2]}")
            print(f"distances(km)=", distarr)
            if nstatsidxs > 0:
                print(f"statsfile_index=", statsidxs)



        KMET = c_K_INTEG_METHOD()

        hs = max(abs(depsrc - deprcv), 1.0)
        KMET.k0 = k0 * np.pi / hs
        KMET.ampk = ampk
        KMET.keps = keps if converg_method == 'none' else 0.0
        KMET.vmin = vmin_ref

        KMET.kcut = filonCut / rmax
        
        KMET.dk = 2.0*np.pi / (Length * rmax)
        
        KMET.applyFIM = filonLength > 0.0
        KMET.filondk = 2.0*np.pi / (filonLength * rmax) if filonLength > 0.0 else 0.0
        
        KMET.applySAFIM = safilonTol > 0.0
        KMET.sa_tol = safilonTol

        KMET.applyDCM = converg_method == 'DCM'
        KMET.applyPTAM = converg_method == 'PTAM'

        # 运行C库函数
        #/////////////////////////////////////////////////////////////////////////////////
        # 计算得到的格林函数的单位：
        #     单力源 HF[ZRT],VF[ZR]                  1e-15 cm/dyne
        #     爆炸源 EX[ZR]                          1e-20 cm/(dyne*cm)
        #     剪切源 DD[ZR],DS[ZRT],SS[ZRT]          1e-20 cm/(dyne*cm)
        #=================================================================================
        C_grt_integ_grn_spec(
            self.c_mod1d, nf1, nf2, c_freqs, nrs, c_rs, wI, keepAllFreq, pointer(KMET),
            print_runtime, calc_upar, 
            c_grnArr, c_grnArr_uiz, c_grnArr_uir,
            c_statsfile, nstatsidxs, c_statsidxs
        )
        #=================================================================================
        #/////////////////////////////////////////////////////////////////////////////////

        # 震源和场点层的物性，写入sac头段变量
        rcv_va = self.c_mod1d.Va[self.ircv]
        rcv_vb = self.c_mod1d.Vb[self.ircv]
        rcv_rho = self.c_mod1d.Rho[self.ircv]
        rcv_qainv = self.c_mod1d.Qainv[self.ircv]
        rcv_qbinv = self.c_mod1d.Qbinv[self.ircv]
        src_va = self.c_mod1d.Va[self.isrc]
        src_vb = self.c_mod1d.Vb[self.isrc]
        src_rho = self.c_mod1d.Rho[self.isrc]
        
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

            for im in range(SRC_M_NUM):
                if(not calc_EX and im==0):
                    continue
                if(not calc_VF and im==1):
                    continue
                if(not calc_HF and im==2):
                    continue
                if(not calc_DC and im>=3):
                    continue

                modr = SRC_M_ORDERS[im]
                sgn = 1
                for c in range(CHANNEL_NUM):
                    if(modr==0 and ZRTchs[c]=='T'):
                        continue
                    
                    sgn = -1 if ZRTchs[c]=='Z'=='Z' else 1
                    stream.append(pygrnLst[ir][im][c].freq2time(delayT, travtP, travtS, sgn ))
                    if(calc_upar):
                        stream.append(pygrnLst_uiz[ir][im][c].freq2time(delayT, travtP, travtS, sgn*(-1) ))
                        stream.append(pygrnLst_uir[ir][im][c].freq2time(delayT, travtP, travtS, sgn  ))


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
        keps:float=-1.0,  
        k0:float=5.0, 
        Length:float=15.0, 
        filonLength:float=0.0,
        safilonTol:float=0.0,
        filonCut:float=0.0,
        converg_method:Union[str,None]=None,
        calc_upar:bool=False,
        statsfile:Union[str,None]=None):

        r"""
            Call the C function to calculate the static Green's functions and return them in a dict

            :param       xarr:          coordinate array in the north direction (km), or a single float.
            :param       yarr:          coordinate array in the east direction (km), or a single float.
            :param       keps:          automatic convergence condition, see (Yao and Harkrider (1983) for more details.
                                        negative value denotes not use.
            :param       k0:            k0 used to define the upper limit :math:`\tilde{k_{max}}=(k_{0}*\pi/hs)^2`, hs=max(abs(depsrc-deprcv),1.0)
            :param       Length:        integration step `dk=2\pi / (L*rmax)`, default L=15
            :param       filonLength:   integration step of Fixed-Interval Filon's Integration Method
            :param       safilonTol:    precision of Self-Adaptive Filon's Integration Method
            :param       filonCut:      The splitting point of DWM and (SA)FIM, k*=<filonCut>/rmax, default is 0
            :param    converg_method:   The method of explicit convergence, you can set "DCM", "PTAM" or "none". Default use "DCM" when abs(depsrc-deprcv) <= 1.0 km
            :param       calc_upar:     whether calculate the spatial derivatives of displacements.
            :param       statsfile:     directory path for saving the statsfile during k integral, used to debug or observe the variations of :math:`F(k,\omega)` and :math:`F(k,\omega)J_m(kr)k` 
            
            :return:
                - **dataDct** -   static Green's function in a dict
        """

        if self.hasLiquid:
            raise NotImplementedError(
                "The feature for calculating static displacements "
                "in a model with liquid layers has not yet been implemented."
            )

        if Length < 0.0:
            raise ValueError(f"Length ({Length}) < 0")
        if filonLength < 0.0:
            raise ValueError(f"filonLength ({filonLength}) < 0") 
        if filonCut < 0.0:
            raise ValueError(f"filonCut ({filonCut}) < 0") 
        if safilonTol < 0.0:
            raise ValueError(f"filonCut ({safilonTol}) < 0") 
        
        # 只能设置一种filon积分方法
        if safilonTol > 0.0 and filonLength > 0.0:
            raise ValueError(f"You should only set one of filonLength and safilonTol.")
        
        # 只能设置规定的收敛方法
        if converg_method is not None and converg_method not in ['DCM', 'PTAM', 'none']:
            raise ValueError(f'Wrong converg_method ({converg_method})')

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
        
        # 设置波数积分间隔
        if Length == 0.0:
            Length = 15.0

        # 若不指定显式收敛方法，则根据情况自动使用PTAM
        if converg_method is None and abs(depsrc - deprcv) <= 1.0:
            converg_method = 'DCM'

        # 积分状态文件
        c_statsfile = None 
        if statsfile is not None:
            os.makedirs(statsfile, exist_ok=True)
            c_statsfile = c_char_p(statsfile.encode('utf-8'))

        # 初始化格林函数
        pygrn = np.zeros((nr, SRC_M_NUM, CHANNEL_NUM), dtype=NPCT_REAL_TYPE, order='C');       c_pygrn = npct.as_ctypes(pygrn)
        pygrn_uiz = np.zeros((nr, SRC_M_NUM, CHANNEL_NUM), dtype=NPCT_REAL_TYPE, order='C');   c_pygrn_uiz = npct.as_ctypes(pygrn_uiz)
        pygrn_uir = np.zeros((nr, SRC_M_NUM, CHANNEL_NUM), dtype=NPCT_REAL_TYPE, order='C');   c_pygrn_uir = npct.as_ctypes(pygrn_uir)

        if not calc_upar:
            c_pygrn_uiz = c_pygrn_uir = None
        

        KMET = c_K_INTEG_METHOD()

        hs = max(abs(depsrc - deprcv), 1.0)
        KMET.k0 = k0 * np.pi / hs
        KMET.keps = keps if converg_method == 'none' else 0.0

        # 最大震中距
        rmax = np.max(rs)

        KMET.kcut = filonCut / rmax
        
        KMET.dk = 2.0*np.pi / (Length * rmax)
        
        KMET.applyFIM = filonLength > 0.0
        KMET.filondk = 2.0*np.pi / (filonLength * rmax) if filonLength > 0.0 else 0.0
        
        KMET.applySAFIM = safilonTol > 0.0
        KMET.sa_tol = safilonTol

        KMET.applyDCM = converg_method == 'DCM'
        KMET.applyPTAM = converg_method == 'PTAM'

        # 运行C库函数
        #/////////////////////////////////////////////////////////////////////////////////
        # 计算得到的格林函数的单位：
        #     单力源 HF[ZRT],VF[ZR]                  1e-15 cm/dyne
        #     爆炸源 EX[ZR]                          1e-20 cm/(dyne*cm)
        #     剪切源 DD[ZR],DS[ZRT],SS[ZRT]          1e-20 cm/(dyne*cm)
        #=================================================================================
        C_grt_integ_static_grn(
            self.c_mod1d, nr, c_rs, pointer(KMET),
            calc_upar, c_pygrn, c_pygrn_uiz, c_pygrn_uir,
            c_statsfile
        )
        #=================================================================================
        #/////////////////////////////////////////////////////////////////////////////////

        # 震源和场点层的物性
        rcv_va = self.c_mod1d.Va[self.ircv]
        rcv_vb = self.c_mod1d.Vb[self.ircv]
        rcv_rho = self.c_mod1d.Rho[self.ircv]
        src_va = self.c_mod1d.Va[self.isrc]
        src_vb = self.c_mod1d.Vb[self.isrc]
        src_rho = self.c_mod1d.Rho[self.isrc]

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
        for isrc in range(SRC_M_NUM):
            src_name = SRC_M_NAME_ABBR[isrc]
            for ic, comp in enumerate(ZRTchs):
                sgn = -1 if comp=='Z' else 1
                dataDct[f'{src_name}{comp}'] = sgn * pygrn[:,isrc,ic].reshape((nx, ny), order='F')
                if calc_upar:
                    dataDct[f'z{src_name}{comp}'] = sgn * pygrn_uiz[:,isrc,ic].reshape((nx, ny), order='F') * (-1)
                    dataDct[f'r{src_name}{comp}'] = sgn * pygrn_uir[:,isrc,ic].reshape((nx, ny), order='F')

        return dataDct