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
from typing import List, Dict, Union, Literal
import tempfile

from time import time
from copy import deepcopy

from ctypes import Array, pointer
from ctypes import _Pointer
from ctypes import cast
from .c_interfaces import *
from .c_structures import *
from .modal import PyDispersion, _normalize_wave_type
from .pygrn import PyGreenFunction

__all__ = [
    "PyModel1D",
]


class PyModel1D:
    def __init__(self, modarr0:np.ndarray, depsrc:float, deprcv:float, allowLiquid:bool=True,
                 topbound:Literal['free', 'rigid', 'halfspace']='free', 
                 botbound:Literal['free', 'rigid', 'halfspace']='halfspace'):
        '''
            Create 1D model instance, and insert the imaginary layer of source and receiver.

            :param    modarr0:    model array, in the format of [thickness(km), Vp(km/s), Vs(km/s), Rho(g/cm^3), Qp, Qs]  
            :param    depsrc:     source depth (km)  
            :param    deprcv:     receiver depth (km)  
            :param    allowLiquid:    (deprecated) unused argument
            :param    topbound:       boundary condition of the top layer
            :param    botbound:       boundary condition of the bottom layer

        '''
        self.depsrc:float = depsrc 
        self.deprcv:float = deprcv 
        self.c_mod1d:c_MODEL1D 
        self.topbound:str = topbound
        self.botbound:str = botbound
        self.hasLiquid = False
        self._modarr0 = np.array(modarr0, dtype=NPCT_REAL_TYPE, copy=True)

        if depsrc < 0.0 or deprcv < 0.0:
            raise ValueError("Negative source depth or receiver depth is not supported.")

        boundDct = {
            'free': 0,
            'rigid': 1,
            'halfspace': 2
        }

        if topbound not in boundDct:
            raise ValueError(f"Unsupported topbound={topbound}.")
        if botbound not in boundDct:
            raise ValueError(f"Unsupported botbound={botbound}.")

        # 将modarr写入临时数组
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmpfile:
            ncol = modarr0.shape[-1]
            np.savetxt(tmpfile, modarr0.reshape(-1, ncol), "%.15e")
            tmp_path = tmpfile.name  # 获取临时文件路径

        try:
            c_mod1d_ptr = C_grt_read_mod1d_from_file(tmp_path.encode("utf-8"), depsrc, deprcv, True)
            self.c_mod1d = c_mod1d_ptr.contents  # 这部分内存在C中申请，需由C函数释放。占用不多，这里跳过
            C_grt_set_mod1d_boundary(self.c_mod1d, boundDct[topbound], boundDct[botbound])
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
        
        # 设置边界条件
        self.c_mod1d.topbound = boundDct[topbound]
        self.c_mod1d.botbound = boundDct[botbound]
        
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


    def _read_modal_model(self):
        '''
            读取未插入震源和台站虚拟层的模型，供面波相关 C 函数使用。
        '''
        boundDct = {
            'free': 0,
            'rigid': 1,
            'halfspace': 2
        }

        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmpfile:
            ncol = self._modarr0.shape[-1]
            np.savetxt(tmpfile, self._modarr0.reshape(-1, ncol), "%.15e")
            tmp_path = tmpfile.name

        try:
            c_mod1d_ptr = C_grt_read_mod1d_from_file(tmp_path.encode("utf-8"), -1.0, -1.0, True)
            C_grt_set_mod1d_boundary(c_mod1d_ptr, boundDct[self.topbound], boundDct[self.botbound])
            return c_mod1d_ptr
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)


    @staticmethod
    def _resolve_dispersion_freqs(freqs, freq_range):
        if freqs is not None and freq_range is not None:
            raise ValueError("Set only one of freqs and freq_range.")
        if freqs is None and freq_range is None:
            raise ValueError("Set one of freqs and freq_range.")

        if freq_range is not None:
            if len(freq_range) != 3:
                raise ValueError("freq_range must be (f1, f2, df).")
            f1, f2, df = [float(x) for x in freq_range]
            if df <= 0.0:
                raise ValueError("freq_range df must be positive.")
            if f1 < 0.0 or f2 <= 0.0:
                raise ValueError("freq_range f1/f2 must be nonnegative/positive.")
            if f1 > f2:
                raise ValueError("freq_range f1 must be <= f2.")
            f1 = max(f1, df)
            f2 = max(f1, f2)
            nfreq = int(np.floor((f2 - f1)/df + 1e-8)) + 1
            freqs = f1 + np.arange(nfreq, dtype=NPCT_REAL_TYPE)*df
        else:
            freqs = np.asarray(freqs, dtype=NPCT_REAL_TYPE)
            if freqs.ndim != 1:
                raise ValueError("freqs must be a 1-D array.")
            freqs = np.sort(freqs.copy())

        if freqs.size == 0:
            raise ValueError("freqs must not be empty.")
        if np.any(freqs <= 0.0):
            raise ValueError("freqs must be positive.")

        return freqs.astype(NPCT_REAL_TYPE, copy=False)


    @staticmethod
    def _modal_wave_type_to_int(wave_type):
        wave_type = _normalize_wave_type(wave_type)
        if wave_type == "Rayleigh":
            return wave_type, GRT_DISPERSION_RAYL
        return wave_type, GRT_DISPERSION_LOVE


    @staticmethod
    def _free_eigenv_roots(eigv_arr, nf):
        for iw in range(nf):
            eigv = eigv_arr[iw]
            if bool(eigv.c_roots):
                C_grt_free(cast(eigv.c_roots, c_void_p))
            if bool(eigv.u_roots):
                C_grt_free(cast(eigv.u_roots, c_void_p))
            if bool(eigv.c_roots_iref):
                C_grt_free(cast(eigv.c_roots_iref, c_void_p))


    def _set_modal_model_source_receiver(self, c_mod1d):
        c_mod1d.depsrc = self.depsrc
        c_mod1d.deprcv = self.deprcv
        c_mod1d.ircvup = self.depsrc >= self.deprcv

        deps = npct.as_array(c_mod1d.Dep, (c_mod1d.n,))
        isrc = 0
        for i in range(c_mod1d.n - 1):
            isrc = i
            if deps[i + 1] >= self.depsrc:
                break
        ircv = 0
        for i in range(c_mod1d.n - 1):
            ircv = i
            if deps[i + 1] >= self.deprcv:
                break
        c_mod1d.isrc = isrc
        c_mod1d.ircv = ircv


    def _default_modal_crange(self, c_mod1d):
        va = npct.as_array(c_mod1d.Va, (c_mod1d.n,))
        vb = npct.as_array(c_mod1d.Vb, (c_mod1d.n,))
        is_liquid = npct.as_array(c_mod1d.isLiquid, (c_mod1d.n,))

        vb_for_min = np.where(vb > 0.0, vb, va)
        vbmin = np.min(vb_for_min)
        vbn = vb[-1] if vb[-1] > 0.0 else va[-1]

        for i in range(c_mod1d.n - 1):
            if (not is_liquid[i]) and is_liquid[i + 1]:
                vbmin *= 1e-3
                break

        if vbmin >= vbn:
            raise ValueError(f"minimum velocity ({vbmin}) >= halfspace velocity ({vbn}).")

        return float(vbmin), float(vbn)


    def compute_dispersion(
        self,
        freqs=None,
        freq_range=None,
        wave_type:Literal['Rayleigh', 'Love']='Rayleigh',
        nmode:Union[int,None]=0,
        cmin:Union[float,None]=None,
        cmax:Union[float,None]=None,
        iref:Union[int,None]=None,
        satol:float=1e-3,
        cgap:float=1e-7,
        rtol:float=1e-4,
        vgap:float=1e-7,
        uniform_dc:float=0.0,
        print_progress:bool=False):
        r'''
            Compute phase-velocity dispersion curves of Rayleigh or Love waves.

            :param    freqs:       frequency points in Hz
            :param    freq_range:  tuple ``(f1, f2, df)``, matching CLI ``-Ff1/f2/df``
            :param    wave_type:   ``"Rayleigh"`` or ``"Love"``
            :param    nmode:       maximum mode order, 0 for fundamental mode, None for all roots

            :return:
                - **dispersion** - :class:`pygrt.modal.PyDispersion`
        '''
        if self.topbound != 'free' or self.botbound != 'halfspace':
            raise NotImplementedError(
                "Surface-wave dispersion currently only supports "
                "topbound='free' and botbound='halfspace'."
            )

        wave_type, wtype = self._modal_wave_type_to_int(wave_type)
        freqs = self._resolve_dispersion_freqs(freqs, freq_range)

        if nmode is not None:
            if int(nmode) < 0:
                raise ValueError("nmode must be nonnegative or None.")
            c_nmode = int(nmode) + 1
        else:
            c_nmode = 0

        if (cmin is None) != (cmax is None):
            raise ValueError("cmin and cmax must be set together.")
        if satol <= 0.0 or cgap <= 0.0 or rtol <= 0.0 or vgap <= 0.0:
            raise ValueError("satol, cgap, rtol and vgap must be positive.")
        if uniform_dc < 0.0:
            raise ValueError("uniform_dc must be nonnegative.")

        c_freqs = npct.as_ctypes(freqs)
        eigv_arr = (c_EIGENV*len(freqs))()
        eigmet = c_EIGENV_INFO()
        eigmet.nf = len(freqs)
        eigmet.freqs = c_freqs
        eigmet.nmode = c_nmode
        eigmet.modes = None
        eigmet.wtype = wtype
        eigmet.print_sec = False
        eigmet.satol = satol
        eigmet.uniform_dc = uniform_dc
        eigmet.cgap = cgap
        eigmet.rtol = rtol
        eigmet.vgap = vgap
        eigmet.eigv = eigv_arr

        c_mod1d_ptr = self._read_modal_model()
        try:
            c_mod1d = c_mod1d_ptr.contents
            if cmin is None:
                eigmet.cmin, eigmet.cmax = self._default_modal_crange(c_mod1d)
                eigmet.manual_crange = False
            else:
                if cmin <= 0.0 or cmax <= cmin:
                    raise ValueError("Require 0 < cmin < cmax.")
                eigmet.cmin = cmin
                eigmet.cmax = cmax
                eigmet.manual_crange = True

            if iref is not None:
                if int(iref) < 0 or int(iref) >= c_mod1d.n:
                    raise ValueError("iref is out of model layer range.")
                if wtype == GRT_DISPERSION_LOVE and c_mod1d.isLiquid[int(iref)]:
                    raise ValueError("Love wave iref cannot be in a liquid layer.")
                eigmet.iref = int(iref)
                eigmet.manual_iref = True

            C_grt_get_secular_roots(c_mod1d_ptr, pointer(eigmet), print_progress)

            max_roots = c_nmode if c_nmode > 0 else max([eigv_arr[i].n for i in range(len(freqs))] + [0])
            c = np.full((len(freqs), max_roots), np.nan, dtype=NPCT_REAL_TYPE)
            ciref = np.full((len(freqs), max_roots), -1, dtype=np.int64)
            cnum = np.zeros((len(freqs),), dtype=np.int64)

            for iw in range(len(freqs)):
                nroot = min(eigv_arr[iw].n, max_roots)
                cnum[iw] = eigv_arr[iw].n
                if nroot > 0:
                    c[iw, :nroot] = npct.as_array(eigv_arr[iw].c_roots, (eigv_arr[iw].n,))[:nroot]
                    ciref[iw, :nroot] = npct.as_array(eigv_arr[iw].c_roots_iref, (eigv_arr[iw].n,))[:nroot]
        finally:
            self._free_eigenv_roots(eigv_arr, len(freqs))
            C_grt_free_mod1d(c_mod1d_ptr)

        return PyDispersion(wave_type, freqs, c, ciref, cnum)


    @staticmethod
    def _validate_modal_fft_freqs(freqs):
        if len(freqs) < 2:
            raise ValueError("Modal summation requires at least two frequency points.")
        df = freqs[1] - freqs[0]
        if df <= 0.0:
            raise ValueError("dispersion frequencies must be increasing.")
        if not np.allclose(np.diff(freqs), df, rtol=1e-8, atol=1e-10):
            raise ValueError("Modal summation requires uniformly sampled frequencies.")
        if not np.isclose(freqs[0], df, rtol=1e-8, atol=1e-10):
            raise ValueError(
                "Modal summation requires dispersion generated like "
                "freq_range=(0, fmax, df)."
            )
        return float(df)


    def _build_eigenfn_info_for_modsum(self, dispersion, modes, freqband):
        all_freqs = dispersion.freqs.astype(NPCT_REAL_TYPE, copy=False)
        if freqband is None:
            freq_mask = np.ones(all_freqs.shape, dtype=bool)
        else:
            if len(freqband) != 2:
                raise ValueError("freqband must be [f1, f2].")
            f1, f2 = [float(x) for x in freqband]
            if f1 > f2:
                raise ValueError("freqband f1 must be <= f2.")
            freq_mask = (all_freqs >= f1) & (all_freqs <= f2)
        freq_idxs = np.nonzero(freq_mask)[0]
        if freq_idxs.size == 0:
            raise ValueError("freqband selects no frequencies.")

        if modes is None:
            modes = np.arange(dispersion.c.shape[1], dtype=np.int64)
        else:
            modes = np.asarray(modes, dtype=np.int64)
            if modes.ndim != 1:
                raise ValueError("modes must be a 1-D array.")
            if np.any(modes < 0) or np.any(modes >= dispersion.c.shape[1]):
                raise ValueError("modes contains out-of-range mode indexes.")

        nf = len(freq_idxs)
        nmode = len(modes)
        eigv_arr = (c_EIGENV*nf)()
        eigfn_rows = (POINTER(c_EIGENFN)*nf)()
        eigfn_storage = []
        c_roots_storage = []
        u_roots_storage = []
        ciref_storage = []

        for out_iw, src_iw in enumerate(freq_idxs):
            valid_modes = [int(mode) for mode in modes if mode < dispersion.cnum[src_iw]]
            nroot = len(valid_modes)
            c_vals = np.array([dispersion.c[src_iw, mode] for mode in valid_modes], dtype=NPCT_REAL_TYPE)
            iref_vals = np.array([dispersion.ciref[src_iw, mode] for mode in valid_modes], dtype=np.uint64)
            u_vals = np.zeros((nroot,), dtype=NPCT_REAL_TYPE)

            c_roots_storage.append(npct.as_ctypes(c_vals))
            u_roots_storage.append(npct.as_ctypes(u_vals))
            ciref_storage.append(npct.as_ctypes(iref_vals))

            eigv_arr[out_iw].c_roots = c_roots_storage[-1]
            eigv_arr[out_iw].u_roots = u_roots_storage[-1]
            eigv_arr[out_iw].c_roots_iref = ciref_storage[-1]
            eigv_arr[out_iw].n = nroot

            eigfn_storage.append((c_EIGENFN*nroot)())
            for ic in range(nroot):
                eigfn_storage[-1][ic].eigenC = c_vals[ic]
            eigfn_rows[out_iw] = eigfn_storage[-1]

        c_eigfn_freqs = npct.as_ctypes(all_freqs[freq_idxs].astype(NPCT_REAL_TYPE, copy=True))
        c_modes = npct.as_ctypes(modes.astype(np.uint64, copy=True))

        eigfnmet = c_EIGENFN_INFO()
        eigfnmet.nf = nf
        eigfnmet.freqs = c_eigfn_freqs
        eigfnmet.nmode = nmode
        eigfnmet.modes = c_modes
        eigfnmet.wtype = self._modal_wave_type_to_int(dispersion.wave_type)[1]
        eigfnmet.cpar_nz = 0
        eigfnmet.eigv = eigv_arr
        eigfnmet.eigfn = eigfn_rows

        keepalive = (
            c_eigfn_freqs, c_modes, eigv_arr, eigfn_rows,
            eigfn_storage, c_roots_storage, u_roots_storage, ciref_storage,
        )
        return eigfnmet, keepalive


    def compute_surface_grn(
        self,
        dispersion:PyDispersion,
        distarr:Union[np.ndarray,List[float],float],
        modes:Union[np.ndarray,List[int],None]=None,
        freqband:Union[np.ndarray,List[float],None]=None,
        upsampling_n:int=1,
        delayT0:float=0.0,
        delayV0:float=0.0,
        skipImagComps:bool=False,
        calc_upar:bool=False,
        gf_source=['EX', 'VF', 'HF', 'DC']):
        r'''
            Compute surface-wave Green's functions using modal summation.

            :param    dispersion:    :class:`pygrt.modal.PyDispersion`
            :param    distarr:       epicentral distances (km), or a single float
            :param    modes:         mode indexes to sum, default all available columns
            :param    freqband:      optional frequency range ``[f1, f2]`` in Hz

            :return:
                - **dataLst** - Green's Functions at multiple distances, in a list of :class:`obspy.Stream`
        '''
        if not isinstance(dispersion, PyDispersion):
            raise TypeError("dispersion must be a PyDispersion instance.")
        if self.topbound != 'free' or self.botbound != 'halfspace':
            raise NotImplementedError(
                "Surface-wave modal summation currently only supports "
                "topbound='free' and botbound='halfspace'."
            )
        if upsampling_n <= 0:
            raise ValueError("upsampling_n must be positive.")

        df = self._validate_modal_fft_freqs(dispersion.freqs)
        fft_nf = len(dispersion.freqs) + 1
        fft_nt = 2*(fft_nf - 1)
        fft_dt = 0.5 / dispersion.freqs[-1]
        fft_freqs = (np.arange(fft_nf)*df).astype(NPCT_REAL_TYPE)
        c_fft_freqs = npct.as_ctypes(fft_freqs)

        if isinstance(distarr, float) or isinstance(distarr, int):
            distarr = np.array([distarr*1.0])
        distarr = np.array(distarr, dtype=NPCT_REAL_TYPE, copy=True)
        if np.any(distarr < 0.0):
            raise ValueError("distarr < 0")
        distarr[distarr == 0.0] = 1e-5

        pygrnLst, c_grnArr = self._init_grn(distarr, fft_nt, fft_dt, upsampling_n, fft_freqs, 0.0, '')
        pygrnLst_uiz = []
        c_grnArr_uiz = None
        pygrnLst_uir = []
        c_grnArr_uir = None
        if calc_upar:
            pygrnLst_uiz, c_grnArr_uiz = self._init_grn(distarr, fft_nt, fft_dt, upsampling_n, fft_freqs, 0.0, 'z')
            pygrnLst_uir, c_grnArr_uir = self._init_grn(distarr, fft_nt, fft_dt, upsampling_n, fft_freqs, 0.0, 'r')

        c_rs = npct.as_ctypes(distarr.astype(NPCT_REAL_TYPE, copy=True))
        grn = c_GRNSPEC()
        grn.nf = fft_nf
        grn.freqs = c_fft_freqs
        grn.nf1 = 0
        grn.nf2 = fft_nf - 1
        grn.nr = len(distarr)
        grn.rs = c_rs
        grn.wI = 0.0
        grn.keepAllFreq = True
        grn.calc_upar = calc_upar
        grn.u = c_grnArr
        grn.uiz = c_grnArr_uiz
        grn.uir = c_grnArr_uir
        grn.statsstr = None
        grn.nstatsidxs = 0
        grn.statsidxs = None

        eigfnmet, keepalive = self._build_eigenfn_info_for_modsum(dispersion, modes, freqband)
        _ = keepalive

        c_mod1d_ptr = self._read_modal_model()
        try:
            c_mod1d = c_mod1d_ptr.contents
            self._set_modal_model_source_receiver(c_mod1d)
            eigfnmet.cpar_nz = c_mod1d.n
            C_grt_modsum_grn_spec(
                c_mod1d_ptr,
                self._modal_wave_type_to_int(dispersion.wave_type)[1],
                pointer(eigfnmet),
                pointer(grn),
            )
        finally:
            C_grt_free_mod1d(c_mod1d_ptr)

        return self._get_stream_from_grn_spectra(
            distarr, pygrnLst, pygrnLst_uiz, pygrnLst_uir,
            delayT0, delayV0, skipImagComps, calc_upar, gf_source
        )

    
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
    
    def _get_grn_spectra(
        self, 
        distarr:np.ndarray, 
        nt:int, 
        dt:float, 
        upsampling_n:int = 1,
        freqband:Union[np.ndarray,List[float]]=[-1,-1],
        zeta:float=0.8, 
        keepAllFreq:bool=False,
        vmin_ref:float=0.0,
        keps:float=-1.0,  
        ampk:float=1.15,
        k0:float=50.0, 
        k0_is_fixed:bool=False,
        Length:float=0.0, 
        filonLength:float=0.0,
        safilonTol:float=0.0,
        filonCut:float=0.0,
        converg_method:Literal['AUTO', 'NONE', 'DCM', 'PTAM']='AUTO',
        delayT0:float=0.0,
        delayV0:float=0.0,
        calc_upar:bool=False,
        statsfile:Union[str,None]=None, 
        statsidxs:Union[np.ndarray,List[int],None]=None, 
        print_runtime:bool=True
    ):
        
        depsrc = self.depsrc
        deprcv = self.deprcv

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


        # ====================================================================
        KPROC = c_K_INTEG_PROCESS()
        hs = max(abs(depsrc - deprcv), 1.0)
        KPROC.k0 = k0 * np.pi / hs
        KPROC.k0_is_fixed = k0_is_fixed
        KPROC.ampk = ampk
        KPROC.keps = keps if converg_method.upper() != 'AUTO' else 0.0
        KPROC.vmin = vmin_ref

        KPROC.kcut = filonCut / rmax
        
        KPROC.dk = 2.0*np.pi / (Length * rmax)
        
        KPROC.applyFIM = filonLength > 0.0
        KPROC.filondk = 2.0*np.pi / (filonLength * rmax) if filonLength > 0.0 else 0.0
        
        KPROC.applySAFIM = safilonTol > 0.0
        KPROC.sa_tol = safilonTol

        KPROC.cvgmet = K_INTEG_CVGMET_DICT[converg_method.upper()]
        # ====================================================================


        # ====================================================================
        grn = c_GRNSPEC()
        grn.nf = nf
        grn.freqs = c_freqs
        grn.nf1 = nf1
        grn.nf2 = nf2
        grn.nr = nrs
        grn.rs = c_rs
        grn.wI = wI
        grn.keepAllFreq = keepAllFreq
        grn.calc_upar = calc_upar
        grn.u = c_grnArr
        grn.uiz = c_grnArr_uiz
        grn.uir = c_grnArr_uir
        grn.statsstr = c_statsfile
        grn.nstatsidxs = nstatsidxs
        grn.statsidxs = c_statsidxs
        # ====================================================================

        # 运行C库函数
        #/////////////////////////////////////////////////////////////////////////////////
        # 计算得到的格林函数的单位：
        #     单力源 HF[ZRT],VF[ZR]                  1e-15 cm/dyne
        #     爆炸源 EX[ZR]                          1e-20 cm/(dyne*cm)
        #     剪切源 DD[ZR],DS[ZRT],SS[ZRT]          1e-20 cm/(dyne*cm)
        #=================================================================================
        C_grt_integ_grn_spec(self.c_mod1d, pointer(KPROC), pointer(grn), print_runtime)
        #=================================================================================
        #/////////////////////////////////////////////////////////////////////////////////

        return pygrnLst, pygrnLst_uiz, pygrnLst_uir

    def _get_stream_from_grn_spectra(
        self, distarr, pygrnLst, pygrnLst_uiz, pygrnLst_uir, 
        delayT0:float=0.0,
        delayV0:float=0.0,
        skipImagComps:bool=False,
        calc_upar:bool=False,
        gf_source=['EX', 'VF', 'HF', 'DC']
    ):
        depsrc = self.depsrc
        deprcv = self.deprcv
        
        calc_EX:bool = 'EX' in gf_source
        calc_VF:bool = 'VF' in gf_source
        calc_HF:bool = 'HF' in gf_source
        calc_DC:bool = 'DC' in gf_source

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
        for ir in range(len(distarr)):
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
                    stream.append(pygrnLst[ir][im][c].freq2time(delayT, travtP, travtS, sgn, skipImagComps))
                    if(calc_upar):
                        stream.append(pygrnLst_uiz[ir][im][c].freq2time(delayT, travtP, travtS, sgn*(-1), skipImagComps))
                        stream.append(pygrnLst_uir[ir][im][c].freq2time(delayT, travtP, travtS, sgn     , skipImagComps))


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
        k0:float=50.0, 
        k0_is_fixed:bool=False,
        Length:float=0.0, 
        filonLength:float=0.0,
        safilonTol:float=0.0,
        filonCut:float=0.0,
        converg_method:Literal['AUTO', 'NONE', 'DCM', 'PTAM']='AUTO',
        delayT0:float=0.0,
        delayV0:float=0.0,
        skipImagComps:bool=False,
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
            :param    vmin_ref:      minimum reference velocity (km/s). the default vmin=max(minimum velocity, 0.1), used to define the upper bound of k integral
            :param    keps:          automatic convergence condition, see Yao and Harkrider (1983) for more details.
                                     negative value denotes not use.
            :param    ampk:          The factor that affect the upper bound of the k integral, see below.
            :param    k0:            k0 used to define the maximum offset of upper bound :math:`\tilde{k_{max}}=\sqrt{(k_{0}*\pi/hs)^2 + (ampk*w/vmin_{ref})^2}` , hs=max(abs(depsrc-deprcv),1.0)
            :param    k0_is_fixed:      directly use k0, rather than choosing a proper offset in [0, k0]
            :param    Length:        integration step `dk=2\pi / (L*rmax)`, see Bouchon (1981) and 张海明 (2021) for the criterion, default set automatically.
            :param    filonLength:   integration step of Fixed-Interval Filon's Integration Method
            :param    safilonTol:    precision of Self-Adaptive Filon's Integration Method
            :param    filonCut:      The splitting point of DWM and (SA)FIM, k*=<filonCut>/rmax, default is 0
            :param    converg_method:   The method of explicit convergence, you can set "AUTO", "NONE", "DCM" or "PTAM". Default use "AUTO".
            :param    skipImagComps:    skip the amplitude compensation from imaginary frequency.
            :param    calc_upar:     whether calculate the spatial derivatives of displacements.
            :param    gf_source:     The source type to be calculated
            :param    statsfile:     directory path for saving the statsfile during k integral, used to debug or observe the variations of :math:`F(k,\omega)` and :math:`F(k,\omega)J_m(kr)k`    
            :param    statsidxs:     only output the statsfile at specific frequency indexes. It is recommended to specify the indexes; 
                                     otherwise, by default, statsfiles of all frequency will be output, which probably occupy a lot of disk space
            :param    print_runtime: whether print runtime and some other infomation.

            :return:
                - **dataLst** -   Green's Functions at multiple distances, in a list of :class:`obspy.Stream`
                
        '''

        if isinstance(distarr, float) or isinstance(distarr, int):
            distarr = np.array([distarr*1.0]) 

        distarr = np.array(distarr)
        distarr = distarr.copy().astype(NPCT_REAL_TYPE)

        pygrnLst, pygrnLst_uiz, pygrnLst_uir = self._get_grn_spectra(
            distarr, nt, dt, upsampling_n, freqband, zeta, keepAllFreq, 
            vmin_ref, keps, ampk, k0, k0_is_fixed, Length, filonLength, safilonTol, filonCut, converg_method,
            delayT0, delayV0, calc_upar,
            statsfile, statsidxs, print_runtime
        )

        dataLst = self._get_stream_from_grn_spectra(
            distarr, pygrnLst, pygrnLst_uiz, pygrnLst_uir,
            delayT0, delayV0, skipImagComps, calc_upar, gf_source
        )
        
        return dataLst  

    

    def compute_static_grn(
        self,
        xarr:Union[np.ndarray,List[float],float,None]=None, 
        yarr:Union[np.ndarray,List[float],float,None]=None, 
        distarr:Union[np.ndarray,List[float],float,None]=None, 
        keps:float=-1.0,  
        k0:float=50.0, 
        k0_is_fixed:bool=False,
        Length:float=15.0, 
        filonLength:float=0.0,
        safilonTol:float=0.0,
        filonCut:float=0.0,
        converg_method:Literal['AUTO', 'NONE', 'DCM', 'PTAM']='AUTO',
        calc_upar:bool=False,
        statsfile:Union[str,None]=None):

        r"""
            Call the C function to calculate the static Green's functions and return them in a dict.
            There're two ways to define the "epicentral distances":
            1. set both "xarr" and "yarr" to define a XY grid in advance.
            2. simply set "distarr", which equal to "xarr=[0.0], yarr=distarr".

            :param       xarr:          coordinate array in the north direction (km), or a single float.
            :param       yarr:          coordinate array in the east direction (km), or a single float.
            :param    distarr:          equal to "xarr=[0.0], yarr=distarr"
            :param       keps:          automatic convergence condition, see (Yao and Harkrider (1983) for more details.
                                        negative value denotes not use.
            :param       k0:            k0 used to define the maximum offset of upper bound :math:`\tilde{k_{max}}=(k_{0}*\pi/hs)^2`, hs=max(abs(depsrc-deprcv),1.0)
            :param       k0_is_fixed:      directly use k0, rather than choosing a proper offset in [0, k0]
            :param       Length:        integration step `dk=2\pi / (L*rmax)`, default L=15
            :param       filonLength:   integration step of Fixed-Interval Filon's Integration Method
            :param       safilonTol:    precision of Self-Adaptive Filon's Integration Method
            :param       filonCut:      The splitting point of DWM and (SA)FIM, k*=<filonCut>/rmax, default is 0
            :param    converg_method:   The method of explicit convergence, you can set "AUTO", "NONE", "DCM" or "PTAM". Default use "AUTO".
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

        depsrc = self.depsrc
        deprcv = self.deprcv

        if distarr is not None:
            if isinstance(distarr, float) or isinstance(distarr, int):
                distarr = np.array([distarr*1.0])
            distarr = np.array(distarr)
    
            xarr = np.array([0.0])
            yarr = distarr.copy()
            if np.any(yarr < 0.0):
                raise ValueError("distances can't be negative.")
        
        if xarr is None or yarr is None:
            raise ValueError("you need to set xarr and yarr or distarr.")

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
        

        # ====================================================================
        KPROC = c_K_INTEG_PROCESS()
        hs = max(abs(depsrc - deprcv), 1.0)
        KPROC.k0 = k0 * np.pi / hs
        KPROC.k0_is_fixed = k0_is_fixed
        KPROC.keps = keps if converg_method.upper() != 'AUTO' else 0.0

        # 最大震中距
        rmax = np.max(rs)
        KPROC.kcut = filonCut / rmax
        KPROC.dk = 2.0*np.pi / (Length * rmax)
        
        KPROC.applyFIM = filonLength > 0.0
        KPROC.filondk = 2.0*np.pi / (filonLength * rmax) if filonLength > 0.0 else 0.0
        
        KPROC.applySAFIM = safilonTol > 0.0
        KPROC.sa_tol = safilonTol

        KPROC.cvgmet = K_INTEG_CVGMET_DICT[converg_method.upper()]
        # ====================================================================



        # 运行C库函数
        #/////////////////////////////////////////////////////////////////////////////////
        # 计算得到的格林函数的单位：
        #     单力源 HF[ZRT],VF[ZR]                  1e-15 cm/dyne
        #     爆炸源 EX[ZR]                          1e-20 cm/(dyne*cm)
        #     剪切源 DD[ZR],DS[ZRT],SS[ZRT]          1e-20 cm/(dyne*cm)
        #=================================================================================
        C_grt_integ_static_grn(
            self.c_mod1d, nr, c_rs, pointer(KPROC),
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
