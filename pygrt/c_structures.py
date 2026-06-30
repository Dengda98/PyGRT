"""
    :file:     c_structures.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包括  
        1、模型结构体的C接口 c_PyModel1D  
        2、格林函数结构体的C接口 c_GRN  

"""


from ctypes import *

__all__ = [
    "CHANNEL_NUM",
    "QWV_NUM",
    "INTEG_NUM",
    "SRC_M_NUM",
    "SRC_M_ORDERS",
    "SRC_M_NAME_ABBR",
    "ZRTchs",
    "ZNEchs",
    "qwvchs",
    "MECHANISM_NUM",
    "GRT_DISPERSION_RAYL",
    "GRT_DISPERSION_LOVE",
    "GRT_RAYL_DIM",
    "GRT_LOVE_DIM",
    "GRT_EGYINTS_MAX",
    "GRT_SNSTVTY_MAX",

    "NPCT_REAL_TYPE",
    "NPCT_CMPLX_TYPE",

    "K_INTEG_CVGMET_DICT",

    "REAL",
    "PREAL",
    "CPLX",
    "PCPLX",

    "c_MODEL1D",
    "c_K_INTEG_PROCESS",
    "c_GRNSPEC",
    "c_EIGENV",
    "c_EIGENV_INFO",
    "c_EIGENFN",
    "c_EIGENFN_INFO",
]


CHANNEL_NUM = 3
QWV_NUM = 3
INTEG_NUM = 4
SRC_M_NUM = 6
SRC_M_ORDERS = [0, 0, 1, 0, 1, 2]
SRC_M_NAME_ABBR = ["EX", "VF", "HF", "DD", "DS", "SS"]
ZRTchs = ['Z', 'R', 'T']
ZNEchs = ['Z', 'N', 'E']
qwvchs = ['q', 'w', 'v']
MECHANISM_NUM = 6
GRT_DISPERSION_RAYL = 0
GRT_DISPERSION_LOVE = 1
GRT_RAYL_DIM = 4
GRT_LOVE_DIM = 2
GRT_EGYINTS_MAX = 4
GRT_SNSTVTY_MAX = 3

NPCT_REAL_TYPE = 'f8'
NPCT_CMPLX_TYPE = 'c16'

K_INTEG_CVGMET_DICT = {
    'AUTO': 0,
    'NONE': 1,
    'DCM': 2,
    'PTAM': 3
}


class CPLX(Structure):
    _fields_ = [
        ('real', c_double),
        ('imag', c_double),
    ]


REAL = c_double
PREAL = POINTER(REAL)
PCPLX = POINTER(CPLX)


class c_MODEL1D(Structure):
    """
    和C结构体 MODEL1D 作匹配
    """
    
    _fields_ = [
        ('n', c_size_t), 
        ("depsrc", REAL),
        ("deprcv", REAL),
        ('isrc', c_size_t),
        ('ircv', c_size_t),
        ('ircvup', c_bool),
        ('io_depth', c_bool),
        ('srcrcv_isInserted', c_bool),
        ('omgref', CPLX),

        ('Thk', PREAL),
        ('Dep', PREAL),
        ('Va', PREAL),
        ('Vb', PREAL),
        ('Rho', PREAL),
        ('Qa', PREAL),
        ('Qb', PREAL),
        ('Qainv', PREAL),
        ('Qbinv', PREAL),
        ('isLiquid', POINTER(c_bool)),

        ('topbound', c_int),
        ('botbound', c_int),
    ]


class c_K_INTEG_PROCESS(Structure):
    """
    和C结构体 K_INTEG_PROCESS 作匹配
    """

    _fields_ = [
        ('k0', REAL),
        ('k0_is_fixed', c_bool),
        ('ampk', REAL),
        ('keps', REAL),
        ('vmin', REAL),

        ('kcut', REAL),
        ('kmax', REAL),

        ('dk', REAL),

        ('applyFIM', c_bool),
        ('filondk', REAL),
        
        ('applySAFIM', c_bool),
        ('sa_tol', REAL),

        ('cvgmet', c_int),

        ('fstats', c_void_p),
        ('ptam_fstatsnr', c_void_p),
    ]

class c_GRNSPEC(Structure):
    """
    和 C 结构体 GRNSPEC 作匹配
    """

    _fields_ = [
        ('nf', c_size_t),
        ('freqs', PREAL),
        ('nf1', c_size_t),
        ('nf2', c_size_t),
        ('nr', c_size_t),
        ('rs', PREAL),
        ('wI', REAL),
        ('keepAllFreq', c_bool),
        ('calc_upar', c_bool),
        
        ('u', POINTER((PCPLX*CHANNEL_NUM)*SRC_M_NUM)),
        ('uiz', POINTER((PCPLX*CHANNEL_NUM)*SRC_M_NUM)),
        ('uir', POINTER((PCPLX*CHANNEL_NUM)*SRC_M_NUM)),

        ('statsstr', c_char_p),
        ('nstatsidxs', c_size_t),
        ('statsidxs', POINTER(c_size_t)),
    ]


class c_EIGENV(Structure):
    """
    和 C 结构体 EIGENV 作匹配
    """

    _fields_ = [
        ('c_roots', PREAL),
        ('u_roots', PREAL),
        ('c_roots_iref', POINTER(c_size_t)),
        ('n', c_size_t),
    ]


class c_EIGENV_INFO(Structure):
    """
    和 C 结构体 EIGENV_INFO 作匹配
    """

    _fields_ = [
        ('nf', c_size_t),
        ('freqs', PREAL),
        ('nmode', c_size_t),
        ('modes', POINTER(c_size_t)),
        ('wtype', c_int),
        ('print_sec', c_bool),

        ('cmin', REAL),
        ('cmax', REAL),
        ('manual_crange', c_bool),
        ('iref', c_size_t),
        ('manual_iref', c_bool),

        ('satol', REAL),
        ('uniform_dc', REAL),

        ('cgap', REAL),
        ('rtol', REAL),
        ('vgap', REAL),

        ('neval', c_size_t),

        ('eigv', POINTER(c_EIGENV)),
    ]


class c_EIGENFN(Structure):
    """
    和 C 结构体 EIGENFN 作匹配
    """

    _fields_ = [
        ('eigenC', REAL),
        ('fn', POINTER(CPLX*4)),
        ('egyint', CPLX*GRT_EGYINTS_MAX),
        ('csens', POINTER(CPLX*GRT_SNSTVTY_MAX)),
        ('usens', POINTER(CPLX*GRT_SNSTVTY_MAX)),
    ]


class c_EIGENFN_INFO(Structure):
    """
    和 C 结构体 EIGENFN_INFO 作匹配
    """

    _fields_ = [
        ('nf', c_size_t),
        ('freqs', PREAL),
        ('nmode', c_size_t),
        ('modes', POINTER(c_size_t)),
        ('wtype', c_int),

        ('nz', c_size_t),
        ('zs', PREAL),
        ('z_irefs', POINTER(c_size_t)),

        ('cpar_nz', c_size_t),
        ('cpar_zs', PREAL),
        ('cpar_z_irefs', POINTER(c_size_t)),

        ('eigv', POINTER(c_EIGENV)),
        ('eigfn', POINTER(POINTER(c_EIGENFN))),
    ]
