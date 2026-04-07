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

    "NPCT_REAL_TYPE",
    "NPCT_CMPLX_TYPE",

    "REAL",
    "PREAL",
    "CPLX",
    "PCPLX",

    "c_MODEL1D",
    "c_K_INTEG_METHOD",
    "c_GRNSPEC"
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

NPCT_REAL_TYPE = 'f8'
NPCT_CMPLX_TYPE = 'c16'


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


class c_K_INTEG_METHOD(Structure):
    """
    和C结构体 K_INTEG_METHOD 作匹配
    """

    _fields_ = [
        ('k0', REAL),
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

        ('applyDCM', c_bool),
        ('applyPTAM', c_bool),

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