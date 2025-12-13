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

    "NPCT_REAL_TYPE",
    "NPCT_CMPLX_TYPE",

    "REAL",
    "PREAL",
    "PCPLX",

    "c_GRT_MODEL1D",
    "c_K_INTEG_METHOD"
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


NPCT_REAL_TYPE = 'f8'
NPCT_CMPLX_TYPE = 'c16'



REAL = c_double
CPLX = REAL*2
PREAL = POINTER(REAL)
PCPLX = POINTER(CPLX)

class c_RT_MATRIX(Structure):
    """
    和C结构体 RT_MATRIX 匹配
    """

    _fields_ = [
        ('RD', CPLX * 4),
        ('RU', CPLX * 4),
        ('TD', CPLX * 4),
        ('TU', CPLX * 4),
        ('RDL', CPLX),
        ('RUL', CPLX),
        ('TDL', CPLX),
        ('TUL', CPLX),
        ('invT', CPLX * 4),
        ('invTL', CPLX),
        ('stats', c_int)
    ]


class c_GRT_MODEL1D(Structure):
    """
    和C结构体 GRT_MODEL1D 作匹配
    """
    
    _fields_ = [
        ('n', c_size_t), 
        ("depsrc", REAL),
        ("deprcv", REAL),
        ('isrc', c_size_t),
        ('ircv', c_size_t),
        ('ircvup', c_bool),
        ('io_depth', c_bool),

        ('omega', CPLX),
        ('k', REAL),
        ('c_phase', CPLX),

        ('Thk', PREAL),
        ('Dep', PREAL),
        ('Va', PREAL),
        ('Vb', PREAL),
        ('Rho', PREAL),
        ('Qa', PREAL),
        ('Qb', PREAL),
        ('Qainv', PREAL),
        ('Qbinv', PREAL),

        ('mu', PCPLX),
        ('lambda', PCPLX),
        ('delta', PCPLX),
        ('atna', PCPLX),
        ('atnb', PCPLX),
        ('xa', PCPLX),
        ('xb', PCPLX),
        ('caca', PCPLX),
        ('cbcb', PCPLX),

        ('stats', c_int),

        ('M_AL', c_RT_MATRIX),
        ('M_BL', c_RT_MATRIX),
        ('M_RS', c_RT_MATRIX),
        ('M_FA', c_RT_MATRIX),
        ('M_FB', c_RT_MATRIX),

        ('M_top', c_RT_MATRIX),

        ('R_EV', CPLX * 4),
        ('R_EVL', CPLX),
        ('uiz_R_EV', CPLX * 4),
        ('uiz_R_EVL', CPLX),

        ('src_coef', CPLX * SRC_M_NUM * QWV_NUM * 2)
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