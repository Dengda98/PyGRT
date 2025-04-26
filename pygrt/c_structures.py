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
    "USE_FLOAT",
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

    "c_PyModel1D",
]


USE_FLOAT = False
CHANNEL_NUM = 3
QWV_NUM = 3
INTEG_NUM = 4
SRC_M_NUM = 6
SRC_M_ORDERS = [0, 0, 1, 0, 1, 2]
SRC_M_NAME_ABBR = ["EX", "VF", "HF", "DD", "DS", "SS"]
ZRTchs = ['Z', 'R', 'T']
ZNEchs = ['Z', 'N', 'E']
qwvchs = ['q', 'w', 'v']


NPCT_REAL_TYPE = 'f4' if USE_FLOAT else 'f8'
NPCT_CMPLX_TYPE = f'c{int(NPCT_REAL_TYPE[1:])*2}'



REAL = c_float if USE_FLOAT else c_double
PREAL = POINTER(REAL)
PCPLX = POINTER(REAL*2)

class c_PyModel1D(Structure):
    """
    和C结构体PYMODEL1D作匹配

    :field n:        层数
    :filed depsrc:   震源深度 km
    :filed deprcv:   接收点深度 km
    :field isrc:     震源所在层位
    :field ircv:     台站所在层位
    :field ircvup:   台站层位是否高于震源 
    :field thk:      数组, 每层层厚(km)
    :field Va:       数组, 每层P波速度(km/s)
    :field Vb:       数组, 每层S波速度(km/s)
    :field Rho:      数组, 每层密度(g/cm^3)
    :field Qa:       数组, 每层P波品质因子Q_P
    :field Qb:       数组, 每层S波品质因子Q_S

    """
    _fields_ = [
        ('n', c_int), 
        ("depsrc", REAL),
        ("deprcv", REAL),
        ('isrc', c_int),
        ('ircv', c_int),
        ('ircvup', c_bool),

        ('Thk', PREAL),
        ('Va', PREAL),
        ('Vb', PREAL),
        ('Rho', PREAL),
        ('Qa', PREAL),
        ('Qb', PREAL),
    ]
