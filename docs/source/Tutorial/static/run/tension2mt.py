"""
    Zhu Dengda
    
    将张位错转为矩张量表示
"""
from math import *
import numpy as np
from scipy.io import netcdf_file
import sys

strike = radians(float(sys.argv[1]))
dip = radians(float(sys.argv[2]))
ncfile = sys.argv[3]

nvec = np.array([
    - sin(dip) * sin(strike),
      sin(dip) * cos(strike),
      cos(dip),
])

vvec = nvec.copy()

# 从 nc 文件中读取 src_va, src_vb
with netcdf_file(ncfile, mmap=False) as f:
    src_va = f._attributes['src_va']
    src_vb = f._attributes['src_vb']

M = ((src_va/src_vb)**2 - 2) * np.eye(3) * np.sum(nvec * vvec) + (np.einsum('i,j', nvec, vvec) + np.einsum('i,j', vvec, nvec))

# 输出为 mrr   mtt   mff   mrt   mrf   mtf
Mrtf = [
    M[2,2], M[0,0], M[1,1], M[0,2], -M[1,2], -M[0,1]
]
print(*Mrtf, sep=' ', end='')
