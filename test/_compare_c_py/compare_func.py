import numpy as np
from obspy import *
from scipy.io import netcdf_file
from typing import List
import matplotlib.pyplot as plt


def compare3(st_py:Stream, c_prefix:str, ZNE:bool=False, dim2:bool=False):
    """return average relative error"""
    
    if ZNE:
        pattern = "[ZNE]"
    else:
        pattern = "[ZRT]"

    if dim2:
        pattern = pattern*2

    st_c = read(f"{c_prefix}{pattern}.sac")
    
    print(f"{c_prefix}{pattern}")

    error = 0.0
    nerr = 0
    for tr_c in st_c:
        tr_py = st_py.select(channel=tr_c.stats.channel)[0]
        if np.all(tr_c.data == 0.0) and np.all(tr_py.data == 0.0):
            continue

        rerr = np.sum(np.abs(tr_c.data - tr_py.data)) / np.mean(np.abs(tr_py.data))
        print(tr_c.stats.channel, rerr)

        error += rerr
        nerr += 1 

    
    # fig, axs = plt.subplots(len(st_c), 1, figsize=(10, 1.5*len(st_c)))
    # for i, tr_c in enumerate(st_c):
    #     ax = axs[i]
    #     t = np.arange(0, tr_c.stats.npts)*tr_c.stats.delta
    #     tr_py = st_py.select(channel=tr_c.stats.channel)[0]
    #     ax.plot(t, tr_c.data, label="c_"+tr_c.stats.channel)
    #     ax.plot(t, tr_py.data, label="py_"+tr_py.stats.channel)
    #     # 总误差水平线
    #     ax.hlines(np.sum(np.abs(tr_c.data - tr_py.data)), *ax.get_xlim())
    #     ax.legend()
    #     ax.grid()

    # plt.show()

    return error/nerr

def update_dict(resDct:dict, Dct0:dict, prefix:str):
    keys = resDct.keys()
    for k in Dct0.keys():
        if k in keys:
            continue

        resDct.update({f"{prefix}{k}": Dct0[k]})

def static_compare3(resDct:dict, c_prefix:str):
    """return average relative error"""

    print(c_prefix)

    # read nc
    with netcdf_file(c_prefix, mmap=False) as f:
        error = 0.0
        nerr = 0

        for k in f.variables:
            if k == 'north' or k == 'east':
                continue
            val1 = resDct[k]
            val2 = f.variables[k][:]
            
            if np.all(val1 == 0.0) and np.all(val2 == 0.0):
                continue

            rerr = np.sum(np.abs(val1 - val2)) / np.mean(np.abs(val2))
            print(k, rerr)

            error += rerr
            nerr += 1 

    # plot_static(resDct, c_prefix)

    return error/nerr

def plot_static(resDct:dict, c_prefix:str):
    n = len([k for k in resDct.keys() if k[0] != '_'])
    fig, axs = plt.subplots(n, 3, figsize=(8, 3*n))
    with netcdf_file(c_prefix, mmap=False) as f:
        norths = resDct['_xarr']
        easts  = resDct['_yarr']

        keys = f.variables
        keys.pop('north')
        keys.pop('east')

        for i, k in enumerate(keys):
            val1 = resDct[k]
            val2 = f.variables[k][:]
            vmin = np.min(val1)
            vmax = np.max(val1)

            norm = np.abs(val2)
            norm[norm == 0.0] = 1.0

            print(k, np.mean(np.abs(val1 - val2) / norm), np.max(np.abs(val1 - val2) / norm))

            pcm0 = axs[i, 0].pcolorfast(easts, norths, val1.T, vmin=vmin, vmax=vmax)
            pcm1 = axs[i, 1].pcolorfast(easts, norths, val2.T, vmin=vmin, vmax=vmax)
            pcm2 = axs[i, 2].pcolorfast(easts, norths, ((val1 - val2) / norm).T)

            axs[i, 0].set_title("py_" + k)
            axs[i, 1].set_title("c_" + k)

            fig.colorbar(pcm0)
            fig.colorbar(pcm1)
            fig.colorbar(pcm2)

    import os
    fig.savefig(os.path.basename(c_prefix)[:-3] + ".pdf", bbox_inches='tight')
