import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf_file

def plot_disp(fLst):
    fig, axs = plt.subplots(3, 3, figsize=(12, 8))
    for i, f in enumerate(fLst):
        ax3 = axs[i]
        chs = ['E', 'N', 'Z']
        for k in range(3):
            ch = chs[k]
            ax = ax3[k]
            val = f.variables[ch][:, 0]
            Y = f.variables['north'][:]
            ax.plot(Y, val, 'bo-', ms=3, label=chs[k])
            ax.grid()
            ax.legend(loc='upper right')

    axs[0, 0].set_ylabel("Strike-slip", fontsize=12)
    axs[1, 0].set_ylabel("Dip-slip", fontsize=12)
    axs[2, 0].set_ylabel("Tensile", fontsize=12)
    return fig, axs

def plot_stress(fLst):
    fig, axs = plt.subplots(2, 3, figsize=(12, 7), gridspec_kw=dict(hspace=0.2))
    for i, f in enumerate(fLst):
        ax = axs[0, i]
        Y = f.variables['north'][:]
        ax.plot(Y, f.variables['stress_ZZ'][:,0], 'bo-', ms=3, label='ZZ')
        ax.plot(Y, f.variables['stress_NN'][:,0], 'ro-', ms=3, label='NN')
        ax.plot(Y, f.variables['stress_EE'][:,0], 'ko-', ms=3, label='EE')
        ax.grid()

        ax = axs[1, i]
        ax.plot(Y, f.variables['stress_ZN'][:,0], 'bo-', ms=3, label='ZN')
        ax.plot(Y, f.variables['stress_ZE'][:,0], 'ro-', ms=3, label='ZE')
        ax.plot(Y, f.variables['stress_NE'][:,0], 'ko-', ms=3, label='NE')
        ax.grid()

    axs[0, 0].set_title("Strike-slip", fontsize=12)
    axs[0, 1].set_title("Dip-slip", fontsize=12)
    axs[0, 2].set_title("Tensile", fontsize=12)

    axs[0, 0].set_ylabel("Normal", fontsize=12)
    axs[1, 0].set_ylabel("Shear", fontsize=12)

    axs[0, 2].legend(bbox_to_anchor=(1.3, 1))
    axs[1, 2].legend(bbox_to_anchor=(1.3, 1))
    return fig, axs


# ======================================================================
with netcdf_file("stsyn_1_stkslip.nc", mmap=False) as f1, \
     netcdf_file("stsyn_1_dipslip.nc", mmap=False) as f2, \
     netcdf_file("stsyn_1_tensile.nc", mmap=False) as f3:
    fig, axs = plot_disp([f1, f2, f3])
    fig.savefig("disp1.svg", bbox_inches='tight')
    fig, axs = plot_stress([f1, f2, f3])
    fig.savefig("stress1.svg", bbox_inches='tight')

with netcdf_file("stsyn_2_stkslip.nc", mmap=False) as f1, \
     netcdf_file("stsyn_2_dipslip.nc", mmap=False) as f2, \
     netcdf_file("stsyn_2_tensile.nc", mmap=False) as f3:
    fig, axs = plot_disp([f1, f2, f3])
    fig.savefig("disp2.svg", bbox_inches='tight')
    fig, axs = plot_stress([f1, f2, f3])
    fig.savefig("stress2.svg", bbox_inches='tight')

with netcdf_file("stsyn_3_stkslip.nc", mmap=False) as f1, \
     netcdf_file("stsyn_3_dipslip.nc", mmap=False) as f2, \
     netcdf_file("stsyn_3_tensile.nc", mmap=False) as f3:
    fig, axs = plot_disp([f1, f2, f3])
    fig.savefig("disp3.svg", bbox_inches='tight')
    fig, axs = plot_stress([f1, f2, f3])
    fig.savefig("stress3.svg", bbox_inches='tight')