# -----------------------------------------------------------------
# BEGIN GRN
import numpy as np
import matplotlib.pyplot as plt
from typing import Union
import pygrt 

modarr = np.loadtxt("mod1")

pymod = pygrt.PyModel1D(modarr, depsrc=0.01, deprcv=0.0)

# 不指定statsidx，默认输出全部频率点的积分过程文件
# vmin_ref 显式给定参考速度（用于定义波数积分上限），避免使用PTAM
# Length 给定波数积分间隔dk
_ = pymod.compute_grn(distarr=[1], nt=500, dt=0.02, vmin_ref=0.1, Length=20, statsfile="pygrtstats")
# END GRN
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# BEGIN read
# 指定待采样的速度数组
vels = np.arange(0.1, 0.6, 0.001)

# 读取所有频率的核函数，并插值到vels
# 不指定ktypes，默认返回全部核函数，均以2D数组的形式保存，shape=(nfreqs, nvels)
kerDct = pygrt.utils.read_kernels_freqs("pygrtstats", vels)
print(kerDct.keys())
# dict_keys(['_vels', '_freqs', 'EXP_q0', 'EXP_w0', 'VF_q0', 'VF_w0', 'HF_q1', 'HF_w1', 'HF_v1', 'DC_q0', 'DC_w0', 'DC_q1', 'DC_w1', 'DC_v1', 'DC_q2', 'DC_w2', 'DC_v2'])
# END read
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# BEGIN plot
# 绘制图像
def plot_kernel(kerDct:dict, RorI:bool, out:Union[str,None]=None):
    funcRorI = np.real if RorI else np.imag

    ktypes = []
    for key in kerDct:
        if key[0] == '_':
            continue
        ktypes.append(key)

    srctypes = ['EX', 'VF', 'HF', 'DD', 'DS', 'SS']

    vels = kerDct['_vels']
    freqs = kerDct['_freqs']

    fig = plt.figure(figsize=(12, 3*len(srctypes)))
    gs = fig.add_gridspec(len(srctypes), 3)
    qwvLst = ['q', 'w', 'v']
    for ik, key in enumerate(ktypes):
        srctype = key.split("_")[0]
        qwv = key.split("_")[1][0]
        
        ax = fig.add_subplot(gs[srctypes.index(srctype), qwvLst.index(qwv)])

        # 对不同速度间取归一化
        data = kerDct[key].copy()
        data[...] = data/np.max(np.abs(data), axis=1)[:,None]

        pcm = ax.pcolormesh(freqs, vels, np.abs(funcRorI(data)).T, vmin=0, vmax=1, shading='nearest')
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Velocity (km/s)")
        ax.set_title(key)
        fig.colorbar(pcm, ax=ax)

    if RorI:
        fig.suptitle("Real parts of Kernels", fontsize=20, x=0.5, y=0.99)
    else:
        fig.suptitle("Imag parts of Kernels", fontsize=20, x=0.5, y=0.99)

    if out is not None:
        fig.tight_layout()
        fig.savefig(out, dpi=100)


plot_kernel(kerDct, False, "imag.png")
plot_kernel(kerDct, True, "real.png")
# END plot
# -----------------------------------------------------------------
