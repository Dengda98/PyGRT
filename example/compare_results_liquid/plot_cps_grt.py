import numpy as np 
import matplotlib.pyplot as plt 
from obspy import * 

def plot(st_grt:Stream, st_cps:Stream, out:str):
    nt = st_grt[0].stats.npts
    dt = st_grt[0].stats.delta
    t = np.arange(0, nt)*dt

    fig, axs = plt.subplots(6, 3, figsize=(12, 15), gridspec_kw=dict(hspace=0.0, wspace=0.1)) # 
    axs = axs.ravel()

    itr = 0
    for src, src2 in zip(['EX', 'VF', 'HF', 'DD', 'DS', 'SS'],
                        ['Explosion', 'Vertical\nForce', 'Horizontal\nForce', '45°\nDip Slip', '90°\nDip Slip', 'Strike Slip']):
        for comp, comp2 in zip(['Z', 'R', 'T'], ['Vertical', 'Radial', 'Transverse']):
            ax = axs[itr]
            
            for (_,spine) in ax.spines.items():
                spine.set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])

            if itr < 3:
                ax.set_title(comp2, fontsize=18, pad=10)

            itr += 1
            try:
                tr = st_grt.select(channel=f'{src}{comp}')[0]
                tr_cps = st_cps.select(channel=f'{comp}{src}')[0]
                data1 = tr.data.copy()
                data2 = tr_cps.data.copy()
            except:
                # data1 = data2 = np.zeros_like(t)
                continue
            
            
            norm = np.linalg.norm(data1, np.inf)+1e-10

            # adjust sign
            if comp == 'T':
                data2 *= -1
            if src == 'DS':
                data2 *= -1


            ax.plot(t, data2/norm, c='k', ls='-', lw=1.8)
            ax.plot(t, data1/norm, c='r', ls='-', lw=0.4)

            ax.set_xlim([10, 60])
            ax.set_ylim([-1.1,1.1])

            ylen = ax.get_ylim()[1] - ax.get_ylim()[0]

            if comp == 'Z':
                ax.text(0.8, 0, src2, fontsize=18, clip_on=False, ha='center', va='center', )


    for i in range(1, 4):        
        ax = axs[-i]
        ax.spines['bottom'].set_visible(True)
        ax.set_xticks(np.arange(10, 70, 10))
        ax.set_xlabel('Time (s)', fontsize=12)

    fig.savefig(out)




if __name__=='__main__':
    st_grt = read("GRN/liquid_20_8.2_100/*")
    st_cps = read("liquid_sdep20_rdep8.2/GRN/*")

    plot(st_grt, st_cps, "compare_cps_grt.pdf")