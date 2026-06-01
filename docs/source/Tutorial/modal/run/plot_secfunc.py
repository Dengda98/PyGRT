import numpy as np
import matplotlib.pyplot as plt

def plot_secfunc(path, title, ax):

    data = np.loadtxt(path)

    carr = data[:,0]

    ax.plot(carr, data[:,1], 'k-', lw=1, marker='o', ms=0, label='Real')
    ax.plot(carr, data[:,2], 'k--', lw=1, marker='o', ms=0, label='Imaginary')

    ax.set_xmargin(0)
    ax.set_xlabel("Phase velocity (km/s)")
    ax.set_ylabel("Secular function")
    ax.set_title(title)
    ax.legend(loc='upper right', handletextpad=0.3, borderaxespad=0.3, labelspacing=0.2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3), layout='constrained')

plot_secfunc("secfunc_R", "Rayleigh", ax1)
plot_secfunc("secfunc_L", "Love", ax2)

fig.savefig("secfunc.svg")
