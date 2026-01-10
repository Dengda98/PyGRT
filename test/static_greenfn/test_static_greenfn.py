import numpy as np
import pygrt

depsrc=2
deprcv=0
xarr = np.arange(-3., 3.1, 0.2)
yarr = np.arange(-2., 2.1, 0.2)

modname="../milrow"

modarr = np.loadtxt(modname)

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)
stgrn = pymod.compute_static_grn(xarr, yarr)
stgrn = pymod.compute_static_grn(xarr, yarr, calc_upar=True)
stgrn = pymod.compute_static_grn(xarr, yarr, Length=20)

stgrn = pymod.compute_static_grn(xarr, yarr, Length=20, converg_method='DCM')
stgrn = pymod.compute_static_grn(xarr, yarr, Length=20, converg_method='PTAM')
stgrn = pymod.compute_static_grn(xarr, yarr, Length=20, converg_method='none')

stgrn = pymod.compute_static_grn(xarr, yarr, k0=4, keps=1e-3)
stgrn = pymod.compute_static_grn(xarr, yarr, k0=4, statsfile="stgrt_stats")

# boundary condition
pymod = pygrt.PyModel1D(modarr, depsrc, deprcv, topbound='free', botbound='free')
stgrn = pymod.compute_static_grn(xarr, yarr)
pymod = pygrt.PyModel1D(modarr, depsrc, deprcv, topbound='halfspace', botbound='free')
stgrn = pymod.compute_static_grn(xarr, yarr)
pymod = pygrt.PyModel1D(modarr, depsrc, deprcv, topbound='rigid', botbound='rigid')
stgrn = pymod.compute_static_grn(xarr, yarr)