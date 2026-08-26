import numpy as np
import pygrt

pymod = pygrt.PyModel1D(stgrn="stgrn.nc")

# 参考点经纬度
lat0 = -8.3101000
lon0 = 121.3517000

# 合成
synout = "stsyn_ff.nc"
pymod.compute_static_syn(
    norths=(-200, 200, 10), easts=(-210, 210, 10),
    src_fault="coulomb-fault.inp", deprcv=0.0,
    output_path=synout, calc_upar=True, zne=True
)
pygrt.utils.compute_stress(synout)
pygrt.utils.compute_sproj(synout, strike=90, dip=32, rake=90)
pygrt.utils.compute_coulomb(synout, friction=0.4)
# 坐标转换
pygrt.utils.xy2geo(ingrid=synout, lat0=lat0, lon0=lon0, outgrid=f"{synout}.geo")

# okada 解作为对比
synout = "okada_ff.nc"
pygrt.utils.compute_okada(
    modelparams=(6, 3.464, 2.7),
    norths=(-200, 200, 10), easts=(-210, 210, 10),
    src_fault="coulomb-fault.inp", deprcv=0.0,
    output_path=synout, calc_upar=True, zne=True
)
pygrt.utils.compute_stress(synout)
pygrt.utils.compute_sproj(synout, strike=90, dip=32, rake=90)
pygrt.utils.compute_coulomb(synout, friction=0.4)
# 坐标转换
pygrt.utils.xy2geo(ingrid=synout, lat0=lat0, lon0=lon0, outgrid=f"{synout}.geo")

