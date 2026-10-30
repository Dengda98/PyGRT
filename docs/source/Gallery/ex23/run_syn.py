import pygrt

pymod = pygrt.PyModel1D(grn="GRN")

# 不再额外细分子断层，合成地表接收点的位移及其空间导数
pymod.syn(
    src_fault="fault.inr", src_fault_size=(0.0, 0.0),
    dist=200.0, azimuth=10.0, deprcv=0.0,
    integrate_order=1, zne=True, calc_upar=True, output_path="syn"
)

# 计算应力、接收断层上的投影应力和库伦应力
pygrt.tensor.stress("syn")
pygrt.tensor.sproj("syn", strike=90, dip=32, rake=90)
pygrt.tensor.coulomb("syn", friction=0.4)
