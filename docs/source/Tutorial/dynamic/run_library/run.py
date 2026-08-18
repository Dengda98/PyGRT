from pathlib import Path

import pygrt


MODEL = Path("milrow")


# BEGIN GRN
pymod = pygrt.PyModel1D(grn="GRN", modelpath=MODEL)
pymod.compute_grn(
    depsrc=[2.0, 4.0],
    deprcv=[0.0, 2.0],
    dists=[5.0, 8.0, 10.0],
    nt=256,
    dt=0.02,
)
# END GRN


# BEGIN SYN
pymod = pygrt.PyModel1D(grn="GRN")
pymod.compute_syn(
    depsrc=4.0,
    deprcv=2.0,
    dist=8.0,
    azimuth=30.0,
    scale=1e24,
    output_path="syn_python",
)
# END SYN
