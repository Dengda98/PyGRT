import pygrt

# BEGIN PYTHON RCVFN
model = pygrt.PyModel1D(modelpath="mod1")

model.rcvfn(
    wtype="P",
    rayp=0.03,
    nt=500,
    dt=0.1,
    output_path="PY_P",
    alp=5.0,
    delay=-10.0,
    write_components=True,
    upsampling_n=5,
)

model.rcvfn(
    wtype="S",
    inca=10.0,
    idx=0,
    nt=500,
    dt=0.1,
    output_path="PY_S",
    alp=5.0,
    delay=-10.0,
    write_components=True,
    upsampling_n=5,
)
# END PYTHON RCVFN
