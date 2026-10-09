import shutil
from unittest import TestCase

import pygrt

raises = TestCase().assertRaises

model = pygrt.PyModel1D(modelpath="../milrow")
model.rcvfn(wtype="P", rayp=0.12, nt=32, dt=0.1, output_path="receiver")
model.rcvfn(wtype="S", inca=20, idx=2, nt=32, dt=0.1, output_path="receiver", zeta=0.9, upsampling_n=2,
            keepAllFreq=True, skipImagComps=True, alp=0.8, delay=0.5, write_components=True)

# 水平射线参数不可大到使 P 波无法传播
with raises(RuntimeError):
    model.rcvfn(wtype="P", rayp=1, nt=32, dt=0.1, output_path="receiver")

# 不可同时设置水平射线参数和入射角
with raises(ValueError):
    model.rcvfn(wtype="P", rayp=0.12, inca=20, nt=32, dt=0.1, output_path="receiver")
# 未指定水平射线参数或入射角
with raises(ValueError):
    model.rcvfn(wtype="P", nt=32, dt=0.1, output_path="receiver")
# 入射角不可达到或超过 90 度
with raises(ValueError):
    model.rcvfn(wtype="P", inca=90, nt=32, dt=0.1, output_path="receiver")

shutil.rmtree("receiver")
