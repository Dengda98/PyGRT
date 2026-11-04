"""读取已有源台文件，合成 ZNE 背景噪声"""

# BEGIN
import pygrt

model = pygrt.PyModel1D(modelpath="model4")
model.synnoise(source_file="sources.txt", receiver_file="stations.txt",
               output_path="ambient_run", modes=(0, 10))
# END
