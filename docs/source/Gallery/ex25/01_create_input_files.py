"""生成本例的震源文件和台站文件"""

import numpy as np

from pygrt import noise


# 输出文件
SOURCE_FILE  = "sources.txt"
STATION_FILE = "stations.txt"

# 本例的源台分布及采样参数
STATION_COUNT   = 20                        # 台站数
STATION_SPACING = 2.0                       # 台站间距，km
DT              = 0.1                       # 采样间隔，s
LENGTH          = 36000.0                   # 记录长度，s
SOURCE_COUNT    = 10000                     # 震源数，每个源只激发一次
RADIUS_RANGE    = (60.0, 100.0)             # 圆环源区的内外半径，km

# 随机单力参数
FORCE_Z         = (-1.0, 1.0)               # 单位力的向下分量范围
FORCE_AZIMUTH   = (0.0, 360.0)              # 单力水平角范围，从北顺时针计，度
STRENGTH        = (1e15, 2e15)              # 单力强度范围，dyne
SEED            = 20261031                  # 随机种子


# 以下生成本例的源台坐标，实际使用时可以替换为自己的坐标或输入文件
# 坐标依次为北、东、向下深度，单位为 km

# 台站沿东西向布设，台阵中点位于原点
stations = np.zeros((STATION_COUNT, 3))
stations[:, 1] = (np.arange(STATION_COUNT) - (STATION_COUNT - 1) / 2) * STATION_SPACING

# 震源在以台阵中点为圆心的圆环内按面积均匀分布，各源等间隔激发一次
sources = noise.random_points(SOURCE_COUNT, radius_range=RADIUS_RANGE, seed=SEED)

noise.write_sources(SOURCE_FILE, sources, interval=LENGTH / SOURCE_COUNT, force_z_range=FORCE_Z,
                    force_azimuth_range=FORCE_AZIMUTH, strength_range=STRENGTH, seed=SEED, fmt=".8g")
noise.write_receivers(STATION_FILE, stations, dt=DT, length=LENGTH, fmt=".8g")
