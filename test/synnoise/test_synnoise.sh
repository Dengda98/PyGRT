#!/bin/bash

set -euo pipefail
source ../common.sh

cat > noise_sources.txt <<'INPUT'
FORCE_Z        -1 1
FORCE_AZIMUTH  0 360
STRENGTH       1e15 2e15
INTERVAL       0.5
SEED           0

# ID     north    east   depth
SRC001       0      -1     0.0
SRC002       1       2     0.5
INPUT

cat > noise_receivers.txt <<'INPUT'
DT      0.2
LENGTH  2

# ID     north    east   depth
STA001       0       0     0.0
STA002       0       1     0.5
INPUT

grt synnoise -h
grt synnoise -M../milrow -Cnoise_sources.txt -Unoise_receivers.txt -Onoise_cli -P2
# 短记录只保留首个事件，同时覆盖阶数范围和静默选项
sed 's/LENGTH  2/LENGTH  0.4/' noise_receivers.txt > noise_short.txt
grt synnoise -M../milrow -Cnoise_sources.txt -Unoise_short.txt -Onoise_clipped -N0/1 -P1 -s

# 未指定模型
expect_fail grt synnoise -Cnoise_sources.txt -Unoise_receivers.txt -Oinvalid
# 起始阶数不可大于结束阶数
expect_fail grt synnoise -M../milrow -Cnoise_sources.txt -Unoise_receivers.txt -Oinvalid -N1/0
# 记录长度必须对应整数采样点
sed 's/LENGTH  2/LENGTH  2.1/' noise_receivers.txt > noise_bad.txt
expect_fail grt synnoise -M../milrow -Cnoise_sources.txt -Unoise_bad.txt -Oinvalid
# 台站参数头不可缺失
sed '1,2d' noise_receivers.txt > noise_bad.txt
expect_fail grt synnoise -M../milrow -Cnoise_sources.txt -Unoise_bad.txt -Oinvalid
# 台站参数头必须按 DT、LENGTH 顺序给出
sed '1s/.*/LENGTH  2/;2s/.*/DT      0.2/' noise_receivers.txt > noise_bad.txt
expect_fail grt synnoise -M../milrow -Cnoise_sources.txt -Unoise_bad.txt -Oinvalid
# 采样间隔必须为正
sed 's/DT      0.2/DT      0/' noise_receivers.txt > noise_bad.txt
expect_fail grt synnoise -M../milrow -Cnoise_sources.txt -Unoise_bad.txt -Oinvalid
# 记录长度必须为正
sed 's/LENGTH  2/LENGTH  0/' noise_receivers.txt > noise_bad.txt
expect_fail grt synnoise -M../milrow -Cnoise_sources.txt -Unoise_bad.txt -Oinvalid
# 采样参数须通过台站文件指定
expect_fail grt synnoise -M../milrow -Cnoise_sources.txt -Unoise_receivers.txt -S0.2/2 -Oinvalid
# 并发进程数必须为正
expect_fail grt synnoise -M../milrow -Cnoise_sources.txt -Unoise_receivers.txt -Oinvalid -P0
# 输出路径不可包含空格
expect_fail grt synnoise -M../milrow -Cnoise_sources.txt -Unoise_receivers.txt '-Onoise output'
# 输出目录必须为空
expect_fail grt synnoise -M../milrow -Cnoise_sources.txt -Unoise_receivers.txt -Onoise_cli
# 参数头不可缺失
sed '1,5d' noise_sources.txt > noise_bad.txt
expect_fail grt synnoise -M../milrow -Cnoise_bad.txt -Unoise_receivers.txt -Oinvalid
# 激发间隔必须为正
sed 's/INTERVAL       0.5/INTERVAL       0/' noise_sources.txt > noise_bad.txt
expect_fail grt synnoise -M../milrow -Cnoise_bad.txt -Unoise_receivers.txt -Oinvalid

python -u test_synnoise.py

rm -rf noise_cli noise_clipped invalid
rm -f noise_sources.txt noise_receivers.txt noise_short.txt noise_bad.txt
