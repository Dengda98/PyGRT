#!/bin/bash

set -euo pipefail

python plot_geometry.py

mkdir -p inputs results

# BEGIN FAULTDIV
# 从两种完整断层母版生成六种尺度的剖分，并设置滑动速率与子断层中心的破裂延迟
for dx in 10 5 2.5 2 1 0.5; do
    grt faultdiv -Fstrike_slip.inr -I${dx}/${dx} -R3.5+d0  -Dt/0/0.5/0 -Oinputs/strike_slip_dx${dx}.inp
    grt faultdiv -Fthrust.inr      -I${dx}/${dx} -R3.5+d90 -Dt/0/0.5/0 -Oinputs/thrust_dx${dx}.inp
done
# END FAULTDIV

python plot_faultdiv.py
bash plot_delay.sh

# BEGIN LAMB
# +i0/0 保留 faultdiv 的剖分，采用各行末尾的时间函数与破裂延迟
# -I1 对矩形滑动速率的响应积分，得到斜坡滑动对应的位移
for dx in 10 5 2.5 2 1 0.5; do
    grt lamb -H8/4.62/3.3 -N200/0.06 -Cinputs/strike_slip_dx${dx}.inp+i0/0 -R20 -A0 -Dr2 -I1 -Oresults/strike_slip_dx${dx}_z2
    grt lamb -H8/4.62/3.3 -N200/0.06 -Cinputs/thrust_dx${dx}.inp+i0/0      -R20 -A0 -Dr2 -I1 -Oresults/thrust_dx${dx}_z2
done

# 最细网格的地表接收点，用于绘制灰色对照曲线
grt lamb -H8/4.62/3.3 -N200/0.06 -Cinputs/strike_slip_dx0.5.inp+i0/0 -R20 -A0 -Dr0 -I1 -Oresults/strike_slip_dx0.5_z0
grt lamb -H8/4.62/3.3 -N200/0.06 -Cinputs/thrust_dx0.5.inp+i0/0      -R20 -A0 -Dr0 -I1 -Oresults/thrust_dx0.5_z0
# END LAMB

python plot.py

cp fig_5_8_3.svg cover.svg


ex=$(basename "$(pwd)")
rm -f "$ex.tar.gz"
cd .. && tar -czvf "$ex.tar.gz" "$ex" && mv "$ex.tar.gz" "$ex" && cd -

# 把输出和结果都打包，然后再删除
rm -rf inputs results
