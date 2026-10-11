:author: 朱邓达
:date: 2026-10-10

.. include:: common_OPTs.rst_


faultdiv
========================

:简介: 剖分 |simplified-Coulomb| 格式的有限断层文件中的有限断层，并按需设置时间函数和破裂延迟


语法
-----------

**grt faultdiv**
|-F|\ *fault*
|-I|\ *dL/dW*
|-O|\ *outfile*
[ |-R|\ *vel*\ **+o**\ *along/down* | |-R|\ *vel*\ **+d**\ *degree* ]
[ |-D|\ *tftype*\ [/*tfparams*][**+d**\ *delay*] ]
[ **-h** ]


描述
--------

:doc:`faultdiv` 将 |simplified-Coulomb| 格式的有限断层文件中的矩形断层剖分为子断层，
输出仍采用同一格式，可直接用于合成。

仅支持 **Kode=100/200/300** 的矩形记录，保持原有滑动分量和开度；
滑动角/总滑动量输入会转换为等价的右旋/逆冲分量。
不接受 **Kode=400/500** 点源。每条输入断层均按整数块数等分，
子断层按先沿走向、再沿下倾排列。
**输入文件中 "bot" 列之后的内容始终全部忽略**。

仅指定 |-F|、|-I| 和 |-O| 时，输出每行只有 11 个数值列。
指定 |-R| 或 |-D| 后，每行末尾追加一个完整的 **-D** 选项，
其中延迟等于 |-D| 给出的整体偏移加上 |-R| 计算的传播时间。

该文件可直接用于 :doc:`syn` 或 :doc:`lamb`，合成时使用 **+i0/0** 保留已有剖分，
并省略全局 **-D**，以采用各行的时间函数。


必选选项
----------

.. _-F:

**-F**\ *fault*
    输入 |simplified-Coulomb| 格式的有限断层文件。

.. _-I:

**-I**\ *dL/dW*
    目标剖分尺寸，两个尺寸均须为正数。
    *dL* 为沿走向目标尺寸，单位为 km。
    *dW* 为沿断层面向下倾的目标尺寸，单位为 km。
    程序内为调整 *dL/dW* 使形成整数份的等距剖分。

.. _-O:

**-O**\ *outfile*
    输出剖分后的 |simplified-Coulomb| 格式的有限断层文件路径。


可选选项
----------

.. _-R:

**-R**\ *vel*\ **+o**\ *along/down* | **-R**\ *vel*\ **+d**\ *degree*
    恒定破裂速度 *vel*，单位为 km/s，须为正数。
    **+o** 和 **+d** 必须且只能设置一个：

    * **+o**\ *along/down*：从局部起裂点向四周传播，原点为断层顶边起点。
      **此时输入断层文件必须且只能包含一条断层记录**。
      *along* 为沿正走向的距离，单位为 km。
      *down* 为沿断层面向下倾的距离，单位为 km。
      起裂点须位于输入断层的断层面内。
      子断层传播时间为其中心到起裂点的断层面距离除以 *vel*。
    * **+d**\ *degree*：平面前沿沿给定方向传播。
      角度采用类似 rake 的定义，范围为 **−180°～180°**：
      0° 沿正走向，90° 向上倾，−90° 向下倾，±180° 沿反走向。
      零时刻前沿接触完整断层的上游边界，传播时间为中心沿该方向距前沿的投影距离除以 *vel*。
      允许多条输入记录，各断层按同一方向独立设置破裂前沿与传播时间。

    破裂从零时刻开始传播。
    设置 |-R| 而省略 |-D| 时，默认使用脉冲时间函数 **-Di**。

.. _-D:

**-D**\ *tftype*\ [/*tfparams*][**+d**\ *delay*]
    指定子断层的时间函数，语法与 :doc:`syn` 的 **-D** 相同。


示例
----------

设置平面破裂和点破裂的延迟，用 :doc:`fault2nc` 转为断层面网格，
再通过 |GMT| 绘制延迟等值线，详见以下脚本。

.. gmtplot::
    :language: bash
    :caption: 平面破裂（左）与点破裂（右）的延迟等值线，横轴沿走向，纵轴向下倾，时间单位为 s

    # 给出一个有限断层母版，沿走向长 20 km，沿倾向长 10 km
    cat > fault.inr <<'EOF'
    #  X-start  Y-start  X-fin  Y-fin  Kode  rake  slip  dip  top  bot
    #  -------  -------  -----  -----  ----  ----  ----  ---  ---  ---
    1  0        0        20     0      100   0     1     90   0    10
    EOF

    # 剖分为多个子断层，每个子断层设置不同的破裂延迟
    # 平面前沿沿 -30 度方向传播
    grt faultdiv -Ffault.inr -R3.0+d-30  -I0.5/0.5 -Odelay_plane.inp
    # 点破裂从局部坐标 (8, 3) km 开始
    grt faultdiv -Ffault.inr -R3.0+o8/3 -I0.5/0.5 -Odelay_point.inp

    # 规则剖分的断层可以将一些信息导出到 nc 文件
    grt fault2nc -Fdelay_plane.inp -Odelay_plane.nc
    grt fault2nc -Fdelay_point.inp -Odelay_point.nc

    # 绘图
    gmt begin delay png
        gmt set FONT_ANNOT_PRIMARY 11p FONT_LABEL 12p FONT_TITLE 14p MAP_FRAME_TYPE plain
        gmt subplot begin 1x2 -Fs10c/5c -M0.8c -R0/20/0/10 -JX10c/-5c
            # 横轴为走向距离，纵轴为下倾距离，保持两个方向的距离比例一致
            gmt subplot set 0
            gmt basemap -Bxa5f1+l"Along strike (km)" -Bya2f1+l"Down dip (km)" -BWSen+t"Plane rupture: delay (s)"
            gmt grdcontour "delay_plane.nc?delay" -C0.25 -A0.5+f10p -Gn1 -Wc0.6p,gray40 -Wa1p,black

            # 每 0.25 s 绘制等值线，每 0.5 s 标注延迟时间，+号标出起裂点
            gmt subplot set 1
            gmt basemap -Bxa5f1+l"Along strike (km)" -Bya2f1+l"Down dip (km)" -BWSen+t"Point rupture: delay (s)"
            gmt grdcontour "delay_point.nc?delay" -C0.25 -A0.5+f10p -Gn1 -Wc0.6p,gray40 -Wa1p,black
            echo "8 3" | gmt plot -S+0.25c -Gred -W0.5p,black
        gmt subplot end
    gmt end
