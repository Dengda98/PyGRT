.. include:: common_OPTs.rst_


ker2asc
==================

:简介: 将波数积分过程中的核函数文件转为 ASCII 文件


语法
-----------

**grt ker2asc** *kernelfile*


描述
--------

**ker2asc** 模块输出结果打印到标准输出，第一列为波数值，
后面每两列的的命名方式为  ``{srcType}_{q/w/v}``，
与 :doc:`/Tutorial/dynamic/gfunc` 部分介绍的积分公式中的核函数 :math:`q_m, w_m, v_m` 保持一致。
两列分别表示复数的实部和虚部。

示例
------------

以下将 :doc:`greenfn` 模块输出的某个频率的核函数文件转为文本文件，
再使用 GMT 绘制其中核函数 *EX_q* 的变化趋势。

.. gmtplot::

    # 示例模型
    cat > mod <<EOF
    0.0  4.0 2.7 2.5
    EOF

    grt greenfn -Mmod -D0/0 -N500/0.02 -OGRN -R10 -S50,100
    grt ker2asc GRN_grtstats/mod_0_0/K_0050_5.00000e+00 > stats

    gmt begin test png
        gmt plot stats -JX12c/3c -Bafg -BWSen -i0,1 -W0.5p,blue -l"Real"
        gmt plot stats -i0,2 -W0.5p,blue,-- -l"Imag"
        gmt legend -DjTL+w1.5c+o0.1c -F+gwhite+c0.2p+p1p
    gmt end show


