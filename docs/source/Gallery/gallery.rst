:author: 朱邓达
:date: 2025-12-19

计算示例 Gallery
====================

以下为 PyGRT 的一些计算示例。绘图部分使用
`Matplotlib <https://matplotlib.org/>`_ 或 `GMT <https://www.generic-mapping-tools.org/>`_ 。
计算部分使用 PyGRT 的 Python 接口或可执行文件 :command:`grt` 。

.. grid:: 1 2 3 3

    .. jinja::

        {% for i in range(1, 10) %}
        {% set i = '%02d' % i %}
        .. grid-item-card:: :doc:`ex{{i}}/ex{{i}}`
            :padding: 1
            :link-type: doc
            :link: ex{{i}}/ex{{i}}

            .. figure:: ex{{i}}/cover.*
        {% endfor %}

.. toctree::
   :hidden:
   :glob:

   ex*/ex*
