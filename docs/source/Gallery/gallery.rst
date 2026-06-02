:author: 朱邓达
:date: 2025-12-19

计算示例 Gallery
====================

以下为 **PyGRT** 的一些计算示例。绘图部分使用 |Matplotlib| 或 |GMT| 。
计算部分使用 **PyGRT** 的 Python 接口或可执行文件 :command:`grt` 。
所有示例均在代码主页的 `docs/source/Gallery <https://github.com/Dengda98/PyGRT/tree/main/docs/source/Gallery>`_
路径下。

.. jinja::

    {% macro card(idx) %}
    .. grid-item-card:: :doc:`ex{{ idx }}/ex{{ idx }}`
        :class-item: cropped-grid
        :padding: 1
        :link-type: doc
        :link: ex{{ idx }}/ex{{ idx }}

        .. figure:: ex{{ idx }}/cover.*
    {% endmacro %}
    
    动态解
    ------------

    .. grid:: 4

        {{ card('01')  | indent(8, first=true) }}  

        {{ card('02')  | indent(8, first=true) }}

        {{ card('03')  | indent(8, first=true) }}

        {{ card('04')  | indent(8, first=true) }}

        {{ card('05')  | indent(8, first=true) }}

        {{ card('06')  | indent(8, first=true) }}

        {{ card('07')  | indent(8, first=true) }}

        {{ card('08')  | indent(8, first=true) }}

        {{ card('09')  | indent(8, first=true) }}

        {{ card('10')  | indent(8, first=true) }}
        
        {{ card('11')  | indent(8, first=true) }}

        {{ card('14')  | indent(8, first=true) }}

        {{ card('16')  | indent(8, first=true) }}

        {{ card('17')  | indent(8, first=true) }}

    静态解
    ------------

    .. grid:: 4

        {{ card('12')  | indent(8, first=true) }}

        {{ card('13')  | indent(8, first=true) }}

        {{ card('15')  | indent(8, first=true) }}


.. toctree::
   :hidden:
   :glob:

   ex*/ex*
