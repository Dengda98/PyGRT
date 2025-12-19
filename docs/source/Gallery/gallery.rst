Gallery
====================

.. grid:: 2 3 3 4

    .. jinja::

        {% for i in range(1, 3) %}
        {% set i = '%02d' % i %}
        .. grid-item-card:: :doc:`ex{{i}}/ex{{i}}`
            :padding: 1
            :link-type: doc
            :link: ex{{i}}/ex{{i}}

            .. figure:: ex{{i}}/cover.png
        {% endfor %}

.. toctree::
   :hidden:
   :glob:

   ex*/ex*
