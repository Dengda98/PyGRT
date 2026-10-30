"""
    :file:     utils.py
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
    :date:     2024-07-24

    旧 utils 模块的弃用报错入口

"""

def __getattr__(name):
    if name.startswith("__"):
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    raise ImportError(
        "pygrt.utils is deprecated and has been split into multiple modules. "
        "Please see the Python API documentation: "
        "https://pygrt.readthedocs.io/zh-cn/latest/API/py_api.html"
    )
