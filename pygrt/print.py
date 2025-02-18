"""
    :file:     print.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2025-2-10

    打印安装包所在位置

"""

import os

MAIN_DIR = os.path.dirname(__file__)

if __name__=='__main__':
    print("PyGRT installation directory: ", MAIN_DIR)
    print("PyGRT executable file directory: ", os.path.join(MAIN_DIR, 'C_extension/bin'))
    print("PyGRT dynamic library directory: ", os.path.join(MAIN_DIR, 'C_extension/lib'))