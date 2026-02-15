from setuptools import setup, find_packages
import os


# 读取版本号
def read_version():
    version_file = os.path.join('pygrt', '_version.py')
    with open(version_file) as f:
        exec(f.read())
    return 'v'+locals()['__version__']

def read_readme():
    with open("README.md", encoding="utf-8") as f:
        return f.read()

setup(
    name='pygrt-kit',
    version=read_version(),
    description='An Efficient and Integrated Python Package for Computing Synthetic Seismograms in a Layered Half-Space Model',
    author='Zhu Dengda',
    author_email='zhudengda@mail.iggcas.ac.cn',
    long_description=read_readme(),  
    long_description_content_type="text/markdown", 
    url="https://github.com/Dengda98/PyGRT",
    project_urls={
        "Documentation": "https://pygrt.readthedocs.io/zh-cn/latest/",
        "Source Code": "https://github.com/Dengda98/PyGRT",
    },
    license='GPL-3.0',

    packages=find_packages(),
    package_data={'pygrt': ['./C_extension/*']},
    include_package_data=True,

    install_requires=[
        'numpy>=1.20, <2.0',
        'scipy>=1.10',
        'matplotlib>=3.5',
        'obspy>=1.4',
        'setuptools'
    ],
    python_requires='>=3.6',
    zip_safe=False,  # not compress the binary file
)
