# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys, os, pathlib
import subprocess 
import setuptools_scm
sys.path.insert(0, os.path.abspath('../../'))  # 指向pygrt模块的根目录

html_last_updated_fmt = '%b %d, %Y'


def setup(app):
    app.add_css_file('my_theme.css')  





project = 'PyGRT'
copyright = '2025, Zhu Dengda'
author = 'Zhu Dengda'
# 引入版本信息
version = setuptools_scm.get_version(root='../..', relative_to=__file__)
release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    'sphinx.ext.viewcode',
    "recommonmark",
    "sphinx_markdown_tables",
    "sphinx_copybutton",
    "sphinx.ext.intersphinx",
    "breathe",
    'sphinx.ext.mathjax',
    "nbsphinx",
    'sphinx_tabs.tabs',
    'sphinxcontrib.mermaid',
    'sphinx.ext.imgconverter',
] 

sphinx_tabs_valid_builders = ['linkcheck']
sphinx_tabs_disable_tab_closing = True   # 选项卡不关闭

nbsphinx_allow_errors = True  # 在构建文档时允许 Notebook 中的错误
nbsphinx_execute = 'never'   # 不执行 Notebook，只是展示内容
nbsphinx_requirejs_path = ''  # 加上该选项，否则nbsphinx插件会导致mermaid无法使用zoom

# 使用mermaid绘制简易流程图
# mermaid_output_format = 'raw'
# mermaid_version = 'latest'  
# mermaid_d3_zoom = True

source_suffix = {
    '.rst': 'restructuredtext',
    # '.txt': 'markdown',
    '.md': 'markdown',
} 

myst_enable_extensions = [
    "tasklist",
    "deflist",
    "dollarmath",
]

# Breathe configuration
breathe_projects = {
    "h_PyGRT": "../doxygen_h/xml",
}
breathe_default_project = "h_PyGRT"

templates_path = ['_templates']
exclude_patterns = []

language = 'zh_CN'
locale_dirs = ['../locales/']  # 存放翻译文件的目录
gettext_uuid = False
gettext_compact = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
master_doc = 'index'  
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'analytics_anonymize_ip': False,
    'logo_only': False,  # True
    'version_selector': True,
    'language_selector': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}

html_logo = "../../figs/logo_transparent.png"
html_static_path = ['_static']
html_css_files = ['my_theme.css']
html_js_files = [
    'my_custom.js',
]


autodoc_default_options = {
    'member-order': 'bysource',
    'special-members': '__init__',
}


intersphinx_mapping = {
    'obspy': ('https://docs.obspy.org/', None),
    'numpy': ('https://numpy.org/doc/stable/', None)
}

# 别名
rst_epilog = f"""
.. |Stream| replace:: :class:`~obspy.core.stream.Stream`
.. |Trace| replace:: :class:`~obspy.core.trace.Trace`
"""

 