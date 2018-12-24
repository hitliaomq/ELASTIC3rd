#!/usr/bin/env python
#

from setuptools import setup

try:
    import numpy as np
    import scipy
except ImportError :
    raise ImportError("Numpy and Scipy must be installed.")

files = ["energy/energyrun"]

setup(
    name = "ELASTIC3rd",
    version = "1.1",
    description = "Claculate 3rd order elastic constant",
    author = "Mingqing Liao",
    author_email = "liaomq1900127@163.com",
    url = 'https://github.com/hitliaomq',
    license = "GPL3",
    packages = ['elastic3rd', 'elastic3rd.crystal', 'elastic3rd.energy', 'elastic3rd.post', 'elastic3rd.symmetry'],
    install_requires = ['numpy', 'scipy', 'matplotlib'],
    #long_description = read('README.md'),
    package_data = {'elastic3rd' : files},
    classifiers = [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Education",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    zip_safe = False
    )