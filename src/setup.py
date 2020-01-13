#!/usr/bin/env python
#

from setuptools import setup

try:
    import numpy as np
    import scipy
except ImportError :
    raise ImportError("Numpy and Scipy must be installed.")

#files_erun = ["energy/energyrun"]

setup(
    name = "ELASTIC3rd",
    version = "2.4.2",
    description = "Claculate 3rd order elastic constants for crystals",
    author = "Mingqing Liao",
    author_email = "liaomq1900127@163.com",
    url = 'https://github.com/hitliaomq/ELASTIC3rd',
    download_url = 'https://github.com/hitliaomq/ELASTIC3rd',
    license = "GPL3",
    platforms = ['linux', 'windows'],
    keywords = ['physics', 'materials', 'elastic constants'],
    packages = ['elastic3rd', 'elastic3rd.crystal', 'elastic3rd.energy', 'elastic3rd.post', 'elastic3rd.symmetry'],
    package_dir = {'elastic3rd': 'elastic3rd', 'elastic3rd.energy' : 'elastic3rd/energy', 
                   'elastic3rd.symmetry' : 'elastic3rd/symmetry', 'elastic3rd.crystal' : 'elastic3rd/crystal', 
                   'elastic3rd.post' : 'elastic3rd/post'},
    package_data = {'elastic3rd.energy' : ['energyrun']},
    install_requires = ['numpy', 'scipy', 'matplotlib'],
    #long_description = read('README.md'),    
    classifiers = [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Education",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    zip_safe = False
    )