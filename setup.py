#!/usr/bin/env python
#

PACKAGE = "elastic3rd"
NAME = "elastic3rd"
DESCRIPTION = "Claculate nth order elastic constants for crystals."
AUTHOR = "Mingqing Liao"
AUTHOR_EMAIL = "liaomq1900127@163.com"
URL = "https://github.com/hitliaomq/elastic3rd"
VERSION = __import__(PACKAGE).__version__

from setuptools import setup, find_packages
import os
from io import open

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()



data_file_coef = ""

setup(
    name = NAME,
    version = VERSION,
    description = DESCRIPTION,
    long_description = LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url = URL,
    author = AUTHOR,
    author_email = "liaomq1900127@163.com",  
    license = "GPL3",

    packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    install_requires=['numpy', 'scipy', 'matplotlib', 'pyyaml'],

    package_data = {
        'elastic3rd.energy' : ['energyrun', 'atomcfg.yaml'], 
        'elastic3rd.symmetry' : ['COEFS/*/*/*', 'COEFS/*/*', 'SymData/*/*'],
    },
    
    entry_points={
        'console_scripts': [
            'elastic3rd = elastic3rd.scripts.run_elastic3rd:run_e3rd',
        ]
    },

    classifiers = [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Education",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.7',
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords = ['physics', 'materials', 'elastic constants'],
    project_urls={
        'Documentation': 'https://elastic3rd.readthedocs.io/en/latest',
        'Author scholar hompage': 'https://scholar.google.com/citations?user=rfEuIUUAAAAJ',
    },

    zip_safe = False,
)