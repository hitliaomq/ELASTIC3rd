#!/usr/bin/env python
# Metadata is declared in pyproject.toml. This shim keeps legacy entry points
# such as "python setup.py install" working.
from setuptools import setup

if __name__ == "__main__":
    setup()
