#!python
#
import os
import sys

__name__ = "elastic3rd"
__author__ = "Mingqing Liao"
__version__ = "7.0.2"
__codedate__ = "2026-10-15"

CurPath = os.path.dirname(__file__)
CurPath = CurPath.strip("\n")
if not CurPath in sys.path:
    sys.path.append(CurPath)
#print(CurPath)