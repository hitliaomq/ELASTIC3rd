#!python
#
import os
import sys

__name__ = "elastic3rd"
__author__ = "Mingqing Liao"
__version__ = "2.5.0"
__codedate__ = "2020-10-16"

CurPath = os.path.dirname(__file__)
CurPath = CurPath.strip("\n")
if not CurPath in sys.path:
    sys.path.append(CurPath)
#print(CurPath)