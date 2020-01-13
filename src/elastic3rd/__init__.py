#!python
#
import os
import sys

CurPath = os.path.dirname(__file__)
CurPath = CurPath.strip("\n")
if not CurPath in sys.path:
    sys.path.append(CurPath)
#print(CurPath)