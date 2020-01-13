#!python

import os
import sys

CurPath = os.path.dirname(__file__)
CurPath = CurPath.split('\\')
#print CurPath
ParPath = "\\".join(CurPath[0:-1])
if not ParPath in sys.path:
    sys.path.append(ParPath)