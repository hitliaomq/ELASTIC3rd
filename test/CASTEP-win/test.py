#!python
#
import os
import sys
import numpy as np

CurPath = os.path.dirname(__file__)
CurPath = CurPath.split('\\')
#print CurPath
ParPath = "\\".join(CurPath[0:-2])
if not ParPath in sys.path:
    sys.path.append(ParPath)

import elastic3rd
import elastic3rd.elastic

elastic3rd.elastic.elastic3()

