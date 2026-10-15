#!python
#
#
import elastic3rd.esutils as esutils

def test_write_input():
    ParaIn = esutils.read_input('INPUT.e3rd')
    esutils.write_input(ParaIn, FileName = "INPUT-NEW")

test_write_input()