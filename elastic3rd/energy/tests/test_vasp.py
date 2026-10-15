#!python
#

import elastic3rd.energy.vasp as eglue
import numpy as np

def test_get_pos():
    pos = eglue.get_pos()
    print(pos)

def test_parse_poscar():
    poscar = eglue.parse_poscar(filename='POSCAR')
    print(poscar)

def test_write_poscar():
    poscar = eglue.parse_poscar(filename='POSCAR')
    eglue.write_poscar(poscar, filename='POSCAR-NEW')

def test_add_perturbation():
    eglue.add_perturbation(tor=2E-10, filename='POSCAR-NEW')

test_add_perturbation()
