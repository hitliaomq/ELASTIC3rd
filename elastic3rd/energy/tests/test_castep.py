#!python
#
#Test for castep.py
#
#This is part of elastic3rd
#

import elastic3rd.energy.castep as eglue
import numpy as np


def test_parse_param():
    param_file = 'C.param'
    params = eglue.parse_param(param_file)
    print(params)

def test_write_param():
    param_file = 'C.param'
    param_new_file = 'C-NEW.param'
    params = eglue.parse_param(param_file)
    eglue.write_param(param_new_file, params)

def test_check_param():
    BaseName = 'C-NEW'
    Settings = {'task' : 'SinglePoint', 'calculate_stress' : 'False'}
    eglue.check_param(BaseName, Settings=Settings)

def test_parse_cell():
    cell_file = 'C.cell'
    cell_param = eglue.parse_cell(cell_file)
    print(cell_param)

def test_write_cell():
    cell_file = 'C.cell'
    cell_new_file = 'C-NEW.cell'
    cell_param = eglue.parse_cell(cell_file)
    eglue.write_cell(cell_new_file, cell_param)

def test_get_base_vec():
    BaseName = 'C'
    BaseVec = eglue.get_base_vec(BaseName)
    print(BaseVec)

def test_remove_comments():
    commemt_line = ['#This is a coments\n', ' ! This is a comments', 'COMMENT This is a comments', ' ;This is a comments']
    for linei in commemt_line:
        print(eglue.remove_comments(linei))
        assert(eglue.remove_comments(linei) is None)
    enmpty_line = ['  ', '\n', '']
    for linei in enmpty_line:
        print(eglue.remove_comments(linei))
        assert(eglue.remove_comments(linei) is None)
    normal_line = ['FIX_ALL_CELL : true', 'KPOINTS_MP_GRID 8 8 8\n', '%ENDBLOCK POSITIONS_FRAC', 'ANYTHINGS']
    for linei in normal_line:
        print(eglue.remove_comments(linei))
    normal_line_with_comments = ['%ENDBLOCK POSITIONS_FRAC#This is comments', 'KPOINTS_MP_GRID 8 8 8! This is comments', 'FIX_ALL_CELL : true ; this is comments']
    for linei in normal_line_with_comments:
        print(eglue.remove_comments(linei))


def test_write_basevec():
    BaseName = 'C-NEW'
    BaseVec = np.array([[1, 0, 3], [0, 2, 0], [0, 0, 1]], dtype=np.float64)
    BaseVec = [[1.2, 2.3, 3.4], [3.4, 4.5, 5.6], [6.7, 7.8, 8.9]]
    eglue.write_base_vec(BaseName, BaseVec)

def test_parse_geom():
    geom_file = 'C.geom'
    geom_param = eglue.parse_geom(geom_file)
    for propertyi in geom_param:
        print(propertyi['Energy'])
        #print(propertyi['Lattice'])
        #print(propertyi['Stress'])
        #print(propertyi['Position'])

def test_get_stress():
    BaseName = 'C'
    Stress0 = eglue.get_stress(BaseName)
    Stress1 = eglue.get_stress2(BaseName)
    print(Stress0)
    print(Stress1)

def test_get_energy():
    BaseName = 'C'
    Energy0 = eglue.get_energy(BaseName)
    Energy1 = eglue.get_energy2(BaseName)
    print(Energy0)
    print(Energy1)

def test_get_pos():
    BaseName = 'C'
    pos0 = eglue.get_pos(BaseName)
    (atoms, pos1) = eglue.get_pos2(BaseName)
    print(pos0)
    print(pos1)
    print(atoms)

def test_get_cell_type():
    cell_value0 = 12
    (cell_type0, cell_value0) = eglue.get_cell_type(cell_value0)
    print(cell_type0)
    assert(cell_type0 == 'SINGLE')
    cell_value1 = 'test'
    (cell_type1, cell_value1) = eglue.get_cell_type(cell_value1)
    assert(cell_type1 == 'SINGLE')
    cell_value2 = [12, 23]
    (cell_type2, cell_value2) = eglue.get_cell_type(cell_value2)
    assert(cell_type2 == 'SINGLE')
    cell_value3 = [[2, 3], [3, 4]]
    (cell_type3, cell_value3) = eglue.get_cell_type(cell_value3)
    assert(cell_type3 == 'BLOCK')
    cell_value4 = np.array(cell_value2)
    (cell_type4, cell_value4) = eglue.get_cell_type(cell_value4)
    assert(cell_type4 == 'SINGLE')
    cell_value5 = np.array(cell_value3)
    (cell_type5, cell_value5) = eglue.get_cell_type(cell_value5)
    assert(cell_type5 == 'BLOCK')

def test_check_cell_type():
    cell_dict = {'TEST1': {'type': 'SINGLE', 'value': 12}}
    cell_dict = eglue.check_cell_type(cell_dict)
    print(cell_dict)
    cell_dict = {'TEST1': {'value': 12}}
    cell_dict = eglue.check_cell_type(cell_dict)
    print(cell_dict)
    cell_dict = {'TEST1': {'value': [[12]]}}
    cell_dict = eglue.check_cell_type(cell_dict)
    print(cell_dict)
    cell_dict = {'TEST1': {'type': 'SINGLE', 'value': [[12]]}}
    cell_dict = eglue.check_cell_type(cell_dict)
    print(cell_dict)
    cell_dict = {'TEST1': {'type': 'BLOCK', 'value': 12}}
    cell_dict = eglue.check_cell_type(cell_dict)
    print(cell_dict)
    cell_dict = {'TEST1': 12}
    cell_dict = eglue.check_cell_type(cell_dict)
    print(cell_dict)
    cell_dict = {'TEST1': [[12, 13], [23, 24]]}
    cell_dict = eglue.check_cell_type(cell_dict)
    print(cell_dict)

def test_check_cell():
    BaseName = 'C-NEW'
    Settings={'FIX_ALL_CELL': ['ture'], 'CELL_CONSTRAINTS': {},
        'KPOINT_MP_GRID': {'type': 'SINGLE', 'value': np.array([16, 16, 16])}, 'KPOINTS_LIST': {}}
    eglue.check_cell2(BaseName, Settings)

def test_update_input_for_geomopt():
    BaseName='C'
    FlagEG='E'
    eglue.update_input_for_geomopt(BaseName, FlagEG)

test_parse_cell()