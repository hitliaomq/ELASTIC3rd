#!python
#
from elastic3rd.crystal.lattice import *
import elastic3rd.energy.castep as escastep
from castepkit.parse.inputs import Cell

def test_switchLatticAxis():
    BaseVec = [[0.000000000000000,       2.649562000000000,       0.000000000000000],
           [5.651847000000000,       0.000000000000000,       0.000000000000000],
           [0.000000000000000,       0.000000000000000,       2.733606000000000]]

    positions = [[0.7305050000000000,  0.3739310000000001,      0.8527239999999991],
                 [-0.7305050000000000, -0.3739310000000001,      1.3527239999999991],
                 [0.7305050000000000, -0.3739310000000001,     -0.3527239999999993],
                 [-0.7305050000000000,  0.3739310000000001,     -0.8527239999999991],
                 [0.5000000000000000,  0.1904570000000000,      0.4999999999999998],
                 [-0.5000000000000000, -0.1904570000000000,      0.9999999999999999],
                 [0.8395260000000000,  1.0000000000000000,      0.2499999999999999],
                 [-0.8395260000000000, -1.0000000000000000,      0.7499999999999999]]

    swax = [0, 2]
    (BaseVec, pos) = switchLatticAxis(swax, BaseVec, positions)
    print(BaseVec)
    print(pos)

'''
import shutil
new='y'
BaseName = 'C8-3'
BaseName_new = '{}-{}'.format(BaseName, new)
shutil.copyfile('{}.cell'.format(BaseName), '{}.cell'.format(BaseName_new))
cell = Cell.from_file('{}.cell'.format(BaseName_new))
print(cell)
basevec = cell.lattice
pos = cell.positions
swax = [1, 2]
(BaseVec, positions) = switchLatticAxis(swax, basevec, pos)
print(cell.cell_dict)
cell_dict = cell.cell_dict
#cell_dict.update({'LATTICE_CART': BaseVec, })
escastep.update_pos(BaseName_new, positions, atoms=None)
escastep.write_base_vec(BaseName_new, BaseVec)
'''

abc = [5.10977,  5.83587,  5.86552, 110.8101, 115.7885,  86.2458]
basevec = ABC2base_vec(abc)
print(basevec)