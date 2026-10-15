#!python
from pymatgen.io.pwscf import PWInput, PWOutput
from pymatgen.core.structure import Structure
import elastic3rd.energy.glue as glue

import numpy as np

def get_base_vec(BaseName=None):
    #Get the base vector
    scfin = PWInput.from_file('{}.scf'.format(BaseName))
    structure = scfin.structure
    basevec = structure.lattice.matrix
    return basevec

def write_base_vec(BaseName=None, BaseVec=None, F=None):
    scfin = PWInput.from_file('{}.scf'.format(BaseName))
    scfin.structure.lattice = BaseVec
    scfin.write_file(filename='{}.scf'.format(BaseName))

def add_perturbation(BaseName=None, tor=2.0E-10, filename='POSCAR', PertMethod='random', PertType='relative'):
    pass

def get_energy(BaseName=None):
    outobj = PWOutput(filename='{}.out'.format(BaseName))
    energy = outobj.final_energy
    Energy = glue.multi_energy(energy)    
    return Energy

def get_stress(BaseName=None):
    pass

def copy_files(BaseName=None, Path=None):
    pass

def update_input_for_geomopt(BaseName=None, FlagEG='G', KeepSym=False):
    if FlagEG == 'G':
        incar_setting = {"ISIF": "3", "NSW": "100"}
    elif FlagEG == 'E':
        incar_setting = {"ISIF": "2", "NSW": "100"}
    elif FlagEG == 'GS':
        incar_setting = {'ISIF': '7', "NSW": "100"}
    elif FlagEG == 'GV':
        incar_setting = {'ISIF': '4', "NSW": "100"}
    elif FlagEG == 'GI':
        incar_setting = {'ISIF': '6', "NSW": "100"}
    elif FlagEG == 'EO':
        incar_setting = {'ISIF': '2', "NSW": "0", 'IBRION': '-1'}

def update_structure(BaseName=None, dstpath="./"):
    pass

def rem_files(BaseName=None):
    pass

'''
get_base_vec()

basevec = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
write_base_vec(BaseName=None, BaseVec=basevec)

get_stress(BaseName=None)
'''