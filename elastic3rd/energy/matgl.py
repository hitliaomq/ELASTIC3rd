#!python
#for matgl: https://github.com/materialsvirtuallab/matgl
import os, re
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from monty.serialization import loadfn
from ase.units import GPa

from elastic3rd.energy.ml_ff import (
    get_struct_from_file, 
    get_base_vec, 
    write_base_vec, 
    get_energy, 
    get_stress,
    copy_files,
    update_structure,
    update_input_for_geomopt,
    add_perturbation,
    rem_files
)

import matgl
from matcalc import RelaxCalc
from matgl.ext.ase import PESCalculator
from ase.filters import FrechetCellFilter, ExpCellFilter, StrainFilter, UnitCellFilter

filter_dict = {'ExpCellFilter': ExpCellFilter, 'FrechetCellFilter': FrechetCellFilter,
               'StrainFilter': StrainFilter, 'UnitCellFilter': UnitCellFilter}

def relax_struct(BaseName=None):
    #device = "cuda" if torch.cuda.is_available() else "cpu"
    if os.path.exists('settings.yaml'):
        settings = loadfn('settings.yaml')
    else:
        settings = {}
    optimizer = settings.get('optimizer', 'BFGS')
    ase_filter = settings.get('ase_filter', 'ExpCellFilter')
    ase_filter = filter_dict[ase_filter]
    constrain_symmetry = settings.get('constrain_symmetry', True)
    model = settings.get('model', 'TensorNet-MatPES-PBE-v2025.1-PES')
    steps = settings.get('seps', 800)
    flageg = settings.get('flageg', 'G')
    (struct, filename) = get_struct_from_file(BaseName=BaseName)

    atoms = AseAtomsAdaptor.get_atoms(struct)
    if flageg == 'E':
        relax_cell = False
        relax_atoms = True
    elif flageg == 'G':
        relax_atoms = True
        relax_cell = True
    elif flageg == 'EO':
        relax_atoms = False
        relax_cell = False
    else:
        relax_cell = False
        relax_atoms = True

    pot = matgl.load_model(model)
    ase_calc = PESCalculator(pot)
    calculator = RelaxCalc(ase_calc, optimizer=optimizer, max_steps=steps,
        relax_atoms=relax_atoms, relax_cell=relax_cell, cell_filter=ase_filter)
    results = calculator.calc(struct)

    return (results, filename)

def run(NP=None, BaseName=None):
    '''
    Run string for vasp. 
    Parameters
    ----------
        BaseName: str
            Any string, just for consist with castep
        NP: int
            NP is the total cores
    Return
    ------
        RunStr: str
            The string for calling first principles code
    '''
    (results, filename) = relax_struct(BaseName=BaseName)
    energy = results['energy']
    stress = results['stress']
    stress = np.array(stress) #/ GPa
    if len(stress.shape) == 1:
        stress = np.array([[stress[0], stress[5], stress[4]], 
                           [stress[5], stress[1], stress[3]], 
                           [stress[4], stress[3], stress[2]]])
    np.savetxt('ENERGY', [energy])
    np.savetxt('STRESS', stress)
    struct = results['final_structure']
    flnm_lst = re.split(r'\.', filename)
    if len(flnm_lst) == 1:
        #POSCAR
        flnm_lst.append('relax')
    else:
        flnm_lst.insert(-1, 'relax')
    filename = '.'.join(flnm_lst)
    struct.to(filename=filename)
