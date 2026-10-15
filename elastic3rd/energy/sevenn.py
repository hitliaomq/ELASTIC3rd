#!python
#for svenn: https://github.com/MDIL-SNU/SEVENNet
#
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
    rem_files,
    Relaxer
)

from sevenn.calculator import SevenNetCalculator
    
def relax_struct(BaseName=None):
    #device = "cuda" if torch.cuda.is_available() else "cpu"
    if os.path.exists('settings.yaml'):
        settings = loadfn('settings.yaml')
    else:
        settings = {}
    optimizer = settings.get('optimizer', 'BFGS')
    device = settings.get('device', 'cpu')
    ase_filter = settings.get('ase_filter', 'ExpCellFilter')
    constrain_symmetry = settings.get('constrain_symmetry', True)
    model = settings.get('model', '7net-mf-ompa')
    tag = model.split('-')[-1]
    if tag.lower() == 'ompa':
        modal_defalt = 'omat24'
    else:
        modal_defalt = None
    modal = settings.get('modal', modal_defalt)
    steps = settings.get('seps', 800)
    flageg = settings.get('flageg', 'E')
    (struct, filename) = get_struct_from_file(BaseName=BaseName)
    atoms = AseAtomsAdaptor.get_atoms(struct)
    if flageg == 'E':
        fix_axis = True
    elif flageg == 'G':
        fix_axis = False
    else:
        fix_axis = False
    atoms.calc = SevenNetCalculator(model=model, modal=modal)
    if flageg == 'EO':
        return (True, atoms, filename)
    else:
        relaxer = Relaxer(
            optimizer=optimizer,
            filter=ase_filter,
            constrain_symmetry=constrain_symmetry,
            fix_axis=fix_axis)
        (flag_relax, atoms_relax) = relaxer.relax(atoms, steps=steps)
        return (flag_relax, atoms_relax, filename)

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
    (flag_relax, atoms_relax, filename) = relax_struct(BaseName=BaseName)
    if not flag_relax:
        print('WARNING: Relax failed.')
    energy = atoms_relax.get_potential_energy()
    stress = atoms_relax.get_stress(voigt=False)
    stress = np.array(stress) / GPa
    np.savetxt('ENERGY', [energy])
    np.savetxt('STRESS', stress)
    struct = AseAtomsAdaptor.get_structure(atoms_relax)
    flnm_lst = re.split(r'\.', filename)
    if len(flnm_lst) == 1:
        #POSCAR
        flnm_lst.append('relax')
    else:
        flnm_lst.insert(-1, 'relax')
    filename = '.'.join(flnm_lst)
    struct.to(filename=filename)