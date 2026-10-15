#!python
#for mace: https://github.com/ACEsuit/mace
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

#from mace.calculators import mace_mp
from matris.applications import MatRISCalculator, StructOptimizer
    
def relax_struct(BaseName=None):
    #device = "cuda" if torch.cuda.is_available() else "cpu"
    if os.path.exists('settings.yaml'):
        settings = loadfn('settings.yaml')
    else:
        settings = {}
    optimizer = settings.get('optimizer', 'BFGS')
    device = settings.get('device', 'cpu')
    fmax = settings.get('fmax', 0.05)
    task = settings.get('task', 'efsm')
    ase_filter = settings.get('ase_filter', 'ExpCellFilter')
    constrain_symmetry = settings.get('constrain_symmetry', True)
    model = settings.get('model', 'matris_10m_oam')
    dtype = settings.get('dtype', 'float32')
    level = settings.get('level', 0)
    steps = settings.get('seps', 800)
    flageg = settings.get('flageg', 'E')
    (struct, filename) = get_struct_from_file(BaseName=BaseName)
    atoms = AseAtomsAdaptor.get_atoms(struct)
    if flageg == 'E':
        fix_axis = True
    elif flageg == 'G':
        fix_axis = False
    elif flageg == 'C':
        fix_axis = False
        ase_filter = 'StrainFilter'
    else:
        fix_axis = False
    
    if flageg == 'EO':
        matris_calc = MatRISCalculator(
            model = model,
            task = task,
            device = device
        )
        matris_calc.calculate(
            atoms = atoms,
            properties = None,
            system_changes = None
        )
        results = matris_calc.results
        energy = results['energy']
        stress = results['stress']
        final_structure = atoms
        return (energy, stress, final_structure, filename)
    else:
        matris_opt = StructOptimizer(
            model = model, 
            task = task,
            optimizer = optimizer, # FIRE, BFGS ...
            device=device
        )
        opt_result = matris_opt.relax(
            atoms = atoms,
            verbose = True,
            steps = steps,
            fmax = fmax,
            relax_cell = not fix_axis,
            ase_filter = ase_filter
        )
        trajectory = opt_result['trajectory']
        energy = trajectory.energies[-1]
        stress = trajectory.stresses[-1]
        final_structure = opt_result['final_structure']

    return (energy, stress, final_structure, filename)

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
    (energy, stress, struct, filename) = relax_struct(BaseName=BaseName)
    stress = np.array(stress)
    if len(stress.shape) == 1:
        (s1, s2, s3, s4, s5, s6) = stress
        stress = np.array([[s1, s6, s5], [s6, s2, s4], [s5, s4, s3]])    
    stress = np.array(stress) / GPa
    np.savetxt('ENERGY', [energy])
    np.savetxt('STRESS', stress)
    #struct = AseAtomsAdaptor.get_structure(final_structure)
    flnm_lst = re.split(r'\.', filename)
    if len(flnm_lst) == 1:
        #POSCAR
        flnm_lst.append('relax')
    else:
        flnm_lst.insert(-1, 'relax')
    filename = '.'.join(flnm_lst)
    #struct.to(filename=filename)