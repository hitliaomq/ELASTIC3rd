#!python
#for chgnet: https://github.com/CederGroupHub/chgnet
import re, os, shutil
import numpy as np
import elastic3rd.energy.glue as glue
import elastic3rd.esutils as esutils
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from chgnet.model.model import CHGNet
from chgnet.model import StructOptimizer
#from mattersim.applications.relax import Relaxer
#from mattersim.forcefield.potential import MatterSimCalculator
from ase.units import GPa
from monty.serialization import loadfn, dumpfn

def get_struct_from_file(BaseName=None):
    file_exts = ['cif', 'poscar', 'json', 'mcsqs', 'yaml']
    filenames = ['{}.{}'.format(BaseName, ext) for ext in file_exts]
    filenames.append('POSCAR')
    for filenamei in filenames:
        if os.path.exists(filenamei):
            struct = Structure.from_file(filenamei)
            filename = filenamei
    return (struct, filename)

def get_base_vec(BaseName=None):
    '''
    Get the base vector from POSCAR
    Parameter
    --------
        BaseName: str
            Optional, not used. Just for consist with castep
    Return
        BaseVec: np.ndarray
            The crystal vector
    '''
    (struct, filename) = get_struct_from_file(BaseName=BaseName)
    BaseVec = struct.lattice.matrix
    return BaseVec

def write_base_vec(BaseName=None, BaseVec=None, F=None):
    '''
    Write Base vector to POSCAR
    Parameter
    --------
        BaseName: str
            Any string, just for consist with castep
        BaseVec: np.ndarray
            The crystal vector
    Return
        None
    '''
    (struct, filename) = get_struct_from_file(BaseName=BaseName)
    struct.lattice = BaseVec
    struct.to(filename=filename)
    
def static(struct):
    chgnet = CHGNet.load()
    prediction = chgnet.predict_structure(struct)
    energy = prediction['energy'][0]    #eV/atom
    stress = prediction['stress'][0]   #GPa
    return (energy, stress, struct)

def relax_struct(BaseName):
    (struct, filename) = get_struct_from_file(BaseName=BaseName)
    if os.path.exists('settings.yaml'):
        settings = loadfn('settings.yaml')
    else:
        settings = {}
    model = settings.get('model', None)
    optimizer = settings.get('optimizer', 'FIRE')
    stress_weight = settings.get('stress_weight', 1.0)
    fmax = settings.get('fmax', 0.1)
    steps = settings.get('steps', 500)
    ase_filter = settings.get('ase_filter', 'FrechetCellFilter')
    verbose = settings.get('verbose', True)
    flageg = settings.get('flageg', 'E')
    if flageg == 'E':
        relax_cell = False
    elif flageg == 'G':
        relax_cell = True
    else:
        relax_cell = False
    if flageg == 'EO':
        (energy, stress, struct) = static(struct, model_name='0.3.0')
    else:  
        relaxer = StructOptimizer(model=model, optimizer_class=optimizer, stress_weight=stress_weight)
        result = relaxer.relax(struct, fmax=fmax, steps=steps, relax_cell=relax_cell,
            ase_filter=ase_filter, verbose=verbose)
        energy = result['trajectory'].energies[-1]
        stress = result['trajectory'].stresses[-1]
        struct = result['final_structure']
    stress = voigt2mat_stress(stress)
    return (energy, stress, struct, filename)

def voigt2mat_stress(stress_v):
    s1, s2, s3, s4, s5, s6 = stress_v
    stress_mat = np.array([[s1, s6, s5], 
        [s6, s2, s4], 
        [s5, s4, s3]])
    return stress_mat

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
    np.savetxt('ENERGY', [energy])
    np.savetxt('STRESS', stress)
    flnm_lst = re.split(r'\.', filename)
    if len(flnm_lst) == 1:
        #POSCAR
        flnm_lst.append('relax')
    else:
        flnm_lst.insert(-1, 'relax')
    filename = '.'.join(flnm_lst)
    struct.to(filename=filename)

def get_energy(BaseName=None):
    '''
    Get the energy from OSZICAR
    Parameter
    ---------
        BaseName: str
            Optional, any string, just for consist with castep
    Return
    ------
        energy: float
            energy in multi unit
    '''
    energy = np.loadtxt('ENERGY')
    # eV, Hartree, kJ/mol, kcal/mol
    Energy = glue.multi_energy(energy)
    return Energy

def get_stress(BaseName=None):
    '''
    Get the stress from OUTCAR
    Parameter
    ---------
        BaseName: str
            Optional, any string, just for consist with castep
    Return
        Stress: np.ndarray
            Stress matrix
    '''
    Stress = np.loadtxt('STRESS')
    return Stress

def copy_files(BaseName=None, Path=None):
    '''
    Copy required files(INCAR, POSCAR, POTCAR, KPOINTS) for run VASP
    Parameters
    ----------
        BaseName: str
            Any string, just for consist with castep
        Path: path-like str
            The destination folder to store the files
    Return
    ------
        None
    '''
    #TODO: need to be improved
    (struct, filename) = get_struct_from_file(BaseName=BaseName)
    shutil.copyfile(filename, os.path.join(Path, filename))
    if os.path.exists('settings.yaml'):
        shutil.copyfile('settings.yaml', os.path.join(Path, 'settings.yaml'))

def update_structure(BaseName=None, dstpath="./"):
    '''
    Update POSCAR for geom opt
    Parameters
    ----------
        BaseName: str
            Any string, just for consist with castep
        dstpath: path-like str
            The destination folder to store the files
    Return
    ------
        None
    '''
    (struct, filename) = get_struct_from_file(BaseName=BaseName)
    flnm_lst = re.split(r'\.', filename)
    flnm_lst_relax = re.split(r'\.', filename)
    if len(flnm_lst) == 1:
        #POSCAR
        flnm_lst.append('org')
        flnm_lst_relax.append('relax')
    else:
        flnm_lst.insert(-1, 'org')
        flnm_lst_relax.insert(-1, 'relax')
    filename_org = '.'.join(flnm_lst)
    filaname_relax = '.'.join(flnm_lst_relax)
    os.rename(os.path.join(dstpath, filename), os.path.join(dstpath, filename_org))
    shutil.copyfile(filaname_relax, os.path.join(dstpath, filename))

def update_input_for_geomopt(BaseName=None, FlagEG='G', KeepSym=False):
    '''
    Update the incar file for GeomOpt
    Parameters
    ----------
        BaseName: str
            Any string, just for consist with castep
        FlagEG: str
            'G' or 'E' for full relax and single-point calculations
            E: only relax the ions
            G: full relaxiation (ions and cell)
            GS: fix shape
            GV: fix volume
            GI: fix ions
            EO: energy only (fix all)
    Return
    ------
        None
    '''
    if os.path.exists('settings.yaml'):
        settings = loadfn('settings.yaml')
    else:
        settings = {}
    settings.update({'flageg': FlagEG, 'constrain_symmetry': KeepSym})
    dumpfn(settings, 'settings.yaml')

def add_perturbation(BaseName, tor, PertMethod):
    #Not implemented in mattersim
    pass

def rem_files(BaseName):
    #For mattersim, nothing to be removed
    pass

#..............TEST...................#
