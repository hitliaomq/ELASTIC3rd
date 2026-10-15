#!python
import re, os, copy, shutil
import linecache
import numpy as np
from monty.serialization import loadfn
import elastic3rd.energy.glue as glue
import elastic3rd.esutils as esutils
from pymatgen.core.units import ang_to_bohr, bohr_to_ang
from abacuspy import Input, Stru

def get_base_vec(BaseName, FileName='STRU'):
    struct = Stru.from_file(filename=FileName)
    basevec = bohr_to_ang * np.array(struct.latvec)
    return basevec

def write_base_vec(BaseName=None, BaseVec=None, F=None):
    basevec = np.array(BaseVec) * ang_to_bohr
    struct = Stru.from_file('STRU')
    struct.latvec = basevec.tolist()
    struct.to_file(filename='STRU')

def run(NP, BaseName=None):
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
    RunStr = "mpirun -np " + str(int(NP)) + " abacus"
    return RunStr

def get_energy(BaseName=None):
    FileName = 'abacus.json'
    abacus_out = loadfn(FileName)['output']
    energy = abacus_out[-1]['energy']
    Energy = glue.multi_energy(energy)
    return Energy

def get_stress(BaseName=None):
    FileName = 'abacus.json'
    abacus_out = loadfn(FileName)['output']
    stress = abacus_out[-1]['stress']
    stress = 0.1 * np.array(stress)
    return -stress

def copy_files(BaseName=None, Path=None):
    '''
    Copy required files(STRU INPUT KPT) for run VASP
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
    shutil.copyfile("STRU", Path + "/" + "STRU")
    shutil.copyfile("KPT", Path + "/" + "KPT")
    shutil.copyfile("INPUT", Path + "/" + "INPUT")

def update_structure(BaseName=None, dstpath="./"):
    '''
    Update POSCAR according to CONTCAR
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
    os.rename(os.path.join(dstpath, "STRU"), os.path.join(dstpath, "STRU-org"))
    shutil.copyfile(os.path.join('OUT.{}'.format(BaseName), "STRU_ION_D"), dstpath + "STRU")

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
    base_setting = {'calculation': 'cell-relax', 'relax_nmax': '100', 'cal_stress': '1'}
    if FlagEG == 'G':
        input_setting = copy.deepcopy(base_setting)
    elif FlagEG == 'E':
        input_setting = copy.deepcopy(base_setting)
        input_setting.update({'calculation': 'relax'})
    elif FlagEG == 'GS':
        input_setting = {'fixed_axes': 'shape', 'fixed_atoms': True}
        input_setting.update(base_setting)
    elif FlagEG == 'GV':
        input_setting = {'fixed_axes': 'volume'}
        input_setting.update(base_setting)
    elif FlagEG == 'GI':
        input_setting = {'fixed_atoms': True}
        input_setting.update(base_setting)
    elif FlagEG == 'EO':
        input_setting = {'calculation': 'scf', 'cal_stress': '1', 'relax_nmax': 1}
    input_file = Input.from_file('INPUT')
    input_file.update(input_setting)
    input_file.to_file(filename='INPUT')

def add_perturbation(BaseName=None, tor=2.0E-10, filename='STRU', PertMethod='random', PertType='relative'):
    #Currently not implemented in abacus
    pass

def rem_files(BaseName):
    #For mattersim, nothing to be removed
    pass