#!python
#
#This is the new interface to VASP (by pymatgen)
#In this version, only the POSCAR is required.
#(with correct pymatgen configuration)
#
from pymatgen.io.vasp.inputs import Poscar, Incar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.core.structure import Structure
import elastic3rd.energy.glue as glue
import numpy as np

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
    RunStr = "mpirun -np {%d} vasp_std".format(NP)
    return RunStr

def get_base_vec(BaseName=None):
    #Get the base vector
    structure = Structure.from_file('POSCAR')
    basevec = structure.lattice.matrix
    return basevec

def write_base_vec(BaseName=None, BaseVec=None, F=None):
    structure = Structure.from_file('POSCAR')
    structure.lattice = BaseVec
    structure.to(fmt='poscar', filename='POSCAR')

def add_perturbation(BaseName=None, tor=2.0E-10, filename='POSCAR', PertMethod='random', PertType='relative'):
    pass

def get_energy(BaseName=None):
    outcar = Outcar(filename='OUTCAR')
    energy = outcar.final_energy
    Energy = glue.multi_energy(energy)
    return Energy

def get_stress(BaseName=None):
    vasprun = Vasprun(filename='vasprun.xml')
    ionic_steps = vasprun.ionic_steps
    print(ionic_steps[-1]['stress'])

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
    os.rename(os.path.join(dstpath, "POSCAR"), os.path.join(dstpath, "POSCAR-org"))
    shutil.copyfile("CONTCAR", dstpath + "POSCAR")

def rem_files(BaseName=None):
    pass

'''
get_base_vec()

basevec = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
write_base_vec(BaseName=None, BaseVec=basevec)

get_stress(BaseName=None)
'''