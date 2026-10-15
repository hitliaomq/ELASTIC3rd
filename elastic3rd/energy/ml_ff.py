#!python
#for mattersim: https://github.com/microsoft/mattersim
import re, os, shutil, warnings
import numpy as np
from typing import Iterable, List, Tuple, Union

import elastic3rd.energy.glue as glue
import elastic3rd.esutils as esutils
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase import Atoms
from ase.constraints import FixSymmetry
from ase.filters import Filter
from ase.filters import ExpCellFilter, FrechetCellFilter, StrainFilter, UnitCellFilter
from ase.optimize import BFGS, FIRE, FIRE2, LBFGS
from ase.optimize.optimize import Optimizer
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

class Relaxer(object):
    """Relaxer is a class for structural relaxation with fixed volume."""

    SUPPORTED_OPTIMIZERS = {"BFGS": BFGS, "FIRE": FIRE, "LBFGS": LBFGS, "FIRE2": FIRE2}
    SUPPORTED_FILTERS = {
        "EXPCELLFILTER": ExpCellFilter,
        "FRECHETCELLFILTER": FrechetCellFilter,
        "STRAINFILTER": StrainFilter,
        "UNITCELLFILTER": UnitCellFilter
    }

    def __init__(
        self,
        optimizer: Union[Optimizer, str] = "FIRE",
        filter: Union[Filter, str, None] = None,
        constrain_symmetry: bool = True,
        fix_axis: Union[bool, Iterable[bool]] = False,
    ) -> None:
        """
        Args:
            optimizer (Union[Optimizer, str]): The optimizer to use.
            filter (Union[Filter, str, None]): The filter to use.
            constrain_symmetry (bool): Whether to constrain the symmetry.
            fix_axis (Union[bool, Iterable[bool]]): Whether to fix the axis.
        """
        self.optimizer = (
            self.SUPPORTED_OPTIMIZERS[optimizer.upper()]
            if isinstance(optimizer, str)
            else optimizer
        )
        self.relax_cell = filter is not None
        if filter is not None:
            self.filter = (
                self.SUPPORTED_FILTERS[filter.upper()]
                if isinstance(filter, str)
                else filter
            )
        self.constrain_symmetry = constrain_symmetry
        self.fix_axis = fix_axis

    def relax(
        self,
        atoms: Atoms,
        steps: int = 500,
        fmax: float = 0.01,
        params_filter: dict = {},
        **kwargs,
    ) -> Tuple[bool, Atoms]:
        """
        Relax the atoms object.

        Args:
            atoms (Atoms): The atoms object to relax.
            steps (int): The maximum number of steps to take.
            fmax (float): The maximum force allowed.
            params_filter (dict): The parameters for the filter.
            kwargs: Additional keyword arguments for the optimizer.
        """

        if atoms.calc is None:
            raise ValueError("Atoms object must have a calculator.")

        if self.constrain_symmetry:
            atoms.set_constraint(FixSymmetry(atoms))

        if self.relax_cell:
            # Set the mask for the fixed axis
            if isinstance(self.fix_axis, bool):
                mask = [not self.fix_axis for i in range(6)]
            else:
                assert (
                    len(self.fix_axis) == 6
                ), "The length of fix_axis list not equal 6."
                mask = [not elem for elem in self.fix_axis]

            # check if the scalar_pressure is provided
            if (
                "scalar_pressure" in params_filter
                and params_filter["scalar_pressure"] > 1
            ):
                warnings.warn(
                    "The scalar_pressure used in ExpCellFilter assumes "
                    "eV/A^3 unit and 1 eV/A^3 is already 160 GPa. "
                    "Please make sure you have converted your pressure "
                    "from GPa to eV/A^3 by dividing by 160.21766208."
                )
            ecf = self.filter(atoms, mask=mask, **params_filter)
        else:
            ecf = atoms
        optimizer = self.optimizer(ecf, **kwargs)
        optimizer.run(fmax=fmax, steps=steps)

        converged = optimizer.get_number_of_steps() < steps

        if self.constrain_symmetry:
            atoms.set_constraint(None)

        return converged, atoms

    @classmethod
    def relax_structures(
        cls,
        atoms: Union[Atoms, Iterable[Atoms]],
        optimizer: Union[Optimizer, str] = "FIRE",
        filter: Union[Filter, str, None] = None,
        constrain_symmetry: bool = False,
        fix_axis: Union[bool, Iterable[bool]] = False,
        pressure_in_GPa: Union[float, None] = None,
        **kwargs,
    ) -> Union[Tuple[bool, Atoms], Tuple[List[bool], List[Atoms]]]:
        """
        Args:
            atoms: (Union[Atoms, Iterable[Atoms]]):
                The Atoms object or an iterable of Atoms objects to relax.
            optimizer (Union[Optimizer, str]): The optimizer to use.
            filter (Union[Filter, str, None]): The filter to use.
            constrain_symmetry (bool): Whether to constrain the symmetry.
            fix_axis (Union[bool, Iterable[bool]]): Whether to fix the axis.
            **kwargs: Additional keyword arguments for the relax method.
        Returns:
            converged (Union[bool, List[bool]]):
                Whether the relaxation converged or a list of them
            Atoms (Union[Atoms, List[Atoms]]):
                The relaxed atoms object or a list of them
        """
        params_filter = {}

        if filter is None and pressure_in_GPa is None:
            pass
        elif filter is None and pressure_in_GPa is not None:
            filter = "ExpCellFilter"
            params_filter["scalar_pressure"] = (
                pressure_in_GPa * GPa
            )  # GPa = 1 / 160.21766208
        elif filter is not None and pressure_in_GPa is None:
            params_filter["scalar_pressure"] = 0.0
        else:
            params_filter["scalar_pressure"] = (
                pressure_in_GPa * GPa
            )  # GPa = / 160.21766208

        relaxer = Relaxer(
            optimizer=optimizer,
            filter=filter,
            constrain_symmetry=constrain_symmetry,
            fix_axis=fix_axis,
        )

        if isinstance(atoms, (list, tuple)):
            relaxed_results = relaxed_results = [
                relaxer.relax(atom, params_filter=params_filter, **kwargs)
                for atom in atoms
            ]
            converged, relaxed_atoms = zip(*relaxed_results)
            return list(converged), list(relaxed_atoms)
        else:
            return relaxer.relax(atoms, params_filter=params_filter, **kwargs)

