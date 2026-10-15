#!python
#
import re
import os
import copy
import shutil
import platform
import warnings
import numpy as np
import elastic3rd.energy.glue as glue
import elastic3rd.esutils as esutils
from scipy.constants import physical_constants
#from elastic3rd.tools.IdealStrength_stress import get_smat_by_eig, , update_input_for_is, getPrinciple
import elastic3rd.tools.IdealStrength_stress as iss

def get_base_vec(BaseName):
    '''
    Get the base vector from cell file
    Parameter
    --------
        BaseName: str
            The filename of the cell file (without ext)
    Return
        BaseVec: np.ndarray
            The crystal vector
    '''
    FileName = BaseName + ".cell"
    cell_param = parse_cell(FileName)
    basevec = cell_param['LATTICE_CART']['value']
    return np.array(basevec, dtype=np.float64)

def get_base_vec_by_geom_file(BaseName):
    '''
    Get the base vector from geom file
    Parameter
    --------
        BaseName: str
            The filename of the geom file (without ext)
    Return
        BaseVec: np.ndarray
            The crystal vector
    '''
    geom_file = BaseName + ".geom"
    geom_param = parse_geom(geom_file)
    return np.array(geom_param[-1]['Lattice'], dtype=np.float64)

def write_base_vec(BaseName, BaseVec, F=None):
    '''
    Write Base vector to cell file(update cell file by base vector)
    Parameter
    --------
        BaseName: str
            The filename of the cell file (without ext)
        BaseVec: np.ndarray
            The crystal vector
    Return
        None
    '''
    FileName = BaseName + ".cell"
    cell_param = parse_cell(FileName)
    BaseVec = [[format(float(baseveci), '.15f') for baseveci in BaseVeci] for BaseVeci in BaseVec]
    cell_param['LATTICE_CART'] = {'type': 'BLOCK', 'value': BaseVec}
    write_cell(FileName, cell_param)

def run(NP, BaseName):
    '''
    Run string for CASTEP. 
        Note: this string is for Materials Studio's version
    Parameters
    ----------
        BaseName: str
            The filename of the cell and param file (without ext)
        NP: int
            NP is the total cores
    Return
    ------
        RunStr: str
            The string for calling first principles code
    '''
    RunStr = "mpirun -np " + str(int(NP)) + " castep " + BaseName
    return RunStr

def copy_files(BaseName, Path):
    '''
    Copy required files(cell, param, RunCASTEP.sh/bat file) for run CASTEP
    Parameters
    ----------
        BaseName: str
            The filename of the cell and param file (without ext)
        Path: path-like str
            The destination folder to store the files
    Return
    ------
        None
    '''
    #TODO: need to be improved
    plat = platform.platform().split("-")[0].lower()
    shutil.copyfile(BaseName + ".cell", Path + "/" + BaseName + ".cell")
    shutil.copyfile(BaseName + ".param", Path + "/" + BaseName + ".param")

#TODO: add a parameter to control the remove
def rem_files(BaseName):
    '''
    Remove large files, e.g. check file
    Parameter
    --------
        BaseName: str
            The filename of the cell and param file (without ext)
    Return
    ------
        None
    '''
    #os.remove(path)
    os.system("rm *.check *.pid *.castep_bin *.bands")

def get_energy(BaseName):
    '''
    Get the energy from geom file
    Parameter
    ---------
        BaseName: str
            The filename of the cell and param file (without ext)
    Return
    ------
        energy: float
            energy in multi unit
    '''
    FileName = BaseName + '.geom'
    if os.path.exists(FileName):
        geom_param = parse_geom(FileName)
        Energy = geom_param[-1]['Energy'][0][0]
    else:
        Energy = get_energy_from_castep(BaseName)
    return glue.multi_energy(float(Energy))

def get_stress(BaseName, tovec=False):
    '''
    Get the stress from geom file, if the STRESS is not in geom file, then try to get the stress from castep file
    Parameter
    ---------
        BaseName: str
            The filename of the cell and param file (without ext)
    Return
        Stress: np.ndarray
            Stress matrix
    '''
    FileName = BaseName + '.geom'
    FlagFromCASTEP = False
    if not os.path.exists(FileName):
        FlagFromCASTEP = True
    else:
        geom_param = parse_geom(FileName)
        Stress = np.array(geom_param[-1]['Stress'])
        if Stress.size == 0:
            FlagFromCASTEP = True
    if FlagFromCASTEP:
        Stress = get_stress_from_castep(BaseName)
    if tovec:
        Stress = [Stress[0][0], Stress[1][1], Stress[2][2], Stress[1][2], Stress[0][2], Stress[0][1]]
    return Stress

def get_stress_from_castep(BaseName):
    '''
    Get stress from castep file
    Parameter
    ---------
        BaseName: str
            The filename of the cell and param file (without ext)
    Return
        Stress: np.ndarray
            Stress matrix
    '''
    Stress = np.zeros((3, 3))
    FileName = BaseName + ".castep"
    flag = 0
    count = 1
    i = 0
    with open(FileName, 'r') as fopen:
        for eachline in fopen:
            eachline = eachline.strip();
            eachlinei = eachline.split("*")
            if ' Stress Tensor ' in eachlinei:
                flag = 1
            if flag :
                count = count + 1
            if count < 8:
                continue
            else:
                Stressi = re.split(r'\s+', eachlinei[1])
                Stress[i, :] = Stressi[2:5]
                i = i + 1
                if count > 9:
                    flag = 0
                    count = 1
                    i = 0
    if Stress.shape[0] == 0:
        print("Please set calculate_stress as ture.")
    return Stress

def get_strain(BaseName=None, basevec0=None, tovec=False):
    #This function used after geom optimization, hence get from geom file
    if basevec0 is None:
        strain = np.zeros((3, 3))
    else:
        new_basevec = get_base_vec_by_geom_file(BaseName)
        strain = np.dot(np.linalg.inv(basevec0), new_basevec) - np.eye(3)
    if tovec:
        strain = [strain[0, 0], strain[1, 1], strain[2, 2], strain[1, 2], strain[0, 2], strain[0, 1]]
    return strain

def get_pos(BaseName):
    '''
    Get atom's position from geom file
    Parameter
    ---------
        BaseName: str
            The filename of the cell and param file (without ext)
    Return
        pos: list
            list of atomic position
        atoms: list
            list of atoms
    '''
    FileName = BaseName + '.geom'
    geom_param = parse_geom(FileName)
    pos = geom_param[-1]['Position']
    atoms = geom_param[-1]['Atom']
    return (pos, atoms)

def get_pos_vca():
    pass

def normalize_composition(x):
    '''
    Normalize composition
    Parameter
    ---------
        x: dict
            e.g. x = {'Nb': 1, 'Ti': 1, 'V': 1, 'Zr': 1}
    '''
    totalx = sum([x[ele] for ele in x])
    x = {ele: x[ele]/totalx for ele in x}
    return x

def get_com_vca(BaseName):
    celldict = parse_cell('{}.cell'.format(BaseName))['POSITIONS_FRAC']['value']
    eles = celldict['ELEMENTS']
    mixs = celldict['MIXTURE']
    if 'ISMIX' in celldict and celldict['ISMIX']:
        com = [1 if len(mix) == 0 else mix[1] for mix in mixs]
    else:
        com = [1 for ele in eles]

    composition = {}
    for i, ele in enumerate(eles):
        if ele in composition:
            composition[ele] = composition[ele] + float(com[i])
        else:
            composition[ele] = float(com[i])

    return normalize_composition(composition)

def update_pos(BaseName, NewPos, atoms=None):
    '''
    Update the cell file according to atom's position, 
        if atoms not provided, try to parse from cell file in current folder
        TODO: add support for VCA's cell file
    Parameter
    --------
        BaseName: str
            The filename of the cell and param file (without ext)
        NewPos: list
            The atomic position
        atoms: list
            The list of atom
    '''
    cell_param = parse_cell(BaseName + ".cell")
    pos_Settings = cell_param['POSITIONS_FRAC']
    if not atoms:
        cell_param = parse_cell(BaseName + ".cell")
        atoms = pos_Settings['value']['ELEMENTS']
    pos = []
    for i, atomi in enumerate(atoms):
        pos.append([format(float(item), '.15f') for item in NewPos[i]])
    pos_Settings['value']['POSITIONS'] = pos
    Settings = {'POSITIONS_FRAC': pos_Settings}
    check_cell(BaseName, Settings)

def add_perturbation(BaseName, tor=2.0E-10, PertMethod='random'):
    '''
    Add a small perturbation in the position in the cell file

    Parameter
    ---------
        BaseName: str
            The name of the *cell file
        tor: float
            The amplitude of the pertubation, usually a small value
        PertMethod: str
            The method for perturbation, random/uniform, only the first letter works
    Return
    ------
    '''
    cell_param = parse_cell(BaseName + ".cell")
    pos_and_atom = cell_param['POSITIONS_FRAC']['value']
    atoms = pos_and_atom['ELEMENTS']
    pos = pos_and_atom['POSITIONS']
    #atoms = [atom[0] for atom in pos_and_atom]
    #pos = np.array([posi[1:] for posi in pos_and_atom], dtype=np.float64)
    pert = esutils.gen_perturbation(pos, tor=tor, PertMethod=PertMethod)
    new_pos = np.array(pos).astype('float64') + np.array(pert).astype('float64')
    update_pos(BaseName, new_pos, atoms=atoms)    

def parse_param(param_file):
    '''
    Parse the param file, return the settings(dict)
    Rules: http://www.tcm.phy.cam.ac.uk/castep/documentation/WebHelp/content/modules/castep/expcastepfileparam.htm
    Parameter
    --------
        param_file: str
            The filename of param file (with ext)
    Return
    ------
        parameters: dict
            The parameters dict
    '''
    parameters = {}
    with open(param_file, 'r') as fopen:
        for eachline in fopen:
            linei = remove_comments(eachline)
            if linei:
                #Not empty
                if len(linei) > 1 :
                    key = linei[0].lower()  #case insensitive
                    val = linei[1]
                    parameters[key] = val
    return parameters

def write_param(param_file, param_dict):
    '''
    Write out param file according to param_dict
    Parameter
    --------
        param_file: str
            The filename of param file (with ext)
        param_dict: dict
            The parameters dict
    Return
    ------
        None
    '''
    with open(param_file, 'w+') as fid:
        for key in param_dict:
            if not ((param_dict[key] is None) or (param_dict[key] == '')):
                fid.write('{} : {}\n'.format(key, param_dict[key]))

def check_param(BaseName, Settings={'task': 'GeometryOptimization', 'calculate_stress': 'true'}):
    #Check the param(update param file according to Settings)
    #Checking task==GeometryOptimization, calculate_stress==true, cut_off_energy=360
    '''
    Update param file by Settings
    Parameter
    --------
        BaseName: str
            The filename of the cell and param file (without ext)
        Settings: dict
            The parameters dict
    Return
    ------
        None
    '''
    FileName = BaseName + '.param'
    parameters = parse_param(FileName)
    empty_key = []
    for key in Settings:
        value = Settings[key]
        if (value is None) or (value == ''):
            empty_key.append(key)            
    for key in empty_key:
        Settings.pop(key)
        if key in parameters:
                parameters.pop(key)
    parameters.update(Settings)
    write_param(FileName, parameters)

def update_input_for_geomopt(BaseName, FlagEG='G', KeepSym=False):
    '''
    Update cell and param file for geom controled by FlagEG
    Param
    -----
        BaseName: str
            The base name of cell and param file
        FlagEG: str
            E: only relax the ions
            G: full relaxiation (ions and cell)
            GS: fix shape
            GV: fix volume
            GI: fix ions
            EO: energy only (fix all)
    '''
    GEOM_PARAM = ['GEOM_CONVERGENCE_WIN', 'GEOM_DISP_TOL', 'GEOM_ENERGY_TOL', 'GEOM_FORCE_TOL', 'GEOM_FREQUENCY_EST', 'GEOM_MAX_ITER',
                  'geom_method', 'GEOM_MODULUS_EST', 'GEOM_PRECONDITIONER', 'GEOM_SPIN_FIX', 'GEOM_STRESS_TOL']
    param_settings = {'task' : 'GeometryOptimization', 'calculate_stress' : 'true'}
    FlagEG = FlagEG.upper()
    if FlagEG == 'G':
        cell_settings = {'FIX_ALL_CELL' : {}, 'CELL_CONSTRAINTS': {'type': 'BLOCK', 'value': [['1', '2', '3'],['4', '5', '6']]}}
    elif FlagEG == 'E':
        cell_settings = {'FIX_ALL_CELL' : {'type': 'SINGLE', 'value': 'true'}, 'CELL_CONSTRAINTS': {}}
        param_settings.update({'geom_method': 'Delocalized'})
    elif FlagEG == 'GS':
        #Fix shape
        cell_settings = {'FIX_ALL_CELL' : {}, 'CELL_CONSTRAINTS': {'type': 'BLOCK', 'value': [['1', '1', '1'],['0', '0', '0']]}}
    elif FlagEG == 'GV':
        #Fix volume
        cell_settings = {'FIX_ALL_CELL' : {}, 'CELL_CONSTRAINTS': {}, 'FIX_VOL': {'type': 'SINGLE', 'value': 'true'}}
    elif FlagEG == 'GI':
        cell_settings = {'FIX_ALL_CELL' : {}, 'CELL_CONSTRAINTS': {}, 'FIX_ALL_IONS': {'type': 'SINGLE', 'value': 'true'}}
    elif FlagEG == 'EO':
        param_settings = {'task' : 'SinglePoint', 'calculate_stress' : 'true'}
        for key in GEOM_PARAM:
            param_settings[key.lower()] = None
        cell_settings = {'FIX_ALL_CELL' : {'type': 'SINGLE', 'value': 'true'}, 'CELL_CONSTRAINTS': {}}

    if KeepSym:
        cell_settings['SYMMETRY_GENERATE'] = {'type': 'SINGLE', 'value': ['']}
    else:
        cell_settings['SYMMETRY_GENERATE'] = {}

    check_param(BaseName, param_settings)
    check_cell(BaseName, cell_settings)

def remove_comments(line_str):
    '''
    Remove comments(#, ;, !, COMMENT, comment)
    If it's comments line or empty line, return None
    If there are some comment behind the line, remove it, return a list

    Parameter
    ---------
        line_str: str
            A single line string
    Return
    ------
        line_list: list
            list that generated by split with : and space
    '''
    line_str = line_str.strip()
    if line_str.startswith(('#', ';', '!', 'COMMENT', 'comment')):
        return None
    else:
        if line_str:
            line_list = re.split(r'(?:#|!|;)', line_str)
            #The string after # or ! or ; are comments. So only the first one is the content need to be kept
            line_list = line_list[0]
            #split with : or space or = or ( or )
            line_list = re.split(r'(?:\s*:\s*|\s*=\s*|\s*\(\s*|\s*\)\s*|\s+)', line_list)
            #remove empty elements
            line_list = [line_list_i for line_list_i in line_list if line_list_i]
            return line_list
        else:
            #Empty
            return None

#Not used now, to be deleted
def remove_comments_back(line_str):
    '''
    Remove comments(#, ;, !, COMMENT, comment)
    If it's comments line or empty line, return None
    If there are some comment behind the line, remove it, return a list

    Parameter
    ---------
        line_str: str
            A single line string
    Return
    ------
        line_list: list
            list that generated by split with : and space
    '''
    line_str = line_str.strip()
    if line_str.startswith(('#', ';', '!', 'COMMENT', 'comment')):
        return None
    else:
        if line_str:
            #split with : or space or = or ( or )
            line_list = re.split(r'(?:\s*:\s*|\s*=\s*|\s*\(\s*|\s*\)\s*|\s+)', line_str)
            #remove empty elements
            line_list = [line_list_i for line_list_i in line_list if line_list_i]
            for i, linei in enumerate(line_list):
                if linei.startswith(('#', ';', '!')):
                    line_list = line_list[0:i]
                    break
                else:
                    linei_j = re.split(r'(?:#|!|;)', linei)
                    if len(linei_j) > 1:
                        line_list[i] = linei_j[0]
                        line_list = line_list[0:i+1]
                        break
            return line_list
        else:
            #Empty
            return None

def parse_cell(cell_file):
    '''
    Parse the cell file, return the parameters (dict)
    Rules: http://www.tcm.phy.cam.ac.uk/castep/documentation/WebHelp/content/modules/castep/expcastepfilecell.htm

    Parameter
    ---------
        cell_file: str (path-like)
            The filename of cell file
    Return
    ------
        cell_parameters: dict
            cell_parameters = {key1: {'type': 'SINGLE'/'BLOCK', 'value': keyValue}, key2: {}}
            SINGLE for the value and key in the same line
            BLOCK for the value was located between %BLOCK KEY and %ENDBLOCK KEY
            the keyValue should be a 1-d list (for SINGLE) or 2-d list (for BLOCK)
    '''
    cell_parameters = {}
    flag_block = 0
    with open(cell_file, 'r') as fopen:
        for eachline in fopen:
            linei = remove_comments(eachline)
            if linei:
                #Not empty
                tmpkey = linei[0].upper()
                if tmpkey == '%BLOCK':
                    key = linei[1].upper()
                    val = {'type': 'BLOCK', 'value': []}
                    flag_block = 1
                elif tmpkey == '%ENDBLOCK':
                    flag_block = 0
                    cell_parameters[key] = val
                else:
                    if flag_block:
                        if key == 'SPECIES_POT':
                            val['value'].append(re.split(r'\s+', eachline.strip()))
                        else:
                            val['value'].append(linei)   #Note: linei is a list
                    else:
                        val = {'type': 'SINGLE'}
                        n_linei = len(linei)
                        key = tmpkey
                        if n_linei == 1:
                            val['value'] = ['']
                        else:
                            val['value'] = linei[1:]
                        cell_parameters[key] = val
    cell_parameters['POSITIONS_FRAC']['value'] = split_pos_vca(cell_parameters['POSITIONS_FRAC']['value'])
    return cell_parameters

def split_pos_vca(vca_pos):
    '''
    Split atomic position for VCA or with SPIN
    example:
        vca_pos = [['Ti', '0.0', '0.0', '0.0', 'MIXTURE', '1', '0.250000'], 
                   ['Ti', '0.5', '0.5', '0.5', 'MIXTURE', '2', '0.250000'], 
                   ['V', '0.0', '0.0', '0.0', 'MIXTURE', '1', '0.250000'], 
                   ['V', '0.5', '0.5', '0.5', 'MIXTURE', '2', '0.250000'], 
                   ['Zr', '0.0', '0.0', '0.0', 'MIXTURE', '1', '0.250000'], 
                   ['Zr', '0.5', '0.5', '0.5', 'MIXTURE', '2', '0.250000'], 
                   ['Nb', '0.0', '0.0', '0.0', 'MIXTURE', '1', '0.250000'], 
                   ['Nb', '0.5', '0.5', '0.5', 'MIXTURE', '2', '0.250000']]
    return:
        new_pos = {'ELEMENTS': ['Ti', 'Ti', 'V', 'V', 'Zr', 'Zr', 'Nb', 'Nb'],
                   'POSITIONS': [['0.0', '0.0', '0.0'], ['0.5', '0.5', '0.5'], 
                                 ['0.0', '0.0', '0.0'], ['0.5', '0.5', '0.5'],
                                 ['0.0', '0.0', '0.0'], ['0.5', '0.5', '0.5'],
                                 ['0.0', '0.0', '0.0'], ['0.5', '0.5', '0.5']], 
                   'MIXTURE': [['1', '0.25'], ['2', '0.25'], ['1', '0.25'], ['2', '0.25'], 
                               ['1', '0.25'], ['2', '0.25'], ['1', '0.25'], ['2', '0.25']]}
    if NO MIXTURE, then MIXTURE=None
    '''
    new_pos = {'ELEMENTS': [], 'POSITIONS': [], 'MIXTURE': [], 'ISMIX': False, 'SPIN': [], 'MAGMOM': []}
    for posi in vca_pos:
        new_pos['ELEMENTS'].append(posi[0])
        new_pos['POSITIONS'].append([posij for posij in posi[1:4]])
        if 'MIXTURE' in posi:
            new_pos['ISMIX'] = True
            mix_ind = posi.index('MIXTURE')
            new_pos['MIXTURE'].append([posi[mix_ind + 1], posi[mix_ind + 2]])
        else:
            new_pos['MIXTURE'].append([])
        if 'SPIN' in posi:
            spin_index = posi.index('SPIN')
            new_pos['SPIN'].append(posi[spin_index + 1])
        else:
            new_pos['SPIN'].append(0)
        if 'MAGMOM' in posi:
            mag_index = posi.index('MAGMOM')
            new_pos['MAGMOM'].append(posi[mag_index + 1])
        else:
            new_pos['MAGMOM'].append(0)
    return new_pos

def write_cell(cell_file, cell_dict):
    '''
    Write out the cell file according to the cell_dict
    Parameter
    ---------
        cell_file: str
            The filename of cell file (with ext)
        cell_dict: dict
            The cell dict, and the format is the same with parse_cell function
            If the type is not given, it will determine it automatically
    Return
        None
    '''
    species_keys = ['SPECIES_MASS', 'SPECIES_POT', 'SPECIES_LCAO_STATES']
    cell_dict = check_cell_type(cell_dict)
    elements = cell_dict['POSITIONS_FRAC']['value']['ELEMENTS']
    with open(cell_file, 'w+') as fid:
        for key in cell_dict:
            val_type = cell_dict[key]['type']
            val_value = cell_dict[key]['value']
            if val_type == 'BLOCK':
                fid.write('%BLOCK {}\n'.format(key))
                if key == 'POSITIONS_FRAC':
                    for i, elei in enumerate(val_value['ELEMENTS']):
                        valueForWrite = val_value['POSITIONS'][i]
                        valueForWrite = [format(float(valuei), '.15f') for valuei in valueForWrite]
                        linei = "  {}  {}".format(elei, "  ".join(valueForWrite))
                        if val_value['ISMIX']:
                            if val_value['MIXTURE'][i]:
                                #Not empty
                                linei = "{}  MIXTURE:(  {})".format(linei, "  ".join(val_value['MIXTURE'][i]))
                        if 'SPIN' in val_value:
                            if val_value['SPIN'][i]:
                                linei = "{}  SPIN=  {}".format(linei, format(float(val_value['SPIN'][i]), '.10f'))
                        fid.write(linei)
                        fid.write('\n')
                elif key in species_keys:
                    for linei in val_value:
                        if linei[0] in elements:
                            linei = [str(item) for item in linei]
                            if len(linei) > 0:
                                fid.write('    {}\n'.format('    '.join(linei)))
                elif key == 'EXTERNAL_PRESSURE':
                    val_value = [[float(valii) for valii in valp] for valp in val_value]
                    space0 = ' '*16
                    if val_value:
                        fid.write('{:>16.10f}{:>16.10f}{:>16.10f}\n'.format(*val_value[0]))
                        fid.write('{}{:>16.10f}{:>16.10f}\n'.format(space0, *val_value[1]))
                        fid.write('{}{}{:>16.10f}\n'.format(space0, space0, *val_value[2]))
                else:
                    if val_value:
                        for linei in val_value:
                            linei = [str(item) for item in linei]
                            if len(linei) > 0:
                                fid.write('    {}\n'.format('    '.join(linei)))
                fid.write('%ENDBLOCK {}\n\n'.format(key))
            elif val_type == 'SINGLE':
                val_value = [str(item) for item in val_value]
                fid.write('{}  {}\n'.format(key, ' '.join(val_value)))
                #if len(val_value) == 1:
                #    fid.write('{} : {}\n'.format(key, ' '.join(val_value)))
                #else:
                #    fid.write('{} {}\n'.format(key, ' '.join(val_value)))

def check_cell_type(cell_dict):
    '''
    Check if the type matches the value, and turn uncorrected type into {'type': cell_type, 'value': cell_value}
    Parameter
    ---------
        cell_dict: dict
            1. {'item1': cell_value, 'item2': ...}
            2. {'ITEM1':{'value': cell_value}, 'item2':{'value': ...}}
            3. {'ITEM1':{'type': cell_type, 'value': cell_value}, 'item2':{'type': ...}}
    Return
    ------
        cell_dict: dict
            {'ITEM1':{'type': cell_type, 'value': cell_value}, 'item2':{'type': ...}}

    '''
    for item in cell_dict:
        cell_dict_i = cell_dict[item]
        if isinstance(cell_dict_i, dict):
            cell_type = cell_dict_i.get('type', 'UNKNOW')
            cell_value = cell_dict_i.get('value', [])
        else:
            cell_type = 'UNKNOW'
            cell_value = cell_dict_i

        (cell_type0, cell_value0) = get_cell_type(cell_value)
        if cell_type == cell_type0:
            cell_value = cell_value0
        elif (cell_type == 'BLOCK') and (cell_type0 == 'SINGLE'):
            cell_value = [cell_value0]
        elif (cell_type == 'SINGLE') and (cell_type0 == 'BLOCK'):
            cell_value = cell_value[0]       
            warnings.warn('The type is {}, the value should be a 1-dimension list, not 2-dimension. \
                Here the first elements is adopted.'.format(cell_type))
        else:
            #Type is not SINGLE or BLOCK
            cell_type = cell_type0
            cell_value = cell_value0
        cell_dict[item] = {'type': cell_type, 'value': cell_value}
    return cell_dict

def get_cell_type(cell_value):
    '''
    Get the cell type by the cell value
    Parameter
    ---------
        cell_value: int/str/float/list/tuple/np.ndarray
    Return
    ------
        cell_type:
            SINGLE: one dimension list
                The str or int or float will convert into list
            BLOCK: two dimension list
        cell_value: list
            The corrected cell value
    '''
    if isinstance(cell_value, np.ndarray):
        cell_value = cell_value.tolist()
    if isinstance(cell_value, (int, float, str)):
        cell_type = 'SINGLE'
        cell_value = [cell_value]
    elif isinstance(cell_value, (list, tuple, dict)):
        #the max dimension is 2
        if len(cell_value) == 0:
            cell_type = 'SINGLE'
            #cell_value = []
        else:
            if isinstance(cell_value, dict):
                #The dict must be atomic position
                cell_type = 'BLOCK'
            else:
                cell_value0 = cell_value[0]
                if isinstance(cell_value0, (int, float, str)):
                    cell_type = 'SINGLE'
                elif isinstance(cell_value0, (list, tuple)) and isinstance(cell_value0[0], (int, float, str)):
                    cell_type = 'BLOCK'
                else:
                    raise TypeError('Unsupported type or the dimension is larger than 2.')
    else:
        raise TypeError('Unsupported type. The type must be one of (int, float, str, list, tuple or np.ndarray)')
    return (cell_type, cell_value)

def parse_geom(geom_file):
    '''
    Parse geom file, return energy, position, stress tensor
    Rules: http://www.tcm.phy.cam.ac.uk/castep/documentation/WebHelp/content/modules/castep/expcastepfiletask.htm
    Note: All in atomic unit
        E - total energy and enthalpy
        T - temperature
        R - atomic coordinates (Cartesian)
        h - cell matrix
        hv - derivatives (velocities) of the cell matrix
        F - forces
        V - velocities
        S - stress tensor
        P - pressure

    Parameter
    ---------
        geom_file: file name
            the file name of .geom
    Return
    ------
        geom_param: dict
            geom_param = {'Energy': [], 'Position': [], 'Lattice': [], 'Force': [], 'Stress': [], 'Atom': []}
    '''
    eV2Ha = physical_constants['electron volt-hartree relationship'][0]
    Ha2J = physical_constants['Hartree energy'][0]
    bohr2m = physical_constants['Bohr radius'][0]
    bohr2ang = bohr2m * 1E10
    #atomic_pressure = Ha/bohr^3
    atomic_pressure = Ha2J / (bohr2m ** 3)
    atomic_pressure2GPa = atomic_pressure * 1E-9
    atomic_force = physical_constants['atomic unit of force'][0]
    atomic_velocity = physical_constants['atomic unit of velocity']

    geom_param = []

    #'c' here is the count
    property_map = {'E': 'Energy', 'T': 'Temperature', 'R': 'Position', 'h': 'Lattice', 'F': 'Force',
                    'V': 'Velocity', 'S': 'Stress', 'P': 'Pressure'}
    unit_map = {'E': 1./eV2Ha, 'T': 1, 'R': bohr2ang, 'h': bohr2ang, 'F': atomic_force, 
                'V': atomic_velocity, 'S': atomic_pressure2GPa, 'P': atomic_pressure2GPa, 'hv': 1}
    atomic_prop = ['R', 'F', 'V']
    #nonatomic_prop
    not_used = ['T', 'P', 'hv', 'V']
    
    with open(geom_file, 'r') as fopen:
        for eachline in fopen:
            eachline = eachline.strip()
            linei = re.split(r'(?:\s*<--\s*)', eachline)

            first_ele = linei[0]
            flag = linei[-1]

            if first_ele:
                #Not empty
                first_ele_list = re.split(r'\s+', first_ele)
                if first_ele_list[0].isdigit():
                    if int(first_ele_list[0]) == 0:
                        pass
                    else:
                        geom_param.append(property_val)
                    property_val = {'Energy': [], 'Position': [], 'Lattice': [], 'Force': [], 'Stress': [], 'Atom': [], 'Velocity': []}
                else:
                    #print(flag)
                    if flag in unit_map:
                        if flag in not_used:
                            #Do not parse it
                            pass
                        elif flag in atomic_prop:
                            values = re.split(r'\s+', linei[0])
                            num_value = [float(item)*unit_map[flag] for item in values[2:]]
                            property_val[property_map[flag]].append(num_value)
                            if flag == 'R':
                                property_val['Atom'].append(values[0])
                        else:
                            num_value = [float(item)*unit_map[flag] for item in re.split(r'\s+', linei[0])]
                            property_val[property_map[flag]].append(num_value)
    geom_param.append(property_val)
    for i, property_vali in enumerate(geom_param):
        geom_param[i]['Position'] = glue.cart2frac(property_vali['Lattice'], property_vali['Position'])
    return geom_param

def parse_castep(castep_file):
    '''
    '''
    #TODO
    cell_file = '{}.cell'.format(castep_file.split('.')[0])
    cell_parameters = parse_cell(cell_file)

    #NOTE: for mixture, it's the number of real atoms
    n_atom = len(cell_parameters['POSITIONS_FRAC']['value']['ELEMENTS'])
    distct_map = {'Energy': {'flag': 'NB est. 0K energy', 'n_empty': 0, 'n': 1}, 
                  'Position': {'flag': 'Fractional', 'n_empty': 2, 'n': n_atom}, 
                  'Lattice': {'flag': 'Unit Cell', 'n_empty': 2, 'n': 3},
                  'Force': {'flag': 'Forces', 'n_empty': 5, 'n': n_atom}, 
                  'Stress': {'flag': 'Stress Tensor', 'n_empty': 5, 'n': 3}}
    castep_param = []
    property_val = {'Energy': [], 'Position': [], 'Lattice': [], 'Force': [], 'Stress': [], 'Atom': [], 'Velocity': []}
    #The following is copy from parse_stress_from_castep
    with open(castep_file, 'r') as fopen:
        for eachline in fopen:
            eachline = eachline.strip();
            eachlinei = eachline.split("*")
            if ' Stress Tensor ' in eachlinei:
                flag = 1
            if flag :
                count = count + 1
            if count < 8:
                continue
            else:
                Stressi = re.split(r'\s+', eachlinei[1])
                Stress[i, :] = Stressi[2:5]
                i = i + 1
                if count > 9:
                    flag = 0
                    count = 1
                    i = 0
    if Stress.shape[0] == 0:
        print("Please set calculate_stress as ture.")

def get_energy_from_castep(BaseName):
    '''
    Get the energy from castep file
    Parameter
    ---------
        BaseName: str
            The filename of the cell and param file (without ext)
    Return
    ------
        energy: float
            energy in multi unit
    '''
    FileName = BaseName + ".castep"
    fopen = open(FileName, 'r')
    for eachline in fopen:
        linei = eachline.split("=")
        flag = linei[0].strip()
        flag0 = flag.split(",")[0]
        if flag0 == "Final energy":
            energy = linei[1].strip().split(' ')[0]
    fopen.close()
    Energy = float(energy)
    return Energy

def check_cell(BaseName, Settings={'FIX_ALL_CELL': {'type': 'SINGLE', 'value': ['true']}, 'CELL_CONSTRAINTS': {},
    'KPOINT_MP_GRID': {'type': 'SINGLE', 'value': np.array([16, 16, 16])}, 'KPOINTS_LIST': {}}):
    '''
    Update cell file by settings
    Parameter
    ---------
        BaseName: str
            The base name of the file
        Settings: dict
            Settings = {'CELL_KEYS': {'type': 'SINGLE/BLOCK', 'value': CELL_VALUE}}
    Return
    ------
        None
    '''
    FileName = BaseName + '.cell'
    cell_param = parse_cell(FileName)
    rem_keys = []
    for key in Settings:
        cell_value = Settings[key]
        if not cell_value:
            rem_keys.append(key)
        else:
            if isinstance(cell_value, dict):
                cell_value1 = cell_value.get('value', None)
            else:
                cell_value1 = cell_value
            if isinstance(cell_value1, np.ndarray):
                cell_value1 = cell_value1.tolist()
            if not cell_value1:
                rem_keys.append(key)
    for key in rem_keys:
        if key in cell_param:
            cell_param.pop(key)
        Settings.pop(key)
    Settings = check_cell_type(Settings)
    cell_param.update(Settings)
    write_cell(FileName, cell_param)

def check_files(BaseName, para_sets, cell_sets):
    '''
    Update cell file and param file according to sets
    Parameter
    ---------
        BaseName: str
            The base name of the file
        para_sets: dict
            The dict of parameters
        cell_sets: dict
            The dict of cell settings
    '''
    check_param(BaseName, para_sets)
    check_cell(BaseName, cell_sets)

def update_structure(BaseName, dstpath="./"):
    '''
    Update the cell file according to .geom file
    Parameter
    ---------
        BaseName: str
            The base name of the file
        dstpath: str
            Destination path
    '''
    shutil.copyfile(os.path.join(dstpath, BaseName+".cell"), os.path.join(dstpath, BaseName+"-orig.cell"))
    (position, atoms) = get_pos(BaseName)
    update_pos(BaseName, position, atoms)
    BaseVec = get_base_vec_by_geom_file(BaseName)
    write_base_vec(BaseName, BaseVec)
    shutil.copyfile(BaseName + ".cell", dstpath + "/" + BaseName + ".cell")

def update_structure_after_geom(BaseName, backupFix='bk'):
    #celldict = parse_cell(BaseName)
    shutil.copyfile("{}.cell".format(BaseName), "{}-{}.cell".format(BaseName, backupFix))
    (position, atoms) = get_pos(BaseName)
    update_pos(BaseName, position, atoms)
    BaseVec = get_base_vec_by_geom_file(BaseName)
    write_base_vec(BaseName, BaseVec)

def optimize_with_fix_direction(BaseName, sval, direction, direction2=None, flag_type='t', basevec0=None,
    NP=1, ismode='strain', tor_stress0=5.e-1, verbose=True, Nmax=5):
    '''
    Optimize
    '''
    print('----------------------Begin the optimization-------------------------')
    print('Fixed: {} {} {}'.format(*direction))
    print(type(sval))
    print('Value: {:8.4f}'.format(float(sval)))
    max_sval = abs(sval)
    if basevec0 is None:
        basevec0 = np.eye(3)

    #TODO: here only one direction and one value
    #Create and enter the folder and copy the input file for castep to the folder
    foldername = 'GeomOpt-Fix{}-Val{:3.2f}'.format(''.join(map(str, direction)), sval)
    esutils.creat_folder(foldername)
    copy_files(BaseName, Path=foldername)
    os.chdir(foldername)

    RunStr = run(NP, BaseName)
    tor_stress = 10.
    n = 1
    k = 0.9
    strains, stresses, sqrt_stress = [], [], []
    smat = iss.get_strain_or_stressBydirection(sval, direction, flag=flag_type, direction2=direction2, basevec=basevec0)
    #(eigval0, eigvec0) = iss.getPrinciple(smat, flag_ss='strain')
    eigvec0 = esutils.direction2dcm(direction, base_vec=basevec0).T

    coef = [[0., 0.]]
    stress_tmp = []
    while (tor_stress > tor_stress0):
        print('----------------Begin the {}th step'.format(n))
        strains.append(smat)
        print('The {} matrix is :'.format(ismode))
        print('{}\n'.format('\n'.join(['\t'.join(['{:<8.4f}'.format(smatii) for smatii in smati]) for smati in smat])))
        #Create and enter the folder
        folderi = 'smat{}'.format(n)
        esutils.creat_folder(folderi)
        copy_files(BaseName, Path=folderi)
        os.chdir(folderi)

        #Update the input file and run
        # Fix cell, optimize the position
        update_input_for_geomopt(BaseName, FlagEG='E')
        iss.update_input_for_is(smat, BaseName, basevec=basevec0, ismode=ismode)
        os.system(RunStr)

        #Get the stress and strain
        stress_opt = get_stress(BaseName)
        stresses.append(stress_opt)
        (PriStress, DirStress) = iss.getPrinciple(stress_opt, flag_ss='stress')
        PriStress_abs = np.abs(PriStress)
        DirStress_sort = esutils.sort_x_by_y(DirStress, PriStress_abs)
        PriStress_sort = esutils.sort_x_by_y(PriStress, PriStress_abs)
        PriStress_abs = sorted(PriStress_abs)
        if n == 1:
            stress_tmp.append(PriStress_sort[0:2])
        else:
            stress_tmp = np.append(stress_tmp, [PriStress_sort[0:2]], axis=0)
        slope_stress = max_sval/PriStress_abs[-1]

        sqrt_stress.append(PriStress_abs[0]**2 + PriStress_abs[1]**2)
        
        if n == 1:
            k1 = -slope_stress * stress_tmp[-1][0]
            k2 = -slope_stress * stress_tmp[-1][1]
        else:
            p1 = np.polyfit(np.array(coef)[:,0], stress_tmp[:,0], 1)
            p2 = np.polyfit(np.array(coef)[:,1], stress_tmp[:,1], 1)
            k1 = -p1[1]/p1[0]
            k1 = -p2[1]/p2[0]
            k1 = coef[-1][0] - (coef[-1][0] - coef[-2][0])/(stress_tmp[-1][0] - stress_tmp[-2][0]) * stress_tmp[-1][0]
            k2 = coef[-1][1] - (coef[-1][1] - coef[-2][1])/(stress_tmp[-1][1] - stress_tmp[-2][1]) * stress_tmp[-1][1]
        coef = np.append(coef, [[k1, k2]], axis=0)
        eigval = [k1, k2 ,sval]
        if verbose:
            print('Strains:')
            print(np.array(eigval))
            print('Direction:')
            print(np.array(DirStress_sort))

        #smat = iss.get_smat_by_eig(eigval=eigval, eigvec=DirStress_sort)
        smat = iss.get_smat_by_eig(eigval=eigval, eigvec=eigvec0)

        tor_stress = max(abs(PriStress[[0, 1]]))
        print('tor_stress: {}'.format(tor_stress))
        os.chdir('..')
        if n > Nmax:
            min_stress = min(sqrt_stress)
            ind_min = sqrt_stress.index(min_stress)
            print(stresses)
            print(ind_min)
            stress_opt = stresses[ind_min]
            smat = strains[ind_min]
            break
        n = n + 1

    os.chdir('..')
    print('----------------------End of the optimization-------------------------\n')
    #stress and strain
    return (stress_opt, smat)

