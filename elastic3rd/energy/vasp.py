#!python
#
import re
import os
import shutil
import linecache
import numpy as np
import elastic3rd.energy.glue as glue
import elastic3rd.esutils as esutils

def get_base_vec(BaseName=None, FileName='POSCAR'):
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
    BaseVec = np.zeros((3, 3))
    #FileName = "POSCAR"
    fopen = open(FileName, 'r')
    count = 0
    for eachline in fopen:
        if count > 1:
            linei = eachline.strip('\n').strip()
            linei = re.split(r'\s+', linei)
            for j in range(0, 3):
                BaseVec[count - 2][j] = float(linei[j])
        if count > 3:
            return BaseVec
        count = count + 1
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
    FileName = "POSCAR"
    fopen = open(FileName, 'r')
    tmpopen = open('tmpfile', 'a')
    count = 1
    lines = []
    for eachline in fopen:
        if (count > 2) and (count < 6):
            for j in range(0, 3):
                tmpopen.write("      ")
                tmpopen.write("%.15f" % BaseVec[count - 3][j])
            tmpopen.write("\n")
        else:
            tmpopen.write(eachline)
        count = count + 1
    fopen.close()
    tmpopen.close()
    os.remove(FileName)
    os.rename('tmpfile', FileName)

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
    RunStr = "mpirun -np " + str(int(NP)) + " vasp_std"
    return RunStr

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
    FileName = "OSZICAR"
    fopen = open(FileName, 'r')
    for eachline in fopen:
        linei = eachline.split("=")
        if len(linei) > 2:
            energy = linei[2].strip().split(" ")[0]
    fopen.close()
    energy = float(energy)
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
    Stress = np.zeros((3, 3))
    FileName = "OUTCAR"
    fopen = open(FileName, 'r')
    flag = 0
    count = 1
    i = 0
    for eachline in fopen:
        eachline = eachline.strip('\n')
        eachline = eachline.strip();
        if 'in kB' in eachline:
            flag = 1
            stress_vec = re.split(r'\s+', eachline)[2:]
    fopen.close()
    if flag == 0:
        print("ISIF should be 2.")
    else:
        Stress[0, 0] = -float(stress_vec[0])/10.0
        Stress[1, 1] = -float(stress_vec[1])/10.0
        Stress[2, 2] = -float(stress_vec[2])/10.0
        Stress[0, 1] = -float(stress_vec[3])/10.0
        Stress[1, 0] = Stress[0, 1]
        Stress[1, 2] = -float(stress_vec[4])/10.0
        Stress[2, 1] = Stress[1, 2]
        Stress[0, 2] = -float(stress_vec[5])/10.0
        Stress[2, 0] = Stress[0, 2]
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
    shutil.copyfile("INCAR", Path + "/" + "INCAR")
    shutil.copyfile("POSCAR", Path + "/" + "POSCAR")
    shutil.copyfile("POTCAR", Path + "/" + "POTCAR")
    shutil.copyfile("KPOINTS", Path + "/" + "KPOINTS")

def parse_poscar(filename='POSCAR'):
    '''
    Parse the POSCAR
        Note: there should no extra space line and not support MD
    Parameter
    ---------
        filename: str
            The filename of POSCAR like
    Return
    ------
        poscar: dict
            The parameters poscar, it has following keys
            {'SYSTEM': 'test', 'SCARLE_FACTOR': 1.0, 'LATTICE': [], 'ATOMS': [], 
              'POSITION': [], 'ELE': [], 'ELE_NUM': [], 'POS_TYPE': 'D'}
    '''
    count_row = 0
    poscar = {'SYSTEM': 'test', 'SCARLE_FACTOR': 1.0, 'LATTICE': [], 'ATOMS': [], 
              'POSITION': [], 'ELE': [], 'ELE_NUM': [], 'POS_TYPE': 'D'}
    with open(filename, 'r') as fopen:
        for eachline in fopen:
            eachline = eachline.strip()
            if count_row == 0:
                poscar['SYSTEM'] = eachline
            elif count_row == 1:
                poscar['SCARLE_FACTOR'] = float(eachline)
            elif count_row > 1 and count_row < 5:
                poscar['LATTICE'].append(re.split(r'\s+', eachline)[0:3])
            elif count_row == 5:
                ele = re.split(r'\s+', eachline)
                poscar['ELE'] = ele
            elif count_row == 6:
                ele_num = [int(item) for item in re.split(r'\s+', eachline)]
                poscar['ELE_NUM'] = ele_num
                ni_ele = 0
                for i in range(0, len(ele_num)):
                    ni_ele0 = ni_ele
                    ni_ele = ni_ele + ele_num[i]
                    for j in range(ni_ele0, ni_ele):
                        poscar['ATOMS'].append(ele[i])
            elif count_row == 7:
                poscar['POS_TYPE'] = eachline.upper()
            elif count_row > 7 and count_row < 8 + len(poscar['ATOMS']):
                poscar['POSITION'].append(re.split(r'\s+', eachline)[0:3])
            count_row = count_row + 1
    if poscar['POS_TYPE'].startswith(('C', 'K')):
        poscar['POSITION'] = glue.cart2frac(poscar['LATTICE'], poscar['POSITION'])
        poscar['POS_TYPE'] = 'Direct'
    return poscar

def write_poscar(poscar, filename='POSCAR'):
    '''
    Write out POSCAR
    Parameter
    ---------
        poscar: dict
            The dict of poscar, the format is the same the the result of parse_poscar
        filename: str
            The filename of POSCAR-like file
    Return
    ------
        None
    '''
    key_order = ['SYSTEM', 'SCARLE_FACTOR', 'LATTICE', 'ELE', 'ELE_NUM', 'POS_TYPE', 'POSITION']
    with open(filename, 'w+') as fopen:
        for key in key_order:
            value = poscar[key]
            if isinstance(value, (str, int, float)):
                fopen.write("{}\n".format(value))
            elif isinstance(value, (list, tuple)):
                flag_list = True
                for i, valuei in enumerate(value):
                    if isinstance(valuei, (str, int, float)):
                        flag_list = True
                        value[i] = str(valuei)
                    elif isinstance(valuei, (list, tuple)):
                        flag_list = False
                        fopen.write("  {}\n".format(' '.join(valuei)))
                if flag_list:
                    fopen.write("  {}\n".format(' '.join(value)))
            elif isinstance(value, np.ndarray):
                value_shape = value.shape
                value = value.tolist()
                if len(value_shape) == 1:
                    fopen.write("  {}\n".format(' '.join(value)))
                elif len(value_shape) == 2:
                    for valuei in value:
                        fopen.write("  {}\n".format(' '.join(valuei)))

def get_pos(filename='POSCAR'):
    '''
    Get the position, the formation is below:
      [[ele1, pos1, pos2, pos3],
       [ele2, pos1, pos2, pos3], 
       [], [], ..]
    Parameter
    ---------
        filename: str
            The filename of POSCAR-like file
    Return
    ------
        pos: list
            The list of atomic position
    '''
    position_tmp = []
    ele_list = []
    ni_ele = 0
    flag = 0
    #flag_be  the flag of begin or end of the position
    count_row = 0
    count_ele = 0
    with open(filename, 'r') as fopen:
        for eachline in fopen:
            eachline = eachline.strip()
            #eachlinei = eachline.split("*")
            if count_row == 5:
                ele = re.split(r'\s+', eachline)
            if count_row == 6:
                ele_num = re.split(r'\s+', eachline)
                for i in range(0, len(ele_num)):
                    ni_ele0 = ni_ele
                    ni_ele = ni_ele + int(ele_num[i])
                    for j in range(ni_ele0, ni_ele):
                        ele_list.append(ele[i])
            if 'Direct' in eachline:
                flag = 1
                continue
            if flag :
                if not eachline:
                    # is empty
                    flag = 0
                else:
                    pos = re.split(r'\s+', eachline)
                    #print(pos)
                    posi = [ele_list[count_ele]]
                    posi.extend(pos)
                    position_tmp.append(posi)
                    position = position_tmp                
                    count_ele = count_ele + 1
            count_row = count_row + 1
    return position

def add_perturbation(BaseName=None, tor=2.0E-10, filename='POSCAR', PertMethod='random', PertType='relative'):
    '''
    Add a small perturbation in the position in the POSCAR file

    Parameters
    ----------
        PertType: str
            The type of perturbation. Only first letter works, case insensitive
                abs: absolute perturbation, in unit of A
                relative: relative perturbation, in unit of a(lattice constant)
    Return
    ------
        None
    '''
    poscar = parse_poscar(filename=filename)
    #the position is converted into fration when parsing
    pos = np.array(poscar['POSITION'], dtype=np.float64)
    pert = esutils.gen_perturbation(pos, tor=tor, PertMethod=PertMethod)
    if PertType.lower().startswith('a'):
        pert = glue.cart2frac(np.array(poscar['LATTICE'], dtype=np.float64) * float(poscar['SCARLE_FACTOR']), pert)
    elif PertType.lower().startswith('r'):
        pass
    else:
        raise ValueError('Unsupported type of perturbation. Currently support the abs or relative')
    pos = pos + pert
    new_pos = [[format(float(posij), '.15f') for posij in posi] for posi in pos]
    poscar['POSITION'] = new_pos
    write_poscar(poscar, filename=filename)

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
    os.rename(os.path.join(dstpath, "POSCAR"), os.path.join(dstpath, "POSCAR-org"))
    shutil.copyfile("CONTCAR", dstpath + "POSCAR")

def check_cell(BaseName, flag = "E", Settings = {'KPOINTS' : '16 16 16'}):
    pass

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
    check_incar(incar_setting)

def check_kpoints(kpoints, kpoint_folder= "."):
    '''
    Update KPOINTS
        support two mode, 1.Only one number, a means a x a x a
                          2. three numbers: a, b, c, means a x b x c
    Parameter
    ---------
        kpoints: int or list
            the kpoint parameters
        kpoint_folder: path-like str
            the path of KPOINTS file 
    Return
    ------
        None
    '''
    if (type(kpoints) is int) or (type(kpoints) is float):
        kpoints = str(int(kpoints))
        kpoints = [kpoints, kpoints, kpoints]
    elif type(kpoints) is str :
        try:
            kpoints = str(int(float(kpoints)))
            kpoints = [kpoints, kpoints, kpoints]
        except Exception as e:
            raise e
    elif type(kpoints) is list and len(kpoints) == 3:
        pass
    elif type(kpoints) is np.ndarray and len(kpoints) == 3:
        pass
    else:
        raise IOError('The input of the type of length of kpoints is not correct')
    kp_template = open(kpoint_folder + "/KPOINTS", 'r')
    kp_file = open(kpoint_folder + "/KPOINTStmp", 'w+')
    k_count = 0
    for eachline in kp_template:
        if k_count == 3:
            kp_file.write("  %d  %d  %d\n" % (int(kpoints[0]), int(kpoints[1]), int(kpoints[2])))
        else:
            kp_file.write(eachline)
        k_count = k_count + 1
    kp_template.close()
    kp_file.close()
    os.remove(kpoint_folder + "/KPOINTS")
    os.rename(kpoint_folder + "/KPOINTStmp", kpoint_folder +  "/KPOINTS")

def check_incar(Settings = {'ISIF' : '3'}):
    '''
    Update INCAR
    Parameter
    ---------
        Settings: dict
            The dict of INCAR parameters
    Return
    ------
        None
    '''
    #Check the param
    #Checking task==GeometryOptimization, calculate_stress==true, cut_off_energy=360
    (INCAR_dict, key_order) = incar_parser(INCAR = "./INCAR")
    if INCAR_dict['ISIF'] == "3":
        INCAR_dict["ENCUT"] = str(1.3 * float(INCAR_dict["ENCUT"]))
    for key in Settings:
        INCAR_dict[key] = Settings[key]
    incar_write(INCAR_dict, key_order, dst_folder = ".")

def check_poscar(pos_sets = {}):
    pass

def incar_parser(INCAR = "./template/INCAR"):
    '''
    The test INCAR file is taken from https://cms.mpi.univie.ac.at/vasp/guide/node91.html
    Input: INCAR
        The filename of INCAR 
    Output :
        INCAR_dict : A dict contains all the parameters(key) in the INCAR file and corresponding settings(value)
        key_order : A list contains all the parameters as the order of the original INCAR file
    '''
    INCAR_dict = {}
    key_order = []
    fopen = open(INCAR, 'r')
    for eachline in fopen:
        eachline = eachline.strip('\n')
        #print(eachline)
        if not is_blank_comments(eachline):
            eachline = eachline.split("#")[0].strip()  #the content after # is neglected
            linei = eachline.split(";")
            for tags in linei:
                #neglect lines not tag-value pairs, such as : Start parameter for this Run:
                tagsi = tags.split("=")
                if len(tagsi) > 1:
                    incar_key = re.split(r'\s+', tagsi[0].strip())[0]
                    incar_val = tagsi[1].strip()
                    #incar_val = re.split('\s+', tagsi[1].strip())[0]
                    INCAR_dict[incar_key] = incar_val
                    key_order.append(incar_key)          
    return INCAR_dict, key_order

def incar_write(INCAR_dict, key_order, dst_folder = "."):
    """
    Write INCAR using INCAR_dict
    Parameter
    ---------
        INCAR_dict: dict
            The dict of the INCAR, if the value if "" or None, don't write this key (delete)
        key_order: list
            The list containing the key of the INCAR
    Return
    ------
        None
    """
    incar = open(dst_folder + '/INCAR', 'w+')
    for keys in key_order:
        if INCAR_dict[keys]:
            #Not empty
            incar.write("%s  =  %s\n" % (keys, INCAR_dict[keys]))
    incar.close()

def is_blank_comments(line):
    '''
    Judge if the line is comment of not
    Parameter
    ---------
        line: str
            The line in INCAR file
    Return
    ------
        flag: int
          1 for blank or comments
          0 for not
    '''
    flag = 0
    if len(line) < 1:   
        #blank line
        flag = 1
    elif line[0] == "#":    
        #comment line
        flag = 1
    return flag

def check_files(BaseName, incar_sets, kpoints_sets):
    '''
    Update INCAR and KPOINTS
    Parameter
    ---------
        BaseName: str
            Any string, just for consist with castep
        incar_sets: dict
            The dict of INCAR parameters
        kpoints_sets: dict
            The dict for KPOINTS setting
            For consistent with CASTEP, should contain KPOINTS key in kpoints_sets dict
    Return
    ------
        None
    '''
    check_incar(incar_sets)
    if 'KPOINTS' in kpoints_sets:
        kpoints = kpoints_sets['KPOINTS']
        check_kpoints(kpoints)
    #check_poscar(pos_sets)

def rem_files(BaseName):
    '''
    Remove some large files, e.g. WAVECAR
    '''
    #os.remove(path)
    os.system("rm WAVECAR CHG* vasprun.xml IBZ* EIG* POT*")

