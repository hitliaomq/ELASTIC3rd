#!python
#
import logging
import numpy as np
from scipy.constants import pi
import time
import json
import yaml
import math
import re
import os
VERSION = __import__('elastic3rd').__version__
CODEDATE = __import__('elastic3rd').__codedate__

def init_log():
    logger = logging.getLogger(name='e3rd')
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(message)s")

    #To file
    fh = logging.FileHandler('OUTPUT.e3rd')
    fh.setFormatter(fmt)
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger
logger = init_log()

def read_input(INPUT = "INPUT.e3rd"):
    '''
    Read the input file and set default value for Input parameters
    Parameters
    ----------
        INPUT: str
            The file name of the input file
                For normal text, The key and value in the INPUT file should be divided by space
                In addition, it support json and yaml file
                The key is capital sensitive and the value is not
    Return
    ------
        ParaIn: dict
            The dict of input parameters
    '''
    ParaIn = {'Continue': 0, 'EnergyCode' : 'CASTEP', 'EnergyRun' : 1, 'CrystalType': 'C1', 'Ord': 3, 
        'MaxStrain' : 5.0, 'NP' : 1, 'STEPS' : 5, 'BaseName' : 'Si', 'CalMode': 'W', 'FlagSE': 'E', 
        'GenSM' : 1, 'CoefMethod': 'fit', 'Acc': 4, 'Geom': 0, 'AddPert': 0, 'PertMethod': 'random',
        'StableEst': 1, 'StrainType': 'Uniform'}
    INPUTtype = {'Continue': 'int', 'EnergyCode' : 'str', 'EnergyRun' : 'int', 'CrystalType': 'str', 'Ord': 'int', 
        'MaxStrain' : 'float', 'NP' : 'int', 'STEPS' : 'int', 'BaseName' : 'str', 'CalMode': 'str', 'FlagSE': 'str', 
        'GenSM' : 'int', 'CoefMethod': 'str', 'Acc': 'int', 'MinStrain': 'float', 'Geom': 'int', 'AddPert': 'int',
        'PertMethod': 'str', 'StableEst': 'int', 'StrainType': 'str', 'ISModel': 'str', 'ISStep': 'float',
        'ISMethod': 'str'}
    with open(INPUT, 'r') as fopen:
        if INPUT.lower().endswith('.json'):
            ParaIn = json.load(fopen)
        elif INPUT.lower().endswith('.yaml'):
            ParaIn = yaml.load(fopen)
        else:
            for eachline in fopen:
                eachline = eachline.strip()
                linei = re.split(r'\s+', eachline)
                if len(linei) > 1 :
                    flag = linei[0]
                    value = linei[1]
                    ParaIn[flag] = value
    for item in ParaIn:
        key_type = INPUTtype[item]
        if key_type == 'str':
            ParaIn[item] = str(ParaIn[item])
        elif key_type == 'int':
            ParaIn[item] = int(ParaIn[item])
        elif key_type == 'float':
            ParaIn[item] = float(ParaIn[item])
    ParaIn['EnergyCode'] = ParaIn['EnergyCode'].lower()
    ParaIn['CrystalType'] = crystaltypemap(ParaIn['CrystalType'])
    if 'MinStrain' not in ParaIn:
        ParaIn['MinStrain'] = -ParaIn['MaxStrain']
    return ParaIn

def crystaltypemap(crystaltype):
    crystaltype = crystaltype.lower()
    crystaltype_dict = {"isotropic": "i", "i" : "i", "cubic1" : "c1", "c1" : "c1", 
                        "cubic2" : "c2", "c2" : "c2", "hexagonal1" : "h1", "h1" : "h1", 
                        "hexagonal2" : "h2", "h2" : "h2", "rhombohedral1" : "r1", "r1" : "r1",
                        "rhombohedral2" : "r2", "r2" : "r2", "tetragonal1" : "t1", "t1" : "t1", 
                        "tetragonal2" : "t2", "t2" : "t2", "orthorhombic" : "o", "o" : "o", 
                        "monoclinic" : "m", "m" : "m", "triclinic" : "n", "n" : "n"}
    crystaltype = crystaltype_dict[crystaltype]
    return crystaltype

def crystaltype2lattice(crystaltype, lat):
    '''
    Convert lat(independent lattice constants) into lattice
    '''
    crystaltype = crystaltypemap(crystaltype)
    if crystaltype in ['c1', 'c2', 'i']:
        lattice = [lat[0], lat[0], lat[0], 90, 90, 90]
    elif crystaltype in ['h1', 'h2']:
        lattice = [lat[0], lat[0], lat[1], 90, 90, 120]
    elif crystaltype in ['t1', 't2']:
        lattice = [lat[0], lat[0], lat[1], 90, 90, 90]
    elif crystaltype in ['r1', 'r2']:
        lattice = [lat[0], lat[0], lat[0], lat[1], lat[1], lat[1]]
    elif crystaltype in ['o']:
        lattice = [lat[0], lat[1], lat[2], 90, 90, 90]
    elif crystaltype in ['m']:
        lattice = [lat[0], lat[1], lat[2], 90, lat[3], 90]
    else:
        lattice = lat
    return lattice

def write_input(ParaIn, FileName = "INPUT.e3rd"):
    '''
    Write INPUT file according ParaIn dict
    Parameters
    ----------
        ParaIn: dict
            The dict of INPUT parameters
        FileName: str
            The file name of the INPUT file
    Return
    ------
        None, but it will write out an input file
    '''
    with open(FileName, 'w+') as fopen:
        for key in ParaIn:
            fopen.write('{}  {}\n'.format(key, ParaIn[key]))

def check_input(INPUT='INPUT.e3rd', NewParaIn={}):
    '''
    Update current INPUT file according New settings
    '''
    ParaIn = read_input(INPUT=INPUT)
    ParaIn.update(NewParaIn)
    write_input(ParaIn, FileName=INPUT)

def read_strainmode(STRAINMODE = "STRAINMODE"):
    '''
    Read strain mode file
    Parameters
    ----------
        STRAINMODE: str
            The filename of strain mode file, STRAINMODE by default
    Return
    ------
        strainmode: np.ndarray
            strain modes, 2D matrix, e.g. [[1, 2, -0.5, 1, 1, 0]]
    '''
    strainmodetmp = np.loadtxt(STRAINMODE)
    m = strainmodetmp.shape
    if len(m) == 1:
        strainmode = np.zeros((1, 6))
        strainmode[0, :] = strainmodetmp
    else:
        strainmode = strainmodetmp
    return strainmode

def print_parain(ParaIn):
    '''
    Print the input parameters
    Parameters
    ----------
        ParaIn: dict
            The dict of input parameters
    Return
    ------
        None
    '''
    if ParaIn['Continue']:
        RunMode = "CONTINUE"
    else:
        RunMode = "NEW"
    if ParaIn['FlagSE'].lower() == 'e':
        Method = "Strain-Energy method"
    else:
        Method = "Strain-Stress method"   
    logger.info("The Runing mode is : " + RunMode)
    logger.info("The BaseName is : " + ParaIn['BaseName'])
    logger.info("The Method used here is: " + Method)
    if ParaIn['CalMode'].lower() == 's':
        logger.info("The strain mode(s) is(are) specified.")
    logger.info("The crystal symmetry is : " + ParaIn['CrystalType'])
    logger.info("The order of elastic is : " + str(ParaIn['Ord']))
    logger.info("The code for calculating the energy : " + ParaIn['EnergyCode'].upper())
    logger.info("The maximum of the strain : " + str(ParaIn['MaxStrain']))
    logger.info("The steps run in every mode : " + str(ParaIn['STEPS']))
    logger.info("The core used in energy calculate : " + str(ParaIn['NP']) + "\n")
    #print("\n")

def creat_folder(path):
    '''
    Create folder if not exist, if exist, skip
    Parameters
    ----------
        path: str
            The path of the folder
    Return
    ------
        None
    '''
    if os.path.exists(path):
        logger.info("Warning: The folder(" + path + ") exists!\n")
        return 0
    else:
        os.mkdir(path)
        return 1

def print_e(E):
    '''
    Print the energy
    Parameters
    ----------
        E: list
            The energy in different unit, it is generated by elastic3rd.energy.glue.multi_energy function
    Return
    ------
        None
    '''
    logger.info("       eV             Hartree            kJ/mol            kcal/mol")
    E_str = []
    for i in range(0, 4):
        strparam = "%.10f"
        E_str.append(strparam % E[i])
    logger.info("  ".join(E_str))

def gen_perturbation(pos, tor=2.0E-10, PertMethod='random'):
    '''
    Generate perturbation
    Parameter
    ---------
        pos: list or np.ndarray
            The name of the *cell file
        tor: float
            The amplitude of the pertubation, usually a small value
        PertMethod: str
            The method for perturbation, random/uniform, only the first letter works
    Return
    ------
        pert: np.ndarray
            The perturbation matrix
    '''
    pos = np.array(pos)
    if os.path.exists('PERTUR_MATRIX'):
        pert = np.loadtxt('PERTUR_MATRIX')
    else:
        if PertMethod.lower().startswith('r'):
            pert = np.random.random(pos.shape) * tor - tor/2.
        elif PertMethod.lower().startswith('u'):
            pert = np.tile(np.random.random(pos.shape[1])*tor, (pos.shape[0], 1))
        np.savetxt('PERTUR_MATRIX', pert)
    return pert

def gen_strain_list(ParaIn):
    if os.path.exists('STRAINAMP'):
        StrainList = np.loadtxt('STRAINAMP')
    else:
        Step = int(ParaIn['STEPS'])
        MaxStrain = float(ParaIn['MaxStrain'])
        if "MinStrain" in ParaIn:
            MinStrain = float(ParaIn['MinStrain'])
        else:
            MinStrain = - MaxStrain
        #if Step % 2 == 0:
        #    Step = Step + 1
        if Step < 3:
            logger.error("ERROR:The STEPS is too small, at least 3!")
            return
        #StepStrain = round((MaxStrain - MinStrain)/(Step - 1), 2)
        #EndStrain = MaxStrain + StepStrain
        #StrainList1 = np.arange(0, MaxStrain, StepStrain)
        #StrainList2 = np.arange(MinStrain, 0, StepStrain)
        #StrainList = np.arange(MinStrain, MaxStrain+StepStrain, StepStrain)
        StrainType = ParaIn['StrainType'].upper()
        if StrainType.startswith('C'):
            StrainList = [MaxStrain * np.cos((2*i + 1)*pi/(2*Step)) for i in range(Step)]
        else:
            StrainList = np.linspace(MinStrain, MaxStrain, Step)
        StrainList = np.array([round(StrainListi, 2) for StrainListi in StrainList])
    return StrainList

def iscontinue(Path, ModeOrStrain, flag_se = "e"):
    #In mode step, the energy file named as Energy_Mode
    #In strain step, the energy file named as Energy_Strain
    #ModeOrStrain equal to Mode or Strain
    if flag_se == "e":
        filename = Path + "/" + "Energy_" + ModeOrStrain + ".txt"
    elif flag_se == "s":
        filename = Path + "/" + "Stress_" + ModeOrStrain + ".txt"
    if os.path.exists(filename):
        flag = 1
    else:
        flag = 0
    return flag

def sort_x_by_y(x, y):
    """
    Sort a list of x in the order of sorting y
    [This function is taken from dfttk]
    """
    return [xx for _, xx in sorted(zip(y, x), key=lambda pair: pair[0])]

def print_logo():
    logger.info("+=============================================================+")
    logger.info("| EEEEE L       A   SSSSS TTTTT IIIII CCCCC 33333 RRRRR DDD   |")
    logger.info("| E     L      A A  S       T     I   C         3 R   R D  D  |")
    logger.info("| EEEE  L     A   A SSSSS   T     I   C      3333 RRRRR D   D |")
    logger.info("| E     L     AAAAA     S   T     I   C         3 R  R  D  D  |")
    logger.info("| EEEEE LLLLL A   A SSSSS   T   IIIII CCCCC 33333 R   R DDD   |")
    logger.info("+-------------------------------------------------------------+")
    logger.info("|             Version: ELASTIC3RD  version {}              |".format(VERSION))
    logger.info("|                Date: {}                             |".format(CODEDATE))
    logger.info("|              Author: Liao Mingqing                          |")
    logger.info("|              E_mail: liaomq1900127@163.com                  |")
    logger.info("|     FGMS @ Harbin Institute of Technology(HIT)              |")
    logger.info("|  Jiangsu University of Science and Technolgy(JUST)          |")
    logger.info("|    URL:https://github.com/hitliaomq/ELASTIC3rd              |")
    logger.info("|Please Cite:                                                 |")
    logger.info("|M.Q. Liao, et al., Elastic3rd: A tool for calculating third- |")
    logger.info("|  order elastic constants from first-principles calculations.|")
    logger.info("|  Computer Physics Communications 261 (2021) 107777          |")
    logger.info("|If the strain-stress method is used, please also cite:       |")
    logger.info("|M.Q. Liao, et al., A high-efficient strain-stress method for |")
    logger.info("|  calculating higher- order elastic constants from first-    |")
    logger.info("|  principles.                                                |")
    logger.info("|  Computer Physics Communications 280 (2022) 108478          |")
    logger.info("+=============================================================+")

def rotationByEuler(theta, euler_type='zxz', axis='inner', clockwise=False, degree=False):
    '''
    Rotation Matrix according to Euler angles
    Parameters
    ----------
        theta: float or list
            if the size is 1x2, the final one is zeros
        euler_type: str
            'zxz' means :rotation along z, and then along x and finally along z
        axis: 'inner' or 'outer'
            'inner': Rotation around inner axis (multify at right side)
            'outer': Rotation around outer(fixed) axis (multify at left side)
        clockwise: bool
            True for clock wise, False for inverse clock wise
        degree: bool
            True for degree in theta, and False for radian in theta
    '''
    R = []
    euler_type = euler_type.lower()
    if isinstance(theta, (int, float)):
        theta = [float(theta), 0., 0.]
    elif isinstance(theta, (list, tuple)) and len(theta) == 2:
        theta = list(theta).append(0.)
    if degree:
        theta = [math.radians(item) for item in theta]
    if clockwise:
        theta = [-item for item in theta]
    
    for i, flagi in enumerate(euler_type):
        R.append(rotationMatrix(theta[i], flag=flagi))

    if axis.lower().startswith('i'):
        R_Euler = np.dot(R[2], np.dot(R[1], R[0]))
    elif axis.lower().startswith('o'):
        R_Euler = np.dot(R[0], np.dot(R[1], R[2]))
    else:
        raise TypeError('The axis should start with i or o.')
    
    return R_Euler

def rotationMatrix(theta, flag='x', degree=False):
    if degree:
        theta = math.radians(theta)
    if flag == 'x':
        R = np.array([[1, 0, 0],
                      [0, np.cos(theta), np.sin(theta)], 
                      [0, -np.sin(theta), np.cos(theta)]])
    elif flag == 'y':
        R = np.array([[np.cos(theta), 0, -np.sin(theta)],
                      [0, 1, 0],
                      [np.sin(theta), 0, np.cos(theta)]])
    elif flag == 'z':
        R = np.array([[np.cos(theta), np.sin(theta), 0],
                      [-np.sin(theta), np.cos(theta), 0],
                      [0, 0, 1]])
    else:
        raise ValueError('The flag must be x or y or z')
    return R

def norm_vector(vector):
    return np.array(vector) / np.linalg.norm(vector)

def direction2dcm(direction, direction2=None, base_vec=None, degree=False, tor=1e-7):
    '''
    Convert direction to dcm (direction cosine matrix)
    Parameter
    ---------
    Return
    ------
        dcm: np.array
            NOTE: this is a column matrix, 
    '''
    theta = direction2Euler(direction, direction2=direction2, base_vec=base_vec, degree=degree, tor=tor)
    dcm = rotationByEuler(theta, euler_type='zxz', axis='inner', clockwise=False, degree=degree)
    return dcm

def direction2Euler(direction, direction2=None, base_vec=None, degree=False, tor=1e-7):
    '''
    Get the Euler angle from directions (Rotate current axis to the directions)
        direction is the z axis
        direction2 is the x axis
    Note: The Euler type is 'zxz'
    Parameter
    ---------
        direction: list
            The z axis
        direction2: None or list
            None: the final angle is zero
            list: the x axis
        base_vec: 3x3 list
            The current crystal vector in which the directions are defined
        tor: float
            The torlerance
    Return
    ------
    '''
    if base_vec is None:
        base_vec = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
    #convert to cart coord
    direction = np.dot(direction, base_vec)

    direction = norm_vector(direction)
    (x, y, z) = direction
    theta1 = np.arccos(z)
    if abs(z) > 1. - tor:
        theta0 = 0.
    else:
        theta0 = np.arccos(-y / np.sqrt(1 - z ** 2))
    if direction[0] < 0:
        theta0 = 2. * np.pi - theta0

    if direction2 is None:
        theta2 = 0.
    else:
        direction2 = np.dot(direction2, base_vec)
        direction2 = norm_vector(direction2)
        newBaseVec = np.dot(base_vec, rotationByEuler(theta=[theta0, theta1, 0.], euler_type='zxz', axis='inner'))
        newX = norm_vector(newBaseVec[0])
        theta2 = np.arccos(np.dot(direction2, newX))
    theta = [theta0, theta1, theta2]
    if degree:
        theta = [math.degrees(item) for item in theta]
    return theta

def get_cos_mat(direction, direction2, base_vec=None):
    '''
    Parameter
    ---------
        direction: list
            z axis
        direction2: list
            x axis
    '''
    if base_vec is None:
        base_vec = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
    direction = norm_vector(direction)
    direction2 = norm_vector(direction2)
    y = np.cross(direction, direction2)
    logger.info(y)
    newAxes = np.array([direction2, y, direction])
    dcm = [[np.dot(basevi, newaxesi) for newaxesi in newAxes] for basevi in base_vec]
    return np.array(dcm).T

def GetNowTime():
    logger.info(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(time.time())))
    return time.time()

def rotate_around_arb_axes(point, direction, alpha, point0=None, degree=False):
    '''
    Rotation point around any axes
    Ref. https://blog.csdn.net/maple_2014/article/details/104443928
    '''
    if degree:
        alpha = math.radians(alpha)
    if point0 is None:
        (x0, y0, z0) = (0., 0., 0.)

    direction = norm_vector(direction)
    (nx, ny, nz) = direction
    cosa = math.cos(alpha)
    sina = math.sin(alpha)
    k = 1. - cosa
    m = nx * x0 + ny * y0 + nz * z0
    RM = [[k*nx**2+cosa, k*nx*ny-nz*sina, k*nx*nz+ny*sina, k*(x0-nx*m)+(nz*y0-ny*z0)*sina],
          [k*nx*ny+nz*sina, k*ny**2+cosa, k*ny*nz-nx*sina, k*(y0-ny*m)+(nx*z0-nz*x0)*sina],
          [k*nx*nz-ny*sina, k*ny*nz+nx*sina, k*nz**2+cosa, k*(z0-nz*m)+(ny*x0-nx*y0)*sina],
          [0., 0., 0., 1.]]
    point_tmp = np.array(point).tolist() + [1.]

    point_new = np.dot(RM, point_tmp)
    return point_new[0:3]

def mat2vec(mat):
    mat = np.array(mat)
    return [mat[0, 0], mat[1, 1], mat[2, 2], mat[1, 2], mat[0, 2], mat[0, 1]]

def vec2mat(vec):
    mat = np.array([[vec[0], vec[5], vec[4]],
                    [vec[5], vec[1], vec[3]],
                    [vec[4], vec[3], vec[2]]])
    return mat


'''
direction = [1, 1, 0]  #[1, 1, 0]
direction2 = [1, -1, 0]
theta = direction2Euler(direction, direction2=direction2, base_vec=None, degree=True)
print(theta)
print('')


flag = 'x'
theta = [90, 90, 90]
base_vec = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
RM = rotationByEuler(theta, euler_type='zxz', axis='i', clockwise=False, degree=True)

new_bv = np.dot(base_vec, RM)
print(new_bv)
print('')

direction = [1, 1, 1]
direction2 = [1, -1, 0]
dcm = get_cos_mat(direction, direction2)
print(dcm)

theta = direction2Euler(direction, direction2=direction2)
print(theta)
RM = rotationByEuler(theta, euler_type='zxz', axis='i')
print(RM)
'''

'''
point =[[0.5790674639430065,    0.5793349865922877,    0.2090882713534079],
        [0.5792071164814082,    0.2090787976170879,    0.5792159412739355],
        [0.2091508592243757,    0.5791456447669229,    0.5792209509198858],
        [0.7908499468446257,    0.4207457365466742,    0.4208879608769175],
        [0.4209055972495299,    0.7909237156398796,    0.4206711851719773],
        [0.4208204543152991,    0.4207768400004352,    0.7909088330090087],
        [0.3507312423217488,    0.3508416736866657,    0.1144113709055920],
        [0.1143596955316159,    0.3507161913987690,    0.3508329601852644],
        [0.3508510228920914,    0.1144195436759206,    0.3507471462823143],
        [0.8856398593687324,    0.6491730714356058,    0.6492788828174293],
        [0.6492610664228414,    0.8855783892924517,    0.6491452273331835],
        [0.6491556754047221,    0.6492654093472984,    0.8855912698710875]]

alpha = 60
direction = [1., 1., 1.]
for pointi in point:
    newpoint = rotate_around_arb_axes(pointi, direction, alpha, degree=True)
    print(newpoint)
'''

