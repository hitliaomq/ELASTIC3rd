# -*- coding: utf-8 -*-
# This is the script run elastic3rd in the terminal
# Simply run     elastic3rd
import sys
import os
import argparse
import importlib
from elastic3rd import __version__
import elastic3rd.elastic
import elastic3rd.esutils as esutils
import elastic3rd.post.post as espost
import elastic3rd.crystal.lattice as crylat
import elastic3rd.symmetry.symmetry as essym

def run(args):
    INPUT = args.INPUT
    elastic3rd.elastic.elastic3(INPUT=INPUT)

def post(args):
    INPUT = args.INPUT
    N = args.N
    STRAINMODE = args.STRAINMODE
    FIG = args.FIG
    VOLUME = args.VOLUME
    ENERGY = args.ENERGY

    ParaIn = esutils.read_input(INPUT)
    BaseName = ParaIn['BaseName']
    EnergyCode = ParaIn['EnergyCode'].lower()
    CrystalType = ParaIn['CrystalType'].lower()
    if ENERGY is None:
        if os.path.exists('EEnergy.txt'):
            ENERGY = 'EEnergy.txt'
        elif os.path.exists(os.path.join(BaseName, 'EEnergy.txt')):
            ENERGY = os.path.join(BaseName, 'EEnergy.txt')
        else:
            raise FileNotFoundError('The energy file not exists in both current folder or the BaseName folder, please specify it by -e parameter.')
    if VOLUME is None:
        if ParaIn['FlagSE'].lower() == 's':
            VOLUME = 1.0
        else:
            eglue = importlib.import_module('elastic3rd.energy.{}'.format(EnergyCode))
            if os.path.exists(os.path.join(BaseName, 'Mode0')):
                os.chdir(os.path.join(BaseName, 'Mode0'))
                BaseVec = eglue.get_base_vec(BaseName)
                os.chdir('../../')
            else:
                try:
                    BaseVec = eglue.get_base_vec(BaseName)
                except Exception as e:
                    raise FileNotFoundError('The structure file not found, pleas provide the strucutre file in current folder or provide the volume by -v parameter.')
            VOLUME = crylat.print_lattice(BaseVec)

    if STRAINMODE:
        #Not None
        (C2, C3) = espost.post_mode(V0=VOLUME, Flag_Fig=FIG, Flag_Ord=N, EEnergy=ENERGY, INPUT=INPUT, STRAINMODE=STRAINMODE)
    else:
        (C2, C3) = espost.post(V0=VOLUME, Flag_Fig=FIG, Flag_Ord=N, EEnergy=ENERGY, INPUT=INPUT)

    essym.print_cijk(CrystalType, 2)
    print(C2)
    essym.print_cijk(CrystalType, 3)
    print(C3)

def run_e3rd():
    print("ELASTIC3RD version: " + __version__)
    print("Copyright \u00a9 FGMS @ HIT\n")

    parser = argparse.ArgumentParser(description='Run elastic3rd.')
    
    subparsers = parser.add_subparsers()

    #SUB-PROCESS: run
    prun = subparsers.add_parser("run", help="Run elastic3rd.")
    prun.add_argument("-i", "--input", dest="INPUT", type=str, default="INPUT",
                      help="The input file for Elastic3rd. Default: 'INPUT'\n")
    prun.set_defaults(func=run)

    ppost = subparsers.add_parser('post', help='post process of elastic3rd.')
    ppost.add_argument('-n', '--n', dest='N', type=int, default=3,
                      help='The order of elastic constant taken into consideration in the post process. Default: 3\n')
    ppost.add_argument("-i", "--input", dest="INPUT", type=str, default="INPUT",
                      help="The input file for Elastic3rd. Default: 'INPUT'\n")
    ppost.add_argument("-sm", "--strainmode", dest="STRAINMODE", type=str,
                      help="The strainmode file for Elastic3rd, None means using the code to generate it. Default: None\n")
    ppost.add_argument("-fig", "--fig", dest="FIG", type=int, default=0,
                      help="Show the fitting fig. Default: 0\n")
    ppost.add_argument("-v", "--volume", dest="VOLUME", type=float,
                      help="The volume of initial structure. Default: None\n")
    ppost.add_argument("-e", "--energy", dest="ENERGY", type=str,
                      help="The file of the energy result. If None, it will try to find in current folder or in the BaseName folder. Default: None\n")
    ppost.set_defaults(func=post)

    args = parser.parse_args()

    try:
        a = getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)

if __name__ == '__main__':
    run_e3rd()