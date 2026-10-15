# -*- coding: utf-8 -*-
# This is the script run elastic3rd in the terminal
# Simply run     elastic3rd
import sys
import os
import json
import argparse
import importlib
import elastic3rd
import elastic3rd.elastic
from elastic3rd import __version__
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
    if STRAINMODE is None:
        STRAINMODE = 'STRAINMODE'
    if ~os.path.exists(ENERGY):
        ENERGY = None
    if ENERGY is None:
        if os.path.exists('EEnergy.txt'):
            ENERGY = 'EEnergy.txt'
        elif os.path.exists(os.path.join(BaseName, 'EEnergy.txt')):
            ENERGY = os.path.join(BaseName, 'EEnergy.txt')
        else:
            raise ValueError('The energy file not exists in both current folder or the BaseName folder, please specify it by -e parameter.')
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
                    raise ValueError('The structure file not found, pleas provide the strucutre file in current folder or provide the volume by -v parameter.')
            VOLUME = crylat.print_lattice(BaseVec)

    C = espost.post_mode(V0=VOLUME, Flag_Sep=0, Flag_Fig=FIG, Flag_Data=0, Flag_Ord=N, EEnergy=ENERGY, 
            INPUT=INPUT, STRAINMODE=STRAINMODE, FigName='diamond-3rd', DataFileName='diamond-3rd')
    #if STRAINMODE:
    #    #Not None
    #    (C2, C3) = espost.post_mode(V0=VOLUME, Flag_Fig=FIG, Flag_Ord=N, EEnergy=ENERGY, INPUT=INPUT, STRAINMODE=STRAINMODE)
    #else:
    #    (C2, C3) = espost.post(V0=VOLUME, Flag_Fig=FIG, Flag_Ord=N, EEnergy=ENERGY, INPUT=INPUT)

    essym.print_cijk(CrystalType, 2)
    print(C.C2)
    essym.print_cijk(CrystalType, 3)
    print(C.C3)

def nd(args):
    import elastic3rd.tools.sep_nd as esnd
    ParaIn = esutils.read_input(args.INPUT)
    # Align the energy module with the given INPUT file (module names are lowercase)
    esnd.eglue = importlib.import_module('elastic3rd.energy.{}'.format(ParaIn['EnergyCode'].lower()))

    BaseName = ParaIn['BaseName']
    CrystalType = ParaIn['CrystalType']
    Ord = args.N if args.N is not None else ParaIn['Ord']
    NP = args.NP if args.NP is not None else ParaIn['NP']
    acc = args.ACC if args.ACC is not None else ParaIn['Acc']
    StrainAmp = args.STRAIN if args.STRAIN is not None else [ParaIn['MaxStrain']]

    if Ord not in (2, 3):
        raise ValueError('The separate numerical differentiation method supports 2nd and 3rd order elastic constants.')
    if acc not in (2, 4):
        raise ValueError('The accuracy of the numerical differentiation should be 2 or 4.')

    esnd.calc_ec(StrainAmp, CrystalType=CrystalType, Ord=Ord, BaseName=BaseName, NP=NP, acc=acc)

    ECfile = os.path.join('{}-E3rd-ND'.format(BaseName), 'EC_ave.json')
    if os.path.exists(ECfile):
        with open(ECfile) as fopen:
            EC = json.load(fopen)
        print('\n=======Averaged elastic constants (GPa)=======')
        for key in sorted(EC, key=lambda k: int(k)):
            print('{}  {}'.format(key, EC[key]))

def run_e3rd():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(errors='backslashreplace')
    print("ELASTIC3RD version: " + __version__)
    print("Copyright \u00a9 FGMS@HIT & JUST\n")

    parser = argparse.ArgumentParser(description='Run elastic3rd.')
    
    subparsers = parser.add_subparsers()

    #SUB-PROCESS: run
    prun = subparsers.add_parser("run", help="Run elastic3rd.")
    prun.add_argument("-i", "--input", dest="INPUT", type=str, default="INPUT.e3rd",
                      help="The input file for Elastic3rd. Default: 'INPUT.e3rd'\n")
    prun.set_defaults(func=run)

    ppost = subparsers.add_parser('post', help='post process of elastic3rd.')
    ppost.add_argument('-n', '--n', dest='N', type=int, default=3,
                      help='The order of elastic constant taken into consideration in the post process. Default: 3\n')
    ppost.add_argument("-i", "--input", dest="INPUT", type=str, default="INPUT.e3rd",
                      help="The input file for Elastic3rd. Default: 'INPUT.e3rd'\n")
    ppost.add_argument("-sm", "--strainmode", dest="STRAINMODE", type=str, default="STRAINMODE",
                      help="The filename of strainmodes. Default: 'STRAINMODE'\n")
    ppost.add_argument("-fig", "--fig", dest="FIG", type=int, default=0,
                      help="If show the figure or not. Default: 0\n")
    ppost.add_argument("-v", "--volume", dest="VOLUME", type=float, default=1.,
                      help="The volume for strain-energy method.")
    ppost.add_argument("-e", "--energy", dest="ENERGY", type=str, default="EEnergy.txt",
                      help="The energy file for post")
    ppost.set_defaults(func=post)

    pnd = subparsers.add_parser('nd', help='Calculate elastic constants by the separate numerical '
                                'differentiation method (PRL 121, 216001), each Cij component from '
                                'dedicated strain-stress finite differences.')
    pnd.add_argument("-i", "--input", dest="INPUT", type=str, default="INPUT.e3rd",
                      help="The input file for Elastic3rd. Default: 'INPUT.e3rd'\n")
    pnd.add_argument('-n', '--n', dest='N', type=int, default=None,
                      help='The order of elastic constants, 2 or 3. Default: the Ord in the input file\n')
    pnd.add_argument('-s', '--strain', dest='STRAIN', type=float, nargs='+', default=None,
                      help='The strain amplitudes, in unit of %%. Multiple values give a mean over '
                           'them. Default: the MaxStrain in the input file\n')
    pnd.add_argument('-a', '--acc', dest='ACC', type=int, default=None,
                      help='The accuracy of the differentiation, 2 or 4. Default: the Acc in the input file\n')
    pnd.add_argument('--np', dest='NP', type=int, default=None,
                      help='The number of cores used in the energy/stress calculation. Default: the NP in the input file\n')
    pnd.set_defaults(func=nd)

    args = parser.parse_args()

    try:
        a = getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)

if __name__ == '__main__':
    run_e3rd()