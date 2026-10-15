#!python
#
#
from elastic3rd.elastic import elastic3
from elastic3rd.esutils import write_input
import numpy as np
import os, re
from monty.serialization import loadfn, dumpfn
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

def create_folder(fld):
    if os.path.exists(fld):
        pass
    else:
        os.mkdir(fld)

def build_structure(latparam, lat, ele):
    lat = lat.lower()
    if lat == 'bcc':
        sg = 229
        coords = [[0., 0., 0.]]
        lattice = Lattice.cubic(latparam)
        strainmode = [[1., 0., 0., 2., 0., 0.]]
    elif lat == 'fcc':
        sg = 225
        coords = [[0., 0., 0.]]
        lattice = Lattice.cubic(latparam)
        strainmode = [[1., 0., 0., 2., 0., 0.]]
    elif lat == 'hcp':
        sg = 194
        coords = [[1./3., 2./3., 0.25]]
        lattice = Lattice.hexagonal(a=latparam[0], c=latparam[1]*latparam[0])
        strainmode = [[1., 0., 0., 0., 0., 0.], [0., 0., 1., 2., 0., 0.]]
    struct = Structure.from_spacegroup(sg=sg, lattice=lattice, species=[ele], coords=coords)
    #struct.to(filename='POSCAR')
    return (struct, strainmode)


pure_file = 'PureEle.json'
puredb = loadfn(pure_file)
skip_ele = ['Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'S', 'I', 'P', 'C', 'Si']
energycode = 'sevenn'
settings = {'model': 'MatterSim-v1.0.0-5M'}

ParaIn = {'Continue': 1, 'EnergyCode' : energycode, 'EnergyRun' : 1, 'CrystalType': 'C1', 'Ord': 2, 
        'MaxStrain' : 1.0, 'NP' : 1, 'STEPS' : 9, 'BaseName' : 'Si', 'CalMode': 'S', 'FlagSE': 'S', 
        'GenSM' : 0, 'CoefMethod': 'fit', 'Acc': 4, 'Geom': 1, 'AddPert': 0, 'PertMethod': 'random',
        'StableEst': 1, 'StrainType': 'Uniform'}

lats = ['bcc', 'fcc', 'hcp']
cwd = os.getcwd()
for lat in lats:
    create_folder(lat)
    os.chdir(lat)
    results = {}
    for ele, eledb in puredb.items():
        if ele in skip_ele:
            continue
        else:
            create_folder(ele)
            os.chdir(ele)
            if os.path.exists('{}-{}-elastic.json'.format(lat, ele)):
                C = loadfn('{}-{}-elastic.json'.format(lat, ele))
            else:
                latparam = eledb['a_{}'.format(lat)]
                (struct, strainmode) = build_structure(latparam, lat, ele)
                ParaIn.update({'BaseName': ele})
                if lat in ['bcc', 'fcc', 'dia']:
                    ParaIn.update({'CrystalType': 'C1'})
                elif lat in ['hcp']:
                    ParaIn.update({'CrystalType': 'H1'})
                write_input(ParaIn, FileName='INPUT.e3rd')
                struct.to(filename='POSCAR')
                np.savetxt('STRAINMODE', strainmode)
                if settings:
                    dumpfn(settings, 'settings.yaml')
                C = elastic3()
                dumpfn(C, '{}-{}-elastic.json'.format(lat, ele), indent=4)
            results.update({ele: C})
            os.chdir('..')
    dumpfn(results, '{}-elastic.json'.format(lat), indent=4)
    os.chdir('..')

