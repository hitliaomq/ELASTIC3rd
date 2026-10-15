#!python
#
import os, re
from monty.serialization import loadfn, dumpfn
import elastic3rd.post.post2 as espost2
import elastic3rd.post.post as espost
from elastic3rd.esutils import gen_strain_list, read_strainmode, read_input, write_input
import elastic3rd.symmetry.symmetry_stress as symstr

def post_elastic3rd(filenameout, strainmode='STRAINMODE', INPUT='INPUT.e3rd', EEnergy='EEnergy.txt'):
    StrainIn = read_strainmode(STRAINMODE=strainmode)
    ParaIn = read_input(INPUT=INPUT)
    CrystalType = ParaIn['CrystalType']
    Ord = ParaIn['Ord']
    flag_se = ParaIn['FlagSe']
    StrainList = gen_strain_list(ParaIn)
    Cij_mode, coef_e, StrainMode = symstr.CoefForSingleMode(CrystalType, Ord, StrainIn)
    E = espost.read_e(EEnergy=EEnergy)
    V0 = 1.
    coef_fit = espost.get_coef(StrainList/100., E, V0, flag_se, 3)
    C = espost2.get_cij(coef_fit, coef_e, flag_se, Ord=Ord)
    dumpfn(C, filenameout, indent=4)

def interleave_lists(list1, list2):
    combined = [item for pair in zip(list1, list2) for item in pair]
    return combined + list1[len(list2):] + list2[len(list1):]

cwd = os.getcwd()
lats = ['bcc', 'fcc', 'hcp']
pure_file = 'PureEle.json'
puredb = loadfn(pure_file)
for lat in lats:
    dirs = os.listdir(lat)
    results = []
    for diri in dirs:
        fulldir = os.path.join(lat, diri)
        if os.path.isdir(fulldir):
            fullname = os.path.join(fulldir, '{}-{}-elastic.json'.format(lat, diri))
            if os.path.exists(fullname):
                if lat == 'hcp':
                    path_cwd = os.path.join(cwd, lat, diri, diri)
                    STRAINMODE = os.path.join(path_cwd, 'STRAINMODE')
                    INPUT = os.path.join(path_cwd, 'INPUT.e3rd')
                    ParaIn = read_input(INPUT=INPUT)
                    ParaIn.update({'CrystalType': 'H1'})
                    write_input(ParaIn, FileName=INPUT)
                    EEnergy = os.path.join(path_cwd, 'EEnergy.txt')
                    post_elastic3rd(fullname, strainmode=STRAINMODE, INPUT=INPUT, EEnergy=EEnergy)
                data = loadfn(fullname)[0]
                data_dft = puredb[diri]['CIJ_{}'.format(lat)]
                resulti = [diri] + interleave_lists(data_dft, data)
                results.append(resulti)
    filename_result = '{}-result.txt'.format(lat)
    with open(filename_result, 'w+') as fid:
        for linei in results:
            linei_lst = [str(item) for item in linei]
            linei_str = '\t'.join(linei_lst)
            fid.write(linei_str)
            fid.write('\n')