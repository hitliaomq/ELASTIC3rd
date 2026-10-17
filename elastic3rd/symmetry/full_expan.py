#!python
#
from elastic3rd.symmetry.symmetry import Num2IJK, cijk2num, num2cijk
from elastic3rd.symmetry.symdata import coef_ind, coef_crystal
from elastic3rd.esutils import sort_x_by_y
import fractions as frct
import numpy as np
import json
import yaml
import math
import logging
import copy
#import latexify


class FullExpan(object):
    """docstring for FullExpan
        {11: MultiEta, 12: MultiEta, 44: MultiEta}
    """
    def __init__(self, CrystalType='C1', Ord=2, values=None, latex=False):
        super(FullExpan, self).__init__()
        self.CrystalType = CrystalType
        self.Ord = Ord
        self.values = {}
        self.latex = latex
        pre_coef_meta = math.factorial(self.Ord)
        self.coef = 1. / pre_coef_meta
        if values is None:
            self.values = self.full_expan()
        elif isinstance(values, dict):
            for keyi, valuei in values.items():
                if isinstance(valuei, dict):
                    self.values[keyi] = MultiEta.from_dict(valuei)
                elif isinstance(valuei, MultiEta):
                    self.values[keyi] = valuei
                elif isinstance(valuei, list) or isinstance(valuei, Eta):
                    self.values[keyi] = MultiEta.from_etas(valuei)
                elif isinstance(valuei, (int, float)):
                    #self.values[keyi] = valuei
                    if float(valuei) == 0.:
                        pass
                else:
                    logging.warning('The values should be dict, the code generate it by full_expan.')
                    self.values = self.full_expan()
        else:
            logging.warning('The values should be dict, the code generate it by full_expan.')
            self.values = self.full_expan()

    def __str__(self, latex=False):
        pre_coef_meta = math.factorial(self.Ord)
        strings = ''
        count_f = 0
        for keyi, valuei in self.values.items():
            if latex:
                stringi = '\\frac{{1}}{{{}}}C_{{{}}}({})'.format(pre_coef_meta, keyi, valuei.__str__(latex=latex)) #valuei.string)
            else:
                stringi = '1/{} * C{} * ({})'.format(pre_coef_meta, keyi, valuei.__str__(latex=latex)) #valuei.string)
            if count_f == 0:
                strings = strings + stringi
            else:
                strings = strings + ' + ' + stringi
            count_f = count_f + 1
        return strings

    def full_expan(self):
        CrystalType = self.CrystalType
        Ord = self.Ord
        #TODO: add the coefficinets
        pre_coef_meta = 1. / float(math.factorial(Ord))
        result = {}
        CoefInd = coef_ind(CrystalType=CrystalType, Ord=Ord)
        CoefCoef = coef_crystal(CrystalType=CrystalType, Ord=Ord)
        for i_tennum in range(3 ** (2 * Ord)):
            k = i_tennum + 1
            indexk = Num2IJK(Num=k, N=3, Ord=2*Ord)

            etai = 1.
            cij_voigt_list = []
            for i_eta in range(Ord):
                eta_ind_list = [indexk[2*i_eta], indexk[2*i_eta + 1]]
                eta_ind = cijlist2ind(eta_ind_list)
                etai = Eta.from_dict({eta_ind: [1, 1]}) * etai

                cij_voigt_list.append(ten2voigt(eta_ind))
            
            cij_num = cijk2num(cijlist2ind(sorted(cij_voigt_list)), Ord=Ord)
            cij_num_true = CoefInd[cij_num[0]]
            coefcoef_true = CoefCoef[cij_num[0]]
            if isinstance(cij_num_true, list):
                for i_num, numi in enumerate(cij_num_true):
                    cij_voigt = int(num2cijk(num=numi, Ord=Ord))
                    if cij_voigt in result:
                        result[cij_voigt] = result[cij_voigt] + copy.deepcopy(etai) * coefcoef_true[i_num]
                    else:
                        result[cij_voigt] = copy.deepcopy(etai) * coefcoef_true[i_num]
            else:
                if cij_num_true == 0:
                    pass
                elif isinstance(cij_num_true, int):
                    cij_voigt = int(num2cijk(num=cij_num_true, Ord=Ord))
                    if cij_voigt in result:
                        result[cij_voigt] = result[cij_voigt] + copy.deepcopy(etai) * coefcoef_true
                    else:
                        result[cij_voigt] = copy.deepcopy(etai) * coefcoef_true
                else:
                    raise TypeError
        return result
    
    def derivation(self, d_eta_ind=11, CrystalType='C1', Ord=2):
        dmeta = {}
        for keyi in self.values:
            metai = self.values[keyi]
            #print('C{} : {}'.format(keyi, self.values[keyi]))
            #print(metai)
            #print(type(metai))
            if isinstance(metai, Eta):
                metai = MultiEta.from_etas(metai)
            dmetai = metai.derivation(d_eta_ind=d_eta_ind)
            dmeta[keyi] = dmetai
        dmeta = FullExpan(CrystalType=CrystalType, Ord=Ord, values=dmeta)
        return dmeta

    def to_dict(self):
        values = self.values
        value_dict = {}
        for keyi in values:
            value_dict[keyi] = values[keyi].to_dict()
        return value_dict

    @classmethod
    def from_dict(cls, value_dict):
        if isinstance(value_dict, dict):
            CrystalType = value_dict.pop('CrystalType', 'C1')
            Ord = value_dict.pop('Ord', 2)
            return FullExpan(CrystalType=CrystalType, Ord=Ord, values=value_dict)
        else:
            raise TypeError('Please provide dict.')

    @classmethod
    def from_file(cls, filename):
        with open(filename, 'r') as fid:
            if filename.lower().endswith('.yaml'):
                value_dict = yaml.load(fid)
            elif filename.lower().endswith('.json'):
                value_dict = json.load(fid)
            else:
                raise ValueError
        return FullExpan.from_dict(value_dict=value_dict)

    def write2file(self, filename=None, filetype='yaml', **kwargs):
        value_dict = self.to_dict()
        value_dict['CrystalType'] = self.CrystalType
        value_dict['Ord'] = self.Ord
        if filename is None:
            filename = '{}-Ord{}.{}'.format(self.CrystalType, self.Ord, filetype)
        with open(filename, 'w+') as fid:
            if filetype.lower().startswith('j'):
                json.dump(value_dict, fid, **kwargs)
            elif filetype.lower().startswith('y'):
                yaml.dump(value_dict, fid, **kwargs)
            else:
                raise ValueError

    def to_beta(self):
        val = self.values
        beta = 0
        for keyi, vali in val.items():
            vali_beta = vali.to_beta()
            val_dict = vali_beta.values
            betadict = {}
            for key_beta, val_beta in val_dict.items():
                val_tmp = []
                for i, coefi in enumerate(val_beta):
                    if float(coefi) == 0.:
                        val_tmp.append(0)
                    else:
                        val_tmp.append({keyi: coefi})
                betadict[key_beta] = val_tmp
            betai = Beta(betadict=betadict)
            beta = beta + betai
        return beta

    def fullexpression(self, tind=11, withPrecoef=False, Ord=None):
        '''
        Parameters
            withPrecoef: bool
                If False, using tylor to replace (1+2*eta)^(-1/2)
                If True, keep the coefficient
        '''
        etalist0 = [1., -1., 3./2., -5./2., 35./8., -63./8., 231./16., -429./16., 6435./128., -12155./128.,
                    46189./256., -88179./256., 676039./1024., -1300075./1024., 5014575./2048.]
        if Ord is None:
            Ord = self.Ord + 1
        etalist0 = etalist0[0:Ord]
        precoef_eta = SingleEta(etalist=etalist0, latex=self.latex)
        betadict0 = {0: [1], 11: [-1], 22: [-1], 33:[-1]}
        precoef_beta = Beta(betadict=betadict0, latex=self.latex)
        if withPrecoef:
            prebeta = precoef_beta
        else:
            prebeta = precoef_eta * precoef_beta

        tind_list = [int(i) for i in str(tind)]
        result = 0
        for i in range(1, 4):
            beta_ind1 = int(10 * i + tind_list[0])
            if i == tind_list[0]:
                beta1 = Beta(betadict={0: [1], beta_ind1: [1]})
            else:
                beta1 = Beta(betadict={beta_ind1: [1]})
            for j in range(1, 4):
                beta_ind2 = int(10 * j + tind_list[1])
                if j == tind_list[1]:
                    beta2 = Beta(betadict={0: [1], beta_ind2: [1]})
                else:
                    beta2 = Beta(betadict={beta_ind2: [1]})
                beta_tmp = beta1 * beta2

                eta_ind = int(''.join([str(i), str(j)]))
                d_eta_ij = self.derivation(d_eta_ind=eta_ind, CrystalType=self.CrystalType, Ord=self.Ord)
                
                d_eta_ij_beta = d_eta_ij.to_beta()

                if i == 1 and j == 1:
                    P = Beta(betadict={0: d_eta_ij_beta.values[0]}) * (-1)
                    #P = {0: d_eta_ij.values[0]}
                result0 = beta_tmp * d_eta_ij_beta
                result = result + result0
        result = prebeta * result
        #P = prebeta * P
        pre_coef_meta = 1. / math.factorial(self.Ord)
        result = result * float(pre_coef_meta)
        P = P * float(pre_coef_meta)
        print('P : {}\n'.format(P))
        result_val = result.values
        for keyi, vali in result_val.items():
            if keyi != 0:
                print('{} : {}\n'.format(keyi, Beta(betadict={keyi: vali})))
        #print(type(result))
        return result

        #print(prebeta)
        #print(type(prebeta))

def cijlist2ind(cijlist):
    return int(''.join([str(i) for i in cijlist]))

def ten2voigt(ten_ind):
    voigtmap = {11: 1, 12: 6, 21: 6, 13: 5, 31: 5, 22: 2, 23: 4, 32: 4, 33: 3}
    return voigtmap[ten_ind]

def voigt2ten(voigt_ind):
    tenmap = {1: 11, 2: 22, 3: 33, 4: 23, 5: 13, 6: 12}
    return tenmap[voigt_ind]


class MultiEta(object):
    """docstring for MultiEta
        [{11: [1, 1], 12: [1, 2]}, {11: [2, 2], 13: [1, 1]}]
        {'coef': [1, 2],
         'eta':[{11: 1, 12: 1}, {11: 2, 13: 1}]}
        eta_11 * eta_12^2 + 2 * eta_11^2 * eta_13

    """
    def __init__(self, eta_inds, pre_coefs, exponents, latex=False):
        super(MultiEta, self).__init__()
        self.eta_inds = eta_inds
        self.pre_coefs = pre_coefs
        self.exponents = exponents
        self.latex = latex
        if not ((len(eta_inds) == len(pre_coefs)) and (len(eta_inds) == len(exponents))):
            raise ValueError('The length of eta_inds, pre_coefs and exponents should match.')        
        self.values = self.to_dict()
        self.string = self.__str__(latex=self.latex)

    def __str__(self, latex=False):
        count_p = 0     #For plus
        for i, valuei in enumerate(self.values['eta']):
            count_m = 0     #For multipy
            precoef_m = float(self.values['coef'][i])
            for keyi in valuei:
                etai = Eta.from_dict({keyi: [1., valuei[keyi]]})
                #print(etai)
                etai_str = etai.__str__(latex=latex) #etai.string
                if not isinstance(etai_str, (int, float)):
                    if count_m == 0:
                        string_m = etai_str
                    else:
                        if latex:
                            string_m = '{} {}'.format(string_m, etai_str)
                        else:
                            string_m = '{} * {}'.format(string_m, etai_str)
                    count_m = count_m + 1
            if count_m == 0:
                   string_m = etai_str

            if isinstance(string_m, (float, int)):
                if str(string_m) == '0.0':
                    continue
                else:
                    string_m = precoef_m * float(string_m)
            elif precoef_m != 1.:
                if latex:
                    string_m = '{} {}'.format(precoef_m, string_m)
                else:
                    string_m = '{} * {}'.format(precoef_m, string_m)
            if count_p == 0:
                string_p = string_m
            else:
                if precoef_m < 0:
                    string_p = '{} {}'.format(string_p, string_m)
                else:
                    string_p = '{} + {}'.format(string_p, string_m)
            count_p = count_p + 1
        return string_p

    def to_etas(self):
        Etas = []
        for i, valuei in enumerate(self.values['eta']):
            Etai = []
            count_m = 0
            for keyi in valuei:
                if count_m == 0:
                    etai = Eta.from_dict({keyi: [self.values['coef'][i], valuei[keyi]]})
                else:
                    etai = Eta.from_dict({keyi: [1., valuei[keyi]]})
                Etai.append(etai)
                count_m = count_m + 1
            Etas.append(Etai)   
        return Etas

    def to_beta(self):
        val = self.to_dict()
        coefs = val['coef']
        etas = val['eta']
        beta = 0
        for i, etai in enumerate(etas):
            betai = coefs[i]
            for keyi, vali in etai.items():
                betai = betai * SingleBeta(beta_ind=keyi) ** vali
            beta = beta + betai
        return beta

    def to_dict(self):
        eta_inds = self.eta_inds
        pre_coefs = self.pre_coefs
        exponents = self.exponents
        
        values = {'coef': pre_coefs, 'eta': []}
        for i, eta_indi in enumerate(eta_inds):
            valuei = {}
            for j, eta_indij in enumerate(eta_indi):
                if eta_indij in valuei:
                    valuei[eta_indij] = valuei[eta_indij] + exponents[i][j]
                else:
                    valuei[eta_indij] = exponents[i][j]
            values['eta'].append(valuei)
        return values

    @classmethod
    def from_dict(cls, etasdict):
        # etalist format: [{eta_ind1: [pre_coef1, exponent1], eta_ind2: [...]}, {...}]
        eta_inds = [[keyi for keyi in etai] for etai in etasdict['eta']]
        pre_coefs = etasdict['coef']
        exponents = [[etai[keyi] for keyi in etai] for etai in etasdict['eta']]
        return MultiEta(eta_inds=eta_inds, pre_coefs=pre_coefs, exponents=exponents)

    @classmethod
    def from_etas(cls, etas):
        #etas format: [[Eta1, Eta2], [Eta3, Eta4], Eta5, ...]
        if isinstance(etas, Eta):
            eta_inds = [[etas.eta_ind]]
            pre_coefs = [etas.pre_coef]
            exponents = [[etas.exponent]]
        elif isinstance(etas, list):
            eta_inds, pre_coefs, exponents = [], [], []
            for etai in etas:
                if isinstance(etai, Eta):
                    eta_inds.append([etai.eta_ind])
                    pre_coefs.append(etai.pre_coef)
                    exponents.append([etai.exponent])
                elif isinstance(etai, (list, tuple)):
                    eta_indsi, exponentsi = [], []
                    pre_coefsi = 1.
                    for etaij in etai:
                        if isinstance(etaij, Eta):
                            eta_indsi.append(etaij.eta_ind)
                            exponentsi.append(etaij.exponent)
                            pre_coefsi = pre_coefsi * etaij.pre_coef
                        else:
                            raise TypeError('etas must be Eta class, or list of Eta class.')
                    eta_inds.append(eta_indsi)
                    pre_coefs.append(pre_coefsi)
                    exponents.append(exponentsi)
        else:
            raise TypeError('etas must be Eta class, or list of Eta class.')
        return MultiEta(eta_inds=eta_inds, pre_coefs=pre_coefs, exponents=exponents)

    @staticmethod
    def sort_etas_by_ind(etas):
        '''
        Sort etas according to eta_ind
        Note: Only work for [eta1, eta2, ...]
        '''
        inds = [etai.eta_ind for etai in etas]                        
        return sort_x_by_y(etas, inds)

    def __add__(self, newMultiEta):
        if isinstance(newMultiEta, Eta):
            return newMultiEta + self
        if isinstance(newMultiEta, MultiEta):
            pass
        elif isinstance(newMultiEta, (int, float)):
            return self + Eta.from_dict({11: [newMultiEta, 0]})
        else:
            raise TypeError('The input for * shoulb be int or float or MultiEta class.')

        etas_new = newMultiEta.to_etas()
        etas_self = self.to_etas()
        for i_new, etai_new in enumerate(etas_new):
            flag_add = False
            for i_self, etai_self in enumerate(etas_self):
                if isinstance(etai_new, Eta) and isinstance(etai_self, Eta):
                    eta_tmp = etai_new + etai_self
                    if isinstance(eta_tmp, Eta):
                        etas_self[i_self] = eta_tmp
                        flag_add = True
                        break
                elif isinstance(etai_new, list) and isinstance(etai_self, list):
                    flag_multipy = True
                    if len(etai_new) == len(etai_self):
                        etai_new = self.sort_etas_by_ind(etai_new)
                        meta_new = MultiEta.from_etas([etai_new])

                        etai_self = self.sort_etas_by_ind(etai_self)
                        meta_self = MultiEta.from_etas([etai_self])
                        for i_etai_new, etai_newi in enumerate(etai_new):
                            if etai_newi.eta_ind == etai_self[i_etai_new].eta_ind:
                                if etai_newi.exponent != etai_self[i_etai_new].exponent:
                                    flag_multipy = False
                                    break
                            else:
                                flag_multipy = False
                                break
                    else:
                        flag_multipy = False
                    if flag_multipy:
                        meta_tmp = meta_self * (1./meta_self.pre_coefs[0]) * (meta_self.pre_coefs[0] + meta_new.pre_coefs[0])
                        etas_self[i_self] = meta_tmp.to_etas()[0]
                        flag_add = True
                else:
                    raise TypeError
            if not flag_add:
                etas_self.append(etai_new)
        return MultiEta.from_etas(etas_self)

    def __radd__(self, newMultiEta):
        return self.__add__(self, newMultiEta)

    def __mul__(self, newMultiEta):
        if isinstance(newMultiEta, (int, float)):
            pre_coefs = [float(newMultiEta) * pre_coefi for pre_coefi in self.pre_coefs]
            return MultiEta(eta_inds=self.eta_inds, pre_coefs=pre_coefs, exponents=self.exponents)
        elif isinstance(newMultiEta, Eta):
            return self * MultiEta.from_etas([newMultiEta])
        elif isinstance(newMultiEta, MultiEta):
            pass
        else:
            raise TypeError('The input for * shoulb be int or float or MultiEta class.')

        etas_new = newMultiEta.to_etas()
        etas_self = self.to_etas()
        count_m = 0
        for i_new, etai_new in enumerate(etas_new):
            flag_add = False
            for i_self, etai_self in enumerate(etas_self):
                etai_self_tmp = copy.deepcopy(etai_self)
                for etai_n in etai_new:
                    flag_m = True
                    for i_s, etai_s in enumerate(etai_self_tmp):
                        eta_tmp = etai_n * etai_s
                        if isinstance(eta_tmp, Eta):
                            etai_self_tmp[i_s] = eta_tmp
                            flag_m = False
                            break
                    if flag_m:
                        etai_self_tmp.append(etai_n)
                if count_m == 0:
                    resultMeta = MultiEta.from_etas([etai_self_tmp])
                else:
                    resultMeta = resultMeta + MultiEta.from_etas([etai_self_tmp])
                count_m = count_m + 1
        return resultMeta

    def __rmul__(self, newMultiEta):
        return self.__mul__(newMultiEta)

    def derivation1(self, d_eta_ind=11):
        #derivation for one term (no add)
        etas_dict = self.to_dict()
        etas_inds = etas_dict['eta'][0]
        if d_eta_ind in etas_inds:
            etas_dict['coef'][0] = etas_dict['coef'][0] * etas_dict['eta'][0][d_eta_ind]
            etas_dict['eta'][0][d_eta_ind] = etas_dict['eta'][0][d_eta_ind] - 1
            return MultiEta.from_dict(etas_dict)
        else:
            return 0

    def derivation(self, d_eta_ind=11):
        etas_list = self.to_etas()
        dmeta = 0
        count_d = 0
        for etai in etas_list:
            metai = MultiEta.from_etas([etai])
            dmetai = metai.derivation1(d_eta_ind=d_eta_ind)
            if isinstance(dmetai, MultiEta):
                if count_d == 0:
                    dmeta = copy.deepcopy(dmetai)
                else:
                    dmeta = dmeta + dmetai
                count_d = count_d + 1
        return dmeta


class Eta(object):
    """docstring for Eta
    The class of strain component
        It has the following format
            {11, [1, 2]} means eta_{11}^2
            11: is the strain component index
            1: is the pre-coefficient
            2: is the exponent
    """
    def __init__(self, eta_ind=11, pre_coef=1, exponent=1, isvoigt=False, latex=False):
        super(Eta, self).__init__()
        self.isvoigt = isvoigt
        self.latex = latex
        if isinstance(eta_ind, (int, float)):
            eta_ind = int(''.join([i for i in str(eta_ind)]))
        elif isinstance(eta_ind, list):
            eta_ind = [str(i) for i in eta_ind]
            eta_ind = int(''.join(eta_ind))
        if eta_ind < 7 and eta_ind > 0:
            if not isvoigt:
                eta_ind = voigt2ten(eta_ind)
                pre_coef = 0.5 * pre_coef
        elif eta_ind in [11, 12, 13, 21, 22, 23, 31, 32, 33]:
            if isvoigt:
                eta_ind = ten2voigt(eta_ind)
                pre_coef = 2. * pre_coef
        else:
            raise ValueError('The eta_ind must be the index of strain component in tensor format or voigt format')
        self.eta_ind = eta_ind
        self.pre_coef = float(pre_coef)
        self.exponent = exponent
        self.eta = self.to_dict()
        self.string = self.__str__(latex=self.latex)

    @classmethod
    def from_dict(cls, etadict):
        #TODO: dict in dict
        if len(etadict) > 1:
            logging.warning('The length of the dict is larger than 1, only one item is used.')
        for keyi, valuei in etadict.items():
            if isinstance(valuei, (int, float)):
                return Eta(eta_ind=keyi, pre_coef=valuei)
            elif isinstance(valuei, (list, tuple)):
                len_v = len(valuei)
                if len_v == 1:
                    return Eta(eta_ind=keyi, pre_coef=valuei[0])
                elif len_v == 2:
                    return Eta(eta_ind=keyi, pre_coef=valuei[0], exponent=valuei[1])
                elif len_v == 3:
                    return Eta(eta_ind=keyi, pre_coef=valuei[0], exponent=valuei[1], isvoigt=valuei[2])
                else:
                    raise ValueError
            else:
                raise TypeError

    def to_dict(self):
        return {self.eta_ind: [self.pre_coef, self.exponent]}

    def __str__(self, latex=False):
        if latex:
            str0 = '\\eta_{{{}}}'.format(self.eta_ind)
        else:
            str0 = '{}_{}'.format(chr(951), self.eta_ind)
        if self.pre_coef == 1.:
            if self.exponent == 1:
                return str0
            elif self.exponent == 0:
                return 1.
            else:
                if latex:
                    return '{}^{{{}}}'.format(str0, self.exponent)
                else:
                    return '{}^{}'.format(str0, self.exponent)
        else:
            if self.exponent == 1:
                if latex:
                    return '{} {}'.format(self.pre_coef, str0)
                else:
                    return '{}*{}'.format(self.pre_coef, str0)
            elif self.exponent == 0:
                return self.pre_coef
            else:
                if latex:
                    return '{} {}^{{{}}}'.format(self.pre_coef, str0, self.exponent)
                else:
                    return '{}*{}^{}'.format(self.pre_coef, str0, self.exponent)

    def to_beta(self):
        #eta: {eta_ind: [pre_coef, exponent]}
        #beta: {0: [0, 1], 11: [1, 2, {11: 1, 12: 2}]}
        beta = SingleBeta(beta_ind=self.eta_ind) ** self.exponent
        beta = beta * self.pre_coef
        return beta

    def derivation(self):
        #return {self.eta_ind: [self.pre_coef * self.exponent, self.exponent - 1]}
        return self.update(pre_coef=self.pre_coef * self.exponent, exponent=self.exponent-1)

    def voigt(self):
        return self.update(isvoigt=True)

    def full(self):
        return self.update(isvoigt=False)

    def sort_ind(self):
        eta_ind = int(''.join(sorted(str(self.eta_ind))))
        return self.update(eta_ind=eta_ind)

    def update(self, **kwargs):
        eta_ind = kwargs.get('eta_ind', self.eta_ind)
        pre_coef = kwargs.get('pre_coef', self.pre_coef)
        exponent = kwargs.get('exponent', self.exponent)
        isvoigt = kwargs.get('isvoigt', self.isvoigt)
        return Eta(eta_ind=eta_ind, pre_coef=pre_coef, exponent=exponent, isvoigt=isvoigt)

    def __add__(self, newEta):
        #if self.exponent != newEta.exponent:
        #    #If use correct, this case should not happen
        #   raise TypeError('The eta of this two eta is not the same.')
        if isinstance(newEta, Eta):
            pass
        elif isinstance(newEta, MultiEta):
            return newEta + MultiEta.from_etas([self])
        elif isinstance(newEta, (int, float)):
            return self + Eta.from_dict({11: [newEta, 0]})
        else:
            raise TypeError('The input data should be Eta class.')
        if self.eta_ind == newEta.eta_ind:
            return self.update(pre_coef=self.pre_coef + newEta.pre_coef)
        else:
            return MultiEta.from_etas([self, newEta])

    def __radd__(self, newEta):
        return self.__add__(newEta)
    
    def negtive(self):
        return self.update(pre_coef=-self.pre_coef)

    def __sub__(self, newEta):
        newEta = newEta.negtive()
        return self.__add__(newEta)

    def __rsub__(self, newEta):
        return self.__sub__(newEta)

    def __mul__(self, newEta):
        if isinstance(newEta, (int, float)):
            return self.update(pre_coef=self.pre_coef * float(newEta))
        elif isinstance(newEta, Eta):
            if self.eta_ind == newEta.eta_ind:
                return self.update(pre_coef=self.pre_coef * newEta.pre_coef, exponent=self.exponent + newEta.exponent)
            else:
                return MultiEta.from_etas([[self, newEta]])
        elif isinstance(newEta, MultiEta):
            return newEta * self
        else:
            raise TypeError('Only suport numerical data of Eta class')

    def __rmul__(self, newEta):
        return self.__add__(newEta)

    def __div__(self, newEta):
        if self.eta_ind != newEta.eta_ind:
            raise TypeError('This two eta is not the same index.')
        return self.update(pre_coef=self.pre_coef / newEta.pre_coef, exponent=self.exponent - newEta.exponent)


class SingleEta(object):
    """docstring for SingleEta
    [k0, k1, k2, k3, ...] means [k0*eta^0, k1*eta^1, ...]
    for ki, it should be a dict or (int, float)
        if ki is int or float
            k0 is just the constants coefficients
        elif ki is dict
            ki = {0: constants, Cij: Cij_coef}
    """
    def __init__(self, etalist, latex=False):
        super(SingleEta, self).__init__()
        if isinstance(etalist, (int, float)):
            etalist = [etalist]
        self.values = etalist
        self.latex = latex
        self.string = self.__str__(latex=latex)
        self.n = len(etalist)
        for etai in etalist:
            if isinstance(etai, (int, float)) and float(etai) == 0.:
                self.n = self.n - 1

    def __str__(self, latex=False):
        etalist = self.values
        strings = ''
        count_eta = 0
        for i, etai in enumerate(etalist):
            if i == 0:
                str0 = '1.'
            else:
                if latex:
                    if i == 1:
                        str0 = '\\eta'
                    else:
                        str0 = '\\eta^{{{}}}'.format(i)
                else:
                    if i == 1:
                        str0 = '{}'.format(chr(951))
                    else:
                        str0 = '{}^{}'.format(chr(951), i)
            if isinstance(etai, (int, float)):
                if float(etai) == 0.:
                    continue
                if count_eta == 0:
                    if i == 0:
                        strings = str(etai)
                    else:
                        if float(etai) == 1.:
                            strings = '{}'.format(str0)
                        else:
                            if latex:
                                strings = '{}{}'.format(etai, str0)
                            else:
                                strings = '{} * {}'.format(etai, str0)
                else:
                    if i == 0:
                        strings = '{} + {}'.format(strings, etai)
                    elif float(etai) == 1.:
                        strings = '{} + {}'.format(strings, str0)
                    else:
                        if latex:
                            strings = '{} + {}{}'.format(strings, etai, str0)
                        else:
                            strings = '{} + {} * {}'.format(strings, etai, str0)
                count_eta = count_eta + 1
            elif isinstance(etai, dict):
                count_d = 0
                n_dict = len(etai)
                for eta_keyi, eta_vali in etai.items():
                    if eta_vali == 0:
                        n_dict = n_dict - 1
                        continue
                    if eta_keyi == 0:
                        #constants
                        string_etai = str(eta_vali)
                    else:
                        #Cij
                        if latex:
                            if float(eta_vali) == 1.:
                                string_etai = 'C_{{{}}}'.format(eta_keyi)
                            else:
                                string_etai = '{}C_{{{}}}'.format(eta_vali, eta_keyi)
                        else:
                            if float(eta_vali) == 1.:
                                string_etai = 'C{}'.format(eta_keyi)
                            else:
                                string_etai = '{} * C{}'.format(eta_vali, eta_keyi)
                    if count_d == 0:
                        if n_dict == 1:
                            string_etas = '{}'.format(string_etai)
                        else:
                            string_etas = '({}'.format(string_etai)
                    else:
                        string_etas = '{} + {}'.format(string_etas, string_etai)
                    count_d = count_d + 1
                if n_dict == 1:
                    string_etas = '{}'.format(string_etas)
                else:
                    string_etas = '{})'.format(string_etas)
                if count_eta == 0:
                    if i == 0:
                        strings= string_etas
                    else:
                        if latex:
                            strings= '{}{}'.format(string_etas, str0)
                        else:
                            strings= '{} * {}'.format(string_etas, str0)
                else:
                    if i == 0:
                        strings= '{} + {}'.format(strings, string_etas)
                    else:
                        if latex:
                            strings= '{} + {}{}'.format(strings, string_etas, str0)
                        else:
                            strings= '{} + {} * {}'.format(strings, string_etas, str0)
                count_eta = count_eta + 1
            else:
                raise TypeError('The etai should be int or float or dict.')
            count_eta = count_eta + 1
        return strings

    def __add__(self, newSingleEta):
        etalist_self = self.values
        if isinstance(newSingleEta, (int, float)):            
            etalist_self[0] = etalist_self[0] + newSingleEta
            return SingleEta(etalist=etalist_self, latex=self.latex)
        elif isinstance(newSingleEta, SingleEta):
            etalist_new = newSingleEta.values
            n_self = len(etalist_self)
            n_new = len(etalist_new)
            n_min = min(n_new, n_self)
            if n_self < n_new:
                etalist = copy.deepcopy(etalist_new)
            else:
                etalist = copy.deepcopy(etalist_self)
            for i in range(n_min):
                etai_self = etalist_self[i]
                etai_new = etalist_new[i]
                if isinstance(etai_self, (int, float)) and isinstance(etai_new, (int, float)):
                    etalist[i] = etai_self + etai_new
                elif isinstance(etai_self, (int, float)) and isinstance(etai_new, dict):
                    if 0 in etai_new:
                        etai_new[0] = etai_new[0] + etai_self
                    else:
                        etai_new[0] = etai_self
                    etalist[i] = etai_new
                elif isinstance(etai_self, dict) and isinstance(etai_new, (int, float)):
                    etalist[i][0] = etalist_self[i][0] + etai_new
                elif isinstance(etai_self, dict) and isinstance(etai_new, dict):
                    etai_self_tmp = copy.deepcopy(etai_self)
                    for key_new, val_new in etai_new.items():
                        if key_new in etai_self_tmp:
                            etai_self_tmp[key_new] = etai_self_tmp[key_new] + val_new
                        else:
                            etai_self_tmp[key_new] = val_new
                    etalist[i] = etai_self_tmp
                else:
                    raise TypeError
            return SingleEta(etalist=etalist, latex=self.latex)
        elif isinstance(newSingleEta, (SingleBeta, Beta)):
            return newSingleEta + self
        else:
            raise TypeError('The newSingleEta should be int or float or SingleEta.')
    
    def __radd__(self, newSingleEta):
        return self.__add__(newSingleEta)

    def __mul__(self, newSingleEta):
        val_self = copy.deepcopy(self.values)
        #print(val_self)
        n_self = len(val_self)
        try:
            n_new = len(newSingleEta.values)
        except Exception as e:
            n_new = 1
        result = [0 for i in range(n_self + n_new - 1)]
        result_se = SingleEta(etalist=result, latex=self.latex)
        if isinstance(newSingleEta, (int, float)):
            newSingleEta = SingleEta(etalist=[newSingleEta], latex=self.latex)
            return self * newSingleEta
        elif isinstance(newSingleEta, list):
            newSingleEta = SingleEta(etalist=newSingleEta, latex=self.latex)
            return self * newSingleEta
        elif isinstance(newSingleEta, Beta):
            return newSingleEta * self
        elif isinstance(newSingleEta, SingleEta):
            val_new = newSingleEta.values
            for i_self, vali_self in enumerate(val_self):
                for i_new, vali_new in enumerate(val_new):
                    i_result = i_self + i_new
                    if isinstance(vali_self, (int, float)) and isinstance(vali_new, (int, float)):
                        result_tmp = vali_self * vali_new
                        #result[i_result] = result[i_result] + vali_self * vali_new
                    elif isinstance(vali_self, (int, float)) and isinstance(vali_new, dict):
                        result_tmp = {keyi: vali*vali_self for keyi, vali in vali_new.items()}
                    elif isinstance(vali_self, dict) and isinstance(vali_new, (int, float)):
                        result_tmp = {keyi: vali*vali_new for keyi, vali in vali_self.items()}
                    elif isinstance(vali_self, dict) and isinstance(vali_new, dict):
                        logging.warning('* for both dicts skipped.')
                        continue
                    else:
                        raise TypeError('* in SingleEta only support for numerical by numerical and numerical by dict.')

                    result0 = copy.deepcopy(result)
                    result0[i_result] = result_tmp
                    result_SingleEta = SingleEta(etalist=result0, latex=self.latex)
                    result_se = result_se + result_SingleEta
        else:
            raise TypeError('* for SingleEta only support int or float or list or SingleEta.')
        return result_se

    def __rmul__(self, newSingleEta):
        return self.__mul__(newSingleEta)


class Beta(object):
    """docstring for Beta
        {0: [0, 1, {11: 2, 12: 3}], 11: []}
    """
    def __init__(self, betadict, latex=False):
        super(Beta, self).__init__()
        self.betadict = betadict
        self.values = betadict
        self.latex = latex
        self.string = self.__str__(latex=latex)

    def __str__(self, latex=False):
        values = self.values
        count_b = 0
        string_out = '0'
        for keyi, valuei in values.items():
            
            if isinstance(valuei, (int, float)):
                valuei = [valuei]
                return self.__str__(latex=latex)
            elif isinstance(valuei, list):
                seta = SingleEta(etalist=valuei, latex=latex)
                str0 = seta.__str__(latex=latex)
                n_i = seta.n
                if n_i == 0:
                    continue
                else:
                    str0 = seta.__str__(latex=latex)

                if keyi == 0:
                    string_beta = '{}'.format(str0)
                else:
                    if n_i == 1:
                        if latex:
                            string_beta = '{}\\beta_{{{}}}'.format(str0, keyi)
                        else:
                            string_beta = '{} * {}_{}'.format(str0, chr(946), keyi)
                    else:
                        if latex:
                            string_beta = '({})\\beta_{{{}}}'.format(str0, keyi)
                        else:
                            string_beta = '({}) * {}_{}'.format(str0, chr(946), keyi)

                if count_b == 0:
                    string_out = string_beta
                else:
                    string_out = '{} + {}'.format(string_out, string_beta)

            else:
                raise TypeError('Only support int or float or list')
            count_b = count_b + 1
        return string_out

    def __add__(self, newBeta):
        if isinstance(newBeta, (int, float)):
            beta_new = Beta(betadict={0: [newBeta]}, latex=self.latex)
            return self + beta_new
        elif isinstance(newBeta, SingleEta):
            beta_new = Beta(betadict={0: newBeta.values}, latex=self.latex)
            return self + beta_new
        elif isinstance(newBeta, (Beta, SingleBeta)):
            beta_new_val = newBeta.values
            result = copy.deepcopy(self.values)
            for keyi, vali in beta_new_val.items():
                if keyi in result:
                    val_self_list = result[keyi]
                    if len(val_self_list) > len(vali):
                        result_list_tmp = copy.deepcopy(val_self_list)
                        n_min = len(vali)
                    else:
                        result_list_tmp = copy.deepcopy(vali)
                        n_min = len(val_self_list)

                    for i in range(n_min):
                        if isinstance(val_self_list[i], (int, float)) and isinstance(vali[i], (int, float)):
                            result_list_tmp[i] = val_self_list[i] + vali[i]
                        elif isinstance(val_self_list[i], (int, float)) and isinstance(vali[i], dict):
                            dct_tmp = copy.deepcopy(vali[i])
                            dct_val_tmp = val_self_list[i]
                            if 0 in dct_tmp:
                                dct_tmp[0] = dct_tmp[0] + dct_val_tmp
                            else:
                                dct_tmp[0] = dct_val_tmp
                            result_list_tmp[i] = dct_tmp
                        elif isinstance(val_self_list[i], dict) and isinstance(vali[i], (int, float)):
                            dct_tmp = copy.deepcopy(val_self_list[i])
                            dct_val_tmp = vali[i]
                            if 0 in dct_tmp:
                                dct_tmp[0] = dct_tmp[0] + dct_val_tmp
                            else:
                                dct_tmp[0] = dct_val_tmp
                            result_list_tmp[i] = dct_tmp
                        elif isinstance(val_self_list[i], dict) and isinstance(vali[i], dict):
                            dct_tmp = copy.deepcopy(val_self_list[i])
                            for key_val, val_val in vali[i].items():
                                if key_val in dct_tmp:
                                    dct_tmp[key_val] = dct_tmp[key_val] + val_val
                                else:
                                    dct_tmp[key_val] = val_val
                            result_list_tmp[i] = dct_tmp
                        else:
                            raise TypeError('Should be dict or int or float.')
                    result[keyi] = result_list_tmp
                else:
                    result[keyi] = vali
            return Beta(betadict=result, latex=self.latex)

    def __radd__(self, newBeta):
        return self.__add__(newBeta)

    def __mul__(self, newBeta):
        val_self = self.values
        #print(val_self)
        if isinstance(newBeta, (int, float)):
            #print(val_self)
            result = {}
            for keyi, vali in val_self.items():
                result[keyi] = []
                for i, valii in enumerate(vali):
                    if isinstance(valii, (int, float)):
                        result[keyi].append(valii * newBeta) 
                    elif isinstance(valii, dict):
                        result[keyi].append({key_tmp: newBeta * val_tmp for key_tmp, val_tmp in valii.items()})
            #result = {keyi: [valii * newBeta for valii in vali] for keyi, vali in val_self.items()}
        elif isinstance(newBeta, SingleEta):
            #result = {0: newBeta.values}
            return self * Beta(betadict={0: newBeta.values}, latex=self.latex)
        elif isinstance(newBeta, (SingleBeta, Beta)):
            #print()
            val_new = newBeta.values
            result = {}
            for keyi_self, vali_self in val_self.items():
                for keyi_new, vali_new in val_new.items():
                    result_tmp = SingleEta(etalist=vali_self) * SingleEta(etalist=vali_new)
                    if (keyi_self == 0) and (keyi_new == 0):
                        result_ind = keyi_self
                        #result[keyi_self] = result_tmp.values
                    elif keyi_self == 0:
                        result_ind = keyi_new
                    elif keyi_new == 0:
                        result_ind = keyi_self
                    else:
                        continue

                    if result_ind in result:
                        result_tmp = result_tmp + SingleEta(etalist=result[result_ind])
                    result[result_ind] = result_tmp.values
        else:
            raise TypeError('* for SingleBeta should be int or float or SingleBeta or Beta.')
        return Beta(betadict=result, latex=self.latex)

    def __rmul__(self, newBeta):
        return self.__mul__(newBeta)


    @classmethod
    def from_SingleEta(cls, beta_ind):
        sbeta = SingleBeta(beta_ind=beta_ind)
        return Beta(betadict=sbeta.values)
        

class SingleBeta(object):
    """docstring for SingleBeta
    eta * delta_ij + (1 + 2*eta) * beta_ij
    """
    def __init__(self, beta_ind, latex=False):
        super(SingleBeta, self).__init__()
        self.beta_ind = beta_ind
        self.latex = latex
        ind_list = [int(i) for i in str(beta_ind)]
        if len(ind_list) < 2:
            raise ValueError('The length of beta_ind should equal to 2.')
        elif len(ind_list) > 2:
            logging.warning('The length of beta_ind larger than 2 and only the first two element are taken into consideration.')
            beta_ind = int(str(beta_ind)[0:2])
        else:
            pass
        if ind_list[0] == ind_list[1]:
            self.values = {0: [0, 1], beta_ind: [1, 2]}
        else:
            self.values = {0: [0], beta_ind: [1, 2]}
        self.string = self.__str__(latex=self.latex)

    def __str__(self, latex=False):
        values = self.values
        return Beta(betadict=self.values, latex=latex).string
    
    def __add__(self, newSigleBeta):
        beta_self = Beta(betadict=self.values, latex=self.latex)
        return beta_self + newSigleBeta

    def __radd__(self, newSigleBeta):
        return self.__add__(newSigleBeta)

    def __sub__(self, newSigleBeta):
        pass

    def __rsub__(self, newSigleBeta):
        return self.__sub__(newSigleBeta)

    def __mul__(self, newSigleBeta):
        val_self = self.values
        #print(val_self)
        if isinstance(newSigleBeta, (int, float)):
            result = {keyi: [valii * newSigleBeta for valii in vali] for keyi, vali in val_self.items()}
            #print(result)
        elif isinstance(newSigleBeta, SingleEta):
            result = {0: newSigleBeta.values}
        elif isinstance(newSigleBeta, (SingleBeta, Beta)):
            #print()
            val_new = newSigleBeta.values
            result = {}
            for keyi_self, vali_self in val_self.items():
                for keyi_new, vali_new in val_new.items():
                    result_tmp = SingleEta(etalist=vali_self) * SingleEta(etalist=vali_new)
                    if (keyi_self == 0) and (keyi_new == 0):
                        result_ind = keyi_self
                        #result[keyi_self] = result_tmp.values
                    elif keyi_self == 0:
                        result_ind = keyi_new
                    elif keyi_new == 0:
                        result_ind = keyi_self
                    else:
                        continue

                    if result_ind in result:
                        result_tmp = result_tmp + SingleEta(etalist=result[result_ind])
                    result[result_ind] = result_tmp.values
        else:
            raise TypeError('* for SingleBeta should be int or float or SingleBeta or Beta.')
        return Beta(betadict=result, latex=self.latex)

    def __rmul__(self, newSigleBeta):
        return self.__mul__(newSigleBeta)

    def __pow__(self, n):
        self_bk = copy.deepcopy(self)
        if n == 0:
            return 1.
        elif n < 0:
            raise ValueError('n should not less than 0 in pow of SingleBeta.')
        for i in range(n-1):
            self_bk = self_bk * self
        return self_bk

    '''        
    def __add__(self, newSigleBeta):
        result = copy.deepcopy(self).values #dict
        if isinstance(newSigleBeta, (int, float)):
            result[0][0] = result[0][0] + newSigleBeta
            return SingleBeta(beta_ind, latex)
        elif isinstance(newSigleBeta, SingleEta):
            result[0]
        elif isinstance(newSigleBeta, SingleBeta):
            pass
    '''

    @classmethod
    def from_SingleEta(cls, setaslist):
        pass



CrystalType = 'O'
Ord = 2
fe1 = FullExpan(CrystalType=CrystalType, Ord=Ord)
#fe1 = FullExpan.from_file(filename='C1-Ord5.yaml')
print(fe1)
fexp = fe1.fullexpression(tind=11, withPrecoef=False, Ord=7)


#for i in [0.1, 0.25, 0.38472, 0.19999999999999999]:
#    print('{} = {}'.format(i, frct.Fraction(i).limit_denominator()))


'''
CrystalType = 'C1'
Ord = 5
fe1 = FullExpan(CrystalType=CrystalType, Ord=Ord)
print('Full expansion for {} and {} :'.format(CrystalType, Ord))
print(fe1)
fe1.write2file()

beta_fe1 = fe1.to_beta()
print(beta_fe1)
print('')

print('Derivations:')
d_eta_ind = [11, 12]
for d_eta_indi in d_eta_ind:
    print('C{}'.format(d_eta_indi))
    dfe1 = fe1.derivation(d_eta_ind=d_eta_indi)
    print(dfe1)
    beta_dfe1 = dfe1.to_beta()
    print(beta_dfe1)
    dfe1.write2file(filename='{}-Ord{}-d{}.yaml'.format(CrystalType, Ord, d_eta_indi))
    print('')
'''


