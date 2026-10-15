#!python
#
from elastic3rd.symmetry.full_expan import Eta, MultiEta, SingleEta, SingleBeta, Beta, FullExpan

def test_Eta():
    eta_ind = 11
    pre_coef = 2
    exponent = 1
    eta1 = Eta(eta_ind=eta_ind, pre_coef=pre_coef, exponent=exponent)
    #test for __str__
    print('eta1: {}'.format(eta1))
    print(eta1.eta)
    # 2.0*η_11
    
    #test for from_dict
    etadict = {11: [1, 1]}
    eta2 = Eta.from_dict(etadict)
    print('eta2: {}'.format(eta2))

    #test for +
    eta3 = eta1 + eta2
    print('eta1 + eta2 : {}'.format(eta3))

    #test for -
    eta4 = eta1 - eta2
    print('eta1 - eta2 : {}'.format(eta4))

    #test for *
    eta5 = eta1 * eta2
    print('eta1 * eta2 : {}'.format(eta5))

    #TODO
    #eta6 = eta1 / eta2
    #print('eta1 / eta2 : {}'.format(eta6))
    
    eta7 = Eta.from_dict({12: [2, 2]})
    print('eta7 : {}'.format(eta7))

    #test for + for different eta_ind
    meta1 = eta1 + eta7
    print('eta1 + eta7 : {}'.format(meta1))
    print(type(meta1))

    #test for * for different eta_ind
    meta2 = eta1 * eta7
    print('eta1 * eta7 : {}'.format(meta2))

    eta8 = Eta(eta_ind=5, pre_coef=1, exponent=1)
    print('eta8 : {}'.format(eta8))

    eta9 = eta8.voigt()
    print('eta9 : {}'.format(eta9))

    eta10 = eta8.update(pre_coef=2)
    print('eta10: {}'.format(eta10))

    print('eta10 ({}) to dict is {}'.format(eta10, eta10.to_dict()))

    eta11 = eta10 + 0

    print(type(eta11))
    print(eta11.to_dict())
    print('eta10 + 1 : {}'.format(eta11))


def test_MultiEta():
    eta_inds = [[11, 12], [22]]
    pre_coefs = [2, 4]
    exponents = [[2, 1], [3]]
    meta1 = MultiEta(eta_inds=eta_inds, pre_coefs=pre_coefs, exponents=exponents)
    print('meta1 : {}'.format(meta1))

    eta1 = Eta(eta_ind=11, pre_coef=2, exponent=2)
    eta2 = Eta(eta_ind=21, pre_coef=1, exponent=1)
    eta3 = Eta(eta_ind=12, pre_coef=2., exponent=1)
    meta2 = MultiEta.from_etas(etas=[[eta1, eta3], eta2])
    print('meta2 : {}'.format(meta2))

    etadict = {'coef': [2, 2], 'eta': [{11: 2, 21: 2}, {21: 4}]}
    meta3 = MultiEta.from_dict(etadict)
    print('meta3 : {}'.format(meta3))

    meta4 = meta1 + meta2
    print('meta1 + meta2 : {}'.format(meta4))

    meta5 = meta2 * meta3
    print('meta2 * meta3 : {}'.format(meta5))

    print('meta2 * meta3 to_dict is {}'.format(meta5.to_dict()))

    print('meta1 * eta1 : {}'.format(meta1 * eta1))

    print('eta1 * meta1 : {}'.format(eta1 * meta1))

    print('eta1 + meta1 : {}'.format(eta1 + meta1))

    print('meta1 + eta1 : {}'.format(meta1 + eta1))

    meta6 = MultiEta.from_dict({'coef': [1, 0.5, 0.5], 'eta': [{11: 2}, {12: 2}, {21: 1, 12: 1}]})
    meta7 = MultiEta.from_dict({'coef': [1], 'eta': [{12: 1, 21: 1}]})
    print('meta6 : {}'.format(meta6))
    print('meta7 : {}'.format(meta7))
    print('meta6 + meta7 : {}'.format(meta6 + meta7))

def test_SingleEta():
    etalist1 = [0, 1, {0: 2, 11: 2, 12: 1}, 1, {0: 1, 23: 1, 33: 2}]
    seta1 = SingleEta(etalist=etalist1)
    print('seta1: {}'.format(seta1))

    etalist2 = [1, 2, {0: 3, 11: 2, 13: 2}, 1, {0: 1, 22: 1, 33: 2}, {0: 0, 23: 1}]
    seta2 = SingleEta(etalist=etalist2)
    print('seta2: {}'.format(seta2))

    seta3 = seta1 + seta2
    print('seta1 + seta2 : {}'.format(seta3))
    print('seta1 + 2 : {}'.format(seta1 + 2))
    print('2 + seta1 : {}'.format(2 + seta1))

    print('seta1: {}'.format(seta1))
    print('seta1 * 2 : {}'.format(seta1 * 2))

    seta4 = SingleEta(etalist=[1, 1, 3])
    print('seta4: {}'.format(seta4))

    seta5 = SingleEta(etalist=[1, 1, {0: 3, 11: 2, 12: 3}])
    print('seta5: {}'.format(seta5))

    print('seta4 * seta5: {}'.format(seta4 * seta5))

    seta6 = SingleEta(etalist=[1, {0: 2, 12: 3}, 3])
    print('seta6: {}'.format(seta6))
    print('seta5: {}'.format(seta5))
    print('seta5 * seta6 : {}'.format(seta5 * seta6))

def test_SingleBeta():
    sbeta1 = SingleBeta(beta_ind=11)
    print('sbeta1: {}'.format(sbeta1))
    sbeta2 = SingleBeta(beta_ind=12)
    print('sbeta2: {}'.format(sbeta2))
    sbeta3 = SingleBeta(beta_ind=22)
    print('sbeta3: {}'.format(sbeta3))

    sbeta4 = sbeta1 + sbeta2
    print('sbeta1 + sbeta2 : {}'.format(sbeta4))

    beta1 = Beta(betadict={0: [1, 2, 3], 11: [0, 1, 2], 13: [1, 2]})
    print('beta1 : {}'.format(beta1))
    beta2 = sbeta1 + beta1
    print('sbeta1 + beta1 : {}'.format(beta2))
    print('1 + sbeta1 : {}'.format(1 + sbeta1))
    print('sbeta1 + 1 : {}'.format(sbeta1 + 1))

    beta3 = sbeta1 * sbeta2
    print('sbeta1 * sbeta2 : {}'.format(beta3))
    print('sbeta1 * 3 : {}'.format(sbeta1 * 3))
    print('3 * sbeta1 : {}'.format(3 * sbeta1))

    beta4 = sbeta1 ** 3
    print('sbeta1 ** 3 : {}'.format(beta4))
    print('sbeta2 ** 2 : {}'.format(sbeta2 ** 2))

def test_Beta():
    betadict = {0: [0, 1], 11: [1, 2, 3], 22: [0, 1, 1]}
    beta1 = Beta(betadict=betadict)
    print('beta1 : {}'.format(beta1))

    betadict2 = {0: [0], 12: [1, 2]}
    beta2 = Beta(betadict=betadict2)
    print('beta2 : {}'.format(beta2))

    betadict3 = {0: [0, 1], 11: [1, 2], 12: [1, 3]}
    beta3 = Beta(betadict=betadict3)
    print('beta3 : {}'.format(beta3))

    beta4 = beta1 + beta3
    print('beta1 + beta3 : {}'.format(beta4))
    print('beta1 + 1 : {}'.format(beta1 + 1))
    print('1 + beta1 : {}'.format(1 + beta1))

    beta5 = beta1 * beta3
    print('beta1 * beta3 : {}'.format(beta5))

    #beta5 = beta1 * beta3
    #print('beta1 * beta3 : {}'.format(beta5))
    
    betadict6 = {0: [0, 0, 1], 11: [1, 2, 3], 22: [0, 1, 1]}
    beta6 = Beta(betadict=betadict6)
    print('beta6 : {}'.format(beta6))

def test_FullExpan():
    CrystalType = 'C1'
    Ord=2
    fe1 = FullExpan(CrystalType=CrystalType, Ord=Ord)
    print('fe1 : {}'.format(fe1))
    for keyi, vali in fe1.values.items():
        print('{} : {}'.format(keyi, vali))
    print(fe1.values)

test_Beta()

#test_Eta()
#print('')
#test_MultiEta()
#

'''
sbeta = SingleBeta(beta_ind=11, latex=False)
print(sbeta)
sbeta2 = SingleBeta(beta_ind=22, latex=False)
print(sbeta2)

print('\n')
betadict = {0: [0, 1], 11: [1, 2, 3], 22: [0, 1, 1]}
beta1 = Beta(betadict=betadict, latex=latex)
print(beta1)



print(seta2)
beta4 = beta1 + seta2
print(beta4)

print('')
print(seta1)
seta4 = SingleEta(etalist=[1, 2], latex=False)
print(seta4)
seta3 = seta1 * 1
print(seta3)

print('')
print(sbeta)
print(sbeta2)
sbeta3 = sbeta * sbeta2
print(sbeta3)

print(beta1)
beta5 = sbeta * beta1
print(beta5)

#sbeta4 = sbeta * 2
#print(sbeta4)

print(sbeta)
beta6 = sbeta ** 2
print(beta6)
beta7 = sbeta ** 4
print(beta7)

'''

#'''
#full_expan()
'''
fe1 = FullExpan(Ord=2)
result = fe1.full_expan()
for keyi in result:
    metai = result[keyi]
    print('C{} : {}'.format(keyi, result[keyi]))
    dmetai = metai.derivation(d_eta_ind=11)
    print(dmetai)
    dmetai12 = metai.derivation(d_eta_ind=12)
    print(dmetai12)

print(fe1.to_dict())
fe1.write2file()
'''

'''
with open('C1-Ord2.yaml', 'r') as fid:
    a = yaml.load(fid)
print(a)

fe2 = FullExpan.from_dict(value_dict=a)
dfe2 = fe2.derivation(d_eta_ind=11)

fe3 = FullExpan.from_file(filename='C1-Ord2.yaml')
print(fe3)
'''

#print(result[111])

#'''

'''
meta6 = MultiEta.from_dict({'coef': [1, 0.5, 0.5], 'eta': [{11: 2}, {12: 2}, {21: 1, 12: 1}]})
meta7 = MultiEta.from_dict({'coef': [0.5], 'eta': [{12: 1, 21: 1}]})
print('meta6 : {}'.format(meta6))
print('meta7 : {}'.format(meta7))
print('meta6 + meta7 : {}'.format(meta6 + meta7))
'''

'''
eta1 = Eta.from_dict({11: [2, 2]})
print(eta1)

print(type(eta1))

beta1 = eta1.to_beta()
print(beta1)

metadict = {'coef': [1, 2, 3], 'eta': [{11: 2, 22: 1, 12: 1}, {11: 1, 13: 1}, {22: 1, 11: 2}]}
meta = MultiEta.from_dict(etasdict=metadict)
print(meta)

beta = meta.to_beta()
print(beta)

fe1 = FullExpan(CrystalType='C1', Ord=2)
print('')
print(fe1)
beta = fe1.to_beta()
print(beta)

print('')
dfe1 = fe1.derivation()
print(dfe1)
dfe1_beta = dfe1.to_beta()
print(dfe1_beta)
'''

'''
print(dfe1_beta.values)

print(dfe1.to_dict())
print(fe1.to_dict())
'''
