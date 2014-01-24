import pickle
metal = 'Pd'
ads_dict = {'CO':'CO',
            'CHO':'CHO',
            'CH2O':'CH2O',
            'OCH3':'OCH3',
            'OH':'OH',
            'CO-CO':'COCO_OCCO_TS',
            'CHO-CO':'CHOCO_TS',
            'CH2O-CO':'COCH2O_OCCH2O_TS',
            'CHO-CHO':'CHOCHO_OCHCHO_TS',
            'CH2O-CHO':'CH2OCHO_OCH2CHO_TS',
            'CH2O-CH2O':'CH2OCH2O_OCH2CH2O_TS',
            'OCCO':'OCCO',
            'OCCHO':'OCCHO',
            'OCHCHO':'OCHCHO',
            'OCH2CHO':'OCH2CHO',
            'OCH2CH2O':'OCH2CH2O',
            'COOH':'COOH',
            #'OCH2CH2OH':'OCH2CH2OH',
            'C':'C',
            'OCHO':'OCHO',
            'COH':'COH',
            'CHOH':'CHOH',
            'CH':'CH',
            'CH2':'CH2',
            'CH3':'CH3',
            'CH2OH':'CH2OH'
           }

#References

CH4 = pickle.load(open('CH4'))['electronic energy']
H2O = pickle.load(open('H2O'))['electronic energy']
H2 = pickle.load(open('H2'))['electronic energy']
slab = pickle.load(open(metal+'-fcc211_4x3'))['electronic energy']
refs = {'CO': CH4+H2O-3.*H2,
            'CHO': CH4+H2O - 5./2. * H2,
            'CH2O':CH4+H2O - 2.*H2,
            'OCH3':CH4+H2O - 3./2.*H2,
        'OH':H2O - 1./2.*H2,
        'C':CH4 - 2.*H2,
        'OCHO':CH4+2*H2O - 7./2.*H2,
        'CHOH':CH4 + H2O - 2 * H2,
        'CH': CH4 - 3./2. * H2,
        'CH2':CH4 - H2,
        'CH3':CH4 - 1./2. * H2,
        'CH2OH':CH4 + H2O - 3./2. * H2
        }
refs['COH'] = refs['CHO']
refs['CO-CO'] = 2.*(CH4+H2O-3.*H2)
refs['CHO-CO'] = refs['CHO']+refs['CO']
refs['CH2O-CO'] = refs['CO']+refs['CH2O']
refs['CHO-CHO'] = refs['CHO']*2.
refs['CH2O-CHO'] = refs['CHO']+refs['CH2O']
refs['CH2O-CH2O'] = refs['CH2O']*2
refs['OCCO'] = refs['CO-CO']
refs['OCCHO'] = refs['CHO-CO']
refs['OCHCHO'] = refs['CHO-CHO']
refs['OCH2CHO'] = refs['CH2O-CHO']
refs['OCH2CH2O'] = refs['CH2O-CH2O']
refs['COOH'] = refs['CO']+refs['OH']
refs['OCH2CH2OH'] = refs['OCH2CH2O'] + 1./2.*H2

f = open('ads_jhm.csv','w')
for ads in ads_dict.keys():
    f.write(metal+',')
    f.write(ads+'*,')
    ads_e = pickle.load(open(ads_dict[ads]+'_'+metal+'-fcc211_4x3'))['electronic energy'] -\
                                                              refs[ads] - slab
    f.write(str(ads_e)+',')
    f.write('fcc,')
    f.write('211,')
    f.write('Peterson/Montoya (Unpublished 2012)\n')
f.close()
