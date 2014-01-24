from ase.all import *
dirs_pref = '/nfs/slac/g/suncatfs/montoyjh/MAIN/PAR_NEB/'
rxns = {'COCO':'COCO_OCCO_TS',
        'CHOCO':'CHOCO_TS', 
        'CHOCHO':'CHOCHO_OCHCHO_TS', 
        'COCH2O':'COCH2O_OCCH2O_TS', 
        'CH2OCHO':'CH2OCHO_OCH2CHO_TS', 
        'CH2OCH2O':'CH2OCH2O_OCH2CH2O_TS'}
finals = {'COCO':'OCCO',
        'CHOCO':'OCCHO', 
        'CHOCHO':'OCHCHO', 
        'COCH2O':'OCCH2O', 
        'CH2OCHO':'OCH2CHO', 
        'CH2OCH2O':'OCH2CH2O'}

metals = ['Au', 'Ag', 'Pt', 'Pd','Zn', 'Ir']

import os, pickle

for rxn in rxns.keys():
    for metal in metals:
        # Extract barrier from highest fbl bond length
        inter_es = []
        for i in range(1,8):
            if rxn =='CH2OCHO' and metal=='Ag' and (i==1 or i==6 or i==7):
                continue
            if rxn =='CH2OCH2O' and metal=='Ag' and (i==6 or i==7):
                continue
            if rxn=='CH2OCH2O' and metal=='Pt' and (i==1 or i==2 or i==3 or i==4):
                continue
            if rxn=='CH2OCH2O' and metal=='Pd' and (i==2 or i==3 or i==4):
                continue
            if rxn=='COCH2O' and metal=='Ir' and (i>=4):
                continue
            inter_es.append(float(open(dirs_pref+rxn+'/metals/fbl/'+metal+'/inter_'+str(i)+'/out.energy').read()))
        barrier = max(inter_es)
        # Extract final state energy
        f_e =\
        float(open(dirs_pref+rxn+'/metals/runs/'+metal+'/final/out.energy').read())
        d = {'electronic energy': barrier,
             'author': 'JHM',
             'note': 'calculated using fixed bond length constrants near the adsorbate transition state geometry calculated using NEB on Cu-211',
             'path':dirs_pref+rxn+'/metals/fbl/'+metal+'/inter_'+str(i)+'/'}
        f = {'electronic energy': f_e,
             'author': 'JHM',
             'path': dirs_pref+rxn+'/metals/runs/'+metal+'/final/'}
        file_barrier = open(rxns[rxn]+'_'+metal+'-fcc211_4x3','w')
        file_final = open(finals[rxn]+'_'+metal+'-fcc211_4x3','w')
        pickle.dump(d,file_barrier)
        file_barrier.close()
        pickle.dump(f,file_final)
        file_final.close()
