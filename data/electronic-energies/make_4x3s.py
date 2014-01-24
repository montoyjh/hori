import pickle
from glob import glob
import sys

metal = sys.argv[1]
cu_pkl_list = glob('*_'+metal+'-fcc211')
ftlist = glob('*_'+metal+'-fcc211_4x3')
ftlistnew = []
for a in ftlist:
    ftlistnew.append(a[:-4])
convertlist=[]

for pkl in cu_pkl_list:
    if pkl in ftlistnew:
        print pkl+' already converted'
    else:
        convertlist.append(pkl)

e_slab_33 = pickle.load(open(metal+'-fcc211'))['electronic energy']
e_slab_43 = pickle.load(open(metal+'-fcc211_4x3'))['electronic energy']

for name in convertlist:
    f = open(name+'_4x3','w')
    e_ads_33 = pickle.load(open(name))['electronic energy']
    e_ads_43 = e_ads_33-e_slab_33+e_slab_43
    d = {'electronic energy': e_ads_43,
         'author' : 'JHM',
         'note' : 'calculated by taking 3x3 ee, subtracting 3x3 slab ee, \
         and adding 4x3 slab ee'}
    pickle.dump(d,f)
    f.close()

print convertlist
