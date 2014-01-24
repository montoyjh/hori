import pickle
from glob import glob
import os
newdir = '/nfs/slac/g/suncatfs/montoyjh/MAIN/PAR_NEB/adsorbates_beef/'
try:
    os.mkdir(newdir)
except:
    print newdir, ' failed'
startjobs = False
filelist = glob('*_Cu-fcc211')
for filename in filelist:
    data = pickle.load(open(filename))
    #print data['path'][10:]
    try:
        os.mkdir(newdir+filename.split('_')[0])
    except:
        print filename, ' failed'
    os.system('cp /nfs/slac/g/suncatfs/aap/niflheim'+data['path'][10:]+'/qn.traj ' \
             +newdir+filename.split('_')[0])
    os.system('cp /nfs/slac/g/suncatfs/montoyjh/SCRIPTS/beefRELAX.py '+\
              newdir+filename.split('_')[0])
    os.chdir(newdir+filename.split('_')[0])
    if startjobs:
        os.system('beefsub beefRELAX.py')
    os.chdir('/nfs/slac/g/suncatfs/montoyjh/usr/pylib/pylib2/hori/data/electronic-energies/')
    #print os.getcwd()


