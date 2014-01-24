import pickle
from glob import glob

A = glob('*4x3')
for filename in A:
    try:
        B = pickle.load(open(filename))
        d = {}
        d['path'] = B['path']
        d['electronic energy'] = B['electronic   energy']
        d['calculation script'] = B['calculation script']
        d['author'] = B['author']
        pickle.dump(d,open(filename,'w'))
    except KeyError:
        print filename+' is correct already!'
