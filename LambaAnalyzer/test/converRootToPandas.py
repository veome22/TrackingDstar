from root_pandas import read_root
import sys

cols = ['isSharedHit','trackPt','trackEta','trackPhi']
for i in xrange(16*16):
    cols.append('pixel_%i' % i)
    
df = read_root(sys.argv[1], columns=cols)
print df
df.to_hdf("pixelTrain.h5",key='df',mode='w')
print "done"
