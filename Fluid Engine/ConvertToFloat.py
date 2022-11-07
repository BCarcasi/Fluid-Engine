# magick convert *.png simsph.gif
import os
import sys
import numpy as np
import glob


outpath = 'processed1_test'
if (os.path.exists(outpath) == False):
    os.mkdir(outpath)

path = 'sph_sim_output'
files = sorted(glob.glob(path + "/*"))

for i,f in enumerate(files):
    print (f)
    df = np.loadtxt(f)    

    floats = df.astype(float)
    np.save(outpath + '/outputFloat-%02d.npy' % i, floats)

    print("running" + str(i))
