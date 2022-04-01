#!/usr/bin/env python3

import numpy as np
import time as tm

import pymirrors as mrr

from os.path import exists
from os import mkdir
import matplotlib.pyplot as pl
pl.ioff()

#
nthread = 8
#
#
#
#
idealc = True
cdust = False
plotsec = False
printout = False
savedata = True
overwrite = False
#
odirseed = 'results_article/'
#
#
#
telescope = 'gtc'
#
# Checked:
methods = ['equal']
#
# Orientation:
alpha = 0.
x_alpha = 0.
y_alpha = 0.



#
# Number of rays to consider in the radial and proxi for the azimuth directions:
nums=[4096]
#
# Wavelengths:
lambs = np.linspace(2100.,24000.,219+1)
#
# Timing control parameters:
deltat=1
cleandust=1.
tstep=1.
tlong=1
mltch=1
#

################################################################################

if (not exists(odirseed)):
  mkdir(odirseed)

for method in methods:

  for lamb in lambs:
  
    for num in nums:
    
      outdir = '%s/tel%s_%s_alpha%.2f_%08i' % (odirseed \
          , telescope.upper(), method, alpha, num)
      if (not exists(outdir)):
        mkdir(outdir)
      
      oname = '%s/w%.2f_Dt%.2f_St%.2f_%08i_%010.2f_%i' \
            % (outdir, lamb, deltat, tstep, tlong, cleandust, mltch, )
      
      if ( (not(exists(oname))) \
          or (overwrite==True) ):
      #
        t0 = tm.time()
        teles = mrr.init_telescope(telescope, num)
        secondary = mrr.secondary(teles)
        beam = mrr.init_beam(alpha, x_alpha, y_alpha,degrees=True)
        
        t1 = tm.time()
        segments = mrr.primary(teles, method, tstep, tlong, cleandust \
            , deltat=deltat, ideal=idealc, multiplechange=mltch)
        
        t2 = tm.time()
        materials = mrr.dirt_timing(lamb, segments.time_map)
        
        t3 = tm.time()
        mrr.get_geometry(teles, beam, secondary, segments, osecondary=plotsec \
            , nthreads=nthread)
      
        t4 = tm.time()
        avg_mat,mat,npts_seg = mrr.get_mueller_time(segments, materials \
            , cdust=cdust, nthreads=nthread)
        
        t5 = tm.time()
        
        print('   ***   %.3f Sec. ===   ' % (t5-t0,))
        #
        if (savedata==True):
          np.savez(oname, mat=mat, avg_mat=avg_mat, npts_seg=npts_seg)

        if (printout==True):
          print(avg_mat)

      else:
        print('exists!')
      
