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

idealc = False
cdust = False
plotsec = False
printout = False
savedata = True
overwrite = False

odirseed = 'results_article/'

#
telescope = 'gtc'
#
methods = []
methods.append('azimuthal')
methods.append('lineal')
methods.append('random')
methods.append('symmetric')
#methods = ['lineal']
#
print('Telescope chosen!')
alpha = 0.
x_alpha = 0.
y_alpha = 0.
print('Orientation chosen!')
print('Number of rays per segment chosen!')

if (telescope=='eelt'):
  #
  #
  nums = [4096]
  #
  #
  lambs = [5000.]
  #
  #
  deltat = 1.
  tstep = 1. # SIMULATION TIME STEP   $NUMBER OF TIMESTEPS FOR EACH MIRROR
  tlong = 798//2+1
  mltch = 2

  if (cdust==True):
    cleandust = 30.      # DUST CLEANING TIME FREQUENCY [in days?]
  else:
    cleandust = 1.      # DUST CLEANING TIME FREQUENCY [in days?]

elif (telescope=='gtc'):
  #
  #
  nums = [4096]
  #
  #
  lambs = [5000.]
  #
  #
  deltat = 10. # Days between segment exchange
  tstep = 10. # SIMULATION TIME STEP   $NUMBER OF TIMESTEPS FOR EACH MIRROR
  tlong = 370 # Simulated days
  mltch = 1

  if (cdust==True):
    cleandust = 30.      # DUST CLEANING TIME FREQUENCY [in days?]
  else:
    cleandust = 1.      # DUST CLEANING TIME FREQUENCY [in days?]



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
      
      if ( (not exists(oname)) \
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
      
