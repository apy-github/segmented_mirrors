#!/usr/bin/env python3

import numpy as np
from pdb import set_trace as stop
import time as tm

import apy_utils as uts
import pymirrors as mrr

#TO BE MOVED WHEN SATISFACTORY RESULTS ARE OBTAINED
from os.path import exists
from os import mkdir
import matplotlib.pyplot as pl
pl.ioff()

nthread = 8

idealc = False
cdust = True
plotsec = False
printout = True
savedata = True
overwrite = True

odirseed = 'results_2022_v05/'

#
telescope = 'gtc'
#
methods = []
methods.append('azimuthal')
methods.append('lineal')
methods.append('random')
methods.append('symmetric')
#
print('Telescope chosen!')
alpha = 0.
x_alpha = 0.
y_alpha = 0.
print('Orientation chosen!')
print('Number of rays per segment chosen!')

if (telescope=='eelt'):
  deltat = 1.
  lambs = [5000.]
  nums = [720]
  tstep = 1. # SIMULATION TIME STEP   $NUMBER OF TIMESTEPS FOR EACH MIRROR
  tlong = 798//2+1
  mltch = 2

  if (cdust==True):
    cleandust = 30.      # DUST CLEANING TIME FREQUENCY [in days?]
  else:
    cleandust = 1.      # DUST CLEANING TIME FREQUENCY [in days?]

elif (telescope=='gtc'):
  deltat = 10. # Days between segment exchange
  tstep = 10. # SIMULATION TIME STEP   $NUMBER OF TIMESTEPS FOR EACH MIRROR
  lambs = [5000.]
  nums = [720]
  tlong = 370 # Simulated days
  mltch = 1

  if (cdust==True):
    cleandust = 30.      # DUST CLEANING TIME FREQUENCY [in days?]
  else:
    cleandust = 1.      # DUST CLEANING TIME FREQUENCY [in days?]

nums = [4096]

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
      
      oname1 = '%s/full_mueller_matrix_w%.2f_Dt%.2f_St%.2f_%08i_%010.2f_%i.fits' \
            % (outdir, lamb, deltat, tstep, tlong, cleandust, mltch, )
      oname2 = '%s/avg_mueller_matrix_w%.2f_Dt%.2f_St%.2f_%08i_%010.2f_%i.fits' \
            % (outdir, lamb, deltat, tstep, tlong, cleandust, mltch, )
      oname3 = '%s/npts_mueller_matrix_w%.2f_Dt%.2f_St%.2f_%08i_%010.2f_%i.fits' \
            % (outdir, lamb, deltat, tstep, tlong, cleandust, mltch, )
      
      if ( (not ( exists(oname1) & exists(oname2) & exists(oname3) ) ) \
          or (overwrite==True) ):
      #
        print('Initialize telescope:')
        t0 = tm.time()
        teles = mrr.init_telescope(telescope, num)
        print('Initialize secondary:')
        secondary = mrr.secondary(teles)
        print('Initialize beam:')
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
        
        print('   ***   ')
        print((t1-t0))
        print((t2-t1))
        print((t3-t2))
        print((t4-t3))
        print((t5-t4))
        print('   ---   ')
        #
        if (savedata==True):
          uts.writefits_v3(mat, oname1, overwrite=1)
          uts.writefits_v3(avg_mat, oname2, overwrite=1)
          uts.writefits_v3(npts_seg, oname3, overwrite=1)

        if (printout==True):
          print(avg_mat)

      else:
        print('exists!')
      
