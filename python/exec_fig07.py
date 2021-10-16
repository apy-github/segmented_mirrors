#!/usr/bin/env python3.8

import numpy as np
from pdb import set_trace as stop
import time as tm

import apy_utils as uts
import mirror_lib_v03 as mrr

#TO BE MOVED WHEN SATISFACTORY RESULTS ARE OBTAINED
from os.path import exists
from os import mkdir
import matplotlib.pyplot as pl
pl.ioff()

#
#
telescope = 'gtc'
#telescope = 'anular'
#telescope = 'sgtc'
telescope = 'eelt'
#telescope = 'seelt'
#telescope = 'aeelt'
#
# Checked:
#method = 'azimuthal'
#method = 'lineal'
#method = 'random'
#method = 'symmetric'
#methods = ['equal']
methods = []
#methods.append('azimuthal')
#methods.append('lineal')
#methods.append('random')
methods.append('symmetric')
#
print('Telescope chosen!')
alpha = 0.
x_alpha = 0.
y_alpha = 0.
print('Orientation chosen!')
print('Number of rays per segment chosen!')


deltat=1
cleandust=1.
tstep=1.#4.4
tlong=798//2+1
lambs = [5000.]
nums=[720]

mltch=2
nums=[4096]



idealc=False#True
cdust=False#True
plotsec=False#True#False
printout=True#True
savedata=True#False#False
overwrite=False#False

odirseed = 'results_202110_v03/'
if (not exists(odirseed)):
  mkdir(odirseed)

#True

################################################################################

for method in methods:

  for lamb in lambs:
  
    for num in nums:
    
      outdir = '%s/tel%s_%s_alpha%.2f_%08i' % (odirseed \
          , telescope.upper(), method, alpha, num)
      if (not exists(outdir)):
        mkdir(outdir)
      
      # Posibles necesidades: tstep, tlong, cleansdust, period, tlong
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
      #, period=period)
        #><import matplotlib.pyplot as pl
        #><pl.ion()
        #><pl.figure(1)
        #><pl.clf() ; pl.imshow(segments.time_map[:,:,0].T) ; pl.show()
        #><#stop()
        #><pl.figure(2)
        #><pl.clf() ; pl.imshow(segments.time_map[:,:,1].T) ; pl.show()
        #><stop()
        
        t2 = tm.time()
        materials = mrr.dirt_timing(lamb, segments.time_map)
        
        t3 = tm.time()
        #mrr.get_geometry(teles, beam, secondary, segments, osecondary=False)
        mrr.get_geometry(teles, beam, secondary, segments, osecondary=plotsec)
      
        t4 = tm.time()
        avg_mat,mat,npts_seg = mrr.get_mueller_time(segments, materials, cdust=cdust)
        to_save_mat = mat * 1.
        
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
        to_save_mat = mat * 1.
        if (printout==True):
          print(avg_mat)

      else:
        print('exists!')
      
