#!/usr/bin/env python3

import numpy as np
from pdb import set_trace as stop
import time as tm

import apy_utils as uts
import mirror_lib_v03 as mrr

#TO BE MOVED WHEN SATISFACTORY RESULTS ARE OBTAINED
from os.path import exists
from os import mkdir
import matplotlib.pyplot as pl
pl.ion()

imfact = 1./2. * 2.

pl.rcParams.update({'font.size':15*imfact})
pl.rcParams.update({'xtick.labelsize':14})  
pl.rcParams.update({'ytick.labelsize':14})  
pl.rcParams.update({'text.usetex':True}) 
pl.rcParams['xtick.minor.visible'] = True 
pl.rcParams['ytick.minor.visible'] = True


#---------------------------------------------------------------------
#
def auto_label_subplots(iax, label,px=None,py=None):
  #
  if (px==None):
    px=0.9
  if (py==None):
    py=0.9
#  xpts = iax.get_position().get_points()[:,0]
#  ypts = iax.get_position().get_points()[:,1]
  #
  iax.text(0.95, 0.95, label, horizontalalignment='center'\
     ,verticalalignment='center', transform=iax.transAxes\
     , bbox={'facecolor': 'white', 'alpha': 0.8, 'pad': 0.2\
     , 'edgecolor':'k', 'boxstyle':'round'})
  #
  return
#
#---------------------------------------------------------------------

quadrature = 'polar'
#quadrature = 'cartesian'


#
telescopes = ['gtc', 'eelt']#[::-1]
#telescopes = ['sgtc', 'gtc']#[::-1]
#
# Checked:
method = 'equal'
#
print('Telescope chosen!')
alpha = 0.
x_alpha = 0.
y_alpha = 0.
print('Orientation chosen!')
#
#period = 365 #798.#365.24 * 2.
deltat= 1.#365./36.
cleandust = 1.      # DUST CLEANING TIME FREQUENCY [in days?]
tstep = 1. # SIMULATION TIME STEP   $NUMBER OF TIMESTEPS FOR EACH MIRROR
tlong = 1#365*2+6#798. # Simulated days

print('Time parameters chosen!')
lamb = 5.e3    #ANGSTROM

print('Wavelength chosen!')
num = 4096
num = 25

print('Number of rays per segment chosen!')

################################################################################

idealc=True#False#True
cdust=False#True
plotsec=False#True#False
printout=False#True
savedata=False#False
overwrite=False#True

################################################################################

opath = 'figures_202110_v03'
oname = 'figure_2.pdf'

option=2

xlabels = [r'x [m]', r'x [m]']
ylabels = [r'y [m]', r'y [m]']
pane_labs = [r'a)', r'b)', r'c)']

if (option==1):

  pl.close(2)
  fig, ax = pl.subplots(num=2,nrows=1,ncols=4,figsize=(12,4))
  
   
  for itt, telescope in enumerate(telescopes):
   
    print('Initialize telescope:')
    t0 = tm.time()
    teles = mrr.init_telescope(telescope, num)
    print('Initialize secondary:')
    secondary = mrr.secondary(teles)
    print('Initialize beam:')
    beam = mrr.init_beam(alpha, x_alpha, y_alpha,degrees=True)
    
    segments = mrr.primary(teles, method, tstep, tlong, cleandust \
        , deltat=deltat, ideal=idealc)
    
    materials = mrr.dirt_timing(lamb, segments.time_map)
    
    #mrr.get_geometry(teles, beam, secondary, segments, osecondary=False
    if (itt<1):
      mrr.plot_primary(teles, beam, secondary, segments \
          , pax=ax[itt*len(telescopes)+0], prays=False, fcol='tab10', layout=quadrature)
    else:
      mrr.plot_primary(teles, beam, secondary, segments \
          , pax=ax[itt*len(telescopes)+0], prays=False, layout=quadrature)
    
    ax[itt*len(telescopes)+0].set_aspect('equal')
  
    mrr.plot_primary(teles, beam, secondary, segments \
        , pax=ax[itt*len(telescopes)+1], prays=True, layout=quadrature)
  
    ax[itt*len(telescopes)+1].set_aspect('equal')
  
    mrr.get_geometry(teles, beam, secondary, segments, osecondary=plotsec)

elif (option==2):

  pl.close(2)
  fig, ax = pl.subplots(num=2,nrows=3,ncols=1,figsize=(6*imfact,17*imfact))
  pl.close(3)
  fig2, ax2 = pl.subplots(num=3,nrows=1,ncols=1,figsize=(9,9))
  
   
  for itt, telescope in enumerate(telescopes):
   
    print('Initialize telescope:')
    t0 = tm.time()
    teles = mrr.init_telescope(telescope, num)
    print('Initialize secondary:')
    secondary = mrr.secondary(teles)
    print('Initialize beam:')
    beam = mrr.init_beam(alpha, x_alpha, y_alpha,degrees=True)
    
    segments = mrr.primary(teles, method, tstep, tlong, cleandust \
        , deltat=deltat, ideal=idealc)
    
    materials = mrr.dirt_timing(lamb, segments.time_map)
    
    #mrr.get_geometry(teles, beam, secondary, segments, osecondary=False
    if (itt<1):
      mrr.plot_primary(teles, beam, secondary, segments \
          , pax=ax[itt*len(telescopes)+0], prays=False, fcol='tab10', layout=quadrature)
      ax[itt*len(telescopes)+0].xaxis.set_ticklabels([])
      ax[itt*len(telescopes)+0].set_ylabel(ylabels[itt])
      auto_label_subplots(ax[itt*len(telescopes)+0], pane_labs[itt*len(telescopes)+0])
    else:
      mrr.plot_primary(teles, beam, secondary, segments \
          , pax=ax[itt*len(telescopes)+0], prays=False, layout=quadrature, highlight=np.pi/2.)
      ax[itt*len(telescopes)+0].set_xlabel(xlabels[itt])
      ax[itt*len(telescopes)+0].set_ylabel(ylabels[itt])
      auto_label_subplots(ax[itt*len(telescopes)+0], pane_labs[itt*len(telescopes)+0])

    ax[itt*len(telescopes)+0].set_aspect('equal')
  
    if (itt<1):
      #mrr.plot_primary(teles, beam, secondary, segments \
      mrr.new_plot_primary(teles, beam, secondary, segments \
          , pax=ax[itt+1], prays=True)
  
      ax[itt+1].set_aspect('equal')
      ax[itt+1].set_ylabel(ylabels[itt])
      ax[itt+1].set_xlabel(xlabels[itt])
      auto_label_subplots(ax[itt+1], pane_labs[itt+1])
#      mrr.new_plot_primary(teles, beam, secondary, segments \
#          , pax=ax2, prays=True, layout=quadrature)
#  
#      ax2.set_aspect('equal')
#      xcirc = np.linspace(0., 2.*np.pi, 1001)
#      ycirc = np.cos(xcirc)
#      xcirc = np.sin(xcirc)
#      ax2.plot(xcirc * teles.rad_max, ycirc*teles.rad_max)
#      ax2.set_ylabel(ylabels[itt])
#      ax2.set_xlabel(xlabels[itt])
#      auto_label_subplots(ax2, pane_labs[itt+1])
  
    mrr.get_geometry(teles, beam, secondary, segments, osecondary=plotsec)

  

fig.tight_layout() 

# Finally, set panel b) closer to a):

pos1 = ax[0].get_position()
pos2 = ax[1].get_position()

deltay = (pos1.y0 - pos2.y1)/2.

pos2.y0=pos2.y0+deltay         
pos2.y1=pos2.y1+deltay         
ax[1].set_position([pos2.x0 ,pos2.y0, pos2.width, pos2.height])

#


if (not exists(opath)):
  mkdir(opath)

fig.savefig("%s/%s" % (opath, oname, ), dpi=600, format='pdf')



