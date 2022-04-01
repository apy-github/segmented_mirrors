#!/usr/bin/env python2.7

import numpy as np
from pdb import set_trace as stop
import time as tm

import apy_utils as uts
import pymirrors as mrr

#TO BE MOVED WHEN SATISFACTORY RESULTS ARE OBTAINED
from os.path import exists
from os import mkdir
import matplotlib.pyplot as pl
from matplotlib import ticker
pl.ioff()
pl.rcParams.update({'font.size':15})
pl.rcParams.update({'xtick.labelsize':14})
pl.rcParams.update({'ytick.labelsize':14})
pl.rcParams.update({'text.usetex':True})
pl.rcParams['xtick.minor.visible'] = True
pl.rcParams['ytick.minor.visible'] = True

pl.rcParams['text.latex.preamble'] = [r"\usepackage{amsfonts}",]

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
  iax.text(px, py, label, horizontalalignment='center'\
     ,verticalalignment='center', transform=iax.transAxes\
     , bbox={'facecolor': 'white', 'alpha': 0.8, 'pad': 0.2\
     , 'edgecolor':'k', 'boxstyle':'round'})
  #
  return
#
#---------------------------------------------------------------------

#
#
# Remove exponential axes: start:
# Taken from: https://stackoverflow.com/questions/33416541/remove-axis-scale#33416790
# Note: some changes have been done.
# Note 2: additional line to each axes:
#

def scientificNotation(value):
    if value == 0:
        return '0'
    else:
        #e = np.floor(np.log10(np.abs(value)))
        e = np.log10(np.abs(value))
        if (np.abs(e)<2):
          return r'${:.2f}$'.format(value)
        ep = np.sign(e)*np.ceil(np.abs(e))
        #m = np.sign(value) * 10 ** (e - int(e))
        m = value * 10 ** (-ep)
        #return r'${:.1f} \cdot 10^{{{:d}}}$'.format(m, int(ep))
        #return r'${:.1f} _{O}^{{{:d}}}$'.format(m, int(ep))
        return r'$%.1f^{%d}$' % (m, int(ep))

formatter = ticker.FuncFormatter(lambda x, p: scientificNotation(x))

#
# Remove exponential axes: end.
#


################################################################################

# Main

opath = 'figures_2022_v05'
odirseed='results_2022_v05/'
oname = 'figure_4.pdf'

#
#
telescope = 'gtc'
#telescope = 'anular'
#telescope = 'gtc'
#telescope = 'eelt'
#
# Checked:
#method = 'azimuthal'
#method = 'lineal'
#method = 'random'
#method = 'symmetric'
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
#2. * 365.24 #TIME NEEDED FOR CHANGING ALL THE MIRRORS IN DAYS
cleandust = 1.      # DUST CLEANING TIME FREQUENCY [in days?]
tstep = 1. # SIMULATION TIME STEP   $NUMBER OF TIMESTEPS FOR EACH MIRROR
#
tlong = 1#365*2+6#798. # Simulated days
print('Time parameters chosen!')
#lamb = 3.e3    #ANGSTROM
lambs = np.linspace(2100.,24000.,219+1)[:]#(np.arange(7)+1.)*1000.+3000.
print('Wavelength chosen!')
nums = [720]#(np.arange(100,dtype=np.int32)+1)*10
nums = [4096]#(np.arange(100,dtype=np.int32)+1)*10
print('Number of rays per segment chosen!')

mltch=1

#odirseed='~/scratch/apy/segmented_mirrors/results/'

################################################################################


pl.close(2)
fig, ax = pl.subplots(num=2,nrows=4,ncols=4,figsize=(12,8),sharex=True)

res = np.zeros((len(nums), len(lambs), 4, 4))

for itn, num in enumerate(nums):
  
  for itl, lamb in enumerate(lambs):

    outdir = '%s/tel%s_%s_alpha%.2f_%08i' % (odirseed\
        , telescope.upper(), method, alpha, num)
    
    # Posibles necesidades: tstep, tlong, cleansdust, period, tlong
    oname2 = '%s/avg_mueller_matrix_w%.2f_Dt%.2f_St%.2f_%08i_%010.2f_%i.fits' \
        % (outdir, lamb, deltat, tstep, tlong, cleandust, mltch, )

    if (not exists(oname2)):
      continue
    
    int_mat = uts.readfits_v3(oname2)[1]
    res[itn,itl,:,:] = int_mat[0,:,:] * 1.
  
  nt, nx, ny = int_mat.shape

  for it_xxx in range(nx):
    for it_yyy in range(ny):
  
      if ((it_xxx==0)&(it_yyy==0)):
        ax[it_xxx, it_yyy].plot(lambs*1.e-4 \
            , res[itn,:, it_xxx, it_yyy], linewidth=1., color='k')
      else:
        ax[it_xxx, it_yyy].plot(lambs*1.e-4 \
            , res[itn,:, it_xxx, it_yyy] / res[itn,:, 0,0] \
            , linewidth=1., color='k')

      ax[it_xxx,it_yyy].yaxis.set_major_formatter(formatter)

      if ( (it_xxx==it_yyy) \
            | ( (it_xxx==2) & (it_yyy==1) ) \
            | ( (it_xxx==1) & (it_yyy==0) ) \
            | ( (it_xxx==0) & (it_yyy==1) ) \
            | ( (it_xxx==0) & (it_yyy==3) ) \
            ):
        ipx=0.15
        ipy=0.85
      elif ( (it_xxx==2) & (it_yyy==3) ):
        ipx=0.85
        ipy=0.75
      else:
        ipx=0.85
        ipy=0.85
      if ((it_xxx == it_yyy)&(it_xxx==0)):
        elseed = r'\mathfrak{m}'
      else:
        elseed = r'\tilde{\mathfrak{m}}'
      auto_label_subplots(ax[it_xxx,it_yyy] \
          , r'$%s_{\rm %i,%i}$' % (elseed, it_xxx+1,it_yyy+1) \
          , px=ipx, py=ipy)


      if (itl==0):
        ax[it_xxx, it_yyy].set_yscale('log')
        ax[it_xxx, it_yyy].set_xscale('log')
        if (it_xxx!=it_yyy):
          if ( (it_xxx<2 ) & (it_yyy>1) ):
            ax[it_xxx, it_yyy].set_ylim(1.e-22,1.e-12)
          elif ( (it_yyy<2 ) & (it_xxx>1) ):
            ax[it_xxx, it_yyy].set_ylim(1.e-22,1.e-12)
          else:
            ax[it_xxx, it_yyy].set_ylim(1.e-22,1.e-3)
        else:
          ax[it_xxx, it_yyy].set_ylim(7.e-1,1.e0)
    
      if (it_xxx==3):
        ax[it_xxx, it_yyy].set_xlabel(r'$\lambda$ [$\mu$m]')

if (not exists(opath)):
  mkdir(opath)

pl.savefig("%s/%s" % (opath, oname, ), dpi=600, format='pdf')


