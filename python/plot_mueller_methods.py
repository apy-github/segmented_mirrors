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
from matplotlib import ticker 
pl.ioff()

pl.rcParams.update({'font.size':15})  
pl.rcParams.update({'xtick.labelsize':14})  
pl.rcParams.update({'ytick.labelsize':14})  
pl.rcParams.update({'text.usetex':True}) 
pl.rcParams.update({'xtick.minor.visible':True })
pl.rcParams.update({'ytick.minor.visible':True})

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
          return r'${:.3f}$'.format(value)
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

opath = 'figures_202110_v03'
oname = 'figure_'
#odirseed = 'results/'
odirseed='results_202110_v03/'

# Main

dd='neelt'
dd='eelt'
relative = False#True#False#True


cdust=False

#
#
if (dd=='eelt'):
  telescopes = ['seelt', 'eelt']
  telescopes = ['eelt']
  oname = "%s7.pdf" % (oname, )
  ylims = [0.810, 0.830]
else:
  telescopes = ['sgtc', 'gtc']
  telescopes = ['gtc']
  oname = "%s5.pdf" % (oname, )
  ylims = [0.824, 0.830]
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
methods = []
methods.append('azimuthal')
methods.append('lineal')
methods.append('random')
methods.append('symmetric')

omethod = [r'Azimuth', r'Linear', r'Random', r'Symmetry']

#
print('Telescope chosen!')
alpha = 0.
x_alpha = 0.
y_alpha = 0.
print('Orientation chosen!')
#
#period = 365 #798.#365.24 * 2.

if (dd=='eelt'):
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

else:
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
#


print('Time parameters chosen!')
#lamb = 3.e3    #ANGSTROM
print('Wavelength chosen!')
print('Number of rays per segment chosen!')


x = np.arange(tlong//tstep) * tstep

for itt, telescope in enumerate(telescopes):
  pl.close(20+itt)
  fig, ax = pl.subplots(num=20+itt,nrows=4,ncols=4,figsize=(12,8),sharex=True)
  pl.subplots_adjust(left=0.07,bottom=0.07,right=0.99 \
      ,top=0.91,wspace=0.44,hspace=0.25)
  
  for itl, lamb in enumerate(lambs):
    
    for itn, num in enumerate(nums):
  
      for itm, method in enumerate(methods):
  
        outdir = '%s/tel%s_%s_alpha%.2f_%08i' % (odirseed, telescope.upper() \
            , method, alpha, num)
        
        # Posibles necesidades: tstep, tlong, cleansdust, period, tlong
        oname2 = '%s/avg_mueller_matrix_w%.2f_Dt%.2f_St%.2f_%08i_%010.2f_%i.fits' % \
              (outdir, lamb, deltat, tstep, tlong, cleandust, mltch, )
        print(oname2)
        int_mat = uts.readfits_v3(oname2, path='./')[1]
        
        nt, nx, ny = int_mat.shape
      
        for it_xxx in range(nx):
          for it_yyy in range(ny):
            if (relative==False):
              ax[it_xxx, it_yyy].plot(x,int_mat[:, it_xxx, it_yyy], linewidth=1. \
                  , label=omethod[itm])
            else:
              ax[it_xxx, it_yyy].plot(x,int_mat[:, it_xxx, it_yyy]\
                  / int_mat[0,0,0]*100., linewidth=1. \
                  , label=omethod[itm])
      
 
  for it_xxx in range(nx):
    for it_yyy in range(ny):
      ax[it_xxx,it_yyy].yaxis.set_major_formatter(formatter)

   #   if ( (it_xxx==it_yyy) \
   #         | ( (it_xxx==2) & (it_yyy==1) ) \
   #         | ( (it_xxx==1) & (it_yyy==0) ) \
   #         | ( (it_xxx==0) & (it_yyy==1) ) \
   #         | ( (it_xxx==0) & (it_yyy==3) ) \
   #         ):
   #     auto_label_subplots(ax[it_xxx,it_yyy] \
   #         , r'M$_{\rm %i,%i}$' % (it_xxx+1,it_yyy+1) \
   #         , px=0.15, py=0.85)
   #   elif ( (it_xxx==2) & (it_yyy==3) ):
   #     auto_label_subplots(ax[it_xxx,it_yyy] \
   #         , r'M$_{\rm %i,%i}$' % (it_xxx+1,it_yyy+1) \
   #         , px=0.85, py=0.75)
   #   else:
   #     auto_label_subplots(ax[it_xxx,it_yyy] \
   #         , r'M$_{\rm %i,%i}$' % (it_xxx+1,it_yyy+1) \
   #         , px=0.85, py=0.85)
      auto_label_subplots(ax[it_xxx,it_yyy] \
          , r'$\mathfrak{m}_{\rm %i,%i}$' % (it_xxx+1,it_yyy+1) \
          , px=0.85, py=0.85)

      if (it_xxx == it_yyy):
        ax[it_xxx, it_yyy].set_ylim(ylims[0], ylims[1])
      if (it_xxx == 3):
        ax[it_xxx,it_yyy].set_xlabel(r't [days]')

  ax[0,3].legend(loc='lower center',bbox_to_anchor=(0.5,0.92) \
      , fancybox=True, shadow=True,ncol=4,  bbox_transform=fig.transFigure)
  #fig.canvas.set_window_title('Telescope: %s' % (telescope, ))
    


if (not exists(opath)):
  mkdir(opath)

pl.savefig("%s/%s" % (opath, oname, ), dpi=600, format='pdf')


















