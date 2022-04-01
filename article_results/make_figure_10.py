#!/usr/bin/env python3.8

import numpy as np
import time as tm

import pymirrors as mrr

#TO BE MOVED WHEN SATISFACTORY RESULTS ARE OBTAINED
from os.path import exists
from os import mkdir
import matplotlib.pyplot as pl
from matplotlib import ticker 
pl.ioff()

pl.rcParams.update({'font.size':15/2})  
pl.rcParams.update({'xtick.labelsize':14/2})
pl.rcParams.update({'ytick.labelsize':14/2})
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

# Main

opath = 'figures_article'
odirseed='results_article/'
fname = 'figure_10.pdf'


dd='neelt'
relative = False

#
#
if (dd=='eelt'):
  telescopes = ['seelt', 'eelt']
  telescopes = ['eelt']
else:
  telescopes = ['sgtc', 'gtc']
  telescopes = ['gtc']


methods = []
methods.append('azimuthal')
methods.append('lineal')
methods.append('random')
methods.append('symmetric')

omethod = [r'Azimuth', r'Linear', r'Random', r'Symmetry']

#
alpha = 0.
x_alpha = 0.
y_alpha = 0.
#

if (dd=='eelt'):
  deltat = 1.
  lambs = [5000.]
  nums = [720]
  tstep = 1. # SIMULATION TIME STEP   $NUMBER OF TIMESTEPS FOR EACH MIRROR
  tlong = 798//2+1
  mltch = 2
  cleandust = 30.      # DUST CLEANING TIME FREQUENCY [in days?]
else:
  deltat = 10. # Days between segment exchange
  tstep = 10. # SIMULATION TIME STEP   $NUMBER OF TIMESTEPS FOR EACH MIRROR
  lambs = [5000.]
  nums = [720]
  tlong = 370 # Simulated days
  mltch = 1
  cleandust = 30.      # DUST CLEANING TIME FREQUENCY [in days?]

nums = [4096]

#


x = np.arange(tlong//tstep) * tstep

for itt, telescope in enumerate(telescopes):
  pl.close(20+itt)
  fig, ax = pl.subplots(num=20+itt,nrows=4,ncols=4,figsize=(12/2,8/2),sharex=True)
  pl.subplots_adjust(left=0.07,bottom=0.07,right=0.99 \
      ,top=0.91,wspace=0.44,hspace=0.25)
  
  for itl, lamb in enumerate(lambs):
    
    for itn, num in enumerate(nums):
  
      for itm, method in enumerate(methods):
  
        outdir = '%s/tel%s_%s_alpha%.2f_%08i' % (odirseed, telescope.upper() \
            , method, alpha, num)
        
        # Posibles necesidades: tstep, tlong, cleansdust, period, tlong
        oname = '%s/w%.2f_Dt%.2f_St%.2f_%08i_%010.2f_%i.npz' % \
              (outdir, lamb, deltat, tstep, tlong, cleandust, mltch, )
        
        data = np.load(oname)
        int_mat = data['avg_mat']
 
        nt, nx, ny = int_mat.shape
      
        for it_xxx in range(nx):
          for it_yyy in range(ny):
            if (relative==False):
              if ((it_xxx==0)&(it_yyy==0)):
                ax[it_xxx, it_yyy].plot(x,int_mat[:, it_xxx, it_yyy], linewidth=1. \
                    , label=omethod[itm])
              else:
                ax[it_xxx, it_yyy].plot(x,int_mat[:, it_xxx, it_yyy] \
                  / int_mat[:, 0, 0], linewidth=1. \
                  , label=omethod[itm])
            else:
              ax[it_xxx, it_yyy].plot(x,int_mat[:, it_xxx, it_yyy]\
                  / int_mat[0,0,0]*100., linewidth=1. \
                  , label=omethod[itm])
      
 
  for it_xxx in range(nx):
    for it_yyy in range(ny):
      ax[it_xxx,it_yyy].yaxis.set_major_formatter(formatter)

      if ((it_xxx == it_yyy)&(it_xxx==0)):
        elseed = r'\mathfrak{m}'
      else:
        elseed = r'\tilde{\mathfrak{m}}'
      auto_label_subplots(ax[it_xxx,it_yyy] \
          , r'$%s_{\rm %i,%i}$' % (elseed, it_xxx+1,it_yyy+1) \
          , px=0.85, py=0.85)

      if (it_xxx == 3):
        ax[it_xxx,it_yyy].set_xlabel(r't [days]')

      if ((it_xxx == it_yyy)&(it_xxx==0)):
        ax[it_xxx, it_yyy].set_ylim(0.765, 0.825)

  ax[0,3].legend(loc='lower center',bbox_to_anchor=(0.5,0.92) \
      , fancybox=True, shadow=True,ncol=4,  bbox_transform=fig.transFigure)
    





if (not exists(opath)):
  mkdir(opath)

pl.savefig("%s/%s" % (opath, fname, ), dpi=600, format='pdf')

















