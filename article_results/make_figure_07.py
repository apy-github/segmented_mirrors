#!/usr/bin/env python3

import numpy as np
import time as tm

#import pymirrors as mrr

from os.path import exists
from os import mkdir
import matplotlib.pyplot as pl
from matplotlib import ticker 
pl.ioff()

pl.rcParams.update({'font.size':15/2.})  
pl.rcParams.update({'xtick.labelsize':14/2.})  
pl.rcParams.update({'ytick.labelsize':14/2.})  
pl.rcParams.update({'text.usetex':False}) 
pl.rcParams.update({'xtick.minor.visible':True })
pl.rcParams.update({'ytick.minor.visible':True})

if (pl.rcParams['text.usetex']):
  pl.rcParams['text.latex.preamble'] = [r"\usepackage{amsfonts}",]

#---------------------------------------------------------------------
#
def auto_label_subplots(iax, label,px=None,py=None):
  #
  if (px==None):
    px=0.9
  if (py==None):
    py=0.9
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


opath = 'figures_article'
odirseed='results_article/'
oname = 'figure_'

# Main

dd='neelt'
#dd='eelt'
relative = False


cdust=False

#
#
if (dd=='eelt'):
  telescopes = ['eelt']
  fname = "%s9.pdf" % (oname, )
  ylims = [0.824, 0.827]
else:
  telescopes = ['gtc']
  fname = "%s7.pdf" % (oname, )
  ylims = [0.824, 0.830]
#
# Checked:
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

  if (cdust==True):
    cleandust = 30.      # DUST CLEANING TIME FREQUENCY [in days?]
  else:
    cleandust = 1.      # DUST CLEANING TIME FREQUENCY [in days?]

  nums = [4096*4]
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
        #int_mat = uts.readfits_v3(oname2, path='./')[1]
        
        nt, nx, ny = int_mat.shape
      
        for it_xxx in range(nx):
          for it_yyy in range(ny):
            if (relative==False):
              if ((it_xxx==0)&(it_yyy==0)):
                ax[it_xxx, it_yyy].plot(x,int_mat[:, it_xxx, it_yyy], linewidth=1. \
                    , label=omethod[itm])
              else:
                ax[it_xxx, it_yyy].plot(x,int_mat[:, it_xxx, it_yyy] \
                    / int_mat[:, 0, 0] \
                    , linewidth=1., label=omethod[itm])
            else:
              ax[it_xxx, it_yyy].plot(x,int_mat[:, it_xxx, it_yyy]\
                  / int_mat[0,0,0]*100., linewidth=1. \
                  , label=omethod[itm])
      
 
  for it_xxx in range(nx):
    for it_yyy in range(ny):
      ax[it_xxx,it_yyy].yaxis.set_major_formatter(formatter)

      if ((it_xxx == it_yyy)&(it_xxx==0)):
        if (pl.rcParams['text.usetex']):
          elseed = r'\mathfrak{m}'
        else:
          elseed = r'm'
      else:
        if (pl.rcParams['text.usetex']):
          elseed = r'\tilde{\mathfrak{m}}'
        else:
          elseed = r'\tilde{m}'

      auto_label_subplots(ax[it_xxx,it_yyy] \
          , r'$%s_{\rm %i,%i}$' % (elseed, it_xxx+1,it_yyy+1) \
          , px=0.85, py=0.85)


      if (it_xxx == it_yyy):
        if ((it_xxx==it_yyy)&(it_xxx!=0)):
          ax[it_xxx, it_yyy].set_ylim(ylims[0]/int_mat[:, 0, 0].max() \
              , ylims[1]/int_mat[:, 0, 0].min())
        else:
          ax[it_xxx, it_yyy].set_ylim(ylims[0], ylims[1])
      if (it_xxx == 3):
        ax[it_xxx,it_yyy].set_xlabel(r't [days]')

  ax[0,3].legend(loc='lower center',bbox_to_anchor=(0.5,0.92) \
      , fancybox=True, shadow=True,ncol=4,  bbox_transform=fig.transFigure)


if (not exists(opath)):
  mkdir(opath)

pl.savefig("%s/%s" % (opath, fname, ), dpi=600, format='pdf')


















