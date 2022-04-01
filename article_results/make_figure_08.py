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

opath = 'figures_article'
odirseed='results_article/'
fname = "figure_8.pdf"

dd='neelt'
relative = False

#
#
if (dd=='eelt'):
  telescopes = ['eelt']
else:
  telescopes = ['gtc']
#
methods = []
methods.append('azimuthal')
methods.append('lineal')
methods.append('random')
methods.append('symmetric')
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
  #nums = [120]
  tlong = 799#365*2+6#798. # Simulated days



  deltat=1
  cleandust=1.
  tstep=1.
  tlong=798//2+1
  lambs = [5000.,10000., 24000.]
  nums=[720]
  mltch=1 

  mltch=2

  cleandust = 10000.      # DUST CLEANING TIME FREQUENCY [in days?]
  cleandust = 1.      # DUST CLEANING TIME FREQUENCY [in days?]


else:
  deltat = 30.
  lambs = [5000.]
  nums = [720]
  tstep = 1. # SIMULATION TIME STEP   $NUMBER OF TIMESTEPS FOR EACH MIRROR
  tlong = 370 # Simulated days
  mltch=4 

  deltat=10
  cleandust=1.
  tstep=10.
  tlong=370
  lambs = [5000.]
  lambs = [6500.]
  nums=[720]
  mltch=1 


  deltat=10
  cleandust=1.
  tstep=10.
  tlong=370
  lambs = np.linspace(2100.,24000.,219+1)[:]#(np.arange(7)+1.)*1000.+3000.
#  lambs = np.linspace(3500.,24000.,205+1)[:]#(np.arange(7)+1.)*1000.+3000.
  nums=[720]
  mltch=1 


nums = [4096]
#



omethod = [r'Azimuth', r'Linear', r'Random', r'Symmetry']
x = np.arange(tlong/tstep) * tstep

toplot = np.zeros((4,4,len(lambs), len(methods))) - np.nan

for itt, telescope in enumerate(telescopes):
  
  for itl, lamb in enumerate(lambs):
    
    for itn, num in enumerate(nums):
  
      for itm, method in enumerate(methods):
  
        outdir = '%s/tel%s_%s_alpha%.2f_%08i' % (odirseed, telescope.upper() \
            , method, alpha, num)
        
        # Posibles necesidades: tstep, tlong, cleansdust, period, tlong
        oname = '%s/w%.2f_Dt%.2f_St%.2f_%08i_%010.2f_%i.npz' \
              % (outdir, lamb, deltat, tstep, tlong, cleandust, mltch, )

        if (not exists(oname)):
          continue

        data = np.load(oname)
        int_mat = data['avg_mat']

        tmp = int_mat[:,0,0] * 1.

        int_mat = int_mat / tmp[:,None,None]
        int_mat[:,0,0] = tmp * 1.

        nt, nx, ny = int_mat.shape

        toplot[:,:,itl,itm] = (int_mat.max(0)-int_mat.min(0))
      
    
pl.close(17)
fig, ax = pl.subplots(num=17,nrows=4,ncols=4,figsize=(12/2,8/2),sharex=True)

for it_xxx in range(nx):
  for it_yyy in range(ny):
    for itm, method in enumerate(methods):
      ax[it_xxx, it_yyy].plot(lambs*1.e-4, toplot[it_xxx, it_yyy, :, itm] \
          , linewidth=1., label=omethod[itm])
 
for it_xxx in range(nx):
  for it_yyy in range(ny):
    ax[it_xxx,it_yyy].yaxis.set_major_formatter(formatter)

    ax[it_xxx,it_yyy].set_yscale("log")
    if ((it_xxx == it_yyy)&(it_xxx==0)):
      elseed = r'\Delta\mathfrak{m}'
    else:
      elseed = r'\Delta\tilde{\mathfrak{m}}'
    auto_label_subplots(ax[it_xxx,it_yyy] \
        , r'$%s_{\rm %i,%i}$' % (elseed, it_xxx+1,it_yyy+1) \
        , px=0.85, py=0.85)


    if (it_xxx == 3):
      ax[it_xxx, it_yyy].set_xlabel(r'$\lambda$ [$\mu$m]')

ax[0,3].legend(loc='lower center',bbox_to_anchor=(0.5,0.92) \
    , fancybox=True, shadow=True,ncol=4,  bbox_transform=fig.transFigure)

pl.subplots_adjust(left=0.07,bottom=0.07,right=0.99 \
    ,top=0.91,wspace=0.44,hspace=0.25)


if (not exists(opath)):
  mkdir(opath)

pl.savefig("%s/%s" % (opath, fname, ), dpi=600, format='pdf')




