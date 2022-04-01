#/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import matplotlib.pyplot as pl
#pl.ion()
pl.ioff()
pl.rcParams['hatch.linewidth']=0.4
import numpy as np
from scipy.interpolate import interp1d as interpol
from os.path import exists
from os import mkdir

pl.rcParams.update({'font.size':15})
pl.rcParams.update({'text.usetex':True}) 
pl.rcParams['xtick.minor.visible'] = True 
pl.rcParams['ytick.minor.visible'] = True

################################################################################

opath = 'figures_202110_v03'
oname = 'figure_1.pdf'

option = 2
pattern = 1

if (option == 2):

  pl.close(2)
  fg,ax=pl.subplots(num=2,nrows=1,ncols=1,figsize=(6,5))
  pl.subplots_adjust(left=-0.02, bottom=0.01, wspace=0., hspace=0., top=1., right=0.68)


  num = 21
  yl0 = 2.2
  yl1 = 5.4
  #
  xlim = np.linspace(yl0,yl1,num)
  ylim0 = xlim * 0. + 10.
  ylim3 = ylim0 * 0. + yl0
  ylim4 = ylim0 * 0. + 2
  ylim5 = ylim0 * 0. + 1
  ylim6 = ylim0 * 0. - 2.

  
  # Zerodur
  ax.plot(xlim,ylim5, color='k')
  if (pattern > 0):
    ax.fill_between(xlim, ylim6, ylim5, hatch='\\', facecolor='white')
  # Conductor:
  ax.plot(xlim,ylim4, color='k')
  if (pattern > 0):
    ax.fill_between(xlim, ylim5, ylim4, hatch='/', facecolor='white')

  # Oxide:
  ax.plot(xlim,ylim3, color='k')
  if (pattern > 0):
    ax.fill_between(xlim, ylim4, ylim3, hatch='x', facecolor='white')

  # Air:
  if (pattern > 0):
    ax.fill_between(xlim, ylim3, ylim0, hatch='+', facecolor='white')

  # Dust:
  #ax.plot(xlim, ylim2, color='k')
  #ax.plot(xlim,ylim1, color='k')
  #if (pattern > 0):
  #  ax.fill_between(xlim, ylim2, ylim1, hatch='x', facecolor='white')
  #
  # Dust: two cubes with random angles:
  #
  np.random.seed(5)
  ylimd = []
  xl0 = xlim[xlim.size//5]
  xl1 = xlim[xlim.size//7*4]
  y1 = yl1+np.random.rand(num)*0.2
  y2 = yl0+np.random.rand(num)*0.2
  x = xlim / (xlim[-1]-xlim[0]) * (yl1-yl0)
  sgm = 2.e-2
  for itn in range(2):

    xsq = np.hstack([x,np.hstack([0, sgm*np.random.randn(num-2), 0]) \
        + x[-1],x[::-1],np.hstack([0, sgm*np.random.randn(num-2),0]) \
        + x[0]])
    ysq = np.hstack([y1,np.linspace(y1[-1],y2[-1],num) \
        ,y2[::-1],np.linspace(y2[0],y1[0],num)])

    xa = xsq.mean()
    ya = ysq.mean()

    ang = [0.1, -np.pi/7][itn]#np.random.rand() * 2. * np.pi
    ca = np.cos(ang)
    sa = np.sin(ang)
    rot = np.array([ca, sa, -sa, ca]).reshape(2,2)
    res = np.dot(rot, np.vstack([xsq-xa, ysq-ya]))

    xsq = res[0,:] + xa
    ysq = res[1,:] + ya

    if (itn==0):
      xsq += (xl0 - xsq.max())
    else:
      xsq += (xl1 - xsq.min())

    ww = np.where( (xsq>xlim[0]) & (xsq<xlim[-1]) )[0]
    ww = np.where(xsq==xsq)

    ax.plot(xsq[ww], ysq[ww] - ysq[ww].min() + yl0, color='k')
    if (pattern > 0):
      ax.fill_between(xsq, ysq - ysq[ww].min() + yl0, hatch='x', facecolor='white')
  np.random.seed()

  # Labels:
  ax.text(xlim[-1]+0.05, 3.8, r'Dust: $p_{0}(t)$')
  ax.text(xlim[-1]+0.05, 1.9, r'Al. Oxide: $d_{2}(t)$, $n_{2}$')
  ax.text(xlim[-1]+0.05, 1.4, r'Aluminum: $d_{3}(t)$, $\hat{n}_{3}$')
  ax.text(xlim[-1]+0.05, -0.2, r'Zerodur: $d_{4}$, $n_{4}$')
  ax.text(xlim[-1]+0.05, 7.5, r'air: $n_{1}$')

#//  ax.text(1.05, 3.8, r'Dust: $d_{\rm dust}(t)$, $n_{\rm dust}$')
#//  ax.text(1.05, 1.9, r'Al. Oxide: $d_{\rm ox}(t)$, $n_{\rm ox}$')
#//  ax.text(1.05, 1.4, r'Aluminum: $d_{\rm al}(t)$, $\hat{n}_{\rm al}$')
#//  ax.text(1.05, -0.2, r'Zerodur: $d_{\rm zer}$, $n_{\rm zer}$')
#//  ax.text(1.05, 7.5, r'air: $n_{\rm air}$')

#  ax.arrow(0.3, 8.8, 0.2, -2.4 \
#      , length_includes_head=False, width=0.01, color='k', )

  pl.ylim(-1., 9.)
  pl.xlim(xlim[0], xlim[-1])

  ax.set_axis_off()

  #pl.adjust()

  #pl.tight_layout()


if (not exists(opath)):
  mkdir(opath)

pl.savefig("%s/%s" % (opath, oname, ), dpi=600, format='pdf')


  
