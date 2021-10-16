
import matplotlib.pyplot as pl
pl.ion()
pl.rcParams.update({'font.size':12\
    , 'xtick.labelsize':12\
    , 'ytick.labelsize':12\
    })

import numpy as np
#
#
#
def get_labs(it, ia, ib, ic, ee):

  ea = "%s^{%s}%s" % (ia, it, ee,)
  eb = "%s^{%s}%s" % (ib, it, ee,)
  ec = "%s^{%s}%s" % (ic, it, ee,)

  return ea, eb, ec
#
def plot_eaxes(ax, x, y, kws1, kws2, kws3, pkwe):

  plot_arrow(ax, x, y, arrowhead=0.2 \
      , plotkwargs=pkwe, **kws1)

  plot_arrow(ax, x, y, arrowhead=0.2 \
      , plotkwargs=pkwe, **kws2)

  plot_arrow(ax, x, y, arrowhead=0.2 \
      , plotkwargs=pkwe, **kws3)

  return
#
def get_intersection_plane(cdir, x0, z0, z1):

  a = cdir[2] / cdir[0]

  b = z0 - a * x0

  x1 = (z1 - b) / a

  return x1, z1
#
def get_intersection(cdir, x1, z1, r, k, offset):

  m = cdir[2] / cdir[0]

  b = z1 - m * x1

  cte1 = -2. * r
  cte2 = 1. + k

  ctea = -cte1 * offset + cte2 * offset**2
  cteb = cte1 - 2. * cte2 * offset

  cteq = 1. + cte2 * m**2
  ctel = 2. * cte2 * m * b + cteb * m
  ctei = cte2 * b**2 + cteb * b + ctea

  xs2 = np.sqrt(ctel**2-4.*cteq*ctei)
  xs1 = (-ctel + xs2) / (2. * cteq)
  xs2 = (-ctel - xs2) / (2. * cteq)

  zs1 = get_conicNd(xs1,k,r, offset=offset)
  zs2 = get_conicNd(xs2,k,r, offset=offset)

  if (abs(zs1-offset)<abs(zs2-offset)):
    return xs1, 0., zs1
  else:
    return xs2, 0., zs2

  return

def get_conicNd(x, k, r, offset = 0.):

  dims = np.shape(x)
  if (len(dims)==0):
    x = np.array([x])
  x1d = x.flatten()
  res = x1d * 0.
  for itn in range(x1d.size):
    _, res[itn] = get_conic(x1d[itn], k, r, offset=offset)

  return res.reshape(dims)

def get_conic(x, k, r, offset=0.):
  #
  def get_conic_sols(ia,ib,ic):
    is1=(-ib+np.sqrt(ib**2-4.*ia*ic))/(2.*ia)
    is2=(-ib-np.sqrt(ib**2-4.*ia*ic))/(2.*ia)
    return is1, is2
  #
  c1 = k + 1
  c2 = -2. * r
  c3 = x**2
  #
  if (k==-1):
    s = -c3 / c2
  else:
    s1, s2 = get_conic_sols(c1, c2, c3)
    t1, t2 = get_conic_sols(c1, c2, 0)

    if (abs(t1)<abs(t2)):
      cs = s1 * 1.
    else:
      cs = s2 * 1.
  #
    zs = cs + offset
  #
  return cs, zs

def rotate(dx, dy, dz, rang, axis, rot=1):

  ca = np.cos(rang)
  sa = np.sin(rang)

  mat = np.zeros((3,3,))
  for its in range(3):
    for itf in range(3):
      val = 0.
      if ( (its==itf) & (its==axis) ):
          val = 1.
      if ( (its==itf) & (its!=axis) & (itf!=axis)):
        val = ca * 1.
      if ( (its==itf) & (its!=axis) & (itf!=axis)):
        val = ca * 1.
      if ( (its<itf) & (its!=axis) & (itf!=axis)):
        val = - sa * rot/np.abs(rot)
      if ( (its>itf) & (its!=axis) & (itf!=axis)):
        val = sa * rot/np.abs(rot)
      mat[its, itf] = val * 1.

  res = np.dot(mat, np.array([dx, dy, dz]))

  return res

def get_cosines(dx, dy, dz):

  dd = np.sqrt(dx**2+dy**2+dz**2)

  dx /= dd
  dy /= dd
  dz /= dd

  return dx, dy, dz

def get_normal(x, z, k, r, offset=0., norm=False):
  #
  c1 = k + 1
  c2 = -2. * r
  c3 = x**2
  #
  cs = z - offset
  #
  nx = - 2. * x
  ny = - 0.
  nz = - (c2 + 2 * c1 * cs)
  #
  if (norm):
    nx, ny, nz = get_cosines(nx, ny, nz)
  #
  return nx, ny, nz

def get_circle(x, y, r):

  x = np.linspace(-r, r, 1001)
  xc = r * np.cos(2. * np.pi * x / (2.*r))
  yc = r * np.sin(2. * np.pi * x / (2.*r))

  return xc, yc

def plot_arrow(ax, x0, y0 \
    , x1=np.nan, y1=np.nan, dx=np.nan, dy=np.nan, a0=np.nan, l0=np.nan, dh=np.nan\
    , dz=np.nan\
    , arrowhead=0.02, label='', labelposition='tr'\
    , plotkwargs={}, textkwargs={}):

  if (dz==dz):
    adz = np.abs(dz)
    xcirc, ycirc = get_circle(x0, y0, adz)
    ax.plot(xcirc+x0, ycirc+y0, **plotkwargs)
    if (dz<0):
      dd = adz * np.sin(45./180.*np.pi)
      ax.plot([x0-dd, x0+dd], [y0+dd, y0-dd], **plotkwargs)
      ax.plot([x0-dd, x0+dd], [y0-dd, y0+dd], **plotkwargs)
    else:
      ax.plot(x0, y0, marker='.', markersize=5., **plotkwargs)

    xmax = xcirc.max() + x0
    ymax = ycirc.max() + y0
    xmin = xcirc.min() + x0
    ymin = ycirc.min() + y0

  else:
    #
    #
    # Check inputs:
    if ( (x1!=x1) | (y1!=y1) ):
      if ( (l0!=l0) | (a0!=a0) ):
        if ( (dx!=dx) | (dy!=dy) ):
          raise Exception("Either x1 and y1 are given, dx and dy are given or l0 and a0 are given!")
        else:
          print(dx, dy)
          l0 = np.sqrt(dx**2+dy**2)
          a0 = np.arctan2(dy, dx)
          x1 = x0 + dx
          y1 = y0 + dy
      else:
        dx = l0 * np.cos(a0)
        dy = l0 * np.sin(a0)
        x1 = x0 + dx
        y1 = y0 + dy
    else:
      dx = x1 - x0
      dy = y1 - y0
      l0 = np.sqrt(dx**2+dy**2)
      a0 = np.arctan2(dy, dx)
    #
    #
    # Segment:
    ax.plot([x0, x1],[y0, y1], **plotkwargs)
    #
    #
    # Arrow head:
    if (dh!=dh):
      dh = l0 * arrowhead
    dh2 = dh / 2.
    x2, y2 = x1 + dh2 * np.cos(a0+np.pi/2.), y1 + dh2 * np.sin(a0+np.pi/2.)
    x3, y3 = x1 + dh * np.cos(a0), y1 + dh * np.sin(a0)
    x4, y4 = x1 + dh2 * np.cos(a0-np.pi/2.), y1 + dh2 * np.sin(a0-np.pi/2.)

    ax.plot([x1,x2,x3,x4,x1], [y1,y2,y3,y4,y1], **plotkwargs)
    #
    #
    xmax = np.max([x1,x2,x3,x4,x1])
    ymax = np.max([y1,y2,y3,y4,y1])
    xmin = np.min([x1,x2,x3,x4,x1])
    ymin = np.min([y1,y2,y3,y4,y1])

  # Label:
  if (len(label)!=0):
    if (labelposition=='tr'):
      va='bottom'
      ha='left'
      xt = xmax
      yt = ymax
    elif (labelposition=='tl'):
      va='bottom'
      ha='right'
      xt = xmin
      yt = ymax
    elif (labelposition=='br'):
      va='top'
      ha='left'
      xt = xmax
      yt = ymin
    elif (labelposition=='bl'):
      va='top'
      ha='right'
      xt = xmin
      yt = ymin

    kws = ['verticalalignment', 'horizontalalignment']
    kvs = [va, ha]
    for itk, itv in zip(kws, kvs):
      textkwargs[itk] = itv
    ax.text(xt, yt, label, **textkwargs)
  #
  #
  #
  return

#
#
#
#
eperplab = r'$\hat{{\rm e}}_{\perp}$'
eparlab = r'$\hat{{\rm e}}_{\parallel}$'
elonlab = r'$\hat{{\rm e}}_{k}$'
iperplab = r'$\hat{{\rm e}}_{\perp}'
iparlab = r'$\hat{{\rm e}}_{\parallel}'
ilonlab = r'$\hat{{\rm e}}_{k}'
eendlab = r'$'
#
#
#
#
#
#
x0 = 4.7
z0 = 14.
epos1 = 4.
eposa = 7.
epos2 = 10.

#
# Telescope properties:
#
r1 = 33. #m
k1 = -1.002250 #
#
dm = 14.739410 #m
#
r2 = 3.899678 #m
k2 = -1.504835
#
rf = 1.792969 #m
#
m1r = 6.
m2r = 1.1726/2.
mfr = 1.2363 #m (0.3956 #m ; unvignetted)
zf = -3.4
#
m1n = 171
m1a = np.linspace(-m1r, m1r, m1n)
m1a = m1a[np.abs(m1a) >= 0.936]  * 1.
m2n = 19
m2a = np.linspace(-m2r, m2r, m2n)
#
m1z = get_conicNd(m1a, k1, r1)
m2z = get_conicNd(m2a, k2, r2, offset=dm)
#
#
px0 = -8.   # m
px1 = 8.    # m
pxn = 11
px = np.linspace(px0, px1, pxn)

py0 = zf-0.5   # m
py1 = dm+2.   # m
pyn = 11
py = np.linspace(py0, py1, pyn)
#
#
pkwaxes = {'color':'k','linewidth':0.8,'linestyle':'dashed'}
pkwrays = {'color':'k','linewidth':0.8,'linestyle':'solid'}

pkwe0 = {'color':'darkmagenta','linewidth':1.3,'linestyle':'solid'}
pkwe1 = {'color':'orange','linewidth':1.3,'linestyle':'solid'}
pkwe2 = {'color':'r','linewidth':1.3,'linestyle':'solid'}
pkwe3 = {'color':'g','linewidth':1.3,'linestyle':'solid'}
pkwe4 = {'color':'b','linewidth':1.3,'linestyle':'solid'}
pkwe5 = {'color':'darkmagenta','linewidth':1.3,'linestyle':'solid'}

# Axes: a1, z
#
#fg, ax = pl.subplots(num=4, clear=True, squeeze=False, ncols=1, nrows=3)
#tax = ax[0,0]
#sax = ax[1,0]
#fax = ax[2,0]

a4width = 8.+1./4.
a4height = 11.+17./24.

fnum = 4
ffrac = 0.9

if (0):
  fsize = (a4width/2.*ffrac, a4height*ffrac)

  pl.close(fnum)
  fg = pl.figure(num=fnum, figsize=fsize)
  figdh = 0.0#4
  figx0 = 0.02
  figw = 0.92
  #figfx0 = 
  figfy0 = figdh
  figfh = 0.18
  #figsx0 = 
  figsy0 = figfy0 + figfh + figdh
  figsh = 0.58
  #figtx0 =
  figty0 = figsy0 + figsh + figdh
  figth = 0.25

  tax = pl.axes((figx0,figty0,figw,figth))
  sax = pl.axes((figx0,figsy0,figw,figsh))
  fax = pl.axes((figx0,figfy0,figw,figfh))

else:
  fsize = (a4width*ffrac, a4height/3.*2.*ffrac)

  pl.close(fnum)
  fg = pl.figure(num=fnum, figsize=fsize)
  figdh = 0.04#4
  figw = 0.5
# Focal view
  figfx0 = 0.
  figfy0 = 0.
  figfh = 0.5 - 2. * figdh
# Side view:
  figsx0 = 0.5
  figsy0 = 0
  figsh = 1.
# Top view:
  figtx0 = 0
  figty0 = figfy0 + figfh + 2.*figdh
  figth = 0.5 - 2. * figdh

  tax = pl.axes((figtx0,figty0,figw,figth))
  sax = pl.axes((figsx0,figsy0,figw,figsh))
  fax = pl.axes((figfx0,figfy0,figw,figfh))

  tax.text(0.5, 1.07, 'Optical system top view' \
        , verticalalignment='bottom', horizontalalignment='center'\
        , transform=tax.transAxes)

  sax.text(0.5, 0.955, 'Side view' \
        , verticalalignment='bottom', horizontalalignment='center'\
        , transform=sax.transAxes)

  fax.text(0.5, 1.07, 'Focal plane top view' \
        , verticalalignment='bottom', horizontalalignment='center'\
        , transform=fax.transAxes)


#bfg, ax = pl.subplots(num=1, clear=True, squeeze=False, ncols=1, nrows=1)
#sax = ax[0,0]
#tfg, ax = pl.subplots(num=2, clear=True, squeeze=False, ncols=1, nrows=1)
#tax = ax[0,0]
#ffg, ax = pl.subplots(num=3, clear=True, squeeze=False, ncols=1, nrows=1)
#fax = ax[0,0]


y0 = np.sqrt(-x0**2+m1r**2) * 0.83
#
# Angle:
ang = np.arctan2(y0, x0)
#
m1x, m1y = get_circle(0., 0., m1r)
m2x, m2y = get_circle(0., 0., m2r)


xaxlab = r'$\hat{\rm x}$'
yaxlab = r'$\hat{\rm y}$'
zaxlab = r'$\hat{\rm z}$'
aaxlab = r'$\hat{\rm a}$'


#
# 
# b) frame:
plot_arrow(sax, 0., 0., dx=m1x.max()+0.1, dy=0., dh=0.3, label=aaxlab, labelposition='br' \
    , plotkwargs=pkwaxes)
plot_arrow(sax, 0., 0., dx=0., dy=dm+0.4, dh=0.3, label=zaxlab, labelposition='tl' \
    , plotkwargs=pkwaxes)

sax.set_xlim([px0,px1])
sax.set_xlim([-0.5,px1])
sax.set_ylim([py0,py1])
sax.set_aspect('equal')

#
#
# Mirrors:
sax.plot(m1a[m1a>0], m1z[m1a>0], color='k')
sax.plot(m1a[m1a<0], m1z[m1a<0], color='k')
sax.plot(m2a, m2z, color='k')





#
#
# Incident ray:
#
#
# From above: a
# Mirrors:
tax.plot(m1x, m1y, color='k', linewidth=0.8, linestyle='solid')
tax.plot(m2x, m2y, color='k', linewidth=0.8, linestyle='solid')
tax.set_xlim([min(0., x0)-0.7,max(x0,m1x.max())+1.5])
tax.set_ylim([min(0., y0)-0.7,max(y0,m1y.max())+1.5])
tax.set_aspect('equal')

# Axes:
plot_eaxes(tax, 0, 0 \
    , {'dx':m1r+0.2,'dy':0,'dh':0.3, 'label':xaxlab, 'labelposition':'tr'} \
    , {'dx':0.,'dy':m1r+0.2,'dh':0.3, 'label':yaxlab, 'labelposition':'tr'} \
    , {'dz':0.3,'dh':0.3, 'label':zaxlab, 'labelposition':'bl'} \
    , pkwaxes)

#


d0x = 1.
d0y = 0.
d1x, d1y, _ = rotate(d0x, d0y, 0., np.pi/2., 2)
#><# Perp. (0):
#><plot_arrow(tax, x0, y0, dx=d0x, dy=d0y, arrowhead=0.2 \
#><    , plotkwargs=pkwe0, label=eperplab)
#><# Prop. (0):
#><plot_arrow(tax, x0, y0, dz=-0.2, arrowhead=0.2 \
#><    , plotkwargs=pkwe0, label=elonlab)
#><# Par. (0):
#><plot_arrow(tax, x0, y0, dx=d1x, dy=d1y, arrowhead=0.2 \
#><    , plotkwargs=pkwe0, label=eparlab)

ea, eb, ec = get_labs('0', iperplab, ilonlab, iparlab, eendlab)

# Axis through the point containing the normal and the incident beam
plot_arrow(tax, 0., 0., a0=ang, l0=abs(m1r)*1.4, dh=0.3, label=aaxlab, labelposition='tr' \
    , plotkwargs=pkwaxes)

# Plot angle alpha1:
plot_arrow(tax, x0, y0, a0=np.pi/2., l0=3., arrowhead=0. \
    , plotkwargs={'color':(0.6,0.6,0.6),'linewidth':0.8,'linestyle':'dashed'})
plot_arrow(tax, x0, y0, a0=np.pi+ang, l0=3., arrowhead=0. \
    , plotkwargs={'color':(0.6,0.6,0.6),'linewidth':0.8,'linestyle':'dashed'})

# Arc:
ang0 = np.pi/2.
ang1 = np.pi+ang
larc = 2.7
angs = np.linspace(ang0, ang1, 1001)
xarc = x0 + np.cos(angs) * larc
yarc = y0 + np.sin(angs) * larc
tax.plot(xarc, yarc, color=(0.6,0.6,0.6))
# Text:
frac = 1.2
tax.text(xarc.mean()*(1.+(1-frac)), yarc.mean()*frac, r'$\alpha_{1}$'\
        , verticalalignment='center', horizontalalignment='center'\
        , rotation=180+90.+angs.mean()*180./np.pi \
    )
# Arrow:
arrowhead = 0.3
dang = np.arctan2(arrowhead/2., arrowhead)
dx1 = arrowhead * np.cos(ang1-dang)
dy1 = arrowhead * np.sin(ang1-dang)
plot_arrow(tax, xarc[-1], yarc[-1] \
    , a0=angs[-1]-np.pi/2.+dang, l0=arrowhead, arrowhead=0. \
    , plotkwargs={'color':(0.6,0.6,0.6),'linewidth':0.8,'linestyle':'solid'})
plot_arrow(tax, xarc[-1], yarc[-1] \
    , a0=angs[-1]-np.pi/2.-dang, l0=arrowhead, arrowhead=0. \
    , plotkwargs={'color':(0.6,0.6,0.6),'linewidth':0.8,'linestyle':'solid'})
#

#


plot_eaxes(tax, x0, y0 \
    , {'dx':d0x,'dy':d0y, 'label':ea, 'labelposition':'tr'} \
    , {'dz':-0.2, 'label':eb, 'labelposition':'tr'} \
    , {'dx':d1x, 'dy':d1y, 'label':ec, 'labelposition':'tr'} \
    , pkwe0)

#
#plot_eaxes(tax, x0, y0 \
#    , {'dx':d0x,'dy':d0y, 'label':eperplab, 'labelposition':'tr'} \
#    , {'dz':-0.2, 'label':elonlab, 'labelposition':'tr'} \
#    , {'dx':d1x, 'dy':d1y, 'label':eparlab, 'labelposition':'tr'} \
#    , pkwe0)



# Rotate system:
rang = ang + np.pi/2.
ed0a, ed0b, _ = rotate(d0x, d0y, 0., rang, 2)
ed1a, ed1b, _ = rotate(ed0a, ed0b, 0., np.pi/2., 2)

ea, eb, ec = get_labs('1', iperplab, ilonlab, iparlab, eendlab)
plot_eaxes(tax, x0, y0 \
    , {'dx':ed0a,'dy':ed0b, 'label':ea, 'labelposition':'bl'} \
    , {'dz':-0.25, 'label':eb, 'labelposition':'br'} \
    , {'dx':ed1a, 'dy':ed1b, 'label':ec, 'labelposition':'br'} \
    , pkwe1)
#
#plot_eaxes(tax, x0, y0 \
#    , {'dx':ed0a,'dy':ed0b, 'label':eperplab, 'labelposition':'bl'} \
#    , {'dz':-0.25, 'label':elonlab, 'labelposition':'br'} \
#    , {'dx':ed1a, 'dy':ed1b, 'label':eparlab, 'labelposition':'br'} \
#    , pkwe1)

tax.set_axis_off()



#















a0 = np.sqrt(x0**2+y0**2)
# Primary:
a1 = a0
z1 = get_conicNd(a1, k1, r1)
# Secondary:

# Incident beam:
plot_arrow(sax, a0, z0, x1=a1, y1=z1, arrowhead=0\
    , plotkwargs=pkwrays)



e1a = a0
e1z = 12.

# Normal:
n1a, _, n1z = get_normal(a1, z1, k1, r1, norm=True)
d1a, _, d1z = get_cosines(a1 - a0, 0., z1 - z0)


# e Axes:
ed1a, _, ed1z = rotate(d1a, 0., d1z, -np.pi/2., 1)
ea, eb, ec = get_labs('1', iparlab, iperplab, ilonlab, eendlab)
plot_eaxes(sax, e1a, e1z \
    , {'dx':ed1a,'dy':ed1z, 'label':ea} \
    , {'dz':0.2, 'label':eb} \
    , {'dx':d1a, 'dy':d1z, 'label':ec} \
    , pkwe1)
#plot_eaxes(sax, e1a, e1z \
#    , {'dx':ed1a,'dy':ed1z, 'label':eparlab} \
#    , {'dz':0.2, 'label':eperplab} \
#    , {'dx':d1a, 'dy':d1z, 'label':elonlab} \
#    , pkwe1)




# Incident ray is inversed in sign!
th1 = np.arccos(n1a*(-d1a)+n1z*(-d1z))

plot_arrow(sax, a1, z1, dx=2.*n1a, dy=2.*n1z, arrowhead=0\
    , plotkwargs=pkwaxes)

# Outgoing beam (1):
d2a, _, d2z = rotate(-d1a, 0., -d1z, 2.*th1, 1)


a2, _, z2 = get_intersection([d2a, 0., d2z], a1, z1, r2, k2, dm)

plot_arrow(sax, a1, z1, x1=a2, y1=z2, arrowhead=0\
    , plotkwargs=pkwrays)


# Perp (2):
e2a = a1 + epos1 * d2a
e2z = z1 + epos1 * d2z
ed2a, _, ed2z = rotate(d2a, 0., d2z, -np.pi/2., 1)

ea, eb, ec = get_labs('2', iparlab, iperplab, ilonlab, eendlab)
plot_eaxes(sax, e2a, e2z \
    , {'dx':ed2a,'dy':ed2z, 'label':ea} \
    , {'dz':-0.2, 'label':eb,'labelposition':'bl'} \
    , {'dx':d2a, 'dy':d2z, 'label':ec} \
    , pkwe2)
#
#plot_eaxes(sax, e2a, e2z \
#    , {'dx':ed2a,'dy':ed2z, 'label':eparlab} \
#    , {'dz':-0.2, 'label':eperplab,'labelposition':'bl'} \
#    , {'dx':d2a, 'dy':d2z, 'label':elonlab} \
#    , pkwe2)


# Plot angle alpha2:
eaa = a1 + eposa * d2a
eaz = z1 + eposa * d2z
# Arc:
ang0 = ang+0.1
ang1 = np.pi+(ang-0.1)
larc = 0.7
angs = np.linspace(ang0, ang1, 1001)
xarc = eaa + np.cos(angs) * larc
yarc = eaz + np.sin(angs) * larc
sax.plot(xarc, yarc, color=(0.6,0.6,0.6))
# Arrow:
arrowhead = 0.3
dang = np.arctan2(arrowhead/2., arrowhead)
dx1 = arrowhead * np.cos(ang1-dang)
dy1 = arrowhead * np.sin(ang1-dang)
plot_arrow(sax, xarc[-1], yarc[-1] \
    , a0=angs[-1]-np.pi/2.+dang, l0=arrowhead, arrowhead=0. \
    , plotkwargs={'color':(0.6,0.6,0.6),'linewidth':0.8,'linestyle':'solid'})
plot_arrow(sax, xarc[-1], yarc[-1] \
    , a0=angs[-1]-np.pi/2.-dang, l0=arrowhead, arrowhead=0. \
    , plotkwargs={'color':(0.6,0.6,0.6),'linewidth':0.8,'linestyle':'solid'})
#
# Arc (b):
ang0 = np.pi+(ang+0.1)
ang1 = 2.*np.pi + (ang-0.1)
larc = 0.7
angs = np.linspace(ang0, ang1, 1001)
xarc = eaa + np.cos(angs) * larc
yarc = eaz + np.sin(angs) * larc
sax.plot(xarc, yarc, color=(0.6,0.6,0.6))
# Arrow:
arrowhead = 0.3
dang = np.arctan2(arrowhead/2., arrowhead)
dx1 = arrowhead * np.cos(ang1-dang)
dy1 = arrowhead * np.sin(ang1-dang)
plot_arrow(sax, xarc[-1], yarc[-1] \
    , a0=angs[-1]-np.pi/2.+dang, l0=arrowhead, arrowhead=0. \
    , plotkwargs={'color':(0.6,0.6,0.6),'linewidth':0.8,'linestyle':'solid'})
plot_arrow(sax, xarc[-1], yarc[-1] \
    , a0=angs[-1]-np.pi/2.-dang, l0=arrowhead, arrowhead=0. \
    , plotkwargs={'color':(0.6,0.6,0.6),'linewidth':0.8,'linestyle':'solid'})
#

# Text:
frac = 1.13
sax.text(eaa*frac, eaz*frac, r'$\alpha_{2}$'\
        , verticalalignment='center', horizontalalignment='left'\
        , rotation=0 \
    )















e3a = a1 + epos2 * d2a
e3z = z1 + epos2 * d2z
ed3a, _, ed3z = rotate(d2a, 0., d2z, np.pi/2., 1)

ea, eb, ec = get_labs('3', iparlab, iperplab, ilonlab, eendlab)
plot_eaxes(sax, e3a, e3z \
    , {'dx':ed3a,'dy':ed3z, 'label':ea} \
    , {'dz':0.2, 'label':eb} \
    , {'dx':d2a, 'dy':d2z, 'label':ec} \
    , pkwe3)


n2a, _, n2z = get_normal(a2, z2, k2, r2, offset=dm, norm=True)

plot_arrow(sax, a2, z2, dx=-2.*n2a, dy=-2.*n2z, arrowhead=0\
    , plotkwargs=pkwaxes)

# Incident (2) and normal (2) ray are inversed in sign!
th2 = np.arccos((-n2a)*(-d2a)+(-n2z)*(-d2z))

# Outgoing beam (2):
d3a, _, d3z = rotate(-d2a, 0., -d2z, -2.*th2, 1)
#
# On the focal plane:
a3, z3 = get_intersection_plane([d3a, 0., d3z], a2, z2, zf)

#plot_arrow(sax, a2, z2, dx=100.*d3a, dy=100.*d3z, arrowhead=0\
plot_arrow(sax, a2, z2, x1=a3, y1=z3, arrowhead=0\
    , plotkwargs=pkwrays)


e4d = 10.
e4a = a2 + e4d * d3a
e4z = z2 + e4d * d3z
ed4a, _, ed4z = rotate(d3a, 0., d3z, np.pi/2., 1)

ea, eb, ec = get_labs('4', iparlab, iperplab, ilonlab, eendlab)
plot_eaxes(sax, e4a, e4z \
    , {'dx':ed4a,'dy':ed4z, 'label':ea} \
    , {'dz':0.2, 'label':eb} \
    , {'dx':d3a, 'dy':d3z, 'label':ec} \
    , pkwe4)

sax.set_axis_off()
#
#
#
#
#
#
# From above c)

xf, yf, _ = rotate(a3, 0., 0., ang, 2)
fax.set_aspect('equal')

dfa = -1.
df1x, df1y, _ = rotate(dfa, 0., 0., rang, 2)
df2x, df2y, _ = rotate(df1x, df1y, 0., np.pi/2., 2)

ea, eb, ec = get_labs('4', iperplab, ilonlab, iparlab, eendlab)
plot_eaxes(fax, xf, yf \
    , {'dx':df1x,'dy':df1y, 'label':ea} \
    , {'dz':-0.2, 'label':eb} \
    , {'dx':df2x, 'dy':df2y, 'label':ec} \
    , pkwe4)


# Return to original setting:


dfc1x, dfc1y, _ = rotate(df1x, df1y, 0., np.pi/2.-ang, 2)


# Arc:
ang0 = ang
ang1 = np.pi / 2.
# Plot angle alpha2:
plot_arrow(fax, xf, yf, a0=ang0, l0=3., arrowhead=0. \
    , plotkwargs={'color':(0.6,0.6,0.6),'linewidth':0.8,'linestyle':'dashed'})
plot_arrow(fax, xf, yf, a0=ang1, l0=3., arrowhead=0. \
    , plotkwargs={'color':(0.6,0.6,0.6),'linewidth':0.8,'linestyle':'dashed'})
larc = 2.7
angs = np.linspace(ang0, ang1, 1001)
xarc = xf + np.cos(angs) * larc
yarc = yf + np.sin(angs) * larc
fax.plot(xarc, yarc, color=(0.6,0.6,0.6))
# Text:
frac = 1.2
fax.text(xarc.mean()*frac, yarc.mean()*frac, r'$\alpha_{3}$'\
        , verticalalignment='top', horizontalalignment='center'\
        , rotation=180+90.+angs.mean()*180./np.pi \
    )
# Arrow:
arrowhead = 0.3
dang = np.arctan2(arrowhead/2., arrowhead)
dx1 = arrowhead * np.cos(ang1-dang)
dy1 = arrowhead * np.sin(ang1-dang)
plot_arrow(fax, xarc[-1], yarc[-1] \
    , a0=angs[-1]-np.pi/2.+dang, l0=arrowhead, arrowhead=0. \
    , plotkwargs={'color':(0.6,0.6,0.6),'linewidth':0.8,'linestyle':'solid'})
plot_arrow(fax, xarc[-1], yarc[-1] \
    , a0=angs[-1]-np.pi/2.-dang, l0=arrowhead, arrowhead=0. \
    , plotkwargs={'color':(0.6,0.6,0.6),'linewidth':0.8,'linestyle':'solid'})
#


dfc2x, dfc2y, _ = rotate(dfc1x, dfc1y, 0., np.pi/2., 2)

ea, eb, ec = get_labs('5', iperplab, ilonlab, iparlab, eendlab)
plot_eaxes(fax, xf, yf \
    , {'dx':dfc1x,'dy':dfc1y, 'label':ea} \
    , {'dz':-0.2, 'label':eb, 'labelposition':'tl'} \
    , {'dx':dfc2x, 'dy':dfc2y, 'label':ec} \
    , pkwe5)


# Axes:
plot_eaxes(fax, 0, 0 \
    , {'dx':2.5,'dy':0,'dh':0.3, 'label':xaxlab, 'labelposition':'tr'} \
    , {'dx':0.,'dy':2.5,'dh':0.3, 'label':yaxlab, 'labelposition':'tr'} \
    , {'dz':0.3,'dh':0.3, 'label':zaxlab, 'labelposition':'bl'} \
    , pkwaxes)
fax.set_axis_off()

