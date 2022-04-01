#
import time as tm
import numpy as np
from pdb import set_trace as stop
import matplotlib.pyplot as pl
pl.ion()
import sys
#
import cython
from numpy cimport ndarray as ar
#
cdef extern from "mirror_calculation.hpp":
    cdef void eval_segment "mrr::eval_mueller_segment<long,int,double>"( \
        long nrays, int nlayers, double* i1, double* th, double* i2 \
        , double lamb, double* ln, double* lk, double* ld, int nthreads \
        , double* dat)

    cdef void inside_circle "mrr::eval_pts_circle<long,double>"( \
        int nthreads, int inside, double radius, long nrays \
        , double* xc, double* yc, long* result)

    cdef void inside_polygon "mrr::eval_pts_polygon<int,long,double>"( \
        int nthreads, int inside, int nverts, double* xverts \
        , double* yverts, long nrays, double* xc, double* yc, long* result)
#
def init_telescope(name, num):
#
  class telescope(object):

    def __init__(self, name, num):
      if (name.upper() == 'GTC'):
        self.ID = 'GTC'
        self.rc_1=33.
        self.k_1=-1.002250
        self.df_1=self.rc_1/2.
                #distancia focal en metros
        self.e_1=np.sqrt(-self.k_1)

        self.rc_2=-3.899678
        self.nrc_2=3.899678
        self.k_2=-1.504835
        self.df_2=self.rc_2/2.
        self.ndf_2=self.nrc_2/2.
        self.e_2=np.sqrt(-self.k_2)

        self.dmm=14.739410
#FOR PRIMAR MIRROR:
        self.num=num
        self.radius=0.936
        self.num_esp=1*36
#AQO
        self.dim_x=6.
        self.dim_y=6.
        self.rad_max=6.
        self.rad_min=0.2
        self.primary_shape='hexagonal'
#FOR SECONDARY MIRROR:
        self.sec_num = 1024
        self.sec_num_esp=36
        self.sec_ext_radius = 0.096714416e0                    #(1172.6/(14.e0*np.cos(30.e0/1.8e2*np.pi)))*1.e-3
        self.sec_int_radius = 0.1388e0
        self.secondary_shape='hexagonal'
        #self.sec_num_esp=1
        #self.sec_ext_radius = 0.65#0.590#19#0.096714416e0                    #(1172.6/(14.e0*np.cos(30.e0/1.8e2*np.pi)))*1.e-3
        #self.sec_int_radius = 0.09#0.1388e0
        #self.secondary_shape='anular'

      elif (name.upper() == 'SGTC'):
        self.ID = 'SGTC'
        self.rc_1 = 33.
        self.k_1 = -1.002250
        self.df_1 = self.rc_1/2.
                #distancia focal en metros
        self.e_1=np.sqrt(-self.k_1)

        self.rc_2=-3.899678
        self.nrc_2=3.899678
        self.k_2=-1.504835
        self.df_2=self.rc_2/2.
        self.ndf_2=self.nrc_2/2.
        self.e_2=np.sqrt(-self.k_2)

        self.dmm=14.739410
#FOR PRIMAR MIRROR:
        self.num=num
        self.radius=0.936
        #self.radius=5.2
        self.num_esp=36
        #self.num_esp=1
        self.dim_x=6.
        self.dim_y=6.
        self.rad_max=6.
        self.rad_min=0.2
        self.primary_shape='hexagonal'
#FOR SECONDARY MIRROR:
        self.sec_num_esp=1
        self.sec_num = 1024
        #self.sec_ext_radius = 0.8#0.500#19#0.096714416e0                    #(1172.6/(14.e0*np.cos(30.e0/1.8e2*np.pi)))*1.e-3
        self.sec_ext_radius = 0.5#0.500#19#0.096714416e0                    #(1172.6/(14.e0*np.cos(30.e0/1.8e2*np.pi)))*1.e-3
        self.sec_int_radius = 0.12#0.11#0.1388e0
        self.secondary_shape='anular'

      elif (name.upper() == 'AGTC'):
        self.ID = 'AGTC'
        self.rc_1 = 33.
        self.k_1 = -1.002250
        self.df_1 = self.rc_1/2.
                #distancia focal en metros
        self.e_1=np.sqrt(-self.k_1)

        self.rc_2=-3.899678
        self.nrc_2=3.899678
        self.k_2=-1.504835
        self.df_2=self.rc_2/2.
        self.ndf_2=self.nrc_2/2.
        self.e_2=np.sqrt(-self.k_2)

        self.dmm=14.739410
#FOR PRIMAR MIRROR:
        self.num=num
        self.radius=0.936
        #self.radius=5.2
        self.num_esp=1
        #self.num_esp=1
        self.dim_x=6.
        self.dim_y=6.
        self.rad_max=5.2
        self.rad_min=0.2
        self.primary_shape='anular'
#FOR SECONDARY MIRROR:
        self.sec_num_esp=1
        self.sec_num = 1024
        self.sec_ext_radius = 0.519#0.096714416e0                    #(1172.6/(14.e0*np.cos(30.e0/1.8e2*np.pi)))*1.e-3
        self.sec_int_radius = 0.11#0.1388e0
        self.secondary_shape='anular'

      elif (name.upper() == 'EELT'):
        self.ID = 'EELT'
        diam_1 = 39.3
        fr_1 = 0.88
        self.rc_1 = diam_1*fr_1*2. # 69 m comes from page 94 in eelt construction proposal
        self.k_1 = -0.993295
        self.df_1 = self.rc_1/2.                  #distancia focal en metros
        self.e_1 = np.sqrt(-self.k_1)

        self.rc_2 = -9.313
        self.nrc_2 = 9.313
        self.k_2 = -2.28962
        self.df_2 = self.rc_2/2.
        self.ndf_2 = self.nrc_2/2.
        self.e_2 = np.sqrt(-self.k_2)

        self.dmm = 31.415
                  #37.742d0/(42.d0*0.93)*39.3d0*0.88d0 
        #self.dmm = 37.742/(42.0*0.93)*39.3*0.88 
#FOR PRIMAR MIRROR:
        self.num = num
        self.radius = 1.45/2.
        self.num_esp = 798
        self.dim_x = 20.
        self.dim_y = 20.
        self.rad_max = 20.
        self.rad_min = 5.4
        self.primary_shape='hexagonal'
#FOR SECONDARY MIRROR:
        self.sec_num_esp=1
        self.sec_num = 451
        self.sec_int_radius = 1.091/2.#0.891e0/2.e0
        self.sec_ext_radius = 4.0906/2.#4.250e0/2.e0
        self.secondary_shape='anular'

      elif (name.upper() == 'SEELT'):
        self.ID = 'SEELT'
        diam_1 = 39.3
        fr_1 = 0.88
        self.rc_1 = diam_1*fr_1*2.
        self.k_1 = -0.993295
        self.df_1 = self.rc_1/2.                  #distancia focal en metros
        self.e_1 = np.sqrt(-self.k_1)

        self.rc_2 = -9.313
        self.k_2 = -2.28962
        self.df_2 = self.rc_2/2.
        self.e_2 = np.sqrt(-self.k_2)

        self.dmm = 31.415
                  #37.742d0/(42.d0*0.93)*39.3d0*0.88d0 
#FOR PRIMAR MIRROR:
        self.num = num
        self.radius = 1.45/2.
        self.num_esp = 798
        self.dim_x = 20.
        self.dim_y = 20.
        self.rad_max = 20.
        self.rad_min = 5.4
        self.primary_shape='hexagonal'
#FOR SECONDARY MIRROR:
        self.sec_num_esp=1
        self.sec_num = 4096#451
        self.sec_int_radius = 1.091/2.#0.891e0/2.e0
        self.sec_ext_radius = 4.0906/2.-0.1#4.250e0/2.e0
        self.secondary_shape='anular'

      elif (name.upper() == 'AEELT'):
        self.ID = 'AEELT'
        diam_1 = 39.3
        fr_1 = 0.88
        self.rc_1 = diam_1*fr_1*2.
        self.k_1 = -0.993295
        self.df_1 = self.rc_1/2.                  #distancia focal en metros
        self.e_1 = np.sqrt(-self.k_1)

        self.rc_2 = -9.313
        self.k_2 = -2.28962
        self.df_2 = self.rc_2/2.
        self.e_2 = np.sqrt(-self.k_2)

        self.dmm = 31.415
                  #37.742d0/(42.d0*0.93)*39.3d0*0.88d0 
#FOR PRIMAR MIRROR:
        self.num = num
        self.radius = diam_1/2. #1.45/2.
        self.num_esp = 1 #798
        self.dim_x = 20.
        self.dim_y = 20.
        self.rad_max = 20.
        self.rad_min = 5.4
        self.primary_shape='anular'
#FOR SECONDARY MIRROR:
        self.sec_num_esp=1
        self.sec_num = 4096#451
        self.sec_int_radius = 1.091/2.#0.891e0/2.e0
        self.sec_ext_radius = 4.0906/2.-0.1#4.250e0/2.e0
        self.secondary_shape='anular'

      elif (name.upper() == 'SYMMETRY'):
        print('Not yet')
      else:
        print("only GTC or EELT are available!")
        sys.exit()

      return

  avail = ['EELT', 'GTC', 'SGTC', 'AGTC', 'SEELT', 'AEELT']

  if (not (name.upper() in avail)):
    print("only GTC or EELT are available!")
    sys.exit()

  teles = telescope(name, num)

  return teles
#
#                                                                              #
#
def init_beam(alpha,x_alpha,y_alpha,degrees=False):

  class beam(object):

    def __init__(self, alp, x_alp, y_alp):

      self.alp = alp * 1.
      self.a_alp = np.arctan2(y_alp,x_alp)
      self.vx_alp = np.sin(alp)*np.cos(self.a_alp)
      self.vy_alp = np.sin(alp)*np.sin(self.a_alp)
      self.vz_alp = np.cos(alp)
      self.l_inc = - self.vx_alp * 1.
      self.m_inc = - self.vy_alp * 1.
      self.n_inc = - self.vz_alp * 1.

      return

  if (degrees != False):
    alp = alpha / 180. * np.pi
  else:
    alp = alpha * 1.

  beams = beam(alp,x_alpha,y_alpha)

  return beams
#
#                                                                              #
#
def primary(telescope, order, tstep, tlong, cleandust \
    , period=0., deltat=0. ,SimulatedTime=-np.inf \
    , ideal=False, multiplechange=1):
#
  class primary_obj(object):

    def __init__(self,telescope,tstep,tlong,cleandust,period,deltat \
            , mltch):

        def get_primary_geometry(inum_esp \
            , iradius, irad_min, irad_max, iprimary_shape):

          if (iprimary_shape=='hexagonal'):
            #
            #GET AN INITIAL ESTIMATION OF MIRRORS:
            #
              ones = np.ones(inum_esp, dtype=np.float64)
              x_odds = np.arange(inum_esp) - inum_esp / 2
              _, x_odds = np.meshgrid(ones, x_odds)
              y_odds = np.arange(inum_esp) - inum_esp / 2
              y_odds, _ = np.meshgrid(y_odds, ones)
              x_even = (np.arange(inum_esp) - \
                  inum_esp / 2) * 2. - 1.
              _, x_even = np.meshgrid(ones, x_even)
              y_even = (np.arange(inum_esp) - \
                  inum_esp / 2) * 2. - 1.
              y_even, _ = np.meshgrid(y_even, ones)
            #
            #GIVE SPATIAL DIMENSIONS TO THE POSITIONS:
            #
              ap = iradius * np.cos(30./180.*np.pi)
              x_fact_odds = 3. * iradius
              y_fact_odds = 2. * ap
              x_fact_even = 1.5 * iradius
              y_fact_even = 1.* ap
              x_odds = x_odds.flatten() * x_fact_odds
              y_odds = y_odds.flatten() * y_fact_odds
              x_even = x_even.flatten() * x_fact_even
              y_even = y_even.flatten() * y_fact_even
            #
              rads_odds = np.sqrt(x_odds**2 + y_odds**2)
              rads_even = np.sqrt(x_even**2 + y_even**2)
            #
            #SORT:
            #
              sort = np.argsort(rads_odds)
              x_odds = x_odds[sort]
              y_odds = y_odds[sort]
              rads_odds = rads_odds[sort]
              sort = np.argsort(rads_even)
              x_even = x_even[sort]
              y_even = y_even[sort]
              rads_even = rads_even[sort]
            #
            #CONCATENATE:
            #
              x = np.concatenate([x_odds, x_even])
              y = np.concatenate([y_odds, y_even])
              rads = np.concatenate([rads_odds, rads_even])
              sort =  np.argsort(rads)
              x = x[sort]
              y = y[sort]
              rads = rads[sort]
            #
            #  ww = np.where(rads >= irad_min)[0]
            #  rads = rads[ww]
            #  x = x[ww]
            #  y = y[ww]
            #
              ww = np.where( (rads >= irad_min) & (rads < irad_max) )[0]
              rads = rads[ww]
              x = x[ww]
              y = y[ww]
#
              if (inum_esp!= rads.size):
                #print("")
                #print("\t Warning!, number of segments inside maximum radius is too big")
                #print(" %i -> %i" % (inum_esp,rads.size,))
                #print("")

                if (inum_esp<rads.size):
                  ww = np.argsort(rads)
                  rads = rads[ww[0:inum_esp]]
                  x = x[ww[0:inum_esp]]
                  y = y[ww[0:inum_esp]]
                else:
                  inum_esp=np.min([inum_esp, rads.size])
  

              # 
              # Area:
              # 
              area = x * 0.
              b = iradius * 1.
              h = ap * 1.
              area = area + b * h / 2. * 6.
              print(" Area covered by the primary mirror: %.4f sqr. meters." % (area.sum(),))

              #
              #FAMILIES:
              #
              family = np.int64(x * 0)
              condit = 0
              num = 0
              step = 1
              Nperfamily = []
              while (condit == 0):
                diff = rads - rads[num]
                wwn = np.where(diff > 0.001)[0]
                #wwn = np.where(rads > rads[num])[0]
                ww = np.where(np.abs(diff) < 0.001)[0]
                #ww = np.where(rads == rads[num])[0]
      # This is not right. Now, once we now what segments are at the same distance...
      # ... we have to look for the azimuth values:
                azis = (np.arctan2(y,x)[ww]*180./np.pi) % 60.
                wws = np.argsort(azis)
                azis = azis[wws]
      
                condit2 = 0
                num2 = 0
                stcnt = 0
                step2 = np.zeros((azis.size, ), dtype=np.int32)
                while (condit2 == 0):
      
                  diff2 = np.abs(azis - azis[num2])
                  wwn2 = np.where( (diff2 > 0.001) & (diff2 < 59.) )[0]
                  ww2 = np.where( (diff2 < 0.001) | (diff2 > 59.) )[0]
      
                  num2 = num2 + np.size(ww2)
      
                  step2[ww2] = stcnt * 1
                  stcnt = stcnt + 1
      
                  if (num2 >= azis.size):
                    condit2 = 1
      
                num = num + np.size(ww)
                for it_nnn in range(np.max(step2)+1):
                  ww2 = np.where(np.abs(step2-it_nnn) < 1.e-3)[0]
                  family[ww[wws][ww2]] = step+it_nnn
                  Nperfamily.append(np.size(ww[wws][ww2]))
                  if (np.size(ww[wws][ww2])==0):
                    stop()
      
                step = step + 1 + np.max(step2)
      
                if (np.size(wwn) == 0):
                  condit = 1
          elif (iprimary_shape=='anular'):
              x = np.zeros((1,), dtype=np.float64)
              y = np.zeros((1,), dtype=np.float64)
              rads = np.zeros((1,)) + iradius
              family = np.zeros((1,),dtype=np.int64)
              Nperfamily = [1]

          else:
            print('primary:')
            print('\tUnknown primary mirror shape!')
            print('')
            stop()

          return x, y, rads, family, Nperfamily, inum_esp


        x,y,rads,family,Nperfamily,telescope.num_esp=\
            get_primary_geometry(telescope.num_esp \
            , telescope.radius, telescope.rad_min, telescope.rad_max \
            , telescope.primary_shape)
      #
        self.xpos = x * 1.
        self.ypos = y * 1.
        self.dist = rads * 1.
        self.azimuth = np.arctan2(y, x)#+self.dist/np.max(self.dist)*0.01
        self.families = family * 1
        self.Nfamily = np.max(family)
        self.Nperfamily = np.array(Nperfamily) * 1
        #
        # Initialize various variables:
        self.sort = self.xpos * 0. - np.nan
        self.N = np.int64(np.size(x))
        #
        if (type(mltch)!=int):
          print("")
          print("\tMultiplechange MUST be an integer!")
          print("")
          sys.exit()
        #
        # May be I should change this limitation:
        if (mltch<1):
          print("")
          print("\tMultiplechange MUST be >=1!")
          print("")
          sys.exit()
        if (self.N%mltch!=0):
          print("")
          print("\tMultiplechange MUST give MODULUS(N,MLT)==0")
          print("")
          sys.exit()
        self.mltch = mltch * 1
        #
# If we change the way of initializing incident rays, it might not be useful

        self.i = {}
        self.i2 = {}
        self.th1 = {}
        self.area = {}
        self.rad = {}

        for i in range(telescope.num_esp):
          self.i['%i' % (i,)]=np.array([])
          self.i2['%i' % (i,)]=np.array([])
          self.th1['%i' % (i,)]=np.array([])
          self.area['%i' % (i,)]=np.array([])
          self.rad['%i' % (i,)]=np.array([])

        self.d11d = 0.
        self.d21d = 0.

      #TIMING:
        # We are given the following:
        #
        #   tstep: elapsed time (days) in each simulated step
        #          if 21, it means, the simulation would have
        #          one point every 21 days.
        #
        #   tsim: total amount (days) of simulated time.
        #          if smaller than tstep, only the first time
        #          step is simulated. In general, the number
        #          of simulated points would be period/tstep
        #
        #   period/deltat:
        #          -period: time (days) in which all the mirrors
        #           have to be changed
        #          -deltat: we change a mirror every this amount
        #           of time (days) 
        #
        #   self.mltch: It says how many mirrors we do change...
        #          ... everytime we change one
        #          Hints:
        #               * It is as if the mirror has N/mltch effective...
        #                 ... segments in order to estimate the period
        #
        #
        #
        #
        if (period!=0):
          #ideltat = (1.*period)/(1.*self.N)
# NEW; MLTCH::
          ideltat = (1.*period)/(1.*self.N/(1.*self.mltch))
          iperiod = period * 1.
        if (deltat!=0):
          iperiod = deltat * (self.N * 1.) / (1.*self.mltch)
          ideltat = deltat * 1.
        #



#AQI






        self.tstep = tstep * 1. # SIMULATION STEP (IN DAYS)
        self.period = iperiod * 1.  # NUMBER OF DAYS FOR SUBSTITUTING ALL MIRRORS
        self.tlong = tlong * 1.  # NUMBER OF DAYS TO BE SIMULATED
        self.deltat = ideltat * 1. # Days in between mirror exchanges
        #
        # Number of simulated steps:
        self.Ntimes = np.int64(np.ceil(self.tlong / self.tstep))
        #
        self.days_between_change = self.deltat * 1.#self.period / (1. * self.N)
        self.days_for_each_timestep = self.deltat / self.tstep #self.tstep * 1.















        # Caution, this name is missleading, it is  not days for each timestep...
        # ... but rather simulation steps in between mirror changes!
        if (self.days_for_each_timestep < 1.):
          print('At least one time step per change is needed, returning nan')
          sys.exit()
        #
        self.SimulatedTimes = self.Ntimes * 1
        self.SimulatedTime = self.tlong * 1
        print('\tNumber of simulated steps: %i' % (self.Ntimes, ))
        print('\tNumber of simulated steps in between mirror changes: %i' \
            % (self.days_for_each_timestep, ))

# DOES IT WORK AT ALL?
        self.time_map = np.zeros((self.N, self.SimulatedTimes, 3), dtype=np.float64)
        # Caution, this variable (cleandustcadence) comes in days...
        # ... and it has to be transformed to simultion steps:
        self.cleandustcadence = cleandust / self.tstep

      #
        return
      #
    def get_order(self, order):

      if (order.upper() == 'RANDOM'):
        sort = np.zeros((self.N,), dtype=np.int32)
        np.random.seed(1)
        fsort = np.argsort(np.random.randn(self.Nfamily))+1
        for it_nnn in range(self.Nfamily):
          ww=np.where(np.abs(self.families-fsort[it_nnn])<1.e-3)[0]
          sort[ww] = (np.argsort(np.random.randn(6))) * self.Nfamily + it_nnn
        sort=np.argsort(sort)
        np.random.seed()

      elif (order.upper() == 'SYMMETRIC'):
        sort = np.int16(np.array([]))
        ref_pos = np.int16(np.arange(self.N))
        for it_nnn in range(self.Nfamily):
          ww = np.where(self.families == (it_nnn + 1))[0]
          it_x = self.xpos[ww]
          it_y = self.ypos[ww]
          it_a = self.azimuth[ww]
          it_sort = np.argsort(it_a)
          #
          it_index = list(ref_pos[ww][it_sort])
          azi360 = list(it_a[it_sort]*180./np.pi)
          #
          if (it_nnn == 0):
            init_val = 0
          else:
            dif = (azi360 - ref_val + 360) % 360
            adif = np.abs(dif - 120.)
            init_val = np.where(adif == np.min(adif))[0][0]
          #
          iff_azi = np.array([azi360[init_val]])
          itf_sort = np.array([it_index[init_val]])
          it_index.remove(it_index[init_val])
          azi360.remove(azi360[init_val])
          #
          for it_nn2 in range(np.size(ww)-1):
            it2_a = iff_azi[-1]
            dif = (azi360 - it2_a + 360) % 360
            adif = np.abs(dif - 120.)
            ww2 = np.where(adif == np.min(adif))[0][0]
            itf_sort = np.hstack([itf_sort, it_index[ww2]])
            iff_azi = np.hstack([iff_azi, azi360[ww2]])
            #
            it_index.remove(it_index[ww2])
            azi360.remove(azi360[ww2])
          #
          ref_val = iff_azi[-1] * 1.
          #
          sort = np.hstack([sort, itf_sort])

      elif (order.upper() == 'LINEAL'):
        
        sort = np.int32(np.array([]))
        ref_pos = np.int32(np.arange(self.N))
        for it_nnn in range(self.Nfamily):
          ww = np.where(self.families == (it_nnn + 1))[0]
          it_x = self.xpos[ww]
          it_y = self.ypos[ww]
          it_a = self.azimuth[ww]
          it_sort = np.argsort(it_a)
          sort = np.hstack([sort, ref_pos[ww][it_sort]])

      elif (order.upper() == 'AZIMUTHAL'):
        sort = np.argsort(self.azimuth)[::-1]

        off = 0
        touse = self.azimuth[sort]
        rdist = self.dist[sort] * 1.
        while (off<sort.size):
          ref = touse[off]
          ww = np.where(np.abs(touse-ref)<1.e-6)[0]
          if (np.size(ww)>1):
            sort[ww]=sort[ww][np.argsort(rdist[ww])[::-1]]
          off = off+ww.size

      elif (order.upper() == 'EQUAL'):
        sort = np.int32(self.azimuth * 0)
      else:
        print("THE ORDER MUST BE ONE OF THE FOLLOWINGS:")
        print("\tRANDOM")
        print("\tSYMMETRIC")
        print("\tLINEAL")
        print("\tAZIMUTHAL")
        sys.exit()

    #
      self.sort = sort * 1
    #
      return
    #
    def get_times(self, ideal):
    #

      if (ideal!=True):

        to_be_rolled = np.linspace(0.,self.SimulatedTime \
            , self.SimulatedTimes) % self.period

# Dust:
        frac = 0.099999

        cnt = 0
        while (cnt<self.N):

          for i in range(self.mltch):

            offset = ( (self.sort[cnt]//self.mltch) * self.days_between_change ) \
                / self.tstep

            self.time_map[cnt,:,0] = (np.arange(self.Ntimes \
                , dtype=np.float64) * self.tstep + offset * self.tstep) % self.period
            # in days.

            self.time_map[cnt,:,1] = np.min(np.vstack([self.time_map[cnt \
                ,:,0], (np.arange(self.Ntimes,dtype=np.float64) \
                %self.cleandustcadence)*self.tstep]),axis=0)

            cnt = cnt+1

#

        initial_num_cln_since_last_sch = (self.time_map[:,0,0] \
            / self.tstep) / self.cleandustcadence

        inital_amount_dust_in_sim_steps = initial_num_cln_since_last_sch \
            - np.floor(initial_num_cln_since_last_sch)

        first_acc_amount_dust_in_sim_steps \
            = inital_amount_dust_in_sim_steps * frac 

        #
        # Accumulate initial dust for every segment:
        aux = np.floor(initial_num_cln_since_last_sch).astype("i4")
        ginit_dust_val = aux * 0. + first_acc_amount_dust_in_sim_steps \
            * self.tstep

        bu_ginit_dust_val = ginit_dust_val * 1.

        while (np.max(aux)>0):

          ww = np.where(aux > 0)[0]

          ginit_dust_val[ww] = (ginit_dust_val[ww] + self.cleandustcadence \
              * self.tstep) * frac

          aux = aux - 1

        it_sgm = 0
        for it_sgm in range(self.N):

  
          for it_ttt in range(self.Ntimes):
  
            # First time:
            if (it_ttt == 0):
              #print(" Segment %i dust value: %.4f" % (it_sgm, ginit_dust_val[it_sgm]))
              self.time_map[it_sgm, it_ttt, 2] = ginit_dust_val[it_sgm] * 1.
              continue
  
            # If we have changed the segment in this last sim. time:
            aux_seg_chg = self.time_map[it_sgm, it_ttt, 0] / self.tstep
            if (aux_seg_chg < 1):
              if (aux_seg_chg != 0.):
                excess_sim_step = (it_ttt % aux_seg_chg) * self.tstep
              else:
                excess_sim_step = 0.
  
              init_dust_val = excess_sim_step * self.tstep
              self.time_map[it_sgm, it_ttt, 2] = init_dust_val * 1.
  
            else:
  
              # If we have cleaned the segment in this last sim. time:
              if (it_ttt % self.cleandustcadence < 1):
    
                excess_sim_step = (it_ttt % self.cleandustcadence ) \
                    * self.tstep
    
                last_dust_val = self.time_map[it_sgm, it_ttt-1, 2] \
                    + self.tstep * (1. - excess_sim_step)
    
                init_dust_val = excess_sim_step * self.tstep
                self.time_map[it_sgm, it_ttt, 2] = init_dust_val \
                    + last_dust_val * frac
    
              else:
              # Otherwise we continue to accumulate dust:
                # The amount of dust is updated:
                self.time_map[it_sgm, it_ttt, 2] = self.time_map[it_sgm \
                    , it_ttt-1, 2] + self.tstep
    
  
            dum = self.time_map[it_sgm, it_ttt, 2] * 1.
            if (dum!=dum):
              stop()

      #
      return
    #
    #
  #
  if ( (deltat==0.) & (period==0.) ):
    print('deltat and period are 0, one of them must not be 0!')
    return np.nan
  if ( (deltat!=0.) & (period!=0.) ):
    print('deltat and period are non 0, one of them must be 0!')
    return np.nan
  #
  #
  segments = primary_obj(telescope, tstep, tlong, cleandust \
      , period, deltat, multiplechange)
  segments.get_order(order)
  segments.get_times(ideal)
#
  return segments
#

def cinside(\
    ar [double, ndim=1] xir \
    , ar [double, ndim=1] yir \
    , ar [double, ndim=1] xer \
    , ar [double, ndim=1] yer \
    , ar [double, ndim=1] xc \
    , ar [double, ndim=1] yc \
    , int nthreads \
    , int vectors=1 \
    , int complementary=0 \
    , int verbose=0):

  cdef long cnrays = xc.size
  cdef ar[long, ndim=1] res = np.ones((cnrays,), dtype=long)

  cdef int nrefs

  #
  # Internal perimeter:
  # 
  nrefs = xir.size
  if (nrefs==1):
    # Anular mirror:
    # xir: circle radius
    # yir: circle radius (only xr is used)
    inside_circle(nthreads, 0, xir[0], cnrays, <double*>xc.data \
        , <double*>yc.data, <long*>res.data)
  else:
    # Polygon-shaped mirror:
    inside_polygon(nthreads, 0, nrefs, <double*> xir.data \
        , <double*> yir.data, cnrays, <double*>xc.data \
        , <double*>yc.data, <long*>res.data)

  #
  # External perimeter:
  # 
  nrefs = xer.size
  if (nrefs==1):
    # Anular mirror:
    # xer: circle radius
    # yer: circle radius (only xr is used)
    inside_circle(nthreads, 1, xer[0], cnrays, <double*>xc.data \
        , <double*>yc.data, <long*>res.data)
  else:
    # Polygon-shaped mirror:
    inside_polygon(nthreads, 1, nrefs, <double*> xer.data \
        , <double*> yer.data, cnrays, <double*>xc.data \
        , <double*>yc.data, <long*>res.data)

  if (complementary):
    ww = res < 0.5
  else:
    ww = res > 0.5

  if (vectors):
    return xc[ww], yc[ww]
  else:
    return ww



def secondary(telescope):
#
  class mirror_obj(object):

    def __init__(self,telescope):
      vmax = telescope.sec_ext_radius * np.sqrt(telescope.sec_num_esp) * 1.1
      if (telescope.ID == 'GTC'):
        self.ID = 'GTC'
      if (telescope.ID == 'SGTC'):
        self.ID = 'SGTC'
      if (telescope.ID == 'AGTC'):
        self.ID = 'AGTC'
      elif (telescope.ID == 'EELT'):
        self.ID = 'EELT'
      elif (telescope.ID == 'SEELT'):
        self.ID = 'SEELT'
      elif (telescope.ID == 'AEELT'):
        self.ID = 'AEELT'

      return
#
    def get_secondary_geometry(self,inum_esp, isec_ext_radius \
        , isec_int_radius,ishape):

      if (ishape=='hexagonal'):
        #
        #GET AN INITIAL ESTIMATION OF MIRRORS:
        #
        ones = np.ones(inum_esp, dtype=np.float64)
        x_odds = np.arange(inum_esp) - inum_esp / 2
        _, x_odds = np.meshgrid(ones, x_odds)
        y_odds = np.arange(inum_esp) - inum_esp / 2
        y_odds, _ = np.meshgrid(y_odds, ones)
        x_even = (np.arange(inum_esp) - \
            inum_esp / 2) * 2. - 1.
        _, x_even = np.meshgrid(ones, x_even)
        y_even = (np.arange(inum_esp) - \
            inum_esp / 2) * 2. - 1.
        y_even, _ = np.meshgrid(y_even, ones)
        #
        #GIVE SPATIAL DIMENSIONS TO THE POSITIONS:
        #
        ap = isec_ext_radius * np.cos(30./180.*np.pi)
        x_fact_odds = 3. * isec_ext_radius
        y_fact_odds = 2. * ap
        x_fact_even = 1.5 * isec_ext_radius
        y_fact_even = 1.* ap
        x_odds = x_odds.flatten() * x_fact_odds
        y_odds = y_odds.flatten() * y_fact_odds
        x_even = x_even.flatten() * x_fact_even
        y_even = y_even.flatten() * y_fact_even
        #
        rads_odds = np.sqrt(x_odds**2 + y_odds**2)
        rads_even = np.sqrt(x_even**2 + y_even**2)
        #
        #SORT:
        #
        sort = np.argsort(rads_odds)
        x_odds = x_odds[sort]
        y_odds = y_odds[sort]
        rads_odds = rads_odds[sort]
        sort = np.argsort(rads_even)
        x_even = x_even[sort]
        y_even = y_even[sort]
        rads_even = rads_even[sort]
        #
        #CONCATENATE:
        #
        x = np.concatenate([x_odds, x_even])
        y = np.concatenate([y_odds, y_even])
        rads = np.concatenate([rads_odds, rads_even])
        sort =  np.argsort(rads)
        x = x[sort]
        y = y[sort]
        rads = rads[sort]
        #
        ww = np.where(rads >= isec_ext_radius)[0]
        rads = rads[ww]
        x = x[ww]
        y = y[ww]
        #
        #LAST MIRROR
        #
        last = rads[inum_esp-1]
        ww = np.where(rads == last)[0]
        if ( np.max(ww) > inum_esp-1 ):
          sys.exit()
        rads = rads[0:inum_esp]
        x = x[0:inum_esp]
        y = y[0:inum_esp]
#
#
#
#
#
#
#
        sangs = np.arange(6) * np.pi / 3.
        srads = np.sqrt(x**2 + y**2)

        ww = np.where(srads+isec_ext_radius>srads.max())[0]

        sdx = isec_ext_radius * np.cos(sangs)
        sdy = isec_ext_radius * np.sin(sangs)

        sx = x[ww][:,None] + sdx[None,:]
        sy = y[ww][:,None] + sdy[None,:]

        srads = np.sqrt(x[ww]**2 + y[ww]**2)
        srs = np.sqrt(sx**2 + sy**2)

        ww = np.where( (srs-srads[:,None]) > 0 )
        sx = sx[ww[0], ww[1]]
        sy = sy[ww[0], ww[1]]

        sazis = np.arctan2(sy.flatten(), sx.flatten())

        sx = sx.flatten()[np.argsort(sazis)]
        sy = sy.flatten()[np.argsort(sazis)]

        sx = np.concatenate([sx, [sx[0],]])
        sy = np.concatenate([sy, [sy[0],]])

        # Remove duplicated points:
        diffs = np.abs(sx[1:]-sx[0:-1]) + np.abs(sy[1:]-sy[0:-1])
        ww = np.where(diffs > 1.e-5)

        self.esx = sx[ww[0]+1]
        self.esy = sy[ww[0]+1]
#
#
#
        self.isx = isec_int_radius * np.cos(sangs)
        self.isy = isec_int_radius * np.sin(sangs)

      elif (ishape=='anular'):

        self.isx = np.ones((1,)) * isec_int_radius
        self.isy = np.ones((1,)) * isec_int_radius

        self.esx = np.ones((1,)) * isec_ext_radius
        self.esy = np.ones((1,)) * isec_ext_radius

      else:
        print('Secondary:')
        print('\tShape unknown!')
        print('')
        stop()

      return

    def get_mirror(self, telescope):
#
      self.get_secondary_geometry(telescope.sec_num_esp \
          , telescope.sec_ext_radius, telescope.sec_int_radius \
          , telescope.secondary_shape)

      return
    #
#
  secondary = mirror_obj(telescope)
  secondary.get_mirror(telescope)

  return secondary
#
#
#
# ---------------------------------------------------------------------
#

#
def dirt_timing(lamb, time):
#
  def zerodux(lamb, angstrom=True):
  #
    if (angstrom == True):
      lamb_mu = lamb * 1.e-4
    else:
      lamb_mu = lamb * 1.e0
  #
    b=[1.3182408e0,2.44e-2,1.08915181]
    c=[8.79e-3,6.09e-2,1.1e2]
    n=np.sqrt(b[0] * lamb_mu**2 / (lamb_mu**2 - c[0]) + b[1] * lamb_mu**2 / \
        (lamb_mu**2 - c[1]) + b[2] * lamb_mu**2 / (lamb_mu**2 - c[2]) + 1.e0)
  #
    return n
#
  def indexes(elem, lam_in):
#la longitud de onda entra en amstrons
    if (elem.upper() == 'DUST'):
      lam = lam_in * 1.e-4
      n2 = 1.28604141e0 + (1.07044083e0 * lam**2) / \
          (lam**2 - 1.00585997e-2) + (1.10202242e0 * lam**2.) / \
          (lam**2 - 100.e0) 
      n_out = np.sqrt(n2)
      return n_out, np.nan
    if (elem.upper() == 'OX'):
      lam = np.array([2.480000 ,2.066000 ,1.771000 ,1.550000 ,1.378000 ,\
          1.240000 ,1.127000 ,1.033000 ,0.953700 ,0.885600 ,0.826600 ,\
          0.774900 ,0.729300 ,0.688800 ,0.652500 ,0.619900 ,0.590400 ,\
          0.563600 ,0.539100 ,0.516600 ,0.495900 ,0.476900 ,0.459200 ,\
          0.442800 ,0.427500 ,0.413300 ,0.399900 ,0.387500 ,0.375700 ,\
          0.364700 ,0.354200 ,0.344400 ,0.335100 ,0.326300 ,0.317900 ,\
          0.310000 ,0.302400 ,0.295200 ,0.288300 ,0.281800 ,0.275500 ,\
          0.269500 ,0.263800 ,0.258300 ,0.253000 ,0.248000 ,0.243100 ,\
          0.238400 ,0.233900 ,0.229600 ,0.225400 ,0.221400 ,0.217500 ,\
          0.213800 ,0.210100 ,0.206600]) * 1.e4
      n = np.array([1.726000 ,1.736000 ,1.742000 ,1.746000 ,1.749000 ,\
          1.751000 ,1.752600 ,1.754200 ,1.755800 ,1.757400 ,1.759000 ,\
          1.760600 ,1.762200 ,1.763800 ,1.765400 ,1.767000 ,1.768568 ,\
          1.770144 ,1.771736 ,1.773352 ,1.775000 ,1.776656 ,1.778368 ,\
          1.780152 ,1.782024 ,1.784000 ,1.786192 ,1.788496 ,1.790904 ,\
          1.793408 ,1.796000 ,1.798608 ,1.801304 ,1.804096 ,1.806992 ,\
          1.810000 ,1.813576 ,1.817168 ,1.820672 ,1.823984 ,1.827000 ,\
          1.829104 ,1.830832 ,1.832208 ,1.833256 ,1.834000 ,1.834336 ,\
          1.834448 ,1.834392 ,1.834224 ,1.834000 ,1.833776 ,1.833608 ,\
          1.833552 ,1.833664 ,1.834000])
      from scipy.interpolate import interp1d as interpol
      sort = np.argsort(lam)
      lam = lam[sort]
      n = n[sort]
      f = interpol(lam, n)
      #print(lam[0], lam_in, lam[-1])
      n_out = f(lam_in)
      k_out = 0.
      return n_out, k_out

    elif (elem.upper() == 'AL'):
      lam = np.array([200.000000, 177.119995, 153.850006, 137.759995, \
          125.000000, 99.996002, 80.000999, 66.666000, 57.144001, 50.000000, \
          44.444000, 40.000000, 33.333000, 30.996000, 27.552000, 24.797001, \
          22.542999, 20.664000, 19.075001, 17.712000, 15.498000, 13.776000, \
          12.399000, 10.332000, 8.856100, 7.749100, 6.888100, 6.199300, \
          5.635700, 5.166000, 4.768700, 4.428000, 4.132800, 3.874500, \
          3.646600, 3.444000, 3.262800, 3.099600, 2.755200, 2.479700, \
          2.066400, 1.771200, 1.549800, 1.377600, 1.239900, 1.127100, \
          1.033200, 0.999880, 0.968630, 0.939280, 0.911660, 0.885610, \
          0.837740, 0.815690, 0.794780, 0.774910, 0.729320, 0.688810, \
          0.652250, 0.619930, 0.563570, 0.516600, 0.476870, 0.442800, \
          0.413280, 0.364660, 0.326280, 0.309960, 0.247970, 0.206640, \
          0.177120, 0.154980, 0.137760, 0.123990, 0.112710, 0.103320, \
          0.095373, 0.088561, 0.086101, 0.084921, 0.083774, 0.082657, \
          0.082109, 0.081569, 0.081036, 0.080510, 0.079990, 0.079478, \
          0.078472, 0.077491, 0.072932, 0.068881, 0.065255, 0.061993, \
          0.049594, 0.041328, 0.035424, 0.030996, 0.027552, 0.024797, \
          0.022543, 0.020664, 0.019998, 0.019373, 0.018786, 0.018233, \
          0.017712, 0.017463, 0.017220, 0.017208, 0.017196, 0.017172, \
          0.017149, 0.017125, 0.017101, 0.017078, 0.017054, 0.017031, \
          0.017008, 0.016984, 0.016961, 0.016938, 0.016915, 0.016892, \
          0.016869, 0.016846, 0.016823, 0.016800, 0.016777, 0.016755, \
          0.016642, 0.016531, 0.016422, 0.016314, 0.016102, 0.015896, \
          0.015694, 0.015498, 0.015120, 0.014760, 0.014417, 0.014089, \
          0.013776, 0.013477, 0.013190, 0.012915, 0.012652, 0.012399, \
          0.011808, 0.011271, 0.010781, 0.010332, 0.009919, 0.009537, \
          0.009184, 0.008856, 0.008551, 0.008266, 0.007749, 0.007293, \
          0.006888, 0.006526, 0.006199, 0.004959, 0.004133, 0.003100, \
          0.002480, 0.002066, 0.001771, 0.001550, 0.001378, 0.001240, \
          0.001127, 0.001033, 0.000954, 0.000886, 0.000855, 0.000827, \
          0.000816, 0.000805, 0.000800, 0.000799, 0.000798, 0.000796, \
          0.000795, 0.000795, 0.000794, 0.000793, 0.000785, 0.000777, \
          0.000740, 0.000729, 0.000689, 0.000653, 0.000620, 0.000496, \
          0.000413, 0.000354, 0.000310, 0.000248, 0.000207, 0.000177, \
          0.000155, 0.000138, 0.000131, 0.000124]) * 1.e4
      n=np.array([423.959991, 397.929993, 364.040009, 339.619995, \
          318.809998, 274.380005, 233.559998, 202.630005, 177.929993, \
          157.300003, 140.050003, 125.139999, 102.099998, 91.955002, \
          75.748001, 68.535004, 63.554001, 58.580002, 54.412998, 50.951000, \
          43.775002, 38.460999, 33.519001, 26.216000, 20.837000, 16.754999, \
          14.088000, 12.195000, 10.742000, 9.558000, 8.588100, 7.775700, \
          7.079600, 6.480800, 5.956400, 5.490300, 5.073500, 4.709700, \
          3.938000, 3.337200, 2.473800, 1.920500, 1.578200, 1.389900, \
          1.315700, 1.328100, 1.399800, 1.435900, 1.486700, 1.678400, \
          1.973900, 2.280200, 2.694500, 2.766800, 2.767500, 2.615400, \
          2.160600, 1.830100, 1.572400, 1.366000, 1.072800, 0.873400, \
          0.727800, 0.607900, 0.521350, 0.398770, 0.314740, 0.280030, \
          0.181370, 0.126770, 0.094236, 0.072505, 0.057167, 0.046304, \
          0.038468, 0.035753, 0.036437, 0.044168, 0.054863, 0.067041, \
          0.094517, 0.150650, 0.179430, 0.205690, 0.233440, 0.259360, \
          0.282710, 0.303730, 0.340310, 0.371970, 0.491310, 0.572510, \
          0.632420, 0.679120, 0.815120, 0.880130, 0.918020, 0.941890, \
          0.958340, 0.970480, 0.979980, 0.988270, 0.991430, 0.994570, \
          0.997910, 1.001900, 1.007000, 1.010800, 1.016900, 1.017400, \
          1.017900, 1.019200, 1.020600, 1.022600, 1.024900, 1.030500, \
          1.034900, 1.030500, 1.025500, 1.024600, 1.026200, 1.025900, \
          1.021900, 1.019400, 1.018100, 1.017300, 1.016700, 1.016100, \
          1.015600, 1.015100, 1.013200, 1.011800, 1.011000, 1.010600, \
          1.009500, 1.007800, 1.007500, 1.007500, 1.007700, 1.007400, \
          1.006500, 1.006000, 1.005800, 1.004100, 1.001200, 0.996520, \
          0.992650, 0.991230, 0.992850, 0.994150, 0.992330, 0.991390, \
          0.989410, 0.987610, 0.987930, 0.988830, 0.989340, 0.989660, \
          0.989120, 0.989090, 0.990070, 0.990540, 0.991110, 0.993130, \
          0.994800, 0.996940, 0.997970, 0.998600, 0.998980, 0.999240, \
          0.999400, 0.999530, 0.999610, 0.999680, 0.999740, 0.999790, \
          0.999810, 0.999830, 0.999850, 0.999870, 0.999890, 0.999910, \
          0.999910, 0.999920, 0.999920, 0.999930, 0.999920, 0.999900, \
          0.999870, 0.999860, 0.999840, 0.999840, 0.999850, 0.999860, \
          0.999870, 0.999910, 0.999940, 0.999950, 0.999970, 0.999980, \
          0.999984, 0.999989, 0.999991, 0.999993, 0.999994, 0.999995])#*2.#*0.
      #n=1.38888498e-4*lam-1.62620406e-1
      k = np.array([4.837000e+02 ,4.585000e+02 ,4.296200e+02 ,4.089200e+02 ,\
          3.917100e+02 ,3.543500e+02 ,3.210800e+02 ,2.954200e+02 ,\
          2.753400e+02 ,2.582600e+02 ,2.434300e+02 ,2.301900e+02 ,\
          2.081000e+02 ,1.999900e+02 ,1.817800e+02 ,1.648100e+02 ,\
          1.534500e+02 ,1.442300e+02 ,1.360900e+02 ,1.294900e+02 ,\
          1.183900e+02 ,1.089600e+02 ,1.012800e+02 ,8.819700e+01 ,\
          7.827400e+01 ,6.985700e+01 ,6.284100e+01 ,5.715600e+01 ,\
          5.251800e+01 ,4.859300e+01 ,4.525700e+01 ,4.236700e+01 ,\
          3.982600e+01 ,3.759500e+01 ,3.560800e+01 ,3.381400e+01 ,\
          3.218300e+01 ,3.073700e+01 ,2.758000e+01 ,2.500400e+01 ,\
          2.098200e+01 ,1.799100e+01 ,1.565600e+01 ,1.378400e+01 ,\
          1.224500e+01 ,1.096900e+01 ,9.891400e+00 ,9.493900e+00 ,\
          9.065500e+00 ,8.597000e+00 ,8.305800e+00 ,8.113400e+00 ,\
          8.187800e+00 ,8.257300e+00 ,8.386600e+00 ,8.491400e+00 ,\
          8.356500e+00 ,8.060100e+00 ,7.735400e+00 ,7.405200e+00 ,\
          6.783900e+00 ,6.241800e+00 ,5.778100e+00 ,5.367600e+00 ,\
          5.000800e+00 ,4.395700e+00 ,3.916500e+00 ,3.708100e+00 ,\
          2.902900e+00 ,2.356300e+00 ,1.951900e+00 ,1.636600e+00 ,\
          1.377500e+00 ,1.155500e+00 ,9.567700e-01 ,7.716300e-01 ,\
          5.908600e-01 ,3.911500e-01 ,2.929300e-01 ,2.342000e-01 ,\
          1.658900e-01 ,1.104100e-01 ,9.422300e-02 ,7.995900e-02 ,\
          6.834800e-02 ,6.140700e-02 ,5.669700e-02 ,5.334900e-02 ,\
          4.832000e-02 ,4.420200e-02 ,3.240900e-02 ,2.768100e-02 ,\
          2.477000e-02 ,2.234000e-02 ,1.589400e-02 ,1.165100e-02 ,\
          9.312100e-03 ,7.846600e-03 ,6.619100e-03 ,5.746900e-03 ,\
          5.000400e-03 ,4.369600e-03 ,4.239700e-03 ,4.109200e-03 ,\
          3.892600e-03 ,3.673000e-03 ,3.542500e-03 ,3.495700e-03 ,\
          3.524900e-03 ,3.510800e-03 ,3.487700e-03 ,3.476700e-03 ,\
          3.440000e-03 ,3.414100e-03 ,3.621800e-03 ,4.116400e-03 ,\
          1.247600e-02 ,2.007200e-02 ,2.001200e-02 ,1.914500e-02 ,\
          1.956400e-02 ,2.422700e-02 ,2.601800e-02 ,2.543200e-02 ,\
          2.483100e-02 ,2.449900e-02 ,2.437500e-02 ,2.428200e-02 ,\
          2.422800e-02 ,2.418400e-02 ,2.434300e-02 ,2.402000e-02 ,\
          2.395500e-02 ,2.385300e-02 ,2.518000e-02 ,2.475700e-02 ,\
          2.450100e-02 ,2.447600e-02 ,2.546000e-02 ,2.682600e-02 ,\
          2.823200e-02 ,2.895600e-02 ,3.091800e-02 ,3.339200e-02 ,\
          3.531100e-02 ,3.588300e-02 ,3.306100e-02 ,2.992000e-02 ,\
          2.441500e-02 ,2.545200e-02 ,2.492800e-02 ,2.406300e-02 ,\
          2.342100e-02 ,2.060600e-02 ,1.776500e-02 ,1.630400e-02 ,\
          1.543700e-02 ,1.477300e-02 ,1.372800e-02 ,1.098700e-02 ,\
          9.651700e-03 ,8.471600e-03 ,7.509900e-03 ,4.186300e-03 ,\
          2.349200e-03 ,9.688700e-04 ,4.350300e-04 ,2.200100e-04 ,\
          1.181100e-04 ,7.299200e-05 ,4.677100e-05 ,3.111500e-05 ,\
          2.158600e-05 ,1.562800e-05 ,1.163600e-05 ,8.940000e-06 ,\
          7.842200e-06 ,6.804900e-06 ,6.481800e-06 ,6.190100e-06 ,\
          6.048400e-06 ,1.102600e-05 ,2.368700e-05 ,3.291100e-05 ,\
          4.113100e-05 ,4.999800e-05 ,7.184700e-05 ,8.541800e-05 ,\
          7.848000e-05 ,7.486200e-05 ,6.163300e-05 ,5.797400e-05 ,\
          4.713800e-05 ,3.849700e-05 ,3.207900e-05 ,1.364200e-05 ,\
          7.140800e-06 ,4.048400e-06 ,2.484300e-06 ,1.072600e-06 ,\
          5.510200e-07 ,3.146300e-07 ,1.917300e-07 ,1.248800e-07 ,\
          1.272000e-07 ,8.241000e-08])#*2.#*0.
      #k=0.00099066*lam+0.4526242
      from scipy.interpolate import interp1d as interpol
      sort = np.argsort(lam)
      lam = lam[sort]
      n = n[sort]
      k = k[sort]
      f = interpol(lam, n)
      n_out = f(lam_in)
      f = interpol(lam, k)
      k_out = f(lam_in)

      return n_out,k_out
#
  def thick(th_max, tau, t):
#
    th = th_max * (1.e0 - np.exp( -t / tau))
    return th
#
  def Dthick(alpha, t):
#
    th = alpha * t
    
    return th
#
  class materials(object):
#
    def __init__(self, n_s, n_in, n_d, n_c, n_ox, d_ox, d_c, d_s, d_d, lamb):
      self.n_s = n_s + 0j
      self.n_in = n_in + 0j
      self.n_d = n_d + 0j
      self.n_c = n_c * 1.
      self.n_ox = n_ox + 0j

      self.d_ox = d_ox * 1.
      self.d_c = d_c * 1.
      self.d_s = d_s * 1.
      self.d_d = d_d * 1.

      self.lamb = lamb * 1.

      return
#
  #FIRST APPROACH: ONE DUST PARTICLE THICK PER DAY
  # AlPhA
# Dust:
  alpha = 0.003

  n_s = zerodux(lamb * 1.e-4)
  n_in = 1.e0
  n_d = 1.52e0
  elem_esp = 'al'
  n_g, k_g = indexes(elem_esp,lamb)
  n_c = n_g + k_g * 1j
  elem_ox = 'ox'
  n_ox, k_ox = indexes(elem_ox,lamb)
  n_ox = n_ox + k_ox * 1j
  
  tau_o = 1392.1484e0
  tau_d = 1.8e2
  th_m_o = 1.2e3
  th_m_d = 1.2e4
  d_mir = 1.2e3   #angstroms

  d_ox = 1. * thick( th_m_o , tau_o , time[:,:,0] ) 
  d_c = (d_mir - d_ox)  #el espesor del conductor sera el original (d_mir) menos lo que se ha oxidado.
  d_d = 1. * Dthick ( alpha, time[:,:,2] )

  d_s = d_d.real * 0. + 1.e4

  material = materials(n_s, n_in, n_d, n_c, n_ox, d_ox, d_c, d_s, d_d, lamb)

  return material
#
#
#
# ---------------------------------------------------------------------
#
def get_intersection_focus(dir_cos, x, z, f):

  a = dir_cos[2,:]/dir_cos[0,:]

  b = z - a * x

  xf = (f - b) / a

  return xf

def get_intersection(dir_cos, x, z, r, k, offset):

  a = dir_cos[2,:]/dir_cos[0,:]

  b = z - a * x

  cte1 = -2. * r
  cte2 = 1. + k

  ctea = -cte1 * offset + cte2 * offset**2
  cteb = cte1 - 2. * cte2 * offset

  cteq = 1. + cte2 * a**2
  ctel = 2. * cte2 * a * b + cteb * a
  ctei = cte2 * b**2 + cteb * b + ctea

  xs2 = np.sqrt(ctel**2-4.*cteq*ctei)
  xs1 = (-ctel + xs2) / (2. * cteq)
  xs2 = (-ctel - xs2) / (2. * cteq)

  dir_cos_n1, zs1 = new2_get_gcurv(xs1,k=k,r=r, ind=offset,norm=True)
  dir_cos_n2, zs2 = new2_get_gcurv(xs2,k=k,r=r, ind=offset, norm=True)

  x = xs1 * 1.
  z = zs1 * 1.
  dir_cos = dir_cos_n1 * 1.

  adiff1 = np.abs(zs1-offset)
  adiff2 = np.abs(zs2-offset)

  ww = np.where(adiff2<adiff1)

  x[ww] = xs2[ww] * 1.
  z[ww] = zs2[ww] * 1.
  dir_cos[:,ww] = dir_cos_n2[:,ww] * 1.

  return x,z,dir_cos

def rotate_2d(vec, ang):

  cs = np.cos(ang)
  sn = np.sin(ang)

  res = vec * 0.
  res[0,:]=vec[0,:]*cs+vec[2,:]*sn
  res[2,:]=vec[2,:]*cs-vec[0,:]*sn

  return res


def get_direction_cosines(x,y,z):

  dd = np.sqrt(x**2 + y**2 + z**2)
  l = x / dd
  m = y / dd
  n = z / dd

  return np.array([l,m,n])

def new2_get_gcurv(x,k=0.,r=0.,f=0.,e=0.,norm=False,ind=0.):
  #
  # Following: https://en.wikipedia.org/wiki/Conic_constant
  # n_ are the derivatives of the equation with respect to...
  # ... each variable.

  def get_conic_sols(ia,ib,ic):
    isol1=(-ib+np.sqrt(ib**2-4.*ia*ic))/(2.*ia)
    isol2=(-ib-np.sqrt(ib**2-4.*ia*ic))/(2.*ia)
    return isol1, isol2

  if ( ( (k==0) & (e==0) ) | ( (r==0) & (f==0) ) ):
    print('new2_get_gcurv requires either k/e and r/f')
    sys.exit()

  if (k==0):
    k = -e**2
  if (r==0):
    r = 2. * f

  b=-2.*r
  a=k+1.

  #
  # If x=0,y=0 I want z to be at ind position:
  #cte = b * ind + a * ind**2

  c=x**2#+y**2

  n_x = 2.e0 * x
  n_y = 0.e0 * x # We rotate the system so that only x and z are relevant

  if (k==-1):

    sol = -c/b

  else:

    sol1, sol2 = get_conic_sols(a,b,c)

    t1, t2 = get_conic_sols(a,b,0)#-cte)
    #
    # Assuming a mirror for a telescope is given by the solution closest to the...
    # ... optical axis, we can easily solve the ambiguity:
    if (abs(t1)<abs(t2)):
      sol = sol1 * 1.
    else:
      sol = sol2 * 1.
    nsol = b + 2 * a * sol
    sol = sol + ind

  if (norm==True):
    dir_cos = get_direction_cosines(n_x,n_y,nsol)

    return dir_cos, sol

  return np.array([n_x, n_y, nsol]), sol

#
#
# ---------------------------------------------------------------------
#
def get_primary_rays(layout, telescope):

  if ( (layout!='polar') & (layout!='cartesian') ):
     print('')
     print('\t Unknown layout: %s for ray tracing!' % (layout,))
     print('\t Allowed options are: "polar", "cartesian"')
     print('')
     exit(1)

  factor=1.001
  # Polar:
  if (layout=='polar'):
    # set rays to be considered in a polar layout:
    azi_val = 24 #24
    rad_to_azi_rat = 24
    irads = np.linspace(0., telescope.dim_x * factor, np.max([3,telescope.num//rad_to_azi_rat]), dtype=np.float64)#[1:]
    irads = np.linspace(0., telescope.dim_x * factor, (telescope.num//2)*2+1, dtype=np.float64)
    drads = irads[1:]-irads[0:-1]
    irads = irads[1:] - drads / 2.

    val = 0.#0.#14.#.0

    iangs = (np.linspace(0.+val, 360.+val, np.max([(rad_to_azi_rat*telescope.num)//azi_val,1])*azi_val+1, dtype=np.float64)) / 180. * np.pi#[0:-1]
    iangs = np.linspace(0.+val, 360.+val, (telescope.num//2)*2+1, dtype=np.float64) / 180. * np.pi
    dangs = iangs[1:]-iangs[0:-1]
    iangs = iangs[1:] - dangs / 2.

    #print(np.min(iangs))

    # Transform to x and y coordinates:
    
    ixs = irads[:,None] * np.cos(iangs)[None,:]
    iys = irads[:,None] * np.sin(iangs)[None,:]

    d1s = drads.mean()
    d2s = dangs.mean()

    x1d = ixs.flatten()
    y1d = iys.flatten()

    rad = np.sqrt(x1d**2 + y1d**2)
    area = rad * d1s * d2s

#
#
#
#
#
#
#

    ang_ref = (2. * np.pi) / telescope.num
    angs = ang_ref * (irads.min() / irads)
    nds = np.round((2. * np.pi) / angs)

    nds = np.round((nds / nds[nds.size//3]) * telescope.num)

    ww = np.where(nds<telescope.num)
    nds[ww] = telescope.num
    nds = (nds//2) * 2
    nds = (nds//azi_val) * azi_val
    angs = (2. * np.pi) / nds

    totn = np.int32(nds.sum())
    #print(nds.min(), telescope.num**2, totn)

    x1d = np.zeros((totn,), dtype="f8")
    y1d = np.zeros((totn,), dtype="f8")
    rad = np.zeros((totn,), dtype="f8")
    area = np.zeros((totn,), dtype="f8")


    offset = 0
    for itn in range(irads.size):

      itnum = np.int32(nds[itn])

      iangs = np.linspace(0.+val, 360.+val, itnum+1, dtype=np.float64) / 180. * np.pi
      dangs = iangs[1:] - iangs[0:-1]
      iangs = iangs[1:] - dangs / 2.

      x1d[offset:offset+itnum] = irads[itn] * np.cos(iangs)
      y1d[offset:offset+itnum] = irads[itn] * np.sin(iangs)

      rad[offset:offset+itnum] = np.sqrt(x1d[offset:offset+itnum]**2 + y1d[offset:offset+itnum]**2)
      area[offset:offset+itnum] = rad[offset:offset+itnum] * drads[itn] * dangs

      print(" %i  ;  %.4f  ;  %.8e" % (itnum, rad[offset], area[offset], ), end='\r', flush=True)

      offset += itnum

    print("")
#
#
#
#
#
#
#


  # Cartesian:
  elif (layout=='cartesian'):
    ixs = np.linspace(-telescope.dim_x*factor, telescope.dim_x*factor, (telescope.num//2)*2+1, dtype=np.float64)#[1:]
    iys = np.linspace(-telescope.dim_y*factor, telescope.dim_y*factor, (telescope.num//2)*2+1, dtype=np.float64)#[1:]
    ixs = np.linspace(-telescope.dim_x*factor, telescope.dim_x*factor, telescope.num, dtype=np.float64)#[1:]
    iys = np.linspace(-telescope.dim_y*factor, telescope.dim_y*factor, telescope.num, dtype=np.float64)#[1:]
  
    d1s = np.mean(ixs[1:]-ixs[0:-1])
    d2s = np.mean(iys[1:]-iys[0:-1])

    ix1s = ixs[1:]-d1s/2
    iy1s = iys[1:]-d2s/2

    print(ixs.shape)

    ixs = ix1s[:,None] * np.ones(iy1s.size, dtype=np.float64)[None,:]
    iys = iy1s[None,:] * np.ones(ix1s.size, dtype=np.float64)[:,None]

    print(ixs.shape)

    x1d = ixs.flatten()
    y1d = iys.flatten()

    rad = np.sqrt(x1d**2 + y1d**2)
    area = x1d * 0. + d1s * d2s

  else:
    print('Get geometry:')
    print('\tUnknown layout.')
    print('')
    stop()
#
  return x1d, y1d, rad, area

def get_geometry(telescope, beam, secondary, segments \
    , nthreads=1, osecondary=False, layout='polar', check_rays=False):
#
  from pdb import set_trace as stop
  x1d, y1d, rad1d, area1d = get_primary_rays(layout, telescope)
#
  if (osecondary==True):
    print(" !! PLOT !!")
    pl.figure(1)
    pl.clf()
    fg,ax=pl.subplots(ncols=1, nrows=1, num=1)
    ax.contour(secondary.mirror.T, [0.5], extent=[np.min(secondary.x) \
        , np.max(secondary.x), np.max(secondary.y), np.min(secondary.y)],colors='r')

  if (check_rays):

    print(" !! PLOT !!")
    pl.figure(1)
    pl.clf()
    fg,ax=pl.subplots(ncols=1, nrows=1, num=1)
    ax.scatter(x1d,y1d,color='k', s=0.1, edgecolor=None)

# First of all, we remove all the rays that come from behind the secondary:

  ww = cinside(secondary.isx, secondary.isy, secondary.esx, secondary.esy, x1d, y1d, nthreads, complementary=1, verbose=1, vectors=0)

  if (len(ww)>0):
    x1d = x1d[ww] * 1.
    y1d = y1d[ww] * 1.
    rad1d = rad1d[ww] * 1.
    area1d = area1d[ww] * 1.

#

  if (check_rays):
    print(" !! PLOT !!")
    ax.scatter(x1d,y1d,color='r', s=0.1, edgecolor=None)
    pl.figure(2)
    pl.clf()
    fg,ax=pl.subplots(ncols=1, nrows=1, num=2)
    ax.scatter(x1d,y1d,color='k', s=0.1, edgecolor=None)
    stop()

#
  def check_primary_id(iteles, ixpos, iypos, igxs, igys, igrad, igarea):
#
    if (iteles.primary_shape == 'hexagonal'):
      dumgx = igxs - ixpos
      dumgy = igys - iypos
      resgx = igxs - ixpos
      resgy = igys - iypos
#
      a = 30.e0 / 1.8e2 * np.pi
      costhU = np.cos(a) * iteles.radius
# Select pixels that fulfil first criterion:
      ww = np.where(np.abs(dumgy) <= costhU)[0]
      dumgx = dumgx[ww] * 1.
      dumgy = dumgy[ww] * 1.
      resgx = resgx[ww] * 1.
      resgy = resgy[ww] * 1.
      resrad = igrad[ww] * 1.
      resarea = igarea[ww] * 1.
# Rotate system 60 deg:
      a2 = 60.e0 / 1.8e2 * np.pi
      ca = np.cos(a2)
      sa = np.sin(a2)
      rot = np.array([ca, -sa, sa, ca]).reshape(2,2)
#
      out = np.dot(rot, np.vstack([dumgx,dumgy]))
      dumxg = out[0,:] * 1.
      dumyg = out[1,:] * 1.
# Select pixels that fulfil second criterion:
      ww = np.where( np.abs(dumyg) <= costhU)[0]
      dumgx = dumgx[ww] * 1.
      dumgy = dumgy[ww] * 1.
      resgx = resgx[ww] * 1.
      resgy = resgy[ww] * 1.
      resrad = resrad[ww] * 1.
      resarea = resarea[ww] * 1.
#
      rot = np.array([ca, sa, -sa, ca]).reshape(2,2)
#
      out = np.dot(rot, np.vstack([dumgx,dumgy]))
      dumxg = out[0,:] * 1.
      dumyg = out[1,:] * 1.
#
      ww = np.where( np.abs(dumyg) <= costhU)[0]
      dumgx = dumgx[ww] * 1.
      dumgy = dumgy[ww] * 1.
      resgx = resgx[ww] * 1.
      resgy = resgy[ww] * 1.
      resrad = resrad[ww] * 1.
      resarea = resarea[ww] * 1.
#
# Anular:
    elif (iteles.primary_shape == 'anular'):
      dumgx = igxs - ixpos
      dumgy = igys - iypos
      resgx = igxs - ixpos
      resgy = igys - iypos

      dumrad = np.sqrt(dumgx**2+dumgy**2)
      ww = np.where(dumrad <= iteles.rad_max)[0]

      dumrad = dumrad[ww]
      resgx = resgx[ww]
      resgy = resgy[ww]
      resrad = igrad[ww] * 1.
      resarea = igarea[ww] * 1.

      ww = np.where(dumrad >= iteles.rad_min)[0]

      dumrad = dumrad[ww]
      resgx = resgx[ww]
      resgy = resgy[ww]
      resrad = resrad[ww] * 1.
      resarea = resarea[ww] * 1.

# Otherwise:
    else:
      print('Primary shape unknown!')
      print('')
      print('')
      exit()
#
    return resgx + ixpos, resgy + iypos, resrad, resarea

  def check_primary(teles, six, siy, x, y, nthreads):
    if (teles.primary_shape == 'hexagonal'):

      hangs = np.arange(6) * np.pi / 3. + np.pi/6.
      sxprim = telescope.radius * np.sin(hangs)
      syprim = telescope.radius * np.cos(hangs)

      rx = np.ones((1,)) * 0.
      ry = np.ones((1,)) * 0.

      ix = x - six
      iy = y - siy

      wpos = cinside(rx, ry, sxprim, syprim, ix, iy, nthreads, vectors=0)
    elif (teles.primary_shape == 'anular'):

      rix = np.ones((1,)) * teles.rad_min
      riy = np.ones((1,)) * teles.rad_min

      rex = np.ones((1,)) * teles.rad_max
      rey = np.ones((1,)) * teles.rad_max

      ix = x - six
      iy = y - siy

      wpos = cinside(rix, riy, rex, rey, ix, iy, nthreads, vectors=0)

    return wpos


  ntot=0


  if (check_rays):
    print(" !! PLOT !!")
    pl.figure(1)
    pl.clf()
    pl.figure(2)
    pl.clf()
    fg1, ax1 = pl.subplots(num=1,nrows=1,ncols=1)
    fg2, ax2 = pl.subplots(num=2,nrows=1,ncols=1)

    pl.figure(5)
    pl.clf()
    fg,ax=pl.subplots(ncols=1, nrows=1, num=5)
    ax.scatter(x1d,y1d,color='k', s=0.1, edgecolor=None)
    pl.draw()
    print('?')

  for it_sgm in range(segments.N):

    tm1 = tm.time()

    wpos = check_primary(telescope, segments.xpos[it_sgm], segments.ypos[it_sgm], x1d, y1d, nthreads)

    ww = np.where(wpos>0.5)[0]
    if (len(ww)>0):

      wwn = np.where(wpos<0.5)[0]


      x = x1d[ww] * 1.
      y = y1d[ww] * 1.
      rad = rad1d[ww] * 1.
      area = area1d[ww] * 1.
    #  print("")
    #  print(x1d.size, ww.size, wwn.size)
    #  print("")
      x1d = x1d[wwn] * 1.
      y1d = y1d[wwn] * 1.
      rad1d = rad1d[wwn] * 1.
      area1d = area1d[wwn] * 1.













    if (check_rays):
      print(" !! PLOT !!")
      ax1.scatter(x, y, s=5, edgecolor='none')

    #
    # First, we have to rotate the system so that each point is in the plane...
    # ...containing the normal and incident ray

    th = np.arctan2(y, x)

    #
    # Rotate system:
    cs = np.cos(th)
    sn = np.sin(th)
    x1 = x * cs + y * sn
    y1 = x * (-sn) + y * cs

    # Director cosines for input beam:
    dir_cos_i = get_direction_cosines(beam.l_inc, beam.m_inc \
        , beam.n_inc)[:,None]*np.ones(x1.size, dtype=np.float64)[None,:]
    dir_cos_n1, z_1 = new2_get_gcurv(x1,k=telescope.k_1,r=telescope.rc_1 \
        , norm=True)

    i_1 = np.arccos(np.sum(dir_cos_n1*dir_cos_i,axis=0))

    dir_cos_o1 = rotate_2d(dir_cos_i, -2 * i_1 * np.sign(x1))

    x2, z_2, dir_cos_n2 = get_intersection(dir_cos_o1, x1, z_1 \
        , telescope.nrc_2, telescope.k_2, telescope.dmm)

    i_2 = np.arccos(np.sum(dir_cos_n2*dir_cos_o1,axis=0))

# Not necessary here      dir_cos_f = rotate_2d(dir_cos_o1, 2 * i_2 * np.sign(x2))
# Not necessary here
# Not necessary here      zf=-10.
# Not necessary here      xf = get_intersection_focus(dir_cos_f, x2, z_2,zf)

    #
    # Back-rotate system:
    cs = np.cos(-th)
    sn = np.sin(-th)
    y2 = x2 * 0.
    x_sec = x2 * cs + y2 * sn
    y_sec = x2 * (-sn) + y2 * cs

    tm2 = tm.time()
    isin = cinside(secondary.isx, secondary.isy, secondary.esx, secondary.esy, x_sec, y_sec, nthreads, vectors=0)
    tm3 = tm.time()

    ww = np.where(isin > 0.5)[0]

    if (check_rays):
      print(" !! PLOT !!")
      ax1.scatter(x[ww], y[ww], s=5, edgecolor='none')
      print(" !! PLOT !!")
      ax2.scatter(x_sec, y_sec, s=1, edgecolor='none')


    segments.i['%i' % (it_sgm,)]=i_1[ww]
    segments.i2['%i' % (it_sgm,)]=i_2[ww]
    segments.th1['%i' % (it_sgm,)]=th[ww]

    segments.rad['%i' % (it_sgm,)] = rad[ww] * 1.
    segments.area['%i' % (it_sgm,)] = area[ww] * 1.

    if (osecondary==True):
      print(" !! PLOT !!")
      ax.plot(x_sec[ww], y_sec[ww],marker='.' \
          , markersize=1, linestyle='none', alpha=0.5)

    ntot = ntot + ww.size
    tm4 = tm.time()


    txt = " Segment id[#]=%i ; nrays= %i" % (it_sgm, isin.sum())
    txt = txt + " -> it done! %.5f (%.5f sec.)" % ((tm3-tm2) / (tm4-tm1) * 100.,(tm4-tm1),)
    print(txt, end='\r', flush=True)

  #
  if (osecondary==True):
    print(" !! PLOT !!")
    ax.set_xlim(np.min(secondary.x), np.max(secondary.x))
    ax.set_ylim(np.min(secondary.y), np.max(secondary.y))
    ax.set_aspect('equal')
    stop()

  print("")
  print('Total number of rays considered: %i' % ntot)

  if (check_rays):
    print(" !! PLOT !!")
    ax1.scatter(x1d, y1d, s=5, edgecolor='none', color='k')
    fg1.show()
    fg2.show()
    stop()

  return
#
#
#
# ---------------------------------------------------------------------
#

def chex_syst(\
    ar [double, ndim=1] i \
    , ar [double, ndim=1] th_1 \
    , ar [double, ndim=1] i_22 \
    , double lamb \
    , complex n_s \
    , complex n_c \
    , complex n_ox \
    , complex n_d \
    , complex n_in \
    , double d_s \
    , double d_c \
    , double d_ox \
    , double d_d \
    , int nthreads):
#


  cdef long cnrays = i.size
  cdef int cnlayers = 4  # Dust is not considered with thin-film approach

  cdef double clamb = lamb

  cdef ar[double, ndim=1] cnarr = np.zeros((cnlayers), dtype="f8")
  cnarr[3] = n_s.real
  cnarr[2] = n_c.real
  cnarr[1] = n_ox.real
  cnarr[0] = n_in.real

  cdef ar[double, ndim=1] ckarr = np.zeros((cnlayers), dtype="f8")
  ckarr[3] = n_s.imag
  ckarr[2] = n_c.imag
  ckarr[1] = n_ox.imag
  ckarr[0] = n_in.imag

  cdef ar[double, ndim=1] cdarr = np.zeros((cnlayers), dtype="f8")
  cdarr[3] = d_s
  cdarr[2] = d_c
  cdarr[1] = d_ox
  cdarr[0] = -1

  cdef ar[double, ndim=3] result = np.zeros((cnrays,4,4), dtype="f8")
#
  eval_segment(cnrays, cnlayers, <double*>i.data, <double*>th_1.data \
      , <double*>i_22.data, lamb, <double*>cnarr.data, <double*>ckarr.data \
      , <double*>cdarr.data, nthreads, <double*>result.data)
#
#
#  return result
#
#

  return result
#
#
# ---------------------------------------------------------------------
#
#
#
# ---------------------------------------------------------------------
#

#
# ---------------------------------------------------------------------
#

def get_mdust(p0, mean_dust):
#
#
  mat = np.zeros((4,4))
  mat[0,0] = 1. - mean_dust.real
  mat[1,1] = 1. - mean_dust.real
  mat[2,2] = 1. - mean_dust.real
  mat[3,3] = 1. - mean_dust.real
#
  return mat
#
#
# ---------------------------------------------------------------------
#

def get_mueller_time(segments, materials,cdust=False, thetaref=0., nthreads=1):

#
  tarea = 0.
  for it in range(segments.N):
    tarea += np.sum(segments.area['%i' % (it,)])

#RUN TIMES:
  res = np.zeros((segments.SimulatedTimes, segments.N, 4, 4))
  avg_res = np.zeros((segments.SimulatedTimes, 4, 4))
  for it_tms in range(segments.SimulatedTimes):
#RUN SEGMENTS:
    npts_seg = np.ones((segments.N,), dtype=np.float64)
    ntot2=0

    for it_sgm in range(segments.N):
#
      if (segments.i['%i' % (it_sgm,)].size<1):
        npts_seg[it_sgm] = 0.
        res[it_tms,it_sgm,:,:] = res[it_tms,it_sgm,:,:] * 0.
        continue

      mmpa_full = chex_syst(segments.i['%i' % (it_sgm,)], \
          segments.th1['%i' % (it_sgm,)], \
          segments.i2['%i' % (it_sgm,)], \
          materials.lamb, materials.n_s, materials.n_c\
          , materials.n_ox, materials.n_d, materials.n_in, \
          materials.d_s[it_sgm, it_tms], \
          materials.d_c[it_sgm, it_tms]\
          , materials.d_ox[it_sgm, it_tms], materials.d_d[it_sgm, it_tms] \
          , nthreads)

      #GET Mdust
      if (cdust==True):
        mmpd = get_mdust(0.1, materials.d_d[it_sgm, it_tms])
        res_full = mmpa_full.dot(mmpd)
      else:
        res_full = mmpa_full * 1.

      nn,_,_ = res_full.shape
      ntot2=ntot2+nn

      res[it_tms,it_sgm,:,:]=np.sum(res_full[:,:,:] \
          * segments.area['%i' % (it_sgm,)][:,None,None], axis=0) / tarea
      #

    txt = " t=%i(/%i) ; " % (it_tms, segments.SimulatedTimes,)
    txt += " nrays=%i ; " % (ntot2, )
    txt += " area=%.8f" % (tarea, )
    print(txt, end='\r', flush=True)
    

    # Since we consider all the points of each segment, the average must...
    # ... be done accordingly:
    avg_res[it_tms,:,:] = np.sum(res[it_tms,:,:,:], axis=0)

  return avg_res, res, npts_seg


#
# ---------------------------------------------------------------------
#

###><TA::
###><TA::def plot_primary(telescope, beam, secondary, segments \
###><TA::    , pax=None, prays=False, fcol='', highlight=None\
###><TA::    , layout='polar', scol=''):
###><TA::#
###><TA::  import PyColors as PyC
###><TA::
###><TA::  def plot_primary_segment(iteles, ixpos, iypos, i2ax, ifcol, segment_label=-np.inf):
###><TA::#
###><TA::    if (iteles.primary_shape == 'hexagonal'):
###><TA::
###><TA::      dy = np.sqrt(iteles.radius**2-(iteles.radius/2.)**2)
###><TA::
###><TA::      p1x = ixpos + iteles.radius
###><TA::      p1y = iypos + 0.
###><TA::
###><TA::      p2x = ixpos + iteles.radius / 2.
###><TA::      p2y = iypos + dy
###><TA::
###><TA::      p3x = ixpos - iteles.radius / 2.
###><TA::      p3y = iypos + dy
###><TA::
###><TA::      p4x = ixpos - iteles.radius
###><TA::      p4y = iypos + 0.
###><TA::
###><TA::      p5x = ixpos - iteles.radius / 2.
###><TA::      p5y = iypos - dy
###><TA::
###><TA::      p6x = ixpos + iteles.radius / 2.
###><TA::      p6y = iypos - dy
###><TA::
###><TA::      ix = [ p1x, p2x, p3x, p4x, p5x, p6x, p1x ]
###><TA::      iy = [ p1y, p2y, p3y, p4y, p5y, p6y, p1y ]
###><TA::
###><TA::      if (len(ifcol)==0):
###><TA::
###><TA::        im = i2ax.plot(ix,iy,color='k',lw=1.)
###><TA::
###><TA::      else:
###><TA::
###><TA::        #from pdb import set_trace as stop
###><TA::        #stop()
###><TA::
###><TA::        ix1 = [p4x, p3x, p2x, p1x]
###><TA::        iy1 = [p4y, p3y, p2y, p1y]
###><TA::        iy2 = [p4y, p5y, p6y, p1y]
###><TA::
###><TA::        im = i2ax.fill_between(ix1, iy1, iy2, color=ifcol[0])
###><TA::        im = i2ax.plot(ix,iy,color='k',lw=1.)
###><TA::
###><TA::  #
###><TA::#WRONG!::# Anular:
###><TA::#WRONG!::    elif (iteles.primary_shape == 'anular'):
###><TA::#WRONG!::      dumgx = igxs - ixpos
###><TA::#WRONG!::      dumgy = igys - iypos
###><TA::#WRONG!::      resgx = igxs - ixpos
###><TA::#WRONG!::      resgy = igys - iypos
###><TA::#WRONG!::
###><TA::#WRONG!::      dumrad = np.sqrt(dumgx**2+dumgy**2)
###><TA::#WRONG!::      ww = np.where(dumrad <= iteles.rad_max)[0]
###><TA::#WRONG!::
###><TA::#WRONG!::      dumrad = dumrad[ww]
###><TA::#WRONG!::      resgx = resgx[ww]
###><TA::#WRONG!::      resgy = resgy[ww]
###><TA::#WRONG!::
###><TA::#WRONG!::      ww = np.where(dumrad >= iteles.rad_min)[0]
###><TA::#WRONG!::
###><TA::#WRONG!::      dumrad = dumrad[ww]
###><TA::#WRONG!::      resgx = resgx[ww]
###><TA::#WRONG!::      resgy = resgy[ww]
###><TA::#WRONG!::
###><TA::#WRONG!::      im = []
###><TA::
###><TA::# Otherwise:
###><TA::    else:
###><TA::      print('Primary shape unknown!')
###><TA::      print('')
###><TA::      print('')
###><TA::      sys.exit()
###><TA::#
###><TA::    if (segment_label!=-np.inf):
###><TA::      font = {'size':10}
###><TA::      i2ax.text(ixpos, iypos, r"%i" % (segment_label,) \
###><TA::          , ha='center', va='center', fontdict=font \
###><TA::          , bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2', alpha=1.))
###><TA::#
###><TA::    return im
###><TA::
###><TA::
###><TA::  #
###><TA::  # Main: plot_primary
###><TA::  #
###><TA::  x1d, y1d, d1s, d2s = get_primary_rays(layout, telescope)
###><TA::
###><TA::
###><TA::  if (np.any(pax)==None):
###><TA::    pl.figure(1)
###><TA::    pl.clf()
###><TA::    fg,iax=pl.subplots(ncols=1, nrows=1, num=1)
###><TA::  else:
###><TA::    iax=pax
###><TA::
###><TA::
###><TA::  ntot=0
###><TA::
###><TA::  aux = [[],]*segments.N
###><TA::  if (fcol != ''):
###><TA::    tpfcol=PyC.get_colors(segments.Nfamily, cmap=fcol)
###><TA::    aux = []
###><TA::    for itn in range(segments.N):
###><TA::      aux.append([tpfcol[np.int32(segments.families[itn])-1]])
###><TA::  # Put color legend in GTC layout:
###><TA::    for it_fml in range(segments.Nfamily):
###><TA::      iax.plot([],[],linewidth=5, label=r'Fam.: %i' % (it_fml+1,) \
###><TA::          , color=tpfcol[it_fml])
###><TA::    iax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0) \
###><TA::      , fancybox=True, shadow=True,ncol=2)
###><TA::
###><TA::  elif (scol != ''):
###><TA::    print(" I should not be here!")
###><TA::    
###><TA::    tpfcol=PyC.get_colors(segments.N, cmap=scol)
###><TA::    aux = []
###><TA::    for itn in range(segments.N):
###><TA::      aux.append([tpfcol[itn]])
###><TA::#  # Put color legend in GTC layout:
###><TA::#    for it_fml in range(segments.Nfamily):
###><TA::#      iax.plot([],[],linewidth=5, label=r'Fam.: %i' % (it_fml+1,) \
###><TA::#          , color=tpfcol[it_fml])
###><TA::#    iax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0) \
###><TA::#      , fancybox=True, shadow=True,ncol=2)
###><TA::
###><TA::
###><TA::  if (type(highlight)!=type(None)):
###><TA::    for it_sgm in range(segments.N):
###><TA::      tmp = np.arctan2(segments.ypos[it_sgm], segments.xpos[it_sgm])
###><TA::      if (tmp==highlight):
###><TA::        aux[it_sgm] = ['red',]
###><TA::
###><TA::  for it_sgm in range(segments.N):
###><TA::    if (scol != ''):
###><TA::      wi = segments.sort[it_sgm]
###><TA::      itim = plot_primary_segment(telescope\
###><TA::         , segments.xpos[wi], segments.ypos[wi], iax, aux[it_sgm], segment_label=it_sgm+1)
###><TA::    else:
###><TA::      itim = plot_primary_segment(telescope\
###><TA::         , segments.xpos[it_sgm], segments.ypos[it_sgm], iax, aux[it_sgm])
###><TA::
###><TA::
###><TA::  if (prays != False):
###><TA::
###><TA::    ccol = (1.0,0.,0.)
###><TA::    clw = 2.
###><TA::    calp = 0.6
###><TA::
###><TA::#    for rad in irads[1:]:
###><TA::#      var = 2. * np.pi * np.linspace(0.,1.,1001)
###><TA::#      iax.plot(rad*np.sin(var), rad*np.cos(var) \
###><TA::#          , color=ccol, lw=clw, alpha=calp) 
###><TA::#    for ang in iangs[1:]:
###><TA::#      ca = np.cos(ang)
###><TA::#      sa = np.sin(ang)
###><TA::#      ix = [irads[1]*ca, irads[-1]*ca]
###><TA::#      iy = [irads[1]*sa, irads[-1]*sa]
###><TA::#      iax.plot(ix, iy, color=ccol, lw=clw, alpha=calp) 
###><TA::
###><TA::
###><TA::  iax.set_xlim(-telescope.dim_x, telescope.dim_x)
###><TA::  iax.set_ylim(-telescope.dim_y, telescope.dim_y)
###><TA::
###><TA::#
###><TA::# ---------------------------------------------------------------------
###><TA::#
###><TA::
###><TA::def plot_geometry(telescope, beam, secondary, segments):
###><TA::    #
###><TA::    #
###><TA::    #
###><TA::    # This has been copied from get_geometry!
###><TA::    # If get_geometry changes I have to change this...
###><TA::    # ...accordingly
###><TA::    #
###><TA::    #
###><TA::    #
###><TA::#
###><TA::    layout='polar'
###><TA::    #layout='cartesian'
###><TA::
###><TA::    x1d, y1d, d1s, d2s = get_primary_rays(layout, telescope)
###><TA::
###><TA::
###><TA::
###><TA::# First of all, we remove all the rays that come from behind the secondary:
###><TA::    isin=secondary.f2d.ev(x1d, y1d)
###><TA::    ww = np.where(isin < 0.5)[0]
###><TA::    if (len(ww)>0):
###><TA::      x1d = x1d[ww]
###><TA::      y1d = y1d[ww]
###><TA::
###><TA::#    ax.scatter(x1d,y1d,color='r', s=0.1, edgecolor=None)
###><TA::
###><TA::#
###><TA::    segments.d11d = d1s * 1.
###><TA::    segments.d21d = d2s * 1.
###><TA::#
###><TA::    def check_primary_id(iteles, ixpos, iypos, igxs, igys, fill_nan=False):
###><TA::#
###><TA::      if (iteles.primary_shape == 'hexagonal'):
###><TA::        dumgx = igxs - ixpos
###><TA::        dumgy = igys - iypos
###><TA::        resgx = igxs - ixpos
###><TA::        resgy = igys - iypos
###><TA::  #
###><TA::        a = 30.e0 / 1.8e2 * np.pi
###><TA::        costhU = np.cos(a) * iteles.radius
###><TA::  # Select pixels that fulfil first criterion:
###><TA::        if (fill_nan==True):
###><TA::          ww = np.where(np.abs(dumgy) > costhU)[0]
###><TA::          resgx[ww] = np.nan
###><TA::          resgy[ww] = np.nan
###><TA::        else:
###><TA::          ww = np.where(np.abs(dumgy) <= costhU)[0]
###><TA::          dumgx = dumgx[ww] * 1.
###><TA::          dumgy = dumgy[ww] * 1.
###><TA::          resgx = resgx[ww] * 1.
###><TA::          resgy = resgy[ww] * 1.
###><TA::  # Rotate system 60 deg:
###><TA::        a2 = 60.e0 / 1.8e2 * np.pi
###><TA::        ca = np.cos(a2)
###><TA::        sa = np.sin(a2)
###><TA::        rot = np.array([ca, -sa, sa, ca]).reshape(2,2)
###><TA::  #
###><TA::        out = np.dot(rot, np.vstack([dumgx,dumgy]))
###><TA::        dumxg = out[0,:] * 1.
###><TA::        dumyg = out[1,:] * 1.
###><TA::  # Select pixels that fulfil second criterion:
###><TA::        if (fill_nan==True):
###><TA::          ww = np.where(np.abs(dumyg) > costhU)[0]
###><TA::          resgx[ww] = np.nan
###><TA::          resgy[ww] = np.nan
###><TA::        else:
###><TA::          ww = np.where( np.abs(dumyg) <= costhU)[0]
###><TA::          dumgx = dumgx[ww] * 1.
###><TA::          dumgy = dumgy[ww] * 1.
###><TA::          resgx = resgx[ww] * 1.
###><TA::          resgy = resgy[ww] * 1.
###><TA::  #
###><TA::        rot = np.array([ca, sa, -sa, ca]).reshape(2,2)
###><TA::  #
###><TA::        out = np.dot(rot, np.vstack([dumgx,dumgy]))
###><TA::        dumxg = out[0,:] * 1.
###><TA::        dumyg = out[1,:] * 1.
###><TA::  #
###><TA::        if (fill_nan==True):
###><TA::          ww = np.where(np.abs(dumyg) > costhU)[0]
###><TA::          resgx[ww] = np.nan
###><TA::          resgy[ww] = np.nan
###><TA::        else:
###><TA::          ww = np.where( np.abs(dumyg) <= costhU)[0]
###><TA::          dumgx = dumgx[ww] * 1.
###><TA::          dumgy = dumgy[ww] * 1.
###><TA::          resgx = resgx[ww] * 1.
###><TA::          resgy = resgy[ww] * 1.
###><TA::  #
###><TA::# Anular:
###><TA::      elif (iteles.primary_shape == 'anular'):
###><TA::        dumgx = igxs - ixpos
###><TA::        dumgy = igys - iypos
###><TA::        resgx = igxs - ixpos
###><TA::        resgy = igys - iypos
###><TA::
###><TA::        dumrad = np.sqrt(dumgx**2+dumgy**2)
###><TA::        ww = np.where(dumrad <= iteles.rad_max)[0]
###><TA::
###><TA::        dumrad = dumrad[ww]
###><TA::        resgx = resgx[ww]
###><TA::        resgy = resgy[ww]
###><TA::
###><TA::        ww = np.where(dumrad >= iteles.rad_min)[0]
###><TA::
###><TA::        dumrad = dumrad[ww]
###><TA::        resgx = resgx[ww]
###><TA::        resgy = resgy[ww]
###><TA::
###><TA::# Otherwise:
###><TA::      else:
###><TA::        print('Primary shape unknown!')
###><TA::        print('')
###><TA::        print('')
###><TA::        sys.exit()
###><TA::#
###><TA::      return resgx + ixpos, resgy + iypos
###><TA::
###><TA::#
###><TA::# Starting plot_geometry:
###><TA::#
###><TA::#
###><TA::    from mpl_toolkits.mplot3d import Axes3D
###><TA::    from scipy.interpolate import bisplrep, bisplev
###><TA::    from matplotlib import cm
###><TA::
###><TA::    fignum=2
###><TA::    pl.close(fignum)
###><TA::    fg, ax = pl.subplots(num=fignum,subplot_kw={"projection": "3d"})
###><TA::
###><TA::    if ('GTC' in telescope.ID):
###><TA::      im1max=12.
###><TA::      im2max=15
###><TA::    elif ('ELT' in telescope.ID):
###><TA::      im1max=20.
###><TA::      im2max=25.
###><TA::
###><TA::
###><TA::#
###><TA::#
###><TA::
###><TA::##    plt1 = []
###><TA::##    plt2 = []
###><TA::##    plt3 = []
###><TA::    ntot=0
###><TA::    for it_sgm in range(segments.N):
###><TA::
###><TA::      x, y = check_primary_id(telescope\
###><TA::         , segments.xpos[it_sgm], segments.ypos[it_sgm], x1d, y1d)
###><TA::
###><TA::# Plot:
###><TA::      num=21
###><TA::      #f2d=bisplrep(x,y,z_1)
###><TA::      x1 = np.linspace(x.min(), x.max(), num)
###><TA::      y1 = np.linspace(y.min(), y.max(), num)
###><TA::      #z1=bisplev(x1, y1, f2d)
###><TA::
###><TA::      x1_2d = x1[:,None] * np.ones(num)[None,:]
###><TA::      y1_2d = y1[None,:] * np.ones(num)[:,None]
###><TA::
###><TA::      x1_2d = x1_2d.flatten()
###><TA::      y1_2d = y1_2d.flatten()
###><TA::
###><TA::      x, y = check_primary_id(telescope\
###><TA::         , segments.xpos[it_sgm], segments.ypos[it_sgm], x1_2d, y1_2d, fill_nan=True)
###><TA::
###><TA::
###><TA::      #
###><TA::      # First, we have to rotate the system so that each point is in the plane...
###><TA::      # ...containing the normal and incident ray
###><TA::
###><TA::      th = np.arctan2(y, x)
###><TA::
###><TA::      #
###><TA::      # Rotate system:
###><TA::      cs = np.cos(th)
###><TA::      sn = np.sin(th)
###><TA::      x1 = x * cs + y * sn
###><TA::      y1 = x * (-sn) + y * cs
###><TA::
###><TA::      # Director cosines for input beam:
###><TA::      dir_cos_i = get_direction_cosines(beam.l_inc,beam.m_inc,beam.n_inc)[:,None]*np.ones(x1.size)[None,:]
###><TA::      dir_cos_n1, z_1 = new2_get_gcurv(x1,k=telescope.k_1,r=telescope.rc_1,norm=True)
###><TA::
###><TA::      i_1 = np.arccos(np.sum(dir_cos_n1*dir_cos_i,axis=0))
###><TA::
###><TA::      dir_cos_o1 = rotate_2d(dir_cos_i, -2 * i_1 * np.sign(x1))
###><TA::
###><TA::      x2, z_2, dir_cos_n2 = get_intersection(dir_cos_o1, x1, z_1, telescope.nrc_2, telescope.k_2, telescope.dmm)
###><TA::
###><TA::      i_2 = np.arccos(np.sum(dir_cos_n2*dir_cos_o1,axis=0))
###><TA::
###><TA::      dir_cos_f = rotate_2d(dir_cos_o1, 2 * i_2 * np.sign(x2))
###><TA::
###><TA::      zf=-10.
###><TA::      xf = get_intersection_focus(dir_cos_f, x2, z_2,zf)
###><TA::
###><TA::      aux = i_1[i_1==i_1] * 180. / np.pi
###><TA::      print(" Min(i1)=%.2f ; Max(i1)=%.2f ; <i1>=%.2f in deg." % (aux.min(), aux.max(), aux.mean()))
###><TA::      aux = i_2[i_2==i_2] * 180. / np.pi
###><TA::      print(" Min(i2)=%.2f ; Max(i2)=%.2f ; <i2>=%.2f in deg." % (aux.min(), aux.max(), aux.mean()))
###><TA::      del(aux)
###><TA::
###><TA::      # Following: https://pundit.pratt.duke.edu/wiki/Python:Plotting_Surfaces
###><TA::      def get_color(a, clim): return np.nan_to_num(a)*180./np.pi/clim
###><TA::
###><TA::      x1s = x.reshape(num, num)
###><TA::      y1s = y.reshape(num, num)
###><TA::      cmm = cm.viridis
###><TA::      surf1 = ax.plot_surface(x1s, y1s, z_1.reshape(num, num) \
###><TA::          ,linewidth=0, antialiased=False, facecolors \
###><TA::          =cmm(get_color(i_1.reshape(num,num), im1max)) \
###><TA::          , rcount=51, ccount=51,cmap=cmm)
###><TA::
###><TA::      #
###><TA::      # Back-rotate system:
###><TA::      cs = np.cos(-th)
###><TA::      sn = np.sin(-th)
###><TA::      y2 = x2 * 0.
###><TA::      x_sec = x2 * cs + y2 * sn
###><TA::      y_sec = x2 * (-sn) + y2 * cs
###><TA::      yf = xf * 0.
###><TA::      x_f = xf * cs + yf * sn
###><TA::      y_f = xf * (-sn) + yf * cs
###><TA::
###><TA::
###><TA::      x_sec_2d = x_sec.reshape(num, num)
###><TA::      y_sec_2d = y_sec.reshape(num, num)
###><TA::      z_2_2d = z_2.reshape(num, num)
###><TA::      i_2_2d = i_2.reshape(num, num)
###><TA::      x_f_2d = x_f.reshape(num, num)
###><TA::      y_f_2d = y_f.reshape(num, num)
###><TA::
###><TA::      cmm = cm.magma
###><TA::      surf2 = ax.plot_surface(x_sec_2d, y_sec_2d, z_2_2d \
###><TA::          ,linewidth=0, antialiased=False, facecolors=cmm(get_color(i_2_2d,im2max)) \
###><TA::          , rcount=61, ccount=61, cmap=cmm)
###><TA::
###><TA::
###><TA::###      itpt = x.size//2
###><TA::###      plt1.append([[x[itpt], x[itpt]], [y[itpt],y[itpt]], [telescope.dmm,z_1[itpt]]])
###><TA::###      plt2.append([[x[itpt], x_sec[itpt]], [y[itpt],y_sec[itpt]], [z_1[itpt], z_2[itpt]]])
###><TA::###      plt3.append([[x_sec[itpt], x_f[itpt]], [y_sec[itpt],y_f[itpt]], [z_2[itpt], zf]])
###><TA::
###><TA::      #ax.plot([x[itpt], x[itpt]], [y[itpt],y[itpt]], [telescope.dmm,z_1[itpt]], lw=1)
###><TA::      #ax.plot([x[itpt], x_sec[itpt]], [y[itpt],y_sec[itpt]], [z_1[itpt], z_2[itpt]], lw=1)
###><TA::      #ax.plot([x_sec[itpt], x_f[itpt]], [y_sec[itpt],y_f[itpt]], [z_2[itpt], zf], lw=1)
###><TA::
###><TA::      isin=secondary.f2d.ev(x_sec, y_sec)
###><TA::      ww = np.where(isin > 0.5)[0]
###><TA::
###><TA::      segments.i['%i' % (it_sgm,)]=i_1[ww]
###><TA::      segments.i2['%i' % (it_sgm,)]=i_2[ww]
###><TA::      segments.th1['%i' % (it_sgm,)]=np.arctan2(y, x)[ww]#th_1[ww]
###><TA::
###><TA::      segments.rad['%i' % (it_sgm,)] = np.sqrt(x[ww]**2+y[ww]**2)
###><TA::      segments.area['%i' % (it_sgm,)] = segments.d11d * segments.d21d
###><TA::
###><TA::      ntot=ntot+ww.size
###><TA::    #
###><TA::
###><TA::
###><TA::##    for ipl1, ipl2, ipl3 in zip(plt1, plt2, plt3):
###><TA::##      ax.plot(*ipl1, lw=1)
###><TA::##      ax.plot(*ipl2, lw=1)
###><TA::##      ax.plot(*ipl3, lw=1)
###><TA::
###><TA::    print('Total number of rays considered: %i' % ntot)
###><TA::
###><TA::    cax1 = pl.axes([0.05,0.05,0.01,0.4])
###><TA::    cax2 = pl.axes([0.05,0.55,0.01,0.4])
###><TA::
###><TA::    cbar1 = fg.colorbar(surf1,cax=cax1,boundaries=np.linspace(0,im1max,256)\
###><TA::        ,ticks=np.round(np.linspace(0.,im1max,5)))
###><TA::    cbar1.ax.set_ylabel(r'i$_{1}$')
###><TA::    cbar2 = fg.colorbar(surf2,cax=cax2,boundaries=np.linspace(0,im2max,256)\
###><TA::        ,ticks=np.round(np.linspace(0.,im2max,5)))
###><TA::    cbar2.ax.set_ylabel(r'i$_{2}$')
###><TA::
###><TA::    pl.show()
###><TA::    pl.draw()
###><TA::    stop()
###><TA::
###><TA::    return
###><TA::#
###><TA::#
###><TA::#
###><TA::# ---------------------------------------------------------------------
###><TA::#
###><TA::#
###><TA::#
###><TA::# ---------------------------------------------------------------------
###><TA::#
###><TA::
###><TA::def new_plot_primary(telescope, beam, secondary, segments \
###><TA::    , pax=None, prays=False, fcol='', layout='polar'):
###><TA::#
###><TA::  import PyColors as PyC
###><TA::
###><TA::  def plot_primary_segment(iteles, ixpos, iypos, i2ax, ifcol):
###><TA::#
###><TA::    if (iteles.primary_shape == 'hexagonal'):
###><TA::
###><TA::      dy = np.sqrt(iteles.radius**2-(iteles.radius/2.)**2)
###><TA::
###><TA::      p1x = ixpos + iteles.radius
###><TA::      p1y = iypos + 0.
###><TA::
###><TA::      p2x = ixpos + iteles.radius / 2.
###><TA::      p2y = iypos + dy
###><TA::
###><TA::      p3x = ixpos - iteles.radius / 2.
###><TA::      p3y = iypos + dy
###><TA::
###><TA::      p4x = ixpos - iteles.radius
###><TA::      p4y = iypos + 0.
###><TA::
###><TA::      p5x = ixpos - iteles.radius / 2.
###><TA::      p5y = iypos - dy
###><TA::
###><TA::      p6x = ixpos + iteles.radius / 2.
###><TA::      p6y = iypos - dy
###><TA::
###><TA::      ix = [ p1x, p2x, p3x, p4x, p5x, p6x, p1x ]
###><TA::      iy = [ p1y, p2y, p3y, p4y, p5y, p6y, p1y ]
###><TA::
###><TA::      if (len(ifcol)==0):
###><TA::
###><TA::        im = i2ax.plot(ix,iy,color='k',lw=1.)
###><TA::
###><TA::      else:
###><TA::
###><TA::        #stop()
###><TA::
###><TA::        ix1 = [p4x, p3x, p2x, p1x]
###><TA::        iy1 = [p4y, p3y, p2y, p1y]
###><TA::        iy2 = [p4y, p5y, p6y, p1y]
###><TA::
###><TA::        im = i2ax.fill_between(ix1, iy1, iy2, color=ifcol[0])
###><TA::        im = i2ax.plot(ix,iy,color='k',lw=1.)
###><TA::
###><TA::  #
###><TA::# WRONG!::# Anular:
###><TA::# WRONG!::    elif (iteles.primary_shape == 'anular'):
###><TA::# WRONG!::      dumgx = igxs - ixpos
###><TA::# WRONG!::      dumgy = igys - iypos
###><TA::# WRONG!::      resgx = igxs - ixpos
###><TA::# WRONG!::      resgy = igys - iypos
###><TA::# WRONG!::
###><TA::# WRONG!::      dumrad = np.sqrt(dumgx**2+dumgy**2)
###><TA::# WRONG!::      ww = np.where(dumrad <= iteles.rad_max)[0]
###><TA::# WRONG!::
###><TA::# WRONG!::      dumrad = dumrad[ww]
###><TA::# WRONG!::      resgx = resgx[ww]
###><TA::# WRONG!::      resgy = resgy[ww]
###><TA::# WRONG!::
###><TA::# WRONG!::      ww = np.where(dumrad >= iteles.rad_min)[0]
###><TA::# WRONG!::
###><TA::# WRONG!::      dumrad = dumrad[ww]
###><TA::# WRONG!::      resgx = resgx[ww]
###><TA::# WRONG!::      resgy = resgy[ww]
###><TA::# WRONG!::
###><TA::# WRONG!::      im = []
###><TA::
###><TA::# Otherwise:
###><TA::    else:
###><TA::      print('Primary shape unknown!')
###><TA::      print('')
###><TA::      print('')
###><TA::      sys.exit()
###><TA::#
###><TA::    return im
###><TA::
###><TA::
###><TA::  #
###><TA::  # Main: plot_primary
###><TA::  #
###><TA::
###><TA::  #layout='cartesian'
###><TA::
###><TA::  x1d, y1d, d1s, d2s = get_primary_rays(layout, telescope)
###><TA::#######  if ( (layout!='polar') & (layout!='cartesian') ):
###><TA::#######    print('')
###><TA::#######    print('\t Unknown layout: %s for ray tracing!' % (layout,))
###><TA::#######    print('\t Allowed options are: "polar", "cartesian"')
###><TA::#######    print('')
###><TA::#######    exit(1)
###><TA::
###><TA::  if (np.any(pax)==None):
###><TA::    pl.figure(1)
###><TA::    pl.clf()
###><TA::    fg,iax=pl.subplots(ncols=1, nrows=1, num=1)
###><TA::  else:
###><TA::    iax=pax
###><TA::
###><TA::  ntot=0
###><TA::  if (fcol != ''):
###><TA::    tpfcol=PyC.get_colors(segments.Nfamily, cmap=fcol)
###><TA::    aux = []
###><TA::    for itn in range(segments.N):
###><TA::      aux.append([tpfcol[np.int32(segments.families[itn])-1]])
###><TA::  # Put color legend in GTC layout:
###><TA::    for it_fml in range(segments.Nfamily):
###><TA::      iax.plot([],[],linewidth=5, label=r'Fam.: %i' % (it_fml+1,) \
###><TA::          , color=tpfcol[it_fml])
###><TA::    iax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0) \
###><TA::      , fancybox=True, shadow=True,ncol=2)
###><TA::
###><TA::  else:
###><TA::    aux = [[],]*segments.N
###><TA::
###><TA::  for it_sgm in range(segments.N):
###><TA::    itim = plot_primary_segment(telescope\
###><TA::       , segments.xpos[it_sgm], segments.ypos[it_sgm], iax, aux[it_sgm])
###><TA::
###><TA::
###><TA::  if (prays != False):
###><TA::
###><TA::###    import matplotlib.pyplot as pl
###><TA::###    pl.ion()
###><TA::###    pl.close(1)
###><TA::###    fg, ax = pl.subplots(nrows=1,ncols=1,num=1)
###><TA::
###><TA::# First of all, we remove all the rays that come from behind the secondary:
###><TA::    isin=secondary.f2d.ev(x1d, y1d)
###><TA::
###><TA::    azis = np.arctan2(y1d, x1d)
###><TA::    dsts = np.sqrt(y1d**2+x1d**2)
###><TA::
###><TA::    tmp = np.sort(dsts[np.abs(azis-azis[0])<1.e-10])
###><TA::    dr = (tmp[1:]-tmp[0:-1]).mean() / 2.
###><TA::    dtmp = np.array([dsts[0]])
###><TA::    for it in dsts:
###><TA::      if (np.abs(it-dtmp).min() < dr/2.):
###><TA::        continue
###><TA::      dtmp = np.concatenate([dtmp, np.array([it,])])
###><TA::
###><TA::    x = np.linspace(0., 1., 1001)
###><TA::    for itdn, itdv in enumerate(dtmp):
###><TA::      xcirc = (itdv + dr/2.) * np.cos(2. * np.pi * x / 1.)
###><TA::      ycirc = (itdv + dr/2.) * np.sin(2. * np.pi * x / 1.)
###><TA::      iax.plot(xcirc, ycirc, color='magenta')
###><TA::
###><TA::    tmp = np.sort(azis[np.abs(dsts-dsts[0])<1.e-10])
###><TA::    for itdn, itdv in enumerate(tmp):
###><TA::      x0 = dtmp.min()*np.cos(itdv)
###><TA::      x1 = dtmp.max()*np.cos(itdv)
###><TA::      y0 = dtmp.min()*np.sin(itdv)
###><TA::      y1 = dtmp.max()*np.sin(itdv)
###><TA::      iax.plot([x0, x1], [y0, y1], color='magenta')
###><TA::    
###><TA::
###><TA::    #iax.scatter(x1d,y1d,color='r', marker='o', s=10, edgecolor=None)
###><TA::    ww = np.where(isin < 0.5)[0]
###><TA::    if (len(ww)>0):
###><TA::      x1d = x1d[ww]
###><TA::      y1d = y1d[ww]
###><TA::
###><TA::    #iax.scatter(x1d,y1d,color='k', marker='o', s=7, edgecolor=None)
###><TA::
###><TA::
###><TA::    def check_primary_id(iteles, ixpos, iypos, igxs, igys):
###><TA::#
###><TA::      if (iteles.primary_shape == 'hexagonal'):
###><TA::        dumgx = igxs - ixpos
###><TA::        dumgy = igys - iypos
###><TA::        resgx = igxs - ixpos
###><TA::        resgy = igys - iypos
###><TA::  #
###><TA::        a = 30.e0 / 1.8e2 * np.pi
###><TA::        costhU = np.cos(a) * iteles.radius
###><TA::  # Select pixels that fulfil first criterion:
###><TA::        ww = np.where(np.abs(dumgy) <= costhU)[0]
###><TA::        dumgx = dumgx[ww] * 1.
###><TA::        dumgy = dumgy[ww] * 1.
###><TA::        resgx = resgx[ww] * 1.
###><TA::        resgy = resgy[ww] * 1.
###><TA::  # Rotate system 60 deg:
###><TA::        a2 = 60.e0 / 1.8e2 * np.pi
###><TA::        ca = np.cos(a2)
###><TA::        sa = np.sin(a2)
###><TA::        rot = np.array([ca, -sa, sa, ca]).reshape(2,2)
###><TA::  #
###><TA::        out = np.dot(rot, np.vstack([dumgx,dumgy]))
###><TA::        dumxg = out[0,:] * 1.
###><TA::        dumyg = out[1,:] * 1.
###><TA::  # Select pixels that fulfil second criterion:
###><TA::        ww = np.where( np.abs(dumyg) <= costhU)[0]
###><TA::        dumgx = dumgx[ww] * 1.
###><TA::        dumgy = dumgy[ww] * 1.
###><TA::        resgx = resgx[ww] * 1.
###><TA::        resgy = resgy[ww] * 1.
###><TA::  #
###><TA::        rot = np.array([ca, sa, -sa, ca]).reshape(2,2)
###><TA::  #
###><TA::        out = np.dot(rot, np.vstack([dumgx,dumgy]))
###><TA::        dumxg = out[0,:] * 1.
###><TA::        dumyg = out[1,:] * 1.
###><TA::  #
###><TA::        ww = np.where( np.abs(dumyg) <= costhU)[0]
###><TA::        dumgx = dumgx[ww] * 1.
###><TA::        dumgy = dumgy[ww] * 1.
###><TA::        resgx = resgx[ww] * 1.
###><TA::        resgy = resgy[ww] * 1.
###><TA::  #
###><TA::# Anular:
###><TA::      elif (iteles.primary_shape == 'anular'):
###><TA::        dumgx = igxs - ixpos
###><TA::        dumgy = igys - iypos
###><TA::        resgx = igxs - ixpos
###><TA::        resgy = igys - iypos
###><TA::
###><TA::        dumrad = np.sqrt(dumgx**2+dumgy**2)
###><TA::        ww = np.where(dumrad <= iteles.rad_max)[0]
###><TA::
###><TA::        dumrad = dumrad[ww]
###><TA::        resgx = resgx[ww]
###><TA::        resgy = resgy[ww]
###><TA::
###><TA::        ww = np.where(dumrad >= iteles.rad_min)[0]
###><TA::
###><TA::        dumrad = dumrad[ww]
###><TA::        resgx = resgx[ww]
###><TA::        resgy = resgy[ww]
###><TA::
###><TA::# Otherwise:
###><TA::      else:
###><TA::        print('Primary shape unknown!')
###><TA::        print('')
###><TA::        print('')
###><TA::        sys.exit()
###><TA::#
###><TA::      return resgx + ixpos, resgy + iypos
###><TA::
###><TA::#
###><TA::    ntot=0
###><TA::
###><TA::    x1df = x1d * 0.
###><TA::    y1df = y1d * 0.
###><TA::
###><TA::    for it_sgm in range(segments.N):
###><TA::
###><TA::      x, y = check_primary_id(telescope\
###><TA::         , segments.xpos[it_sgm], segments.ypos[it_sgm], x1d, y1d)
###><TA::
###><TA::      #
###><TA::      # First, we have to rotate the system so that each point is in the plane...
###><TA::      # ...containing the normal and incident ray
###><TA::
###><TA::      th = np.arctan2(y, x)
###><TA::
###><TA::      #
###><TA::      # Rotate system:
###><TA::      cs = np.cos(th)
###><TA::      sn = np.sin(th)
###><TA::      x1 = x * cs + y * sn
###><TA::      y1 = x * (-sn) + y * cs
###><TA::
###><TA::      # Director cosines for input beam:
###><TA::      dir_cos_i = get_direction_cosines(beam.l_inc,beam.m_inc,beam.n_inc)[:,None]*np.ones(x1.size)[None,:]
###><TA::      dir_cos_n1, z_1 = new2_get_gcurv(x1,k=telescope.k_1,r=telescope.rc_1,norm=True)
###><TA::
###><TA::      i_1 = np.arccos(np.sum(dir_cos_n1*dir_cos_i,axis=0))
###><TA::
###><TA::      dir_cos_o1 = rotate_2d(dir_cos_i, -2 * i_1 * np.sign(x1))
###><TA::
###><TA::      x2, z_2, dir_cos_n2 = get_intersection(dir_cos_o1, x1, z_1, telescope.nrc_2, telescope.k_2, telescope.dmm)
###><TA::
###><TA::      i_2 = np.arccos(np.sum(dir_cos_n2*dir_cos_o1,axis=0))
###><TA::
###><TA::      #
###><TA::      # Back-rotate system:
###><TA::      cs = np.cos(-th)
###><TA::      sn = np.sin(-th)
###><TA::      y2 = x2 * 0.
###><TA::      x_sec = x2 * cs + y2 * sn
###><TA::      y_sec = x2 * (-sn) + y2 * cs
###><TA::
###><TA::      isin=secondary.f2d.ev(x_sec, y_sec)
###><TA::      ww = np.where(isin > 0.5)[0]
###><TA::
###><TA::      x1df[ntot:ntot+ww.size] = x[ww] * 1.
###><TA::      y1df[ntot:ntot+ww.size] = y[ww] * 1.
###><TA::
###><TA::      ntot=ntot+ww.size
###><TA::
###><TA::    ##iax.scatter(x1df[0:ntot],y1df[0:ntot],color='cyan', marker='o', s=4, edgecolor=None)
###><TA::    print(np.arctan2(x1df[0:ntot],y1df[0:ntot])*180./np.pi)
###><TA::
###><TA::
###><TA::    print('Total number of rays considered: %i' % ntot)
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::
###><TA::  return
###><TA::
###><TA::
###><TA::def plot_mueller_elements(segments, materials,cdust=False):
###><TA::#
###><TA::  tarea = 0.
###><TA::  for it in range(segments.N):
###><TA::    tarea += np.sum(segments.area['%i' % (it,)])
###><TA::
###><TA::#RUN TIMES:
###><TA::  res = np.zeros((segments.SimulatedTimes, segments.N, 4, 4))
###><TA::  avg_res = np.zeros((segments.SimulatedTimes, 4, 4))
###><TA::  for it_tms in range(segments.SimulatedTimes):
###><TA::    print("Timestep: %i of %i" % (it_tms, segments.SimulatedTimes))
###><TA::#RUN SEGMENTS:
###><TA::    npts_seg = np.zeros((segments.N,))
###><TA::    ntot2=0
###><TA::
###><TA::#
###><TA::# Testing:
###><TA::    pl.close(1)
###><TA::    pl.close(2)
###><TA::    pl.close(3)
###><TA::    pl.close(4)
###><TA::    #fg, ax = pl.subplots(num=1, nrows=6, ncols=6)
###><TA::    #ax=ax.flatten()
###><TA::    fg1, ax1 = pl.subplots(num=1, nrows=1, ncols=1)
###><TA::    fg2, ax2 = pl.subplots(num=2, nrows=1, ncols=1)
###><TA::    fg3, ax3 = pl.subplots(num=3, nrows=1, ncols=1)
###><TA::    fg4, ax4 = pl.subplots(num=4, nrows=1, ncols=1)
###><TA::# Testing.
###><TA::#
###><TA::
###><TA::
###><TA::    for it_sgm in range(segments.N):
###><TA::#
###><TA::      if (segments.i['%i' % (it_sgm,)].size<1):
###><TA::        npts_seg[it_sgm] = 0.
###><TA::        res[it_tms,it_sgm,:,:] = res[it_tms,it_sgm,:,:] * 0.
###><TA::        continue
###><TA::
###><TA::      mmpa_full = hex_syst(segments.i['%i' % (it_sgm,)], \
###><TA::          segments.th1['%i' % (it_sgm,)], \
###><TA::          segments.i2['%i' % (it_sgm,)], \
###><TA::          materials.lamb, materials.n_s, materials.n_c\
###><TA::          , materials.n_ox, materials.n_d, materials.n_in, \
###><TA::          materials.d_s[it_sgm, it_tms], \
###><TA::          materials.d_c[it_sgm, it_tms]\
###><TA::          , materials.d_ox[it_sgm, it_tms], materials.d_d[it_sgm, it_tms])
###><TA::
###><TA::      #GET Mdust
###><TA::      if (cdust==True):
###><TA::        ###mmpd = get_mdust(materials.n_d, materials.d_d[it_sgm, it_tms])
###><TA::        mmpd = get_mdust(0.1, materials.d_d[it_sgm, it_tms])
###><TA::        res_full = mmpa_full.dot(mmpd)
###><TA::      else:
###><TA::        res_full = mmpa_full * 1.
###><TA::
###><TA::
###><TA::#
###><TA::# Testing:
###><TA::
###><TA::      xx = segments.rad['%i' % (it_sgm,)]*np.cos(segments.th1['%i' % (it_sgm,)])
###><TA::      yy = segments.rad['%i' % (it_sgm,)]*np.sin(segments.th1['%i' % (it_sgm,)])
###><TA::      im1 = ax1.scatter(xx, yy, marker='s', s=10, c=res_full[:,1,0] * segments.area['%i' % (it_sgm,)] / tarea, cmap='RdGy', vmin=-0.0000001, vmax=0.0000001)
###><TA::      im2 = ax2.scatter(xx, yy, marker='s', s=10, c=res_full[:,2,0] * segments.area['%i' % (it_sgm,)] / tarea, cmap='RdGy', vmin=-0.0000001, vmax=0.0000001)
###><TA::      im3 = ax3.scatter(xx, yy, marker='s', s=10, c=res_full[:,1,0], cmap='RdGy', vmin=-0.0001, vmax=0.0001)
###><TA::      im4 = ax4.scatter(xx, yy, marker='s', s=10, c=res_full[:,2,1], cmap='RdGy', vmin=-0.0001, vmax=0.0001)
###><TA::      #ax[it_sgm].hist(res_full[:,1,0], bins=101)
###><TA::# Testing.
###><TA::#
###><TA::
###><TA::
###><TA::      nn,_,_ = res_full.shape
###><TA::      ntot2=ntot2+nn
###><TA::########      for itn in range(nn):
###><TA::########        res[it_tms,it_sgm,:,:]=res[it_tms,it_sgm,:,:] \
###><TA::########            + res_full[itn,:,:] * segments.rad['%i' % (it_sgm,)][itn]
###><TA::########      res[it_tms,it_sgm,:,:] = res[it_tms,it_sgm,:,:] \
###><TA::########          * segments.area['%i' % (it_sgm,)]
###><TA::########      #
###><TA::########      # Area to normalize:
###><TA::########      npts_seg[it_sgm] = np.sum(segments.rad['%i' % (it_sgm,)] \
###><TA::########          * segments.area['%i' % (it_sgm,)])
###><TA::
###><TA::      for itn in range(nn):
###><TA::        res[it_tms,it_sgm,:,:]=res[it_tms,it_sgm,:,:] \
###><TA::            + res_full[itn,:,:] * segments.area['%i' % (it_sgm,)][itn]
###><TA::      #
###><TA::      # Area to normalize:
###><TA::      npts_seg[it_sgm] = np.sum(segments.area['%i' % (it_sgm,)])
###><TA::
###><TA::
###><TA::    print('Total number of rays considered: %i' % ntot2)
###><TA::    print('Total area considered: %.4f==%.4f?' % (npts_seg.sum(), tarea,))
###><TA::    
###><TA::
###><TA::    # Since we consider all the points of each segment, the average must...
###><TA::    # ... be done accordingly:
###><TA::    avg_res[it_tms,:,:] = np.sum(res[it_tms,:,:,:], axis=0) \
###><TA::        / np.sum(npts_seg)
###><TA::
###><TA::#
###><TA::  cbar1 = pl.colorbar(im1, ax=ax1)
###><TA::  cbar2 = pl.colorbar(im2, ax=ax2)
###><TA::  cbar3 = pl.colorbar(im3, ax=ax3)
###><TA::  cbar4 = pl.colorbar(im4, ax=ax4)
###><TA::  for ax in [ax1, ax2, ax3, ax4]:
###><TA::    ax.set_aspect('equal')
###><TA::    ax.grid()
###><TA::  pl.draw()
###><TA::  stop()
###><TA::  return avg_res, res, npts_seg
###><TA::
###><TA::
###><TA::#
###><TA::# ---------------------------------------------------------------------
###><TA::#
###><TA::
###><TA::
###><TA::
###><TA::##########
