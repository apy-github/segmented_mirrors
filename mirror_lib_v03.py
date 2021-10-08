#!/usr/bin/env python3.8
#
import time as tm
import numpy as np
#                                                                              #
from pdb import set_trace as stop
import matplotlib.pyplot as pl
pl.ion()
#
def init_telescope(name, num):
#
  import numpy as np
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
        import sys
        sys.exit()

      return

  avail = ['EELT', 'GTC', 'SGTC', 'AGTC', 'SEELT', 'AEELT']
  isavail = np.zeros(len(avail))
  for itn, itc in enumerate(avail):
    if (name.upper() == itc):
      isavail[itn]=1

  if ( np.sum(isavail) < 0.5):
#(name.upper() != 'EELT') & (name.upper() != 'GTC')  & (name.upper() != 'SGTC') ):
    print("only GTC or EELT are available!")
    import sys
    sys.exit()

  teles = telescope(name, num)

  return teles
#
#                                                                              #
#
def init_beam(alpha,x_alpha,y_alpha,degrees=False):

  import numpy as np

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
  import numpy as np
  from pdb import set_trace as stop
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
                print("")
                print("\t Warning!, number of segments inside maximum radius is too big")
                print(" %i -> %i" % (inum_esp,rads.size,))
                print("")

                if (inum_esp<rads.size):
                  ww = np.argsort(rads)
                  rads = rads[ww[0:inum_esp]]
                  x = x[ww[0:inum_esp]]
                  y = y[ww[0:inum_esp]]
                else:
                  inum_esp=np.min([inum_esp, rads.size])
       ###??:     #
       ###??:     #LAST MIRROR
       ###??:     #
       ###??:       last = rads[inum_esp-1]
       ###??:       ww = np.where(rads == last)[0]
       ###??:       if ( np.max(ww) > inum_esp-1 ):
       ###??:         import pdb as pdb
       ###??:         pdb.set_trace()
       ###??:         inum_esp = np.max(ww) + 1
       ###??:       rads = rads[0:inum_esp]
       ###??:       x = x[0:inum_esp]
       ###??:       y = y[0:inum_esp]
  

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
          exit()
        #
        # May be I should change this limitation:
        if (mltch<1):
          print("")
          print("\tMultiplechange MUST be >=1!")
          print("")
          exit()
        if (self.N%mltch!=0):
          print("")
          print("\tMultiplechange MUST give MODULUS(N,MLT)==0")
          print("")
          exit()
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
          exit()
        #
        self.SimulatedTimes = self.Ntimes * 1
        self.SimulatedTime = self.tlong * 1
        print('\tNumber of simulated steps: %i' % (self.Ntimes, ))
        print('\tNumber of simulated steps in between mirror changes: %i' \
            % (self.days_for_each_timestep, ))

        #self.time_map = np.zeros((self.N, self.Ntimes, 2))
# DOES IT WORK AT ALL?
        self.time_map = np.zeros((self.N, self.SimulatedTimes, 3), dtype=np.float64)
        #self.cleandustcadence = cleandust * 1.
        # Caution, this variable (cleandustcadence) comes in days...
        # ... and it has to be transformed to simultion steps:
        self.cleandustcadence = cleandust / self.tstep
        #from pdb import set_trace as stop
        #stop()

      #
        return
      #
    def get_order(self, order):
      from pdb import set_trace as stop
      if (order.upper() == 'RANDOM'):
        sort = np.zeros((self.N,), dtype=np.int32)
        np.random.seed(1)
        fsort = np.argsort(np.random.randn(self.Nfamily))+1
        for it_nnn in range(self.Nfamily):
          ww=np.where(np.abs(self.families-fsort[it_nnn])<1.e-3)[0]
          sort[ww] = (np.argsort(np.random.randn(6))) * self.Nfamily + it_nnn
          #print(sort[ww])
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
          #while (np.size(itf_sort) < np.size(ww)):
          for it_nn2 in range(np.size(ww)-1):
            it2_a = iff_azi[-1]
            dif = (azi360 - it2_a + 360) % 360
            #print dif, it2_a
            #print it_index
            adif = np.abs(dif - 120.)
            ww2 = np.where(adif == np.min(adif))[0][0]
            #print ww2
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
          #stop()
          sort = np.hstack([sort, ref_pos[ww][it_sort]])

      elif (order.upper() == 'AZIMUTHAL'):
        #sort = np.argsort(self.azimuth+180./np.pi)[::-1]
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
        #nazim = np.arctan2(self.xpos, self.ypos)
        #sort = np.argsort(nazim)

        #sort = np.arange(self.N)
      elif (order.upper() == 'EQUAL'):
        sort = np.int32(self.azimuth * 0)
      else:
        print("THE ORDER MUST BE ONE OF THE FOLLOWINGS:")
        print("\tRANDOM")
        print("\tSYMMETRIC")
        print("\tLINEAL")
        print("\tAZIMUTHAL")
        import sys as sys
        sys.exit()
      #.stop()
      #.print self.families[sort]
#!><#.      import matplotlib.pyplot as pl
#!><#.      pl.ion()
#!><#.      pl.clf()
#!><#.      for i in range(self.N):
#!><#.        pl.plot([self.xpos[sort][i]], [self.ypos[sort][i]], label='%i' % (i+1), marker='o',color=((i*1.)/(self.N-1.),1.-(i*1.)/(self.N-1.),0.5))
#!><#.        pl.text(self.xpos[sort][i], self.ypos[sort][i], '%i' % (i+1))
#!><#.      pl.legend()
#!><#.      #pl.imshow(self.time_map[:,:,1].T)
#!><#.      from pdb import set_trace as stop
#!><#.      stop()
#AQUI

    #
      self.sort = sort * 1
#>< CHECK SORTING      pl.ion()
#>< CHECK SORTING      pl.clf()
#>< CHECK SORTING      for it_nnn in range(self.ypos.size):
#>< CHECK SORTING        pl.text(self.xpos[self.sort[it_nnn]], self.ypos[self.sort[it_nnn]], '%i' % (it_nnn+1))
#>< CHECK SORTING        print('i=%i ; azi[i]=%.2f ; r[i]=%.2f ; family[i]=%i' % (it_nnn, self.azimuth[self.sort[it_nnn]], self.dist[self.sort[it_nnn]], self.families[self.sort[it_nnn]], ))
#>< CHECK SORTING      pl.axis('equal')
#>< CHECK SORTING      pl.xlim(-40,40)
#>< CHECK SORTING      pl.ylim(-40,40)
#>< CHECK SORTING      pl.draw()
#>< CHECK SORTING      pl.show()
#>< CHECK SORTING      stop()
    #
      return
    #
    def get_times_new(self, ideal):
    #
      import matplotlib.pyplot as pl
      pl.ion()
      from pdb import set_trace as stop

      if (ideal!=True):

        #to_be_rolled = np.linspace(0.,self.period, self.Ntimes)
        to_be_rolled = np.linspace(0.,self.SimulatedTime \
            , self.SimulatedTimes) % self.period

# Dust:
        frac = 0.099999
  #      frac = 0.30

        cnt = 0
        while (cnt<self.N):

          for i in range(self.mltch):

            offset = self.days_between_change \
                * (self.deltat/self.tstep * (self.sort[cnt] // self.mltch))
            #stop()

            self.time_map[cnt,:,0] = (np.arange(self.Ntimes \
                , dtype=np.float64)*self.tstep+offset) % self.period

            self.time_map[cnt,:,1] = np.min(np.vstack([self.time_map[cnt \
                ,:,0], (np.arange(self.Ntimes,dtype=np.float64) \
                %self.cleandustcadence)*self.tstep]),axis=0)

#            self.time_map[cnt,:,2] = ((1.+frac)**np.floor( \
#                (self.time_map[cnt,:,0] / self.tstep) \
#                / (1.*self.cleandustcadence)) - 1.) \
#                * (self.cleandustcadence * self.tstep)
#                = 

            cnt = cnt+1

        import apy_utils as uts

#
# Take into account not complete cleaning after using carbon snow:
        #
        # First: where are the closest to the segment exchange:
#><        mdiff = self.time_map[:,1:,0]-self.time_map[:,0:-1,0]
#><        schanges = []
#><        for its in range(self.N):
#><          schanges.append(np.where(mdiff[its,:] < 0)[0] + 1)
#><        #


#
#
#        #
#        # First: Simulation steps since last segment change:
#        simsteplastsch = self.time_map[0,:,0]/self.tstep
#        #
#        # Second: Number of cleanings since last change:
#        numclnlastsch = simsteplastsch / self.cleandustcadence
#        #
#        # Third: simulation steps since we last clean:
#        dum = np.arange(self.Ntimes) % self.cleandustcadence
#        #
#        # Fourth: days since last cleaning:
#        dum2 = dum * self.tstep
#        #
#
#
#

        initial_num_cln_since_last_sch = (self.time_map[:,0,0] \
            / self.tstep) / self.cleandustcadence

        inital_amount_dust_in_sim_steps = initial_num_cln_since_last_sch \
            - np.floor(initial_num_cln_since_last_sch)

        first_acc_amount_dust_in_sim_steps \
            = inital_amount_dust_in_sim_steps * frac 

        #
        # Accumulate initial dust for every segment:
        aux = np.floor(initial_num_cln_since_last_sch)
        ginit_dust_val = aux * 0. + first_acc_amount_dust_in_sim_steps \
            * self.tstep
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

###          uts.plot(self.time_map[it_sgm, :, 2], vmin = 0)
###          uts.plot(self.time_map[it_sgm, :, 1], noerase=1, color='b')
###  
###          stop()



#
#
#


# PLOT:        uts.tv(self.time_map[:,0:20,0], cmap='viridis', bar=1, num=1)
# PLOT:        uts.tv(self.time_map[:,0:20,1], cmap='viridis', bar=1, num=2)
# PLOT:        uts.tv(self.time_map[:,0:20,2], cmap='viridis', bar=1, num=3)
# PLOT:        stop()


        #><  print('it=%i' % it_ttt)
        #><  print('D(0)= %.2f' % (first_acc_amount_dust_in_sim_steps[0] * self.tstep,))
        #><  print('D(t>t0)= %.2f' % (accum_amount_dust_in_sim_steps[0]* self.tstep,))
        #><  stop()

#        for it_ttt in range(self.Ntimes):
#
#          sim_steps_since_zero = self.time_map[0,it_ttt,0]/self.tstep
#          sim_steps_since_lastcln = dum[it_ttt]
#
#
#
#
#          print('Sim. steps since exchange: %.2f ; Idem since last cleaning: %.2f' % (sim_steps_since_zero, sim_steps_since_lastcln))
#          print(it_ttt==schanges[0][0])


        #

#
# End carbon snow.
#




        #import matplotlib.pyplot as pl
        #pl.ion()
        #pl.close(1)
        #pl.figure(1)
        #im=pl.imshow(self.time_map[:,:,1].T)
        #cbar=pl.colorbar(im)
        #pl.gca(). invert_yaxis()
        #pl.close(2)
        #pl.figure(2)
        #im=pl.imshow(self.time_map[:,:,0].T)
        #cbar=pl.colorbar(im)
        #pl.gca(). invert_yaxis()
        #from pdb import set_trace as stop
        #stop()

    #
      return
    #

    def get_times(self, ideal):
    #
      import matplotlib.pyplot as pl
      pl.ion()
      from pdb import set_trace as stop

      if (ideal!=True):

        #to_be_rolled = np.linspace(0.,self.period, self.Ntimes)
        to_be_rolled = np.linspace(0.,self.SimulatedTime, self.SimulatedTimes) % self.period
        for i in range(self.N):
  
          offset = self.days_between_change * self.sort[i]
          #print offset, i
          self.time_map[i,:,0] = (np.arange(self.Ntimes, dtype=np.float64)*self.tstep+offset) % self.period
  
          #self.time_map[i,:,1] = self.time_map[i,:,0] * 1.
          #self.time_map[i,:,1] = self.time_map[i,:,0] % self.cleandustcadence
          self.time_map[i,:,1] = np.min(np.vstack([self.time_map[i,:,0], (np.arange(self.Ntimes,dtype=np.float64)*self.tstep)%self.cleandustcadence]),axis=0)

      #><  dum=np.min(np.vstack([self.time_map[i,:,0], np.arange(self.Ntimes,dtype=np.float64)%self.cleandustcadence]),axis=0)

      #><  from pdb import set_trace as stop
      #><  import apy_utils as uts
      #><  uts.tv(self.time_map[:,:,0],cmap='viridis',bar=1)
      #><  stop()
      #><  uts.tv(self.time_map[:,:,1],cmap='viridis',bar=1)
      #><  stop()


 #>< #ALREADY CHECKED
 #><       #DUST: CLEANED EVERY CLEANDUSTCADENCE:
 #><       to_be_rolled = np.linspace(0.,self.SimulatedTime, self.SimulatedTimes)
 #><       store = to_be_rolled % self.cleandustcadence
 #><       snow_clean = np.where(store[1:]-store[0:-1] < 0)[0] + 1
 #><       if (np.size(snow_clean) != 0):
 #><         if (snow_clean[0]!=0):
 #><           snow_clean = np.hstack([0, snow_clean])
 #><         if (snow_clean[-1]!=self.Ntimes):
 #><           snow_clean = np.hstack([snow_clean,self.SimulatedTimes])
 #><         for i in range(snow_clean.size-1):
 #><           it=snow_clean[i]*1
 #><           it1=snow_clean[i+1]*1
 #><           to_be_subtracted = self.time_map[:,it,1] * 1. - store[i]
 #><           corrected = self.time_map[:,it:it1,1] - to_be_subtracted[:,None]
 #><           self.time_map[:,it:it1,1] = corrected * 1
 #><         # UP TO HERE IT IS ALMOST RIGHT, THE ONLY THING IS THAT FOR THOSE 
 #><         # TIMES WHEN A MIRROR IS SUBSTITUTED, I GOT TIMES BELOW ZERO
 #><         # I CORRECT THEM BY LOOKING FOR THE MIRROR CHANGING TIMES:
 #><         for i in range(self.N):
 #><           ww = np.where(self.time_map[i,:,1] < 0)[0]
 #><           if (np.size(ww) == 0):
 #><             continue
 #><           self.time_map[i,ww,1] = self.time_map[i,ww,0] * 1.
# It seems to work:        import matplotlib.pyplot as pl
# It seems to work:        pl.ion()
# It seems to work:        pl.close(1)
# It seems to work:        pl.figure(1)
# It seems to work:        im=pl.imshow(self.time_map[:,:,1].T)
# It seems to work:        cbar=pl.colorbar(im)
# It seems to work:        pl.gca(). invert_yaxis()
# It seems to work:        pl.close(2)
# It seems to work:        pl.figure(2)
# It seems to work:        im=pl.imshow(self.time_map[:,:,0].T)
# It seems to work:        cbar=pl.colorbar(im)
# It seems to work:        pl.gca(). invert_yaxis()
# It seems to work:        from pdb import set_trace as stop
# It seems to work:        stop()
    #
      return
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
  #segments.get_times(ideal)
  segments.get_times_new(ideal)
#
  return segments
#

def secondary(telescope):
#
  import numpy as np
#
  def it_hexagon(radius, xx, yy, itf_xpos, itf_ypos):
  #
    import numpy as np
  #
    it_xpos = 0.
    it_ypos = 0.
  #
  #WE ONLY WORK WITH A SMALL PART OF THE WHOLE ARRAY:
  #
    gnum, _ = xx.shape
    #size_to_cut = np.int(np.ceil(radius / (xx[1,1] - xx[0,0]) * 1.2))
    sxx = xx - itf_xpos
    syy = yy - itf_ypos
  #
    num, _ = sxx.shape
  #
    a = 30.e0 / 1.8e2 * np.pi
    costhU = it_ypos + np.cos(a) * radius
    costhL = it_ypos - np.cos(a) * radius
    hex_1 = sxx * 0.
    ww = np.where( (syy < costhU) & (syy > costhL) )
    hex_1[ww[0], ww[1]] = 1.
  #
    a2 = 60.e0 / 1.8e2 * np.pi
    ca = np.cos(a2)
    sa = np.sin(a2)
    rot = np.array([ca, -sa, sa, ca]).reshape(2,2)
  #
    out = np.dot(rot, np.vstack([sxx.reshape(num*num),syy.reshape(num*num)]))
    xxn = out[0,:].reshape(num,num)
    yyn = out[1,:].reshape(num,num)
    it_nxpos, it_nypos = np.dot(rot, np.array([it_xpos, it_ypos]))
  #
    costhU = it_nypos + np.cos(a) * radius
    costhL = it_nypos - np.cos(a) * radius
    hex_2 = xxn * 0.
    ww = np.where( (yyn < costhU) & (yyn > costhL) )
    hex_2[ww[0], ww[1]] = 1.
  #
    rot = np.array([ca, sa, -sa, ca]).reshape(2,2)
  #
    out = np.dot(rot, np.vstack([sxx.reshape(num*num),syy.reshape(num*num)]))
    xxn = out[0,:].reshape(num,num)
    yyn = out[1,:].reshape(num,num)
    it_nxpos, it_nypos = np.dot(rot, np.array([it_xpos, it_ypos]))
  #
    costhU = it_nypos + np.cos(a) * radius
    costhL = it_nypos - np.cos(a) * radius
    hex_3 = xxn * 0.
    ww = np.where( (yyn < costhU) & (yyn > costhL) )
    hex_3[ww[0], ww[1]] = 1.
  #
    hexag = hex_1 * hex_2 * hex_3
  #
    return hexag
#
  class mirror_obj(object):

    def __init__(self,telescope):
      vmax = telescope.sec_ext_radius * np.sqrt(telescope.sec_num_esp) * 1.1
      if (telescope.ID == 'GTC'):
        self.ID = 'GTC'
        self.x = np.linspace(-vmax, vmax, telescope.sec_num)
        self.y = np.linspace(-vmax, vmax, telescope.sec_num)
        self.mirror = np.zeros((telescope.sec_num, telescope.sec_num))
      if (telescope.ID == 'SGTC'):
        self.ID = 'SGTC'
        self.x = np.linspace(-vmax, vmax, telescope.sec_num)
        self.y = np.linspace(-vmax, vmax, telescope.sec_num)
        self.mirror = np.zeros((telescope.sec_num, telescope.sec_num))
      if (telescope.ID == 'AGTC'):
        self.ID = 'AGTC'
        self.x = np.linspace(-vmax, vmax, telescope.sec_num)
        self.y = np.linspace(-vmax, vmax, telescope.sec_num)
        self.mirror = np.zeros((telescope.sec_num, telescope.sec_num))
      elif (telescope.ID == 'EELT'):
        self.ID = 'EELT'
        self.x = np.linspace(-vmax, vmax, telescope.sec_num)
        self.y = np.linspace(-vmax, vmax, telescope.sec_num)
        self.mirror = np.zeros((telescope.sec_num, telescope.sec_num))
      elif (telescope.ID == 'SEELT'):
        self.ID = 'SEELT'
        self.x = np.linspace(-vmax, vmax, telescope.sec_num)
        self.y = np.linspace(-vmax, vmax, telescope.sec_num)
        self.mirror = np.zeros((telescope.sec_num, telescope.sec_num))
      elif (telescope.ID == 'AEELT'):
        self.ID = 'AEELT'
        self.x = np.linspace(-vmax, vmax, telescope.sec_num)
        self.y = np.linspace(-vmax, vmax, telescope.sec_num)
        self.mirror = np.zeros((telescope.sec_num, telescope.sec_num))

      return
#
    def get_secondary_geometry(self,inum_esp, isec_ext_radius \
        , isec_int_radius,ishape, ix, iy):

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
          import sys as sys
          sys.exit()
        rads = rads[0:inum_esp]
        x = x[0:inum_esp]
        y = y[0:inum_esp]
        
        ones = self.x * 0. + 1.
        _, xxax = np.meshgrid(ones, self.x)
        yyax, _ = np.meshgrid(self.y, ones)
        mir = 0.
    
        for it_hxs in range(x.size):
          # THIS FACTOR IS TO NOT HAVE INTERSTICES BETWEEN HEXAGONS
          mir = mir + it_hexagon(isec_ext_radius*1.001, xxax, yyax, x[it_hxs], y[it_hxs])
        
        internal = it_hexagon(isec_int_radius, xxax, yyax, 0., 0.).T
        
        mirror = mir*np.abs(1.-internal)
      #
        ww = np.where(mirror != 0.)
        mirror = mirror * 0.
        mirror[ww[0],ww[1]] = 1.
      elif (ishape=='anular'):
        yy, xx = np.meshgrid(ix, iy)
        mirror = xx * 0.
        radius = np.sqrt(xx**2+yy**2)
        ww = np.where( (radius > isec_int_radius ) & \
            ( radius < isec_ext_radius ) )
        mirror[ww] = 1.

      else:
        print('Secondary:')
        print('\tShape unknown!')
        print('')
        stop()

      return mirror

    def get_mirror(self, telescope):
#

      mirror = self.get_secondary_geometry(telescope.sec_num_esp, telescope.sec_ext_radius \
          , telescope.sec_int_radius, telescope.secondary_shape, self.x, self.y)

      #from pdb import set_trace as stop
      #import matplotlib.pyplot as pl
      #pl.ion()
      #pl.imshow(mirror.T)
      #stop()

      self.mirror = mirror * 1.
      return
    #
    def get_interpolation(self):
      #
      from scipy.interpolate import RectBivariateSpline as interpolate
      self.f2d = interpolate(secondary.x,secondary.y,secondary.mirror)
#
  secondary = mirror_obj(telescope)
  secondary.get_mirror(telescope)

  secondary.get_interpolation()

#  from pdb import set_trace as stop
#  import matplotlib.pyplot as pl
#  pl.ion()
#  pl.imshow(secondary.mirror)
#  stop()

  return secondary
#
#
#
# ---------------------------------------------------------------------
#

#
def dirt_timing(lamb, time):
#
  import numpy as np
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
      print(lam[0], lam_in, lam[-1])
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
#      self.n_s = np.zeros(1, dtype=complex)
#      self.n_s.real = n_s * 1.
      self.n_s = n_s + 0j
#      self.n_in = np.zeros(1, dtype=complex)
#      self.n_in.real = n_in * 1.
      self.n_in = n_in + 0j
#      self.n_d  = np.zeros(1, dtype=complex)
#      self.n_d.real = n_d * 1.
      self.n_d = n_d + 0j
      self.n_c = n_c * 1.
      self.n_ox = n_ox + 0j

      #self.d_ox  = np.zeros(d_ox.shape, dtype=complex)
      #self.d_ox.real = d_ox * 1.
      self.d_ox = d_ox * 1.
      #self.d_c = np.zeros(d_c.shape, dtype=complex)
      #self.d_c.real = d_c * 1.
      self.d_c = d_c * 1.
      #self.d_s = np.zeros(d_s.shape, dtype=complex)
      #self.d_s.real = d_s * 1.
      self.d_s = d_s * 1.
      #self.d_d = np.zeros(d_d.shape, dtype=complex)
      #self.d_d.real = d_d * 1.
      self.d_d = d_d * 1.

      self.lamb = np.zeros(1, dtype=np.complex64)
      self.lamb.real = lamb * 1.
      return
#
  #from pdb import set_trace as stop
  #stop()
  #FIRST APPROACH: ONE DUST PARTICLE THICK PER DAY
  # AlPhA
# Dust:
  #alpha = 0.06
  #alpha = 0.0006
  #alpha = 0.008
  #alpha = 0.0006
  alpha = 0.003#25

  n_s = zerodux(lamb * 1.e-4)
  n_in = 1.e0
  n_d = 1.52e0
  elem_esp = 'al'
  n_g, k_g = indexes(elem_esp,lamb)
  n_c = np.complex( n_g , k_g )
  elem_ox = 'ox'
  n_ox, k_ox = indexes(elem_ox,lamb)
  n_ox = np.complex( n_ox , k_ox )
  
  tau_o = 1392.1484e0
  tau_d = 1.8e2
  th_m_o = 1.2e3
  th_m_d = 1.2e4
  d_mir = 1.2e3   #angstroms

  d_ox = 1. * thick( th_m_o , tau_o , time[:,:,0] ) 
  d_c = (d_mir - d_ox)  #el espesor del conductor sera el original (d_mir) menos lo que se ha oxidado.
  d_d = 1. * Dthick ( alpha, time[:,:,2] )


#  import matplotlib.pyplot as pl
#  pl.ion()
#  import apy_utils as uts
#  from pdb import set_trace as stop
#  for it_sgm in range(36):
#    uts.plot(d_d[it_sgm,:])
#    from pdb import set_trace as stop
#    stop()


  d_s = d_d.real * 0. + 1.e4

  material = materials(n_s, n_in, n_d, n_c, n_ox, d_ox, d_c, d_s, d_d, lamb)

  return material
#
#
#
# ---------------------------------------------------------------------
#

def get_hexagon(num, radius, it_xpos, it_ypos):
#
  import numpy as np
#
  marc = radius * 1.0e0
  x = np.linspace(it_xpos - marc, it_xpos + marc, num)
  ones = x * 0. + 1.
  _, xx = np.meshgrid(ones, x)
  y = np.linspace(it_ypos - marc, it_ypos + marc, num)
  yy, _ = np.meshgrid(y, ones)
#
  a = 30.e0 / 1.8e2 * np.pi
  costhU = it_ypos + np.cos(a) * radius
  costhL = it_ypos - np.cos(a) * radius
  hex_1 = xx * 0.
  ww = np.where( (yy < costhU) & (yy > costhL) )
  hex_1[ww[0], ww[1]] = 1.
#
  a2 = 60.e0 / 1.8e2 * np.pi
  ca = np.cos(a2)
  sa = np.sin(a2)
  rot = np.array([ca, -sa, sa, ca]).reshape(2,2)
#
  out = np.dot(rot, np.vstack([xx.reshape(num*num),yy.reshape(num*num)]))
  xxn = out[0,:].reshape(num,num)
  yyn = out[1,:].reshape(num,num)
  it_nxpos, it_nypos = np.dot(rot, np.array([it_xpos, it_ypos]))
#
  costhU = it_nypos + np.cos(a) * radius
  costhL = it_nypos - np.cos(a) * radius
  hex_2 = xxn * 0.
  ww = np.where( (yyn < costhU) & (yyn > costhL) )
  hex_2[ww[0], ww[1]] = 1.
#
  rot = np.array([ca, sa, -sa, ca]).reshape(2,2)
#
  out = np.dot(rot, np.vstack([xx.reshape(num*num),yy.reshape(num*num)]))
  xxn = out[0,:].reshape(num,num)
  yyn = out[1,:].reshape(num,num)
  it_nxpos, it_nypos = np.dot(rot, np.array([it_xpos, it_ypos]))
#
  costhU = it_nypos + np.cos(a) * radius
  costhL = it_nypos - np.cos(a) * radius
  hex_3 = xxn * 0.
  ww = np.where( (yyn < costhU) & (yyn > costhL) )
  hex_3[ww[0], ww[1]] = 1.
#
  hexag = hex_1 * hex_2 * hex_3
#
  return hexag, xx, yy
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
    exit()

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

  factor=1.1
  # Polar:
  if (layout=='polar'):
    # set rays to be considered in a polar layout:
    azi_val = 24
    rad_to_azi_rat = 24
    #azi_val = 360
    #irads = np.linspace(0., telescope.dim_x * factor, (telescope.num//2)*2)#[1:]
    #iangs = (np.linspace(0., 360., (telescope.num//azi_val)*azi_val+1)) / 180. * np.pi#[0:-1]
    #irads = np.linspace(0., telescope.dim_x * factor, telescope.num+1)#[1:]
    #iangs = np.arange(telescope.num+1) * (2. * np.pi / np.float(telescope.num))#np.linspace(0., 360., telescope.num) / 180. * np.pi#[0:-1]
    ###############irads = np.linspace(0., telescope.dim_x * factor, (telescope.num//2)*2)#[1:]
    irads = np.linspace(0., telescope.dim_x * factor, np.max([3,telescope.num//rad_to_azi_rat]), dtype=np.float64)#[1:]
    irads = np.linspace(0., telescope.dim_x * factor, (telescope.num//2)*2+1, dtype=np.float64)
    drads = irads[1:]-irads[0:-1]
    irads = irads[1:] - drads / 2.

    val = 14.#0.#14.#.0
####    iangs = np.arctan(drads.mean()/irads)*180./np.pi
####
####    from pdb import set_trace as stop
####    stop()

    iangs = (np.linspace(0.+val, 360.+val, np.max([(rad_to_azi_rat*telescope.num)//azi_val,1])*azi_val+1, dtype=np.float64)) / 180. * np.pi#[0:-1]
    iangs = np.linspace(0.+val, 360.+val, (telescope.num//2)*2+1, dtype=np.float64) / 180. * np.pi
    dangs = iangs[1:]-iangs[0:-1]
    iangs = iangs[1:] - dangs / 2.

    print(np.min(iangs))
    ###########from pdb import set_trace as stop
    #iangs = iangs - np.min(iangs)
    #iangs = iangs + val
    #iangs = iangs - np.random.randn() * np.pi
    ###########stop()
######    tpi = 2. * np.pi
######    iangs = (iangs + 10 * tpi) % (tpi)
######    print(np.min(iangs))
######    print(np.min(iangs))
######    print(np.min(iangs))
######    print(np.min(iangs))
######    print(np.min(iangs))
######    print(np.min(iangs))
######    print(np.min(iangs))
######    print(np.min(iangs))
######    print(np.min(iangs))

    ################ww = np.where( (irads<1.3) | (irads>4.5) )
    ################irads = irads[ww]

    # Transform to x and y coordinates:
    
    #ixs = irads[:,None] * np.sin(iangs)[None,:]
    #iys = irads[:,None] * np.cos(iangs)[None,:]
    ixs = irads[:,None] * np.cos(iangs)[None,:]
    iys = irads[:,None] * np.sin(iangs)[None,:]

    d1s = drads.mean()
    d2s = dangs.mean()

  # Cartesian:
  elif (layout=='cartesian'):
    #ixs = np.linspace(-telescope.dim_x*factor, telescope.dim_x*factor, (telescope.num//2)*2)[1:]
    #iys = np.linspace(-telescope.dim_y*factor, telescope.dim_y*factor, (telescope.num//2)*2)[1:]
    ixs = np.linspace(-telescope.dim_x*factor, telescope.dim_x*factor, (telescope.num//2)*2+1, dtype=np.float64)#[1:]
    iys = np.linspace(-telescope.dim_y*factor, telescope.dim_y*factor, (telescope.num//2)*2+1, dtype=np.float64)#[1:]
    ixs = np.linspace(-telescope.dim_x*factor, telescope.dim_x*factor, telescope.num, dtype=np.float64)#[1:]
    iys = np.linspace(-telescope.dim_y*factor, telescope.dim_y*factor, telescope.num, dtype=np.float64)#[1:]
    #ixs = np.linspace(-telescope.dim_x*factor, telescope.dim_x*factor, (telescope.num//2)*2, dtype=np.float64)#[1:]
    #iys = np.linspace(-telescope.dim_y*factor, telescope.dim_y*factor, (telescope.num//2)*2, dtype=np.float64)#[1:]
    #ixs = np.linspace(-telescope.dim_x*factor, telescope.dim_x*factor, telescope.num+1, dtype=np.float64)[1:]
    #iys = np.linspace(-telescope.dim_y*factor, telescope.dim_y*factor, telescope.num+1, dtype=np.float64)[1:]
  
    d1s = np.mean(ixs[1:]-ixs[0:-1])
    d2s = np.mean(iys[1:]-iys[0:-1])

    ix1s = ixs[1:]-d1s/2
    iy1s = iys[1:]-d2s/2
    #ixs = ixs - np.mean(ixs[1:]-ixs[0:-1])/2.
    #iys = iys - np.mean(iys[1:]-iys[0:-1])/2.

    print(ixs.shape)

    #ixs = ixs[:,None] * np.ones(telescope.num, dtype=np.float64)[None,:]
    #iys = iys[None,:] * np.ones(telescope.num, dtype=np.float64)[:,None]
    ixs = ix1s[:,None] * np.ones(iy1s.size, dtype=np.float64)[None,:]
    iys = iy1s[None,:] * np.ones(ix1s.size, dtype=np.float64)[:,None]

    print(ixs.shape)

  else:
    print('Get geometry:')
    print('\tUnknown layout.')
    print('')
    stop()
#
  x1d = ixs.flatten()
  y1d = iys.flatten()

###  import matplotlib.pyplot as pl
###  pl.ion()
###
###  pl.scatter(ixs, iys, edgecolor='none', s=1)
###  pl.gca().set_aspect('equal')
###
###  from pdb import set_trace as stop
###  stop()

  return x1d, y1d, d1s, d2s


def get_geometry(telescope, beam, secondary, segments, osecondary=False, layout='polar', check_rays=False):
#
    import numpy as np

    x1d, y1d, d1s, d2s = get_primary_rays(layout, telescope)
#
    segments.d11d = d1s * 1.
    segments.d21d = d2s * 1.
#
    if (osecondary==True):
      import matplotlib.pyplot as pl
      pl.ion()
      pl.figure(1)
      pl.clf()
      fg,ax=pl.subplots(ncols=1, nrows=1, num=1)
      ax.contour(secondary.mirror.T, [0.5], extent=[np.min(secondary.x) \
          , np.max(secondary.x), np.max(secondary.y), np.min(secondary.y)],colors='r')

    if (check_rays):
      from pdb import set_trace as stop
      import matplotlib.pyplot as pl
      pl.ion()

      pl.figure(1)
      pl.clf()
      fg,ax=pl.subplots(ncols=1, nrows=1, num=1)
      ax.scatter(x1d,y1d,color='k', s=0.1, edgecolor=None)
      stop()

# First of all, we remove all the rays that come from behind the secondary:
    isin=secondary.f2d.ev(x1d, y1d)
    ww = np.where(isin < 0.5)[0]
    if (len(ww)>0):
      x1d = x1d[ww]
      y1d = y1d[ww]

    if (check_rays):
      ax.scatter(x1d,y1d,color='r', s=0.1, edgecolor=None)
      pl.figure(2)
      pl.clf()
      fg,ax=pl.subplots(ncols=1, nrows=1, num=2)
      ax.scatter(x1d,y1d,color='k', s=0.1, edgecolor=None)
      stop()

#
    def check_primary_id(iteles, ixpos, iypos, igxs, igys):
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

        ww = np.where(dumrad >= iteles.rad_min)[0]

        dumrad = dumrad[ww]
        resgx = resgx[ww]
        resgy = resgy[ww]

# Otherwise:
      else:
        print('Primary shape unknown!')
        print('')
        print('')
        exit()
#
      return resgx + ixpos, resgy + iypos

#
    ntot=0


    if (check_rays):
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

####    dum = np.arctan2(y1d, x1d) * 180. / np.pi
####    from pdb import set_trace as stop
####    stop()
    for it_sgm in range(segments.N):

      x, y = check_primary_id(telescope\
         , segments.xpos[it_sgm], segments.ypos[it_sgm], x1d, y1d)

      if (check_rays):
        ax1.scatter(x, y, s=1, edgecolor='none')

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
      dir_cos_i = get_direction_cosines(beam.l_inc,beam.m_inc,beam.n_inc)[:,None]*np.ones(x1.size, dtype=np.float64)[None,:]
      dir_cos_n1, z_1 = new2_get_gcurv(x1,k=telescope.k_1,r=telescope.rc_1,norm=True)

      i_1 = np.arccos(np.sum(dir_cos_n1*dir_cos_i,axis=0))

      dir_cos_o1 = rotate_2d(dir_cos_i, -2 * i_1 * np.sign(x1))

      x2, z_2, dir_cos_n2 = get_intersection(dir_cos_o1, x1, z_1, telescope.nrc_2, telescope.k_2, telescope.dmm)

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

      if (check_rays):
        ax2.scatter(x_sec, y_sec, s=1, edgecolor='none')

      isin=secondary.f2d.ev(x_sec, y_sec)
      ww = np.where(isin > 0.5)[0]

#      ax2.scatter(x_sec[ww], y_sec[ww], s=1, edgecolor='none')

      segments.i['%i' % (it_sgm,)]=i_1[ww]
      segments.i2['%i' % (it_sgm,)]=i_2[ww]
      segments.th1['%i' % (it_sgm,)]=th[ww]

      segments.rad['%i' % (it_sgm,)] = np.sqrt(x[ww]**2+y[ww]**2)

      if (layout=='polar'):
        segments.area['%i' % (it_sgm,)] = segments.rad['%i' % (it_sgm,)] * segments.d11d * segments.d21d
      elif (layout=='cartesian'):
        segments.area['%i' % (it_sgm,)] = segments.rad['%i' % (it_sgm,)] * 0. +  segments.d11d * segments.d21d

      if (osecondary==True):
        ax.plot(x_sec[ww], y_sec[ww],marker='.' \
            , markersize=1, linestyle='none', alpha=0.5)
# COLOR, color=(1./np.float(segments.N)*(it_sgm+1),0.5 \
# COLOR            +0.2/np.float(segments.N)*(it_sgm+1),0.5-0.5/np.float(segments.N)*(it_sgm+1))
      ntot=ntot+ww.size
    #
    if (osecondary==True):
      from pdb import set_trace as stop
      import apy_utils as uts
      ax.set_xlim(np.min(secondary.x), np.max(secondary.x))
      ax.set_ylim(np.min(secondary.y), np.max(secondary.y))
      ax.set_aspect('equal')
      stop()

    print('Total number of rays considered: %i' % ntot)

    if (check_rays):
      fg1.show()
      fg2.show()
      stop()

    return
#
#
#
# ---------------------------------------------------------------------
#

def hex_syst(i, thetaref, th_1, i_22, lamb, n_s, n_c, n_ox, n_d, n_in, d_s, d_c, d_ox, d_d):
#
  import numpy as np
  import time as tm
#
  def reflec_matrix_cnew(i,lamb,n_s,n_c,n_ox,n_d,n_in,d_c,d_ox,d_d):
    #
    # New Muller matrix ellements calculation:
    #
    import time as tm
    #
    def snell(in1,ith1,in2):
    
      return np.arcsin(in1*np.sin(ith1)/in2)
    
    def snell_new(in1,ith1,in2):
    
      if (in2.size!=ith1.size):
        return np.arcsin(in1*np.sin(ith1[None,:])/in2[:,None])
      else:
        return np.arcsin(in1*np.sin(ith1)/in2)
    
    def get_ncosth(in1,ith1,in2,ith2,in3,ith3,pol=None):
    
      if (pol==None):
    #TE wave
        ip1 = in1 * np.cos(ith1)
        ip2 = in2 * np.cos(ith2)
        ip3 = in3 * np.cos(ith3)
      else:
    #TM wave
        ip1 = np.cos(ith1) / in1
        ip2 = np.cos(ith2) / in2
        ip3 = np.cos(ith3) / in3
    
    
      return ip1, ip2, ip3

    def get_ncosth_new(inn,ith,pol=None):
    
      if (pol==None):
    #TE wave
        if (len(ith.shape)==2):
        #if (inn.size!=ith.size):
          ipp = inn[:,None] * np.cos(ith)
        else:
          ipp = inn * np.cos(ith)
          ipp = ipp.reshape(-1,1)
# Comprobar!
      else:
    #TM wave
        #if (inn.size!=ith.size):
        if (len(ith.shape)==2):
          ipp = inn[:,None] / np.cos(ith)
        else:
          ipp = np.cos(ith) / inn
          ipp = ipp.reshape(-1,1)
# Comprobar!
      return ipp
    
    def get_matrix_new(ith,inn,ik,ih,lambda0):
    
      nfilms = len(inn)-2
      npts = ith.size

      res = np.zeros((npts, 12), dtype=np.complex64)

      cnt_pol = -1
      for it_pol in [None, True]:

        cnt_pol += 1
    
        char_mat = np.zeros((npts, 2, 2), dtype=np.complex64)
        char_mat[:,0,0]=1.
        char_mat[:,1,1]=1.

        cn = inn+1j*ik
        th = snell_new(inn[0], ith, cn)
        p = get_ncosth_new(cn,th,pol=it_pol)
        beta = 2. * np.pi / lambda0 * cn[:,None] * ih[:,None] * np.cos(th) 
        cb = np.cos(beta)
        sb = np.sin(beta)
        pos = sb / p
        pts = p * sb

        for it_nnn in range(nfilms):
    
          it_mat = np.zeros((2,2,npts), dtype=np.complex64)
          it_mat[0,0,:]=cb[it_nnn+1]
          it_mat[0,1,:]=-np.complex(0.,1.) * pos[it_nnn+1]
          it_mat[1,0,:]=-np.complex(0.,1.) * pts[it_nnn+1]
          it_mat[1,1,:]=cb[it_nnn+1]
          for it_ith in range(npts):
            #res = np.dot(res, it_mat)
            #char_mat[it_ith,:,:] = char_mat[it_ith,:,:].dot(it_mat)
            char_mat[it_ith,:,:] = np.dot(char_mat[it_ith,:,:],it_mat[:,:,it_ith])
        # Reflection:
        tr1 = (char_mat[:,0,0] + char_mat[:,0,1] * p[-1,:]) * p[0,:]
        tr2 = char_mat[:,1,0] + char_mat[:,1,1] * p[-1,:]
        r = (tr1 - tr2) / (tr1 + tr2)
        res[:,cnt_pol+4]=np.abs(r)
        res[:,cnt_pol+6]=np.arctan2(r.imag,r.real)
        res[:,cnt_pol]=res[:,cnt_pol+4]**2
        del(r)
#>< NOT USED HERE!        #rho.append(np.abs(r))
#>< NOT USED HERE!        #phi.append(np.arctan(r.imag/r.real))
#>< NOT USED HERE!        #rr.append(np.abs(r)**2)
#>< NOT USED HERE!        # Transmission:
#>< NOT USED HERE!        tr1 = 2. * p1
#>< NOT USED HERE!        tr2a = (char_mat[:,0,0] + char_mat[:,0,1] * pl) * p1
#>< NOT USED HERE!        tr2b = char_mat[:,1,0] + char_mat[:,1,1] * pl
#>< NOT USED HERE!    
#>< NOT USED HERE!        t = tr1 / (tr2a + tr2b)
#>< NOT USED HERE!        #tau.append(np.abs(t))
#>< NOT USED HERE!        #chi.append(np.arctan(t.imag/t.real))
#>< NOT USED HERE!        #tt.append(pl/p1*np.abs(t)**2)
#>< NOT USED HERE!        res[:,cnt_pol+8]=np.abs(t)
#>< NOT USED HERE!        res[:,cnt_pol+10]=np.arctan(t.imag/t.real)
#>< NOT USED HERE!        res[:,cnt_pol+2]=pl/p1*res[:,cnt_pol+8]**2
#>< NOT USED HERE!        del(t)
    
      return res

# Last revision:
    nc = [  n_in, n_ox, n_c, n_s]
    n = np.array(nc).real
    k = np.array(nc).imag
    h = np.array([np.inf, d_ox, d_c,d_s])

    check = get_matrix_new(i,n,k,h,lamb)
    rp = np.abs(check[:,5])**2
    x = check[:,4].real / check[:,5].real
    tau = check[:,6].real-check[:,7].real

    return rp, x, tau

# End new calculation
#
  def util(i, n, k):
    f = 0.5e0 * (n**2 - k**2 - np.sin(i)**2 + \
        np.sqrt((n**2 - k**2 - np.sin(i)**2)**2 + 4.e0 * n**2 * k**2))
    g = 0.5e0 * (k**2 - n**2 + np.sin(i)**2 + \
        np.sqrt((n**2 - k**2 - np.sin(i)**2)**2 + 4.e0 * n**2 * k**2))
    #
    r_p1 = f + g - 2.e0 * np.sqrt(f) * np.cos(i) + np.cos(i)**2
    r_p2 = f + g + 2.e0 * np.sqrt(f) * np.cos(i) + np.cos(i)**2
    r_p = r_p1 / r_p2
    #
    s_tau = 2.e0 * np.sqrt(g) * np.sin(i) * np.tan(i)
    c_tau = np.sin(i)**2 * np.tan(i)**2 - (f + g)
    tau = np.arctan(s_tau / c_tau)
    cos_tau = np.cos(tau)
    sin_tau = np.sin(tau)
    #
    x_1 = f + g - 2.e0 * np.sqrt(f) * np.sin(i) * np.tan(i) + np.sin(i)**2 * np.tan(i)**2
    x_2 = f + g + 2.e0 * np.sqrt(f) * np.sin(i) * np.tan(i) + np.sin(i)**2 * np.tan(i)**2
    x_s = x_1 / x_2
    x = np.sqrt(x_s)
    #
    return r_p, x, tau
#
  def get_rot(iang,dims):

    res = np.zeros((dims, 4, 4), dtype=np.float64)
    res[:,0,0] = 1.
    res[:,3,3] = 1.

    c = np.cos(iang*2.0)
    s = np.sin(iang*2.0)
    #print(c[0], s[0])

#    cang = np.zeros(dims, dtype=np.complex64)
#    cang.imag = iang * 2.
#    test = (np.exp(cang)+np.exp(-cang))/2.
#    c = test.real
#    test = (np.exp(cang)-np.exp(-cang))/2.
#    s = test.imag
    #print(c[0], s[0])
    #from pdb import set_trace as stop
    #stop()
    res[:,1,2] = -s * 1.
    res[:,2,1] = s * 1.
    res[:,1,1] = c * 1.
    res[:,2,2] = c * 1.

    return res

  def get_ref_mat(r,x,t,sz):

    mat = np.zeros((sz, 4,4), dtype=np.float64)

    mat[:,0,0]=x**2+1.
    mat[:,1,1]=mat[:,0,0]*1.
    mat[:,0,1]=x**2-1.
    mat[:,1,0]=mat[:,0,1]*1.
    
    mat[:,2,2]=2.*x*np.cos(t)
    mat[:,3,3]=mat[:,2,2]*1.
    mat[:,2,3]=2.*x*np.sin(t)
    mat[:,3,2]=-mat[:,2,3]*1.
    mat = mat / 2. * r[:,None,None]
   
    return mat
 
#
#
#
######## testing:
#######  print(np.size(i))
#######  th_1=np.linspace(0.,360.,361)/180.*np.pi
#######  i = i[0:np.size(th_1)]
#######  i_22 = i_22[0:np.size(th_1)]
######## testing.
#
#TO ACTIVATE ONCE TESTS ARE DONE!!!
#TM  print('__________')
#TM  t0=tm.time()
  #
  nrays = np.size(i)
  #
  r_p1,X1,tau1 = reflec_matrix_cnew(i,lamb,n_s,n_c,n_ox,n_d,n_in,d_c,d_ox,d_d)
  r_p2, X2, tau2 = util(i_22, n_c.real, n_c.imag)

#END TO ACTIVATE
##########???:  theta1 = th_1 * 1
##########???:  theta2 = np.ones(np.size(i)) * np.pi
##########???:  theta3 = np.pi - th_1

  #theta1 = - (np.pi / 2. + th_1)
  #theta1 = np.pi / 2. - th_1
  theta1 = - (np.pi - (th_1-thetaref))
  theta2 = np.ones(nrays, dtype=np.float64) * np.pi
  #theta3 = - (np.pi / 2. - th_1)
  theta3 = theta2 * 1.
  theta4 = -theta1


#END TO ACTIVATE
#
  #
  res = np.zeros((nrays,4,4), dtype=np.float64)
  
  rot1 = get_rot(theta1, nrays)
  rot2 = get_rot(theta2, nrays)
  rot3 = get_rot(theta3, nrays)
  rot4 = get_rot(theta4, nrays)

#  from pdb import set_trace as stop
#  stop()

  ref1 = get_ref_mat(r_p1, X1, tau1, nrays)
  ref2 = get_ref_mat(r_p2, X2, tau2, nrays)
 
  for it_nnn in range(nrays):
  
    res[it_nnn,:,:] = np.dot(rot4[it_nnn,:,:], np.dot(rot3[it_nnn,:,:]\
        , np.dot(ref2[it_nnn,:,:]\
        , np.dot(rot2[it_nnn,:,:]\
        , np.dot(ref1[it_nnn,:,:], rot1[it_nnn,:,:])))))

###    res[it_nnn,:,:] = np.dot(rot2[it_nnn,:,:]\
###        , np.dot(ref1[it_nnn,:,:], rot1[it_nnn,:,:]))

#ref1:    res[it_nnn,:,:] = \
#ref1:        np.dot(rot3[it_nnn,:,:], np.dot(rot2[it_nnn,:,:], np.dot(ref1[it_nnn,:,:], rot1[it_nnn,:,:])))
#ref2:    res[it_nnn,:,:] = \
#ref2:        np.dot(rot4[it_nnn,:,:], np.dot(rot3[it_nnn,:,:], np.dot(ref2[it_nnn,:,:], np.dot(rot2[it_nnn,:,:], rot1[it_nnn,:,:]))))
###    res[it_nnn,:,:] = \
###        np.dot(ref2[it_nnn,:,:], ref1[it_nnn,:,:])
######    res[it_nnn,:,:] = \
######        np.dot(ref2[it_nnn,:,:], np.dot(rot2[it_nnn,:,:], ref1[it_nnn,:,:]))
####    res[it_nnn,:,:] = \
####        rot2[it_nnn,:,:] * 1.
##    res[it_nnn,:,:] = np.dot(rot4[it_nnn,:,:], np.dot(rot3[it_nnn,:,:]\
##        , np.dot(rot2[it_nnn,:,:], rot1[it_nnn,:,:])))

######## testing:
#######  print(res[0,:,:])
##  from pdb import set_trace as stop
##  stop()

#>< IIR? :  el_11 = np.sum(((1.e0 + X1_2)*(1.e0 + X2_2) + (-1.e0 + X1_2)*(-1.e0 + X2_2) * ctheta2)*rp)/auxx
#>< IIR? :  el_12 = np.sum((ctheta1*((-1.e0+X1_2)*(1.e0+X2_2)+(1.e0+X1_2)*(-1.e0+X2_2)*ctheta2)-2.e0*X1*\
#>< IIR? :      (-1.e0+X2_2)*ctau1*stheta1*stheta2)*rp)/auxx
#>< IIR? :  el_13 = np.sum((((-1.e0+X1_2)*(1.e0+X2_2)+(1.e0+X1_2)*(-1.e0+X2_2)*ctheta2)*stheta1+2.e0*X1*\
#>< IIR? :      (-1.e0+X2_2)*ctau1*ctheta1*stheta2)*rp)/auxx
#>< IIR? :  el_14 = np.sum((2.e0*X1*(-1.e0+X2_2)*stau1*stheta2)*rp)/auxx
#>< IIR? :  el_21 = np.sum(((1.e0+X1_2)*(-1.e0+X2_2)*ctheta3+(-1.e0+X1_2)*((1.e0+X2_2)*ctheta2*ctheta3-\
#>< IIR? :      2.e0*X2*ctau2*stheta2*stheta3))*rp)/auxx
#>< IIR? :  el_22 = np.sum((-stheta1*(-4.e0*X1*X2*stau1*stau2*stheta3+2.e0*X1*ctau1*((1.e0+X2_2)*\
#>< IIR? :      ctheta3*stheta2+2.e0*X2*ctau2*ctheta2*stheta3))+ctheta1*((-1.e0+X1_2)*\
#>< IIR? :      (-1.e0+X2_2)*ctheta3+(1.e0+X1_2)*((1.e0+X2_2)*ctheta2*ctheta3-2.e0*X2*ctau2*\
#>< IIR? :      stheta2*stheta3)))*rp)/auxx
#>< IIR? :  el_23 = np.sum((ctheta1*(-4.e0*X1*X2*stau1*stau2*stheta3+2.e0*X1*ctau1*((1.e0+X2_2)*\
#>< IIR? :      ctheta3*stheta2+2.e0*X2*ctau2*ctheta2*stheta3))+stheta1*((-1.e0+X1_2)*\
#>< IIR? :      (-1.e0+X2_2)*ctheta3+(1.e0+X1_2)*((1.e0+X2_2)*ctheta2*ctheta3-2.e0*X2*ctau2*\
#>< IIR? :      stheta2*stheta3)))*rp)/auxx
#>< IIR? :  el_24 = np.sum((2.e0*X1*(1.e0+X2**2)*ctheta3*stau1*stheta2+4.e0*X1*X2*(ctau2*ctheta2*stau1+\
#>< IIR? :      ctau1*stau2)*stheta3)*rp)/auxx
#>< IIR? :  el_31 = np.sum((-(1.e0+X1_2)*(-1.e0+X2_2)*stheta3+(-1.e0+X1_2)*(-2.e0*X2*ctau2*ctheta3*stheta2-\
#>< IIR? :      (1.e0+X2_2)*ctheta2*stheta3))*rp)/auxx
#>< IIR? :  el_32 = np.sum((-2.e0*X2*ctheta3*(-2.e0*X1*stau1*stau2*stheta1+(1.e0+X1_2)*ctau2*ctheta1*\
#>< IIR? :      stheta2)-ctheta1*((-1.e0+X1_2)*(-1.e0+X2_2)+(1.e0+X1_2)*(1.e0+X2_2)*ctheta2)*\
#>< IIR? :      stheta3+2.e0*X1*ctau1*stheta1*(-2.e0*X2*ctau2*ctheta2*ctheta3+\
#>< IIR? :      (1.e0+X2_2)*stheta2*stheta3))*rp)/auxx
#>< IIR? :  el_33 = np.sum((stheta1*(-(-1.e0+X1_2)*(-1.e0+X2_2)*stheta3+(1.e0+X1_2)*(-2.e0*X2*ctau2*ctheta3*\
#>< IIR? :      stheta2-(1.e0+ X2_2)*ctheta2*stheta3))+ctheta1*(-4.e0*X1*X2*ctheta3*\
#>< IIR? :      stau1*stau2-2.e0*X1*ctau1*(-2.e0*X2*ctau2*ctheta2*ctheta3+(1.e0+X2_2)*\
#>< IIR? :      stheta2*stheta3)))*rp)/auxx
#>< IIR? :  el_34 = np.sum((4.e0*X1*X2*ctau1*ctheta3*stau2+2.e0*X1*stau1*(2.e0*X2*ctau2*ctheta2*ctheta3-\
#>< IIR? :      (1.e0+X2_2)*stheta2*stheta3))*rp)/auxx
#>< IIR? :  el_41 = np.sum((2.e0*(-1.e0+X1_2)*X2*stau2*stheta2)*rp)/auxx
#>< IIR? :  el_42 = np.sum((4.e0*X1*X2*(ctau2*stau1+ctau1*ctheta2*stau2)*stheta1+2.e0*(1.e0+X1_2)*X2*\
#>< IIR? :      ctheta1*stau2*stheta2)*rp)/auxx
#>< IIR? :  el_43 = np.sum((-4.e0*X1*X2*ctheta1*(ctau2*stau1+ctau1*ctheta2*stau2)+2.e0*(1.e0+X1_2)*X2*\
#>< IIR? :      stau2*stheta1*stheta2)*rp)/auxx
#>< IIR? :  el_44 = np.sum((4.e0*X1*X2*(ctau1*ctau2-ctheta2*stau1*stau2))*rp)/auxx
#>< IIR? :
#>< IIR? :  stop()
#>< IIR? :
#>< IIR? :  if (el_33<0.):
#>< IIR? :    print el_33
#>< IIR? :    stop()
#>< IIR? :  if (el_22<0.):
#>< IIR? :    print el_22
#>< IIR? :    stop()
#>< IIR? :  if (el_44<0.):
#>< IIR? :    print el_44
#>< IIR? :    stop()
#>< IIR? :  if (el_11<0.):
#>< IIR? :    print el_11
#>< IIR? :    stop()

#TM  t6=tm.time()
#TM  print((t1-t0)/(t6-t0))
#TM  print((t2-t1)/(t6-t0))
#TM  print((t3-t2)/(t6-t0))
#TM  print((t4-t3)/(t6-t0))
#TM  print((t5-t4)/(t6-t0))
#TM  print((t6-t5)/(t6-t0))
#TM  stop()
#
#>< IIR? :  matr = np.array([[el_11,el_12,el_13,el_14],[el_21,el_22,el_23,el_24],[el_31,el_32,el_33,el_34],[el_41,el_42,el_43,el_44]])
#
#>< IIR? :  return matr
  return res
#
#
# ---------------------------------------------------------------------
#

def get_mdust(n, mean_dust):
#
  import numpy as np
#
  def snell(n1, n2, theta1):
  
    return np.arcsin(n1 * np.sin(theta1) / n2)
#
  def tpar(n1, n2, theta1, theta2):
  
    return (2. * n1 * np.cos(theta1)) / (n1 * np.cos(theta1) + n2 * np.cos(theta2))
#
  def tper(n1, n2, theta1, theta2):
  
    return (2. * n1 * np.cos(theta2)) / (n1 * np.cos(theta2) + n2 * np.cos(theta1))
#
  #t_par_1 = 
  #t_par_2 = 
  #t_per_1 = 
  #t_per_2 = 
#
  calculate = 0 # This calculation has to be repeated when changing n of dust
  if (calculate == 1):
    ####upplim = 40.
    rtd = 180. / np.pi
    upplim = 90.#np.arcsin(np.sin(np.pi/2.)*n2.real/n1.real)*rtd
#
    n1 = 1.
    n2 = n * 1.
    n3 = 1.
#
    num = 100001
    dtor = np.pi / 180.

    thetas = np.linspace(0., upplim, num)
    dtheta = (thetas[1] - thetas[0]) * dtor

    int1 = 0.
    int3 = 0.

    for it_th1 in thetas:

      theta1 = it_th1 * dtor
      theta2 = snell(n1, n2, theta1)
  ####    print(it_th1, theta2/dtor)
  ####    if ( theta2/dtor > 45.):
  ####      print(theta2/dtor)
  ####      continue
      theta3 = snell(n2.real, n3.real, theta2)

      tpar1 = tpar(n1.real, n2.real, theta1, theta2)
      tpar2 = tpar(n2.real, n3.real, theta2, theta3)
      tper1 = tper(n1.real, n2.real, theta1, theta2)
      tper2 = tper(n2.real, n3.real, theta2, theta3)
####      theta3 = snell(n2, n3, theta2)
####
####      tpar1 = tpar(n1, n2, theta1, theta2)
####      tpar2 = tpar(n2, n3, theta2, theta3)
####      tper1 = tper(n1, n2, theta1, theta2)
####      tper2 = tper(n2, n3, theta2, theta3)

      it1 = (tpar1*tpar2)**2+(tper1*tper2)**2
      it3 = 2. * tpar1 * tpar2 * tper1 * tper2

      int1 += it1 * np.sin(theta1) * dtheta
      int3 += it3 * np.sin(theta1) * dtheta


    print("a = %.8f / 2. / np.pi" % (int1,))
    print("d = %.8f / 2. / np.pi" % (int3,))

###    print(int1)
###    print(int3)
    from pdb import set_trace as stop
    stop()

  #a = (t_par_2 * t_par_1)**2 + (t_per_2 * t_per_1)**2
  #d = 2. * (np.conj(t_par_2*t_par_1)*(t_per_2*t_per_1)).real
  #a = 1. + ((4. * n**2)/((1.+n**2)**2))**2
  #d = (8. * n**2)/((1.+n**2)**2)
  #a = 0.42751143
  #d = 0.42720953
  a = 1.38942494 / 2. / np.pi
  d = 1.35253306 / 2. / np.pi
#
  #a = a.real[0] / 4. * np.pi
  #d = d.real[0] / 4. * np.pi
  a = a / 4. * np.pi
  d = d / 4. * np.pi
#
  mat = np.zeros((4,4))
  mat[0,0] = np.exp(mean_dust.real * (2. * a - 1.)) 
  mat[1,1] = np.exp(mean_dust.real * ((a + d) - 1.)) 
  mat[2,2] = np.exp(mean_dust.real * ((a + d) - 1.)) 
  mat[3,3] = np.exp(mean_dust.real * (2. * d - 1.)) 
#
  return mat

#
# ---------------------------------------------------------------------
#

def new_get_mdust(p0, mean_dust):
#
  import numpy as np
#
  mat = np.zeros((4,4))
  mat[0,0] = 1. - mean_dust.real# * (1.-p0)
  mat[1,1] = 1. - mean_dust.real# * (1.-p0)
  mat[2,2] = 1. - mean_dust.real# * (1.-p0)
  mat[3,3] = 1. - mean_dust.real# * (1.-p0)
#
  return mat
#
#
# ---------------------------------------------------------------------
#

def get_mueller_time(segments, materials,cdust=False, thetaref=0.):

  import numpy as np
  import time as tm
  from pdb import set_trace as stop
  import apy_utils as uts
#
  tarea = 0.
  for it in range(segments.N):
    tarea += np.sum(segments.area['%i' % (it,)])

#RUN TIMES:
  res = np.zeros((segments.SimulatedTimes, segments.N, 4, 4))
  avg_res = np.zeros((segments.SimulatedTimes, 4, 4))
#  avg_res_test = np.zeros((4, 4))
  for it_tms in range(segments.SimulatedTimes):
    print("Timestep: %i of %i" % (it_tms, segments.SimulatedTimes))
#RUN SEGMENTS:
    npts_seg = np.ones((segments.N,), dtype=np.float64)
    ntot2=0

#
# Testing:
#    pl.close(1)
#    pl.close(2)
#    pl.close(3)
#    pl.close(4)
#    #fg, ax = pl.subplots(num=1, nrows=6, ncols=6)
#    #ax=ax.flatten()
#    fg1, ax1 = pl.subplots(num=1, nrows=1, ncols=1)
#    fg2, ax2 = pl.subplots(num=2, nrows=1, ncols=1)
#    fg3, ax3 = pl.subplots(num=3, nrows=1, ncols=1)
#    fg4, ax4 = pl.subplots(num=4, nrows=1, ncols=1)
# Testing.
#


    for it_sgm in range(segments.N):
#
      if (segments.i['%i' % (it_sgm,)].size<1):
        npts_seg[it_sgm] = 0.
        res[it_tms,it_sgm,:,:] = res[it_tms,it_sgm,:,:] * 0.
        continue

      mmpa_full = hex_syst(segments.i['%i' % (it_sgm,)], thetaref, \
          segments.th1['%i' % (it_sgm,)], \
          segments.i2['%i' % (it_sgm,)], \
          materials.lamb, materials.n_s, materials.n_c\
          , materials.n_ox, materials.n_d, materials.n_in, \
          materials.d_s[it_sgm, it_tms], \
          materials.d_c[it_sgm, it_tms]\
          , materials.d_ox[it_sgm, it_tms], materials.d_d[it_sgm, it_tms])

      #GET Mdust
      if (cdust==True):
        ###mmpd = get_mdust(materials.n_d, materials.d_d[it_sgm, it_tms])
        mmpd = new_get_mdust(0.1, materials.d_d[it_sgm, it_tms])
        res_full = mmpa_full.dot(mmpd)
      else:
        res_full = mmpa_full * 1.


#
# Testing:

#      avg_res_test += np.sum(res_full * segments.area['%i' % (it_sgm,)][:,None,None], axis=0)
#
#      xx = segments.rad['%i' % (it_sgm,)]*np.cos(segments.th1['%i' % (it_sgm,)])
#      yy = segments.rad['%i' % (it_sgm,)]*np.sin(segments.th1['%i' % (it_sgm,)])
#      #im1 = ax1.scatter(xx, yy, marker='s', s=10, c=res_full[:,1,0] * segments.area['%i' % (it_sgm,)] / tarea, cmap='RdGy', vmin=-0.0000001, vmax=0.0000001)
#      #im2 = ax2.scatter(xx, yy, marker='s', s=10, c=res_full[:,2,0] * segments.area['%i' % (it_sgm,)] / tarea, cmap='RdGy', vmin=-0.0000001, vmax=0.0000001)
#      im1 = ax1.scatter(xx, yy, marker='s', s=10, c=res_full[:,1,2], cmap='RdGy', vmin=-0.0001, vmax=0.0001)
#      im2 = ax2.scatter(xx, yy, marker='s', s=10, c=res_full[:,2,1], cmap='RdGy', vmin=-0.0001, vmax=0.0001)
#      im3 = ax3.scatter(xx, yy, marker='s', s=10, c=res_full[:,1,0], cmap='RdGy', vmin=-0.0002, vmax=0.0002)
#      im4 = ax4.scatter(xx, yy, marker='s', s=10, c=res_full[:,2,0], cmap='RdGy', vmin=-0.0002, vmax=0.0002)
#      #ax[it_sgm].hist(res_full[:,1,0], bins=101)
# Testing.
#

      nn,_,_ = res_full.shape
      ntot2=ntot2+nn
########      for itn in range(nn):
########        res[it_tms,it_sgm,:,:]=res[it_tms,it_sgm,:,:] \
########            + res_full[itn,:,:] * segments.rad['%i' % (it_sgm,)][itn]
########      res[it_tms,it_sgm,:,:] = res[it_tms,it_sgm,:,:] \
########          * segments.area['%i' % (it_sgm,)]
########      #
########      # Area to normalize:
########      npts_seg[it_sgm] = np.sum(segments.rad['%i' % (it_sgm,)] \
########          * segments.area['%i' % (it_sgm,)])

  ####    for itn in range(nn):
  ####      res[it_tms,it_sgm,:,:]=res[it_tms,it_sgm,:,:] \
  ####          + res_full[itn,:,:] * segments.area['%i' % (it_sgm,)][itn]
      res[it_tms,it_sgm,:,:]=np.sum(res_full[:,:,:] * segments.area['%i' % (it_sgm,)][:,None,None], axis=0) / tarea
      #
      #
      # Area to normalize:
      #######npts_seg[it_sgm] = np.sum(segments.area['%i' % (it_sgm,)])


    print('Total number of rays considered: %i' % ntot2)
    print('Total area considered: %.4f' % (tarea,))
    

    # Since we consider all the points of each segment, the average must...
    # ... be done accordingly:
    avg_res[it_tms,:,:] = np.sum(res[it_tms,:,:,:], axis=0)#### \
#####        / np.sum(npts_seg)

#
#  cbar1 = pl.colorbar(im1, ax=ax1)
#  cbar2 = pl.colorbar(im2, ax=ax2)
#  cbar3 = pl.colorbar(im3, ax=ax3)
#  cbar4 = pl.colorbar(im4, ax=ax4)
#  for ax in [ax1, ax2, ax3, ax4]:
#    ax.set_aspect('equal')
#    ax.grid()
#  pl.draw()


#  stop()

  return avg_res, res, npts_seg


#
# ---------------------------------------------------------------------
#


def plot_primary(telescope, beam, secondary, segments \
    , pax=None, prays=False, fcol='', layout='polar'):
#
  import numpy as np
  import PyColors as PyC

  def plot_primary_segment(iteles, ixpos, iypos, i2ax, ifcol):
#
    if (iteles.primary_shape == 'hexagonal'):

      dy = np.sqrt(iteles.radius**2-(iteles.radius/2.)**2)

      p1x = ixpos + iteles.radius
      p1y = iypos + 0.

      p2x = ixpos + iteles.radius / 2.
      p2y = iypos + dy

      p3x = ixpos - iteles.radius / 2.
      p3y = iypos + dy

      p4x = ixpos - iteles.radius
      p4y = iypos + 0.

      p5x = ixpos - iteles.radius / 2.
      p5y = iypos - dy

      p6x = ixpos + iteles.radius / 2.
      p6y = iypos - dy

      ix = [ p1x, p2x, p3x, p4x, p5x, p6x, p1x ]
      iy = [ p1y, p2y, p3y, p4y, p5y, p6y, p1y ]

      if (len(ifcol)==0):

        im = i2ax.plot(ix,iy,color='k',lw=1.)

      else:

        from pdb import set_trace as stop
        #stop()

        ix1 = [p4x, p3x, p2x, p1x]
        iy1 = [p4y, p3y, p2y, p1y]
        iy2 = [p4y, p5y, p6y, p1y]

        im = i2ax.fill_between(ix1, iy1, iy2, color=ifcol[0])
        im = i2ax.plot(ix,iy,color='k',lw=1.)

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

      ww = np.where(dumrad >= iteles.rad_min)[0]

      dumrad = dumrad[ww]
      resgx = resgx[ww]
      resgy = resgy[ww]

      im = []

# Otherwise:
    else:
      print('Primary shape unknown!')
      print('')
      print('')
      exit()
#
    return im


  #
  # Main: plot_primary
  #
  x1d, y1d, d1s, d2s = get_primary_rays(layout, telescope)


  if (np.any(pax)==None):
    import matplotlib.pyplot as pl
    pl.ion()
    pl.figure(1)
    pl.clf()
    fg,iax=pl.subplots(ncols=1, nrows=1, num=1)
  else:
    iax=pax


  ntot=0
  from pdb import set_trace as stop
  if (fcol != ''):
    tpfcol=PyC.get_colors(segments.Nfamily, cmap=fcol)
    aux = []
    for itn in range(segments.N):
      aux.append([tpfcol[np.int32(segments.families[itn])-1]])
  # Put color legend in GTC layout:
    for it_fml in range(segments.Nfamily):
      iax.plot([],[],linewidth=5, label=r'Fam.: %i' % (it_fml+1,) \
          , color=tpfcol[it_fml])
    iax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0) \
      , fancybox=True, shadow=True,ncol=2)

  else:
    aux = [[],]*segments.N

  for it_sgm in range(segments.N):
    itim = plot_primary_segment(telescope\
       , segments.xpos[it_sgm], segments.ypos[it_sgm], iax, aux[it_sgm])


  if (prays != False):

    ccol = (1.0,0.,0.)
    clw = 2.
    calp = 0.6

    for rad in irads[1:]:
      var = 2. * np.pi * np.linspace(0.,1.,1001)
      iax.plot(rad*np.sin(var), rad*np.cos(var) \
          , color=ccol, lw=clw, alpha=calp) 
    for ang in iangs[1:]:
      ca = np.cos(ang)
      sa = np.sin(ang)
      ix = [irads[1]*ca, irads[-1]*ca]
      iy = [irads[1]*sa, irads[-1]*sa]
      iax.plot(ix, iy, color=ccol, lw=clw, alpha=calp) 


  iax.set_xlim(-telescope.dim_x, telescope.dim_x)
  iax.set_ylim(-telescope.dim_y, telescope.dim_y)

#
# ---------------------------------------------------------------------
#

def plot_geometry(telescope, beam, secondary, segments):
    #
    #
    #
    # This has been copied from get_geometry!
    # If get_geometry changes I have to change this...
    # ...accordingly
    #
    #
    #
#
    import numpy as np
    layout='polar'
    #layout='cartesian'

    x1d, y1d, d1s, d2s = get_primary_rays(layout, telescope)



# First of all, we remove all the rays that come from behind the secondary:
    isin=secondary.f2d.ev(x1d, y1d)
    ww = np.where(isin < 0.5)[0]
    if (len(ww)>0):
      x1d = x1d[ww]
      y1d = y1d[ww]

#    ax.scatter(x1d,y1d,color='r', s=0.1, edgecolor=None)

#
    segments.d11d = d1s * 1.
    segments.d21d = d2s * 1.
#
    def check_primary_id(iteles, ixpos, iypos, igxs, igys, fill_nan=False):
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
        if (fill_nan==True):
          ww = np.where(np.abs(dumgy) > costhU)[0]
          resgx[ww] = np.nan
          resgy[ww] = np.nan
        else:
          ww = np.where(np.abs(dumgy) <= costhU)[0]
          dumgx = dumgx[ww] * 1.
          dumgy = dumgy[ww] * 1.
          resgx = resgx[ww] * 1.
          resgy = resgy[ww] * 1.
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
        if (fill_nan==True):
          ww = np.where(np.abs(dumyg) > costhU)[0]
          resgx[ww] = np.nan
          resgy[ww] = np.nan
        else:
          ww = np.where( np.abs(dumyg) <= costhU)[0]
          dumgx = dumgx[ww] * 1.
          dumgy = dumgy[ww] * 1.
          resgx = resgx[ww] * 1.
          resgy = resgy[ww] * 1.
  #
        rot = np.array([ca, sa, -sa, ca]).reshape(2,2)
  #
        out = np.dot(rot, np.vstack([dumgx,dumgy]))
        dumxg = out[0,:] * 1.
        dumyg = out[1,:] * 1.
  #
        if (fill_nan==True):
          ww = np.where(np.abs(dumyg) > costhU)[0]
          resgx[ww] = np.nan
          resgy[ww] = np.nan
        else:
          ww = np.where( np.abs(dumyg) <= costhU)[0]
          dumgx = dumgx[ww] * 1.
          dumgy = dumgy[ww] * 1.
          resgx = resgx[ww] * 1.
          resgy = resgy[ww] * 1.
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

        ww = np.where(dumrad >= iteles.rad_min)[0]

        dumrad = dumrad[ww]
        resgx = resgx[ww]
        resgy = resgy[ww]

# Otherwise:
      else:
        print('Primary shape unknown!')
        print('')
        print('')
        exit()
#
      return resgx + ixpos, resgy + iypos

#
# Starting plot_geometry:
#
#
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as pl
    pl.ion()
    from pdb import set_trace as stop
    from scipy.interpolate import bisplrep, bisplev
    from matplotlib import cm

    fignum=2
    pl.close(fignum)
    fg, ax = pl.subplots(num=fignum,subplot_kw={"projection": "3d"})

    if ('GTC' in telescope.ID):
      im1max=12.
      im2max=15
    elif ('ELT' in telescope.ID):
      im1max=20.
      im2max=25.


#
#

##    plt1 = []
##    plt2 = []
##    plt3 = []
    ntot=0
    for it_sgm in range(segments.N):

      x, y = check_primary_id(telescope\
         , segments.xpos[it_sgm], segments.ypos[it_sgm], x1d, y1d)

# Plot:
      num=21
      #f2d=bisplrep(x,y,z_1)
      x1 = np.linspace(x.min(), x.max(), num)
      y1 = np.linspace(y.min(), y.max(), num)
      #z1=bisplev(x1, y1, f2d)

      x1_2d = x1[:,None] * np.ones(num)[None,:]
      y1_2d = y1[None,:] * np.ones(num)[:,None]

      x1_2d = x1_2d.flatten()
      y1_2d = y1_2d.flatten()

      x, y = check_primary_id(telescope\
         , segments.xpos[it_sgm], segments.ypos[it_sgm], x1_2d, y1_2d, fill_nan=True)


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
      dir_cos_i = get_direction_cosines(beam.l_inc,beam.m_inc,beam.n_inc)[:,None]*np.ones(x1.size)[None,:]
      dir_cos_n1, z_1 = new2_get_gcurv(x1,k=telescope.k_1,r=telescope.rc_1,norm=True)

      i_1 = np.arccos(np.sum(dir_cos_n1*dir_cos_i,axis=0))

      dir_cos_o1 = rotate_2d(dir_cos_i, -2 * i_1 * np.sign(x1))

      x2, z_2, dir_cos_n2 = get_intersection(dir_cos_o1, x1, z_1, telescope.nrc_2, telescope.k_2, telescope.dmm)

      i_2 = np.arccos(np.sum(dir_cos_n2*dir_cos_o1,axis=0))

      dir_cos_f = rotate_2d(dir_cos_o1, 2 * i_2 * np.sign(x2))

      zf=-10.
      xf = get_intersection_focus(dir_cos_f, x2, z_2,zf)

      aux = i_1[i_1==i_1] * 180. / np.pi
      print(" Min(i1)=%.2f ; Max(i1)=%.2f ; <i1>=%.2f in deg." % (aux.min(), aux.max(), aux.mean()))
      aux = i_2[i_2==i_2] * 180. / np.pi
      print(" Min(i2)=%.2f ; Max(i2)=%.2f ; <i2>=%.2f in deg." % (aux.min(), aux.max(), aux.mean()))
      del(aux)

      # Following: https://pundit.pratt.duke.edu/wiki/Python:Plotting_Surfaces
      def get_color(a, clim): return np.nan_to_num(a)*180./np.pi/clim

      x1s = x.reshape(num, num)
      y1s = y.reshape(num, num)
      cmm = cm.viridis
      surf1 = ax.plot_surface(x1s, y1s, z_1.reshape(num, num) \
          ,linewidth=0, antialiased=False, facecolors \
          =cmm(get_color(i_1.reshape(num,num), im1max)) \
          , rcount=51, ccount=51,cmap=cmm)

      #
      # Back-rotate system:
      cs = np.cos(-th)
      sn = np.sin(-th)
      y2 = x2 * 0.
      x_sec = x2 * cs + y2 * sn
      y_sec = x2 * (-sn) + y2 * cs
      yf = xf * 0.
      x_f = xf * cs + yf * sn
      y_f = xf * (-sn) + yf * cs


      x_sec_2d = x_sec.reshape(num, num)
      y_sec_2d = y_sec.reshape(num, num)
      z_2_2d = z_2.reshape(num, num)
      i_2_2d = i_2.reshape(num, num)
      x_f_2d = x_f.reshape(num, num)
      y_f_2d = y_f.reshape(num, num)

      cmm = cm.magma
      surf2 = ax.plot_surface(x_sec_2d, y_sec_2d, z_2_2d \
          ,linewidth=0, antialiased=False, facecolors=cmm(get_color(i_2_2d,im2max)) \
          , rcount=61, ccount=61, cmap=cmm)


###      itpt = x.size//2
###      plt1.append([[x[itpt], x[itpt]], [y[itpt],y[itpt]], [telescope.dmm,z_1[itpt]]])
###      plt2.append([[x[itpt], x_sec[itpt]], [y[itpt],y_sec[itpt]], [z_1[itpt], z_2[itpt]]])
###      plt3.append([[x_sec[itpt], x_f[itpt]], [y_sec[itpt],y_f[itpt]], [z_2[itpt], zf]])

      #ax.plot([x[itpt], x[itpt]], [y[itpt],y[itpt]], [telescope.dmm,z_1[itpt]], lw=1)
      #ax.plot([x[itpt], x_sec[itpt]], [y[itpt],y_sec[itpt]], [z_1[itpt], z_2[itpt]], lw=1)
      #ax.plot([x_sec[itpt], x_f[itpt]], [y_sec[itpt],y_f[itpt]], [z_2[itpt], zf], lw=1)

      isin=secondary.f2d.ev(x_sec, y_sec)
      ww = np.where(isin > 0.5)[0]

      segments.i['%i' % (it_sgm,)]=i_1[ww]
      segments.i2['%i' % (it_sgm,)]=i_2[ww]
      segments.th1['%i' % (it_sgm,)]=np.arctan2(y, x)[ww]#th_1[ww]

      segments.rad['%i' % (it_sgm,)] = np.sqrt(x[ww]**2+y[ww]**2)
      segments.area['%i' % (it_sgm,)] = segments.d11d * segments.d21d

      ntot=ntot+ww.size
    #


##    for ipl1, ipl2, ipl3 in zip(plt1, plt2, plt3):
##      ax.plot(*ipl1, lw=1)
##      ax.plot(*ipl2, lw=1)
##      ax.plot(*ipl3, lw=1)
    from pdb import set_trace as stop

    print('Total number of rays considered: %i' % ntot)

    cax1 = pl.axes([0.05,0.05,0.01,0.4])
    cax2 = pl.axes([0.05,0.55,0.01,0.4])

    cbar1 = fg.colorbar(surf1,cax=cax1,boundaries=np.linspace(0,im1max,256)\
        ,ticks=np.round(np.linspace(0.,im1max,5)))
    cbar1.ax.set_ylabel(r'i$_{1}$')
    cbar2 = fg.colorbar(surf2,cax=cax2,boundaries=np.linspace(0,im2max,256)\
        ,ticks=np.round(np.linspace(0.,im2max,5)))
    cbar2.ax.set_ylabel(r'i$_{2}$')

    pl.show()
    pl.draw()
    stop()

    return
#
#
#
# ---------------------------------------------------------------------
#
#
#
# ---------------------------------------------------------------------
#

def new_plot_primary(telescope, beam, secondary, segments \
    , pax=None, prays=False, fcol='', layout='polar'):
#
  import numpy as np
  import PyColors as PyC

  def plot_primary_segment(iteles, ixpos, iypos, i2ax, ifcol):
#
    if (iteles.primary_shape == 'hexagonal'):

      dy = np.sqrt(iteles.radius**2-(iteles.radius/2.)**2)

      p1x = ixpos + iteles.radius
      p1y = iypos + 0.

      p2x = ixpos + iteles.radius / 2.
      p2y = iypos + dy

      p3x = ixpos - iteles.radius / 2.
      p3y = iypos + dy

      p4x = ixpos - iteles.radius
      p4y = iypos + 0.

      p5x = ixpos - iteles.radius / 2.
      p5y = iypos - dy

      p6x = ixpos + iteles.radius / 2.
      p6y = iypos - dy

      ix = [ p1x, p2x, p3x, p4x, p5x, p6x, p1x ]
      iy = [ p1y, p2y, p3y, p4y, p5y, p6y, p1y ]

      if (len(ifcol)==0):

        im = i2ax.plot(ix,iy,color='k',lw=1.)

      else:

        from pdb import set_trace as stop
        #stop()

        ix1 = [p4x, p3x, p2x, p1x]
        iy1 = [p4y, p3y, p2y, p1y]
        iy2 = [p4y, p5y, p6y, p1y]

        im = i2ax.fill_between(ix1, iy1, iy2, color=ifcol[0])
        im = i2ax.plot(ix,iy,color='k',lw=1.)

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

      ww = np.where(dumrad >= iteles.rad_min)[0]

      dumrad = dumrad[ww]
      resgx = resgx[ww]
      resgy = resgy[ww]

      im = []

# Otherwise:
    else:
      print('Primary shape unknown!')
      print('')
      print('')
      exit()
#
    return im


  #
  # Main: plot_primary
  #

  #layout='cartesian'

  x1d, y1d, d1s, d2s = get_primary_rays(layout, telescope)
#######  if ( (layout!='polar') & (layout!='cartesian') ):
#######    print('')
#######    print('\t Unknown layout: %s for ray tracing!' % (layout,))
#######    print('\t Allowed options are: "polar", "cartesian"')
#######    print('')
#######    exit(1)

  if (np.any(pax)==None):
    import matplotlib.pyplot as pl
    pl.ion()
    pl.figure(1)
    pl.clf()
    fg,iax=pl.subplots(ncols=1, nrows=1, num=1)
  else:
    iax=pax

  ntot=0
  if (fcol != ''):
    tpfcol=PyC.get_colors(segments.Nfamily, cmap=fcol)
    aux = []
    for itn in range(segments.N):
      aux.append([tpfcol[np.int32(segments.families[itn])-1]])
  # Put color legend in GTC layout:
    for it_fml in range(segments.Nfamily):
      iax.plot([],[],linewidth=5, label=r'Fam.: %i' % (it_fml+1,) \
          , color=tpfcol[it_fml])
    iax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0) \
      , fancybox=True, shadow=True,ncol=2)

  else:
    aux = [[],]*segments.N

  for it_sgm in range(segments.N):
    itim = plot_primary_segment(telescope\
       , segments.xpos[it_sgm], segments.ypos[it_sgm], iax, aux[it_sgm])


  if (prays != False):

###    import matplotlib.pyplot as pl
###    pl.ion()
###    pl.close(1)
###    fg, ax = pl.subplots(nrows=1,ncols=1,num=1)

# First of all, we remove all the rays that come from behind the secondary:
    isin=secondary.f2d.ev(x1d, y1d)
    iax.scatter(x1d,y1d,color='r', marker='o', s=10, edgecolor=None)
    ww = np.where(isin < 0.5)[0]
    if (len(ww)>0):
      x1d = x1d[ww]
      y1d = y1d[ww]

    iax.scatter(x1d,y1d,color='k', marker='o', s=7, edgecolor=None)


    def check_primary_id(iteles, ixpos, iypos, igxs, igys):
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

        ww = np.where(dumrad >= iteles.rad_min)[0]

        dumrad = dumrad[ww]
        resgx = resgx[ww]
        resgy = resgy[ww]

# Otherwise:
      else:
        print('Primary shape unknown!')
        print('')
        print('')
        exit()
#
      return resgx + ixpos, resgy + iypos

#
    ntot=0

    x1df = x1d * 0.
    y1df = y1d * 0.

    for it_sgm in range(segments.N):

      x, y = check_primary_id(telescope\
         , segments.xpos[it_sgm], segments.ypos[it_sgm], x1d, y1d)

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
      dir_cos_i = get_direction_cosines(beam.l_inc,beam.m_inc,beam.n_inc)[:,None]*np.ones(x1.size)[None,:]
      dir_cos_n1, z_1 = new2_get_gcurv(x1,k=telescope.k_1,r=telescope.rc_1,norm=True)

      i_1 = np.arccos(np.sum(dir_cos_n1*dir_cos_i,axis=0))

      dir_cos_o1 = rotate_2d(dir_cos_i, -2 * i_1 * np.sign(x1))

      x2, z_2, dir_cos_n2 = get_intersection(dir_cos_o1, x1, z_1, telescope.nrc_2, telescope.k_2, telescope.dmm)

      i_2 = np.arccos(np.sum(dir_cos_n2*dir_cos_o1,axis=0))

      #
      # Back-rotate system:
      cs = np.cos(-th)
      sn = np.sin(-th)
      y2 = x2 * 0.
      x_sec = x2 * cs + y2 * sn
      y_sec = x2 * (-sn) + y2 * cs

      isin=secondary.f2d.ev(x_sec, y_sec)
      ww = np.where(isin > 0.5)[0]

      x1df[ntot:ntot+ww.size] = x[ww] * 1.
      y1df[ntot:ntot+ww.size] = y[ww] * 1.

      ntot=ntot+ww.size

    iax.scatter(x1df[0:ntot],y1df[0:ntot],color='cyan', marker='o', s=4, edgecolor=None)
    print(np.arctan2(x1df[0:ntot],y1df[0:ntot])*180./np.pi)


    print('Total number of rays considered: %i' % ntot)




















  return


def plot_mueller_elements(segments, materials,cdust=False):
  import numpy as np
  import time as tm
  from pdb import set_trace as stop
  import apy_utils as uts
#
  tarea = 0.
  for it in range(segments.N):
    tarea += np.sum(segments.area['%i' % (it,)])

#RUN TIMES:
  res = np.zeros((segments.SimulatedTimes, segments.N, 4, 4))
  avg_res = np.zeros((segments.SimulatedTimes, 4, 4))
  for it_tms in range(segments.SimulatedTimes):
    print("Timestep: %i of %i" % (it_tms, segments.SimulatedTimes))
#RUN SEGMENTS:
    npts_seg = np.zeros((segments.N,))
    ntot2=0

#
# Testing:
    pl.close(1)
    pl.close(2)
    pl.close(3)
    pl.close(4)
    #fg, ax = pl.subplots(num=1, nrows=6, ncols=6)
    #ax=ax.flatten()
    fg1, ax1 = pl.subplots(num=1, nrows=1, ncols=1)
    fg2, ax2 = pl.subplots(num=2, nrows=1, ncols=1)
    fg3, ax3 = pl.subplots(num=3, nrows=1, ncols=1)
    fg4, ax4 = pl.subplots(num=4, nrows=1, ncols=1)
# Testing.
#


    for it_sgm in range(segments.N):
#
      if (segments.i['%i' % (it_sgm,)].size<1):
        npts_seg[it_sgm] = 0.
        res[it_tms,it_sgm,:,:] = res[it_tms,it_sgm,:,:] * 0.
        continue

      mmpa_full = hex_syst(segments.i['%i' % (it_sgm,)], \
          segments.th1['%i' % (it_sgm,)], \
          segments.i2['%i' % (it_sgm,)], \
          materials.lamb, materials.n_s, materials.n_c\
          , materials.n_ox, materials.n_d, materials.n_in, \
          materials.d_s[it_sgm, it_tms], \
          materials.d_c[it_sgm, it_tms]\
          , materials.d_ox[it_sgm, it_tms], materials.d_d[it_sgm, it_tms])

      #GET Mdust
      if (cdust==True):
        ###mmpd = get_mdust(materials.n_d, materials.d_d[it_sgm, it_tms])
        mmpd = new_get_mdust(0.1, materials.d_d[it_sgm, it_tms])
        res_full = mmpa_full.dot(mmpd)
      else:
        res_full = mmpa_full * 1.


#
# Testing:

      xx = segments.rad['%i' % (it_sgm,)]*np.cos(segments.th1['%i' % (it_sgm,)])
      yy = segments.rad['%i' % (it_sgm,)]*np.sin(segments.th1['%i' % (it_sgm,)])
      im1 = ax1.scatter(xx, yy, marker='s', s=10, c=res_full[:,1,0] * segments.area['%i' % (it_sgm,)] / tarea, cmap='RdGy', vmin=-0.0000001, vmax=0.0000001)
      im2 = ax2.scatter(xx, yy, marker='s', s=10, c=res_full[:,2,0] * segments.area['%i' % (it_sgm,)] / tarea, cmap='RdGy', vmin=-0.0000001, vmax=0.0000001)
      im3 = ax3.scatter(xx, yy, marker='s', s=10, c=res_full[:,1,0], cmap='RdGy', vmin=-0.0001, vmax=0.0001)
      im4 = ax4.scatter(xx, yy, marker='s', s=10, c=res_full[:,2,1], cmap='RdGy', vmin=-0.0001, vmax=0.0001)
      #ax[it_sgm].hist(res_full[:,1,0], bins=101)
# Testing.
#


      nn,_,_ = res_full.shape
      ntot2=ntot2+nn
########      for itn in range(nn):
########        res[it_tms,it_sgm,:,:]=res[it_tms,it_sgm,:,:] \
########            + res_full[itn,:,:] * segments.rad['%i' % (it_sgm,)][itn]
########      res[it_tms,it_sgm,:,:] = res[it_tms,it_sgm,:,:] \
########          * segments.area['%i' % (it_sgm,)]
########      #
########      # Area to normalize:
########      npts_seg[it_sgm] = np.sum(segments.rad['%i' % (it_sgm,)] \
########          * segments.area['%i' % (it_sgm,)])

      for itn in range(nn):
        res[it_tms,it_sgm,:,:]=res[it_tms,it_sgm,:,:] \
            + res_full[itn,:,:] * segments.area['%i' % (it_sgm,)][itn]
      #
      # Area to normalize:
      npts_seg[it_sgm] = np.sum(segments.area['%i' % (it_sgm,)])


    print('Total number of rays considered: %i' % ntot2)
    print('Total area considered: %.4f==%.4f?' % (npts_seg.sum(), tarea,))
    

    # Since we consider all the points of each segment, the average must...
    # ... be done accordingly:
    avg_res[it_tms,:,:] = np.sum(res[it_tms,:,:,:], axis=0) \
        / np.sum(npts_seg)

#
  cbar1 = pl.colorbar(im1, ax=ax1)
  cbar2 = pl.colorbar(im2, ax=ax2)
  cbar3 = pl.colorbar(im3, ax=ax3)
  cbar4 = pl.colorbar(im4, ax=ax4)
  for ax in [ax1, ax2, ax3, ax4]:
    ax.set_aspect('equal')
    ax.grid()
  pl.draw()
  stop()
  return avg_res, res, npts_seg


#
# ---------------------------------------------------------------------
#



##########
