
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm> 
#include <omp.h>

namespace mrr{
  //
  //constexpr double pi() {return std::atan(1)*4;}
  double pi() {return std::atan(1)*4;}
  //
  //
  template <typename Tdr> void get_rotation(Tdr const &ang, Tdr* mat){

    for (int i=0; i<4*4; ++i){mat[i]=0;}
    
    //  0  1  2  3
    //  4  5  6  7
    //  8  9 10 11
    // 12 13 14 15

    mat[0] = 1.;
    mat[15] = 1.;

    Tdr const ca = std::cos(2.*ang);
    Tdr const sa = std::sin(2.*ang);

    mat[5] = ca;
    mat[6] = -sa;
    mat[9] = sa;
    mat[10] = ca;

  }; // get_rotation;
  //
  template <typename Tdr> void get_reflection( Tdr const &r
      , Tdr const &x, Tdr const &tau, Tdr* mat){

    for (int i=0; i<4*4; ++i){mat[i]=0;}
    
    //  0  1  2  3
    //  4  5  6  7
    //  8  9 10 11
    // 12 13 14 15

    mat[0] = pow(x,2) + 1;
    mat[1] = mat[0] - 2;
    mat[4] = mat[1];
    mat[5] = mat[0];

    Tdr const ct = std::cos(tau);
    Tdr const st = std::sin(tau);

    mat[10] = 2 * x * ct;
    mat[11] = 2 * x * st;
    mat[14] = -mat[11];
    mat[15] = mat[10];

    for (int i=0; i<4*4; ++i){mat[i]=mat[i]*0.5*r;}

  }; // get_reflection.

  //
  template <typename Tnr, typename Tdr> void mat_prod(
      Tnr const &i, Tnr const &n, Tdr* const lt, Tdr* mueller){

    /*
      Offset (elements) on the full mueller matrix
    */
    Tnr sz = n*n;
    Tnr off=i*(sz);
    /*
      Allocate and initialize (0) result vector:
    */
    Tdr* res = new Tdr[sz];
    for (Tnr k=0; k<sz; ++k){res[k]=0;}

    /*
      Matrix product:
    */
    for (Tnr k=0; k<sz; ++k){

      Tnr s = k / n;
      Tnr f = k % n;

      for (Tnr m=0; m<n; ++m){

        Tnr i1 = s*n+m;
        Tnr i2 = m*n+f;

        res[k] = res[k] + lt[i1] * mueller[off+i2];

      }
    }

    /*
      Store result on the full mueller matrix
    */
    for (Tnr k=0; k<sz; ++k){mueller[off+k] = res[k];}
    delete [] res;
  }; // mat_prod.
  //
  template <typename Tdr> Tdr snell( Tdr const &in1
      , Tdr const &ang1, Tdr const &in2){
    return asin( in1 * sin(ang1) / in2 );
  }; // snell.
  //
  template <typename Tdr> Tdr get_ncosi( Tdr const &n
      , Tdr const &i, bool const &perp){

    Tdr res=0;

    if (perp){
      res = cos(i) / n;
    } else {
      res = n * cos(i);
    }

    return res;

  }; // get_ncosth.
  //
  template <typename Tdr, typename Tlr> void get_layered_reflection(
      Tdr const &i, Tlr const &nlayers, Tdr* const ln, Tdr* const lk
      , Tdr* const ld, Tdr const &lambda, Tdr &r, Tdr &x, Tdr &tau){

    std::complex<Tdr>* jones = new std::complex<Tdr>[4]();
    std::complex<Tdr>* itj = new std::complex<Tdr>[4]();

    std::complex<Tdr> const cn_air(ln[0],0);
    std::complex<Tdr> const ci_air(i,0);

    std::complex<Tdr> const cn_subs(ln[nlayers-1],lk[nlayers-1]);
    std::complex<Tdr> const th_subs = snell(cn_air, ci_air, cn_subs);

    std::complex<Tdr> const cc(0,1);

    bool perp = false;

    for (int i=0; i<4; ++i){jones[i]=(i%2 == i/2 ? 1 : 0);}

    /*

    */
    for (int ilay=nlayers-2; ilay>0 ; --ilay){

      std::complex<Tdr> cn(ln[ilay],lk[ilay]);
      std::complex<Tdr> th = snell(cn_air, ci_air, cn);
      std::complex<Tdr> p = get_ncosi(cn, th, perp);

      std::complex<Tdr> beta = 2 * pi() / lambda * cn * ld[ilay] * cos(th);

      /*
        0 1
        2 3
      */
      itj[0] = cos(beta);
      itj[1] = -cc * sin(beta) / p ;
      itj[2] = -cc * sin(beta) * p;
      itj[3] = itj[0];

      mat_prod<int,std::complex<Tdr>>(0,2,itj,jones);

    }; // Loop covering the various layers.

    Tdr* modulus = new Tdr[2]();
    Tdr* phase = new Tdr[2]();

    for (int ipol=0 ; ipol<2; ++ipol){
      /*

      */
      if (ipol>0){perp=true;}
      /*

      */
      std::complex<Tdr> const p_air = get_ncosi(cn_air, ci_air, perp);
      std::complex<Tdr> const p_subs = get_ncosi(cn_subs, th_subs, perp);

      /*
        Once the equivalent Characteristic matrix is retrieved, we can proceed:
      */
      std::complex<Tdr> tr1 = (jones[0] + jones[1] * p_subs) * p_air;
      std::complex<Tdr> tr2 = jones[2] + jones[3] * p_subs;

      std::complex<Tdr> cr = (tr1 - tr2) / (tr1 + tr2);

      phase[ipol] = atan2(std::imag(cr),std::real(cr));
      modulus[ipol] = std::abs(cr);

    }; // polarization components.

    r = pow(modulus[1], 2);
    x = std::real(modulus[0]/modulus[1]);
    tau = std::real(phase[0])-std::real(phase[1]);

    /*
      Clean
    */
    delete [] jones;
    delete [] itj;
    delete [] modulus;
    delete [] phase;

  }; // get_layered_reflection.
  //
  template <typename Tdr> void get_ideal_reflection( Tdr const &i
      , Tdr const &n, Tdr const &k, Tdr &r, Tdr &x, Tdr &tau){

    Tdr const si = sin(i);
    Tdr const ci = cos(i);
    Tdr const ti = tan(i);
    Tdr const si2 = si*si;
    Tdr const ci2 = ci*ci;
    Tdr const n2 = n*n;
    Tdr const k2 = k*k;

    Tdr const tmp = pow(pow(n2 - k2 - si2, 2) + 4 * n2 * k2, 0.5);

    Tdr const f = 0.5 * (n2 - k2 - si2 + tmp);
    Tdr const g = 0.5 * (k2 - n2 + si2 + tmp);

    Tdr const r_p1 = f + g - 2 * pow(f,0.5) * ci + ci2;
    Tdr const r_p2 = f + g + 2 * pow(f,0.5) * ci + ci2;
    r = r_p1 / r_p2;

    Tdr const s_tau = 2 * pow(g,0.5) * si * ti;
    Tdr const c_tau = si2 * ti*ti - (f + g);
    tau = atan(s_tau / c_tau);

    Tdr const x_1 = f + g - 2 * pow(f,0.5) * si * ti + si*si * ti*ti;
    Tdr const x_2 = f + g + 2 * pow(f,0.5) * si * ti + si*si * ti*ti;
    x = pow(x_1 / x_2, 0.5);

  }; // get_ideal_reflection.
  //
  template <typename Tnr, typename Tnl, typename Tdr> void do_ray(Tnr const &i
      , Tdr const &ii1, Tdr const &ith, Tdr const &ii2, Tdr const &lamb
      , Tnl const& nlayers, Tdr* const ln, Tdr* const lk, Tdr* const ld
      , Tdr* result){

    Tdr* mat = new Tdr[16]();

    Tdr r=0;
    Tdr x=0;
    Tdr tau=0;

    Tdr const ang1 = -(ith + pi()/2.);
    Tdr const ang2 = pi();
    Tdr const ang3 = -(pi()/2. - ith);
  
    for (int l=0; l<5; ++l){
  
      if (l==0){
        get_rotation<Tdr>(ang1,mat);
      } else if (l==1){
        get_layered_reflection<Tdr,Tnl>(ii1, nlayers, ln, lk, ld
            , lamb, r, x, tau);
        get_reflection<Tdr>(r, x, tau, mat);
      } else if (l==2){
        get_rotation<Tdr>(ang2,mat);
      } else if (l==3){
        /*
          Secondary mirror:
          ln and lk store complex refraction indexes for: substrate, conductor, ...
        */
        Tnl cond=nlayers-2;
        get_ideal_reflection<Tdr>(ii2, ln[cond], lk[cond], r, x, tau);
        get_reflection<Tdr>(r, x, tau, mat);
      } else if (l==4){
        get_rotation<Tdr>(ang3,mat);
      }
  
      mat_prod<Tnr,Tdr>(i,4,mat,result);

    }; // Linear transformations.

    delete [] mat;

  }; // do_ray.
  //
  template <typename Tnr, typename Tnl, typename Tdr> void eval_mueller_segment(
      Tnr const& nrays, Tnl const& nlayers, Tdr* const i1, Tdr* const th
      , Tdr* const i2, Tdr const& lamb, Tdr* const ln, Tdr* const lk, Tdr* const ld
      , int const nthread, Tdr* result){

    /*
      1: initialize diagonal elements to 1:
    */
    for (Tnr i=0; i<nrays; ++i){
      for (Tnl j=0; j<4; ++j){
        Tnr k=i*(4*4)+(j*4+j);
        result[k]=1;
      }
    }

    /*
      2: Evaluate each linear transformation sequentially:
    */

    Tnr i;
    #pragma omp parallel default(shared) private(i) num_threads(nthread)
    {
      #pragma omp for schedule(guided)
      for (i=0; i<nrays; ++i){

        do_ray<Tnr,Tnl,Tdr>(i,i1[i],th[i],i2[i],lamb,nlayers,ln,lk,ld,result);
 
      }; // nrays
    }; // pragma parallel.

  }; // eval_mueller_segment.
  //
  template <typename Td1, typename Td2> void eval_it_circle(int const &inside
      , Td1 const &rr, Td1 const & sx, Td1 const &sy, Td2 &res){

    Td1 const sr = std::pow(sx,2) + std::pow(sy,2);
    Td2 const in = (sr<rr ? 1 : 0);
    Td2 const out = (sr>rr ? 1 : 0);
    res = res * (inside==1 ? in : out);

  }; // eval_it_circle.
  //
  template <typename Tnr, typename Tdr> void eval_pts_circle(
      int const &nthread, int const &inside, Tdr const &radius
      , Tnr const &nrays, Tdr* const xv, Tdr* const yv, Tnr* result){

    Tdr const rad2=std::pow(radius,2);
    Tnr i;
    #pragma omp parallel default(shared) private(i) num_threads(nthread)
    {
      #pragma omp for schedule(guided)
      for (i=0; i<nrays; ++i){
        eval_it_circle<Tdr, Tnr>(inside, rad2, xv[i], yv[i], result[i]);
      }; // nrays.
    }; // pragma parallel.

  }; // eval_pts_anular.
  //
  template <typename Tnd, typename Ttd> Ttd get_mean(Tnd const &nn
      , Ttd* const vv){
    Ttd avg = 0;
    Ttd fact = 1 / (Ttd)(nn);
    for (Tnd i=0; i<nn; ++i){avg+=vv[i]*fact;}
    return  avg;
  }; //
  //
  template <typename Tn, typename Td> void argsort(Tn const &n, Td* const ang, Tn* w){

    std::vector<Tn> indices(n,0);
    for (Tn i=0; i<n; i++){indices[i]=i;}

    std::stable_sort(indices.begin(), indices.end(),
        [&ang](Tn i1, Tn i2) {return ang[i1] < ang[i2];});

    for (Tn i=0; i<n; i++){w[i]=indices[i];}

  }; // argsort.
  //
  template <typename Tnd, typename Ttd> Ttd get_area_polygon(Tnd const &nn
      , Ttd* const xv, Ttd* const yv){

    Ttd area = 0;
    Ttd const xmean = get_mean(nn, xv);
    Ttd const ymean = get_mean(nn, yv);

    Ttd* xv0 = new Ttd[nn]();
    Ttd* yv0 = new Ttd[nn]();
    Ttd* ang = new Ttd[nn]();
    Tnd* ww = new Tnd[nn]();

    for (Tnd i=0; i<nn; ++i){
      xv0[i] = xv[i] - xmean;
      yv0[i] = yv[i] - ymean;
      ang[i] = std::atan2(yv0[i], xv0[i]);
    }    

    argsort<Tnd, Ttd>(nn, ang, ww);

    for (Tnd i=0; i<nn; ++i){
      xv0[i] = xv[ww[i]] - xmean;
      yv0[i] = yv[ww[i]] - ymean;
    }

    Ttd t1 = 0;
    Ttd t2 = 0;

    for (Tnd i=0; i<nn; ++i){
      t1 += xv0[i] * yv0[(i+1)%nn];
      t2 += yv0[i] * xv0[(i+1)%nn];
    }

    area = t1 - t2;
    area = 0.5 * (area<0 ? -area : area);

    delete [] xv0;
    delete [] yv0;
    delete [] ang;
    delete [] ww;

    return area;

  }; // get_area_polygon.
  //
  template <typename Td0, typename Td1, typename Td2> void eval_it_polygon(
      int const &inside
      , Td0 const &nverts, Td1* const xverts, Td1* const yverts, Td1 const apoly
      , Td1 const & sx, Td1 const &sy, Td2 &res){

    /*
      1- Area of the polygon: apoly
    */

    /*
      2- Extend the polygon to include the target point
    */
    Td0 const nve=nverts+1;
    Td1* xve = new Td1[nve]();
    Td1* yve = new Td1[nve]();

    xve[0] = sx;
    yve[0] = sy;
    for (Td0 i=1; i<nve; i++){
      xve[i] = xverts[i-1];
      yve[i] = yverts[i-1];
    }

    Td1 const earea = get_area_polygon<Td0,Td1>(nve,xve,yve);

    Td2 const in = (earea<apoly ? 1 : 0);
    Td2 const out = (earea>apoly ? 1 : 0);

    res = res * (inside==1 ? in : out);

    delete [] xve;
    delete [] yve;

  }; // eval_it_polygon.
  //
  template <typename Tnp, typename Tnr, typename Tdr> void eval_pts_polygon(
      int const &nthread, int const &inside
      , Tnp const &nverts, Tdr* const xverts, Tdr* const yverts
      , Tnr const &nrays, Tdr* const xv, Tdr* const yv, Tnr* result){

    //std::cout<<" Inside="<<inside<<" ; radius="<<radius<<" , nrays="<<nrays<<std::endl;

    Tdr const parea = get_area_polygon<Tnp,Tdr>(nverts, xverts, yverts);

    Tnr i;
    #pragma omp parallel default(shared) private(i) num_threads(nthread)
    {
      #pragma omp for schedule(guided)
      for (i=0; i<nrays; ++i){
        eval_it_polygon<Tnp, Tdr, Tnr>(inside, nverts, xverts, yverts, parea
            , xv[i], yv[i], result[i]);
      }; // nrays.
    }; // pragma parallel.

  }; // eval_pts_polygon.
  //

}; // mirror namespace
