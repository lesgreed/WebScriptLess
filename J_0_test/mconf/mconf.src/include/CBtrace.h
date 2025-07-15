#ifndef MC_CBtrace_
#define MC_CBtrace_

#include "CVector3d.h"
#include "rkf45.h"

namespace MConf {

#define USE_RKF45 1

static const int nODE=2;
static const double ATOL=1e-10; //1e-8;
static const double RTOL=1e-10; //1e-9;

//***************************************************************************
/// The class is used for B-field line tracing.
/// The following equations are solved: \n
///  dr/dfi=Br*r/Bfi \n
///  dz/dfi=Bz*r/Bfi
///  \anchor usageOfCBtrace
///  \par Usage:
/// \code
/// #include "CStconfig.h"
/// int main() {
///   MConf::CStconfig mConf;
///   mConf.load ("w7x-sc1.bin4");
///   if(!mConf.isOK()) exit(1);   // exit with error
///   mConf.truncate(1e-6);        //truncate spectrum
///   mConf.setAccuracy(1e-3);     //set accuracy of transformation
///
///  //Function object returns vector of magnetic field for point in cylindrical coordinates.
///   class Bfield {
///     MConf::CStconfig *mC;
///   public:
///     Bfield(MConf::CStconfig &magConf) { mC = &magConf; }
///   // The operator returns vector of magnetic field for point in cylindrical coordinates
///     Vector3d operator() ( const Vector3d &cyl ) const {
///       return mC->getBcyl(cyl);
///     }
///   };
///
///   // create function object that returns magnetic field
///   Bfield Bf(magConf);
///   // create B-trace object
///   MConf::CBtrace<Bfield> btrace(Bf,dfi);
///   MConf::Vector3d cyl(6,0,0);              // starting point in cylindrical coordinates
///   double fi_increment = 1*(pi/180); // 1 degree increment in toroidal angle
///   btrace.advance(fi_increment,cyl); // move to new point, cyl is updated
/// }
/// \endcode

template <class Bfield> class CBtrace
{
  void init_(double fiStep=0.0017453,double atol=ATOL,double rtol=RTOL) {
    twopi=6.283185307179586476925286766559;
    this->a_tol=atol;
    this->r_tol=rtol;
    h = fiStep==0?1e-3:fiStep;
    q = 1;
    Bfld = 0;
#if USE_RKF45
    rhs.init(Bfld);
    rkf.init(&rhs,2,atol,rtol);
#endif
  }
public:
  CBtrace() { init_(); }
  /// The constructor creates an instance of tracing object.
  /// @param B is the function object that returns vector of magnetic field
  ///   for given point in cylindrical coordinates,
  /// @param fiStep is the Runge-Kutta step along magnetic field line
  ///   in toroidal direction.
  ///   \note Bfield object have to provide B-vector:
  ///  \code
  ///   MConf::CStconfig magConf
  ///   Bfield B(magConf);
  ///   Vector3d cyl(6,0,0);
  ///   Vector3d Bcyl = B(cyl);
  ///  \endcode
  /// \sa MConf::CStconfig::Bfield and \ref usageOfCBtrace "usage of CBtrace"
  CBtrace(Bfield &B, double fiStep=0.0017453) {
    init_(fiStep);
    this->Bfld=&B;
#if USE_RKF45
    rhs.init(&B);
    rkf.init(&rhs,2);
#endif
  }
  virtual ~CBtrace() {;}
  void init(Bfield &B,double fiStep=0.0017453,double atol=ATOL,double rtol=RTOL) {
    init_(fiStep,atol,rtol);
    this->Bfld=&B;
#if USE_RKF45
    rhs.init(&B);
    rkf.init(&rhs,2,atol,rtol);
#endif
  }
  void init() {
#if USE_RKF45
    rkf.init();
#endif
  }
  void setStep(double fiStep) { h = fiStep==0?1e-3:fiStep;  }
  //***************************************************************************
  // input: dfi -- cylindrical angle increment
  // input: cyl -- cylindrical coordinates of magnetic field line
  // output cyl -- cylindrical coordinates of magnetic field line after dfi
  void advance(double dfi, Vector3d &cyl)
  {
    q = 1;
    double fi;
    y[0] = cyl[0]; // r; // set initial values
    fi   = cyl[1]; // fi
    y[1] = cyl[2]; // z;
    double fiEnd=fi+dfi;
#if USE_RKF45
    rhs.setq(q);
    rkf.advance(y,fi,fiEnd);
#else
    RKutt(fi,fiEnd,h);    // solve from fi to fiEnd
#endif
  //  Adams4(fi,fiEnd,h);    // solve from fi to fiEnd
    cyl[0] = y[0];   // r
    cyl[1] = fiEnd;  // fi
    cyl[2] = y[1];   // z
  }

  //***************************************************************************
  // Tracing along poloidal angle is valid only for straight field line coordinates.
  // advancePoloidally -- use for tokamak only!
  // input: dtheta -- poloidal angle increment
  // input: q -- safety factor
  // input: cyl -- cylindrical coordinates of magnetic field line
  //        cyl[1] is the poloidal angle theta: nothing depend on cyl.angle, so we use cyl[1] for poloidal angle
  // output: cyl -- cylindrical coordinates of magnetic field line after dtheta-step
  //        cyl[1] is the poloidal angle theta: nothing depend on cyl.angle, so we use cyl[1] for poloidal angle
  void advancePoloidally(double dtheta, double q, Vector3d &cyl)
  {
    this->q=q;
    double th;
    y[0] = cyl[0]; // r; // set initial values
    th   = cyl[1]; // theta: nothing depend on cyl.angle, so we use cyl[1] for poloidal angle
    y[1] = cyl[2]; // z;
    double thEnd=th+dtheta;
#if USE_RKF45
    rhs.setq(q);
    rkf.advance(y,th,thEnd);
#else
    RKutt(th,thEnd,h);    // solve from th to thEnd
#endif
  //  Adams4(th,thEnd,h);    // solve from th to thEnd
    cyl[0] = y[0];   // r
    cyl[1] = thEnd;  // theta
    cyl[2] = y[1];   // z
  }

private:
  Bfield *Bfld;  // Function object returns vector of magnetic field for given point in cylindrical coordinates.

  double y[nODE];                // we solve dy/dx = F(y,x)
  double a_tol, r_tol;
  double h;                       // step
  double twopi;
  double q;     // 1/iota;

#if USE_RKF45
  class RightHand: public RightHandSide {
    Bfield *Bfld;
    double q;     // 1/iota;
  public:
    RightHand() { Bfld = 0; q = 1;}
    void init(Bfield *B, double q=1) { Bfld = B; this->q = q;}
    void setq(double q) { this->q = q;}
    void operator() (double fi,const double *y,double *yp) { 
      double w;
      Vector3d c(y[0],fi,y[1]); // c==(r,fi,z)
      Vector3d B = (*Bfld)(c);
      if(B[1]==0)
        w=0;
      else
        w=q*y[0]/B[1];     // q*r/Bfi, where q = 1/iota
      yp[0]=w*B[0];        // dr/dfi = q*r*Br/Bfi
      yp[1]=w*B[2];        // dz/dfi = q*r*Bz/Bfi     
    }
  };

  RightHand rhs;
  rkf45 rkf;
#else
  double y1[nODE],u[nODE],r[nODE];            //
  double f0[nODE],f1[nODE],f2[nODE],f3[nODE]; // f0[nODE], .. -- for Adams method
  //***************************************************************************
  // B-field line tracing
  //  dr/dth=Br*q*r/Bfi
  //  dz/dth=Bz*q*r/Bfi
  //
  //  dy = dfi*F(y,fi), where dy/dfi = F(y,fi)
  //  q is the safity factor.
  //  fi is the cylindrical(toroidal) angle
  //  fi = q*th
  //  th is the poloidal angle
  //  normally q=1, in this case th must be treated as toroidal angle fi. 
  // Tracing along poloidal angle is valid only for straight field line coordinates.
  void RHS(double fi,double dfi,double *y,double *dy)
  {
  //  double r = y[0];
  //  double z = y[1];
    Vector3d c(y[0],fi,y[1]);
    Vector3d B = (*Bfld)(c);
    if(B[1]==0)
      dfi=0;
    else
      dfi *= q*y[0]/B[1];  // q*r/Bfi, where q = 1/iota
    dy[0]=dfi*B[0];        // dr = dfi*q*r*Br/Bfi
    dy[1]=dfi*B[2];        // dz = dfi*q*r*Bz/Bfi
  }

  //***************************************************************************
  // input: y[] at x; dx is the step
  // output y[] at xend
  void RKutt(double x, double xend,double dx)
  {
    int n,N = int(fabs((xend-x)/dx));
    if(N>0)
      dx = (xend-x)/N;
    else {
      dx=xend-x;
      N=1;
    }
    for(n=0;n<N;n++) RKutt1step(x,dx);
  }

  //***************************************************************************
  // one step integrator
  // input:  y[] at x; dx is the step
  // output: y[] at x=x+dx, function updates x
  void RKutt1step(double &x, double dx)
  {
    int i;
    RHS(x,dx, y,u);                 // r1 = u = dx*f(y,x)
    for(i=0; i<nODE; i++)
      y1[i] = y[i] + u[i]/2.;       // y1 = y + r1/2
    RHS(x+dx/2,dx, y1,r);           // r2 = dx*f(y+r1/2,x+dx/2)
    for(i=0; i<nODE; i++)  {
      y1[i] = y[i] + r[i]/2.;       // y1 = y + r2/2
      u[i] += 2.*r[i];              // u = r1 + 2*r2
      }
    RHS(x+dx/2,dx, y1,r);           // r3 = dx*f(y+r2/2,x+dx/2)
    for(i=0; i<nODE; i++)  {
      y1[i] = y[i] + r[i];          // y1 = y + r3
      u[i] += 2.*r[i];              // u = r1 + 2*r2 + 2*r3
      }
    RHS(x+dx,dx, y1,r);             // r4 = dx*f(y+r3,x+dx)
    for(i=0; i<nODE; i++)
      y[i] += (u[i] + r[i])/6;      // y = y + (r1 + 2*r2 + 2*r3 + r4)/6
    x += dx;
  }
#endif

  private: // obsolete or not used
#if 0
  //***************************************************************************
  // input: cyl -- cylindrical coordinates of magnetic field line
  // output cyl -- cylindrical coordinates of magnetic field line after toroidal turn
  void advance2pi_not_used(Vector3d &cyl)
  {
    q = 1;
    double fi;
    y[0] = cyl[0]; // r; // set initial values
    fi   = cyl[1]; // fi
    y[1] = cyl[2]; // z;
    double fiEnd=fi+twopi;
    RKutt(fi,fiEnd,h);    // solve from fi to fiEnd
  //  Adams4(fi,fiEnd,h);    // solve from fi to fiEnd
    cyl[0] = y[0];   // r
  //cyl[1] = fiEnd;  // fi
    cyl[2] = y[1];   // z
  }

  //***************************************************************************
  // input: dfi -- cylindrical angle increment
  // input: cyl -- cylindrical coordinates of magnetic field line
  // output cyl -- cylindrical coordinates of magnetic field line after dfi
  void advance1step_not_used(double dfi, Vector3d &cyl)
  {
    q = 1;
    double fi;
    y[0] = cyl[0]; // r; // set initial values
    fi   = cyl[1]; // fi
    y[1] = cyl[2]; // z;
    RKutt1step(fi,dfi);    // solve from fi to fi+dfi
    cyl[0] = y[0];    // r
    cyl[1] = fi;      // fi
    cyl[2] = y[1];    // z
  }

  //***************************************************************************
  // input: y[] at x; dx is the step
  // output y[] at xend
  void Adams4(double x, double xend,double dx)
  {
    static const double p40 = - 9./24.; //predictor coefficients
    static const double p41 =  37./24.;
    static const double p42 = -59./24.;
    static const double p43 =  55./24.;
    static const double c40 =   1./24.; //corrector coefficients
    static const double c41 =  -5./24.;
    static const double c42 =  19./24.;
    static const double c43 =   9./24.;
    int i;
    int n,N = int(fabs((xend-x)/dx));
    N = N<5?5:N;
    dx = (xend-x)/N;
    double *fp0 = f0;
    double *fp1 = f1;
    double *fp2 = f2;
    double *fp3 = f3;
    RHS(x,dx, y,fp0);
    RKutt1step(x,dx);  //RKutt(x,x+dx,dx/4); x+=dx;
    RHS(x,dx, y,fp1);
    RKutt1step(x,dx);  //RKutt(x,x+dx,dx/4); x+=dx;
    RHS(x,dx, y,fp2);
    RKutt1step(x,dx);  //RKutt(x,x+dx,dx/4); x+=dx;
    for(n=4; n<=N; n++) {
  // next 2 lines is the Adams-Bashforth predictor
      RHS(x,dx, y,fp3);
      for(i=0; i<nODE; i++) y1[i] = y[i]+(p43*fp3[i]+p42*fp2[i]+p41*fp1[i]+p40*fp0[i]);
      x += dx;
  // next 2 lines is the Adams-Moulton corrector
      RHS(x,dx, y1,fp0);
      for(i=0; i<nODE; i++) y[i]+=(c43*fp0[i]+c42*fp3[i]+c41*fp2[i]+c40*fp1[i]);
      double *fpw = fp0;
      fp0 = fp1;
      fp1 = fp2;
      fp2 = fp3;
      fp3 = fpw;
    }
  }

  //***************************************************************************
  // input: y[] at x; dx is the step
  // output y[] at xend
  // note: don't use this method, it has bad accuracy
  void Adams3(double x, double xend,double dx)
  {
    int i;
    int n,N = int(fabs((xend-x)/dx));
    N = N<5?5:N;
    dx = (xend-x)/N;
    static const double c0 =   5./12.;
    static const double c1 = -16./12.;
    static const double c2 =  23./12.;

    double *fp0 = f0;
    double *fp1 = f1;
    double *fp2 = f2;
    RHS(x,dx, y,fp0);
    RKutt1step(x,dx);
    RHS(x,dx, y,fp1);
    RKutt1step(x,dx);
    for(n=3; n<=N; n++) {
      RHS(x,dx, y,fp2);
      for(i=0; i<nODE; i++) y[i]+=(c2*fp2[i]+c1*fp1[i]+c0*fp0[i]);
      double *fpw = fp0;
      fp0 = fp1;
      fp1 = fp2;
      fp2 = fpw;
      x += dx;
    }
  }
#endif
};


};   //namespace MConf

#endif
