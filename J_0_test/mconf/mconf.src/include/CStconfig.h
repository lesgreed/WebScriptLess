#ifndef MC_CSTCONFIG_H
#define MC_CSTCONFIG_H

#ifdef _MSC_VER
// These pragmas are to quiet VC++ about the expanded template identifiers exceeding 255 chars.
// You won't be able to see those variables in a debug session, but the code will run normally
#pragma warning( push )
#pragma warning( disable : 4786  )
#endif

#define NO_MCDB

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <stack>
#include "ngarray.h"
#include "CArray1d.h"
#include "CArray2d.h"
#include "CArray3d.h"
#include "CVector3d.h"
#include "CBtrace.h"
#include "threadtypes.h"
#include "mathf.h"

#undef B0  //to fix the ORACLE-headers problems under LINUX

/// Library namespace.
namespace MConf {

class CEfit;
//****************************************************************************
//                          2004-2021,   Author: yuriy turkin at ipp mpg de
/// Magnetic configuration of a stellarator.
/// Class for geometric calculation.
/// This class provides the coordinate transformation between
/// <a class="el" href="index.html#MagneticCoordinates">magnetic coordinates</a> 
/// and real space coordinates.
/// \par Usage:
/// \code
/// #include "CStconfig.h"
/// int main() {
///   CStconfig mConf;
///   mConf.load ("w7x-sc1.bin4"); //load configuration with Boozer coordinates
///   if(!mConf.isOK()) exit(1);   //exit with error
///   mConf.truncate(1e-5);        //truncate spectrum
///   mConf.setAccuracy(epsA);     //set accuracy of transformation
///   mConf.setB0(2.5);            //set minimum value of B on magn. axis
///
///   Vector3d cyl(6,0.01745,0.5);// R=6,fi=1deg,Z=0.5
/// //  coordinate transformation form cylindrical to magnetic coordinates
///   Vector3d boozer = mConf.cyl2mag(cyl); // boozer is (s,theta,phi)
///   Vector3d B = mConf.getBcyl(cyl);         // B-vector in cylindrical coord.
///
///   double s = boozer[0];               // normalized toroidal flux
///   double reff = mConf.r(s);           // effective r in [m]
/// }
/// \endcode
///
/// \todo
///  add function void setWarningLevel(int level);
///
class CStconfig  {
public:
  enum fileFormat  {  /// file-format types, for internal use
    UNDEFINED, ///< unknown file format
    W7X_FMT,   ///< ascii W7X format(J.Geiger)
    W7X_BIN4,  ///< 4-byte(float)  binary W7X format(Y.Turkin)
    W7X_BIN8,  ///< 8-byte(double) binary W7X format(Y.Turkin)
    LHD,       ///< LHD format
    LHD1,      ///< LHD format
    VMEC,      ///< VMEC wout-format; use without transformation to Boozer coordinates
    VMECNETCDF,///< VMEC netCDF wout-format; use without transformation to Boozer coordinates
    W7X_DB,    ///< W7X format from database
    EFIT,      ///< G EQDSK File from EFIT
    W7X_FMT_TSFC,   ///< ascii W7X format(Y.Turkin) for Tokamak Symmetry Flux Coordinates
    W7X_FMT_VMEC,   ///< ascii W7X format(Y.Turkin) for VMEC Coordinates
    W7X_FMT_VMEC_FP,///< ascii W7X format(J.Geiger) for VMEC Coordinates after FP
    TORBEAM,
    NEMEC
  };

private:
/// Function object returns vector of magnetic field for given point in cylindrical coordinates.
  class Bfield {
    CStconfig *mC;
  public:
    Bfield(CStconfig &magConf) { mC = &magConf; }
/// The operator returns vector of magnetic field for given point in cylindrical coordinates.
    Vector3d operator() ( const Vector3d &cyl ) const {
      return this->mC->getBcyl(cyl);
    }
  };
  CBtrace<Bfield> btrace;

  class sumDat;
  friend class sumDat; // data to speed up summation in cyl2guess

  class __lcms;
  friend class __lcms;
  ///  \brief This class is for internal use only
  class __lcms {
    public:
    void create(const CStconfig * const that) {
      if(!vertex.isOK()) 
        const_cast<CStconfig*>(that)->createLCMS();
    };
    void clear() {
      vertex.clear();
      Bfield.clear();
    };
    CArray2d<Vector3d> vertex;       ///< 2d-Array of Vector3d to hold vertixes of the LCMS
    CArray2d<Vector3d> Bfield;       ///< 2d-Array of Vector3d to hold B-field on the LCMS
    double mRmin, mZmin;      ///< sizes of LCMS in meter
    double mRmax, mZmax;      ///< sizes of LCMS in meter
    double mshift; ///< small angle shift to aviod ray hit into vertex
    double mdPhi;  ///< step along cylindrical angle for tabulating of vertexLCMS
  };

                   
  /// \defgroup epsEffCalculation Internal stuff for eps_eff calculation.
  /// \internal
  //@{
  class __epsilon;
  friend class __epsilon;

  ///  \brief This class is for internal use only
  class __epsilon {
    public:
    void invalidate() {
      Btruncation=0;
    };
    void create(const CStconfig * const that) {
      if(!x.isOK()||Btruncation!=that->truncate()) { // recalculate if truncation has been changed 
        const_cast<CStconfig*>(that)->calculateEpsEff();
      }
    };
    void resize() {
      int n=maxP;
      x.resize(n);
      epsGraz.resize(n);
      epsDkes.resize(n);
      epsVmec.resize(n);
      epsNorm.resize(n);
      epsHat.resize(n);
      GwGraz.resize(n);
      GsGraz.resize(n);
      GvGraz.resize(n);
    };
    void clear() {
      x.clear();
      epsGraz.clear();
      epsDkes.clear();
      epsVmec.clear();
      epsNorm.clear();
      epsHat.clear();
      GwGraz.clear();
      GsGraz.clear();
      GvGraz.clear();
    };
    __epsilon() {
      init();
    }
    ~__epsilon() {
      init();
    }
    void init() {
        Btruncation=0;
        dphi = 0.017453;     // twopi/360 = 0.017453 -> 1degree;     // dphi=0.00698  == 0.4 degree
        doAccuracyTest =true;
        avoidResonances=true;
        useRationalIota=false;
        iotaDenominator=151;
        turns=100;
        numBLevels=2000;
        useBoozerCoord=false; // if true then the properties of boozer coordinates is used in calculations
                              // don't set the following constant to true in normal run.
        maxP = 21;
        xmin = .05;
        xmax = .98;
    };
    double Bmin,Bmax;
    double Btruncation;
    double xmin;
    double xmax;
    double dphi;
    int maxP;
    int iotaDenominator;
    int turns;
    int numBLevels;
    bool doAccuracyTest;
    bool avoidResonances;
    bool useRationalIota;
    bool useBoozerCoord; // if true then use properties of boozer coordinates

    // see 28th EPS, http://epsppd.epfl.ch/Madeira/html/pdf/P5.055.pdf
    CArray1d<double> x;        // x = sqrt(s) is the flux label    
    CArray1d<double> epsGraz;  // using Graz definition of minor radius: dr = dV/S = ds/<|grad(s)|>
    CArray1d<double> epsDkes;  // using Dkes definition of minor radius: r = sqrt( Flux(s)/pi/B_00(s) )
    CArray1d<double> epsVmec;  // using VMEC definition of minor radius: r = a*sqrt(s)
    CArray1d<double> epsNorm;  // eps_eff normalized to the 'classical' particle diffusion flux, see eq.(9)
    CArray1d<double> epsHat;   // eps_eff_hat which corresponds to the Boozer difinition of diffusion coefficient

    CArray1d<double> GwGraz;   // \Gamma_w  eq(21) in Phys.Plasmas 12,112507(2005) http://link.aip.org/link/doi/10.1063/1.2131848
    CArray1d<double> GsGraz;
    CArray1d<double> GvGraz;


    double Graz(double s) { return epsGraz.interp1(sqrt(s), x); }
    double Dkes(double s) { return epsDkes.interp1(sqrt(s), x); }
    double Vmec(double s) { return epsVmec.interp1(sqrt(s), x); }
    double Norm(double s) { return epsNorm.interp1(sqrt(s), x); }
    double Hat (double s) { return epsHat.interp1(sqrt(s), x); }
    double Gw  (double s) { return GwGraz.interp1(sqrt(s), x); }
    double Gs  (double s) { return GsGraz.interp1(sqrt(s), x); }
    double Gv  (double s) { return GvGraz.interp1(sqrt(s), x); }
  };

  ///  \brief This class is for internal use only
  struct __allEps {
    double GwGraz;  // eq(21) in Phys.Plasmas 12,112507(2005)
    double GsGraz;  // eq(26) in Phys.Plasmas 12,112507(2005)
    double GvGraz;  // eq(28) in Phys.Plasmas 12,112507(2005)

    double Graz;  // using Graz definition of minor radius: dr = dV/S = ds/<|grad(s)|>
    double Dkes;  // using Dkes definition of minor radius: r = sqrt( Flux(s)/pi*B_00(s) )
    double Vmec;  // using VMEC definition of minor radius: r = a*sqrt(s)
    double Norm;  // eps_eff normalized to the 'classical' particle diffusion flux, see eq.(9)
    double Hat;   // eps_eff_hat which corresponds to the Boozer difinition of diffusion coefficient
    double gradSmodAverage;   //  <|grad(s)|>
    double vmecFactor;        //  a^2*<|grad(s)|>^2 / (4s)
    double dkesFactor;        //  <|grad(Psi)|>^2 / (4pi*B0*Psi)
    double xi;                //  xi = <|grad(s)|>^2 / (B0^2 <(grad(s)/B)^2 > )
  };

//**********************************************************************************
  class __Fbs;
  friend class __Fbs;

  ///  \brief This class is for internal use only
  class __Fbs {
    public:
    void invalidate() {
      Btruncation=0;
    };
    void create(const CStconfig * const that) {
      if(!x.isOK()||Btruncation!=that->truncate()) { // recalculate if truncation has been changed
        const_cast<CStconfig*>(that)->calculateFbs();
      }
    };
    void resize() {
      int n=maxP;
      x.resize(n);
      fbsGraz.resize(n);
      fbsDkes.resize(n);
      fbsVmec.resize(n);
      lambda_b.resize(n);
      ftrapped.resize(n);
      g2vmec.resize(n);
      g2graz.resize(n);
      g2dkes.resize(n);
      xmm.resize(numMuLevels);
      g2_.resize(n);
      g4_.resize(n,numMuLevels);
      u_.resize(n);   
      uB2_.resize(n); 
      u2B2_.resize(n);
      B2_.resize(n);  
   };
    void clear() {
      x.clear();
      fbsGraz.clear();
      fbsDkes.clear();
      fbsVmec.clear();
      lambda_b.clear();
      ftrapped.clear();
      g2vmec.clear();
      g2graz.clear();
      g2dkes.clear();
      xmm.clear();
      g2_.clear();
      g4_.clear();
      u_.clear();   
      uB2_.clear(); 
      u2B2_.clear();
      B2_.clear();  
      init();
    };
    __Fbs() {
      init();
    }
    void init() {
        Btruncation=0;
        dphi = 0.017453;     // twopi/360 = 0.017453 -> 1degree; 
        doAccuracyTest =true;
        avoidResonances=false;
        useRationalIota=false;
        iotaDenominator=151;
        turns=100;
        xmumax=15;
        numMuLevels=257;
        useBoozerCoord=false; // if true then the properties of boozer coordinates is used in calculations
                              // don't set the following constant to true in normal run.
        maxP = 101;
        xmin = .05;
        xmax = .95;
    };
    double Btruncation;
    double xmin;
    double xmax;
    double xmumax; // mu_max = tanh(xmumax)
    double dphi;
    int maxP;
    int iotaDenominator;
    int turns;
    int numMuLevels;
    bool doAccuracyTest;
    bool avoidResonances;
    bool useRationalIota;
    bool useBoozerCoord; // if true then use properties of boozer coordinates
    // fbsVmec = ftrapped*Gbs*dV/dr_vmec, see  Gbs in N.Nakajima et al, NF, vol.29, 605(1989)
    // see  lambda_b in Nemov, Kalyuzhnyj, Kasilov et al, PPCF, vol.46, 179(2004)
    // fbsVmec = lambda_b*B2Average/(B0*B0)*(dr_vmec/dr_graz);   dr_vmec/dr_graz = <|grad(s)|>* dr_vmec/ds 
    CArray1d<double> x;        // x = sqrt(s) is the flux label    
    CArray1d<double> fbsGraz;  // using Graz definition of minor radius: dr = dV/S = ds/<|grad(s)|>
    CArray1d<double> fbsDkes;  // using Dkes definition of minor radius: r = sqrt( Flux(s)/pi/B_00(s) )
    CArray1d<double> fbsVmec;  // using VMEC definition of minor radius: r = a*sqrt(s)
    CArray1d<double> lambda_b; // fbsGraz**B0*B0/B2Average
    CArray1d<double> ftrapped; // fraction of trapped particles
    CArray1d<double> g2vmec;   // <g2> using VMEC definition of minor radius: r = a*sqrt(s)
    CArray1d<double> g2graz;   // <g2> using Graz ........ 
    CArray1d<double> g2dkes;   // <g2> using Dkes ........ 
    CArray1d<double> g2_;      // <g2>, where g2 is defined as B*grad(g2/b^2) = Bxgrad(s)*grad(1/b^2)  
    CArray2d<double> g4_;      // <g4>, where g4 is defined as B*grad(g4/xi)  = Bxgrad(s)*grad(1/xi) 
    CArray1d<double> xmm;      //  mu = tanh(xmm), where 0<=mu<=1 is the magnetic moment 
    CArray1d<double> u_;       // <u>, where u is defined as B*grad(u) = -Bxgrad(s)*grad(1/b^2)  
    CArray1d<double> uB2_;     // <u*B^2>, where u is defined as B*grad(u) = -Bxgrad(s)*grad(1/b^2)  
    CArray1d<double> u2B2_;    // <u^2*B^2>, where u is defined as B*grad(u) = -Bxgrad(s)*grad(1/b^2)  
    CArray1d<double> B2_;      // <B^2>
    double Graz(double s) { return fbsGraz.interp1(sqrt(s), x); }
    double Dkes(double s) { return fbsDkes.interp1(sqrt(s), x); }
    double Vmec(double s) { return fbsVmec.interp1(sqrt(s), x); }
    double GrazLambda(double s) { return lambda_b.interp1(sqrt(s), x); }
    double ftr(double s) { return ftrapped.interp1(sqrt(s), x); }
    double g2(double s)  { return g2_.interp1(sqrt(s), x); }
    double g2Vmec(double s) { return g2vmec.interp1(sqrt(s), x); }
    double g2Graz(double s) { return g2graz.interp1(sqrt(s), x); }
    double g2Dkes(double s) { return g2dkes.interp1(sqrt(s), x); }
    double u(double s)     { return u_.interp1(sqrt(s), x); }
    double uB2(double s)   { return uB2_.interp1(sqrt(s), x); }
    double u2B2(double s)  { return u2B2_.interp1(sqrt(s), x); }
    double B2(double s)    { return B2_.interp1(sqrt(s), x); }
 };

  ///  \brief This class is for internal use only
  struct __allFbs {
    double Graz;     // Fbs using Graz definition of minor radius: dr = dV/S = ds/<|grad(s)|>
    double Dkes;     // using Dkes definition of minor radius: r = sqrt( Flux(s)/pi*B_00(s) )
    double Vmec;     // using VMEC definition of minor radius: r = a*sqrt(s)
    double lambda_b; // Graz**B0*B0/B2Average
    double ftrapped; // fraction of trapped particles
    double g2vmec;   // g2 using VMEC definition of minor radius: r = a*sqrt(s)
    double g2graz;   // g2 using Graz ........
    double g2dkes;   // g2 using Dkes ........
    double g2_;
    double u_;       // <u>, where u is defined as B*grad(u) = -Bxgrad(s)*grad(1/b^2)  
    double uB2_;     // <u*B^2>, where u is defined as B*grad(u) = -Bxgrad(s)*grad(1/b^2)  
    double u2B2_;    // <u^2*B^2>, where u is defined as B*grad(u) = -Bxgrad(s)*grad(1/b^2)  
    double B2_;      // <B^2>
    CArray1d<double> g4_;
  };
  //@}

//**********************************************************************************

  class __averagedData;
  friend class __averagedData;

  ///  \brief This class is for internal use only
  class __averagedData {
    public:
    void create(const CStconfig * const that) {
      if(!x.isOK()||Btruncation!=that->truncate()) {
        const_cast<CStconfig*>(that)->calculateGradsAvrg();
      }
    };
    void resize(int n) {
      x.resize(n);
      r.resize(n);
      grads.resize(n);
      grads2.resize(n);
      grads2overB2.resize(n);
      rdks.resize(n);
      eB00.resize(n);
    };
    void clear() {
      x.clear();
      r.clear();
      grads.clear();
      grads2.clear();
      grads2overB2.clear();
      rdks.clear();
      eB00.clear();
    };
    double Btruncation;
    CArray1d<double> x;        // x = sqrt(s) is the flux label    
    CArray1d<double> r;        // r = integral (ds/<|grad(s)|>) -- Graz definition   
    CArray1d<double> rdks;     // rdks = sqrt( integral (2x*dx*Flux(1)/(pi*B00)) ) = sqrt(integral(dpsi/(pi*B_00)) ); where psi is the toroidal flux Flux(s)
    CArray1d<double> eB00;     // fabs(Flux(1)/( pi*B00(x) ))
    CArray1d<double> grads;    // <|grad(s)|>/sqrt(s);  <|grad(s)|>/sqrt(s) is proportional to <|grad(r)|>,where r=a*sqrt(s)  
    CArray1d<double> grads2;   // <|grad(s)|^2>/s  
    CArray1d<double> grads2overB2;  // <|grad(s)/B|^2>/s
    double r_graz(double s) { return r.interp1(sqrt(s), x); }
    double r_dks(double s) { return rdks.interp1(sqrt(s), x); }
    double Grads (double s) { return grads.interp1(sqrt(s), x)*sqrt(s); }        // <|grad(s)|>
    double Grads2(double s) { return grads2.interp1(sqrt(s), x)*s; }             // <|grad(s)|^2>
    double Grads2overB2(double s) { return grads2overB2.interp1(sqrt(s), x)*s; } // <|grad(s)/B|^2>
    void set_r_graz() 
    {
      r.rw()[0]=0;
      rdks.rw()[0]=0;
      for(int i=1; i<x.size(); i++) {
        double rg, rdk;
        rInt(x[i-1],x[i], rg, rdk);
        r.rw()[i]   = r[i-1]    + rg;
        rdks.rw()[i]= rdks[i-1] + rdk;
      }
      for(int j=1; j<x.size(); j++) {
        rdks.rw()[j]= sqrt(rdks[j]);
      }
    }
  private:
    double rprime(double x1)  { return 2/grads.interp1(x1, x); }  // @return 2sqrt(s)/<|grad(s)|>
    double rdkprime(double x1)  { return 2*x1*eB00.interp1(x1, x); } // @return 2x*Flux(1)/(pi*B00)
    void rInt(double x1, double x2, double &rg, double &rdk) 
    {
      double sum1=0,sum2=0,dx=0.005;
      double Dx=x2-x1;
      rg = rdk = 0;
      if(Dx==0) return;
      int n = int(Dx/dx);
      if(n<4) n=4;
      n=n+n%2;       //must be even for simpson integration
      dx=Dx/n;
      for(int k=0; k<=n; k++) {  //Simpson integration
        double w = (k==0||k==n)?1:(2*(1+k%2));   // weight for Simpson integration
        sum1 += w*rprime(x1+k*dx);
        sum2 += w*rdkprime(x1+k*dx);
      }
      rg  = sum1*dx/3;
      rdk = sum2*dx/3;
    }
  };

  /// \defgroup vmec2boozerTransrorm Internal stuff for VMEC to Boozer transformation.
  /// \internal
  //@{
// vmec to boozer transformation  
  class __booz;
  friend class __booz;

  ///  \brief This class is for internal use only
  class CComplex {
  public:
    double  re; 
    double  im;
    CComplex() :re(0),im(0) {;}
    CComplex(CComplex const &v ) :re(v.re),im(v.im) {;}
    CComplex(double x,double y=0.) :re(x),im(y) {;}
    double &real() { return  re;}
    double &imag() { return  im;}
    const double &real() const { return  re;}
    const double &imag() const { return  im;}
    /// Simple assignment
    CComplex &operator = ( CComplex const &v ) {
      if( this != &v) {
        re = v.re; im = v.im;
      }
      return( *this );
    }
    /// Returns the conjugate of the complex number
    CComplex conj() {
      return( CComplex(re,-im) );
    }
    /// multiplication assignment
    CComplex operator *= ( CComplex const &v ) {
        double x=re;
        re  = x*v.re-im*v.im;
        im  = x*v.im+im*v.re;
        return( *this );
    }
    /// product of two \b CComplex
    friend CComplex operator * ( CComplex const &v1, CComplex const &v2 ) {
        return CComplex(v1.re*v2.re-v1.im*v2.im,v1.re*v2.im+v1.im*v2.re);
    }
  };

  ///  \brief This class is for internal use only
  class __booz {
    void setMN(int M, int N, int method) {
      Mmax = M;
      Nmax = N;
      if(method==0) { // intergate over equidistant boozer angles
        Npol = Mmax*2;
        Ntor = Nmax*2;
      }
      else {
        Npol = Mmax*4;
        Ntor = Nmax*4;
      }
   }
   public:
    int method; 
    int Ns;
    int Mmax; // # of poloidal mods
    int Nmax; // # of toroidal mods
    int Npol;
    int Ntor;
    CArray1d< CComplex > ctm;
    CArray1d< CComplex > cpn;
    CArray3d<double> Rmn, Zmn, Pmn, Bmn;  // array[s](m,n)
    __booz() {
    }
    __booz(int n,int M, int N, int method) {
      resize(n,M,N,method);
    }
    bool getMethod() {return method==0; } 
    void resize(int n,int M, int N, int method) {
      this->method = method;  // 1-integrate over VMEC angles; 0 - intergate over boozer angles
      setMN(M, N, method); 
      Ns = n;       // number of flux surfaces
      Rmn.resize(n);
      Zmn.resize(n);
      Pmn.resize(n);
      Bmn.resize(n);
      ctm.resize(0,Mmax);
      cpn.resize(-Nmax,Nmax);
    }
    void clear() {
      Rmn.clear();
      Zmn.clear();
      Pmn.clear();
      Bmn.clear();
      ctm.clear();
      cpn.clear();
    }
  };
  //@}

  __lcms lcms;
  __epsilon epsilonEff;
  __Fbs Fbs;
  __averagedData  averagedData;
  __booz bz;


  std::string mHdr1;          // header for bc-file
  std::string mHdr2;          // header for bc-file
  std::string mHdr2a;         // header for bc-file
  std::string mHdr21;         // header for bc-file
  std::string mHdr3;          // header for bc-file
  std::vector<std::string> mComments; ///< To hold comment strings from *.bc file

  class __spline;
  friend class __spline;
  ///  \brief This class is for internal use only
  ///  used for splining of Ip, It, iota, Volume, Vprime,.....(flux surface quantities)
  class __spline : public CArray1d<double> {
    public:
    __spline() :CArray1d<double>() {created=false;};

    bool init(int low,int upper,double val) {
      bool ok = (dynamic_cast<CArray1d<double>*>(this))->resize(low,upper,val); 
      ok     &= _d2.resize(low,upper,val);
      created = false;
      return (ok);
    }
    void clear() {
      (dynamic_cast<CArray1d<double>*>(this))->clear();
      _d2.clear();
      created = false;
    }
    void create(CStconfig *mc) {  // _d2 will be filled
      int ns = mc->mNs;     //     int nnnn = mc->mreff.size()-1;
      mc->splineNRC(ns,mc->mreff.constArray(),this->constArray(),_d2.rw(), 0 ); 
      _d2.rw()[ns-2] = 0; // zero second derivative in order to extrapolate linealy
      created=true;
    }
    // Spline interpolation of Ip, It, iota, Volume, Vprime,.....(flux surface quantities)
    double operator() (double s,const CStconfig *mc) const {
      int i = mc->SearchSindx(s);
    //  int ix = SearchXindx(sqrt(s));
    //  if ( i!=ix ) 
    //    i = ix;
      const CArray1d<double> &mr =  mc->mreff;
      double h = mr[i+1] - mr[i];
      double b = (sqrt(s) - mr[i])/h;
      double a = 1-b;
      double c = a*a*a-a;
      double d = b*b*b-b;
      const double *y = this->constArray();
      const double *y_d2 = _d2.constArray();
      return y[i] + b*(y[i+1]-y[i]) + (c*y_d2[i] + d*y_d2[i+1])*h*h/6;
    }
    // Spline interpolation of Ip, It, iota, Volume, Vprime,.....(flux surface quantities)
    // @return dy/ds
    double prime(double s,const CStconfig *mc) const {
      s = (s<=mc->ms[1])?(mc->ms[1]+0.0001*(mc->ms[2]-mc->ms[1])):s;
      double e2r = 0.5/sqrt(s);
      int i = mc->SearchSindx(s);
      const CArray1d<double> &mr =  mc->mreff;
      double h = mr[i+1] - mr[i];
      double b = (sqrt(s) - mr[i])/h;
      double a = 1-b;
    //double c = a*a*a-a;
    //double d = b*b*b-b;
    //double y = a*y[i] + b*y[i+1] + (c*y_d2[i] + d*y_d2[i+1])*h*h/6;
      const double *y = this->constArray();
      const double *y_d2 = _d2.constArray();
      double yp = (y[i+1]-y[i])/h-( (3*a*a-1)*y_d2[i]-(3*b*b-1)*y_d2[i+1] )*h/6;
      return yp*e2r;
    }

    bool created; // spline coeffs have been created if true
    CArray1d<double> _d2; // spline coefficients
  };

  CArray1d<double> ms;         ///< ms[mNs]--array of normalized toroidal flux
  CArray1d<double> mreff;      ///< mreff[mNs]--array of normalized mreff=sqrt(ms)
  CArray1d<double> msFull;     ///< array[mNs]--normalized toroidal flux, VMEC full mesh
  CArray1d<double> mreffFull;  ///< array[mNs]--array of normalized sqrt(msFull), VMEC full mesh

  __spline msPol;      ///< msPol[mNs]--array of normalized poloidal flux
  __spline miota;      ///< array of iota, length is mNs
  __spline mpres;      ///< array of pressure, length is mNs 
  __spline mIpol;      ///< poloidal current on a period
  __spline mItor;      ///< toroidal current
  __spline mpprime;    ///< Pressure derivative, dp/ds
  __spline mg00;       ///< g00:  dV/ds = |g00(s)*Np|
  __spline mVolume;    ///< Volume[m^3] inside ms[is], array of length mNs

  CArray1d<int> mpM; ///< number of poloidal harmonic after truncation for each flux-surface, array of length mNs
  CArray1d<int> mpN; ///< number of toroidal harmonic after truncation for each flux-surface, array of length mNs
  int mNs;            ///< number of flux surface after adding surfaces near magnetic axis and beyond LCMS
  int mNs0;           ///< number of flux surface stored in the file(original number)
  int mM;             ///< number of poloidal harmonic  0  <= m  <= mM,that was read initially
  int mN;             ///< number of toroidal harmonic -mN <= n  <= mN,that was read initially
  int MminCurrent;
  int NminCurrent;
  int MmaxCurrent;
  int NmaxCurrent;
//harmonics
  CArray1d<double> mCOEF;  ///< Bmn,Rmn,Zmn,Phimn are here, see CStconfig::Bmn(), CStconfig::Rmn()
  int offsetN;    ///<               = 8 -- for accessing Bmn,Rmn,Zmn,Phim
  int offsetM;    ///<        =(2mN+1)*8 -- for accessing Bmn,Rmn,Zmn,Phim
  int offsetS;    ///< =(mM+1)*(2mN+1)*8 -- for accessing Bmn,Rmn,Zmn,Phim

//mCOEF[ns][m+1][2*n+1][8];
//mCOEF[I][M][N][K];
//mCOEF[i][m][n][k] = mCOEF[k+ n*K       + m*N*K     + i*M*N*K  ]
//mCOEF[i][m][n][k] = mCOEF[k+ n*offsetN + m*offsetM + i*offsetS]

  double mFlux;       ///< max Flux [Tm^2]
  double mPolFlux;    ///< max poloidal flux [Tm^2]
  double ma;          ///< minor radius [m]
  double mR0;         ///< major radius [m]
  double mZ0;         ///< Z-coord of mag. axis for tokamak
  double mB0;         ///< saved initial B0, min(B) on axis
  double mB0max;      ///< saved initial B0max, max(B) on axis
  double mBnorm;      ///< scale factor for mFlux,mIpol,mItor, Bmn, BfieldLCMS
  double mmodBtol;    ///< relative tolerance in |B|, this parameter is used to reduce number of harmonics, see CStconfig::truncate(double m_modBtol)
  double mCrossArea;  ///< average cross-section area
  double mB2VolAvrg;  ///< volume averaged B^2
  double sLast;      ///<  = ms[mLCMSindex] -- (last s in bc-file)
  double s1st;       ///<  = ms[m1stSindex] -- (first s in bc-file)
  int mLCMSindex;    ///<  index of the LCMS in ms-array
  int m1stSindex;    ///< number of surfaces that will be added in the vicinity of s=0 or 1st index in ms to read surfaces, see CStconfig::loadascii()
  int mNsAdd;        ///< number of additional surfaces that will be added after the LCMS
  int mDirTheta;     ///< direction of theta 
  int mDirPhi;       ///< direction of phi w.r.t cylindrical angle
  int mSignJac_s;    ///< sign of jacobian in (s,theta,phi)-coordinates
  bool mNewFormat;   ///< if true then new format of ascii bc-file is used
  fileFormat formatOfLoadedFile;


  CArray1d<Vector3d> mBorderMagCoord; ///< array of border points in magnetic coordinates
  CArray1d<Vector3d> mBorder;       ///< array of border points in cylindrical coordinates
  double mBorderRmin, mBorderZmin; ///< sizes of Border in meters
  double mBorderRmax, mBorderZmax; ///< sizes of Border in meters

  Vector3d  Smatrix[3];   ///< Susceptance matrix: S_ij = SMatrix[i][j]; i=1,2; j=1,2 in (torFlux,th,ph)-coordinates

  void * workPointer;

protected:
  bool mOK;              ///< result of the last operation
  std::string mfname;         ///< name of the file, that was read
  double pi;             ///< = 3.14159265358979
  double twopi;
  double twopiinv;
  double degree;         ///< = pi/180;
  double mu0;            ///< = pi*4e-7 [N/A^2 = henry/m]
  int mNp;            ///< number of periods
  double mPeriod;     ///< twopi/mNp
  double mepsA, mepsR;///< absolute and relative accuracy of coordinate transformation

//public:  // for testing
  /// \defgroup NewtonMethod Internal stuff for Newton's method.
  /// \internal
  //@{
  int mNGuessCalls;  
  int mNjac;          ///< number of jacobian calculation
  int mMsum;          ///< number of poloiodal mod 
  int mNsum;          ///< number of toroidal mod
  int mNcyl2mag;      ///< number of calling of cyl2mag-function(Newton's iteration)
  int mNmag2cyl;      ///< number of calling of mag2cyl-function
  int mMsum2;
  int mNsum2;
  int mRoughCalculation; ///< used in cyl2guess for
  ///  \brief  The internal structure holds things after last 
  ///   coordinate transformation/NJacobian calculation.
  struct ctResults {
    ctResults() {
      is = 0;
      s = modB = jac  = jacMixProd = 0;     
      jacVmec = jacBoozer = p=ps=pt=pp = lambda = lambdas = lambdat = lambdap = 0; 
      axisLast = xyzLast = cylLast = mBnormLast = 1e128; 
    }
    void invalidate() { axisLast = xyzLast = cylLast = mBnormLast = 1e128;}
    bool isCoordValid() { return (cylLast[0]>=1e127)?false:true; }
    bool isBnormValid(double mBnorma) { return (mBnormLast==mBnorma)?true:false; }
    bool isAxisValid(const Vector3d &cyl) {
      if(!isCoordValid()) return false;
      return (fabs(axisLast[1]-cyl[1])<1e-8)?true:false;
    }
    bool isCylValid(const Vector3d &cyl,double mBnorma) {
      if(!isCoordValid()) return false;
      if(!isBnormValid(mBnorma)) return false;
      return (cyl==cylLast)?true:false;
//      return (diffCyl2(cyl)<1e-20)?true:false;  // may be this will be better?
    }
    bool isXyzValid(const Vector3d &xyz,double mBnorma) {
      if(!isCoordValid()) return false;
      if(!isBnormValid(mBnorma)) return false;
      return (diffXyz2(xyz)<1e-20)?true:false;
    }
    double diffCyl2(const Vector3d &c) {
      // dR - distance between points: dR=(r-cylLast)
      // @return dR^2
      double dfi = c.y()-cylLast.y();
      double dR2 = cylLast.x()*cylLast.x() + c.x()*c.x() 
                - 2*cylLast.x()*c.x()*::cos(dfi) 
                + (c.z()-cylLast.z())*(c.z()-cylLast.z());
      return dR2;
      ////Vector3d dR(cylLast[0]-c[0]*::cos(dfi),c[0]*::sin(dfi),cylLast[2]-c[2]);
      ////return dR.abs2();
    }
    double diffXyz2(const Vector3d &r) {
      return (xyzLast-r).abs2();
    }
    void saveLastJacobian(double Bnorm) {
      mBnormLast = Bnorm;
      magLast = Jmatr[4];
      cylLast = Jmatr[0];
      xyzLast = cylLast.toCartesian();
      Jlast[0]=Jmatr[0]; 
      Jlast[1]=Jmatr[1]; 
      Jlast[2]=Jmatr[2]; 
      Jlast[3]=Jmatr[3]; 
      Jlast[4]=Jmatr[4];
    }
    int       is;
    double    s;
    double    modB;     // |B|
    double    jac;       // Jacobian in coordinates (torFlux,theta,phi), for example in boozer coordinates the  jac ~ mu0/(4*pi^2)*Ipol/B^2
    double    jacMixProd;// Jacobian in coordinates (torFlux,theta,phi) from mixed product of basis vectors: R/mFlux*(Rs*(Rt^Rp))
    double    jacVmec;    // Vmec Jacobian   in  (s,theta,phi)
    double    jacBoozer;  // Boozer Jacobian in  (s,theta,phi)
    double    p;        // p function for transformation vmec to boozer angle: thetaB = theta + lambda + iota*p, phiB = phi + p
    double    ps;       // d(p)/ds
    double    pt;       // d(p)/dtheta
    double    pp;       // d(p)/dphi
    double    lambda;   // Vmec stream function lambda
    double    lambdas;  // d(lambda)/ds
    double    lambdat;  // d(lambda)/dtheta
    double    lambdap;  // d(lambda)/dphi
    double    mBnormLast; // last norm. factor for |B| for which coord.trans. has been made by cyl2mag() or xyz2mag()
    Vector3d  axisLast;  // last axis position at which coord.trans. has been made by xyz2mag
    Vector3d  xyzLast;  // last cartesian coordinates for which coord.trans. has been made by xyz2mag
    Vector3d  cylLast;  // last cylindrical coordinates for which coord.trans. has been made by cyl2mag
    Vector3d  magLast;  // last magnetic coordinates from cyl2mag
    Vector3d  Jlast[5]; // last Coordinates and jacobian matrix are here
    Vector3d  Bxyz;     // B field in cartesian coordinates 
    Vector3d  Bcyl;     // B field in cyl. coordinates 
    Vector3d  Bcontra;  // B field contravariant components (0,B^theta,B^phi) 
    Vector3d  Bcova;    // B field covariant components (B_s,B_theta,B_phi) in (s,theta,phi) coordinates 
    Vector3d  gradTheta;  // grad(theta)  = (Xp^Xf)/jac in cyl. coord, theta is the poloidal angle
    Vector3d  gradPhi;    // grad(phi)    = (Xf^Xt)/jac in cyl. coord, phi is the toroidal angle
    Vector3d  gradtorFlux;// grad(torFlux)= (Xt^Xp)/jac in cyl. coord, Psi_tor is the toroidal flux
    Vector3d  gradS;      // grad(s)    in cyl. coord
    Vector3d  gradReff;   // grad(reff) in cyl. coord

    Vector3d  Jmatr[5]; // Coordinates and jacobian matrix are here, see CStconfig::NJacobian()
                        //Jmatr[0] = F  == (R,fi,Z)
                        //Jmatr[1] = Fs == d(R,fi,Z)/ds   -- partial derivative on s
                        //Jmatr[2] = Ft == d(R,fi,Z)/dtheta
                        //Jmatr[3] = Fp == d(R,fi,Z)/dphi
                        //Jmatr[4] = u  == (s,theta,phi) coordinates
    Vector3d  gsTensor[3];  // Metric tensor g_ij  in (s,theta,phi) coordinates
    //////Vector3d  gTensor[3];   // Metric tensor g_ij  in (torFlux,theta,phi) coordinates
    //////double    detOfgTensor; // det(g_ij) in (torFluxflux,theta,phi) coordinates
  };
  struct ctStack;
  friend struct ctStack;
  /// \brief The internal fast stack for ctResults.
  struct ctStack : std::vector<ctResults> {
    size_t idx;
    ctResults a;
    ctStack() { idx=0; push_back(a); }
    const ctResults & operator() () const { return at(idx); }
    ctResults & operator() () { return top(); }
    ctResults & top() {
      while(size()<=idx) push_back(a);
      //test if(idx>1) std::cerr<<"stack:"<<idx<<std::endl;
      return at(idx);
    };
    /// \brief The internal class increments an index at given address, the destructor restores the index. 
    /// \internal
    /// \ingroup NewtonMethod
    struct saveState {
      ctStack* pct;
      saveState(ctStack &ct) {
        pct = &ct;
        ct.idx++;  // increment stack counter
        ct.top();  // add empty elements to stack if needed
      }
      ~saveState() {
        pct->top().invalidate();
        pct->idx--; // restore stack index
      }
    };
  };
  
  ctStack ct;

#if 0
  std::stack<ctResults> ctStack;  // stack for saving current transformation results
  ctResults res;
  /// for internal use
  void saveNJacData() {
    ctStack.push(res);
  }
  /// for internal use
  void restoreNJacData() {
    if (!ctStack.empty())
      res = ctStack.top();
  }
#endif

  //@}

public:
  /// The method creates a clone of CStconfig using operator new.
  /// The caller is responsible for deleting the clone.
  CStconfig * clone() {
    CStconfig *mc = new CStconfig;
    *mc = *this; 
    return mc;
  }
  /// The method returns the version of MConf in format year.month
  static const char *getVersion();
  /// The constructor creates an empty instance of CStconfig.
  CStconfig() {init();}
  /// The constructor creates an instance of CStconfig loading file.
  CStconfig(const char * filename) {init(); loadfile(filename,UNDEFINED );}
  /// The destructor is a virtual method.
  virtual ~CStconfig() {
    clear();}
  /// The method frees all internal arrays and makes this object empty.
  virtual void clear();
  /// The method reload the current magnetic configuration to the right oriented (s,theta,phi) coordinate system.
  int setRightCoordinateSystem();
  /// The method reload the current magnetic configuration to the left oriented (s,theta,phi) coordinate system.
  int setLeftCoordinateSystem();
  /// The method reload the current magnetic configuration to opposite oriented coordinate system.
  int reverseCoordinateSystem();
  /// The method reload the current magnetic configuration to opposite oriented 
  /// coordinate system by fliping z-coordinate.
  int flipZcoordinate();
  /// \name Load & save functions
  //@{
  /// The method loads the magnetic configuration stored in the file \b filename.
  /// The method recognizes the format of the file by analyzing its contents;
  /// bc-, bc-binary-, LHD- and VMEC(wout, version>=6.20)-format are supported.
  /// For bc-format please see \ref W7X_Format "W7-X format" or ask
  /// <a class="el" href="mailto:joachim.geiger_at_ipp.mpg.de">J.Geiger.</a>
  /// The binary-format file is produced by the CStconfig::write() method.
  /// Usually you may want to use this format for faster loading of the magnetic configuration;
  /// and probably you don't need to know the exact structure of the binary-format file.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// \code
  ///  #include "CStconfig.h"
  ///  int main() {
  ///    CStconfig mConf;
  ///    mConf.load("w7x-sc1.bc");
  ///    if(!mConf.isOK()) exit(1);   // exit if not OK
  ///    if(!mConf.load("w7x-sc1.bin4")) exit(1); // exit if not OK
  ///    if(!mConf.load("w7x-sc1.bin8")) exit(1); // exit if not OK
  ///  }
  ///  \endcode
  ///
  virtual bool load(const char * filename)  { return loadfile(filename,UNDEFINED );}
  /// The method loads the magnetic configuration stored in the file \b filename.
  /// @param filename is the name of a file to read.
  /// @param filePosition is the file starting position from which the data will be read,
  ///  the function fseek(fp,filePosition,SEEK_SET) is used;
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// \note file must be in text format. 
  bool load(const char * filename, long filePosition)  { return loadfile(filename,UNDEFINED,filePosition);}
  /// The method loads the \ref GEQDSK "G EQDSK file" 
  /// and transforms it to the \ref TokamakSymmetryFluxCoordinates "tokamak-symmetry flux coordinates".
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool loadEFIT(const char * filename)  { return loadfile(filename,EFIT);}
  /// The method loads the magnetic configuration with Boozer coordinates 
  /// from the file in LHD-format. 
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool loadLHD(const char * filename)  { return loadfile(filename,LHD );}
  /// The method loads the \ref VmecCoordinates "VMEC wout-file".
  /// No transformation to Boozer representation is performed.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool loadVMEC(const char * filename)  { return loadfile(filename,VMEC);}
  /// The method loads the \ref VmecCoordinates "VMEC netCDF wout-file".
  /// No transformation to Boozer representation is performed.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
#ifdef NETCDF
  bool loadVMECnetcdf(const char * filename)  { return loadfile(filename,VMECNETCDF);}
#endif
  /// The method loads the \ref VmecCoordinates resulting from FP (function parametrization).
  /// No transformation to Boozer representation is performed.
  /// @param filename is the name of a file.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool loadVMECFP(const char * filename)  { return loadfile(filename,W7X_FMT_VMEC_FP);}
private:
  /// The method loads the "NEMEC wout-file".
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool loadNEMEC(const char * filename)  { return loadfile(filename,NEMEC);}
public:
  /// The method loads Boozer coordinates from Magnetic Configuration Data Base.
  /// @param mcdb -- address of the MCDB control structure, see <a href="html/MCDBAPI.html">MCDB documentation</a>
  /// @param equilibrium_id -- the value of \b id field from the \b equilibrium table.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool loadFromDB(const void *mcdb, int equilibrium_id);
  /// The method deletes magnetic configuration from database.
  /// @param mcdb -- address of the MCDB control structure, see <a href="html/MCDBAPI.html">MCDB documentation</a>
  /// @param equilibrium_id -- the value of \b id field from the \b equilibrium table.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool deleteFromDB(const void *mcdb, int equilibrium_id) const;
  /// The method inserts magnetic configuration Data into database.
  /// @param mcdb -- address of the MCDB control structure, see <a href="html/MCDBAPI.html">MCDB documentation</a>
  /// @param confname -- configuration name.
  /// @param machine_name -- machine name (for example W7-X).
  /// @param Nsmax -- the number of flux surfaces to save or zero if original
  ///  number of surfaces is needed to save.
  /// @param magAxis -- if true then store magnetic axis.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool writeIntoDB(const void *mcdb, const char* confname,const char* machine_name,int Nsmax,bool magAxis) const;
  /// Store current magnetic configuration into file.
  /// The format of the file is derived from the file extension.
  /// The following file extention are possible:
  ///  \li bc    -- write("w.bc");   //save into ascii file, see \ref W7X_Format "W7-X format" 
  ///  \li bin4  -- write("w.bin4"); //save into binary file using float numbers
  ///  \li bin8  -- write("w.bin8"); //save into binary file using double numbers
  ///
  /// @param fname -- filename
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// \attention Current truncation level and value of magnetic field are used.
  bool write(const char * fname) const;
  /// The methods transform the current magnetic configuration to new volume or to new mesh on s.
  /// @param fname is the filename under which the configuration is saved.
  /// @param nsmax is the number of surfaces to store, if zero then the original flux surfaces are stored,
  ///        otherwise store the surfaces using a new mesh which is equidistant on sqrt(s).
  /// @param newVolume is the new volume if not zero.
  /// @param newLastSurface if not zero then this is the value of new last-closed-flux-surface 
  ///   in inits of current s, this parameter is used to reduce volume.
  /// @param varNumberMods if true then variable number of harmonics for each flux surface are stored.
  /// @param magAxis if true then store magnetic axis.
  /// @param circular is applicable only for tokamak, 
  ///   if true then only cosine Fourier series for R coordinate and sine 
  ///   Fourier series for Z coordinate are stored.
  /// @param append if true then append to the file \b fname .
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// \attention Current truncation level and value of magnetic field are used;
  bool writeasciiReduced(const char * fname,int nsmax=0,double newVolume=0,double newLastSurface=1,
                          bool varNumberMods=true,bool magAxis=false,bool circular=false,bool append=false) const;
  /// Write |B|-spectrum on the  surface \b s into file \b fname.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// @attention Current truncation level and value of magnetic field are used.
  bool writeasciiBmod(const char * fname,double s) const;
  /// The method returns \b true if file is loaded without error.
  bool isOK() const {return mOK;}
  bool isValid() const {return mOK;}
  /// The method returns the name of the file loaded.
  const char  *fname() const {return mfname.c_str(); }
  //@}
  /// \name Accuracy, B value
  //@{
  /// Set the value of magnetic field on magnetic axis at given cylindrical angle.
  /// @param B0 is the value[Tesla] of magnetic field to set
  /// @param ficyl is the value[radians] of cylindrical angle
  virtual void setB0(double B0,double ficyl);
  /// Set minimum magnetic field on magnetic axis.
  /// @param B0 is the value[Tesla] of magnetic field to set
  virtual void setB0(double B0);
  /// Set the B<SUB>00</SUB> value of magnetic field on magnetic axis.
  /// @param B00 is the value[Tesla] of magnetic field to set
  virtual void setB00(double B00);
  /// Set normalization factor for value of the magnetic field
  /// @param multiplier is the scaling factor
  virtual void scaleB(double multiplier);
  /// Set the I<SUB>pol</SUB> value of the coil currents at LCMS.
  /// @param I is the value[Ampers] of the current to set
  virtual void setIcoil(double I);
  /// @return the absolute value scale factor for B.
  virtual double getBscaleFactor();
  /// Get the minimum value of the magnetic field on axis
  /// Get the I<SUB>pol</SUB> value of the coil currents at LCMS.
  /// @return the absolute value of I<SUB>pol</SUB> in Ampers.
  double getIcoil();

#if 0
  /// Set average toroidal magnetic field B<SUB>fi</SUB>(cylindrical component) on magnetic axis.
  /// @param Bfi0 is the value[Tesla] of magnetic field to set
  virtual void setAverageBphi0(double Bfi0);
  /// Get the average toroidal magnetic field B<SUB>fi</SUB>(cylindrical component) on magnetic axis.
  double getAverageBphi0();
#endif

  /// Restore initial value of the magnetic field(which is stored in a bc-file)
  virtual void restoreB0(); // restore initial value of the magnetic field on axis
  /// Get the B<SUB>00</SUB> value of the magnetic field on axis.
  /// @return the B<SUB>00</SUB> value of the magnetic field on magnetic axis.
  double getB00();
  /// Get the minimum value of the magnetic field on axis
  /// @return the minimum value of the magnetic field on magnetic axis
  double getB0();
  /// Get the maximum value of the magnetic field on axis
  /// @return the maximum value of the magnetic field on magnetic axis
  double getB0max();
  /// Get the value of the magnetic field on axis at cylindrical angle ficyl
  /// @param ficyl is the value[radians] of cylindrical angle
  /// @return the value of the magnetic field on magnetic axis at cylindrical angle \e ficyl
  double getB0(double ficyl);
  /// The method returns the direction of B with respect
  /// to cylindrical angle of the right handed cylindrical coordinate system.
  int getBdirection();
  /// The method sets the direction of B with respect
  /// to cylindrical angle of the right handed coordinate system.
  /// @param sign is the positive or negative integer number.
  void setBdirection(int sign);
  /// The method sets the original direction of B with respect
  /// to cylindrical angle of the right handed coordinate system.
  /// The original direction is the direction of B stored in the equilibrium file. 
  void restoreBdirection();
  /// Set truncation level to reduce number of harmonics used in summation.
  /// For each flux surface this method calculates M and N
  /// such that B<SUB>mn</SUB>/B<SUB>00</SUB> > level for m<M and |n|<N.
  /// Then only harmonics with m<M and |n|<N are used in summation.
  /// @param level is the new truncation level.
  /// The method returns the value of the old truncation level
  double truncate(double level);
  /// The method returns the current value of truncation level
  double truncate() const {return mmodBtol;}
  /// Set accuracy of coordinate transformation.
  /// @param epsA is the absolute accuracy[m] of the transformation from cylindrycal to magnetic coordinates.
  /// @param epsR is the relative accuracy of the transformation from cylindrycal to magnetic coordinates.
  /// \note the default accuracy \e epsA is 10<sup>-5</sup>m if this function is not used.
  void setAccuracy(double epsA, double epsR=1e-5);
  //@}
  /// \name Info
  //@{
  /// The method returns the tokamak magnetic axis position in [m]
  double Z0() const {return mZ0;}
  /// The method returns the major radius in [m]
  double R0() const {return mR0;}
  /// The method returns the minor radius [m]
  double rminor() const {return r(1.);}
  /// The method returns the number of poloidal harmonics
  int    Mpol() const  {return mM;}
  /// The method returns the number of toroidal harmonics
  int    Ntor() const {return mN;}
  /// The method returns the number of poloidal harmonics after the truncation level was set.
  /// \sa truncate()
  int    MpolMinTruncated() const  {return MminCurrent;}
  /// The method returns the number of poloidal harmonics after the truncation level was set.
  /// \sa truncate()
  int    MpolMaxTruncated() const  {return MmaxCurrent;}
  /// The method returns the number of toroidal harmonics after the truncation level was set.
  /// \sa truncate()
  int    NtorMinTruncated() const {return NminCurrent;}
  /// The method returns the number of toroidal harmonics after the truncation level was set.
  /// \sa truncate()
  int    NtorMaxTruncated() const {return NmaxCurrent;}
  /// The method returns the number of periods
  int    nPeriods() const {return mNp;}
  /// The method returns the number of flux surfaces used.
  /// The method returns the size of internal arrays CStconfig::ms,CStconfig::miota, CStconfig::mreff,
  /// CStconfig::mItor, CStconfig::mIpol, CStconfig::mpprime, see <b>Member Data Documentation</b> here.
  int    nSurfaces() const {return mNs;}
  /// The method returns the number of flux surfaces stored in the magnetic configuration file.
  int    nSurfacesInitial() const {return mNs0;}
  /// The method returns \e true if the opened file has tokamak symmetry.
  bool   isTokamak() const {return tokamakConfig;}
  /// The method returns \e true if the opened file is a Tokamak Symmetry Flux Coordinates file.
  bool   isTokamakSymmetryFluxCoordinates() const {return tokamakSymmetryFluxCoordinates;}
  /// The method returns \e true if the opened file is a VMEC file.
  bool   isVmecCoordinates() const {return vmecCoordinates;}
  /// The method returns \e true if the opened file is a VMEC file 
  /// from  the Function Parametrization procedure.
  bool   isVmecCoordinatesFromFP() const {return fpCoordinates;}
  /// The method returns  VMEC completion code.
  int    vmecReturnCode()    const {return vmecErrorCode;}
  //@}
  /// \name Surface quantity
  //@{
  /// The method returns the effective radius r=a*sqrt(s), 
  /// where a is the plasma minor radius defined as sqrt(averageCrossSectionArea(s=1)/pi),
  /// this definition of a is used by VMEC code.
  /// @param s is the normalized toroidal flux.
  /// @return effective radius in [m]
  /// \sa averageCrossSectionArea()
  double r     (double s) const;  
  /// The method returns the effective radius defined as sqrt(Flux(s)/B<sub>00</sub>(s)/pi).
  /// This definition \f$r_{eff}=\sqrt{\psi(s) / \pi B_{00}(s)}\f$
  /// (where \f$\psi\f$ is the toroidal flux)  is used in DKES code.  
  /// @return effective radius in [m]
  double rdkes(double s) const;
  /// @return effective radius in [m]
  double rdkes2(double s);
  /// The method returns the minor radius dr=ds/\<|grad(s)|\> 
  /// (Graz definition of the radius \f$dr=ds/\left\langle\left|\nabla s\right|\right\rangle\f$).
  double rgraz(double s);
  /// The method returns  dr_graz/dr, 
  /// where r=a*sqrt(s) is the VMEC definition of minor radius.
  double drgraz_dr(double s); // dr_graz/dr_vmec 
  /// The method returns  dr_dkes/dr, 
  /// where r=a*sqrt(s) is the VMEC definition of minor radius.
  double drdkes_dr(double s); 
  double drdkes2_dr(double s); 
  /// The method returns the effective helical ripple for \f$1/\nu\f$ transport, 
  /// using Graz definition of the minor radius dr=ds/\<|grad(s)|\>, 
  /// where s is the normalized toroidal flux.
  /// @param s is the normalized toroidal flux
  /// \sa V.V.Nemov, S.V.Kasilov, W. Kernbichler, and M.F.Heyn, Physics of Plasmas, vol. 6, 4622(1999),
  ///   http://link.aip.org/link/PHPAEN/v6/i12/p4622/s1 
  /// \sa 28th EPS, http://epsppd.epfl.ch/Madeira/html/pdf/P5.055.pdf , EGA Vol.25A (2001) 1985-1988
  /// \note  At first call the method creates a look-up table using multithreading technique to spread 
  /// callculation among several CPU cores; this table is used for interpolating over the flux surface label.
  /// \n Use the following methods to set parameters:
  /// epsEffSetXeffParam(), epsEffSetTracingParam(), epsEffSetIotaParam(), epsEffSetMagnMomentParam()
  double epsEffGraz(double s); 
  /// The method returns the effective helical ripple for \f$1/\nu\f$ transport, 
  /// using DKES definition of minor radius r=sqrt(Flux(s)/(pi*B_00(s))), 
  /// where s is the normalized toroidal flux.  
  /// @param s is the normalized toroidal flux
  /// \sa V.V.Nemov, S.V.Kasilov, W. Kernbichler, and M.F.Heyn, Physics of Plasmas, vol. 6, 4622(1999),
  ///   http://link.aip.org/link/PHPAEN/v6/i12/p4622/s1 
  /// \sa 28th EPS, http://epsppd.epfl.ch/Madeira/html/pdf/P5.055.pdf , EGA Vol.25A (2001) 1985-1988
  /// \sa Note to epsEffGraz()
  double epsEffDkes(double s); 
  /// The method returns the effective helical ripple for \f$1/\nu\f$ transport,
  /// using VMEC definition of minor radius r=a*sqrt(s), 
  /// where s is the normalized toroidal flux.
  /// epsEffVmec^1.5 = epsEffGraz^1.5 * \<|grad(s)|\>^2 * (r/(2s))^2
  /// @param s is the normalized toroidal flux.
  /// \sa V.V.Nemov, S.V.Kasilov, W. Kernbichler, and M.F.Heyn, Physics of Plasmas, vol. 6, 4622(1999),
  ///   http://link.aip.org/link/PHPAEN/v6/i12/p4622/s1 
  /// \sa 28th EPS, http://epsppd.epfl.ch/Madeira/html/pdf/P5.055.pdf , EGA Vol.25A (2001) 1985-1988
  /// \sa Note to epsEffGraz()
  double epsEffVmec(double s); 
  /// The method returns the normalized effective helical ripple for \f$1/\nu\f$ transport,
  /// using normalized diffusion coefficient.
  /// \sa 28th EPS, http://epsppd.epfl.ch/Madeira/html/pdf/P5.055.pdf , 
  ///   EGA Vol.25A (2001) 1985-1988
  /// @param s is the normalized toroidal flux
  /// \sa Note to epsEffGraz()
  double epsEffNorm(double s); 
  /// The method returns the effective helical ripple for \f$1/\nu\f$ transport,
  /// using definition (2), (5) in http://epsppd.epfl.ch/Madeira/html/pdf/P5.055.pdf ,
  ///   28th EPS, EGA Vol.25A (2001) 1985-1988
  /// @param s is the normalized toroidal flux
  /// \sa Note to epsEffGraz()
  double epsEffHat (double s);

  /// eq(21) in Physics of Plasmas 12,112507(2005), http://link.aip.org/link/doi/10.1063/1.2131848
  /// \sa Note to epsEffGraz()
  double Gw(double s);
  /// eq(26) in Physics of Plasmas 12,112507(2005), http://link.aip.org/link/doi/10.1063/1.2131848
  /// \sa Note to epsEffGraz()
  double Gs(double s);
  /// The methods return the bounce-averade gradB drift velocity of trapped particles averaged over flux surface.
  /// See eq(28) in Physics of Plasmas 12,112507(2005), http://link.aip.org/link/doi/10.1063/1.2131848
  /// \sa Note to epsEffGraz()
  double Gv(double s);


  /// The method sets parameters for calculating eps_eff and Gw, Gs, Gv.
  /// @param xmin is the minimal value of sqrt(s), where s is the normalized toroidal flux,
  /// @param xmax is the maximal value of sqrt(s),
  /// @param size is the number of points for creating look-up table
  /// \sa Note to epsEffGraz()
 void epsEffSetXeffParam   (double xmin=.05, double xmax=.95, int size=51);   
  /// The method sets parameters for calculating the eps_eff and Gw, Gs, Gv.
  /// @param turns is the maximal number of toroidal turns for a magnetic field line following,
  /// @param doAccuracyTest -- if true then integrating along magnetic field line is stopped 
  ///                          when relative accuracy is less then 1% 
  /// @param dphi is the toroidal angle step in radian for the magnetic field line following.
  /// \sa Note to epsEffGraz()
  void epsEffSetTracingParam(int turns=100, bool doAccuracyTest = true, double dphi = 0.017453);  // twopi/360 == 0.017453
  /// The method sets parameters for calculating the eps_eff and Gw, Gs, Gv.
  /// @param avoidResonances -- if true then change of flux surface label to avoid 
  ///                           iota resonances is attempted,
  /// @param useRationalIota -- if true then rational iota aproximation is used, 
  /// @param iotaDenominator is the denominator guess for the rational iota aproximation.
  /// \sa Note to epsEffGraz()
  void epsEffSetIotaParam   (bool avoidResonances=false, bool useRationalIota=false, int iotaDenominator=150); 
  /// The method sets parameters for calculating the eps_eff and Gw, Gs, Gv.
  /// @param nBLevels is the number of integration points over magnetic moment, 
  ///   see   V.V.Nemov,S.V.Kasilov,W.Kernbichler,M.F.Heyn, Physics of Plasmas, vol.6,4622(1999),
  ///  http://link.aip.org/link/doi/10.1063/1.873749
  /// \sa Note to epsEffGraz()
  void epsEffSetMagnMomentParam (int nBLevels=513); 
  /// The method calculates and stores the eps_eff in the look-up table.
  /// \sa Note to epsEffGraz()
  void epsEffCreate(); 
  /// The method clears the look-up table with the eps_eff  and Gw, Gs, Gv.
  void epsEffInvalidate(); 
  

  /// The method sets parameters for calculating the bootstrap current geometric factor.
  /// @param smin is the minimal value of s, where s is the normalized toroidal flux,
  /// @param smax is the maximal value of s,
  /// @param size is the number of points for creating look-up table
  void FbsSetSlabelParam(double smin=0.0025, double smax=0.9, int size=101); 
  /// The method sets parameters for calculating the bootstrap current geometric factor.
  /// @param xmin is the minimal value of sqrt(s), where s is the normalized toroidal flux,
  /// @param xmax is the maximal value of sqrt(s),
  /// @param size is the number of points for creating look-up table
  void FbsSetXeffParam   (double xmin=.05, double xmax=.95, int size=101);   
  /// The method sets parameters for calculating the bootstrap current geometric factor.
  /// @param turns is the maximal number of toroidal turns for a magnetic field line following,
  /// @param doAccuracyTest -- if true then integrating along magnetic field line is stopped 
  ///                          when relative accuracy is less then 1% 
  /// @param dphi is the toroidal angle step in radian for the magnetic field line following.
  void FbsSetTracingParam(int turns=100, bool doAccuracyTest=true, double dphi=0.017453);  // twopi/360 == 0.017453
  /// The method sets parameters for calculating the bootstrap current geometric factor.
  /// @param avoidResonances -- if true then change of flux surface label value to avoid 
  ///                           iota resonances is attempted,
  /// @param useRationalIota -- if true then rational iota aproximation is used, 
  /// @param iotaDenominator is the denominator guess for the rational iota aproximation.
  void FbsSetIotaParam   (bool avoidResonances=false, bool useRationalIota=false, int iotaDenominator=150); 
  /// The method sets parameters for calculating the bootstrap current geometric factor.
  /// @param lambdaLevels is the number of integration points over magnetic moment,
  /// @param yMax is the new integration limit after variable change from the magnetic moment
  ///    lambda to y with lambda = tanh(y),
  ///    see N.Nakajima et al, Nuclear Fusion,29,605(1989), 
  ///      http://dx.doi.org/10.1088/0029-5515/29/4/006
  void FbsSetMagnMomentParam (int lambdaLevels=257, double yMax=15); 


  /// The method calculates and stores the bootstrap current geometric factor in the look-up table.
  void FbsCreate();
  /// The method clears the look-up table with the bootstrap current geometric factors.
  void FbsInvalidate(); 
  /// The method returns the bootstrap current geometric factor for \f$1/\nu\f$ transport, 
  /// using VMEC definition of minor radius r=a*sqrt(s), 
  /// where s is the normalized toroidal flux.  
  /// @param s is the normalized toroidal flux
  /// \sa P.Helander et al,Phys.Plasmas 18,092505(2011), http://link.aip.org/link/doi/10.1063/1.3633940
  /// \note  At first call the method creates a look-up table using multithreading technique to spread 
  /// calculation among several CPU cores; this table is used for interpolating over the flux surface label.
  /// \n Use the following methods to set parameters:
  /// FbsSetXeffParam(), FbsSetTracingParam(), FbsSetIotaParam(), FbsSetMagnMomentParam()
  double FbsVmec(double s); 
  /// The method returns the bootstrap current geometric factor for \f$1/\nu\f$ transport, 
  /// using DKES definition of minor radius r=sqrt(Flux(s)/(pi*B_00(s))), 
  /// where s is the normalized toroidal flux.  
  /// @param s is the normalized toroidal flux
  /// \sa Note to FbsVmec()
  double FbsDkes(double s);
  /// The method returns the bootstrap current geometric factor for \f$1/\nu\f$ transport, 
  /// using Graz definition of the minor radius dr=ds/\<|grad(s)|\>, 
  /// where s is the normalized toroidal flux.  
  /// @param s is the normalized toroidal flux
  /// \sa Note to FbsVmec()
  double FbsGraz(double s);  
  /// The method returns the bootstrap current geometric factor \f${\lambda}_b\f$ for \f$1/\nu\f$ transport, 
  /// using Graz definition of the minor radius dr=ds/\<|grad(s)|\>, where s is the normalized toroidal flux.
  /// \sa lambda_b in V.Nemov,V.Kalyuzhnyj,S.Kasilov et al,PPCF,46,179(2004), http://dx.doi.org/10.1088/0741-3335/46/1/011
  /// @param s is the normalized toroidal flux
  /// \sa Note to FbsVmec()
  double GrazLambdab(double s);

  double FbsFtrap(double s); 

  double Fbsg2Vmec(double s);
  double Fbsg2Dkes(double s);
  double Fbsg2Graz(double s);


  /// The method returns flux surface average of g2, 
  /// where g2 is defined as B*grad(g2/B^2) = B x grad(s) * grad(1/B^2),
  /// s is the normalized toroidal flux.  
  /// @param s is the normalized toroidal flux
  /// \sa P.Helander et al,Phys.Plasmas 18,092505(2011), http://link.aip.org/link/doi/10.1063/1.3633940
  /// \note  At first call the method creates a look-up table using multithreading technique to spread 
  /// calculation among several CPU cores; this table is used for interpolating over the flux surface label.
  /// \n Use the following methods to set parameters:
  /// FbsSetXeffParam(), FbsSetTracingParam(), FbsSetIotaParam(), FbsSetMagnMomentParam()
  double Fbsg2(double s);
  /// The method returns flux surface average of u, 
  /// where u is defined as B*grad(u) = -B x grad(s) * grad(1/B^2),
  /// s is the normalized toroidal flux.  
  /// @param s is the normalized toroidal flux
  double Fbsu(double s); 
  /// The method returns flux surface average of u*B^2, 
  /// where u is defined as B*grad(u) = -B x grad(s) * grad(1/B^2),
  /// s is the normalized toroidal flux.  
  /// @param s is the normalized toroidal flux
  double FbsuB2(double s); 
  /// The method returns flux surface average of u^2*B^2, 
  /// where u is defined as B*grad(u) = -B x grad(s) * grad(1/B^2),
  /// s is the normalized toroidal flux.  
  /// @param s is the normalized toroidal flux
  double Fbsu2B2(double s); 
  /// The method returns flux surface average of B^2 
  /// @param s is the normalized toroidal flux
  double FbsB2(double s); 

  /// The method returns flux surface average of g4, 
  /// where g4 is defined as B*grad(g4/xi) = B x grad(s) * grad(1/xi) and
  /// xi = sqrt(1 - mu*B/Bmax),
  /// s is the normalized toroidal flux.  
  /// @param s is the normalized toroidal flux
  /// \sa P.Helander et al,Phys.Plasmas 18,092505(2011), http://link.aip.org/link/doi/10.1063/1.3633940
  /// \note  At first call the method creates a look-up table using multithreading technique to spread 
  /// calculation among several CPU cores; this table is used for interpolating over the flux surface label.
  /// \n Use the following methods to set parameters:
  /// FbsSetXeffParam(), FbsSetTracingParam(), FbsSetIotaParam(), FbsSetMagnMomentParam()
  double Fbsg4(double s, double mu);
  double Fbsg4Vmec(double s, double mu);

  /// The method returns bounce average of \f$\omega_d\f$ according formula
  /// \f$\left( \oint\frac{dl}{\sqrt{1-B/B_m}} \right)^{-1}\oint\frac{\omega_d dl}{\sqrt{1-B/B_m}}\f$ 
  /// for elecrons,
  /// where \f$\omega_d=\textbf{v}_d(\nabla\theta-\iota\nabla\varphi)\f$ , 
  ///  \f$\textbf{v}_d=\textbf{b}\times\left({\bf k}(1-B/B_m)+\nabla{B}/2B_m\right)\f$, \f$\textbf{b}={\bf B}/B\f$,
  ///  \f${\bf k}\f$ is the curvature vector
  /// @param mag0 is the magnetic coordinates of the point 
  /// from which integrating along field line is started and where \f$B_m=B(mag0)\f$
  /// @param dphi is the toroidal step in radians for fild line integration (0.4degree)
  /// @param iotaDenominator if not zero then this number
  ///   is the denominator in iota aproximation by rational number.
  /// @return Vector3d with bounce average of \f$\omega_d\f$ in the first element and 
  ///    bounce time \f$\oint\frac{dl}{\sqrt{1-B/B_m}}\f$ 
  ///    in the second element. 
  Vector3d omega_d_BounceAveraged(const Vector3d &mag0, double dphi=0.00698 , int iotaDenominator=0);
  
  /// The method returns the parallel adiabatic invariant with the integral taken between two turning points according formula
  /// \f$\int \sqrt{1-B/B_t} dl\f$ 
  /// @param mag is the magnetic coordinates of the point 
  /// from which integrating along field line is started and where \f$B_t=B(mag)\f$
  /// @param dphi is the toroidal step in radians for fild line integration (default value is 0.4 degree)
  /// @param iotaDenominator if not zero then this number
  ///   is the denominator in iota aproximation by rational number.
  /// @return std::vector with the parallel adiabatic invariant J ,
  ///   partial derivative dJ/ds ,  partial derivative dJ/dalpha (alpha=theta/iota - phi) 
  ///  with the bounce time \f$\oint\frac{dl}{\sqrt{1-B/B_t}}\f$ in the fourth element, 
  ///    \f$B_t=B(mag)\f$  in the fifth element 
  /// \sa  getMagCoordAtB()
  std::vector<double> Jinvariant(const Vector3d &mag, double dphi=0.00698 , int iotaDenominator=0);

  /// The method returns bounce average of \f$\omega_d\f$ ( see omega_d_BounceAveraged() )
  /// and the parallel adiabatic invariant with the integral taken between two turning points according formula
  /// \f$\int \sqrt{1-B/B_t} dl\f$ (see Jinvariant())
  /// @return std::vector<double> = {J, dJds, dJda, t_bounce, Bma, omegaDaverage} with the parallel adiabatic invariant J ,
  ///   partial derivative dJ/ds ,  partial derivative dJ/dalpha (alpha=theta/iota - phi) 
  ///  with the bounce time \f$\oint\frac{dl}{\sqrt{1-B/B_t}}\f$ in the fourth element, 
  ///    \f$B_t=B(mag)\f$  in the fifth element, and with 
  ///    the bounce average of \f$\omega_d\f$ in the sixth element 
  /// \sa  getMagCoordAtB(), omega_d_BounceAveraged(), Jinvariant()
  std::vector<double> Jinv_omega_d_BounceAveraged(const Vector3d &mag0, double dphi=0.00698 , int iotaDenominator=0);

  /// The method returns the magnetic coordinates of the point (that belongs to the given field line)
  /// at which the magnetic field is equal to the given value.
  /// @param mag defines the field line; 
  /// @param Bt is the magnetic field value, must be within the interval Bmin..Bmax for the considered field line;  
  /// @return Vector3d with the magnetic coordinates at which the magnetic field value is Bt or
  ///         Vector3d with (-1,Bmin,Bmax) if Bt is not within the interval Bmin..Bmax;                 
  /// \note the method seaches the point of interest starting from the point where B=Bmin, 
  ///        Bmin is seached within one machine period: 0<=phi<2*pi/Np
  /// \sa Jinvariant(), Bmin(double s), Bmax(double s)
  Vector3d getMagCoordAtB(const Vector3d & mag, double Bt);

  std::vector<Vector3d> getAllMagCoordAtB(const Vector3d & mag0, double Bt); 
  Vector3d getNextMagCoordAtB(const Vector3d & mag0, double Bt);

  /// The method returns  \<|grad(s)|\>  \f$\left\langle\left|\nabla s\right|\right\rangle\f$
  double GradsAvrg(double s); 
  /// The method returns  \<(grad(s))^2\> \f$\left\langle\left|\nabla s\right|^2\right\rangle\f$
  double Grads2Avrg(double s); 
  /// The method returns  \<(grad(s)/B)^2\> \f$\left\langle\left|\nabla s\right|^2/B^2\right\rangle\f$
  double Grads2overB2Avrg(double s); 
  /// The method returns the normalized poloidal flux, where s is the normalized toroidal flux.
  double sPol  (double s);
  /// The method returns the poloidal flux, where s is the normalized toroidal flux.
  double PoloidalFlux(double s);
  /// The method returns the normalized toroidal flux, where s is the normalized poloidal flux.
  double sToroidal(double sPoloidal);
  /// The method returns the iota, where s is the normalized toroidal flux.
  double iota  (double s) const;
  /// The method returns the pressure in Pa, where s is the normalized toroidal flux.
  double pressure  (double s) const;
  /// The method returns rational aproximation to iota: iota = numerator/denominator
  /// @param[in]  s is the normalized toroidal flux.
  /// @param[in]  denom_guess is the guess value for the denominator.
  /// @param[out] numerator is the numerator.
  /// @param[out] denominator is the denominator.
  /// @return iota = numerator/double(denominator).
  double iotaRational(double s,int denom_guess,int &numerator,int &denominator);
  /// The method returns the diota/ds, where s is the normalized toroidal flux.
  /// @param s is the normalized toroidal flux.
  double iotaPrime(double s) const;
  /// The method returns the poloidal current in [A].
  /// @param s is the normalized toroidal flux.
  double Ip    (double s) const;   // poloidal current [A]
  /// The method returns the toroidal current in [A].
  /// @param s is the normalized toroidal flux.
  double It    (double s) const;   // toroidal current [A]
  /// The method returns the toroidal flux in [Wb].
  /// @param s is the normalized toroidal flux.
  double Flux  (double s) const;   // toroidal flux [Tm^2]
  /// The method returns dp/ds [Pa], where p is the pressure in [Pa]
  /// @param s is the normalized toroidal flux.
  double pp    (double s);   // dp/ds,[Pa]  pressure gradient
  /// The method returns the area of the cross section of the surface \e s at angle \e fi
  /// @param s is the normalized toroidal flux
  /// @param ficyl is the cylindrical angle in radians
  /// @return area in [m]
  double crossSectionArea(double s, double ficyl);
  /// The method returns the average area of the cross section of the surface \e s
  /// @param s is the normalized toroidal flux.
  /// @return average area in [m]
  double averageCrossSectionArea(double s=1.);
  /// The method returns Volume in [m<sup>3</sup>] inside surface \b s.
  /// Volume is negative if (s,theta,phi)-coordinate system is left handed
  /// @param s is the normalized toroidal flux.
  /// @return Volume in [m<sup>3</sup>]
  double Volume(double s);   
  /// The method returns Vprime=dVolume/ds in [m<sup>3</sup>].
  /// Vprime is negative if (s,theta,phi)-coordinate system is left handed
  /// @param s is the normalized toroidal flux.
  double Vprime(double s) const; 
  /// The method calculates <b>dVolume/ds</b> in [m^3], integrating the Jacobian.
  /// Simpson method is used.
  /// @param s is the flux label - the normalized toroidal flux.
  /// @param nth specifies the number of points for poloidal integration.
  /// @param nph specifies the number of points for toroidal integration.
  /// @return dVolume/ds in [m<sup>3</sup>]
  /// \note this function is intended for testing of method Vprime(s);
  /// normally it's better to use  Vprime(s).
  /// dVolume/ds is negative if (s,theta,phi)-coordinate system is left handed
 double VprimeIntJ(double s,int nth=120,int nph=120);// dVolume/ds [m^3], as integral from Jacobian
  /// The method calculates the Susceptance matrix in (Flux,theta,phi)-coordinates.
  /// See definition of the Susceptance matrix in P.I.Strand,W.A.Houlberg, Physics of Plasmas, \b 8, 2782(2001)
  /// @param s is the normalized toroidal flux, s = Flux/Flux_max
  /// @param npol specifies the number of points for poloidal integration.
  /// @param ntor specifies the number of points for toroidal integration.
  /// @return Susceptance matrix S<sub>ij</sub>, see example:
  /// \code
  /// #include "CStconfig.h"
  /// int main() {
  ///   CStconfig mConf;
  ///   if(!mConf.load("w7x-sc1.bc")) exit(1);
  ///   Vector3d *S = mConf.SMatrix(0.5,200,200);
  ///   double S12 =  S[1][2];
  /// }
  /// \endcode
  /// \note The Susceptance matrix in (s,theta,phi)-coordinates can be expressed 
  /// through relation S_ij(Flux,theta,phi)=S_ij(s,theta,phi)* dFlux/ds
  /// \note The method is obsolete, use it for testing only. 
 ///  See S11(), S12(), S21(), S22() for faster and better alternative.
  Vector3d * SMatrix(double s, int npol=180, int ntor=200);
  /// Susceptance matrix: S11 in (Flux,theta,phi)-coordinates.
  /// The method calculates the Susceptance matrix in (Flux,theta,phi)-coordinates.
  /// See definition of the Susceptance matrix in P.I.Strand,W.A.Houlberg, Physics of Plasmas, \b 8, 2782(2001)
  /// @param s is the normalized toroidal flux, s = Flux/Flux_max
  /// @return S<sub>11</sub> element of susceptance matrix.
  double S11(double s); // S11(s)
  /// Susceptance matrix: S12 in (Flux,theta,phi)-coordinates.
  /// @param s is the normalized toroidal flux, s = Flux/Flux_max
  /// @return S<sub>12</sub> element of susceptance matrix.
  double S12(double s); // S12(s)
  /// Susceptance matrix: S21 in (Flux,theta,phi)-coordinates.
  double S21(double s); // S21(s)
  /// Susceptance matrix: S22 in (Flux,theta,phi)-coordinates.
  double S22(double s); // S22(s)
  /// The method calculates the flux average of normal and geodesic curvatures.
  /// @param s is the normalized toroidal flux.
  /// @param nth specifies the number of points for poloidal integration.
  /// @param nph specifies the number of points for toroidal integration.
  /// @return Vector3d k with k[0] = <|kN|> and k[1] = <|kG|>, see example:
  /// \code
  /// #include "CStconfig.h"
  /// int main() {
  ///   CStconfig mConf;
  ///   if(!mConf.load("w7x-sc1.bc")) exit(1);
  ///   Vector3d k = mConf.curvatureAvrg(0.5,200,200);
  ///   double kN =  k[0];
  ///   double kG =  k[1];
  /// }
  /// \endcode
  Vector3d curvatureAvrg(double s, int nth=120, int nph=120);
  //@}
  /// \name Surface quantity for estimation of trapped particle fraction
  //@{
  /// The method returns the position of Bmax on the surface \b s.
  void findBminBmaxPosition(double s, Vector3d &magcrdBmin, Vector3d &magcrdBmax, bool useInputAsGuess=false);
  /// The method returns the minimum value of magnetic field on the surface \b s.
  double Bmin  (double s);
  /// The method returns the maximum value of magnetic field on the surface \b s.
  double Bmax  (double s);
  /// The method returns the minimum value of magnetic field.
  double Bmin  ();
  /// The method returns the maximum value of magnetic field.
  double Bmax  ();
  /// The method returns the \<B\> -- flux average of magnetic field on the surface \b s.
  /// @param s is the normalized toroidal flux.
  double Bavrg (double s);
  /// The method returns the \<B<sup>2</sup>\> -- flux average of magnetic field on the surface \b s.
  /// @param s is the normalized toroidal flux.
  double B2avrg (double s);
  /// The method returns the volume average of B<sup>2</sup>.
  double B2VolAvrg();
  /// The method returns the \<B<sup>-2</sup>\> -- flux average of magnetic field on the surface \b s.
  /// @param s is the normalized toroidal flux.
  double B_2avrg(double s); 
  /// The method returns the \<grads^2*B^-2\> -- flux average on the surface \b s.
  /// @param s is the normalized toroidal flux.
  double grads2B_2avrg(double s); 
  /// The method returns the \<R<sup>2</sup>\> -- flux average of R^2 on the surface \b s.
  /// @param s is the normalized toroidal flux.
  double R2avrg(double s); 
  /// The method returns the \<R<sup>-2</sup>\> -- flux average of R^-2 on the surface \b s.
  /// @param s is the normalized toroidal flux.
  double R_2avrg(double s);
  /// The method returns the contribution of passed particles to the neoclassical polarization on the surface \b s.
  /// @param s is the normalized toroidal flux.
  double NeoPolarizationPass(double s);
  /// The method returns the contribution of traped particles to the neoclassical polarization on the surface \b s.
  /// @param s is the normalized toroidal flux.
  double NeoPolarizationTrap(double s);
  /// The method returns the neoclassical polarization on the surface \b s.
  /// @param s is the normalized toroidal flux.
  double NeoPolarization(double s);
  /// The method returns the value of \<(1-b)<sup>0.5</sup>\>, where  b=B(s,th,ph)/Bmax(s).
  /// @param s is the normalized toroidal flux.
  double db12Avrg(double s);  // <(1-b)^0.5>; b = B(s,th,phi)/Bmax(s)
  /// The method returns the value of \<(1-b)<sup>1.5</sup>\>, where  b=B(s,th,ph)/Bmax(s).
  /// @param s is the normalized toroidal flux.
  double db32Avrg(double s);  // <(1-b)^1.5>
  /// The method returns the fraction of trapped particles on the surface \b s.
  /// @param s is the normalized toroidal flux.
  double ftrapped(double s);
// <ro_paralell^-1>
// f(s,x) = <b/sqrt(1-x*b)>, where 0<=x<=1, 0<=s<=1, and b=B(s,th,ph)/Bmax(s)
  double roParInvAvrg(double s, double x) const;
  /// The method returns the \<p\> -- flux average of pitch angle on the surface \b s.
  /// The method returns the value of function P(s,x) = \<sqrt(1-x*B(s,th,ph)/Bmax(s))\>, where 0<=x<=1, 0<=s<=1,
  /// and definition of x is the following x=(1-p^2)/b, p is the pitch angle, b=B(s,th,ph)/Bmax(s).
  double pavrg(double s, double x) const;
  /// The method returns the Integral(from x to 1) dx'/P(s,x') on the surface \b s.
  /// The method returns the Integral(from x to 1) dx'/P(s,x'), where P(s,x) = \<sqrt(1-x*B(s,th,ph)/Bmax(s))\>,
  //  0<=x<=1, 0<=s<=1  and definition of x is the following x=(1-p^2)/b, p is the pitch angle, b=B(s,th,ph)/Bmax(s).
  double integralInvPavrg(double s, double x) const;
  /// The method returns Bmin and Bmax -- the minimum and maximum values of magnetic field on the surface \b s.
  /// @param s is the flux surface label, 0<=s<=1.
  /// @param Bmin is the reference for returned value of Bmin.
  /// @param Bmax is the reference for returned value of Bmax.
  /// @param npol specifies the number of points in poloidal direction,
  /// which are used for finding of minimum and maximum.
  /// @param ntor specifies the number of points in toroidal direction,
  /// which are used for finding of minimum and maximum.
  /// @return Bmin,Bmax -- the minimum and maximum values of magnetic field on the surface \b s.
  /// \note this function is intended for testing of methods Bmin() and Bmax();
  /// normally it's better to use  Bmin() or Bmax().
  void getBminmax(double s, double &Bmin,double &Bmax,int npol=120,int ntor=120);
  /// The method returns the fraction of trapped particles on the surface \b s.
  /// @param s is the flux surface label, 0<=s<=1.
  /// @param Ftrapped is the reference for returned value of the trapped particles fraction.
  /// @param Bmin is the reference for returned value of Bmin.
  /// @param Bmax is the reference for returned value of Bmax.
  /// @param Bavr is the flux average B on the surface \b s.
  /// @param B2avr is the flux average B^2 on the surface \b s.
  /// @param nth specifies the number of points in poloidal direction,
  /// which are used for integration (surface averaging).
  /// @param nph specifies the number of points in toroidal direction,
  /// which are used for integration.
  /// @return
  ///       \li  Ftrapped -- the fraction of trapped particles on the surface \b s.
  ///       \li  Bmin,Bmax -- the minimum and maximum values of magnetic field on the surface \b s.
  ///       \li  Bavr,B2avr -- are the flux average B and B^2 on the surface \b s.
  /// \note this function is intended for internal use and testing of method ftrapped();
  /// normally it's better to use  ftrapped(),Bmin(),Bmax(),Bavrg(),B2avrg().
  void getFtrappedBminmax(double s,double &Ftrapped, double &Bmin,double &Bmax,
                          double &Bavr,double &B2avr,int nth=120,int nph=120);
  //@}
  /// \name Splined Fourier coefficients
  //@{
  // Splined Fourier coefficients
  /// The method returns Fourier coefficient of |B|-expansion.
  /// Cubic spline interpolation is used.
  /// @param s is the flux surface label, 0<=s<=1.
  /// @param m is the poloidal mode number, 0<=m<= Mpol().
  /// @param n is the toroidal mode number. |n|<=Ntor().
  /// @return Fourier coefficient.
  /// @note extrapolation to s>1 is possible, but is not reliable for big \b s.
  double  Bmn_sp  (double s, int m, int n) const;
  /// The method returns Fourier coefficient of R-expansion, for parameters see Bmn_sp()
  double  Rmn_sp  (double s, int m, int n) const;
  /// The method returns Fourier coefficient of Z-expansion, for parameters see Bmn_sp()
  double  Zmn_sp  (double s, int m, int n) const;
  /// The method returns Fourier coefficient of Phi-expansion, for parameters see Bmn_sp()
  double  Phimn_sp(double s, int m, int n) const;
  /// The method returns Fourier coefficient of VMEC lambda expansion, for parameters see Bmn_sp()
  double  lmn_sp  (double s, int m, int n) const;
  /// The method returns Fourier coefficient of VMEC Jacobian expansion, for parameters see Bmn_sp()
  double  gmn_sp  (double s, int m, int n) const;
  /// The method returns Fourier coefficient of expansion of poloidal contravariant B, for parameters see Bmn_sp()
  double  BcontrPolmn_sp  (double s, int m, int n) const;
  /// The method returns Fourier coefficient of expansion of toroidal contravariant B, for parameters see Bmn_sp()
  double  BcontrTormn_sp  (double s, int m, int n) const;

  double Bmn_sp_prime(double s, int m, int n) const;

  double Boozer2Bmn(double s,int m,int n);
  double BoozerBmn(double s,int m,int n);

  //@}
  /// \name s-grid info
  //@{
  /// The method returns flux surface label \b s(norm.tor.flux)
  /// @param is is the flux surface number, 0<=is< nSurfaces().
  double get_s(int is) const { return ms[is];}
  /// The method returns the \b s value of the first surface stored in a file.
  double get1st_s() const { return s1st;}   // get the first s, which is stored in a bc-file.
  /// The method returns the \b s value of the last surface stored in a file.
  double getLast_s() const { return sLast;}  // get the last  s, which is stored in a bc-file
  /// The method returns the index of the first flux surface in array CStconfig::ms.
  /// The first surface means here the first surface stored in a file.
  int getIdx1st_s() const { return m1stSindex;}  // ms[]-index of the first surface from bc-file
  /// The method returns the index of the last flux surface in array CStconfig::ms.
  /// The last surface means here the last surface stored in a file.
  int getIdxLst_s() const { return mLCMSindex;}  // ms[]-index of the LCMS from bc-file
  /// The method returns the index in array CStconfig::ms, which corresponds to \b s.
  /// The left index of segment ms[i] - ms[i+1] in which \b s lies is returned.
  int SearchSindx(double s, bool halfMesh=true) const;
  //@}
  /// \name |B| and Jacobian
  //@{
  /// The function returns  B=sum_mn B_mn(s)*cos(m*theta-n*Np*phi).
  /// @param magCoord is the vector with magnetic coordinates (s,theta,phi).
  double B(const Vector3d &magCoord) const; // |B|
  /// The method returns |B| on the field line that passes through
  ///   the point given by (s0,theta0,phi0). 
  /// @param[in] magCoord0 is the point of B-line in magnetic coordinates
  /// @param[in] phi is the toroidal angle
  /// @return  B(s0,theta(phi),phi)
  /// \sa magCoordFieldLine()
  double BonFieldLine(const Vector3d &magCoord0, double phi); 
  /// The method returns magnetic coordinates of the point of a field line, 
  /// that passes through the given point (s0,theta0,phi0). 
  /// The method  calculates function
  ///  \f[\theta(\varphi)=\theta_{0}+\iota(s_{0})(\varphi-\varphi_{0}) - [\lambda(s_{0},\theta,\varphi)-\lambda(s_{0},\theta_{0},\varphi_{0})]\f]
  /// where \f$\lambda\f$ is non zero for VMEC coordinates.
  /// @param magCoord0 is the point of B-line in magnetic coordinates
  /// @param phi is the toroidal angle
  /// @param iota is the rational approximation to the iota, the genuine iota is used if this parameter is zero. 
  /// @return  Vector3d(s0,theta(phi),phi) 
  Vector3d magCoordFieldLine(const Vector3d &magCoord0, double phi, double iota=0);
private:
  /// TODO: does not work properly
  Vector3d magCoordFieldLine2(const Vector3d &magCoord0, double theta, double iota=0);
public:
  /// The method returns magnetic coordinates of the point of a field line, 
  /// that passes through the given point (s0,theta0,ficyl0). 
  /// The method  calculates function
  ///  \f[\theta(\varphi)=\theta_{0}+\iota(s_{0})(\varphi-\varphi_{0}) - [\lambda(s_{0},\theta,\varphi)-\lambda(s_{0},\theta_{0},\varphi_{0})]\f]
  /// where \f$\lambda\f$ is non zero for VMEC coordinates.
  /// @param mixCoord0 is the point of B-line in mixed coordinates,
  ///     where s, theta coincide with correspondent magnetic coordinates, 
  ///     but toroidal angle is cylindricyl angle
  /// @param fcyl is the cylindricyl angle
  /// @return  Vector3d(s0,theta(ficyl),phi(ficyl)) 
  /// \note for VMEC coordinats this method is the same as magCoordFieldLine()
  Vector3d mixCoordFieldLine(const Vector3d &mixCoord0, double fcyl);
  /// The method returns grad(|B|).
  /// @param magCoord is the magnetic coordinates
  /// @return grad(|B|) in cylindrical coordinates
  Vector3d  gradB(const Vector3d &magCoord);
  /// The method returns covariant components of  grad(|B|).
  /// @param magCoord is the magnetic coordinates
  /// @return Vector3d with covariant components of grad(|B|),
  ///       Vector3d(dB/ds,dB/dtheta,dB/dphi)
  Vector3d  getCovaGradB(const Vector3d &magCoord) const;
  /// The method returns B and gradients in cylindrical coordinates
  ///  @param[in]  u is the point (s,theta,phi) in magnetic coordinates,
  ///  @param[out]  B  is the magnetic field  in cylindrical coordinates,
  ///  @param[out]  gradB  is the grad(|B|), 
  ///  @param[out]  grads  is the grad(s),
  ///  @param[out]  gradTh is the grad(theta),
  ///  @param[out]  gradPh is the grad(phi),
  void getBandGradients(const Vector3d &u, Vector3d &B, Vector3d &gradB, 
                                  Vector3d &grads, Vector3d &gradTh, Vector3d &gradPh);

  /// Get B and gradients in cartesian coordinates
  ///  @param[in] Vector3d u(s,theta,phi) is a point in magnetic coordinates,
  ///  @param[out]  B  is the magnetic field  in cartesian coordinates
  ///  @param[out]  gradB  is the grad(|B|) in cartesian coordinates
  ///  @param[out]  grads  is the grad(s)
  ///  @param[out]  gradTh is the grad(theta)
  ///  @param[out]  gradPh is the grad(phi) 
  void getBandGradientsxyz(const Vector3d &u, Vector3d &B, Vector3d &gradB, 
                                  Vector3d &grads, Vector3d &gradTh, Vector3d &gradPh);

  /// Get B and basis vectors in cartesian coordinates
  ///  @param[in]   u = (s,theta,phi) is the point magnetic-coordinates,
  ///  @param[out]  xyz  is the point cartesian-coordinates,
  ///  @param[out]  Bxyz  is the magnetic field  in cartesian coordinates
  ///  @param[out]  gradBxyz  is the grad(|B|) in cartesian coordinates
  ///  @param[out]  e1con, e2con, e3con are the contravariant-basis vectors e^i = grad(u^i) in cartesian coordinates
  ///  @param[out]  e1cov, e2cov, e3cov are the covariant-basis vectors e_i = dX/du^i in cartesian coordinates
  void  getBandBasisVectorsxyz(const Vector3d &u, Vector3d &xyz, Vector3d &Bxyz, Vector3d &gradBxyz, 
                            Vector3d &e1con, Vector3d &e2con, Vector3d &e3con,
                            Vector3d &e1cov, Vector3d &e2cov, Vector3d &e3cov);

  /// @return s
  double getdBGradscyl(const Vector3d &cyl,Vector3d &B,
                       Vector3d &dBdr,Vector3d &dBdfir,Vector3d &dBdz,Vector3d &grads);
  /// @return s
  double getdBGradsxyz(const Vector3d &xyz,Vector3d &B,
                       Vector3d &dBdx,Vector3d &dBdy,Vector3d &dBdz,Vector3d &grads );

  /// The method returns B and curl(B) in cylindrical coordinates
  ///  @param[in]  u is the point (s,theta,phi) in magnetic coordinates,
  ///  @param[out]  B  is the magnetic field  in cylindrical coordinates,
  ///  @param[out]  curlB  is the curl(B), 
  ///  @param[out]  grads  is the grad(s)
  ///  @return cylindrical coordinates of the point (s,theta,phi)
  Vector3d getBcurlBGradscyl(const Vector3d &u,Vector3d &B,Vector3d &curlB,Vector3d &grads);

  /// The method returns the local shear.
  /// @param[in] Vector3d magCoord(s,theta,phi) is a point in magnetic coordinates,
  /// @return  local_shear = -H*curl(H), where H = grad(s)/abs(grad(s) ^ B/abs(B), ^ is the cross product
  /// @return cyl, B, grads in cylindrical coordinates
  double getLocalShear(const Vector3d &magCrd,Vector3d &B,Vector3d &cyl,Vector3d &grads);
  double getLocalShear(const Vector3d &magCrd);

  /// If sw==true then the method forces to use Jacobian=dX/dFlux*(dX/dtheta^dX/dphi) for calculating B field.
  void useMixedProductForJacobian(bool sw);

  /// Jacobian: dV =  Jacobian(s,theta,phi)*ds*dtheta*dphi.
  /// @param magCoord is the vector with magnetic coordinates (s,theta,phi).
  /// @return the value of Jacobian at the point (s,theta,phi)
  /// @param[out] *Bmod is the B(s,theta,phi) if Bmod is not NULL
  ///  \note Jacobian in coordinates (s,theta,phi) \n
  ///   Boozer Jacobian = FluxTor'*mu0*(Ip*Np+iota*It)/B^2/twopi^2) \n
  ///   Tokamak Jacobian = FluxTor'*R*R/(mu0*Ipol) \n
  ///   Vmec Jacobian = Sum gmn*cos(....) \n
  ///   Jacobian has sign
  double Jacobian(const Vector3d &magCoord, double *Bmod=0) const; 
  //@}
  /// \name Coordinate transformation from magnetic coordinates to real space
  //@{
  /// The method returns the cartesian coordinates of the point.
  /// @param magCoord is the Vector3d(s,theta,phi) in magnetic coordinates
  /// @return cartesian coordinates.
  Vector3d mag2xyz(const Vector3d &magCoord) const; // return vector (x,y,z) for point (s,theta,phi)
  /// The method makes transformation from magnetic coordinates to cartesian.
  /// \overload
  /// @param s,theta,phi are the magnetic coordinates.
  /// @return cartesian coordinates.
  Vector3d mag2xyz(double s,double theta,double phi) const;
  /// The method makes transformation from magnetic coordinates to cylindrical.
  /// @param magCoord is the Vector3d(s,theta,phi) in magnetic coordinates
  /// @return cylindrical coordinates.
  Vector3d mag2cyl(const Vector3d &magCoord) const; // return vector (R,fi,Z) for point (s,theta,phi)
  /// The method makes transformation from magnetic coordinates to cylindrical.
  /// \overload
  /// @param s,theta,phi are the magnetic coordinates.
  /// @return cylindrical coordinates.
  Vector3d mag2cyl(double s,double theta,double phi) const;
  /// The method returns cylindrical coordinates which correspond to the mixed coordinates.
  /// @param s is the flux surface label
  /// @param theta is the magnetic poloidal angle.
  /// @param fi is the cylindrical angle.
  /// @return cylindrical coordinates.
  /// @param[out] phi if not NULL then this address at return will contain the magnetic toroidal angle.
  /// @param[in] sdat for internal use only.
  Vector3d mixcoord2cyl(double s,double theta,double fi,double* phi=0,
                                 sumDat *sdat=0) const;
  /// The method returns cylindrical coordinates which correspond to the mixed coordinates.
  /// \overload
  /// @param mixcoord are the mixed coordinates.
  /// @param phi if not NULL then this address at return will contain the magnetic toroidal angle.
  Vector3d mixcoord2cyl(const Vector3d &mixcoord, double* phi=NULL) const;
  /// The method returns cartesian coordinates which correspond to the mixed coordinates.
  /// @param mixcoord are the mixed coordinates.
  /// @return cartesian coordinates.
  /// @param phi if not NULL then this address at return will contain the magnetic toroidal angle.
  Vector3d mixcoord2xyz(const Vector3d &mixcoord, double* phi=NULL) const;
  /// The method returns cylindrical coordinates of the magnetic axis.
  /// @param ficyl is the value of cylindrical angle in radians.
  /// @return cylindrical coordinates of the magnetic axis.
  Vector3d Axiscyl(double ficyl) const { return mixcoord2cyl(0,0,ficyl);}
  //@}
  /// \name Coordinate transformation from real to magnetic coordinates
  //@{
  /// The method forces to use the stream function lambda for B-field calculatioins,
  /// the method works only for VMEC representation, otherwise it does nothing.
  void useVmecLambdaOn() {
    useVmecLambda=true;
  }
  /// The method forces to use the VMEC contravariant B in B-field calculatioins,
  /// the method works only for VMEC representation, otherwise it does nothing.
  void useVmecLambdaOff() {
    if(vmecCoordinates) useVmecLambda=false;
  }
//Newton's method
  /// The method calculates Jacobian matrix for Newton's iteration.
  /// @param magCoord is the vector with magnetic coordinates (s,theta,phi).
  /// @return pointer to the array J of 4 vectors calculated in the point \b magCoord
  /// \code
  ///   (R )(Rs )(Rt )(Rp )           (R )       (Rs )
  ///J= (fi)(fis)(fit)(fip), i.e J[0]=(fi), J[1]=(fis), ...
  ///   (Z )(Zs )(Zt )(Zp )           (Z )       (Zs )
  ///  where (R,fi,Z) is the corresponding point
  ///  in cylindrical coordinates;
  ///  Rs  = dR/ds,  Rt  = dR/dtheta,  Rp  = dR/dphi;
  ///  fis = dfi/ds, fit = dfi/dtheta, fip = dfi/dphi;
  /// \endcode
  /// \note Along with the Jacobian matrix calculation the method calculates other quantities,
  /// which are available through methods: getmodB(), getJac(), getJac2(), getJ(),
  /// getJ2(), getBxyz(), getBcyl(), getGradFlux(), getGrads(), getGradreff(),getMetricTensor().
  Vector3d *NJacobian (const Vector3d &magCoord);
  /// The method calculates Jacobian matrix for Newton's iteration using linear interpolation on \em s.
  /// see CStconfig::NJacobian
  /// \note Don't use this method, it is needed only for testing.
  Vector3d *NJacobianL(const Vector3d &magCoord);
  /// The method makes transformation from cylindrical to magnetic coordinates.
  /// The method uses Newton's iteration to find magnetic coordinates, the guess value is
  /// estimated by CStconfig::cyl2guess or is taken from previous call.
  /// @param cyl is the cylindrical coordinates
  /// @return magnetic coordinates.
  /// \note After calling this function the other quantities are available through methods:
  /// getmodB(), getJac(), getJac2(), getJ(), getJ2(), getBxyz(),
  /// getBcyl(), getGradPsi(), getGrads(), getGradreff().
  Vector3d  cyl2mag(const Vector3d &cyl); // return vector (s,theta,phi) for point (R,fi,Z)
  /// The method makes transformation from cylindrical to magnetic coordinates using guess value.
  /// @param cyl is the cylindrical coordinates
  /// @param magCoordGuess is the guess value
  /// @return magnetic coordinates.
  /// \note After calling this function the other quantities are available through methods:
  /// getmodB(), getJac(), getJac2(), getJ(), getJ2(), getBxyz(),
  /// getBcyl(), getGradPsi(), getGrads(), getGradreff().
  /// \note This method is obsolete, use cyl2mag(const Vector3d &cyl)
  Vector3d  cyl2mag(const Vector3d &cyl,const Vector3d &magCoordGuess);
private:
  Vector3d  cyl2mag(const Vector3d &cyl, const Vector3d &magCoordGuess, Vector3d *Jguess);
public:
  /// The method makes transformation from cartesian to magnetic coordinates.
  /// The method uses Newton's iteration to find magnetic coordinates; the guess value is
  /// estimated by CStconfig::cyl2guess or is taken from previous call.
  /// @param xyz is the cartesian coordinates
  /// @return magnetic coordinates.
  /// \note After calling this function the other quantities are available through methods:
  /// getmodB(), getJac(), getJac2(), getJ(), getJ2(), getBxyz(),
  /// getBcyl(), getGradPsi(), getGrads(), getGradreff().
  Vector3d  xyz2mag(const Vector3d &xyz); // return vector (s,theta,phi) for point (x,y,z)
  /// The method finds the flux label \b s which corresponds to the point in space.
  /// @param cyl is the input parameter with the cylindrical coordinates of the point of interest.
  /// @return flux label \b s.
  double   cyl2s(const Vector3d &cyl) {return cyl2mag(cyl)[0];}
  /// The method finds the flux label \b s which corresponds to the point in space.
  /// @param xyz is the input parameter with the cartesian coordinates of the point of interest.
  /// @return flux label \b s.
  double   xyz2s(const Vector3d &xyz) {return xyz2mag(xyz)[0];}
  /// The method makes transformation from cylindrical to magnetic coordinates using
  /// 2d-iteration on R-Z plane.
  /// The algorithm is robust and rather fast, use it for getting 
  /// the guess value for the Newton method implemented in cyl2mag().
  /// @param[in]  cyl is the cylindrical coordinates
  /// @param[in]  eps is the relative accuracy of the transformation, 
  ///             use eps ~ 0.06 (for W7-X plasma, 0.06 is 3cm/50cm ) 
  ///             if you need only guess for Newton method,
  ///             use eps ~ 0.001  if you need exact coordinate transformation.     
  /// @param[in]  smax is the value of the 'extra' surface beyond LCMS, use smax <= 1.2 
  /// @param[out] inside is the address of the flag that is true if cyl is inside smax.
  /// @return magnetic coordinates.
  /// \note the quantities
  /// getmodB(), getJac(), getJac2(), getJ(), getJ2(), getBxyz(),
  /// getBcyl(), getGradPsi(), getGrads(), getGradreff()  <em><b>are not available</b></em>
  /// after calling this function. Call NJacobian() in oder to calculate these quantities.
  Vector3d  cyl2guess (const Vector3d &cyl,double eps, double smax=1, bool *inside=0); 
//The following function must be called only after calling of NJacobian, or xyz2mag, or cyl2mag
  /// The method returns |B|.
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  double     getmodB() const;   // |B|
  /// The method returns |B|.
  double     getmodB(const Vector3d &cyl); 
  /// The method returns Jacobian=mu0*(Ipol+iota*Itor)/(4pi<sup>2</sup>B<sup>2</sup>).
  /// Jacobian for coordinates (Flux,theta,phi)
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  double     getJac() const;  // Jacobian in coordinates (torFlux,theta,phi)
  /// The method returns Jacobian=mu0*(Ipol+iota*Itor)/(4pi<sup>2</sup>B<sup>2</sup>).
  double     getJac(const Vector3d &cyl);  // Jacobian in coordinates (torFlux,theta,phi)
  /// The method returns Jacobian=Flux(1)*mu0*(Ipol+iota*Itor)/(4pi<sup>2</sup>B<sup>2</sup>).
  /// Jacobian for coordinates (s,theta,phi)
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  double     getJ() const;      // Jacobian in coordinates (s,theta,phi)
  /// Jacobian for coordinates (s,theta,phi)
  double     getJ(const Vector3d &cyl);      // Jacobian in coordinates (s,theta,phi)
  /// The method returns Jacobian=dX/dFlux*(dX/dtheta^dX/dphi).
  /// Jacobian for coordinates (Flux,theta,phi)
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  double     getJac2() const;
  /// The method returns Jacobian=dX/dFlux*(dX/dtheta^dX/dphi).
  /// Jacobian for coordinates (Flux,theta,phi)
  double     getJac2(const Vector3d &cyl);
  /// The method returns Jacobian=dX/ds*(dX/dtheta^dX/dphi).
  /// Jacobian for coordinates (s,theta,phi)
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  // Jacobian(s,theta,phi) = R*(R,s*R,t^R,p)) -- jacobian, where ^ is a cross product
  double     getJ2() const;
  /// The method returns Jacobian=dX/ds*(dX/dtheta^dX/dphi).
  /// Jacobian for coordinates (s,theta,phi)
  double     getJ2(const Vector3d &cyl);
  /// The method returns cartesian B-vector.
  /// \note Use this function after calling NJacobian(),cyl2mag(),or xyz2mag().
  ///  \sa NJacobian().
  const Vector3d & getBxyz();   // B-vector in cartesian coord.
  /// The method returns cartesian B-vector at cartesian point xyz.
  const Vector3d & getBxyz(const Vector3d &xyz);
  /// The method returns cylindrical B-vector.
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  const Vector3d & getBcyl() const;   // B-vector in cylindrical coord.
  /// The method returns cylindrical B-vector.
  const Vector3d & getBcyl(const Vector3d &cyl); 
  /// The method returns cylindrical vector grad(theta).
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  const Vector3d & getGradtheta() const;
  /// The method returns cylindrical vector grad(theta).
  const Vector3d & getGradtheta(const Vector3d &cyl); 
  /// The method returns cylindrical vector grad(lambda) for VMEC coordinates,
  /// where lambda is the stream function.
  Vector3d getGradlambda(const Vector3d &cyl);
  /// The method returns cylindrical vector grad(lambda) for VMEC coordinates,
  /// where lambda is the stream function.
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  Vector3d getGradlambda() const;

  /// The method returns cylindrical vector grad(phi).
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  const Vector3d & getGradphi() const;
  /// The method returns cylindrical vector grad(phi).
  const Vector3d & getGradphi(const Vector3d &cyl);
  /// The method returns cylindrical vector grad(Flux).
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  const Vector3d & getGradFlux() const;
  /// The method returns cylindrical vector grad(Flux).
  const Vector3d & getGradFlux(const Vector3d &cyl);
  /// The method returns cylindrical vector grad(s).
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  const Vector3d & getGrads() const; 
  /// The method returns cylindrical vector grad(s).
  const Vector3d & getGrads(const Vector3d &cyl);  // grad(s)
  /// The method returns cartesian vector grad(s).
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  Vector3d  getGradsxyz() const;
  /// The method returns cartesian vector grad(s).
  Vector3d  getGradsxyz(const Vector3d &xyz);
  /// The method returns cylindrical vector grad(r<sub>eff</sub>).
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  const Vector3d & getGradreff() const;// grad(r_eff)
  /// The method returns cylindrical vector grad(r<sub>eff</sub>).
  const Vector3d & getGradreff(const Vector3d &cyl); 
  /// The method retuns pointer to the Jacobian matrix for Newton's iteration 
  /// calculated by CStconfig::NJacobian().
  /// This matrix is calculated during coordinate transformation inside following functions:
  /// cyl2mag(),xyz2mag(),cyl2s(),xyz2s(). \sa CStconfig::NJacobian()
  const Vector3d * getNJmatrix();
  /// The method retuns pointer to the Jacobian matrix for Newton's iteration 
  /// calculated by CStconfig::NJacobian() for point \c cyl.
  const Vector3d * getNJmatrix(const Vector3d &cyl); 
  /// The method retuns pointer to the MetricTensor calculated by CStconfig::NJacobian.
  /// This tensor is for coordinates(s,theta,phi) and calculated also during coordinate
  /// transformation inside following functions:
  /// cyl2mag(),xyz2mag(),cyl2s(),xyz2s(). \sa CStconfig::NJacobian
  /// @return metric tensor g<sub>ij</sub>, where i=(0,1,2)=(s,theta,phi)
  ///see example:
  /// \code
  /// #include "CStconfig.h"
  /// int main() {
  ///   CStconfig mConf;
  ///   if(!mConf.load("w7x-sc1.bc")) exit(1);
  ///   Vector3d magCoord(0.5,1.,2.); // magnetic coordinates
  ///   mConf.NJacobian(magCoord);
  ///   Vector3d *g = mConf.getMetricTensor();
  ///   double g12 =  g[1][2];  // g_theta,phi
  /// }
  /// \endcode
  const Vector3d * getMetricTensor();
  /// The method returns pointer to the metricTensor calculated in the point \c cyl.
  /// @return metric tensor g<sub>ij</sub>, where i=(0,1,2)=(s,theta,phi)
  ///see example:
  /// \code
  /// #include "CStconfig.h"
  /// int main() {
  ///   CStconfig mConf;
  ///   if(!mConf.load("w7x-sc1.bc")) exit(1);
  ///   Vector3d cyl(5.5,2.,0.5); // cylindrical coordinates
  ///   Vector3d *g = mConf.getMetricTensor(cyl);
  ///   double g12 =  g[1][2];    // g_theta,phi
  /// }
  /// \endcode
  const Vector3d * getMetricTensor(const Vector3d &cyl);
  /// The method returns  the cylindrical coordinates of the last point used 
  /// in coordinate transformation.
  /// \sa cyl2mag() or xyz2mag()
  Vector3d  getRcyl() const;
  /// The method returns  the cartesian coordinates of the last point used 
  /// in coordinate transformation.
  /// \sa cyl2mag() or xyz2mag()
  Vector3d  getRxyz() const;
  /// The method returns covariant components of B = (B_s,B_theta,B_phi) in (s,theta,phi) coordinates.
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  const Vector3d & getBcova() const;
  /// The method returns contravariant components of B = (0,B^theta,B^phi).
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  const Vector3d & getBcontra() const;
  /// The method returns covariant-basis vectors e_i = dX/du^i in cylindrical coordinates.
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  void getCovaBasis  (Vector3d &e1, Vector3d &e2, Vector3d &e3) const;
  /// The method returns contravariant-basis vectors e^i = grad(u^i) in cylindrical coordinates.
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  void getContraBasis(Vector3d &e1, Vector3d &e2, Vector3d &e3) const;
  /// The method returns grad(|B|) in cylindrical coordinates.
  /// @return grad(|B|) in cylindrical coordinates.
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  Vector3d  getGradBcyl();
  /// The method returns grad(|B|) in cylindrical coordinates.
  /// @param cyl is the point in cylindrical coordinates
  /// @return grad(|B|) in cylindrical coordinates.
  Vector3d  getGradBcyl(const Vector3d &cyl);
  /// The method returns grad(|B|) in cartesian coordinates.
  /// @return grad(|B|) in cartesian coordinates.
  /// \note Use this function after calling NJacobian(),cyl2mag(),xyz2mag(),cyl2s(),xyz2s().
  /// \sa NJacobian().
  Vector3d  getGradBxyz();
  /// The method returns grad(|B|) in cartesian coordinates.
  /// @param xyz is the point in cartesian coordinates.
  /// @return grad(|B|) in cartesian coordinates.
  Vector3d  getGradBxyz(const Vector3d &xyz);

  /// VMEC coordinates vmecCoord=(s,theta,phi)
  /// @return PEST coordinates
  /** \brief The method transforms the VMEC coordinates to PEST.
  * Transformation formulas
  * \f{eqnarray*}
  * \theta_{pest} & = & \theta_{V}+\lambda(s,\theta_{V},\varphi_{V})\\
  * \varphi_{pest} & = & \varphi_{V}
  * \f}
  * @param vmec is the VMEC coordinates (s,theta_v,phi)
  * @return PEST coordinates (s,theta_pest,phi)
  * \note the input file must be the VMEC file
  * \sa \ref PESTCoordinates
  */
  Vector3d Vmec2Pest(const Vector3d &vmec);
  Vector3d Pest2Vmec(const Vector3d &pest);
  Vector3d SFLcoordOnBlineAtTheta(const Vector3d &pest0, double theta, double iotA=0);
  void SFLcontraBasis(const Vector3d &pest, Vector3d &e1, Vector3d &e2, Vector3d &e3, Vector3d &vmec);
  void BoozerContraBasis(const Vector3d &booz, Vector3d &e1, Vector3d &e2, Vector3d &e3, Vector3d &vmec);
  void Vmec2PestContraBasis(const Vector3d &vmec, Vector3d &e1, Vector3d &e2, Vector3d &e3, Vector3d &pest);
  //@}
  /// \name Transformation between VMEC and Boozer Coordinates
  //@{
  /** \brief The method transforms the VMEC coordinates to Boozer.
  * Transformation formulas
  * \f{eqnarray*}
  * \theta_{B} & = & \theta_{V}+\lambda(s,\theta_{V},\varphi_{V})+
  * \iota h(s,\theta_{V},\varphi_{V})\\
  * \varphi_{B} & = & \varphi_{V}+h(s,\theta_{V},\varphi_{V})
  * \f}
  * @param vmecCoord is the VMEC coordinates (s,theta_v,phi_v)
  * @return Boozer coordinates (s,theta_B,phi_B)
  * \note the input file must be the VMEC file
  * \sa \ref VMEC2Boozer
  */
  Vector3d Vmec2Boozer(const Vector3d &vmecCoord);
  /** \brief The method transforms the Boozer coordinates to VMEC.
  * Transformation formulas
  * \f{eqnarray*}
  * \theta_{B} & = & \theta_{V}+\lambda(s,\theta_{V},\varphi_{V})+
  * \iota h(s,\theta_{V},\varphi_{V})\\
  * \varphi_{B} & = & \varphi_{V}+h(s,\theta_{V},\varphi_{V})
  * \f}
  * @param boozer is the Boozer coordinates (s,theta_B,phi_B)
  * @return VMEC coordinates (s,theta_v,phi_v)
  * \note the input file must be the VMEC file
  * \sa \ref VMEC2Boozer
  */
  Vector3d Boozer2Vmec(const Vector3d &boozer);
  /// The method transforms the VMEC wout-file to Boozer-coordinate data file.
  /// @param fname is the filename under which the configuration is saved.
  /// @param method if 0 then Boozer spectrum are obtained by integrating over Boozer angles,
  ///  otherwise the integrations over VMEC angles are performed;
  ///  the method 0 is faster, because summation (numerical integration)
  ///  is performed using FFT;
  /// @param M is the number of poloidal mods, 
  ///  if M==0 then the number of Fourier mods is three times of that of VMEC expansion;    
  /// @param N is the number of toroidal mods, 
  ///  if M==0 then the number of Fourier mods is three times of that of VMEC expansion:   
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// \code
  /// #include "CStconfig.h"
  /// int main() {
  ///   MConf::CStconfig mConf;
  ///   if(!mConf.load("wout.w7x-sc1.txt")) exit(1);
  ///   mConf.writeVmec2Boozer("w7x-sc1.bc");
  /// }
  /// \endcode
  /// \sa \ref VMEC2Boozer
  bool writeVmec2Boozer(const char * fname, int method=0, int M=0, int N=0);
  /// The method peforms test, that 
  /// d(theta_booz,phi_booz)/d(theta_vmec,phi_vmec) == (Jacobian_vmec/Jacobian_booz).
  /// @param vmecCoord is the VMEC coordinates (s,theta_v,phi_v)
  /// @return ratio { d(theta_booz,phi_booz)/d(theta_vmec,phi_vmec) } / (Jacobian_vmec/Jacobian_booz)
  double Vmec2BoozerTest(const Vector3d &vmecCoord);
  //@}
  /// \name Inquire position
  //@{
  /// The method tests whether the point lies inside LCMS
  /// @param cyl -- cylindrical coordinates of the point.
  /// @param distance if not NULL then this address at return
  ///        will contain the distance between cyl and the LCMS,
  ///        distanse<0 if the point lies inside LCMS.
  /// @return  the return value is \b true if the point lies inside LCMS; \b false otherwise.
  bool cylIsInside(const Vector3d &cyl,double * distance=0) const;
  /// The method tests whether the point lies inside LCMS
  /// @param xyz -- cartesian coordinates of the point.
  /// @param distance if not NULL then this address at return
  ///        will contain the distance between xyz and the LCMS,
  ///        distanse<0 if the point lies inside LCMS.
  /// @return  the return value is \b true if the point lies inside LCMS; \b false otherwise.
  bool xyzIsInside(const Vector3d &xyz,double * distance=0) const;
  /// The method tests whether the point lies inside the given border.
  /// @param cyl -- cylindrical coordinates of the point.
  /// @param distance if not NULL then this address at return
  ///        will contain the distance between cyl and the border,
  ///        distanse<0 if the point lies inside the border.
  /// @param point if not NULL then this address at return
  ///        will contain the cylindrical coordinates of nearest point on the border.
  /// @param magCoord if not NULL then this address at return
  ///        will contain the magnetic coordinates of nearest point on the border.
  /// @return true if the point lies inside the border; \b false otherwise.
  ///  \note createBorder() method must be called before using this method.
  bool  cylIsInsideBorder(const Vector3d &cyl,double *distance=0,Vector3d *point=0,Vector3d *magCoord=0) const;
  /// The method creates the border -- 's'-contour at cylindrical angle 'ficyl'.
  /// @param s is the value of flux label
  /// @param ficyl is the cylindrical angle
  /// @param borderSize is the number of point in contour
  /// @return false if no memory
  bool  createBorder(double s,double ficyl,int borderSize=150);
  /// The method returns the number of points in the border created by createBorder() method.
  int   getBorderSize() { return mBorder.size();}
  /// The method returns the cylindrical cooordinates of the k-point of the border.
  /// The border is created by createBorder() method and
  /// 0<=k<getBorderSize().
  const Vector3d & getBorder(int k) const { return mBorder[k];}
  /// The method returns the magnetic cooordinates of the k-point of the border.
  /// The border is created by createBorder() method and
  /// 0<=k<getBorderSize().
  const Vector3d & getBorderMag(int k) const { return mBorderMagCoord[k];}
  /// The method returns the min value of R-cooordinate of the border.
  const double & getBorderRmin() const {return mBorderRmin;}
  /// The method returns the max value of R-cooordinate of the border.
  const double & getBorderRmax() const {return mBorderRmax;}
  const double & getBorderZmin() const {return mBorderZmin;}
  const double & getBorderZmax() const {return mBorderZmax;}
  /// The method traces the Last Closed Magnetic Surface (LCMS).
  /// The method finds the first intersection of the ray with the LCMS, where the
  /// ray is R(t)=r0+rd*t,  and t>=0.
  /// @param r0 is the Cartesian coordinates of the ray origin, must be outside with respect to the LCMS,
  /// @param rd is the direction of the ray in Cartesian coordinates,
  /// @param entryPoint this vector at return contains the Cartesian coordinates of the ray entry point into the plasma,
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool getRayEntryPoint(const Vector3d &r0,const Vector3d &rd,Vector3d &entryPoint);
  /// The method traces the Last Closed Magnetic Surface (LCMS).
  /// The method finds intersections of the ray with the LCMS, where the
  /// ray is R(t)=r0+rd*t,  and t>=0.
  /// @param r0 is the Cartesian coordinates of the ray origin, must be outside with respect to the LCMS,
  /// @param rd is the direction of the ray in Cartesian coordinates,
  /// @param entryPoint this vector at return contains the Cartesian coordinates of the ray entry point into the plasma,
  /// @param exitPoint this vector at return contains the Cartesian coordinates of the ray exit point from the plasma,
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool getRayIntersectionPoints(const Vector3d &r0,const Vector3d &rd,Vector3d &entryPoint,Vector3d &exitPoint);
  Vector3d xyzAtB(double Bres,const Vector3d &r0,const Vector3d &rd, const Vector3d &entryPoint);

  ///****************************************************************************
  /// The method returns the coefficients needed for 
  /// the poloidal flux equation in Astra transport code
  /// @param[in] sqrts is the sqrt(s), where s is the normalized toroidal flux 
  /// @param[out] r = a*sqrt(s) is the effective plasma radius, where a is the minor radius.
  /// @param[out] gradr2Avr is the flux surface average of (grad(r))^2, \<(grad(r))^2\> \f$\left\langle\left|\nabla r\right|^2\right\rangle\f$
  /// @param[out] J = Ipol/Ipol(a)
  /// @param[out] G2 = h*S11 > 0 is the G2 coefficient in the poloidal flux equation,
  ///       h is +1 for the right handed (s,theta,phi) coordinate system and is -1 for the left handed one.
  /// @param[out] hV' = h*dV/dr > 0
  /// @param[out] B0 = h*TorFlux(a) / (pi*a^2)
  /// @param[out] R0 = mu0*Ipol(a) / (2pi*B0)
  void getCoeffForAstraCode(double sqrts, double &r, double &gradr2, double &J, double &G2, double &hV, double &B0, double &R0, double &h);

  //@}
  /// \name Statistics
  //@{
  /// Clear statistics
  void zeroStatistics();
  /// The method returns the number of NJacobian calls
  ///  \note for benchmarking and timing.
  int & nJac() { return mNjac;} // number of NJacobian calls
  /// The method returns the number of coordinate transformations from real space to magnetic coordinates.
  ///  \note for benchmarking and timing.
  int & n2mag() { return mNcyl2mag;} // number of cyl2bmag calls
  /// The method returns the average number of poloidal harmonics used in summation during
  ///  coordinate transformations from real space to magnetic coordinates.
  ///  \note for benchmarking and timing.
  double MAverage() const;
  /// The method returns the average number of toroidal harmonics used in summation during
  ///  coordinate transformations from real space to magnetic coordinates.
  ///  \note for benchmarking and timing.
  double NAverage() const;
  /// The method returns the number of NJacobian calls.
  ///  \note for benchmarking and timing.
  double NNjacCalls() const;
  /// The method returns the number of mag2cyl calls.
  ///  \note for benchmarking and timing.
  double Nmag2cylCalls() const;
  /// The method returns the number of cyl2mag calls.
  ///  \note for benchmarking and timing.
  double Ncyl2magCalls() const;
  double NGuessCalls() const;

  //@}
  /// \name Misc
  //@{
// Find magnetic toroidal angle phi for point (s,theta,ficyl), ficyl is a cylindrical angle
  /// The method returns the magnetic toroidal angle which corresponds to the mixed coordinates.
  /// @param s is the flux surface label
  /// @param theta is the magnetic poloidal angle.
  /// @param ficyl is the cylindrical angle.
  /// @return magnetic toroidal angle.
  double phimag(double s,double theta,double ficyl) const;
  /// The method returns the magnetic toroidal angle which corresponds to the mixed coordinates.
  /// \overload
  /// @param mixcoord is the Vector3d(s,theta,ficyl), see phimag(double s,double theta,double ficyl)
  /// @return magnetic toroidal angle.
  double phimag(const Vector3d &mixcoord) const;
  /// The method returns the sizes of a rectangle circumscribed about the cross
  /// section of the surface \b s at cylindrical angle \b fi.
  void getExtent(double &Rmin,double &Rmax,double &Zmin,double &Zmax,double s,double fi) const;
  /// The method returns the sizes of a rectangle circumscribed about the cross section of the surface \b s.
  void getExtent(double &Rmin,double &Rmax,double &Zmin,double &Zmax,double s) const;
  /// The method returns the sizes of a rectangle circumscribed about the cross section of the LCMS.
  void getExtentLCMS(double &Rmin, double &Rmax, double &Zmin, double &Zmax) const;
  /// The method returns the sizes of parallelepiped circumscribed about the segment f1-f2<pi
  /// of the 's' surface, 'f' is the cylindrical angle
  void getExtent(Vector3d &rmin, Vector3d &rmax, double s, double f1, double f2) const;
  /// The method returns the sign of the Jacobian calculating the mixed product dX/ds*(dX/dtheta^dX/dphi).
  /// Where ^ is a cross product signs,and the (s,theta,phi) are the magnetic coordinates:
  /// (flux label, poloidal, and toroidal angles); X is the cartesian vector.
  int getJacobianSign();
  /// The method returns the direction of the poloidal angle.
  /// Plus sign corresponds to counterclockwise direction on the R-Z plane
  /// when looking on it in fi-direction,where fi is the cylindrical angle.
  /// Plus sign also means that increasing theta increases angle = atan((Z(theta)-Zaxis)/(R(theta)-Raxis)).
  int getThetaDirection();
  /// The method returns the direction of the toroidal angle.
  /// Plus sign corresponds to the direction on the cylindrical angle.
  int getPhiDirection();
  /// The method returns the direction of the toroidal current j^phi (in the right-handed
  /// cylindrical coordinates system), that increases the absolute value of the iota.
  /// @return 1 or -1 depending on the direction of the toroidal current.
  int getSignOfCurrentThatIncreasesAbsIota() const {return (mSignJac_s*Ip(0.25)*iota(0.25) >= 0)?1:-1;}
  /// The method returns the sign of iota increment due to non-inductive current \<j*B\>/\<|B|\>
  /// @param Itor is the current that is equal to integral_0^r \<j*B\>/\<|B|\>  *(2*pi*R0)^-1 *dV
  /// @return 1 or -1 depending on sign of the iota increment;  0 if current is zero.
  int getDeltaIotaSign(double Itor) const {return Itor==0?0:( (mSignJac_s*Itor > 0)?1:-1 );}
#if 0 
  /// The method transform the current from cylindrical coordinates system to magnetic coordinates.
  /// @param Ifi is the toroidal current in cylindrical system
  /// @return toroidal current in magnetic coordinates system.
  double Icyl2Itor(double Ifi) const {return mSignJac_s*Ifi;}
  /// The method transform the current from magnetic coordinates system to cylindrical coordinates.
  /// @param Itor is the toroidal current in magnetic coordinates system
  /// @return toroidal current in the right-handed cylindrical coordinates system.
  double Itor2Icyl(double Itor) const {return mSignJac_s*Itor;}
#endif
  /// The method returns the pointer to the 2d-array of vertexes of LCMS.
  /// @return the pointer to the 2d-array of vertexes of LCMS, the vertexes are in cartesian coordinates.
  const CArray2d<Vector3d> * LCMS();
  /// The method returns the pointer to the 2d-array of the magnetic field value at the LCMS.
  const CArray2d<Vector3d> * BLCMS();
  /// The method reduces \em x to the interval 0..2pi.
  double mod2pi   (double x) const { return (x>=0)?fmod(x,twopi):(fmod(x,twopi)+twopi);}
  /// The method reduces \em x to the interval 0..2pi.
  ///  if 0=<x<=twopi then x is returned,
  ///  otherwise (x<0)?(fmod(x,twopi)+twopi):fmod(x,twopi)   
  double mod2piNoWrap   (double x) const 
  {
    if(x<0) {
      x = fmod(x,twopi)+twopi;
      ////int ia=int(x*twopiinv); return x-(--ia)*twopi;
    }
    else if(x>twopi) {
      x = fmod(x,twopi);
      //////int ia=int(x*twopiinv); return x-ia*twopi;
    }
    return x;
  }
  /// The method reduces \em x to the interval 0..2pi/mNp, where mNp is the number of stellarator periods.
  double modPeriod(double x) const { return (x>=0)?fmod(x,mPeriod):(fmod(x,mPeriod)+mPeriod);}
  /// The method sorts angles, so that f1<f2 and f2-f1<pi
  void sortAngles(double &f1, double &f2) const;
  /// The method creates array of vertexes of the surface s.
  bool createSurface(double s,int nTheta,int nPhi,CArray2d<Vector3d> &vertex,CArray2d<Vector3d> *Bfld=0, bool cylPhi=true);
  /// The method creates array of normals to the surface Vertexes.
  /// @param Vertexes -- array of vertexes of a surface.
  /// @param sign -- direction of normals
  /// @param normals -- array of normals
  /// @return normals -- array of normals
  void createNormals(const CArray2d<Vector3d> &Vertexes, CArray2d<Vector3d> &normals, int sign=1) const;
  /// The method calculates curvature at given point.
  /// @param magCoord is the magnetic coordinates of the point of interest.
  /// @return Vector3d with component:
  ///  \li [0] - normal curvature
  ///  \li [1] - geodesic curvature
  ///  \li [2] - Jacobian in coordinates (s,theta,phi)
  Vector3d curvature(const Vector3d &magCoord);
  //@}
  //   0.01227 = twopi/512
  //Vector3d bscurrentFactor(double s, double dphi=0.01227,int numBlevels=513,int numTurns=100);

private:
  uiMutex mConfMutex;
  typedef double REAL;
  bool getRayIntersection(Vector3d r, Vector3d dr, Vector3d &intersectionPoint);
  bool isRationalIota(double s, double eps=2e-2);
  void advanceDataForEpsEff02(Vector3d &mag, Vector3d &r1, Vector3d &B1, Vector3d &n1, 
                             REAL &Bmod, REAL &grads, REAL &kG, REAL &dl, Vector3d *curvature,
                             const double &dphi, const double &iotaRatnl);
  void advanceDataForEpsEff(Vector3d &mag, Vector3d &r, double &B, double &gradReff, double &kG, const double &dphi);
  void advanceDataForJinvariant(Vector3d &mag, Vector3d &r1,REAL &Bmod,Vector3d &covaGradB, REAL &dl, const double &dphi, const double &iotaRatnl);
  double getOmega_d(Vector3d & mag1, double Bma, Vector3d & curvature1,Vector3d &dBdu1);

  // dphi=0.00698  == 0.4 degree
  void epsEff(double s, __allEps &eps);
  void calculateEpsEff();
  void calculateEpsEffExe(CStconfig * objR, int low, int upper, bool setup);

  void calculateFbs();
  void calculateFbsExe(CStconfig * objR, int low, int upper, bool setup);
  void bscurrentFactor(double s, __allFbs &retValue);

  void calculateGradsAvrg(int size=34);
  void calculateGradsAvrgExe(CStconfig * objR, int low, int upper, bool setup);

  void calculateSMatrix(int size = 30);
  void calculateSMatrixExe(CStconfig * objR, int low, int upper, bool setup);

  /// The methods calculates min and max B on surfaces.
  void calculateBminmax();
  void calculateBminmaxExe(CStconfig * objR, int low, int upper, bool setup);

  void calculateFtrappedBminmax();
  void calculateFtrappedBminmaxExe(CStconfig * objR, int low, int upper, bool setup);

  void calculateIpolItorFromSij(int ntheads=-1);
  void calculateIpolItorFromSijExe(CStconfig * objR, int low, int upper, bool setup);

  void calculateIpolItor(int ntheads=-1);
  void calculateIpolItorExe(CStconfig * objR, int low, int upper, bool setup);
  double calculateItor(int is, int Ntheta=180);
  double calculateIpol(int is, int Ntphi=180);


  void calculateVp(int ntheads=-1);
  void calculateVpExe(CStconfig * objR, int low, int upper, bool setup);
//template <class TR> void calculateVpExe(TR * objR, int low, int upper, bool setup);

  void calcEfitFourierMods(void * efit);
  void calcEfitFourierModsExe(CEfit* ef,int low, int upper, bool setup);
  //template <class efit> void calcEfitFourierModsExe(efit* ef,int low, int upper, bool setup);

  /// Get memory
  /// @return  \em true if the method succeeds; \em false otherwise.
  bool resize(int ns,int m,int n);
  /// function returns fi -- cylindrical angle and dfi/dphi for given point in magnetic coordinates.
  /// fi = phi - (2*pi/Np)*sum_mn Phi_mn(is)*sin(m*theta-n*Np*phi)
  /// @param[in] s,theta,phi -- magnetic coordinates.
  /// @param[out] dfi_dphi is the returned value of dfi/dphi
  /// @param[in] onFildLine derivative along field line if true
  /// @return dfi_dphi -- the value of dfi/dphi
  /// @return cylindrical angle.
  double fi(double s,double theta,double phi,double &dfi_dphi,bool onFildLine=false) const;
  /// \overload
  double fi(const Vector3d &magCrd,double &dfi_dphi,bool onFildLine=false) const;
  /// \overload
  double fi(double s,double costheta,double sintheta,double cosNpPhi,double snNpPhi,
    double phi,double &dfi_dphi,bool onFildLine=false) const;
  /// Ascertain the file format.
  /// @return  fileFormat, see CStconfig::fileFormat.
  fileFormat getfileformat(const char * fullname,const char *nam,long fPostn) const;
protected:
  /// helper function for CStconfig::load()
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool loadfile (const char * filename,fileFormat t, long fPostn=0);
private:
  /// helper function for CStconfig::load()
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool loadascii(const char * filename,fileFormat fmt, long fPostn=0);
  /// helper function for CStconfig::write()
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool writeascii(const char * filename) const;

  bool writeBoozerCoord(const void *mcdb, int magneticgeometryrunlogid, int Nsmax,bool magAxis) const;

  /// The method loads VMEC wout-file and transform it to boozer.
  /// helper function for CStconfig::load()
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// \note must be revised
  bool loadwout (const char * filename);
  /// The method loads VMEC wout-file.
  /// helper function for CStconfig::load()
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool loadvmec (const char * filename,long fPostn=0);
  bool loadnemec(const char * filename,long fPostn=0);
  /// The method loads VMEC netCDF wout-file.
  /// helper function for CStconfig::load()
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
#ifdef NETCDF
  bool loadvmecnetcdf(const char * filename,long fPostn=0);
#endif
  /// The method loads LHD boozer-file.
  /// helper function for CStconfig::load()
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool loadlhd (const char * filename);
  /// The method loads magnetic configuration stored in the EQDSK file \c fname.
  /// The method understands format of a file by analyzing its contents.
  /// \sa the detailed description of CStconfig::load().
  /// The method import equilibrium from the EQDSK File into straight-field-line coordinates 
  /// system (s,theta,phi) with the cylindrical angle phi as the toroidal angle. 
  /// The right handed cylindrical coodinate system (r,phi,z) is used. The system (s,theta,phi) can be left or right handed; 
  ///  the direction of theta and correspondingly the Jacobian sign depends on Bpol, Btor directions and the sign of safety factor, see below.
  ///
  ///  \n The value of the following parameters BpolScale, BtorScale, BpolSign, BtorSign, qSign, psiOverTwopi 
  ///  can be put in the first line of the EQDSK file, for example 
  /// \code
  /// RC-ITER...             06-DEC-02                  3 129 129  EQDSK   BpolScale=0 BtorScale=0 BpolSign=1 BtorSign=1 qSign=1 psi/2pi=no    
  /// \endcode
  /// @param[in] fname is the name of G EQDSK file
  /// @param[in] BpolScale is the scale factor for the poloidal field, the factor sign is ignored, 
  ///            BpolScale=0 means that the scale factor is taken from the the first line of the EQDSK file 
  ///            if it found, otherwise no scaling is performed.
  /// @param[in] BtorScale is the scale factor for the toroidal field, the same logic as for poloidal field scaling.

  /// @param[in] BpolSign sets the sign of the poloidal field (by changing sign of the poloidal flux psi):
  ///    \li     +1 means positive Bpol, positive Bpol is created by the positive currrent density j_phi flowing in direction of the angle phi 
  ///             in the right handed cylindrical coodinate system (r,phi,z);
  ///    \li     -1 means negative Bpol created by the negative currrent density flowing in counter direction of the cylindrical angle phi; 
  ///    \li      0 means that Btor direction is taken from the BpolSign=n (where n is one from =-1,0,1) if it stored in the first line of the EQDSK file,
  ///    \li     if BpolSign not found in the file or BpolSign=0 then Bpol direction is defined by the psirz array stored in the EQDSK file.
  ///    \li       B = (-dpsirz/dz, Fpol, dpsirz/dr)/r 

  /// @param[in] BtorSign sets the sign of the toroidal field:
  ///    \li     +1 means positive Btor in direction of the angle phi in the right handed cylindrical coodinate system (r,phi,z),
  ///    \li     -1 means negative Btor in counter direction of the cylindrical angle phi, 
  ///    \li      0 means that Btor direction is taken from the BtorSign=n (where n is one from =-1,0,1) if it stored in the first line of the EQDSK file,
  ///    \li     if BtorSign not found in the file or BtorSign=0 then Btor direction is defined by the poloidal current function stored in the EQDSK file.
  ///    \li       B = (-dpsirz/dz, Fpol, dpsirz/dr)/r 

  /// @param[in] qSign sets the sign of the safety factor, qSign is taken from the EQDSK file if this parameter is zero, 
  ///            if qSign not found in the file or qSign=0 in it then the safety factor is defined by the q values stored in the EQDSK file.

  /// @param[in] psiOverTwopi defines whether the poloidal flux psi is divided by 2pi:
  ///    \li      -1 means no, i.e psi is not devided over 2pi, this corresponds to COCOS>=11.
  ///    \li      0,1 assume that the flux psi stored in the file is already divided by 2pi.
  ///    \li      in case psiOverTwopi is 0 the psiOverTwopi is derived from the string psi/2pi=val stored in the EQDSK file hader, where val is yes or no 

  /// @param[in] mconf if it not zero then the circular tokamak with the same R0, a, and iota as in the mconf is created. 
  /// @param[in] fPostn is the position in the file fname from where to start reading EQDSK staff. 
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.

  bool loadEfit_p(const char * fname, double scaleBpol=0, double scaleBtor=0,int signBpol=0,int signBtor=0, 
                  int signQ=0, int psiOverTwopi=0, const CStconfig *mconf=0, long fPostn=0);
public:
  /// The method loads the \ref GEQDSK "G EQDSK file" 
  /// and transforms it to the \ref TokamakSymmetryFluxCoordinates "tokamak-symmetry flux coordinates".
  /// @param[in] fname is the name of the file to load.
  /// @param[in] scaleBpol is the scale factor for the poloidal field.
  /// @param[in] scaleBtor is the scale factor for the toroidal field.
  /// @return  if the function succeeds, the return value is \b non-zero, \b zero otherwise.
  /// \sa loadEfit_p()
  bool loadEfitAndScale(const char * fname, double scaleBpol=0, double scaleBtor=0,int signBpol=0,int signBtor=0,int signQ=0,int psiOverTwopi=0);
  /// The method creates the circular tokamak with the same R0, a, and iota as in the mconf.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool loadTokamak(const CStconfig *mconf);
protected:
  /// The method loads torbeam equilibria file.
  /// helper function for CStconfig::load(), not implemented yet.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool loadTopfile(const char * fname) {fname; return mOK=false;}
  /// Adjust number of integration points and steps.
  void adjustNumbersForSimpson(int &Npol, int &Ntor, double &dpol, double &dtor) const;
  /// The method calculates currents and dV/ds if needed.
  void calculateCurrents(); 
  /// The method calculates dV/ds and then Volume integrating Jacobian  
  void calculateVprime();
  void calculateFluxAverages1(int npol=180, int ntor=200);
  /// The method tabulate arrays for finding Sij(s).
  /// \note This method is obsolete.
  void calculateSMatrixArrays(int npol=180, int ntor=200);
  /// ajust Fourier coefficients and create radial meshes
  void ajustFourierCoeffMeshes();
  /// The method adds surfaces near magnetic axis and calculates spline coefficients
  void addSurfacesNearAxis();
  void addSurfacesAfterLCMS();
  void createSplineCoeff();
  /// The method searches nearest to \b s segment in array ms
  /// @return  nearest to 's' index in array 'ms' and
  /// numbers of poloidal \b M and toroidal \b N harmonics for this surface.
  /// \b M, \b N are based on truncation level, see truncate(double level)
  int  getSsegMN(double s, int &M, int &N) const;
  /// The method returns the Volume inside surface \b s integrating Jacobian.
  /// The Volume is negative if (s,theta,phi)-coordinate system is left handed.
  double VolumeInt(double s) const;              //Volume(s), [m^3]
  /// The method returns the volume between surfaces \b s1 and \b s2 integrating Jacobian.
  double VolumeInt(double s1, double s2) const;  //Volume(s2)-Volume(s1)
  /// The method creates array of vertexes of LCMS.
  bool createLCMS();
  /// The method calculates the sign of Jacobian and angle directions
  void signOfJacobianAndAngles();
  /// The method is used to create flux label based on normalized poloidal flux.
  double FluxPoloidal(double s1, double s2) const;

  /// \internal
  //@{
private:
  void V2Btrans(Vector3d &boozer, double &R, double &Z, double &p, double &B);
  void V2Btrans(const Vector3d &vmec, Vector3d &boozer, double &R, double &Z, double &p, double &B, double &jv_jb);
  void FourierUsingBoozerGrid(const CComplex &ctm, const CComplex &cpn, 
                const CArray2d<double> &R, const CArray2d<double> &Z, 
                const CArray2d<double> &P, const CArray2d<double> &B,                      
                double &Rmn,double &Zmn,double &Pmn,double &Bmn); 
  void FourierUsingVmecGrid(int m, int n, const CArray2d<Vector3d> &booz, 
                 const CArray2d<double> &R, const CArray2d<double> &Z, 
                 const CArray2d<double> &P, const CArray2d<double> &B,                      
                 double &Rmn,double &Zmn,double &Pmn,double &Bmn); 
  void createSurfaceDataOnBoozerGrid(double s, 
                         CArray2d<double> &R, CArray2d<double> &Z, 
                         CArray2d<double> &P, CArray2d<double> &B); 
  void createSurfaceDataOnVmecGrid(double s, CArray2d<Vector3d> &booz, 
                         CArray2d<double> &R, CArray2d<double> &Z, 
                         CArray2d<double> &P, CArray2d<double> &B); 
  void calcBoozFourierMods(int Ns, int M, int N, int method);
  void calcBoozFourierModsExe(CStconfig * objR, int low, int upper, bool setup);


  Vector3d *NJacobian1L(const Vector3d &magCoord);
  Vector3d *NJacobian1Q(const Vector3d &magCoord);
  void      NJacobian2 (int typeOfInterp);
  //@}
  /// \name Internal things for accesing Fourier coefficients and spline coefficients
  //@{
  // Fourier coefficients
  /// The method returns Fourier coefficient of |B|-expansion.
  /// @param is is the flux surface number, 0<=is< nSurfaces().
  /// @param m is the poloidal mode number, 0<=m<= Mpol().
  /// @param n is the toroidal mode number. |n|<=Ntor().
  /// @return Fourier coefficient B_mn(is).
  double  Bmn(int is, int m, int n) const {return mCOEF[(n+mN)*offsetN + m*offsetM + is*offsetS]*fabs(mBnorm); }

#define FC(Cmn,i)   double  Cmn(int is,int m,int n) const {return mCOEF[i+(n+mN)*offsetN + m*offsetM + is*offsetS];} \
                    double& Cmn(int is,int m,int n)  {return mCOEF.rw()[i+(n+mN)*offsetN + m*offsetM + is*offsetS];}

#define FC_d2(Cmn,i)   double  Cmn##_d2(int is,int m,int n) const {return mCOEF[i+(n+mN)*offsetN + m*offsetM + is*offsetS];} \
                       double& Cmn##_d2(int is,int m,int n)  {return mCOEF.rw()[i+(n+mN)*offsetN + m*offsetM + is*offsetS];}

  FC(bmn,0) FC_d2(bmn,1)
  /// The method returns Fourier coefficient of R-expansion, for parameters see Bmn()
  FC(Rmn,2) FC_d2(Rmn,3)
  /// The method returns Fourier coefficient of Z-expansion, for parameters see Bmn()
  FC(Zmn,4) FC_d2(Zmn,5);
  /// The method returns Fourier coefficient of Phi-expansion, for parameters see Bmn()
  FC(Phimn,6) FC_d2(Phimn,7)
  /// The methods return Fourier coefficient of ... on half-mesh, for parameters see Bmn().  
  /// VMEC
  FC(lmn,8) FC_d2(lmn,9)
  FC(gmn,10) FC_d2(gmn,11)
  FC(bContraTheta_mn,12) FC_d2(bContraTheta_mn,13)
  FC(bContraFi_mn,14)   FC_d2(bContraFi_mn,15)

 // VMEC additional fourier coefficients in asymmetrical mode  if vmecAsym is true

  FC(bmns,16) FC_d2(bmns,17)
  /// The method returns Fourier coefficient of R-expansion, for parameters see Bmn()
  FC(Rmns,18) FC_d2(Rmns,19)
  /// The method returns Fourier coefficient of Z-expansion, for parameters see Bmn()
  FC(Zmnc,20) FC_d2(Zmnc,21);
  /// The method returns Fourier coefficient of Phi-expansion, for parameters see Bmn()
  FC(Phimns,22) FC_d2(Phimns,23)
  /// The methods return Fourier coefficient of ... on half-mesh, for parameters see Bmn().  
  /// VMEC
  FC(lmnc,24) FC_d2(lmnc,25)
  FC(gmns,26) FC_d2(gmns,27)
  FC(bContraTheta_mns,28) FC_d2(bContraTheta_mns,29)
  FC(bContraFi_mns,30)    FC_d2(bContraFi_mns,31)


  void calcVmecLgB(const Vector3d &u, Vector3d * Bcontrav=0, bool boozTrans=false);
  double  JacobianVmec(const Vector3d &magCoord);
  int vmecErrorCode;
  bool useVmecLambda;
  //@}
  /// \name Internal stuff for extrapolation, interpolation, splining.
  //@{
  double extrapolate0(double u, const double *x, const double *y) const;
  double extrapolateLCMS(double u, const double *x, const double *y,int n,int h=1) const;
  //double S_interp  (double s, const double *y, const double *y_d2) const;
  //double S_interp_D(double s, const double *y, const double *y_d2) const;
  double Linterp(const double *y)  const;  // Linear interpolation
  double interp2(double u, const CArray1d<double> &x, const CArray1d<double> &y) const;
  double interp1(double u, const CArray1d<double> &x, const CArray1d<double> &y) const;

   /// internal structure holds things for splining
  struct splineCoeffs {
    double e2rcurr;    //= 0.5/sqrt(s)
    double eds01,ds0,ds0eds01; // data for Linterp (linear interpolation)
    int is0;                   // current segment in ms[] or in msFull[]
    // coeffs for spline follow
    double eh;
    double a;
    double b;
    double c;
    double d;
    double c1;
    double d1;
  };

  splineCoeffs halfMesh;
  splineCoeffs fullMesh;

  ///  \brief This class is for internal use only
  class sumDat { // data to speed up summation in cyl2guess
  public:
    double s;
    double fi;
    double csNp;
    double snNp;
    const double * fCh;
    const double * fCf;
    int M,N;
    splineCoeffs halfMesh;
    splineCoeffs fullMesh;
    splineCoeffs * hMesh;
    splineCoeffs * RZmesh;
    sumDat() { };
    sumDat(CStconfig &thiz, double s, double fi) {set(thiz,s,fi); };
    void set(CStconfig &thiz, double s, double fi) {
      this->fi = fi;
      set(thiz, s);
      csNp = ::cos(thiz.mNp*fi);
      snNp = ::sin(thiz.mNp*fi);
    }
    void set(CStconfig &thiz, double s) {
      this->s = s;
      thiz.getSsegMN(s, M, N); // search nearest segment in array ms;
      halfMesh = thiz.halfMesh;
      fullMesh = thiz.fullMesh;
      fCh = thiz.mCOEF.constArray()+thiz.mN*thiz.offsetN+ halfMesh.is0*thiz.offsetS;
      fCf = thiz.mCOEF.constArray()+thiz.mN*thiz.offsetN+ fullMesh.is0*thiz.offsetS;
      hMesh  = &halfMesh;
      RZmesh = thiz.vmecCoordinates? &fullMesh:&halfMesh;
   }
  };

  //@}
//arrays for flux surface averages
  int mNsK;           ///< number of flux surface used for tabulating of mBmin, Bmax,...
  CArray1d<double> mreK;       ///< mreK[mNsK]--array of normalized r=sqrt(s),
  CArray1d<double> mBmin;      ///< Bmin on a surface, array of length mNsK
  CArray1d<double> mBmax;      ///< Bmax on a surface, array of length mNsK
  CArray1d<double> mFtrap;     ///< fraction of trapped particles, array of length mNsK
  CArray1d<double> mBavrg;     ///< B average on a surface, array of length mNsK
  CArray1d<double> mB2avrg;    ///< B^2 average on a surface, array of length mNsK
  CArray1d<double> mR2avrg;    ///< R^2 average on a surface, array of length mNsK
  CArray1d<double> mB_2avrg;   ///< B^-2 average on a surface, array of length mNsK
  CArray1d<double> mR_2avrg;   ///< R^-2 average on a surface, array of length mNsK
  CArray1d<double> mGrads2B_2avrg;   ///< grads^2*B^-2 average on a surface, array of length mNsK
  CArray1d<double> mNeoPolarztnPass;   ///< for Mishchenko, array of length mNsK
  CArray1d<double> mNeoPolarztnTrap;   ///< for Mishchenko, array of length mNsK
  CArray1d<double> mb12;    // <(1-b)^0.5>  , array of length mNsK
  CArray1d<double> mb32;    // <(1-b)^1.5>  , array of length mNsK

  double mBminGlobal;           ///< min of mBmin
  double mBmaxGlobal;           ///< max of mBmax

  int mNxK;           // must be odd!!
  CArray1d<double> mxtr;   // = x[i] for trapped particles fraction calculation,array[mNxK]
  CArray1d<double> mpavr2; // = (<sqrt(1-x[i]*b)>)^2 for given surface, array[mNxK]
  CArray1d<double> mpavrI; // = integral(from x to 1) (dx/<sqrt(1-x[i]*b)>) for given surface,array[mNxK]
  CArray1d<double> mpavrb2;// = (<b/sqrt(1-x[i]*b)>)^-2 for given surface,array[mNxK], for Mishchenko
#ifdef __INTEL_COMPILER
// These pragmas are to quiet icc about "virtual function override intended?"
#pragma warning( push )
#pragma warning( disable : 1125 )
#endif
  CArray2d<double> mpavr_sx;  // f(s,x)= square(<sqrt(1-x*b(s))>), array[mNsK][mNxK], see also CStconfig::pavrg
  CArray2d<double> mpavrI_sx; // integral(from x to 1) dx/<sqrt(1-x*b(s))> , array[mNsK][mNxK]
  CArray2d<double> mpavrb_sx; // f(s,x)= (<b/sqrt(1-x*b(s))>)^-2, array[mNsK][mNxK], see also CStconfig::roParInvAvrg, for Mishchenko
#ifdef __INTEL_COMPILER
#pragma warning( pop )
#endif

//arrays for Smatrix
  CArray1d<double> msK;    ///< msK --array of normalized tor. flux s
  CArray1d<double> mS11;   ///< array of length msK.size()
  CArray1d<double> mS12;   ///< array of length msK.size()
  CArray1d<double> mS21;   ///< array of length msK.size()
  CArray1d<double> mS22;   ///< array of length msK.size()

Vector3d getGradsAverages(double s,int ithmax=180,int ipmax=200);

void getFluxAverages(double s,double &Ftrapped,
                                  double &Bmin,double &Bmax,
                                  double &Bavr,
                                  double &B2avr,  // <B^2>
                                  CArray1d<double> &mpavrI,
                                  CArray1d<double> &mpavr2,
                                  int itmax,int ipmax);
// for A.Mishchenko
void getFluxAverages1(double s,
                                  double &B_2avr, // <B^-2>
                                  double &R2avr,  // <R^2>
                                  double &R_2avr,  // <R^-2>
                                  double &g2B_2avr, // <grads^2*B^-2>
                                  double &b_12,   // <(1-b)^0.5>
                                  double &b_32,   // <(1-b)^1.5>
                                  double &neoPolarPass,   // for Mishchenko
                                  double &neoPolarTrap,   // for Mishchenko
                                  int itmax,int ipmax);
protected:
  bool jacobianIsMixedProduct;
  bool Averages1Ready;
  bool ftrappedReady;
  bool BminmaxReady;
  bool SMatrixReady;
  bool tokamakSymmetryFluxCoordinates;
  bool vmecCoordinates;   ///< VMEC Coordinates if true
  bool vmecAsym;          ///< VMEC additional fourier coefficients in asymmetrical mode  if true
  bool fpCoordinates;     ///< VMEC Coordinates from FP if true
  bool fpCoordinates2;    
  bool fpCoordinates5;    
  bool tokamakEQDSK;      ///< same as tokamakSymmetryFluxCoordinates
  bool tokamakConfig;     ///< mN==0 && mNp==1
  bool polCurrentNeeded;  ///< signal to load() function that currents must be calculated
  bool torCurrentNeeded;  ///< signal to load() function that currents must be calculated
  bool vprimeNeeded;      ///< signal to load() function that dV/ds must be calculated
  /// init everything before loading magnetic configuration
  void init();
  /// the method is called when the configuration is reloaded 
  /// to another coordinate system (L.H.S.  to/from  R.H.S) 
  void clearImportantFlags();
  /// Free all memory used.
  void freeConf();
  /// write identification and comments into bc-file.   
  void writeIdentificationComment(FILE * fp) const;
  /// return the sign of B-normalization factor
  int getSignOfBnormFactor() const;
  /// helper function for write(), see also C3dMesh::writeMeshbin();
  template <class T> bool writebin(const char * fname, T szType, FILE *ff=NULL) const;
  /// helper function for write(), see also C3dMesh::writeMeshbin();
  bool writebin4(const char * fname, FILE *ff) const;
  /// helper function for write(), see also C3dMesh::writeMeshbin();
  bool writebin8(const char * fname, FILE *ff) const;
  /// helper function for load(), see also C3dMesh::loadMeshbin();
  template <class T> bool loadbin (const char * fname, T szType, FILE *ff=NULL);
  /// helper function for load(), see also C3dMesh::loadMeshbin();
  bool loadbin4(const char * fname, FILE *ff);
  /// helper function for load(), see also C3dMesh::loadMeshbin();
  bool loadbin8(const char * fname, FILE *ff);
  /// calculate spline coeffitients and add surfaces near axis.
  bool initAfterLoad(fileFormat fmt=UNDEFINED);
  /// test whether point \b ri lies inside triangle \b r
  int isInsideTriangle(const Vector3d *r,const Vector3d &n, Vector3d &ri) const;
  int isTriangleHitted(const Vector3d *r,const Vector3d &rayR0) const;
  /// Returns a pointer to the last period in buf
  const char *fname_ext(const char *buf) const;
  /// returns a pointer to the filename, or null if name ends with '/'
  const char *fname_name(const char *name) const;
  /// The method tests whether the point lies inside LCMS
  /// @param cyl -- cylindrical coordinates of the point, fi=cyl[1] must be given on a period.
  /// @param vertex
  /// @param dFi  -- toroidal angle step, which was used for tabulating array vertex
  /// @param distance if not NULL then this address at return
  ///        will contain the distance between cyl and the LCMS.
  /// @return  the return value is \b true if the point lies inside LCMS; \b false otherwise.
  bool cylIsInsideEx(const Vector3d &cyl,const CArray2d<Vector3d> &vertex, double dFi,double * distance) const;
  
  template <class T> void Swap  (T &v1, T &v2) const {T w = v1;  v1 = v2;  v2 = w;}
  template <class T> void SwapN (T &v1, T &v2) const {T w = -v1;  v1 = -v2;  v2 = w;} // swap and change sign
  template <class T> T square(const T x) const  { return( x*x );  }
  template <class T> T cube  (const T x) const  { return( x*x*x );  }
  template <class T> T mmin ( T x,  T y) const { return (x<y)?x:y; }
  template <class T> T mmax ( T x,  T y) const { return (x>y)?x:y; }
  template <class T> T abs(T x) const { return (x>=0)?x:-x; }
  template <class T> T setSign(T target,T sign) const { return (sign>=0)?(abs(target)):(-abs(target)); }
  template <class T> int sign(T x) const { return x>=0?1:-1; }
  template <class T> double mantissa(T x) { return(x-floor(x)); }
  //*********************************************************************
  // Heaviside function H = H(x)
  static inline double Heaviside(double x) { return x>0?1:0; }
  //*********************************************************************
  // "smooth" Heaviside function H = H(x)
  // H(x)=(1+tanh(4*x/w))/2 = 1/(1+exp(-8*x/w))
  // |smoothHeaviside(x)-H(x)|<0.02 if |x|>w/2
  //    where H(x) is the Heaviside step function
  double smoothHeaviside(double x, double w=1) { return w==0?Heaviside(x):(1/(1+exp(-8*x/w))); }

private:
  /// Cubic interpolating spline.
  ///  Adapted from the text:
  ///  Forsythe G.E., Malcolm M.A., and Moler C.B.(1977),
  ///  "Computer Methods for Mathematical Computations".
  /// @return
  ///       \li = 0 normal return
  ///       \li = 1 less than two data points; cannot interpolate
  ///       \li = 2 x[] are not in ascending order
  int spline(int n,const double x[],const double y[],double b[],double c[],double d[],
            int end1=0,int end2=0,double slope1=0,double slope2=0);
  ///  Evaluate the cubic spline function and the derivative of the cubic spline function.
  double seval (int n,double u,const double x[],const double y[],
                      const double b[],const double c[],const double d[],double &yp);

  /// Cubic interpolating spline.
  ///  Idea from the text:
  ///  "Numerical Recipes in C"
  /// @return
  ///       \li = 0 normal return
  ///       \li = 1 less than three data points; cannot interpolate
  int splineNRC(int n, const double *x, const double *y, double *y2, double yp1=1e30, double ypn=1e30);
  ///  Evaluate the cubic spline function and the derivative of the cubic spline function.
  double sevalNRC(int n, const double *xa, const double *ya, const double *y2a, double x, double &yp);
  
  int SearchXindx(double x, bool halfMesh=true) const;

private:
  /// \name Downhill simplex method
  //@{

/// Function object returns the magnetic field for given point in magnetic coordinates.
  class amoebaFunc {
    double s_; 
    int sign_;
    int ncall;
    CStconfig *mC;
  public:
    amoebaFunc(CStconfig *magConf, double s, int sign) { mC = magConf; s_=s; sign_=sign;ncall=0;}
/// The operator returns the magnetic field for given point in magnetic coordinates.
    double operator() ( const double *x ) {
      ncall++;
//////04jun2011
////mC->NJacobian(Vector3d(s_,x[0],x[1]));
////Vector3d Bvec  = mC->getBcyl();
////return sign_*Bvec.abs(); 
//////04jun2011
      return sign_*mC->B(Vector3d(s_,x[0],x[1])); // sign_*B(s,theta, phi)
    }
    int ncalls() {return ncall;}
  };
  void   amoeba(double **p,double *y,int ndim,double ftol,amoebaFunc &func);
  double amotry(double **p,double *y,double *psum,int ndim,int ihi,double fac,double *ptry,amoebaFunc &func);
  //@}

private: // these are former public functions I disabled
  /// \name These are the former public functions I've disabled.
  //@{
  /// The method makes transformation from cartesian to magnetic coordinates using guess value.
  /// @param xyz is the cartesian coordinates
  /// @param guess is the guess value
  /// @return magnetic coordinates.
  /// \note After calling this function the other quantities are available through methods:
  /// getmodB(), getJac(), getJac2(), getJ(), getJ2(), getBxyz(),
  /// getBcyl(), getGradPsi(), getGrads(), getGradreff().
  Vector3d  xyz2mag(const Vector3d &xyz,const Vector3d &guess);
  // g00 :  dV/ds = |g00(s)*Np|
  double g00   (double s) const;   
  // Following function returns 'true' if the point lies inside LCMS
  bool  xyzIsInside1(Vector3d &xyz) const;
  bool  cylIsInside1(Vector3d &cyl) const;
  //@}

/** \example 1stexample.cpp
 * This is a minimal example of how to use the CStconfig class.
 */
/** \example 2ndexample.cpp
 * This is an example of how to use the CStconfig class.
 */
/** \example remesh.cpp
 * This is an example of how to use the CStconfig::writeasciiReduced() method.
 * The program truncates Fourier spectrum to given accuracy and transforms
 * everything on new s-mesh that is equidistant on r<sub>eff</sub>.
 * The error distribution is printed.
 */
/** \example accuracytest.cpp
 * This is an example of how to use the CStconfig class.
 * The accuracy of coordinate transformation is verified and statistics is printed.
 */
/** \example cpp2for.cpp
 * This is an example of how to write interface C-function in order to call
 * C++ function from within FORTRAN programs.
 */
/** \example eeff_example.cpp
 * This is an example of how to calculate and output the effective helical ripple for \f$1/\nu\f$ transport.<br>
 * <a class="el" href="http://link.aip.org/link/doi/10.1063/1.873749">&nbsp;&nbsp;V.V.Nemov,S.V.Kasilov,W.Kernbichler,M.F.Heyn, Physics of Plasmas, vol.6,4622(1999).</a><br>
 * <a class="el" href="http://epsppd.epfl.ch/Madeira/html/pdf/P5.055.pdf">&nbsp;&nbsp;V.V.Nemov,S.V.Kasilov,W.Kernbichler,M.F.Heyn, 28th EPS, Madeira.</a>
 */
/** \example Fbs_example.cpp
 * This is an example of how to calculate and output the bootstrap current geometric factor for \f$1/\nu\f$ transport.<br>
 * <a class="el" href="http://link.aip.org/link/doi/10.1063/1.3633940">&nbsp;&nbsp;P.Helander et al, Physics of Plasmas vol.18, 092505(2011)</a><br>
 * <a class="el" href="http://dx.doi.org/10.1088/0029-5515/29/4/006">&nbsp;&nbsp;N.Nakajima et al, Nuclear Fusion, vol.29, 605(1989)</a><br>
 * <a class="el" href="http://dx.doi.org/10.1088/0741-3335/46/1/011">&nbsp;&nbsp;V.Nemov et al, Plasma Pysics and Control Fusion, vol.46, 179(2004)</a>
 */
/** \example Gv_example.cpp
 * This is an example of how to calculate and output the gradB drift velocity of trapped particles in stellarators.<br>
 *  <a class="el" href="http://link.aip.org/link/doi/10.1063/1.2131848">&nbsp;&nbsp;The gradB drift velocity of trapped particles in stellarators. Physics of Plasmas vol.12,112507(2005).</a>
 */
  // Gamma_v defined by eq.(28)
};

#ifdef _MSC_VER
#pragma warning( pop )
#endif

};   //namespace MConf

#endif //MC_CSTCONFIG_H


/* usage example

#include "CStconfig.h"

int main() {
  CStconfig mConf;
  mConf.load ("w7x-sc1.bin4");      //load configuration with Boozer coordinates
  mConf.truncate(1e-5);             //truncate spectrum
  mConf.setAccuracy(epsA);          //set accuracy of transformation
  mConf.setB0(2.5);                 //set minimum value of B on magn. axis

  Vector3d cyl(6,1*degree,0.5);     // R=6,fi=1deg,Z=0.5
  //  coordinate transformation form cylindrical to magnetic coordinates
  Vector3d boozer = mConf.cyl2mag(cyl); // boozer is (s,theta,phi)
  Vector3d B = mConf.getBcyl(cyl);         // B-vector in cylindrical coord.

  double s = boozer[0];               // norm. toroidal flux
  double reff = mConf.r(s);           // effective r
}

*/
