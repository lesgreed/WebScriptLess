#include "../include/CStconfig.h"
#include "../include/threads.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <exception>
#include <time.h>

#define CRAIG_BMN_TEST_FIELD 0 

// undefine this if you don't want debug print
//#define _TESTPRINT

#if (defined( _DEBUG ) || defined(_TESTPRINT) )
 #define  TESTPRINT(x)  {(x);}
#else 
 #define  TESTPRINT(x) {;}
#endif

#define  TESTPRINT1(x) {(x);}

namespace MConf { 

static const bool expandLCMS = false;  // do not set this option to true, it needs more time to tune
static const bool extrapolateUsingSpline = true;

static const char *Version = "21.08";

const char * CStconfig::getVersion() {return Version;}

#if defined( _DEBUG )
  static const int  nThreads = 16;
#else 
  static const int  nThreads = 16;
#endif

#define boozerComment    "Boozer-coordinate data file" 
#define tokamakComment1  "Tokamak Symmetry Flux Coordinates"   //PEST coord. see also W.D.D'haeseleer, W.N.G.Hitchon,W.I.van Rij,S.P.Hirshman,and J.L.Shohet,Flux Coordinates and Magnetic Field Structure (Springer-Verlag,New York,1991).
#define tokamakComment2  "tokamak configuration, use phi=pi/2 for summation"
#define vmecComment      "VMEC coordinates, symmetrical mode"
#define fpComment        "Configuration resulting from FP"    // bc-like-vmec-output from Function Parametrization

  #define h_bmn(i)      fCoeff[i    ]      // fc[i].bmn
  #define h_bmn_d2(i)   fCoeff[i + 1]
//Rmn, Zmn in VMEC coordinates are defined on full-mesh
//Rmn, Zmn in Boozer coordinates are defined on half-mesh
//full-mesh or half-mesh depending on the value of the fCoeffPointer
  #define h_Rmn(i)      fCoeffPointer[i + 2]
  #define h_Rmn_d2(i)   fCoeffPointer[i + 3]
  #define h_Zmn(i)      fCoeffPointer[i + 4]
  #define h_Zmn_d2(i)   fCoeffPointer[i + 5]

  #define h_Phimn(i)           fCoeff[i + 6]
  #define h_Phimn_d2(i)        fCoeff[i + 7]

  #define h_lmn(i)                fCoeff[i + 8]
  #define h_lmn_d2(i)             fCoeff[i + 9]
  #define h_gmn(i)                fCoeff[i +10]
  #define h_gmn_d2(i)             fCoeff[i +11]
  #define h_bContraTheta_mn(i)    fCoeff[i +12]
  #define h_bContraTheta_mn_d2(i) fCoeff[i +13]
  #define h_bContraFi_mn(i)       fCoeff[i +14]
  #define h_bContraFi_mn_d2(i)    fCoeff[i +15]

 // VMEC additional fourier coefficients in asymmetrical mode  if vmecAsym is true

  #define h_bmns(i)      fCoeff[i + 16]      // fc[i].bmn
  #define h_bmns_d2(i)   fCoeff[i + 17]
//Rmn, Zmn in VMEC coordinates are defined on full-mesh
//Rmn, Zmn in Boozer coordinates are defined on half-mesh
//full-mesh or half-mesh depending on the value of the fCoeffPointer
  #define h_Rmns(i)        fCoeffPointer[i + 18]
  #define h_Rmns_d2(i)     fCoeffPointer[i + 19]
  #define h_Zmnc(i)        fCoeffPointer[i + 20]
  #define h_Zmnc_d2(i)     fCoeffPointer[i + 21]

  #define h_Phimns(i)             fCoeff[i + 22]
  #define h_Phimns_d2(i)          fCoeff[i + 23]

  #define h_lmnc(i)                fCoeff[i +24]
  #define h_lmnc_d2(i)             fCoeff[i +25]
  #define h_gmns(i)                fCoeff[i +26]
  #define h_gmns_d2(i)             fCoeff[i +27]
  #define h_bContraTheta_mns(i)    fCoeff[i +28]
  #define h_bContraTheta_mns_d2(i) fCoeff[i +29]
  #define h_bContraFi_mns(i)       fCoeff[i +30]
  #define h_bContraFi_mns_d2(i)    fCoeff[i +31]


//****************************************************************************
// remove 0x0d, 0x0a and append "\n"
static char *correctNewlineChar(char *buf) 
{
  if(buf) {
    int i = (int)strlen(buf);
    if(i>0)
      if(buf[i-1]==0x0d||buf[i-1]==0x0a) buf[i-1]=0;
    if(i>1) 
      if(buf[i-2]==0x0d||buf[i-2]==0x0a) buf[i-2]=0;
    strcat(buf,"\n");
  }
  return buf;
}

//****************************************************************************
// remove 0x0d, 0x0a
static char *remove0D0A(char *buf) 
{
  if(buf) {
    int i = (int)strlen(buf);
    if(i>0)
      if(buf[i-1]==0x0d||buf[i-1]==0x0a) buf[i-1]=0;
    if(i>1) 
      if(buf[i-2]==0x0d||buf[i-2]==0x0a) buf[i-2]=0;
  }
  return buf;
}

//****************************************************************************
static void skipRestOfLine(const char *buf,FILE *fp) 
{
  if(buf) {
    int i = (int)strlen(buf);
    if(i>0) { 
      if(buf[i-1]!='\n') 
        for(int c=fgetc(fp);(c!='\n'&&c!=EOF);c=fgetc(fp)) {;}
    }
  }
}

//----------------------------------------------------------------------------

#define INTERPOLATION_ORDER 3  // must be 3  for spline interpolation

#define CORRECT_SIGNES  1      // if 1 then correct the sign of mFlux

/*
// ns == mNs, so that ms[mNs] is the last array element
#define CALCOEFF(y,ns) {splineNRC(ns,mreff.constArray(),y.constArray(),y._d2.rw(), 0 ); y._d2.rw()[ns-2]=0; y.created=true;}

//****************************************************************************
// Spline interpolation function, 
//  used for splining of Ip, It, iota, Volume, Vprime,.....(flux surface quantities)
// @return y(s)
#define SPLINEy__(s,y)  S_interp(s,y.constArray(),y._d2.constArray())
inline double CStconfig::S_interp(double s, const double *y, const double *y_d2) const
{
  int i = SearchSindx(s);
//  int ix = SearchXindx(sqrt(s));
//  if ( i!=ix ) 
//    i = ix;
  double h = mreff[i+1]-mreff[i];
  double b = (sqrt(s)-mreff[i])/h;
  double a = 1-b;
  double c = a*a*a-a;
  double d = b*b*b-b;
  return y[i] + b*(y[i+1]-y[i]) + (c*y_d2[i] + d*y_d2[i+1])*h*h/6;
}

//****************************************************************************
// Spline interpolation of Ip, It, iota, Volume, Vprime,.....(flux surface quantities)
// @return dy/ds
#define SPLINEy__D(s,y)  S_interp_D(s,y.constArray(),y._d2.constArray())
inline double CStconfig::S_interp_D(double s, const double *y, const double *y_d2) const
{
  s = (s<=ms[1])?(ms[1]+0.0001*(ms[2]-ms[1])):s;
  double e2r = 0.5/sqrt(s);
  int i = SearchSindx(s);
  double h = mreff[i+1]-mreff[i];
  double b = (sqrt(s)-mreff[i])/h;
  double a = 1-b;
//double c = a*a*a-a;
//double d = b*b*b-b;
//double y = a*y[i] + b*y[i+1] + (c*y_d2[i] + d*y_d2[i+1])*h*h/6;
  double yp = (y[i+1]-y[i])/h-( (3*a*a-1)*y_d2[i]-(3*b*b-1)*y_d2[i+1] )*h/6;
  return yp*e2r;
}
*/

/*
//****************************************************************************
// Spline interpolation of harmonics defined on half-mesh or full-mesh
// meshType is mreff.constArray() or mreffFull.constArray()
#define CmnSPLINE(Cmn,meshType) \
  int i = SearchSindx(s, meshType==mreff.constArray());    \
  double h = meshType[i+1]-meshType[i];  \
  double b = (sqrt(s)-meshType[i])/h;  \
  double a = 1-b;                    \
  double c = a*a*a-a;        \
  double d = b*b*b-b;        \
  double C=Cmn(i,m,n) + b*(Cmn(i+1,m,n)-Cmn(i,m,n)) + (c*Cmn##_d2(i,m,n) + d*Cmn##_d2(i+1,m,n))*h*h/6;
//****************************************************************************
// Spline interpolation of harmonics defined on half-mesh or full-mesh
// meshType is mreff.constArray() or mreffFull.constArray()
// @return dCmn/ds
#define CmnSPLINE_D(Cmn,meshType)           \
  s = (s<=ms[1])?(ms[1]+0.0001*(ms[2]-ms[1])):s; \
  double e2r = 0.5/sqrt(s);   \
  int i = SearchSindx(s, meshType==mreff.constArray());    \
  double h = meshType[i+1]-meshType[i];  \
  double b = (sqrt(s)-meshType[i])/h;  \
  double a = 1-b;                    \
  double Cp=(Cmn(i+1,m,n)-Cmn(i,m,n))/h - ((3*a*a-1)*Cmn##_d2(i,m,n) - (3*b*b-1)*Cmn##_d2(i+1,m,n))*h/6; \
  Cp *= e2r;

*/

//****************************************************************************
// Spline interpolation of harmonics defined on half-mesh or full-mesh
// halfMesh is true for mreff otherwise  mreffFull is used
#define CmnSPLINE(Cmn,halfMesh) \
  const double *x = halfMesh?(mreff.constArray()):(mreffFull.constArray()); \
  int i = SearchSindx(s, halfMesh);    \
  double h = x[i+1]-x[i];  \
  double b = (sqrt(s)-x[i])/h;  \
  double a = 1-b;               \
  double c = a*a*a-a;        \
  double d = b*b*b-b;        \
  double C=Cmn(i,m,n) + b*(Cmn(i+1,m,n)-Cmn(i,m,n)) + (c*Cmn##_d2(i,m,n) + d*Cmn##_d2(i+1,m,n))*h*h/6;

//****************************************************************************
// Spline interpolation of harmonics defined on half-mesh or full-mesh
// halfMesh is true for mreff otherwise  mreffFull is used
// @return dCmn/ds
#define CmnSPLINE_D(Cmn,halfMesh)           \
  s = (s<=ms[1])?(ms[1]+0.0001*(ms[2]-ms[1])):s; \
  double e2r = 0.5/sqrt(s);   \
  const double *x = halfMesh?(mreff.constArray()):(mreffFull.constArray()); \
  int i = SearchSindx(s, halfMesh);    \
  double h = x[i+1]-x[i];  \
  double b = (sqrt(s)-x[i])/h;  \
  double a = 1-b;                    \
  double Cp=(Cmn(i+1,m,n)-Cmn(i,m,n))/h - ((3*a*a-1)*Cmn##_d2(i,m,n) - (3*b*b-1)*Cmn##_d2(i+1,m,n))*h/6; \
  Cp *= e2r;

double CStconfig::Rmn_sp  (double s, int m, int n) const
{
  bool halfMesh = vmecCoordinates?false:true;
  CmnSPLINE(Rmn,halfMesh);return C;
};
double CStconfig::Zmn_sp  (double s, int m, int n) const
{
  bool halfMesh = vmecCoordinates?false:true;
  CmnSPLINE(Zmn,halfMesh);return C;
};
double CStconfig::Phimn_sp(double s, int m, int n) const
{
  bool halfMesh = true;
  CmnSPLINE(Phimn,halfMesh);return C;
};
double CStconfig::Bmn_sp  (double s, int m, int n) const
{
  CmnSPLINE(bmn,true); return fabs(mBnorm)*C;
};
double CStconfig::Bmn_sp_prime(double s, int m, int n) const
{
  CmnSPLINE_D(bmn,true); return fabs(mBnorm)*Cp;
};
double CStconfig::lmn_sp  (double s, int m, int n) const
{
  CmnSPLINE(lmn,true); return C;
};
double CStconfig::gmn_sp  (double s, int m, int n) const
{
  CmnSPLINE(gmn,true); return C;
};
double CStconfig::BcontrPolmn_sp  (double s, int m, int n) const
{
  CmnSPLINE(bContraTheta_mn,true); return mBnorm*C;   //return fabs(mBnorm)*C;
};
double CStconfig::BcontrTormn_sp  (double s, int m, int n) const
{
  CmnSPLINE(bContraFi_mn,true); return mBnorm*C;  //return fabs(mBnorm)*C;
};


// Spline interpolation of Ip, It, iota
#define SPLINEy(y) (y[halfMesh.is0] + halfMesh.b*(y[halfMesh.is0+1]-y[halfMesh.is0]) + halfMesh.c*y._d2[halfMesh.is0] + halfMesh.d*y._d2[halfMesh.is0+1])

//****************************************************************************
//****************************************************************************
#define SPLINE(Cmn,i,iCmn,meshType)  \
{ int i1 = i + offsetS; \
  iCmn  = h_##Cmn(i) + meshType->b*(h_##Cmn(i1)-h_##Cmn(i)) + meshType->c*h_##Cmn##_d2(i) +  meshType->d*h_##Cmn##_d2(i1); }

#define SPLINE2(Cmn,i,iCmn,iCm_n,meshType)  \
{ int i1 = i + offsetS; \
  iCmn  = h_##Cmn(i) + meshType->b*(h_##Cmn(i1)-h_##Cmn(i)) + meshType->c*h_##Cmn##_d2(i) +  meshType->d*h_##Cmn##_d2(i1); \
  i1 = -i + offsetS; \
  iCm_n  = h_##Cmn(-i) + meshType->b*(h_##Cmn(i1)-h_##Cmn(-i)) + meshType->c*h_##Cmn##_d2(-i) +  meshType->d*h_##Cmn##_d2(i1); }



//****************************************************************************
// Spline interpolation. iCmn=interpolated(Cmn), iCmns = dCmn/dr, where r = sqrt(s)
#define SPLINED(Cmn,i,iCmn,iCmns,meshType)  \
{ int i1 = i + offsetS; \
 iCmn  = h_##Cmn(i) + meshType->b*(h_##Cmn(i1)-h_##Cmn(i)) + meshType->c*h_##Cmn##_d2(i) +  meshType->d*h_##Cmn##_d2(i1); \
 iCmns = (h_##Cmn(i1)-h_##Cmn(i))*meshType->eh    + meshType->c1*h_##Cmn##_d2(i) + meshType->d1*h_##Cmn##_d2(i1); }

//****************************************************************************
// Linear interpolation. iCmn=interpolated(Cmn), iCmns = dCmn/ds
#define LINTERPD(Cmn,i,iCmn,iCmns,meshType)         \
{ int i1 = i + offsetS; \
      iCmns = (h_##Cmn(i1)-h_##Cmn(i))*meshType->eds01; \
      iCmn = h_##Cmn(i) + iCmns*meshType->ds0; }


//****************************************************************************
//****************************************************************************

//****************************************************************************
// Search nearest segment in array ms
// OUTPUT:
//   In 'M' & 'N' function retuns truncation levels
//   Function retuns nearest to 's' index in array 'ms'.
int CStconfig::getSsegMN(double s1, int &M, int &N) const
{
  splineCoeffs &hMesh = const_cast<CStconfig*>(this)->halfMesh;

  double s = (s1==0)?1e-8:s1;

  int is  = SearchSindx(s);         // perform a binary search of the segment in which s lies
  hMesh.e2rcurr = 0.5/sqrt(s);
  hMesh.is0  = is;
  hMesh.eds01=1./(ms[is+1]-ms[is]); // data for linear interpolations
  hMesh.ds0      = s-ms[is];
  hMesh.ds0eds01 = hMesh.ds0*hMesh.eds01;
  double h = mreff[is+1]-mreff[is];
  hMesh.b  = (sqrt(s)-mreff[is])/h;
  hMesh.eh = 1./h;
  hMesh.a  = 1-hMesh.b;
  hMesh.c  = hMesh.a*(hMesh.a*hMesh.a-1)*h*h/6;
  hMesh.d  = hMesh.b*(hMesh.b*hMesh.b-1)*h*h/6;
  hMesh.c1 = -(3*hMesh.a*hMesh.a-1)*h/6;
  hMesh.d1 =  (3*hMesh.b*hMesh.b-1)*h/6;

  if(vmecCoordinates) { // for VMEC the fullMesh is also needed
    splineCoeffs &fMesh = const_cast<CStconfig*>(this)->fullMesh;
    int is = SearchSindx(s,false);   // perform a binary search of the segment in which s lies
    fMesh.is0  = is;
    fMesh.eds01=1./(msFull[is+1]-msFull[is]); // data for linear interpolations
    fMesh.ds0      = s-msFull[is];
    fMesh.ds0eds01 = fMesh.ds0*fMesh.eds01;
    double h = mreffFull[is+1]-mreffFull[is];
    fMesh.b  = (sqrt(s)-mreffFull[is])/h;
    fMesh.eh = 1./h;
    fMesh.a  = 1-fMesh.b;
    fMesh.c  = fMesh.a*(fMesh.a*fMesh.a-1)*h*h/6;
    fMesh.d  = fMesh.b*(fMesh.b*fMesh.b-1)*h*h/6;
    fMesh.c1 = -(3*fMesh.a*fMesh.a-1)*h/6;
    fMesh.d1 =  (3*fMesh.b*fMesh.b-1)*h/6;
  }

  M=mmax(mpM[is],mpM[is+1]);
  N=mmax(mpN[is],mpN[is+1]);
  
  if(s1==0) M=0;

  if(mRoughCalculation) {
    //use with caution!
    // nstx.bc tokamak has many poloidal mods: M is 25 for eps=0.001  
    M=mmin(s<0.2?8:5, M);
    N=mmin(s<0.2?10:7,N);
  }
  return is;
}

//****************************************************************************
// Linear interpolation on s.
// on return y(s)
// used only in NJacobian2 for intrlolating iota, Ipol, Itor (defined on half-mesh)
inline double CStconfig::Linterp(const double *y) const
{
  return y[0] + (y[1]-y[0])*halfMesh.ds0eds01; // p=(u-x[0])/(x[1]-x[0]),
}

//****************************************************************************
// perform a binary search of the segment in which s lies
//   Function retuns nearest to 's' index in array 'ms'.
int CStconfig::SearchSindx(double s, bool halfMesh) const
{
  int k;
  int L = 0;      // Left bracket
  int R = mNs-1;  // Right bracket

  const double *x = halfMesh?(ms.constArray()):(msFull.constArray());

  if(s<=x[0]) return L;
  else if(s>=x[R]) L=R-1;
  else
    while(L+1 < R) {  // bisection
      k = (L + R)/2;
      if(s < x[k]) R = k;
      else L = k;
    }

  return L;
}

//****************************************************************************
// perform a binary search of the segment in which s lies
//   Function retuns nearest to x=sqrt(s) index in array 'mreff'.
int CStconfig::SearchXindx(double x, bool halfMesh) const
{
  int k;
  int L = 0;      // Left bracket
  int R = mNs-1;  // Right bracket

  const double *r = halfMesh?(mreff.constArray()):(mreffFull.constArray());

  if(x<=r[0]) return L;
  else if(x>=r[R]) L=R-1;
  else
    while(L+1 < R) {  // bisection
      k = (L + R)/2;
      if(x < r[k]) R = k;
      else L = k;
    }

  return L;
}

//****************************************************************************
// Extrapolate linealy to axis
// using one-side aproximation of derivative in the first left point [0]
//    x and y are input arrays with size n>2
// on return the value of function in the point u.
//  use for u<=x[0]
double CStconfig::extrapolate0(double u, const double *x, const double *y) const
{
  int i0=2,i1=1,i2=0;
  double c1=(y[i1]-y[i0])/(x[i1]-x[i0]);
  double c2=(y[i2]-y[i1])/(x[i2]-x[i1]);
  c2=(c2-c1)/(x[i2]-x[i0]);
  double yp = c1 + c2*(x[i2]-x[i1])+c2*(x[i2]-x[i0]); // derivative in point i2
//double yp = (y[i2]-y[i1])/(x[i2]-x[i1]); // derivative in point i2
  return  y[i2] + yp*(u-x[i2]);
}

//****************************************************************************
// Extrapolate linealy beyond the last right point of tabulated function
// using one-side aproximation of derivative in the last right point [n-1]
//    x and y are input arrays with size n>=2*h+1
// on return the value of function in the point u.
//  use this function for u>=x[n-1]
double CStconfig::extrapolateLCMS(double u, const double *x, const double *y,int n, int h) const
{
  int i2=n-1;  // last point
  int i1=i2-h; // midle point
  int i0=i1-h; // first point
  double c1=(y[i1]-y[i0])/(x[i1]-x[i0]);
  double c2=(y[i2]-y[i1])/(x[i2]-x[i1]);
  c2=(c2-c1)/(x[i2]-x[i0]);
  double yp = c1 + c2*(x[i2]-x[i1])+c2*(x[i2]-x[i0]); // derivative in point i2
//double yp = (y[i2]-y[i1])/(x[i2]-x[i1]); // derivative in point i2
  return  y[i2] + yp*(u-x[i2]);
}

//****************************************************************************
// clear all arrays
void CStconfig::clear()
{
  CStconfig::freeConf();
}

//****************************************************************************
void CStconfig::freeConf()
{
  init();

  mCOEF.clear();
  msFull.clear();
  mreffFull.clear();
  ms.clear();
  mreff.clear();
  msPol.clear();
  miota.clear();
  mpres.clear();
  mIpol.clear();
  mItor.clear();
  mpprime.clear();
  mg00.clear();
  mVolume.clear();
  mpM.clear(); 
  mpN.clear(); 
  mpavr_sx.clear();
  mpavrb_sx.clear();
  mpavrI_sx.clear();
  mreK.clear();
  mBmin.clear();
  mBmax.clear();
  mFtrap.clear();
  mBavrg.clear();
  mB2avrg.clear();
  mB_2avrg.clear();
  mR2avrg.clear();
  mR_2avrg.clear();
  mGrads2B_2avrg.clear();
  mNeoPolarztnPass.clear();
  mNeoPolarztnTrap.clear();
  mb12.clear();
  mb32.clear();

  msK.clear();
  mS11.clear();
  mS12.clear();
  mS21.clear();
  mS22.clear();

  mxtr.clear();
  mpavr2.clear();
  mpavrb2.clear();
  mpavrI.clear();
  mBorder.clear();
  mBorderMagCoord.clear();
  lcms.clear();
  epsilonEff.clear();
  Fbs.clear();
  averagedData.clear();
}

//****************************************************************************
bool CStconfig::resize(int ns,int m,int n)
{
  mOK = true;
  mmodBtol = 1e-15;
  mM = m;                // # of poloidal harmonic
  mN = n;                // # of toroidal harmonic
  MminCurrent = m;
  NminCurrent = n;
  MmaxCurrent = m;
  NmaxCurrent = n;
  mNs0 = ns;             // save the number of flux surface
  ns += m1stSindex;      // we will add harmonics near axis, see addSurfaces();
  mLCMSindex = ns-1;     // save the index of the LCMS(the last surface in a bc-file)
  ns += mNsAdd;          // we will add harmonics for s>1, see addSurfaces();
  mNs = ns;              // total # of flux surface

  ns++; // add one radial point, this point is used only for interpolation when s>>1

//After thorough testing I decided to use simple linear array mCOEF in order to achieve performance
  int numElem=vmecCoordinates?16:8;  //Bmn,Rmn,Zmn,Phimn, (lmn,gmn,bContraTheta_mn,bContraFi_mn) and second derivative
  numElem=fpCoordinates?10:numElem;  //Bmn,Rmn,Zmn,Phimn, lmn and second derivative
  numElem=vmecAsym?32:numElem;       // VMEC additional fourier coefficients in asymmetrical mode  if true
// possible addressing
//  offsetM=   ns*(2*n+1)*numElem;
//  offsetN=           ns*numElem;
//  offsetS=              numElem;

//  offsetM=   ns*(2*n+1)*numElem;
//  offsetS=      (2*n+1)*numElem;
//  offsetN=              numElem;

//mCOEF[ns][m+1][2*n+1][numElem];
//mCOEF[I][M][N][K];
//mCOEF[i][m][n][k] = mCOEF[k+ n*K       + m*N*K     + i*M*N*K  ]
//mCOEF[i][m][n][k] = mCOEF[k+ n*offsetN + m*offsetM + i*offsetS]
  offsetS=(m+1)*(2*n+1)*numElem;
  offsetM=      (2*n+1)*numElem;
  offsetN=              numElem;

  mOK &= mCOEF.init(0,ns*(2*n+1)*(m+1)*numElem-1,double(0));

  if(vmecCoordinates) {  //  full mesh for VMEC
    mOK &= msFull.init   (0,ns-1,0.0);
    mOK &= mreffFull.init(0,ns-1,0.0);
  }

  mOK &= ms.init     (0,ns-1,0.0);
  mOK &= mreff.init  (0,ns-1,0.0);
  mOK &= msPol.init  (0,ns-1,0.0);
  mOK &= miota.init  (0,ns-1,0.0);
  mOK &= mpres.init  (0,ns-1,0.0);
  mOK &= mIpol.init  (0,ns-1,0.0);
  mOK &= mItor.init  (0,ns-1,0.0);
  mOK &= mpprime.init(0,ns-1,0.0);
  mOK &= mg00.init   (0,ns-1,0.0);
  mOK &= mVolume.init(0,ns-1,0.0);

  mOK &= mpM.init(0,ns-1,mM); // # of poloidal harmonic
  mOK &= mpN.init(0,ns-1,mN); // # of toroidal harmonic

  if(mNxK%2==0) mNxK++; // must be odd!!

  mOK &= mreK.init(mNsK);
  mOK &= mxtr.init(mNxK);
  mOK &= mpavr2.init(mNxK);
  mOK &= mpavrI.init(mNxK);

  if(!mOK) freeConf();

  return mOK;
}

//****************************************************************************
// from FLTK-lib
// Returns a pointer to the last period in buf
const char *CStconfig::fname_ext(const char *buf) const
{
  const char *q = 0;
  const char *p = buf;
  for (p=buf; *p; p++) {
    if (*p == '/') q = 0;
#if ( defined(_WIN32) && !( defined(__MINGW32__) || defined(__MINGW64__) ) )
    else if (*p == '\\') q = 0;
#endif
    else if (*p == '.') q = p;
  }
  return q ? q : p;
}

//****************************************************************************
// from FLTK-lib
// returns a pointer to the filename, or null if name ends with '/'
// Returns a pointer to the character after the last slash,
// or to the start of the filename if there is none.
const char *CStconfig::fname_name(const char *name) const
{
  const char *p,*q;
  if (!name) return (0);
  q = name;
#if defined(_WIN32)
  if (q[0] && q[1]==':') q += 2; // skip leading drive letter
  for (p=q; *p; p++) if (*p == '/' || *p == '\\') q = p+1;
#else
  for (p=q; *p;) if (*p++ == '/') q = p;
#endif
  return q;
}

//****************************************************************************
// we use this table to find out the format of a file
static struct {
  CStconfig::fileFormat fmt;  const char * str;
}
pattern[] = {
  {CStconfig::W7X_BIN4,"<binary magnetic configuration file in W7X-format>float" },
  {CStconfig::W7X_BIN8,"<binary magnetic configuration file in W7X-format>double"},
  {CStconfig::W7X_FMT_TSFC, tokamakComment1         },
  {CStconfig::W7X_FMT_TSFC, tokamakComment2         },
  {CStconfig::W7X_FMT_VMEC, vmecComment             },
  {CStconfig::W7X_FMT_VMEC_FP, fpComment            },
  {CStconfig::VMEC,    "VMEC VERSION"               },
  {CStconfig::VMECNETCDF,    "CDF"                  },
  {CStconfig::W7X_FMT, boozerComment                },
  {CStconfig::W7X_FMT, "(phib-phi)*nper/twopi"      },
  {CStconfig::LHD,     "nofm, ns, ?, msym, a, R"    },  // H.Maassberg's comments in LHD boozer file
  {CStconfig::LHD1,    "D+00 "                      },  // must be LHD-format
  {CStconfig::EFIT,    "LIUQE"                      },
  {CStconfig::EFIT,    "CHEASE"                     },
  {CStconfig::EFIT,    "EFIT"                       },
  {CStconfig::EFIT,    "EFITD"                      },
  {CStconfig::EFIT,    "JET"                        },
  {CStconfig::EFIT,    "circular tokamak"           },
  {CStconfig::EFIT,    "EQDSK"                      },
  {CStconfig::EFIT,    "eqdsk"                      },
  {CStconfig::EFIT,    "micdu"                      }

//  {CStconfig::TORBEAM, "Rin,Rout"                   },
//  {CStconfig::TORBEAM, "Radial grid coordinates"    }
  };

//****************************************************************************
// Ascertain(find out) the file format
//
CStconfig::fileFormat CStconfig::getfileformat(const char * fullname,const char *nam,long fPostn) const
{
  fileFormat fmt=UNDEFINED;
  FILE * fp = fopen(fullname,"rb");
  if(fp==NULL) return fmt;

  fseek(fp,fPostn,SEEK_SET); //if(fPostn!=0)  fsetpos(fp,&fPostn);    // set position to read

  char temp[512];
  int k=50; // read 50 strings and compare with records in the pattern table
  while(--k && fmt==UNDEFINED) {
    if(fgets(temp,510,fp)==NULL) break;
    skipRestOfLine(temp,fp);
    remove0D0A(temp);
    if(strlen(temp)<3) continue;
    int sz=sizeof(pattern)/sizeof(pattern[0]);
    for(int i=0; i<sz; i++)
      if(strstr(temp,pattern[i].str)) {
        fmt=pattern[i].fmt;
        break;  // OK, we have the format
      }
  }

  fclose(fp);

  if(fmt==LHD1) {
    int n=(int)strlen(nam);
    if(strncmp(nam,"lhd", mmin(n,3))==0)
      fmt=LHD;
    else
      fmt=UNDEFINED;
  }
  if(fmt==VMEC) {
    std::istringstream str(temp+15);
    double version;
    str >> version;
    int iversion = int(1000*version); // int iversion = int(1000*atof(temp+15));
    if (iversion < 6200) fmt=UNDEFINED;
  }
  return fmt;
  //////if     (strncmp(ext,".bin4",mmin(i,5))==0) fmt = W7X_BIN4;
  //////else if(strncmp(ext,".bin8",mmin(i,5))==0) fmt = W7X_BIN8;
  //////else if(strncmp(ext,".bc",  mmin(i,3))==0) fmt = W7X_FMT;
  //////else if(strncmp(nam,"wout", mmin(n,4))==0) fmt = VMEC;}
  //////else fmt = UNDEFINED;
  //////return fmt;
}


//****************************************************************************
//
void CStconfig::init()
{
  mOK = false;  // object is no longer valid 
  mfname = "no file loaded";
  mHdr1 = "";
  mComments.clear();

  m1stSindex   = 12; // # of surfaces that will be added in the vicinity of s=0
  mNsAdd  = 12;      // # of surfaces that will be added after LCMS, must be >2
  pi       = 3.1415926535897932384626433832795;                       
  twopi    = 6.283185307179586476925286766559;
  twopiinv = 0.15915494309189533576888376337251;
  degree = pi/180;
  mu0   = 4e-7*pi; // [N/A^2 = henry/m]
  mNsK  = 30;
  mNxK  = 65;     // must be odd!!


  halfMesh.e2rcurr = 1.;
  mepsA = mepsR = 1e-5;
  mmodBtol = 1e-15;
  useVmecLambda  = true;

  torCurrentNeeded = false;
  polCurrentNeeded = false;
  vprimeNeeded   = false;
  tokamakConfig  = false;
  tokamakEQDSK   = false;
  tokamakSymmetryFluxCoordinates = false;
  vmecCoordinates = false;
  vmecAsym        = false;
  fpCoordinates  = false;
  fpCoordinates2 = false;
  fpCoordinates5 = false;
  vmecErrorCode = 0;
  mR0 = 0;
  ma  = 0;
  ftrappedReady  = false;
  BminmaxReady   = false;
  mIpol.created = false;    
  mItor.created = false;    
  mpres.created = false;    

  formatOfLoadedFile = UNDEFINED;
  jacobianIsMixedProduct = false;

  clearImportantFlags();
  zeroStatistics();
}

//****************************************************************************
void CStconfig::useMixedProductForJacobian(bool sw)
{
  jacobianIsMixedProduct = sw;
}
//****************************************************************************
// the method is called when the configuration is reloaded 
// to another coordinate system (L.H.S. <--> R.H.S) 
void CStconfig::clearImportantFlags()
{
  Averages1Ready = false;
  SMatrixReady   = false;
  mRoughCalculation = 0;
  mBnorm = 1;
  mCrossArea = 0;
  mDirTheta  = 0;
  mDirPhi    = 0;
  mSignJac_s = 0;
  mB0 = 0;       // don't remove, see getB0()
  mZ0 = 0;

  msPol.created = false;    
  miota.created = false;    
  mpprime.created = false;  
  mg00.created = false;     
  mVolume.created = false;  
  //while (!ctStack.empty()) { ctStack.pop(); }
  ct().invalidate();
}

#pragma inline_depth(0)

//****************************************************************************
//
bool CStconfig::loadTokamak(const CStconfig *mconf) 
{
  clock_t tStart = clock();  //timing
  freeConf(); // init() and clear
  tokamakSymmetryFluxCoordinates = loadEfit_p("C_TOKAMAK_EQUIVALENT", 0, 0, 0, 0, 0, 0, mconf);
  mOK = initAfterLoad(EFIT);
  if(!mOK) {
    mfname = "no file loaded";
    freeConf();
  }
//std::cerr <<"loading time="<<double(clock()-tStart)/CLOCKS_PER_SEC<<std::endl;
  return mOK;
}


//****************************************************************************
// The method loads the \ref GEQDSK "G EQDSK file" 
// and transforms it to the \ref TokamakSymmetryFluxCoordinates "tokamak-symmetry flux coordinates".
bool CStconfig::loadEfitAndScale(const char * fname, double scaleBpol, double scaleBtor,int signBpol,int signBtor,int signQ,int psiOverTwopi)
{
  const char * nam = fname_name(fname);
  freeConf();
  fileFormat fmt = getfileformat(fname,nam,0);

  if(fmt==EFIT) tokamakSymmetryFluxCoordinates = loadEfit_p(fname, scaleBpol,scaleBtor,signBpol,signBtor,signQ,psiOverTwopi);
  if(tokamakSymmetryFluxCoordinates) mOK = initAfterLoad(fmt);
  else mOK = false;
  if(!mOK) {
    mfname = "no file loaded";
    freeConf();
  }
  return mOK;
}


//****************************************************************************
//
// load of a file must be done from this method
bool CStconfig::loadfile(const char * fname, fileFormat fmt, long fPostn)
{
  clock_t tStart = clock();  //timing

  const char * ext = fname_ext (fname);
  const char * nam = fname_name(fname);
  int i = (int)strlen(ext);
  int n = (int)strlen(nam);
  
  freeConf(); // init() and clear

  if(fmt==UNDEFINED) fmt = getfileformat(fname,nam,fPostn);

  switch (fmt) {
    case LHD:
      if(fPostn!=0) return false;
      loadlhd(fname);
      break;
    case EFIT:
      tokamakSymmetryFluxCoordinates = loadEfit_p(fname,0,0,0,fPostn);
      break;
    case W7X_BIN4:
      if(fPostn!=0) return false;
      loadbin(fname, static_cast<float>(0));
      break;
    case W7X_BIN8:
      if(fPostn!=0) return false;
      loadbin(fname, static_cast<double>(0));
      break;
    case W7X_FMT:
      loadascii(fname,fmt,fPostn);
      break;
    case W7X_FMT_TSFC:
      tokamakSymmetryFluxCoordinates = loadascii(fname,fmt,fPostn);
      break;
    case W7X_FMT_VMEC:
      vmecCoordinates = loadascii(fname,fmt,fPostn); 
      break;
    case W7X_FMT_VMEC_FP:
      fpCoordinates = loadascii(fname,fmt,fPostn); 
      break;
    case VMEC:
      vmecCoordinates = loadvmec(fname,fPostn);
      break;
    case VMECNETCDF:
#ifdef NETCDF
      vmecCoordinates = loadvmecnetcdf(fname,fPostn);
      break;
#else
      return false;
#endif
    case NEMEC:
      vmecCoordinates = loadnemec(fname,fPostn);
      break;
    case TORBEAM:
      if(fPostn!=0) return false;
      loadTopfile(fname);
      break;
    case UNDEFINED:
      return false;
    default:
      return false;
  }

  mOK = initAfterLoad(fmt);
  if(!mOK) {
    mfname = "no file loaded";
    freeConf();
  }
//std::cerr <<"loading time="<<double(clock()-tStart)/CLOCKS_PER_SEC<<std::endl;

//std::cerr<<mg00[1] <<" "<<getJacobianSign()<<" "<<getJacobianSign()*mFlux<<" "<<getJacobianSign()*mFlux<<" "<<mFlux<<std::endl;
//std::cerr<<Ip(0) <<" "<<It(0.5)<<" "<<directionOfTheta()<<" "<<iota(1)<<" "<<std::endl;
//tStart = clock();  //timing
//std::cerr<<"average area="<<averageCrossSectionArea() <<" pi*a^2="<<pi*rminor()*rminor()<<std::endl;
//std::cerr <<"average area time="<<double(clock()-tStart)/CLOCKS_PER_SEC<<std::endl;
  //////{
  //////  double a = sqrt(averageCrossSectionArea()/pi);
  //////  double R0=fabs( Volume(1)/(2*square(pi*ma)) );
  //////  std::cerr <<"ma="<<ma<<" mR0="<<mR0<<std::endl;
  //////  std::cerr <<" a="<<a <<"  R0="<<R0<<std::endl;
  //////}
  return mOK;
}

//****************************************************************************
//
bool CStconfig::initAfterLoad(fileFormat fmt)
{
  if(!mOK) return mOK;
  mPeriod = twopi/mNp;
  if(fmt!=UNDEFINED) formatOfLoadedFile = fmt; // save file format if normal call 
  else fmt = formatOfLoadedFile;               // else take saved

#if 0
  //* 27.05.2010 disabled. not needed
  if(tokamakSymmetryFluxCoordinates == false)  { // is this a tokamak Symmetry Flux Coordinates ?
    tokamakSymmetryFluxCoordinates = true;       // assume it is true then all Phimn must be 0
    for(int i=0; i<mNs0; i++) {
      int is=m1stSindex+i;
      for(int m=0; m<=mM; m++)
        for(int n=-mN; n<=mN; n++)
          if( Phimn(is,m,n)!=0 ) {
            tokamakSymmetryFluxCoordinates=false;
            goto continueInit;
          }
    }
    continueInit:
    if(vmecCoordinates||mNp>1) tokamakSymmetryFluxCoordinates=false;
  }
#endif
                                      
  tokamakEQDSK  = tokamakSymmetryFluxCoordinates; //14Aug2007 =(mN==1&&mNp==1);
  tokamakConfig |= tokamakEQDSK;   //14Aug2007 tokamakConfig = (mN==0&&mNp==1);
  bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);

#if CRAIG_BMN_TEST_FIELD
  { //test for Craig : model(anti-Sigma)) B_mn for small bootstrap current
  epsilonEff.useBoozerCoord = true;
  Fbs.useBoozerCoord = true;
 #pragma message( "===> //test for Craig" )
   // ma  = 0.601; //sqrt(averageCrossSectionArea()/pi);       // VMEC definition
    //mR0 = Rmn(mLCMSindex,0,0);                      // VMEC definition
    //mR0 = Rmn(m1stSindex,0,0);                     
    for(int i=0; i<mNs0; i++) {
      int is=m1stSindex+i;
      for(int m=0; m<=mM; m++) {
        for(int n=-mN; n<=mN; n++) bmn(is,m,n) = 0;
      }
      bmn(is,0,0) =  1.0;
      bmn(is,1,0) = -0.00750;
      bmn(is,1,1) = -0.00300;
      bmn(is,2,1) = -0.05000;
      bmn(is,3,1) = -0.00300;
      miota.rw()[is] = 71./151.;  
      mIpol.rw()[is] = twopi*1.0*mR0/mu0/mNp;  
      mItor.rw()[is] = 0;  
    }
    mFlux=bmn(mLCMSindex,0,0)*pi*square(ma);        //
  }
#endif
#if HENNING_BMN_TEST_FIELD
  { //test for Henning : remove B_mn 
 #pragma message( "===> //test for Henning: remove B_mn" )
    for(int i=0; i<mNs0; i++) {
      int is=m1stSindex+i;
      for(int m=0; m<=mM; m++) 
        for(int n=-mN; n<=mN; n++) 
          bmn(is,m,n) = 0;
      bmn(is,0,0) = 1;
      bmn(is,1,0) = -0.019024;
      bmn(is,1,1) = -0.043514;
      bmn(is,0,1) = bmn(is,1,0)/2.;
      bmn(is,2,1) = bmn(is,0,1);
    }
  }
#endif

#if 0
  { //for test: remove B_mn m=11 n=2
 #pragma message( "===> //for test: remove B_mn m=11 n=2" )
    for(int i=0; i<mNs0; i++) {
      int is=m1stSindex+i;
      bmn(is,11,2)  = 0;
      Rmn(is,11,2)  = 0;
      Zmn(is,11,2)  = 0;
      Phimn(is,11,2)  = 0;
      bmn(is,22,4)  = 0;
      Rmn(is,22,4)  = 0;
      Zmn(is,22,4)  = 0;
      Phimn(is,22,4)  = 0;
    }
  }
#endif
#if 0
  { //for test: increase size
    double volumeFactor = 8; 
    double linearFactor = pow(volumeFactor,1./3); 
    double areaFactor   = pow(linearFactor,2.); 
    for(int i=0; i<mNs0; i++) {
      int is=m1stSindex+i;
      mIpol.rw()[is] *= linearFactor;  
      mItor.rw()[is] *= linearFactor;  
      mg00.rw()[is]  *= volumeFactor;
      for(int m=0; m<=mM; m++) 
        for(int n=-mN; n<=mN; n++) {
          Rmn(is,m,n)  *= linearFactor;
          Zmn(is,m,n)  *= linearFactor;
        }
    }
    ma *= linearFactor;
    mR0 *= linearFactor;
    mFlux *=areaFactor;
  }
#endif
#if 0
  { //for test: change magnetic coordinate system from Left-handed <--> R-handed
  #pragma message( "===> line 871: change magnetic coordinate system from Left-handed (L.H.S.) <--> R-handed (R.H.S)" )
    for(int i=0; i<mNs0; i++) {
      int is=m1stSindex+i;
      for(int m=0; m<=mM; m++) 
        for(int n=-mN; n<=mN; n++) 
          Zmn(is,m,n)  = -Zmn(is,m,n); 
    }
  }
#endif

  if(expandLCMS) mNs = mLCMSindex+1;
  if(extrapolateUsingSpline) mNs = mLCMSindex+1;

try {  
  ajustFourierCoeffMeshes();

  signOfJacobianAndAngles();  // calculate sign of jacobian in (s,theta,phi)-coordinates
                              //  jac = R*Fs*(Ft^Fp); where ^ is the cross product sign, F=(R,fi,Z) -cyl. coordinates
  addSurfacesNearAxis();

  if(!expandLCMS) addSurfacesAfterLCMS();
  createSplineCoeff();

//std::cerr <<"directionOfTheta "<<directionOfTheta()<<", "<<mSignJac<<std::endl;
  mB0 = 0;
  mB0 = getB0(); // save B0; mB0 must be 0 before calling

  if(ma==0)  ma = sqrt(averageCrossSectionArea()/pi);
  if(mR0==0) mR0 = Rmn_sp(1.,0,0);

  double  truncSav = truncate(1.e-7);

  calculateCurrents();
  calculateVprime();

  switch (fmt) {
    case TORBEAM:
      ma =sqrt(Volume(1.)/(mR0*2*square(pi)));
      mFlux=bmn(mLCMSindex,0,0)*pi*square(ma);  //??? TODO: must be corrected
      mR0=Rmn_sp(0.,0,0);  //Rmn(0,0,0);
      break;
    case LHD:   // correct minor & major radii
      //ma =sqrt(fabs(mFlux/bmn(mLCMSindex,0,0))/pi);  // DKES definition
      //mR0=Rmn_sp(0,0,0);                             // DKES definition
      ma = sqrt(averageCrossSectionArea()/pi);       // VMEC definition
      mR0=Rmn_sp(1.,0,0);                            // VMEC definition
      break;
    case EFIT: 
      ma = sqrt(averageCrossSectionArea()/pi);
      mR0=Rmn_sp(0.,0,0);  //Rmn(0,0,0);
      if(ma/mR0>0.01) mR0=fabs(Volume(1.))/(2*square(pi*ma)); // vprime Needed here
      break;
    default:
      break;
  }

  if(tokamakConfig) mZ0 = mag2cyl(0,0,0)[2];

  //////if(!lcms.vertex.isOK())  { // see loadbin, writebin
  //////  clock_t tStart = clock();  //timing
  //////  std::cerr <<"CStconfig: Calculating LCMS...";
  //////  lcms.create(this);  // will be created if needed
  //////  std::cerr <<"time="<<double(clock()-tStart)/CLOCKS_PER_SEC<<std::endl;
  //////}

  if(vmecCoordinates) {
    if(mHdr1.empty()) {
      mHdr1  = " m0b  n0b nsurf nper flux,[Tm^2]     a,[m]     R,[m]\n";
      mHdr2  = "       s-full        s-half           iota           curr_pol/nper    curr_tor          pprime          sqrt g(0,0)\n";
      mHdr2a = "       s             s                iota           curr_pol/nper    curr_tor          pprime          sqrt g(0,0)\n";
      mHdr21 = "                                                          [A]            [A]            dp/ds,[Pa]      (dV/ds)/nper\n";
      mHdr3  = "    m    n      r,[m]           z,[m]           pmn              Bmn,[T]           lmn            gmn        BcontraPol     BcontraTor\n";
    }
  }
  else if(fpCoordinates) {
    if(mHdr1.empty()) {
      mHdr1 = " m0b  n0b nsurf nper flux,[Tm^2]     a,[m]     R,[m]\n";
      if(mNewFormat)
        mHdr2 = "       s              iota           curr_pol/nper    curr_tor          pprime          sqrt g(0,0)\n";
      else
        mHdr2 = "       s              iota           I_pol/nper,[A]   I_tor,[A]        dp/ds,[Pa]    sqrt g(0,0)=(dV/ds)/nper,[m^3]\n";
      mHdr21 = "                                          [A]            [A]            dp/ds,[Pa]      (dV/ds)/nper\n";
//    mHdr3  = "    m    n      r,[m]           z,[m]           lmn                Bmn,[T]\n";
      mHdr3  = "    m    n      r,[m]           z,[m]           Bmn,[T]            lmn\n";
    }
  }
  else {
    if(mHdr1.empty()) {
      mHdr1 = " m0b  n0b nsurf nper flux,[Tm^2]     a,[m]     R,[m]\n";
      if(mNewFormat)
        mHdr2 = "       s              iota           curr_pol/nper    curr_tor          pprime          sqrt g(0,0)\n";
      else
        mHdr2 = "       s              iota           I_pol/nper,[A]   I_tor,[A]        dp/ds,[Pa]    sqrt g(0,0)=(dV/ds)/nper,[m^3]\n";
      mHdr21 = "                                          [A]            [A]            dp/ds,[Pa]      (dV/ds)/nper\n";
      mHdr3  = "    m    n      r,[m]           z,[m]      (phib-phi)*nper/twopi     bmn,[T]\n";
    }
  }

  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
}
  catch (std::exception& e)
  {
    std::cerr << "initAfterLoad exception: " << e.what() << std::endl;
    mOK=false;
  }
  return mOK;

  if(0) {
    clock_t tStart = clock();  //timing
    std::cerr <<"CStconfig: Calculating Bmin(s), Bmax(s) and fraction of trapped particles..."<<std::endl;
    calculateFtrappedBminmax();
    std::cerr <<"BminBmax,f.trapped time="<<double(clock()-tStart)/CLOCKS_PER_SEC<<std::endl;
  }

/*{
  FILE * fp = fopen("Rmn10.prn","w");
  int m=2,n=-3;
  m=1,n= 0;
  for(i=0; i<=mNs; i++)
     fprintf(fp,"%16.8E %16.8E %16.8E %16.8E %16.8E\n",mreff[i], Rmn(i,m,n),Rmn_b(i,m,n),Rmn_c(i,m,n),Rmn_d(i,m,n));
  fclose(fp);
}*/
  return mOK;
}

//****************************************************************************
//  set Right-handed (s,theta,phi) magnetic coordinate system 
int CStconfig::setRightCoordinateSystem()
{ 
  int jacSgn = getJacobianSign();  // get sign of jacobian in (s,theta,phi)-coordinates
  if(jacSgn<0) jacSgn = reverseCoordinateSystem();
  return jacSgn;
}

//****************************************************************************
//  set Left-handed (s,theta,phi) magnetic coordinate system 
int CStconfig::setLeftCoordinateSystem()
{ 
  int jacSgn = getJacobianSign();  // get sign of jacobian in (s,theta,phi)-coordinates
  if(jacSgn>0) jacSgn = reverseCoordinateSystem();
  return jacSgn;
}

//****************************************************************************
//  change magnetic coordinate system from Left-handed (L.H.S.) <--> R-handed (R.H.S)
int CStconfig::reverseCoordinateSystem()
{ 
  if(!mOK) return 0;
  bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);
  double truncSav = truncate();
  double mBnormSave = mBnorm;
  mBnorm = 1;

  ////ajustFourierCoeffMeshes();
  int jacSgn = getJacobianSign();  // get sign of jacobian in (s,theta,phi)-coordinates

  if(truncSav>2e-15) truncate(1e-15); 
  
  for(int i=0; i<mNs0; i++) {
    int is=m1stSindex+i;
    miota.rw()[is] = -miota.rw()[is];  
    for(int m=0; m<=mM; m++) {
      Zmn(is,m,0)   = -Zmn(is,m,0);       // sin coeff.
      if(boozerCrd) {
        Phimn(is,m,0) = -Phimn(is,m,0);   // sin coeff.
      } else if(vmecCoordinates) {
        gmn(is,m,0)  = -gmn(is,m,0);      //  jacobian must change sign
        Phimn(is,m,0) = -Phimn(is,m,0);   // sin coeff. ????2012jul1
      }
      for(int n=1; n<=mN; n++) {
        Swap (bmn(is,m,n), bmn(is,m,-n)); 
        Swap (Rmn(is,m,n), Rmn(is,m,-n)); 
        SwapN(Zmn(is,m,n),   Zmn(is,m,-n));  //*sin(m*theta-Np*n*phi), swap and change sign
        if(boozerCrd) {
          SwapN(Phimn(is,m,n), Phimn(is,m,-n)); //*sin(m*theta-Np*n*phi), swap and change sign
        } else if(fpCoordinates) {
          Swap (lmn(is,m,n),   lmn(is,m,-n) );  //*sin(m*theta-Np*n*phi); despite of sine-term must not change the sign, see formula  theta^* = theta + lambda(s,theta,fi)
        } else if(vmecCoordinates) {
          Swap (bContraTheta_mn(is,m,n),bContraTheta_mn(is,m,-n) );
          Swap (bContraFi_mn(is,m,n)   ,bContraFi_mn(is,m,-n) );
          Swap (lmn(is,m,n),            lmn(is,m,-n) );  //*sin(m*theta-Np*n*phi); despite of sine-term must not change the sign, see formula  theta^* = theta + lambda(s,theta,fi)
          SwapN(gmn(is,m,n),            gmn(is,m,-n) );  //*cos(m*theta-Np*n*phi); jacobian must change sign
          SwapN(Phimn(is,m,n),          Phimn(is,m,-n)); //*sin(m*theta-Np*n*phi), swap and change sign  ????2012jul1
        }
      }
    }
  }                                   
  clearImportantFlags();
  mOK = initAfterLoad();
  jacSgn = getJacobianSign();  // get sign of jacobian in (s,theta,phi)-coordinates
  TESTPRINT (std::cerr <<"CStconfig: Magnetic coordinate system have been changed to "<<(jacSgn>0?"R.H.S.":"L.H.S.")<<std::endl );
  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
  mBnorm = mBnormSave;
  return jacSgn;
}

//****************************************************************************
//  change magnetic coordinate system from Left-handed (L.H.S.) <--> R-handed (R.H.S) by fliping z-coordinate
int CStconfig::flipZcoordinate()
{ 
  if(!mOK) return 0;
  double truncSav = truncate();
  double mBnormSave = mBnorm;
  mBnorm = 1;

  ////ajustFourierCoeffMeshes();
  int jacSgn = getJacobianSign();  // get sign of jacobian in (s,theta,phi)-coordinates

  if(truncSav>2e-15) truncate(1e-15); 
  
  for(int i=0; i<mNs0; i++) {
    int is=m1stSindex+i;
    for(int m=0; m<=mM; m++) 
      for(int n=-mN; n<=mN; n++) {
        Zmn(is,m,n)  = -Zmn(is,m,n);
        if(vmecCoordinates) {
          gmn(is,m,n)  = -gmn(is,m,n);      //  jacobian must change sign
        }
      }
  }

  clearImportantFlags();
  Fbs.invalidate();
  lcms.clear();
  mOK = initAfterLoad();
  jacSgn = getJacobianSign();  // get sign of jacobian in (s,theta,phi)-coordinates
  TESTPRINT (std::cerr <<"CStconfig: Magnetic coordinate system have been changed to "<<(jacSgn>0?"R.H.S.":"L.H.S.")<<std::endl );
  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
  mBnorm = mBnormSave;
  return jacSgn;

}
//****************************************************************************
//
//
bool CStconfig::write(const char * fname) const
{
  if(!isOK()) return false;
  const_cast<CStconfig*>(this)->calculateCurrents();  
  const_cast<CStconfig*>(this)->calculateVprime();  
  const char * ext = fname_ext(fname);
  int i = (int)strlen(ext);
  if(i==0)  return writeascii(fname);
  if(strncmp(ext,".bin4",mmin(i,5))==0)      return writebin (fname,static_cast<float>(0));
  else if(strncmp(ext,".bin8",mmin(i,5))==0) return writebin (fname,static_cast<double>(0));
  else if(strncmp(ext,".bc",mmin(i,3))==0)   return const_cast<CStconfig*>(this)->writeascii(fname);
  else  return writeascii(fname);
}

//****************************************************************************
// Read bc-file  into object
//
//   here fscanf is used because it works faster than C++ streams >>
bool CStconfig::loadascii(const char * fname, fileFormat fmt, long fPostn)
{
  FILE * fp = fopen(fname,"rb");
  if(fp==NULL) return(mOK=false);
  mfname = fname;

  fseek(fp,fPostn,SEEK_SET); //if(fPostn!=0) fsetpos(fp,&fPostn);           // set position to read

  if(fmt==W7X_FMT_VMEC ) vmecCoordinates=true;  // need these flags to skip... 
  if(fmt==W7X_FMT_VMEC_FP ) fpCoordinates=true; //   ... or read some data from file

  const int bufSz=509;
  char buf[bufSz+3];
  
  mNewFormat = false;
  if(fPostn) {  //skip rest of line 
    for(int c=fgetc(fp);(c!='\n'&&c!=EOF);c=fgetc(fp)) {;} 
  }
  while(1) {
    if(fgets(buf,bufSz,fp)==NULL) return(mOK=false);
    skipRestOfLine(buf,fp);
    int len = (int)strlen(buf);
    if(strncmp(buf,"CC",mmin(len,2))!=0&&strncmp(buf,"cc",mmin(len,2))!=0) break;
    mNewFormat = true;  // "CC" is a signal that file is written in a new format
    mComments.push_back(correctNewlineChar(buf));
  }
  mHdr1 = correctNewlineChar(buf);

  int ns,m,n;
  fscanf(fp,"%d %d %d %d",&m,&n,&ns,&mNp);
  fscanf(fp,"%lf %lf %lf",&mFlux,&ma,&mR0);

//MagAxis*********************************
  int col_=0;
  {  // does the file contain the magnetic axis? 
    fpos_t pos0;
    fgetpos(fp,&pos0); // read position in file
    fgets(buf,bufSz,fp); // ignore rest of line
    fgets(buf,bufSz,fp); // read header "s          iota      curr_pol    curr_tor      pprime  sqrt g(0,0)"
    if(!mNewFormat) { //  file written in a new format can be without lines with "CC"
      fpos_t pos;
      fgetpos(fp,&pos);  // read position in file
      fgets(buf,bufSz,fp); // try to read header "       [A]          [A]   dp/ds,[Pa] (dV/ds)/nper"
      fsetpos(fp,&pos);  // restore position
      if(strstr(buf,"dp/ds,[Pa]")  !=NULL ||
         strstr(buf,"(dV/ds)/nper")!=NULL ) 
        mNewFormat = true;
    }
    if(mNewFormat) {
      fgets(buf,bufSz,fp); // read header "       [A]          [A]   dp/ds,[Pa] (dV/ds)/nper"
    }
    {  // now much numbers we have ?
      fpos_t pos;
      fgetpos(fp,&pos);  // read position in file
      fgets(buf,bufSz,fp); // try to read header "       [A]          [A]   dp/ds,[Pa] (dV/ds)/nper"
      fsetpos(fp,&pos);  // restore position
      for(char *tok=strtok(buf, "\n\t "); tok!=NULL; tok=strtok(NULL, "\n\t ") ) col_++;
    }
    double sfull=1,sHalf=1;
    if(vmecCoordinates) {
      fscanf(fp,"%lf",&sfull);
    }
    fscanf(fp,"%lf",&sHalf);
    if(sHalf==0) m1stSindex = 0;  // file has magnetic axis if true
    fsetpos(fp,&pos0);  // restore position
  }
//MagAxis*********************************

  if(!resize(ns,m,n)) return(mOK=false);
  int nHarm = (mM+1)*(2*mN+1);
  if(mNewFormat) nHarm -= mN;

  for(int i=0; i<mNs0; i++) {
    int is = m1stSindex+i;
    fgets(buf,bufSz,fp); // ignore rest of line
    fgets(buf,bufSz,fp); // read header "s          iota      curr_pol    curr_tor      pprime  sqrt g(0,0)"
    mHdr2 = correctNewlineChar(buf);
    if(mNewFormat) {
      fgets(buf,bufSz,fp); // read header "       [A]          [A]   dp/ds,[Pa] (dV/ds)/nper"
      mHdr21 = correctNewlineChar(buf);
    }
    double s,iot,ipol,itor,pp,g00;
    if(vmecCoordinates) {
      fscanf(fp,"%lf",&s);
      msFull.rw()[is] = s;
    }
    fscanf(fp,"%lf %lf",&s,&iot);
    ms.rw()[is]    = s; 
    miota.rw()[is] = iot;
    if(fpCoordinates) {
      fpCoordinates2 = true;
      if(col_==5) {
        fpCoordinates2 = false;
        fpCoordinates5 = true;
        fscanf(fp,"%lf %lf %lf",&ipol,&itor,&pp);
        mIpol.rw()[is] = ipol;  
        mItor.rw()[is] = itor;
        mpprime.rw()[is] = pp;
      }
    }
    else { 
      fscanf(fp,"%lf %lf %lf %lf",&ipol,&itor,&pp,&g00);
      mIpol.rw()[is] = ipol;  
      mItor.rw()[is] = itor;
      mpprime.rw()[is] = pp;
      mg00.rw()[is]    = g00;
    }

    fgets(buf,bufSz,fp); // ignore rest of line
    fgets(buf,bufSz,fp); // read header "m    n        r/[m]           z/[m] (phib-phi)*nper/twopi     bmn/[T]"
    mHdr3 = correctNewlineChar(buf);
    for(int ih=0; ih<nHarm; ih++) {
      int m,n;
      double rmnc,zmns,pmn,bmn1,lmns,gmnc,bth,bfi; 
      fscanf(fp,"%d %d",&m,&n);
      if(vmecCoordinates) {
        pmn = 0;
        fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf",&rmnc,&zmns,&pmn,&bmn1,&lmns,&gmnc,&bth,&bfi);  //see pmn, 2012jul1
      //fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&rmnc,&zmns,&bmn1,&lmns,&gmnc,&bth,&bfi);
        gmn(is, m, n) = gmnc;            // gmnc   on half-mesh
        lmn(is, m, n) = lmns;            // lmns   on half-mesh
        bContraTheta_mn(is, m, n) = bth;   // contravariant poloidal component of B on half-mesh
        bContraFi_mn   (is, m, n) = bfi;   // contravariant toroidal component of B on half-mesh
        Phimn(is, m, n) = pmn;
      }
      else if(fpCoordinates) {
        fscanf(fp,"%lf %lf %lf %lf",&rmnc,&zmns,&bmn1,&lmns);
        lmn(is, m, n) = lmns;
      }
      else { // boozer coordinates
        fscanf(fp,"%lf %lf %lf %lf",&rmnc,&zmns,&pmn,&bmn1);
        Phimn(is, m, n) = pmn;
      }
      Rmn(is, m, n) = rmnc;   
      Zmn(is, m, n) = zmns;  
      bmn(is, m, n) = bmn1;            // |B|   on half-mesh
//#if VAR_NUM_HARM     // if variable number of harmonics
      fpos_t pos;
      fgetpos(fp,&pos); // read position in file
      fgets(buf,30,fp); // ignore rest of line
      fgets(buf,15,fp); // read line
      if(feof(fp)) goto ret;
      //if(strstr(buf,"EOF")!=0) goto ret;
      fsetpos(fp,&pos); // restore position
      if(strchr(buf,'s')!=NULL)
        break; // break if we read header: "s  iota .."
//#endif
    }
  }
ret:
  if(fpCoordinates) {
    if(!fpCoordinates5){
      polCurrentNeeded = true;
      torCurrentNeeded = true;
    }
    vprimeNeeded     = true;
    useVmecLambda    = true;
  }
  fclose(fp);
  return(mOK=true);
}

//****************************************************************************
void CStconfig::writeIdentificationComment(FILE * fp) const
{
  fprintf(fp,"CC this file has been saved by 'MConf::CStconfig' from the original file '%s'\n",fname_name(mfname.c_str()));
  if(vmecCoordinates)
    fprintf(fp,"CC %s\n",vmecComment);
  else if(fpCoordinates)
    fprintf(fp,"CC %s\n",fpComment); 
  else if(tokamakSymmetryFluxCoordinates) {
    fprintf(fp,"CC %s\n",tokamakComment1);
    fprintf(fp,"CC %s\n",tokamakComment2);
  }
  else {
    fprintf(fp,"CC %s\n",boozerComment); 
  }
}

//****************************************************************************
// Write bc-file
//
bool CStconfig::writeVmec2Boozer(const char * fname, int method, int M, int N)
{
  if(!mOK) return false;
  if(!vmecCoordinates) return false;

  FILE * fp = fopen(fname,"w");
  if(fp==NULL) return(false);

  volatile ctStack::saveState state(ct);
  
  double  truncSav = truncate();
  truncate(1e-15);

  clock_t tStart = clock();  //timing
  TESTPRINT1( std::cerr <<"CStconfig::writeVmec2Boozer...." ); 

  if(M==0) M = 6*mM;
  if(N==0) N = 3*mN;
  if(N==0) N = 1;
  calcBoozFourierMods(mNs0,M,N,method);

  TESTPRINT1( std::cerr <<"time="<<double(clock()-tStart)/CLOCKS_PER_SEC<<std::endl );

  fprintf(fp, "CC Boozer-coordinate data file\n");
  fprintf(fp, "CC original file '%s' saved by 'MConf::CStconfig'\n",fname_name(mfname.c_str()));
//  fprintf(fp, "CC saved by 'mcviewer', truncation level Bmn/B00=%g\n",mmodBtol);
  std::string mHdr1 = " m0b  n0b nsurf nper flux/[Tm^2]     a/[m]     R/[m]\n";
  std::string mHdr2 = "       s              iota           curr_pol/nper    curr_tor          pprime          sqrt g(0,0)\n";
  std::string mHdr21 = "                                         [A]            [A]            dp/ds,[Pa]      (dV/ds)/nper\n";
  std::string mHdr3  = "    m    n      r/[m]           z/[m]      (phib-phi)*nper/twopi     bmn/[T]\n";
  
  int Mmax = bz.Mmax;
  int Nmax = bz.Nmax;
  fputs(mHdr1.c_str(),fp);
  fprintf(fp,"% 4d %4d %4d %4d ",Mmax,Nmax,mNs0,mNp);
  fprintf(fp,"%13.6E %8.5f %8.5f\n",mFlux,ma,mR0);
  for(int i=0; i<mNs0; i++) {
    int is = m1stSindex+i;
    fputs(mHdr2.c_str(), fp);
    fputs(mHdr21.c_str(),fp);
    fprintf(fp,"%17.10E %16.8E %16.8E %16.8E %16.8E %16.8E\n",ms[is],miota[is],mIpol[is],mItor[is],mpprime[is],mg00[is]);
    fputs(mHdr3.c_str(),fp);
    for(int m=0; m<=Mmax; m++) {
      for(int n=(m==0)?0:-Nmax; n<= Nmax; n++) {
        fprintf(fp,"% 5d %4d ", m, n);
        fprintf(fp,"%16.8E %16.8E %16.8E %16.8E\n",bz.Rmn[i](m,n),bz.Zmn[i](m,n),bz.Pmn[i](m,n),bz.Bmn[i](m,n));
      }
    }
  }
  bz.clear();
  fclose(fp);

  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
  return(true);
}

//****************************************************************************
void CStconfig::calcBoozFourierMods(int Ns, int M, int N, int method)
{
  bz.resize(Ns,M,N,method);
  Threads2<CStconfig,CStconfig> threads(this, &CStconfig::calcBoozFourierModsExe,this,"calcBoozFourierMods");  
  int low = 0, upper = Ns-1;
  threads.run(low, upper, nThreads); 
}

void CStconfig::calcBoozFourierModsExe(CStconfig * mc, int low, int upper, bool setup)
{
  if(setup) {  // setup part; it allocates arrays and initialize them
    int size = upper-low+1;
    bool integration_over_Boozer_angles = bz.getMethod();
    if(integration_over_Boozer_angles) {  // if integration_over_Boozer_angles 
      int Mmax = bz.Mmax;
      int Nmax = bz.Nmax;
      int Npol = bz.Npol;
      int Ntor = bz.Ntor;
      double dtheta = twopi/Npol;
      double dphi   = twopi/mNp/Ntor;
      for(int m=0; m<=Mmax; m++) {
        bz.ctm.rw()[m] = CComplex (cos(m*dtheta),   sin(m*dtheta));    // exp(i*m*dtheta), where i is the sqrt(-1)
      }
      for(int n=-Nmax; n<= Nmax; n++) {
        bz.cpn.rw()[n] = CComplex (cos(n*mNp*dphi),-sin(n*mNp*dphi));  // exp(-i*n*mNp*dphi), where i is the sqrt(-1)
      }
    }
  }
  else
  {  // exe part; it will be called from threads to fill arrays
    int Mmax = mc->bz.Mmax;
    int Nmax = mc->bz.Nmax;
    int Npol = mc->bz.Npol;
    int Ntor = mc->bz.Ntor;

    bool integration_over_Boozer_angles = mc->bz.getMethod();

    CArray2d<double> R(Npol,Ntor), B(Npol,Ntor), Z(Npol,Ntor), P(Npol,Ntor);
    CArray2d<Vector3d> Booz;
    if(!integration_over_Boozer_angles) Booz.resize(Npol,Ntor);

    for(int i=low; i<=upper; i++) {
      int is = m1stSindex+i;
      if(integration_over_Boozer_angles) {  // if Fourier integration_over_Boozer_angles 
        mc->createSurfaceDataOnBoozerGrid(mc->ms[is], R, Z, P, B); 
      }
      else {
        mc->createSurfaceDataOnVmecGrid(mc->ms[is],Booz, R, Z, P, B); 
      }
      bz.Rmn[i].resize(0,Mmax, -Nmax,Nmax);
      bz.Zmn[i].resize(0,Mmax, -Nmax,Nmax);
      bz.Pmn[i].resize(0,Mmax, -Nmax,Nmax);
      bz.Bmn[i].resize(0,Mmax, -Nmax,Nmax);
      double *pRmn = bz.Rmn[i].ptr();
      double *pZmn = bz.Zmn[i].ptr();
      double *pPmn = bz.Pmn[i].ptr();
      double *pBmn = bz.Bmn[i].ptr();
      for(int m=0; m<=Mmax; m++) {
        for(int n=(m==0)?0:-Nmax; n<= Nmax; n++) {
          double rmn,zmn,pmn,bmn;        
          if(integration_over_Boozer_angles) {  // if Fourier integration_over_Boozer_angles 
            FourierUsingBoozerGrid(bz.ctm[m], bz.cpn[n], R, Z, P, B, rmn,zmn,pmn,bmn); // ??? mc->FourierUsingBoozerGrid
          }
          else {
            FourierUsingVmecGrid(m,n,Booz, R, Z, P, B, rmn,zmn,pmn,bmn);
          }
          if(m==0&&n==0) {
            rmn *= .5;
            zmn *= .5;
            pmn *= .5;
            bmn *= .5;
          }
          int idx = bz.Rmn[i].index(m,n);  // get linear index
          pRmn[idx] = rmn;    //bz.Rmn[i](m,n) = rmn;
          pZmn[idx] = zmn;    //bz.Zmn[i](m,n) = zmn;
          pPmn[idx] = pmn;    //bz.Pmn[i](m,n) = pmn;
          pBmn[idx] = bmn;    //bz.Bmn[i](m,n) = bmn;
        }
      }
    }
  }
}

//****************************************************************************
// Bmn by integrating over Boozer angles 
double CStconfig::Boozer2Bmn(double s,int m,int n)
{
  bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);
  if(boozerCrd) return Bmn_sp(s, m, n);
  if(m==0&&n<0) return 0;
  if(!vmecCoordinates) return 0;

  volatile CStconfig::ctStack::saveState state(ct);
  double  truncSav = truncate();
  truncate(1e-15);

  int M = (m==0)?1:abs(m);
  int N = (n==0)?1:abs(n);

  if(n==0) N = 1;

  int Npol = M*32;  //4
  int Ntor = N*32;  //4
  CArray2d<double> B(Npol,Ntor);
  double dtheta = twopi/Npol;
  double dphi   = twopi/mNp/Ntor;

  { //  createSurfaceDataOnBoozerGrid(s, R, Z, P, B); // Create data on the Boozer-angles grid
    double *pB = B.ptr();
    for(int k=0, idx=0; k<Npol; k++) {
      for(int j=0; j<Ntor; j++) {
        Vector3d booz(s,k*dtheta,j*dphi);
        double b;
        {   // V2Btrans(booz,r,z,p,b);
          Vector3d vmec = Boozer2Vmec(booz);    
          b = CStconfig::B(vmec);
        }
        //int idx = R.index(k,j);  // get linear index  //int idx = k*Ntor+j;
        pB[idx] = b;    // B(k,j)=z;
        ++idx;
      }
    }
  }
  double bmn=0;        
  { // FourierUsingBoozerGrid
    CComplex ctm(cos(m*dtheta),   sin(m*dtheta));    // exp( i*m*dtheta)
    CComplex cpn(cos(n*mNp*dphi),-sin(n*mNp*dphi));  // exp(-i*n*mNp*dphi)
    CComplex expa(1,0);
    const double *pB = B.cptr();
    for(int k=0, idx=0; k<Npol; k++) {
      for(int j=0; j<Ntor; j++) {
        //double cs = expa.re;  // expa = cexp( i*(m*dtheta*k - n*mNp*dphk*j) )
        //double sn = expa.im;  
        //int idx = R.index(k,j); // get linear index; //int idx = k*Ntor+j;
        bmn += pB[idx]*expa.re;    // B(k,j)*cs;
        ++idx;
        expa *= cpn;
      }
      expa *= ctm;
    }
    bmn *= 2./(Npol*Ntor);
  }
  if(m==0&&n==0) {
    bmn *= .5;
  }
  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
  return bmn;   
}

//****************************************************************************
// Bmn by integrating over VMEC angles 
double CStconfig::BoozerBmn(double s,int m,int n)
{
  bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);
  if(boozerCrd) return Bmn_sp(s, m, n);
  if(m==0&&n<0) return 0;
  if(!vmecCoordinates) return 0;

  volatile CStconfig::ctStack::saveState state(ct);
  double  truncSav = truncate();
  truncate(1e-15);

  int M = (m==0)?1:abs(m);
  int N = (n==0)?1:abs(n);

  if(n==0) N = 1;

  int Npol = M*32;  //4
  int Ntor = N*32;  //4
  CArray2d<double> B(Npol,Ntor);
  CArray2d<Vector3d> booz(Npol,Ntor); 
  { // Create data on the VMEC-angles grid, see also createSurfaceDataOnVmecGrid(s,Booz, R, Z, P, B);
    double dtheta = twopi/Npol;
    double dphi   = twopi/mNp/Ntor;
    double *pB = B.ptr();
    Vector3d *pBz = booz.ptr();
    for(int k=0,idx=0; k<Npol; k++) {
      for(int j=0; j<Ntor; j++) {
        Vector3d vmec(s,k*dtheta,j*dphi);
        //int idx = R.index(k,j);  // get linear index
        //int idx = k*Ntor+j;
        double b,jv_jb;
        {   //   V2Btrans(vmec,pBz[idx],r,z,p,b,jv_jb);   // pBz[idx] is the same as booz(k,j)
          calcVmecLgB(vmec,0,true);
          double iota =  SPLINEy(miota);
          Vector3d &boozer = pBz[idx];
          boozer = vmec;
          boozer[1] += ct().lambda + iota*ct().p;
          boozer[2] += ct().p;
          jv_jb = (1+ct().lambdat)*(1+ct().pp) + ct().pt*(iota-ct().lambdap); // d(theta_booz,phi_booz)/d(theta_vmec,phi_vmec)
          b     = CStconfig::B(vmec);
        }
        pB[idx] = b*jv_jb;    // B(k,j)=z*jv_jb;
        ++idx;
      }
    }
  }

  double bmn=0;        
// Fourier transformation by integrating over VMEC angles
  { // see also FourierUsingVmecGrid(m,n,booz, R, Z, P, B, rmn,zmn,pmn,bmn);
    const double *pB = B.cptr();
    const Vector3d *pBz = booz.cptr();
    for(int k=0,idx=0; k<Npol; k++) {
      for(int j=0; j<Ntor; j++) {
        //int idx = R.index(k,j); // get linear index
        //int idx = k*Ntor+j;
        double a  = m*pBz[idx][1] - n*mNp*pBz[idx][2]; // double a  = m*bz[1] - n*mNp*bz[2];
        double cs = ::cos(a);
        bmn += pB[idx]*cs;    //Bmn += B(k,j)*cs;
        ++idx;
      }
    }
    bmn *= 2./(Npol*Ntor);
  }

  if(m==0&&n==0) {
    bmn *= .5;
  }
    if(truncSav != truncate()) truncate(truncSav); // restore truncation level
  return bmn;

}

//****************************************************************************
// Create data on the Boozer-angles grid
void CStconfig::createSurfaceDataOnBoozerGrid(double s, 
                                  CArray2d<double> &R, CArray2d<double> &Z, 
                                  CArray2d<double> &P, CArray2d<double> &B) 
{
  double mBnormSave = mBnorm;
  mBnorm = 1;
  int Npol = R.size1();
  int Ntor = R.size2();
  double dtheta = twopi/Npol;
  double dphi   = twopi/mNp/Ntor;
  double *pR = R.ptr();
  double *pZ = Z.ptr();
  double *pB = B.ptr();
  double *pP = P.ptr();
  for(int k=0, idx=0; k<Npol; k++) {
    for(int j=0; j<Ntor; j++) {
      //Vector3d vmec(s,-pi+i*dtheta,-pi/mNp + j*dphi);
      Vector3d booz(s,k*dtheta,j*dphi);
      double r,z,p,b;
      V2Btrans(booz,r,z,p,b);
      //int idx = R.index(k,j);  // get linear index
      //int idx = k*Ntor+j;
      pR[idx] = r;    // R(k,j)=r;
      pZ[idx] = z;    // Z(k,j)=z;
      pP[idx] = p;    // P(k,j)=r;
      pB[idx] = b;    // B(k,j)=z;
      ++idx;
    }
  }
  mBnorm = mBnormSave;
}

//****************************************************************************
//Fourier transformation by integrating over Boozer angles
// @param[in] ctm == exp( i*m*dtheta) , where i is the sqrt(-1)
// @param[in] cpn == exp(-i*n*mNp*dphi)
// @param[in] R,Z,P,B 
// @param[out] Rmn,Zmn,Pmn,Bmn 
void CStconfig::FourierUsingBoozerGrid(const CComplex &ctm, const CComplex &cpn, 
                        const CArray2d<double> &R, const CArray2d<double> &Z, 
                        const CArray2d<double> &P, const CArray2d<double> &B,                      
                        double &Rmn,double &Zmn,double &Pmn,double &Bmn) 
{
  int Npol = R.size1();
  int Ntor = R.size2();
  ////double dtheta = twopi/Npol;
  ////double dphi   = twopi/mNp/Ntor;
  ////complex<double> ctm(cos(m*dtheta),   sin(m*dtheta));    // exp( i*m*dtheta)
  ////complex<double> cpn(cos(n*mNp*dphi),-sin(n*mNp*dphi));  // exp(-i*n*mNp*dphi)
  CComplex expa(1,0);
  Rmn = Zmn = Pmn = Bmn = 0;
  const double *pR = R.cptr();
  const double *pZ = Z.cptr();
  const double *pB = B.cptr();
  const double *pP = P.cptr();
  for(int k=0, idx=0; k<Npol; k++) {
    for(int j=0; j<Ntor; j++) {
      //double cs = expa.re;  // expa = cexp( i*(m*dtheta*k - n*mNp*dphk*j) )
      //double sn = expa.im;  
      //int idx = R.index(k,j); // get linear index
      //int idx = k*Ntor+j;
#pragma message( "===> Don't change the class CArray2d, otherwise indexing will be destroyed.")
      Rmn += pR[idx]*expa.re;    // R(k,j)*cs;
      Bmn += pB[idx]*expa.re;    // B(k,j)*cs;
      Zmn += pZ[idx]*expa.im;    // Z(k,j)*cs;
      Pmn += pP[idx]*expa.im;    // P(k,j)*cs;
      ++idx;
      expa *= cpn;
    }
    expa *= ctm;
  }
  double c = 2./(Npol*Ntor);
  Rmn *= c;
  Zmn *= c;
  Bmn *= c;
  Pmn *= c;
  Pmn *= mNp/twopi;  // transform to Geiger's format
}

//****************************************************************************
// Create data on the VMEC-angles grid
void CStconfig::createSurfaceDataOnVmecGrid(double s, CArray2d<Vector3d> &booz,
                                  CArray2d<double> &R, CArray2d<double> &Z, 
                                  CArray2d<double> &P, CArray2d<double> &B) 
{
  double mBnormSave = mBnorm;
  mBnorm = 1;
  int Npol = R.size1();
  int Ntor = R.size2();
  double dtheta = twopi/Npol;
  double dphi   = twopi/mNp/Ntor;
  double *pR = R.ptr();
  double *pZ = Z.ptr();
  double *pB = B.ptr();
  double *pP = P.ptr();
  Vector3d *pBz = booz.ptr();
  for(int k=0,idx=0; k<Npol; k++) {
    for(int j=0; j<Ntor; j++) {
      //Vector3d vmec(s,-pi+i*dtheta,-pi/mNp + j*dphi);
      Vector3d vmec(s,k*dtheta,j*dphi);
      double r,z,p,b,jv_jb;
      //int idx = R.index(k,j);  // get linear index
      //int idx = k*Ntor+j;
      V2Btrans(vmec,pBz[idx],r,z,p,b,jv_jb);   // pBz[idx] is the same as booz(k,j)
      pR[idx] = r*jv_jb;    // R(k,j)=r*jv_jb;
      pZ[idx] = z*jv_jb;    // Z(k,j)=z*jv_jb;
      pP[idx] = p*jv_jb;    // P(k,j)=r*jv_jb;
      pB[idx] = b*jv_jb;    // B(k,j)=z*jv_jb;
      ++idx;
    }
  }
  mBnorm = mBnormSave;
}
//****************************************************************************
//Fourier transformation by integrating over VMEC angles
// @param[in] booz, R,Z,P,B 
// @param[out] Rmn,Zmn,Pmn,Bmn 
void CStconfig::FourierUsingVmecGrid(int m, int n, const CArray2d<Vector3d> &booz, 
                        const CArray2d<double> &R, const CArray2d<double> &Z, 
                        const CArray2d<double> &P, const CArray2d<double> &B,                      
                        double &Rmn,double &Zmn,double &Pmn,double &Bmn) 
{
  int Npol = R.size1();
  int Ntor = R.size2();
  //if(n==0) Ntor/=2;
  Rmn = Bmn = Zmn = Pmn = 0; 
  const double *pR = R.cptr();
  const double *pZ = Z.cptr();
  const double *pB = B.cptr();
  const double *pP = P.cptr();
  const Vector3d *pBz = booz.cptr();
  for(int k=0,idx=0; k<Npol; k++) {
    for(int j=0; j<Ntor; j++) {
      //int idx = R.index(k,j); // get linear index
      //int idx = k*Ntor+j;
//#pragma message( "===> Don't change the class CArray2d, otherwise indexing will be destroyed.")
                                                      //Vector3d bz = pBz[idx]; //booz(k,j);
      double a  = m*pBz[idx][1] - n*mNp*pBz[idx][2]; // double a  = m*bz[1] - n*mNp*bz[2];
      double cs = ::cos(a);
      double sn = ::sin(a);
      Rmn += pR[idx]*cs;    //Rmn += R(k,j)*cs;
      Bmn += pB[idx]*cs;    //Bmn += B(k,j)*cs;
      Zmn += pZ[idx]*sn;    //Zmn += Z(k,j)*sn;
      Pmn += pP[idx]*sn;    //Pmn += P(k,j)*sn;
      ++idx;
    }
  }
  double c = 2./(Npol*Ntor);
  Rmn *= c;
  Zmn *= c;
  Bmn *= c;
  Pmn *= c;
  Pmn *= mNp/twopi;  // transform to Geiger's format
}

//****************************************************************************
// Write bc-file
//
bool CStconfig::writeascii(const char * fname) const
{
  FILE * fp = fopen(fname,"w");
  if(fp==NULL) return(false);

  if(mNewFormat) {
    for(unsigned int i=0; i<mComments.size(); i++)
      fputs(mComments[i].c_str(),fp);
    writeIdentificationComment(fp);
  }

  fprintf(fp, "CC saved by 'mcviewer', truncation level Bmn/B00=%g\n",mmodBtol);

  int Mmax = mmin(1,mM);
  int Nmax = mmin(1,mN);
  int i;
  for(i=0; i<mNs0; i++) {
    int is = m1stSindex+i;
    Mmax = mmax(Mmax,mpM[is]);
    Nmax = mmax(Nmax,mpN[is]);
  }

  fputs(mHdr1.c_str(),fp);
  fprintf(fp,"% 4d %4d %4d %4d ",Mmax,Nmax,mNs0,mNp);
  fprintf(fp,"%13.6E %8.5f %8.5f\n",mFlux,ma,mR0);

  for(i=0; i<mNs0; i++) {
    int is = m1stSindex+i;
    fputs(mHdr2.c_str(), fp);
    if(mNewFormat) fputs(mHdr21.c_str(),fp);
    if(vmecCoordinates) {
      fprintf(fp,"%17.10E ",msFull[is]);
    }
    if(fpCoordinates2) {
      fprintf(fp,"%17.10E %16.8E\n",ms[is],miota[is]);
    }
    else if(fpCoordinates5) {
      fprintf(fp,"%17.10E %16.8E %16.8E %16.8E %16.8E\n", ms[is],miota[is],mIpol[is],mItor[is],mpprime[is]);
    }
    else  {
      fprintf(fp,"%17.10E %16.8E %16.8E %16.8E %16.8E %16.8E\n",ms[is],miota[is],mIpol[is],mItor[is],mpprime[is],mg00[is]);
    }
    fputs(mHdr3.c_str(),fp);
    int m=0,n=0;
    for(m=0; m<=Mmax; m++) {
      int lowN = (mNewFormat&&m==0)?0:-Nmax;
      lowN = tokamakEQDSK?0:lowN;
      for(n=lowN; n<= Nmax; n++) {
        fprintf(fp,"% 5d %4d ", m, n);
        if(vmecCoordinates) {
          //fprintf(fp,"%16.8E %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",
          //   Rmn(is,m,n),Zmn(is,m,n),bmn(is,m,n),lmn(is,m,n),gmn(is,m,n),bContraTheta_mn(is,m,n),bContraFi_mn(is,m,n));
          fprintf(fp,"%16.8E %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",           //see p, 2012jul1
             Rmn(is,m,n),Zmn(is,m,n),Phimn(is,m,n),bmn(is,m,n),lmn(is,m,n),gmn(is,m,n),bContraTheta_mn(is,m,n),bContraFi_mn(is,m,n));
        } else if(fpCoordinates) {
          fprintf(fp,"%16.8E %16.8E %16.8E %16.8E\n",Rmn(is,m,n),Zmn(is,m,n),bmn(is,m,n),lmn(is,m,n));
        }
        else {
          fprintf(fp,"%16.8E %16.8E %16.8E %16.8E\n",Rmn(is,m,n),Zmn(is,m,n),Phimn(is,m,n),bmn(is,m,n));
        }
      }
    }
  }
  fclose(fp);
  return(true);
}

//****************************************************************************
// Remesh and write bc-file with surfaces equidistantly distributed on r_eff.
// Nsmax is the number of surfaces.
// Don't remesh if Nsmax==0, use original s-mesh
// include magnetic axis if magAxis==true
// if varNumberMods==true then store variable number of harmonics for each surface
bool CStconfig::writeasciiReduced(const char * fname,int Nsmax,double newVolume,double newLastSurfaceLabel,
                                  bool varNumberModes,bool magAxis,bool circular,bool append) const
{
  bool remesh = (Nsmax!=0)?true:false; // don't remesh if Nsmax==0
  Nsmax = (Nsmax==0)?mNs0:Nsmax;
  Nsmax = abs(Nsmax);
  //if(newLastSurfaceLabel > 1) newLastSurfaceLabel = 1;
  double sLast = newLastSurfaceLabel==1?this->sLast:newLastSurfaceLabel;

  if(newLastSurfaceLabel!=1) remesh = true;

  FILE * fp = append?fopen(fname,"a+"):fopen(fname,"w");
  if(fp==NULL) return(false);

  std::string hdr2 = mHdr2;
  bool vmecCoord = vmecCoordinates;
  bool fpCoord = fpCoordinates;
  if(remesh&&vmecCoord) hdr2=mHdr2a;

  if(append)
    fprintf(fp, "CC*****************************************\n");

//#define TEST_FP

#ifdef TEST_FP
  vmecCoord = false;
  fpCoord = true;
  remesh = true;
  Nsmax = (Nsmax==0)?mNs0:Nsmax;
  fprintf(fp,"CC original file '%s' saved by 'mcviewer'\n",fname_name(mfname.c_str()));
  fprintf(fp,"CC %s\n",fpComment); 
  hdr2 = "       s              iota\n";
#else
  writeIdentificationComment(fp);
#endif

  // change size
  double vol = const_cast<CStconfig*>(this)->Volume(1);
  double volumeFactor = fabs(newVolume>0? (newVolume/vol) : 1); 
  double linearFactor = pow(volumeFactor,1./3.); 
  double areaFactor   = pow(linearFactor,2.); 

  bool mNewFormatCC=true;
  fprintf(fp, "CC truncated at Bmn/B00=%g\n",mmodBtol);
  if(newVolume>0) {
    fprintf(fp, "CC scaled to new size with linear scaling factor %g\n",linearFactor);
    //fprintf(fp, "CC new volume %gm**3\n",newVolume);
  }
  if(newLastSurfaceLabel!=1) {
    fprintf(fp, "CC new last flux surface is at old s=%g \n",newLastSurfaceLabel);
  }
  if(newVolume>0||newLastSurfaceLabel!=1) {
    double vol1 = fabs( const_cast<CStconfig*>(this)->Volume(newLastSurfaceLabel) );
    vol1 *= volumeFactor;
    double rmax = r(newLastSurfaceLabel)*linearFactor;
    fprintf(fp, "CC new volume %gm**3\n",vol1);
    fprintf(fp, "CC new minor radius %gm\n",rmax);
  }
  fprintf(fp, "CC %s \n", varNumberModes?"variable number of harmonics":" ");

// put also the old header
  fprintf(fp, "CC******below is the header of the original file*************\n");
  if(mNewFormat) {
    for(unsigned int i=0; i<mComments.size(); i++)
      fputs(mComments[i].c_str(),fp);
  }

  int Mmax = mmin(1,mM);
  int Nmax = mmin(1,mN);
  int i,is,m,n;
  for(i=0; i<=mLCMSindex; i++) {
    Mmax = mmax(Mmax,mpM[i]);
    Nmax = mmax(Nmax,mpN[i]);
  }
  int Mmax0 = Mmax;
  int Nmax0 = Nmax;

  double Fluxmax = Flux(newLastSurfaceLabel)*areaFactor;
  double rmax = r(newLastSurfaceLabel)*linearFactor;

  fputs(mHdr1.c_str(),fp);
  fprintf(fp,"%4d %4d %4d %4d ",Mmax0,Nmax0,Nsmax,mNp);
  fprintf(fp,"%13.6E %8.5f %8.5f\n",Fluxmax,rmax,mR0*linearFactor);

  if(varNumberModes) {
    Mmax = mmin(1,mM);
    Nmax = mmin(1,mN);
  }

  double s1st_  = CStconfig::s1st;
  if(magAxis) s1st_ = 0;

  double dr = (sqrt(sLast)-sqrt(s1st_))/(Nsmax-1);
  for(i=0; i<Nsmax; i++) {  // loop over surfaces
    double r = sqrt(s1st_) + i*dr;
    ////if(i==0&&magAxis) r = dr/1024; // skip magnetic axis
    double s = (i==Nsmax-1)?sLast:(r*r);
    double sFull;
    fputs(hdr2.c_str(), fp);
    if(mNewFormatCC) fputs(mHdr21.c_str(),fp);
    double iota_,Ip_,It_,pp_,g00_;
    if(remesh) {
      is = SearchSindx(s); //search of the segment in which s lies, needed for varNumberModes
      sFull = s;
      iota_=iota(s);
      Ip_  =Ip(s)/mNp;
      It_  =It(s);
      pp_  =const_cast<CStconfig*>(this)->pp(s);
      g00_ =g00(s);
    }
    else {
      is = m1stSindex+i;
      s    =ms[is];
      sFull=vmecCoordinates?msFull[is]:s;
      iota_=miota[is];
      Ip_  =mIpol[is]*mBnorm;
      It_  =mItor[is]*mBnorm;
      pp_  =mpprime[is]*mBnorm*mBnorm;
      g00_ =mg00[is];
    }
    Ip_   *= linearFactor;  
    It_   *= linearFactor;  
    ///test to scale rminor g00_  *= areaFactor;
    g00_  *= volumeFactor;
    if(vmecCoord)
      fprintf(fp,"%17.10E ", sFull/newLastSurfaceLabel);
    if(fpCoordinates2) {
      fprintf(fp,"%17.10E %16.8E\n", s/newLastSurfaceLabel,iota_);
    }
    else if(fpCoordinates5) {
      fprintf(fp,"%17.10E %16.8E %16.8E %16.8E %16.8E\n", s/newLastSurfaceLabel,iota_,Ip_,It_,pp_);
    }
    else {
      fprintf(fp,"%17.10E %16.8E %16.8E %16.8E %16.8E %16.8E\n", s/newLastSurfaceLabel,iota_,Ip_,It_,
                                                    pp_*newLastSurfaceLabel,g00_*newLastSurfaceLabel);
    }
    fputs(mHdr3.c_str(),fp);
    if(varNumberModes) {         // variable number of harmonics
      int Mma=mpM[is];
      int Nma=mpN[is];
      Mmax = mmax(Mmax,Mma);    // # of harmonics must be monotonically increasing
      Nmax = mmax(Nmax,Nma);    //  function of surface number;
      Mmax = mmin(Mmax,Mmax0);  //  we need this to avoid ill-behaved spline;
      Nmax = mmin(Nmax,Nmax0);  //  see also calculation of spline coef. in addSurfaces()
    }
    for(m=0; m<=Mmax; m++) {
      int lowN = (mNewFormatCC&&m==0)?0:-Nmax;
      lowN = tokamakEQDSK?0:lowN;
      for(n=lowN; n<= Nmax; n++) {
        double r(0),z(0),p(0),b(0),l(0),g(0),bp(0),bt(0);
        if(remesh) {
          r=Rmn_sp(s,m,n);
          z=Zmn_sp(s,m,n);
          p=(fpCoord|vmecCoord)?0:Phimn_sp(s,m,n);
          l=fpCoord?lmn_sp(s,m,n):0;
          b=Bmn_sp(s,m,n);
          if(vmecCoord) {
            l=lmn_sp(s,m,n);
            g=gmn_sp(s,m,n);
            bp=BcontrPolmn_sp(s,m,n);
            bt=BcontrTormn_sp(s,m,n);
            p=Phimn_sp(s,m,n);
          }
        }
        else {
          r=Rmn(is,m,n);
          z=Zmn(is,m,n);
          p=(fpCoord|vmecCoord)?0:Phimn(is,m,n);
          l=fpCoord?lmn(is,m,n):0;
          b=Bmn(is,m,n);   // Bmn() already has normalizing factor fabs(mBnorm);
          if(vmecCoord) {
            g=gmn  (is,m,n);          // gmnc   on half-mesh
            l=lmn  (is,m,n);          // lmns   on half-mesh
            bp=bContraTheta_mn(is,m,n)*mBnorm; // contravariant poloidal component of B on half-mesh
            bt=bContraFi_mn   (is,m,n)*mBnorm; // contravariant toroidal component of B on half-mesh
            p=Phimn(is,m,n);
          }
        }
#if 0   ///test to scale rminor without scaling major radius
        if(m!=0||n!=0) r *= linearFactor;
        if(m!=0||n!=0) z *= linearFactor;
        g *= areaFactor;
#else
        r *= linearFactor;
        z *= linearFactor;
        g *= volumeFactor;
#endif
        if(varNumberModes) {    // variable number of harmonics
          if(fpCoord&&r==0&&z==0&&l==0&&b==0) continue;
          //else if(vmecCoord&&r==0&&z==0&&l==0&&b==0&&g==0&&bp==0&&bt==0) continue;
          else if(vmecCoord&&r==0&&z==0&&l==0&&b==0&&g==0&&bp==0&&bt==0&&p==0) continue; //see p, 2012jul1
          else if(r==0&&z==0&&p==0&&b==0) continue;
          else if(circular&&n>0) continue;
        }
        fprintf(fp,"% 5d %4d ", m, n);
        if(vmecCoord) {
          fprintf(fp,"%16.8E %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",r,z,p,b,l,g,bp,bt);  //see p, 2012jul1
        }
        else if(fpCoord) {
          fprintf(fp,"%16.8E %16.8E %16.8E %16.8E\n",r,z,b,l);
////todo  fprintf(fp,"%16.8E %16.8E %16.8E %16.8E\n",r,z,l,b);
        }
        else {
          if(circular&&n>0) {r=z=p=b=0;}
          fprintf(fp,"%16.8E %16.8E %16.8E %16.8E\n",r,z,p,b);
        }
      }
    }
  }
  fclose(fp);
  return(true);

  #undef TEST_FP
}

//****************************************************************************
// Extract |B|-spectrum on the  surface 's' into file 'fname'
//
bool CStconfig::writeasciiBmod(const char * fname,double s) const
{
  FILE * fp = fopen(fname,"w");
  if(fp==NULL) return false;

  writeIdentificationComment(fp);

  fprintf(fp,"CC |B|-spectrum extracted from the '%s' by 'MConf::CStconfig'\n",fname_name(mfname.c_str()));
  fprintf(fp,"CC truncation level Bmn/B00=%g\n",mmodBtol);
  fprintf(fp,"CC %s\n",mSignJac_s>0?" Right handed (s,theta,phi) system":"Left handed (s,theta,phi) system");
  fputs("CC\n",fp);

  int m,n,M,N;
  const_cast<CStconfig*>(this)->getSsegMN(s,M,N);

  double Fluxmax = Flux(1.);
  double rmax = r(1.);
  double B00  = Bmn_sp(s,0,0);

  fputs("  m0b  n0b      flux(a)    a,[m]    R,[m]\n",fp);
  fprintf(fp,"% 4d %4d ",M,N);
  fprintf(fp,"%13.6E %8.5f %8.5f\n",Fluxmax,rmax,mR0);

  fprintf(fp,"%13s %13s %13s\n", "s","iota","B00");
  fprintf(fp,"%13.6E %13.6E %13.6E\n", s,iota(s),B00);
  fputs("    m    n    Bmn/B00\n", fp);
  for(m=0; m<=M; m++) {
    int lowN = (m==0)?0:-N;
    for(n=lowN; n<= N; n++) {
      fprintf(fp,"% 5d %4d ", m, n);
      ////double r=Rmn_sp(s,m,n);
      ////double z=Zmn_sp(s,m,n);
      ////double p=(fpCoordinates|vmecCoordinates)?0:Phimn_sp(s,m,n);
      double b=Bmn_sp(s,m,n)/B00;
      fprintf(fp,"%16.8E\n",b);
    }
  }
  fclose(fp);
  return(true);
}


//****************************************************************************
// load configuration from binary file 'fname'
// if( ff!=NULL) then read from opened file 'fname'
// if( ff!=NULL) then don't close file
template <class T> bool CStconfig::loadbin(const char * fname, T szType, FILE *ff)
{
  FILE * fp = (ff!=NULL)?ff:fopen(fname,"rb");
  if(fp==NULL) return(mOK=false);

  mfname = fname;

  mNewFormat = true;

  char temp[256];
  if(fgets(temp,250,fp)==NULL) {
    if(ff==NULL) fclose(fp);
    return(mOK=false);
  }

  if(strstr(temp,tokamakComment1)!=NULL) tokamakSymmetryFluxCoordinates=true;
  if(strstr(temp,vmecComment)!=NULL)     vmecCoordinates=true;
  if(strstr(temp,fpComment)!=NULL)   fpCoordinates=true;
  
  tokamakEQDSK  = tokamakSymmetryFluxCoordinates;
  tokamakConfig = tokamakEQDSK; 

  bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);

  int ws;
  fread(&ws, sizeof(int),1,fp);

  int ns,m,n;
  fread(&m, sizeof(int),1,fp);
  fread(&n, sizeof(int),1,fp);
  fread(&ns,sizeof(int),1,fp);
  fread(&mNp,sizeof(int),1,fp);

  T w[10]; // double or float w[10]
  size_t  wsize = sizeof(w[0]);
  fread(w,wsize,3,fp);
  mFlux = w[0];
  ma    = w[1];
  mR0   = w[2];

//MagAxis*********************************
  {  // does the file contain the magnetic axis? 
    fpos_t pos0;
    fgetpos(fp,&pos0); // read position in file
    int count = vmecCoordinates?7:6;
    double sfull=1,sHalf=1;
    fread(w, wsize,count,fp);
    sHalf = w[0];
    if(vmecCoordinates) sfull  = w[6];
    if(sHalf==0) 
      m1stSindex = 0;   // set to zero if axis 
    fsetpos(fp,&pos0);  // restore position
  }
//MagAxis*********************************

  if(!resize(ns,m,n)) {
    if(ff==NULL) fclose(fp);
    return(mOK=false);
  }

  for(int i=0; i<mNs0; i++) {
    int is = m1stSindex+i;
    int count = vmecCoordinates?7:6;
    fread(w, wsize,count,fp);
    ms.rw()[is]      = w[0];
    miota.rw()[is]   = w[1];
    mIpol.rw()[is]   = w[2];
    mItor.rw()[is]   = w[3];
    mpprime.rw()[is] = w[4];
    mg00.rw()[is]    = w[5];
    if(vmecCoordinates)
      msFull.rw()[is]  = w[6];
    int m,n;
    fread(&m,sizeof(int),1,fp);
    fread(&n,sizeof(int),1,fp);
    mpM.rw()[is] = m;
    mpN.rw()[is] = n;
    for(m=0; m<=mpM[is]; m++)
     for(n=-mpN[is]; n<= mpN[is]; n++) {
        int cnt=0;
        T tmp = 0;
        int ncnt = vmecCoordinates?8:4; 
        ncnt = fpCoordinates?4:ncnt;
        fread(w, wsize,ncnt,fp);
        Rmn(is,m,n)   = w[cnt++];
        Zmn(is,m,n)   = w[cnt++];
        if(boozerCrd||vmecCoordinates||tokamakEQDSK) {
          Phimn(is,m,n) = w[cnt++];
        }
        bmn(is,m,n)   = w[cnt++];
        if(fpCoordinates) {
          lmn(is,m,n) = w[cnt++];
        }
        else if(vmecCoordinates) { 
          lmn(is,m,n)             = w[cnt++];
          gmn(is,m,n)             = w[cnt++];
          bContraTheta_mn(is,m,n) = w[cnt++];
          bContraFi_mn(is,m,n)    = w[cnt++];
        }
     }
  }
#if 1    // read LCMS
  fgets(temp,8,fp);
  if(!feof(fp)&&(strncmp(temp,"<LCMS>",6)==0)) {
    double wd[6];
    fread(wd, sizeof(wd[0]),6,fp);
    lcms.mRmin = wd[0];
    lcms.mRmax = wd[1];
    lcms.mZmin = wd[2];
    lcms.mZmax = wd[3];
    lcms.mshift= wd[4];
    lcms.mdPhi = wd[5];
    if(lcms.vertex.read0(szType,fp) != 0) {
      lcms.Bfield.read0(szType,fp);
    }
    if( !(lcms.vertex.isOK()&&lcms.Bfield.isOK()) ) {
      lcms.vertex.clear();
      lcms.Bfield.clear();
    }
  }
#endif
  if(ff==NULL) fclose(fp);
  return(mOK=true);
}

//****************************************************************************
// write configuration into file 'fname' in binary format
// if( ff!=NULL) then continue writing into file 'fname'
//
template <class T> bool CStconfig::writebin(const char * fname, T szType, FILE *ff) const
{
  T w[10]; // double or float w[10]
  size_t  wsize = sizeof(w[0]);

  FILE * fp = (ff!=NULL)?ff:fopen(fname,"wb");
  if(fp==NULL) return false;
  bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);

  fprintf(fp,"<binary magnetic configuration file in W7X-format>%s",(wsize==sizeof(double))?"double":"float ");
  if(tokamakSymmetryFluxCoordinates)
    fprintf(fp,";%s",tokamakComment1);
  if(vmecCoordinates)
    fprintf(fp,";%s",vmecComment);
  if(fpCoordinates)
    fprintf(fp,";%s",fpComment);
  fprintf(fp,"\n");

  fwrite(&wsize, sizeof(int),1,fp);

  int Mmax = mmin(1,mM);
  int Nmax = mmin(1,mN);
  int Mmax0 = Mmax;
  int Nmax0 = Nmax;
  int i;
  for(i=0; i<mNs0; i++) {
    int is = m1stSindex+i;
    Mmax0 = mmax(Mmax0,mpM[is]);
    Nmax0 = mmax(Nmax0,mpN[is]);
  }

  fwrite(&Mmax0, sizeof(int),1,fp);
  fwrite(&Nmax0, sizeof(int),1,fp);
  fwrite(&mNs0,sizeof(int),1,fp);
  fwrite(&mNp,sizeof(int),1,fp);

  w[0] = T(mFlux);
  w[1] = T(ma);
  w[2] = T(mR0);
  fwrite(w, wsize,3,fp);

  for(i=0; i<mNs0; i++) {
    int is = m1stSindex+i;
    w[0] = T(ms[is]);
    w[1] = T(miota[is]);
    w[2] = T(mIpol[is]);
    w[3] = T(mItor[is]);
    w[4] = T(mpprime[is]);
    w[5] = T(mg00[is]);
    int count=6;
    if(vmecCoordinates) {
      w[count] = T(msFull[is]);
      count++;
    }
    fwrite(w, wsize,count,fp);
    {         // variable number of harmonics
      int Mma=mpM[is];
      int Nma=mpN[is];
      Mmax = mmax(Mmax,Mma);    // # of harmonics must be monotonically increasing
      Nmax = mmax(Nmax,Nma);    //  function of surface number;
      Mmax = mmin(Mmax,Mmax0);  //  we need this to avoid ill-behaved spline;
      Nmax = mmin(Nmax,Nmax0);  //  see also calculation of spline coef. in addSurfaces()
    }
    fwrite(&Mmax,sizeof(int),1,fp);
    fwrite(&Nmax,sizeof(int),1,fp);
    int m,n;
    for(m=0; m<=Mmax; m++)
    for(n=-Nmax; n<= Nmax; n++) {
        int cnt=0; 
        w[cnt++] = T(Rmn(is,m,n));
        w[cnt++] = T(Zmn(is,m,n));
        if(boozerCrd||vmecCoordinates||tokamakEQDSK) {
          w[cnt++] = T(Phimn(is,m,n));
        }
        w[cnt++] = T(bmn(is,m,n));
        if(fpCoordinates) {
          w[cnt++] = T(lmn(is,m,n));
        }
        else if(vmecCoordinates) {    
          w[cnt++] = T(lmn(is,m,n));
          w[cnt++] = T(gmn(is,m,n));
          w[cnt++] = T(bContraTheta_mn(is,m,n));
          w[cnt++] = T(bContraFi_mn(is,m,n));
        }
        fwrite(w, wsize,cnt,fp);
    }
  }
#if 1 // write LCMS
  if(lcms.vertex.isOK()) {
    fprintf(fp,"<LCMS>\n");
    double wd[6];
    wd[0] = lcms.mRmin;
    wd[1] = lcms.mRmax;
    wd[2] = lcms.mZmin;
    wd[3] = lcms.mZmax;
    wd[4] = lcms.mshift;
    wd[5] = lcms.mdPhi;
    fwrite(wd, sizeof(wd[0]),6,fp);
    lcms.vertex.write0(szType,fp);
    lcms.Bfield.write0(szType,fp);
  }
#endif
  if(ff==NULL) fclose(fp);
  return true;
}

/// helper function for write(), see also C3dMesh::writeMeshbin();
bool CStconfig::writebin4(const char * fname, FILE *ff) const {return writebin(fname, static_cast<float>(1), ff);};
/// helper function for write(), see also C3dMesh::writeMeshbin();
bool CStconfig::writebin8(const char * fname, FILE *ff) const {return writebin(fname, static_cast<double>(1),ff);};
/// helper function for load(), see also C3dMesh::loadMeshbin();
bool CStconfig::loadbin4(const char * fname, FILE *ff) {return loadbin(fname, static_cast<float>(1), ff);};
/// helper function for load(), see also C3dMesh::loadMeshbin();
bool CStconfig::loadbin8(const char * fname, FILE *ff) {return loadbin(fname, static_cast<double>(1),ff);};

#pragma inline_depth(255)

//****************************************************************************
// Linear interpolation in u. x and y are input arrays with sizes == 2
// on return y(u)
static inline double interp(double u, const double *x, const double *y)
{
  return(y[0] + (y[1]-y[0])/(x[1]-x[0])*(u-x[0]));
}

//****************************************************************************
// Linear interpolation in u between i0 & i1. x and y are input arrays
// on return y(u)
static inline double interp(double u, const double *x, const double *y, int i0, int i1)
{
  return(y[i0] + (y[i1]-y[i0])/(x[i1]-x[i0])*(u-x[i0]));
}

//****************************************************************************
// ajust Fourier coefficients and create radial meshes
void CStconfig::ajustFourierCoeffMeshes()
{
  if(!mOK) return;
// set harmonics bmn = bmn + bm,-n  for m=0, n>0
// set harmonics bmn = 0 for m=0, n<0
  int i;
  for(i=0; i<mNs0; i++) {
    int is=m1stSindex+i;
    for(int n=1; n<= mN; n++) {
      Rmn(is,0,n)  += Rmn(is,0,-n);  Rmn(is,0,-n)=0;
      Zmn(is,0,n)  -= Zmn(is,0,-n);  Zmn(is,0,-n)=0;
      Phimn(is,0,n)-= Phimn(is,0,-n);Phimn(is,0,-n)=0;
      bmn(is,0,n)  += bmn(is,0,-n);  bmn(is,0,-n)=0;
      if(vmecCoordinates) {
        bContraTheta_mn(is,0,n) += bContraTheta_mn(is,0,-n);  bContraTheta_mn(is,0,-n)=0;
        bContraFi_mn(is,0,n)    += bContraFi_mn(is,0,-n);     bContraFi_mn(is,0,-n)=0;
        gmn(is,0,n)  += gmn(is,0,-n);  gmn(is,0,-n)=0;  //*cos(m*theta-Np*n*phi) - stellarator symmetric part
        lmn(is,0,n)  -= lmn(is,0,-n);  lmn(is,0,-n)=0;  //*sin(m*theta-Np*n*phi) - stellarator symmetric part
      } else if(fpCoordinates) {
        lmn(is,0,n)  -= lmn(is,0,-n);  lmn(is,0,-n)=0;  //*sin(m*theta-Np*n*phi) - stellarator symmetric part
      }
    }
  }
// note: all unassigned harmonics are zero, see resize(...), mCOEF[i]=0

  s1st  = ms[m1stSindex];
  sLast = ms[mLCMSindex];
  if(m1stSindex) {
    double dr = sqrt(s1st)/m1stSindex;
    for(int is=0; is<m1stSindex; is++) ms.rw()[is]=square(is*dr); // add points near axis
  }
  for(int is=mLCMSindex+1; is<mNs; is++) ms.rw()[is]=ms[is-1]+0.3/mNsAdd; // add points after LCMS
  ms.rw()[mNs] = ms[mNs-1]+10;   // extra point, far from last surface, see CStconfig::resize(.., ns++; // add one radial point
  if(mNs==mLCMSindex+1) ms.rw()[mNs] = ms[mNs-1]+0.1; 


  if(expandLCMS) {
    for(i=0; i<mNs; i++) {
      double s = ms[i];
        double mapFunc = 2*smoothHeaviside(s-1,0.15)+1;  // 2/(1+exp(-8*(s-1)/0.15))+1;
      ms.rw()[i]=s*mapFunc;
    }
    if(mNs==mLCMSindex+1) ms.rw()[mNs] = ms[mNs-1]+0.1; 
  }
  
  for(i=0; i<=mNs; i++) mreff.rw()[i]=sqrt(ms[i]);

  if(vmecCoordinates) { // full mesh for VMEC
    double s1st_  = msFull[m1stSindex];
    if(m1stSindex) {
      double dr = sqrt(s1st_)/m1stSindex;
      for(int is=0; is<m1stSindex; is++) msFull.rw()[is]=square(is*dr); // add points near axis
    }
    for(int is=mLCMSindex+1; is<mNs; is++) msFull.rw()[is]=ms[is-1]+0.3/mNsAdd; // add points after LCMS
    msFull.rw()[mNs] = msFull[mNs-1]+10;   // extra point, far from last surface, see CStconfig::resize(.., ns++; // add one radial point
    if(mNs==mLCMSindex+1) msFull.rw()[mNs] = msFull[mNs-1]+0.1; 
    if(expandLCMS) {
      for(i=0; i<mNs; i++) {
        double s = msFull[i];
        double mapFunc = 2*smoothHeaviside(s-1,0.15)+1;  // 2/(1+exp(-8*(s-1)/0.15))+1;
        msFull.rw()[i]=s*mapFunc;
      }
      if(mNs==mLCMSindex+1) msFull.rw()[mNs] = msFull[mNs-1]+0.1; 
    }
    for(int i=0; i<=mNs; i++) mreffFull.rw()[i]=sqrt(msFull[i]);
  }

  if(!expandLCMS) return;


  for(int i=0; i<mNs; i++) {
    double s = msFull[i];
    double mapFunc = 2*smoothHeaviside(s-1,0.15)+1;  // 2/(1+exp(-8*(s-1)/0.15))+1;
    for(int m=0; m<=mM; m++) 
      for(int n=-mN; n<=mN; n++) {
        if(m==0) continue;
        Rmn(i,m,n)  *= mapFunc;
        Zmn(i,m,n)  *= mapFunc;
      }
  }

}

//****************************************************************************
// add surfaces near magnetic axis
// calculate spline coefficients
void CStconfig::addSurfacesNearAxis()
{
  if(!mOK) return;
  int is=0,m,n;

// assign correct signes
#if CORRECT_SIGNES > 0
      //#pragma message( "===> CORRECT_SIGNES VMECcoordinates?")
// we must derive the correct sign of the maxFlux despite what was stored in the file
// see also NJacobian2()
//  jac  = mu0*(Ipol*mNp+iota*Itor)/(modB*modB*twopi*twopi); Boozer Jacobian in coordinates (flux,theta,phi)
  is = m1stSindex;
  double signJac_F=1;
  if(vmecCoordinates||fpCoordinates)
    signJac_F = setSign(1.,mIpol[is]+miota[is]*mItor[is]) ; // get the sign of Jacobian in coordinates (Flux,theta,phi)
  else if(tokamakEQDSK)
    signJac_F = setSign(1.,mIpol[is]);
  else  // Boozer Jacobian
    signJac_F = setSign(1.,mIpol[is]+miota[is]*mItor[is]);   // get the sign of jacobian(flux,theta,phi)

  mFlux = setSign(mFlux,mSignJac_s*signJac_F);
  for(is=0; is<=mNs; is++) mg00.rw()[is]  = mSignJac_s*fabs(mg00[is]);
#endif

// linearly extrapolate to s==0
  mItor.rw()[0]   = 0;  // tor. current must be zero on axis
  for(is=0; is<m1stSindex; is++) {
    miota.rw()[is]   = extrapolate0(ms[is],ms.constArray()+m1stSindex,miota.constArray()+m1stSindex);
    mpres.rw()[is]   = extrapolate0(ms[is],ms.constArray()+m1stSindex,mpres.constArray()+m1stSindex);
    mIpol.rw()[is]   = extrapolate0(ms[is],ms.constArray()+m1stSindex,mIpol.constArray()+m1stSindex);
    mpprime.rw()[is] = extrapolate0(ms[is],ms.constArray()+m1stSindex,mpprime.constArray()+m1stSindex);
    mg00.rw() [is]   = extrapolate0(ms[is],ms.constArray()+m1stSindex,mg00.constArray()+m1stSindex);
    mItor.rw()[is]   = interp(ms[is],ms.constArray(),mItor.constArray(),  0,m1stSindex);
  }

// insert surfaces between 0 and m1stSindex
#if 1
{
  double w[8][4];
  for(is=0; is<m1stSindex; is++) {
    for(m=0; m<=mM; m++) {
      int lowN=(m==0)?0:-mN;
      double mpow=(m>10)?10:m;
      for(n=lowN; n<= mN; n++) {
        for(int i=0; i<3; i++) {          // prepare for linear extrapolation
          double erm  = 1./pow(mreff[i+m1stSindex],mpow);
          double erm1 = vmecCoordinates?(1./pow(mreffFull[i+m1stSindex],mpow)) : erm;
          w[0][i]= Rmn  (i+m1stSindex,m,n)*erm1;
          w[1][i]= Zmn  (i+m1stSindex,m,n)*erm1;
          w[2][i]= Phimn(i+m1stSindex,m,n)*erm;
          w[3][i]= bmn  (i+m1stSindex,m,n)*erm;
          if(vmecCoordinates) {
            w[4][i]= bContraTheta_mn(i+m1stSindex,m,n)*erm;
            w[5][i]= bContraFi_mn   (i+m1stSindex,m,n)*erm;
            w[6][i]= gmn            (i+m1stSindex,m,n)*erm;
            w[7][i]= lmn            (i+m1stSindex,m,n)*erm;
          }
          if(fpCoordinates) {
            w[7][i]= lmn            (i+m1stSindex,m,n)*erm;
          }
        }
        double rm = mreff[is];
        double rm1= vmecCoordinates?mreffFull[is] : rm;
        rm = mpow?pow(rm,mpow) : 1;  // some compilers don't handle correctly pow(0.,0.)
        rm1= mpow?pow(rm1,mpow) : 1;
        const double * sp = vmecCoordinates?(msFull.constArray()):(ms.constArray());
        Rmn(is,m,n)  = extrapolate0(sp[is],sp+m1stSindex,w[0])*rm1;
        Zmn(is,m,n)  = extrapolate0(sp[is],sp+m1stSindex,w[1])*rm1;
        Phimn(is,m,n)= extrapolate0(ms[is],ms.constArray()+m1stSindex,w[2])*rm;
        bmn(is,m,n)  = extrapolate0(ms[is],ms.constArray()+m1stSindex,w[3])*rm;
        if(vmecCoordinates) {
          bContraTheta_mn(is,m,n)  = extrapolate0(ms[is],ms.constArray()+m1stSindex,w[4])*rm;
          bContraFi_mn   (is,m,n)  = extrapolate0(ms[is],ms.constArray()+m1stSindex,w[5])*rm;
          gmn            (is,m,n)  = extrapolate0(ms[is],ms.constArray()+m1stSindex,w[6])*rm;
          lmn            (is,m,n)  = extrapolate0(ms[is],ms.constArray()+m1stSindex,w[7])*rm;
        }
        if(fpCoordinates) {
          lmn            (is,m,n)  = extrapolate0(ms[is],ms.constArray()+m1stSindex,w[7])*rm;
        }
      }
    }
  }
}
#else
  // 14-Feb-08: not used since 2004
//*****BEGIN  Smothing splining
{ // try to smooth
  int splineSM(int n,double x[],double y[],double b[],double c[],double d[]);
  int ns = SearchSindx(0.16);
  int nspl = ns - m1stSindex+1;
  // Create spline coefficients
  double *y  = new double[5*nspl];  if(y==NULL) return;
  double *b = y   + nspl;
  double *c = b   + nspl;
  double *d = c   + nspl;
  double *x = d   + nspl;
  double yp;

#define CALCOEFFmn0(bmn)         \
  for(i=m1stSindex;i< ns;i++) {  \
    y[i-m1stSindex]=bmn(i,m,n);  \
    x[i-m1stSindex]=mreff[i];    \
  }                              \
  splineSM(nspl-1,x,y,b,c,d);    \
  for(is=0; is<ns; is++)         \
    bmn(is,m,n)= seval (nspl-1, mreff[is],x,y,b,c,d, yp);
// end of macro CALCOEFFmn0(bmn)

  m=0;
  for(n=0; n<= mN; n++) {
    CALCOEFFmn0(Rmn);
    CALCOEFFmn0(Zmn);
    CALCOEFFmn0(Phimn);
    CALCOEFFmn0(bmn);
  }


#define CALCOEFFmn1(bmn)           \
  for(i=m1stSindex;i< ns;i++) {    \
    y[i-m1stSindex+1]=bmn(i,m,n);  \
    x[i-m1stSindex+1]=mreff[i];    \
  }                                \
  y[0]=bmn(0,m,n);                 \
  x[0]=0;                          \
  splineSM(nspl,x,y,b,c,d);        \
  for(is=0; is<ns; is++)           \
    bmn(is,m,n)= seval (nspl, mreff[is],x,y,b,c,d, yp);
// end of macro CALCOEFFmn1(bmn)

  for(m=1; m<= mM; m++)
    for(n=-mN; n<= mN; n++) {
      CALCOEFFmn1(Rmn);
      CALCOEFFmn1(Zmn);
      CALCOEFFmn1(Phimn);
      CALCOEFFmn1(bmn);
    }
  delete[] y;
}
//*****END  Smothing splining
#endif
}

//////////////////////////////////////////
// add surfaces after LCMS
void CStconfig::addSurfacesAfterLCMS()
{
  if(!mOK) return;
  int i,is=0,m,n;
// linearly extrapolate to extra point  with s>1
// use linear interpolation on s for all quantity except A_1n

  int nLCMS = mLCMSindex;

  if(mpprime[nLCMS]==0)  //correct error in bc-file
    mpprime.rw()[nLCMS] = extrapolateLCMS(ms[nLCMS],ms.constArray(),mpprime.constArray(), nLCMS);

//use one sided approximation for derivatives in point nLCMS, see extrapolateLCMS
  const double *sp = vmecCoordinates?msFull.constArray():ms.constArray();       //pointer to s-array for Rmn and Zmn extrapolation
  const double *rp = vmecCoordinates?mreffFull.constArray():mreff.constArray(); //pointer to r-array for Rmn and Zmn extrapolation

  const int h=3;   // distance between point for calculating derivative in extrapolateLCMS
  double w[8][3];
  double xs[3];
  double xr[3];
  double xsp[3];
  double xrp[3];
  for(i=0; i<3; i++) {          // prepare for extrapolation
    int idx = nLCMS-2*h+i*h;
    xs [i]= ms[idx];
    xr [i]= mreff[idx];
    xsp[i]= sp[idx];
    xrp[i]= rp[idx];
  }

  for(is=nLCMS+1; is<=mNs; is++) {
    miota.rw()[is]   = extrapolateLCMS(ms[is],ms.constArray(),miota.constArray()   , nLCMS+1,h);
    mpres.rw()[is]   = extrapolateLCMS(ms[is],ms.constArray(),mpres.constArray()   , nLCMS+1,h);
    mIpol.rw()[is]   = extrapolateLCMS(ms[is],ms.constArray(),mIpol.constArray()   , nLCMS+1,h);
    mItor.rw()[is]   = extrapolateLCMS(ms[is],ms.constArray(),mItor.constArray()   , nLCMS+1,h);
    mpprime.rw()[is] = extrapolateLCMS(ms[is],ms.constArray(),mpprime.constArray() , nLCMS+1,h);
    mg00.rw()[is]    = extrapolateLCMS(ms[is],ms.constArray(),mg00.constArray()    , nLCMS+1,h);
    for(m=0; m<= mM; m++) {
      int lowN = (m==0)?0:-mN;
      for(n=lowN; n<= mN; n++) {
        for(i=0; i<3; i++) {          // prepare for extrapolation
          int idx = nLCMS-2*h+i*h;
          w[0][i]= Rmn  (idx,m,n);
          w[1][i]= Zmn  (idx,m,n);
          w[2][i]= Phimn(idx,m,n);
          w[3][i]= bmn  (idx,m,n);
          if(vmecCoordinates) {
            w[4][i]= bContraTheta_mn(idx,m,n);
            w[5][i]= bContraFi_mn   (idx,m,n);
            w[6][i]= gmn            (idx,m,n);
            w[7][i]= lmn            (idx,m,n);
          }
          if(fpCoordinates) {
            w[7][i]= lmn            (idx,m,n);
          }
        }
//???        if(m==1) { // linear extrapolation on reff for A_1n
        if(m<0) {
          Rmn(is,m,n)  = extrapolateLCMS(rp   [is],xrp,w[0],3);
          Zmn(is,m,n)  = extrapolateLCMS(rp   [is],xrp,w[1],3);
          Phimn(is,m,n)= extrapolateLCMS(mreff[is],xr,w[2],3);
          bmn(is,m,n)  = extrapolateLCMS(mreff[is],xr,w[3],3);
          if(vmecCoordinates) {
            bContraTheta_mn(is,m,n) = extrapolateLCMS(mreff[is],xr,w[4],3);
            bContraFi_mn   (is,m,n) = extrapolateLCMS(mreff[is],xr,w[5],3);
            gmn            (is,m,n) = extrapolateLCMS(mreff[is],xr,w[6],3);
            lmn            (is,m,n) = extrapolateLCMS(mreff[is],xr,w[7],3);
          }
          if(fpCoordinates) {
            lmn            (is,m,n) = extrapolateLCMS(mreff[is],xr,w[7],3);
          }
        }
        else {
          Rmn(is,m,n)  = extrapolateLCMS(sp[is],xsp,w[0],3);
          Zmn(is,m,n)  = extrapolateLCMS(sp[is],xsp,w[1],3);
          Phimn(is,m,n)= extrapolateLCMS(ms[is],xs,w[2],3);
          bmn(is,m,n)  = extrapolateLCMS(ms[is],xs,w[3],3);
          if(vmecCoordinates) {
            bContraTheta_mn(is,m,n) = extrapolateLCMS(ms[is],xr,w[4],3);
            bContraFi_mn   (is,m,n) = extrapolateLCMS(ms[is],xr,w[5],3);
            gmn            (is,m,n) = extrapolateLCMS(ms[is],xr,w[6],3);
            lmn            (is,m,n) = extrapolateLCMS(ms[is],xr,w[7],3);
          }
          if(fpCoordinates) {
            lmn            (is,m,n) = extrapolateLCMS(ms[is],xr,w[7],3);
          }
        }
      }
    }
  }
  // the iota can change sign due to extrapolation to s>1
  double dIota = miota[nLCMS]-miota[nLCMS-1];
  int iotaSign = sign(miota[nLCMS]);
  bool correctionNeeded = (iotaSign>0&&dIota<0)||(iotaSign<0&&dIota>0);
  if(correctionNeeded) {
    int i0 = 0;
    for(is=mNs; is>nLCMS; is--)
      if(iotaSign != sign(miota[is])) i0=is;
    if(i0) {
      //double y1=miota[i0-1];
      //double y2=y1/10;
      //double x1=ms[i0-1];
      //double x2=ms[mNs];
      //double yp=(y2-y1)/(x2-x1);
      //for(is=i0; is<=mNs; is++)
      //  miota[is] = y1 + (ms[is]-x1)*yp;
      //if(i0<mNs)
      //  miota[i0] = 0.5*(miota[i0-1]+ miota[i0+1]);
  i0 = nLCMS+1;
      double x1=ms   [i0-2], x2=ms   [i0-1];
      double y1=miota[i0-2], y2=miota[i0-1];
      double yp=(y2-y1)/(x2-x1);
      double b=-yp/y2;
      double a=y2;
      for(is=i0; is<=mNs; is++)
        miota.rw()[is] = a*::exp(-b*(ms[is]-x2));
      miota.rw()[mNs-1] = miota[mNs-2];
      miota.rw()[mNs]   = miota[mNs-2];
    }
  }
}

void CStconfig::createSplineCoeff()
{
  if(!mOK) return;
  int i,is=0,m,n;

#if INTERPOLATION_ORDER < 3
//  for(is=0; is<=mNs; is++) mVolume[is]=VolumeInt(ms[is]);
  if(!vprimeNeeded)
     for(is=1,mVolume.rw()[0]=0; is<=mNs; is++) mVolume.rw()[is]=mVolume[is-1]+VolumeInt(ms[is-1],ms[is]);

// calculate poloidal flux using iota
// s_poloidal-grid -- normalized poloidal flux
  for(is=1, msPol.rw()[0]=0; is<=mNs; is++) msPol.rw()[is]=msPol[is-1]+FluxPoloidal(ms[is-1],ms[is]);
  double fluxPolMax1 = sPol(1.);
  for(is=0; is<=mNs; is++) msPol.rw()[is] /= fluxPolMax1;

  return;
#endif
  ////clock_t tStart = clock();  //timing

// Create spline coefficients, use r_eff as abscissa
  int nx, ns=mNs; // ns++;

  double *y  = new double[2*ns];  if(y==NULL) return;
  double *y2 = y   + ns;
  double y0Prime=1e31;  // spline Boundary condition at r=0

#define CALCOEFCmn(Rmn,rmeshPtr)  {   \
  for(i=0,nx=-1;i< ns;i++) {          \
    y[i]=Rmn(i,m,n);                  \
    if(y[i]!=0&&nx==-1) nx=i;         \
    y2[i]=0;                          \
  }                                   \
/* nx means # of zero points to be skiped from splining */  \
  if(nx>0) nx--;                      \
  y0Prime = m!=1?0:1e31;              \
  if(nx>=0) /* all y[i]==0  if nx<0*/ \
    splineNRC(ns-nx, rmeshPtr+nx,y+nx,y2+nx,y0Prime); /* TODO spline BC */ \
  for(i=0; i<ns; i++) {                \
/*  Rmn(i,m,n)   = y[i]; if smooth spline*/ \
    Rmn##_d2(i,m,n)=y2[i];               \
  }                                      \
  Rmn##_d2(ns-2,m,n)=0;  /* linear interpolation for the last segment */ \
  Rmn##_d2(ns-1,m,n)=0;} /* linear interpolation for the last segment */
// end of macro CALCOEFCmn(Rmn)

  if(vmecCoordinates)
    for(m=0; m<= mM; m++)
      for(n=-mN; n<= mN; n++) {
        const double *rmeshFull = mreffFull.constArray();
        CALCOEFCmn(Rmn,rmeshFull);
        CALCOEFCmn(Zmn,rmeshFull);
        const double *rmesh = mreff.constArray();
 //28.11.2010       CALCOEFCmn(Phimn,rmesh);
        CALCOEFCmn(Phimn,rmesh);
        CALCOEFCmn(bmn,rmesh);
        CALCOEFCmn(bContraTheta_mn,rmesh);
        CALCOEFCmn(bContraFi_mn,rmesh);
        CALCOEFCmn(gmn,rmesh);
        CALCOEFCmn(lmn,rmesh);
      }
  else if(fpCoordinates) 
    for(m=0; m<= mM; m++)
      for(n=-mN; n<= mN; n++) {
        const double *rmesh = mreff.constArray();
        CALCOEFCmn(Rmn,rmesh);
        CALCOEFCmn(Zmn,rmesh);
 //28.11.2010 CALCOEFCmn(Phimn,rmesh);
        CALCOEFCmn(bmn,rmesh);
        CALCOEFCmn(lmn,rmesh);
      }
  else
    for(m=0; m<= mM; m++)
      for(n=-mN; n<= mN; n++) {
        const double *rmesh = mreff.constArray();
        CALCOEFCmn(Rmn,rmesh);
        CALCOEFCmn(Zmn,rmesh);
        if(!tokamakEQDSK) CALCOEFCmn(Phimn,rmesh);
        CALCOEFCmn(bmn,rmesh);
      }

  //CALCOEFF(miota,ns);
  miota.create(this);
  mpres.create(this);

  if(!polCurrentNeeded) mIpol.create(this);
  if(!torCurrentNeeded) mItor.create(this);
  if(!vprimeNeeded)     mg00.create(this);

  delete[] y;

#undef CALCOEFCmn

  ////std::cerr <<"splining time="<<double(clock()-tStart)/CLOCKS_PER_SEC<<std::endl;
}

//****************************************************************************
void CStconfig::calculateIpolItorFromSij(int nthead)
{
  if(!mOK) return;
  if(nthead<=0) nthead = nThreads;
  Threads2<CStconfig,CStconfig> threads(this, &CStconfig::calculateIpolItorFromSijExe,this,"calculateIpolItorFromSij");  
  int low = (m1stSindex==0)?1:m1stSindex;
  int upper = mLCMSindex;
  threads.run(low, upper, nthead); 
  if(low==1) mItor.w()[0] = 0;
  if(low==1&&polCurrentNeeded) mIpol.w()[0] = mIpol[1];
}

void CStconfig::calculateIpolItorFromSijExe(CStconfig * mc, int low, int upper, bool setup)
{
  if(setup) {  // setup part; it allocates arrays and initialize them
    int size = upper-low+1;
  }
  else
  {  // exe part; it will be called from threads to fill arrays
    ////int npol = tokamakConfig?300:180; // disabled 2011 Nov. 02
    ////int ntor = tokamakConfig?4:180;
    int npol=180, ntor=250;
    for(int i=low; i<=upper; i++) {
      double s = mc->ms[i];
      Vector3d *sm = mc->SMatrix(s,npol,ntor); // S-matrix in (flux,th,ph)-coordinates
      if(torCurrentNeeded) mItor.w()[i]=sm[1][0];     // toroidal current
      if(polCurrentNeeded) mIpol.w()[i]=sm[2][0]/mNp; // poloidal current on a period
    }
  }
}

//****************************************************************************
void CStconfig::calculateIpolItor(int nthead)
{
  if(!mOK) return;

  if(nthead<=0) nthead = nThreads;
  Threads2<CStconfig,CStconfig> threads(this, &CStconfig::calculateIpolItorExe,this,"calculateIpolItor");  
  int low = (m1stSindex==0)?1:m1stSindex;
  int upper = mLCMSindex;
  threads.run(low, upper, nthead); 
  if(low==1) mItor.w()[0] = 0;
  if(low==1&&polCurrentNeeded) mIpol.w()[0] = mIpol[1];
}

void CStconfig::calculateIpolItorExe(CStconfig * mc, int low, int upper, bool setup)
{
  if(setup) {  // setup part; it allocates arrays and initialize them
    int size = upper-low+1;
  }
  else
  {  // exe part; it will be called from threads to fill arrays
    int Ntheta = tokamakConfig?400:180;
    int Nphi   = tokamakConfig?4:180;
    for(int i=low; i<=upper; i++) {
      if(torCurrentNeeded) mItor.w(i) = mc->calculateItor(i, Ntheta);      // toroidal current
      if(polCurrentNeeded) mIpol.w(i) = mc->calculateIpol(i, Nphi)/mNp;    // poloidal current on a period
    }
  }
}

//****************************************************************************
void CStconfig::calculateVp(int nthead)
{
  if(!mOK) return;

  if(nthead<=0) nthead = nThreads;
  Threads2<CStconfig, CStconfig> threads(this, &CStconfig::calculateVpExe, this,"calculateVp");  
  int low = (m1stSindex==0)?1:m1stSindex;
  int upper = mLCMSindex;
  threads.run(low, upper, nthead); 
  if(low==1) mg00.w(0) = mg00[1];
}

//template <class mConf> void CStconfig::calculateVpExe(mConf * mc, int low, int upper, bool setup)
void CStconfig::calculateVpExe(CStconfig * mc, int low, int upper, bool setup)
{
  if(setup) {  // setup part; it allocates arrays and initialize them
    int size = upper-low+1;
  }
  else
  {  // exe part; it will be called from threads to fill arrays
    for(int i=low; i<=upper; i++) {
      int npol=180, ntor=250;
      double s = mc->ms[i];
      double vp = mc->VprimeIntJ(s,npol,ntor);
      mg00.w(i) = vp/mc->mNp;
    }
  }
}


//****************************************************************************
// Calculate currents  Itor = 1/mu0*Integral(dtheta*B_theta); 
//    where B_theta = B*e_theta, e_theta = dX/dtheta; e_theta is the covariant basis vector in poloidal direction
double CStconfig::calculateItor(int is, int Ntheta)
{
  int n = Ntheta;
  if(n<4) n=4;
  n=n+n%2;        //must be even for simpson integration
  double dt = twopi/n;
  double Itor = 0;
  for(int k=0; k<=n; k++) {                  // Simpson integration over poloidal angle
    double w = (k==0||k==n)?1:(2*(1+k%2));   // weight for Simpson integration
    double s = ms[is];
    double theta=k*dt;
    Vector3d mag(s,theta,0);
    NJacobian(mag);
    Vector3d Bcova = getBcova();
    Itor += w*Bcova[1];
  }
  return (Itor*dt/3/mu0);
}

//****************************************************************************
// Calculate currents  Ipol = 1/mu0*Integral(dphi*B_phi); 
//    where B_phi = B*e_phi, e_phi = dX/dphi; e_phi is the covariant basis vector in toroidal direction
double CStconfig::calculateIpol(int is, int Nphi)
{
  int n = Nphi;
  if(n<4) n=4;
  n=n+n%2;        //must be even for simpson integration
  double dphi = twopi/mNp/n;
  double Ipol = 0;
  for(int k=0; k<=n; k++) {                  // Simpson integration over toroidal angle
    double w = (k==0||k==n)?1:(2*(1+k%2));   // weight for Simpson integration
    double s = ms[is];
    if(s==0) s=ms[1]/2;
    double phi=k*dphi;
    Vector3d mag(s,pi,phi);
    NJacobian(mag);
    Vector3d Bcova = getBcova();
    Ipol += w*Bcova[2];
  }
  return (Ipol*dphi/3/mu0)*mNp;
}

//****************************************************************************
// Calculate currents
// file with fileFormat W7X_FMT_VMEC_FP containes neither Vprime nor currents 
// file with fileFormat LHD  containes no Vprime 
//void CStconfig::calculateCurrentsAndVprime()
void CStconfig::calculateCurrents()
{
  if(!mOK||!(torCurrentNeeded|polCurrentNeeded)) return;
  double mBnormSave = mBnorm;
  mBnorm = 1;
  
  TESTPRINT( std::cerr <<"CStconfig: Calculating currents..."<<std::endl ); 

  double  truncSav = truncate();
  if(truncate()<1.e-8) truncate(1.e-8);

  int Ntheta = tokamakConfig?400:180;
  int Nphi   = tokamakConfig?4:180;

  if(torCurrentNeeded|polCurrentNeeded) calculateIpolItor(); // fill mItor and mIpol using multithreading
  //if(nThreads) 
  //  calculateIpolItor(); // fill mItor and mIpol using multithreading
  //else {
  //  int ifirst = (m1stSindex==0)?1:m1stSindex;
  //  if(torCurrentNeeded) for(int is=ifirst; is<=mLCMSindex; is++) mItor.w()[is] = calculateItor(is, Ntheta);      // toroidal current
  //  if(polCurrentNeeded) for(int is=ifirst; is<=mLCMSindex; is++) mIpol.w()[is] = calculateIpol(is, Nphi)/mNp;    // poloidal current on a period
  //  if(m1stSindex==0)  mIpol.w()[0] = mIpol[1];
  //  if(m1stSindex==0)  mg00.w()[0] = mg00[1];
  //}

// linearly extrapolate to s==0
  mItor.w()[0]   = 0;  // tor. current must be zero on axis
  for(int is=0; is<m1stSindex; is++) {
    if(torCurrentNeeded) mItor.w()[is]  = interp(ms[is],ms.constArray(),mItor.constArray(),  0,m1stSindex);
    if(polCurrentNeeded) mIpol.w()[is]  = extrapolate0(ms[is],ms.constArray()+m1stSindex,mIpol.constArray()+m1stSindex);
  }
// add surfaces after LCMS, // linearly extrapolate to extra points with s>1
  const int h=3;            // distance between point for calculating derivative in extrapolateLCMS
  for(int is=mLCMSindex+1; is<=mNs; is++) {
    if(torCurrentNeeded) mItor.w()[is] = extrapolateLCMS(ms[is],ms.constArray(),mItor.constArray()   , mLCMSindex+1,h);
    if(polCurrentNeeded) mIpol.w()[is] = extrapolateLCMS(ms[is],ms.constArray(),mIpol.constArray()   , mLCMSindex+1,h);
  }

#if INTERPOLATION_ORDER == 3
  if(polCurrentNeeded) mIpol.create(this);  // Create spline coefficients
  if(torCurrentNeeded) mItor.create(this);
#endif

  polCurrentNeeded = false;
  torCurrentNeeded = false;

  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
  mBnorm = mBnormSave;
}

//****************************************************************************
// Calculate Vprime
// Vprime is negative if (s,theta,phi)-coordinate system is left handed
void CStconfig::calculateVprime()
{
  if(!mOK||!vprimeNeeded) return;
  vprimeNeeded = false;

  TESTPRINT( std::cerr <<"CStconfig: Calculating dV/ds..."<<std::endl );

  double  truncSav = truncate();
  if(truncate()<1.e-6) truncate(1.e-6);
  
  int ifirst = (m1stSindex==0)?1:m1stSindex;
  if(nThreads>1) 
    calculateVp(); // fill mg00 using multithreading
  else {
    int Ntheta = tokamakConfig?500:180;
    int Nphi   = tokamakConfig?4:250;
    for(int is=ifirst; is<=mLCMSindex; is++) mg00.w()[is] = mSignJac_s*fabs(VprimeIntJ(ms[is],Ntheta,Nphi))/mNp;
    if(ifirst==1)  mg00.w()[0] = mg00[1];
  }

// linearly extrapolate to s==0
  if(ifirst!=1) for(int is=0; is<m1stSindex; is++)  mg00.w()[is]  = extrapolate0(ms[is],ms.constArray()+m1stSindex,mg00.constArray()+m1stSindex);

#if INTERPOLATION_ORDER < 3
  for(int is=1,mVolume.w()[0]=0; is<=mLCMSindex; is++) mVolume.w()[is]=mVolume[is-1]+VolumeInt(ms[is-1],ms[is]);
#else
//use one sided approximation for derivatives in point nLCMS, see extrapolateLCMS
// add surfaces after LCMS, // linearly extrapolate to extra points with s>1
  const int h=3;            // distance between point for calculating derivative in extrapolateLCMS
  for(int is=mLCMSindex+1; is<=mNs; is++) mg00.w()[is] = extrapolateLCMS(ms[is],ms.constArray(),mg00.constArray(), mLCMSindex+1,h);
  mg00.create(this); // Create spline coefficients, use r_eff as abscissa
#endif
  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
}

//****************************************************************************
void CStconfig::calculateSMatrix(int size)
{
  if(!mOK) return;
  double mBnormSave = mBnorm;
  mBnorm = 1;
  // calculate currents if needed
  calculateCurrents();  // current calculation was in initAfterLoad() (immediately after file reading)
                        //  ... so I keep this order
  double truncSav = truncate();
  if(!tokamakConfig) truncate(1.e-4);

  Threads2<CStconfig,CStconfig> threads(this, &CStconfig::calculateSMatrixExe,this,"calculateSMatrix");  
  int low = 0, upper = low+size-1;
  threads.run(low, upper, nThreads); 

  mBnorm = mBnormSave;
  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
  SMatrixReady  = true;
}

void CStconfig::calculateSMatrixExe(CStconfig * mc, int low, int upper, bool setup)
{
  if(!mOK) return;
  if(setup) {  // setup part; it allocates arrays and initialize them
    int size = upper-low+1;
    bool oK = msK.init(size);
    oK &= mS11.init(size); 
    oK &= mS12.init(size);
    oK &= mS21.init(size);
    oK &= mS22.init(size);
    if(!oK) {
      msK.clear();
      mS11.clear();
      mS12.clear();
      mS21.clear();
      mS22.clear();
      return;
    }
    int np = msK.size();   // # of points along s
    double ds = ms[mLCMSindex-3]/(np-4);
    double s3[3];
    s3[0] = ms[mLCMSindex-2];
    s3[1] = (mNs-1==mLCMSindex)?ms[mLCMSindex-1]:ms[mLCMSindex];
    s3[2] = ms[mNs-1];
    for(int i=0; i<np; i++) {
      msK.rw()[i] = ds*i;
      if(i>=np-3) msK.rw()[i] = s3[i-(np-3)];
    }
  }
  else
  {  // exe part; it will be called from threads to fill arrays
    if(msK.empty()) return;
    for(int i=low; i<=upper; i++) {
      Vector3d *sm = mc->SMatrix(msK[i]);
      mS11.w(i)= i==0?0:sm[1][1];
      mS12.w(i)= i==0?0:sm[1][2];
      mS21.w(i)= i==0?0:sm[2][1];
      mS22.w(i)= sm[2][2];
    }
  }
}

//****************************************************************************
// The method tabulate arrays for finding Sij(s).
// nth, nph are given for twopi
// \note This method is obsolete, see calculateSMatrix()
void CStconfig::calculateSMatrixArrays(int nth, int nph)
{
  if(!mOK) return;
  bool oK = msK.init(mNsK);
  oK &= mS11.init(mNsK); 
  oK &= mS12.init(mNsK);
  oK &= mS21.init(mNsK);
  oK &= mS22.init(mNsK);
  if(!oK) return;

  double  truncSav = truncate();
  if((nth<=180||nph<=200)&&!tokamakConfig) truncate(1.e-4);
  double mBnormSave = mBnorm;
  mBnorm = 1;

  int np = msK.size();  // # of points along s

  double ds = ms[mLCMSindex-3]/(np-4);
  double s3[3];
  s3[0] = ms[mLCMSindex-2];
  s3[1] = (mNs-1==mLCMSindex)?ms[mLCMSindex-1]:ms[mLCMSindex];
  s3[2] = ms[mNs-1];
  for(int i=0; i<np; i++) {
    msK.rw()[i] = ds*i;
    if(i>=np-3) msK.rw()[i] = s3[i-(np-3)];
    Vector3d *sm = SMatrix(msK[i],nth,nph);
    mS11.rw()[i]= i==0?0:sm[1][1];
    mS12.rw()[i]= i==0?0:sm[1][2];
    mS21.rw()[i]= i==0?0:sm[2][1];
    mS22.rw()[i]= sm[2][2];
  }

  mBnorm = mBnormSave;
  if(truncSav != truncate()) truncate(truncSav); // restore truncation level

  SMatrixReady  = true;
  return;
}

//****************************************************************************
// Tabulate arrays for finding the values of Bmin,Bmax
//
// nth, nph are given for twopi
void CStconfig::calculateBminmax()
{
  int size = mreK.size(); // # of points along s
  if(!mOK) return;

  double mBnormSave = mBnorm;
  mBnorm = 1;
  double truncSav = truncate();
  if(!tokamakConfig) truncate(1.e-4);

  Threads2<CStconfig,CStconfig> threads(this, &CStconfig::calculateBminmaxExe,this,"calculateBminmax");  
  int low = 0, upper = low+size-1;
  threads.run(low, upper, nThreads); 
  mBminGlobal = 1e15;
  mBmaxGlobal = -1e15;
  for(int i=0; i<size; i++) {
    mBminGlobal = mmin(mBmin[i],mBminGlobal);
    mBmaxGlobal = mmax(mBmax[i],mBmaxGlobal);
  }
  mBnorm = mBnormSave;
  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
  BminmaxReady=true;
}

void CStconfig::calculateBminmaxExe(CStconfig * mc, int low, int upper, bool setup)
{
  if(!mOK) return;
  if(setup) {  // setup part; it allocates arrays and initialize them
    int np = upper-low+1; // # of points along s
    bool oK = mBmin.init(np);
    oK &= mBmax.init(np);
    if(!oK) return;

    double dr = sqrt(ms[mLCMSindex-3])/(np-4);
    double s3[3];
    s3[0] = ms[mLCMSindex-2];
    s3[1] = (mNs-1==mLCMSindex)?ms[mLCMSindex-1]:ms[mLCMSindex];
    s3[2] = ms[mNs-1];
    for(int i=0; i<np; i++) {
      mreK.rw()[i] = dr*i;
      if(i>=np-3) {
        double s = s3[i-(np-3)];
        mreK.rw()[i] = sqrt(s);
      }
    }
  }
  else
  {  // exe part; it will be called from threads to fill arrays
    int nth=180, nph=250;
    for(int i=low; i<=upper; i++) {
      double x = mc->mreK[i];
      double Bmin, Bmax;
      mc->getBminmax(x*x, Bmin, Bmax, nth, nph);
      mBmin.w(i) = Bmin;
      mBmax.w(i) = Bmax;
    }
  }
}

//****************************************************************************
// min and max B on the surface s
void CStconfig::getBminmax(double s, double &Bmin, double &Bmax,int nth,int nph)
{
  int itmax = tokamakConfig?500:nth;  // theta-points (poloidal)
  int ipmax = tokamakConfig?1:nph;    // phi-points (toroidal)

  itmax = mmax(nth,itmax);
  itmax = mmax(10, itmax);

  double dtheta = twopi/itmax;
  double dphi   = mPeriod/ipmax;

  Bmin =  1e15;
  Bmax = -1e15;
  for(int i=0; i<itmax; i++) {
    double theta = i*dtheta;
    for(int k=0; k<ipmax; k++) {
      Vector3d magCoord(s,theta,k*dphi);
      double b = B(magCoord);
      Bmin = mmin(Bmin,b);
      Bmax = mmax(Bmax,b);
    }
  }
}

double CStconfig::rgraz(double s) 
{
  averagedData.create(this); 
  return averagedData.r_graz(s);
}

double CStconfig::drgraz_dr(double s) 
{
  if(s<1e-14) s=1e-14;
  return 2*s/r(s)/GradsAvrg(s);     // dr_graz/dr_vmec
}

double CStconfig::drdkes_dr(double s) 
{
  if(s<1e-14) s=1e-14;
  double B00 = Bmn_sp(s,0,0);
  double B00p= Bmn_sp_prime(s,0,0);
  return rdkes(s)/r(s)*(1-s*B00p/B00); // dr_dkes/dr_vmec
}

double CStconfig::rdkes2(double s) // minor effective radius[m], another DKES definition
{
  averagedData.create(this); 
  return averagedData.r_dks(s);
}  

double CStconfig::drdkes2_dr(double s) 
{
  if(s<1e-14) s=1e-14;
  double B00 = Bmn_sp(s,0,0);
  return r(s)/rdkes2(s)*fabs(Flux(1.)/(pi*ma*ma*B00)); // dr_dkes2/dr_vmec
}

double CStconfig::GradsAvrg(double s) 
{
  averagedData.create(this); 
  return averagedData.Grads(s);
}

double CStconfig::Grads2Avrg(double s) 
{
  averagedData.create(this); 
  return averagedData.Grads2(s);
}

double CStconfig::Grads2overB2Avrg(double s) 
{
  averagedData.create(this); 
  return averagedData.Grads2overB2(s)/(mBnorm*mBnorm);
}

//****************************************************************************
void CStconfig::calculateGradsAvrg(int size)
{
  if(!mOK) return;
  double mBnormSave = mBnorm;
  mBnorm = 1;

  Threads2<CStconfig,CStconfig> threads(this, &CStconfig::calculateGradsAvrgExe,this, "calculateGradsAvrg");  
  int low = 0, upper = low+size-1;
  threads.run(low, upper, nThreads); 
  mBnorm = mBnormSave;
  averagedData.set_r_graz();
}

//****************************************************************************
void CStconfig::calculateGradsAvrgExe(CStconfig * mc, int low, int upper, bool setup)
{
  if(setup) {  // setup part; it allocates arrays and initialize them
    int size = upper-low+1;
    averagedData.resize(size); // detach data if they are shared 
    if(averagedData.x.empty()) return;
    for(int i=1; i<size; i++) {
      double dx = (1-1e-7)/(size-2);
      double x  = 1e-7 + dx*(i-1);   
      averagedData.x.rw()[i] = x;  // x=sqrt(s)
    }
    averagedData.x.rw()[0] = 0;
    averagedData.Btruncation = truncate(); // save truncation level of B
  }
  else 
  {  // exe part; it will be called from threads to fill arrays
    if(averagedData.x.empty()) return;
    for(int i=low; i<=upper; i++) {
      double x = mc->averagedData.x[i];
      double s = x*x;
      averagedData.eB00.w(i) = fabs( mc->Flux(1.)/( pi*mc->Bmn_sp(s,0,0) ) );
      if(i==0) continue;
      Vector3d g = mc->getGradsAverages(s);
      averagedData.grads.w(i)  = g[0]/x;
      averagedData.grads2.w(i) = g[1]/s;
      averagedData.grads2overB2.w(i) = g[2]/s;
    }
    if(low==0) {
      averagedData.grads.w(0)  = averagedData.grads.w(1);
      averagedData.grads2.w(0) = averagedData.grads2.w(1);
      averagedData.grads2overB2.w(0) = averagedData.grads2overB2.w(1);
    }
  }
}

//****************************************************************************
//  <|grads|>
Vector3d CStconfig::getGradsAverages(double s,int ithmax,int ipmax)
{
  volatile ctStack::saveState state(ct);
  if(s<1e-14) s=1e-14;
  double dtheta, dphi;
  adjustNumbersForSimpson(ithmax,ipmax,dtheta,dphi);

  //Simpson integration: flux surface averaging
  double gradsavrg =0;    // <|grads|>  is the flux average
  double grads2avrg=0;    // <grads^2>  is the flux average
  double grads2B2avrg=0;  // <grads^2/B^2>  is the flux average
  double dVds=0;          // dV/ds=integral( jacobian  dthete dphi)
  for(int i=0; i<=ithmax; i++) {
    double theta = i*dtheta;
    double wt = (i==0||i==ithmax)?1:(2*(1+i%2));    // weight for Simpson integration
    for(int k=0; k<=ipmax; k++) {
      double w = (k==0||k==ipmax)?1:(2*(1+k%2));   // weight for Simpson integration
      w *= wt;

      Vector3d magCoord(s,theta,k*dphi);
      NJacobian(magCoord);
      double jac   = fabs(getJ()); //  Jacobian in coordinates (s,theta,phi) ~ Ipol/B^2
      double grads2 = getGrads().abs2();
      double B = getmodB();
      dVds       += w*jac;
      gradsavrg  += w*jac*sqrt(grads2);    
      grads2avrg += w*jac*grads2; 
      grads2B2avrg += w*jac*grads2/square(B);       
    }
  }
  return Vector3d(gradsavrg/dVds,grads2avrg/dVds,grads2B2avrg/dVds);
}

//****************************************************************************
// Tabulate arrays for finding the values of Bmin,Bmax,<B>,<B^2>,
// fraction of trapped particles, and function
// f(s,x) = <sqrt(1-x*B/Bmax)> 0<=x<=1, 0<=s<=1
//
void CStconfig::calculateFtrappedBminmax()
{
  int size = mreK.size(); // # of points along s
  if(!mOK) return;

  double mBnormSave = mBnorm;
  mBnorm = 1;
  double truncSav = truncate();
  if(!tokamakConfig) truncate(1.e-4);

  Threads2<CStconfig,CStconfig> threads(this, &CStconfig::calculateFtrappedBminmaxExe,this,"calculateFtrappedBminmax");  
  int low = 0, upper = low+size-1;
  threads.run(low, upper, nThreads); 

  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
  ftrappedReady=true;
  BminmaxReady=true;
  {
    int n=512; //n must be even for simpson integration
    double sum(0), ds=1./n;
    for(int k=0; k<=n; k++) {  //Simpson integration
      double w = (k==0||k==n)?1:(2*(1+k%2));   // weight for Simpson integration
      sum += w*Vprime(k*ds)*B2avrg(k*ds);
    }
    mB2VolAvrg = (sum*ds/3)/Volume(1.); 
  }
  mBminGlobal = 1e15;
  mBmaxGlobal = -1e15;
  for(int i=0; i<size; i++) {
    mBminGlobal = mmin(mBmin[i],mBminGlobal);
    mBmaxGlobal = mmax(mBmax[i],mBmaxGlobal);
  }
  mBnorm = mBnormSave;
  ////for(int i=0; i<size; i++) std::cerr<<mreK[i]<<std::endl;
}

// nth, nph are given for twopi
void CStconfig::calculateFtrappedBminmaxExe(CStconfig * mc, int low, int upper, bool setup)
{
  if(!mOK) return;  
  if(setup) {  // setup part; it allocates arrays and initialize them
    int np = upper-low+1;   
    bool oK = mFtrap.init(np);
    oK &= mBmin.init(np);
    oK &= mBmax.init(np);
    oK &= mBavrg.init(np);
    oK &= mB2avrg.init(np);
    oK &= mpavr_sx.init (np,mNxK); 
    oK &= mpavrI_sx.init (np,mNxK);
    if(!oK) return;
    
    int Nx1 = mxtr.size()-1; //Nx must be even for simpson integration
    double *x = mxtr.array();
    double dx=1./Nx1;
    for(int j=0; j<=Nx1; j++) x[j]=j*dx;

    double dr = sqrt(ms[mLCMSindex-3])/(np-4);
    double s3[3];
    s3[0] = ms[mLCMSindex-2];
    s3[1] = (mNs-1==mLCMSindex)?ms[mLCMSindex-1]:ms[mLCMSindex];
    s3[2] = ms[mNs-1];
    for(int i=0; i<np; i++) {
      mreK.rw()[i] = dr*i;
      if(i>=np-3) {
        double s = s3[i-(np-3)];
        mreK.rw()[i] = sqrt(s);
      }
    }
  }
  else
  {  // exe part; it will be called from threads to fill arrays
    int Nx = mc->mxtr.size(); //mNxK; // # of points along x
    int nth=180, nph=250;
    CArray1d<double> Ipavrg_1(Nx);
    CArray1d<double> pavrg2(Nx);
    for(int i=low; i<=upper; i++) {
      double x = mc->mreK[i];   // = sqrt(s)
      double Ftrapped,Bmin,Bmax,Bavr,B2avr;
      ////double * pavrg2  =  const_cast<double *>mpavr_sx[i];
      ////double * Ipavrg_1 = const_cast<double *>mpavrI_sx[i];
      mc->getFluxAverages (x*x,Ftrapped,Bmin,Bmax,Bavr,B2avr,Ipavrg_1, pavrg2,nth,nph);                              
      mFtrap.w(i) = square(Ftrapped); // fraction of trapped part. scales as sqrt(reff), see also interpolation
      mBmin.w(i)  = Bmin;
      mBmax.w(i)  = Bmax;
      mBavrg.w(i) = Bavr;
      mB2avrg.w(i) = B2avr;
      for(int j=0; j<Nx; j++) {
        mpavr_sx.w(i,j)  = pavrg2[j];  
        mpavrI_sx.w(i,j) = square(Ipavrg_1[j]);
      }
    }
  }
}


//****************************************************************************
// for A.Mishchenko
void CStconfig::calculateFluxAverages1(int nth, int nph)
{
  if(!mOK) return;
  volatile ctStack::saveState state(ct);

  bool oK = mB_2avrg.init(mNsK);
  oK &= mR2avrg.init(mNsK);
  oK &= mR_2avrg.init(mNsK);
  oK &= mGrads2B_2avrg.init(mNsK);
  oK &= mNeoPolarztnPass.init(mNsK);
  oK &= mNeoPolarztnTrap.init(mNsK);
  oK &= mb12.init(mNsK);    // <(1-b)^0.5>
  oK &= mb32.init(mNsK);    // <(1-b)^1.5>
  oK &= mpavrb_sx.init(mNsK,mNxK);
  oK &= mpavrb2.init(mNxK);
  if(!oK) return;

  double  truncSav = truncate();
  if((nth<=180||nph<=200)&&!tokamakConfig) truncate(1.e-4);
  double mBnormSave = mBnorm;
  mBnorm = 1;

  int Nx = mxtr.size(); //mNxK; // # of points along x
  int np = mreK.size(); //mNsK; // # of points along s
  int Nx1 = mxtr.size()-1; //Nx1 must be even for simpson integration
  double *x = mxtr.array();
  double dx=1./Nx1;
  for(int j=0; j<=Nx1; j++) x[j]=j*dx;

  double dr = sqrt(ms[mLCMSindex-3])/(np-4);
  double s3[3];
  s3[0] = ms[mLCMSindex-2];
  s3[1] = (mNs-1==mLCMSindex)?ms[mLCMSindex-1]:ms[mLCMSindex];
  s3[2] = ms[mNs-1];
  for(int i=0; i<np; i++) {
    mreK.rw()[i] = dr*i;
    double s = mreK[i]*mreK[i];
    if(i>=np-3) {
      s = s3[i-(np-3)];
      mreK.rw()[i] = sqrt(s);
    }
    getFluxAverages1(s,mB_2avrg.rw()[i],mR2avrg.rw()[i],mR_2avrg.rw()[i],mGrads2B_2avrg.rw()[i],
                    mb12.rw()[i],mb32.rw()[i],mNeoPolarztnPass.rw()[i],mNeoPolarztnTrap.rw()[i],nth,nph);
    for(int j=0; j<Nx; j++)
      mpavrb_sx(i,j) = mpavrb2[j];
  }

  mBnorm = mBnormSave;
  if(truncSav != truncate()) truncate(truncSav); // restore truncation level

  Averages1Ready=true;
  return;
}

//****************************************************************************
// Calculate the fraction of trapped particles and min and max B on the surface s
void CStconfig::getFtrappedBminmax(double s,double &Ftrapped,double &Bmin,double &Bmax,
                                   double &Bavr,double &B2avr,
                                   int itmax,int ipmax)
{
   volatile ctStack::saveState state(ct);
   getFluxAverages(s,Ftrapped,Bmin,Bmax,Bavr,B2avr,mpavrI,mpavr2,itmax,ipmax);
}

//****************************************************************************
// Calculate the fraction of trapped particles and min and max B on the surface s
//
//  The averaged fraction of trapped particles on a magnetic flux
//  surface is only defined by the magnetic field strength:
//
//    FT = 1 - 3/4 * <B^2/Bmax^2> * int (0->1) xdx/<g(x)>
//
//  where <g(x)> is the flux surface average of sqrt(1-x*B/Bmax)
//  with Bmax being the maximum of the magnetic field strength and
//  <A> = int int A*jacobian*dth*dph  /  int int jacobian*dth*dph
//  the definition of the flux surface average in magnetic coordinates.
//
// @return Ftrapped is the fraction of trapped particles
//         Bmin  is the Bmin
//         Bmax  is the Bmax
//         Bavr  is the <B>
//         B2avr is the <B^2>
//         mxtr[i]   = x[i]
//         pavrg2[i] = <sqrt(1-x[i]*B/Bmax)>^2
//         Ipavrg_1[i] = integral(from x[i] to 1) (dx/<sqrt(1-x[i]*b)>)
//
// itmax, ipmax are the number of integration points
//  itmax       theta-points (poloidal)
//  ipmax       phi-points (toroidal)
//
void CStconfig::getFluxAverages(double s,double &Ftrapped,
                                  double &Bmin,double &Bmax,
                                  double &Bavr,
                                  double &B2avr,  // <B^2>
                                  CArray1d<double> &Ipavrg_1,
                                  CArray1d<double> &pavrg2,
                                  int ithmax,int ipmax)
{
  if(s<1e-8) s=1e-8;

  int Nx = mxtr.size()-1; //Nx must be even for simpson integration
  const double *x = mxtr.constArray();

  CArray1d<double> mpavr;     // = <sqrt(1-x[i]*b)> for given surface,array[mNxK]
  mpavr.init(0,Nx,0.0); 
  double *I = mpavr.array();

  int i,j,k;
  double dx=1./Nx;

  double dtheta, dphi;
  adjustNumbersForSimpson(ithmax,ipmax,dtheta,dphi);

  pBase::ng2Matrix<double> jacobian(ithmax+1,ipmax+1);
  pBase::ng2Matrix<double> bn(ithmax+1,ipmax+1);
  Bmin= 1e15;
  Bmax=-1e15;
  for(i=0; i<=ithmax; i++) {
    double theta = i*dtheta;
    for(k=0; k<=ipmax; k++) {
      Vector3d magCoord(s,theta,k*dphi);
//1      NJacobian(magCoord);
//1      jacobian[i][k] = fabs(getJ()); //  Jacobian in coordinates (s,theta,phi)
//1      bn[i][k] = getmodB();
      bn[i][k] = B(magCoord);
      if(tokamakEQDSK||vmecCoordinates||fpCoordinates)
        jacobian[i][k] = fabs(Jacobian(magCoord)); //  Jacobian in coordinates (s,theta,phi)
      else
        jacobian[i][k] = 1/square(bn[i][k]); // must be boozer jacobian that is proportional to b^-2
      Bmin = mmin(Bmin,bn[i][k]);
      Bmax = mmax(Bmax,bn[i][k]);
    }
  }

//Simpson integration: flux surface averaging
// calculate  <b>; <b^2>;  I[j]=<sqrt(1-x[j]*b)>
  double bavrg=0;  // <b>    -- flux average, where b=B/Bmax(s)
  double b2avrg=0; // <b^2>  -- flux average, where b=B/Bmax(s)
  double norm=0;   // integral( jacobian  dthete dphi)
  for(i=0; i<=ithmax; i++) {
    double wt = (i==0||i==ithmax)?1:(2*(1+i%2));    // weight for Simpson integration
    for(k=0; k<=ipmax; k++) {
      double w = (k==0||k==ipmax)?1:(2*(1+k%2));   // weight for Simpson integration
      w *= wt;
      double b = bn[i][k]/Bmax;
      double jac = jacobian[i][k];
//12.08.2005 for testing integration    b=1; Bmin=b*Bmax; jac=1;
      bavrg  += w*jac*b;            // =w*jac*b
      b2avrg += w*jac*b*b;          // =w*jac*b^2
      norm   += w*jac;
      for(j=0; j<=Nx; j++) {
        double xb=1-x[j]*b;
        xb = sqrt(xb<0?0:xb);
        I[j] += w*jac*xb; //=w*jac*sqrt(1-x[j]*b)
      }
    }
  }
  //////double dA=mNp*dtheta*dphi/9;
  //////std::cerr <<b2avrg*square(Bmax)*dA-square(twopi)<<std::endl; // here must be zero for boozer coordinates

  for(j=0; j<=Nx; j++) {
    I[j] /= norm;    // I[j]=<sqrt(1-x[j]*B/Bmax)> for the given s
    pavrg2.rw()[j]=square(I[j]); 
  }

  bavrg  /= norm;
  b2avrg /= norm;
  Bavr  = bavrg*Bmax;           // <B>
  B2avr = b2avrg*square(Bmax);  // <B^2>

//-------------12.08.2005------------------------
// revised-----31.10.2011------------------------
// I[Nx] == 0 only if b == 1, that is possible only for tokamak magnetic axis
// calculate f_trapped
// calculate integral(from x to 1) dx'/<sqrt(1-x'*b)>, Ipavrg_1[Nx] is the destination array.
//    for testing: integral(0 to 1)   dx/sqrt(1-x) == 2
//    for testing: integral(0 to 1) x*dx/sqrt(1-x) == 4/3
  bool bIs1 = (I[Nx]==0);
  double eps = bIs1?(1e-4):0;  // to avoid singularity for b==1 (tokamak magnetic axis)
  int NS = (ma/mR0<1e-2||bIs1)?256:32;  //NS must be even for simpson integration
  // integrate from c=x[Nx-1] to d=x[Nx]-eps, Simpson integration over x
  double c=x[Nx-1];
  double d=1-eps; 
  double dxx = (d-c)/NS;
  double Intpav=0;
  double IntFc=0;
  for(j=0; j<=NS; j++) {
    double w = (j==0||j==NS)?1:(2*(1+j%2));   // weight for Simpson integration
    double x1= (j==NS)?d:(c+dxx*j);
    double p = sqrt(interp2(x1, mxtr, pavrg2));  // == <sqrt(1-x1*b)>
    Intpav += w/p;
    IntFc += w*x1/p;
  }
  IntFc *= dxx/3;
  Ipavrg_1.rw()[Nx-1] = Intpav*dxx/3;
  Ipavrg_1.rw()[Nx]   = 0;
  if(bIs1) {  // add analytic contribution 
    Ipavrg_1.rw()[Nx-1] += 2*sqrt(eps);            // +integral(1-eps to 1)   dx/sqrt(1-x)
    IntFc               += 2*sqrt(eps)*(1-eps/3);  // +integral(1-eps to 1) x*dx/sqrt(1-x)
  }
#if 1 
  NS=bIs1?128:16;
// These code improves accuracy of integration; it becomes less then 6e-6 for b=1
  Intpav=0;
  double integral=0;
  c=x[Nx-2];d=x[Nx-1];
  dxx = (d-c)/NS;
  for(j=0; j<=NS; j++) {
    double w = (j==0||j==NS)?1:(2*(1+j%2));   // weight for Simpson integration
    double x1 = c+dxx*j;
    double p = sqrt(interp2(x1, mxtr, pavrg2));  // == w/<sqrt(1-x1*b)>
    Intpav += w/p;
    integral += w*x1/p;
  }
  IntFc += integral*dxx/3;
  Ipavrg_1.rw()[Nx-2] = Ipavrg_1.rw()[Nx-1] + Intpav*dxx/3;
// add the rest of points
  dxx = dx/6;
  for(j=Nx-3; j>=0; j--) {
#else
// add the rest of points
  dxx = dx/6;
  for(j=Nx-2; j>=0; j--) {
#endif
    double c = x[j];
    double d = x[j+1];
    double x1= (c+d)/2;
    double p = sqrt(interp2(x1, mxtr, pavrg2));   // == <sqrt(1-x1*b)>
    IntFc += dxx*(c/I[j] + 4*x1/p + d/I[j+1]);    //integral x*dx/I(x) from x[j] to x[j+1], 3-point Simpson
    Intpav = dxx*(1/I[j] + 4/p + 1/I[j+1]);       //integral dx/I(x) from x[j] to x[j+1], 3-point Simpson
    Ipavrg_1.rw()[j] = Ipavrg_1[j+1] + Intpav;
  }
  double fc = 0.75*b2avrg*IntFc;
  if(tokamakConfig&&s<=1e-8) fc = 1; 
  fc = mmin(1.,fc);
  Ftrapped = 1-fc;                   // fraction of trapped particles
  //double tst = Ipavrg_1[0]-2;    // must be 0 at tokamak magnetic axis
}

// for A.Mishchenko
void CStconfig::getFluxAverages1(double s,
                                  double &B_2avr, // <B^-2>
                                  double &R2avr,  // <R^2>
                                  double &R_2avr,  // <R^-2>
                                  double &g2B_2avr, // <grads^2*B^-2>
                                  double &b_12,   // <(1-b)^0.5>
                                  double &b_32,   // <(1-b)^1.5>
                                  double &neoPolarPass,   // for Mishchenko
                                  double &neoPolarTrap,   // for Mishchenko
                                  int ithmax,int ipmax)
{
  if(s<1e-8) s=1e-8;
  int Nx = mxtr.size()-1; //Nx must be even for simpson integration
  double dx=1./Nx;
  const double *x = mxtr.constArray();

  int i,j,k;

  double dtheta, dphi;
  adjustNumbersForSimpson(ithmax,ipmax,dtheta,dphi);

  pBase::ng2Matrix<double> jacobian(ithmax+1,ipmax+1);
  pBase::ng2Matrix<double> bn(ithmax+1,ipmax+1);
  pBase::ng2Matrix<double> Rn(ithmax+1,ipmax+1);
  pBase::ng2Matrix<double> grads2(ithmax+1,ipmax+1);  //  |grads|^2
  double Bmin= 1e15;
  double Bmax=-1e15;
  for(i=0; i<=ithmax; i++) {
    double theta = i*dtheta;
    for(k=0; k<=ipmax; k++) {
      Vector3d magCoord(s,theta,k*dphi);
      NJacobian(magCoord);
      jacobian[i][k] = fabs(getJ()); //  Jacobian in coordinates (s,theta,phi) ~ Ipol/B^2
      bn[i][k] = getmodB();
      Rn[i][k] = getRcyl()[0];
      grads2[i][k] = getGrads().abs2();
      Bmin = mmin(Bmin,bn[i][k]);
      Bmax = mmax(Bmax,bn[i][k]);
    }
  }

  double xm = Bmax/Bmin;
  pBase::ngArray<double> Ips1(mxtr.size());
  pBase::ngArray<double> Ips2(mxtr.size());
  pBase::ngArray<double> mx4(mxtr.size());
  pBase::ngArray<double> Itr(mxtr.size());
  double dxm=(xm-1)/Nx;
  for(j=0; j<=Nx; j++) { mx4[j]=1+j*dxm;Ips1[j]=Ips2[j]=Itr[j]=0;}

//Simpson integration: flux surface averaging
  double r2avrg=0;  // <r^2>  -- flux average, where b=B/Bmax(s)
  double b_2avrg=0; // <b^-2>  -- flux average, where b=B/Bmax(s)
  double r_2avrg=0; // <r^-2>  -- flux average, where b=B/Bmax(s)
  double g2b_2avrg=0; // <grads^2/b^2>  -- flux average, where b=B/Bmax(s)
  double norm=0;   // integral( jacobian  dthete dphi)
  b_12=0;   // <(1-b)^0.5>
  b_32=0;   // <(1-b)^1.5>
  for(i=0; i<=ithmax; i++) {
    double wt = (i==0||i==ithmax)?1:(2*(1+i%2));    // weight for Simpson integration
    for(k=0; k<=ipmax; k++) {
      double w = (k==0||k==ipmax)?1:(2*(1+k%2));   // weight for Simpson integration
      w *= wt;
      double b = bn[i][k]/Bmax;
      double r2 = square(Rn[i][k]);
      double grd2 = grads2[i][k];
//12/08/05for testing integration; b=1; Bmin=b*Bmax;
      double jac = jacobian[i][k];
      b_2avrg+= w*jac/(b*b);     // =w*jac*b^-2
      g2b_2avrg+= w*jac*grd2/(b*b);  // =w*jac*grads^2*b^-2
      r2avrg += w*jac*r2;        // =w*jac*R^2
      r_2avrg+= w*jac/r2;        // =w*jac*R^-2
      norm   += w*jac;
      double b1 = sqrt(1-b);
      b_12   += w*jac*b1;         // <(1-b)^0.5>
      b_32   += w*jac*(1-b)*b1;   // <(1-b)^1.5>
      for(j=0; j<=Nx; j++) {
        double xb=1-x[j]*b;
        xb = sqrt(xb<0?0:xb);
// following 5 lines for A.Mishchenko
        Ips1[j]+= w*jac*xb/b;
        Ips2[j]+= xb==0?(1e20*w*jac):(w*jac*b/xb);
        xb = 1-mx4[j]*b;
        xb = sqrt(xb<0?0:xb);
        Itr[j]+= w*jac*xb/b;
      }
    }
  }
//Simpson integration over x
  double integralPass=0;
  double integralTrap=0;
  for(j=0; j<=Nx; j++) {
    double w = (j==0||j==Nx)?1:(2*(1+j%2));   // weight for Simpson integration
    integralPass += w*(Ips1[j]/norm-norm/Ips2[j]);
    integralTrap += w*Itr[j]/norm;
  }
  neoPolarPass = 1.5*integralPass*dx/3;
  neoPolarTrap = 1.5*integralTrap*dxm/3;

  B_2avr = b_2avrg/norm/square(Bmax);  // <B^-2>
  g2B_2avr = g2b_2avrg/norm/square(Bmax);  // <grads^2*B^-2>

  R2avr  = r2avrg/norm;    // <R^2>
  R_2avr = r_2avrg/norm;    // <R^2>

  b_12   /= norm;         // <(1-b)^0.5>
  b_32   /= norm;         // <(1-b)^1.5>

  for(j=0; j<=Nx; j++) {
    Ips2[j] /= norm; // =<b/sqrt(1-x[j]*b)> for the given s, b=B/Bmax
    mpavrb2.rw()[j]=1/square(Ips2[j]);
  }
}

//****************************************************************************
// Set accuracy of coordinate transformation.
void CStconfig::setAccuracy(double epsA, double epsR) 
{
  if(mepsA!=epsA||mepsR!=epsR) 
    ct().invalidate();
  mepsA=epsA; 
  mepsR=epsR;
}

//****************************************************************************
// Set truncation level to reduce number of harmonics in summation
double CStconfig::truncate(double truncNew)
{
  double truncOld = mmodBtol; 
  if(truncOld==fabs(truncNew)) return truncOld; 
  if(!mOK) return truncOld;

  mmodBtol=fabs(truncNew);
  ct().invalidate();

  if(mmodBtol<=2e-15) { // restore initial value
    for(int is=0; is<=mNs; is++) {
      mpM.rw()[is] = mM;
      mpN.rw()[is] = mN;
    }
    MminCurrent = mM;
    NminCurrent = mN;
    MmaxCurrent = mM;
    NmaxCurrent = mN;
    return truncOld;
  }

  int Mmin = mmax(1,mM);
  int Nmin = mmax(1,mN);
  int Mmax = mmin(1,mM);
  int Nmax = mmin(1,mN);
  for(int is=0; is<=mNs; is++) {
    mpM.rw()[is] = mmin(1,mM);
    mpN.rw()[is] = mmin(1,mN);
    double B0tol = mmodBtol*fabs(bmn(is,0,0));
    double R0tol = mmodBtol*fabs(Rmn(is,0,0));
    for(int m=0; m<=mM; m++)
      for(int n=-mN; n<= mN; n++) {
        double Bw = fabs(bmn(is,m,n));
        double Rw = fabs(Rmn(is,m,n));
        double Zw = fabs(Zmn(is,m,n));
        if(Bw>B0tol || Rw>R0tol || Zw>R0tol) {
          mpM.rw()[is] = mmax(m ,mpM[is]);
          mpN.rw()[is] = mmax(abs(n),mpN[is]);
        }
      }
    Mmin = mmin(Mmin,mpM[is]);
    Nmin = mmin(Nmin,mpN[is]);
    Mmax = mmax(Mmax,mpM[is]);
    Nmax = mmax(Nmax,mpN[is]);
  }
  MminCurrent = Mmin;
  NminCurrent = Nmin;
  MmaxCurrent = Mmax;
  NmaxCurrent = Nmax;
  return truncOld;
}

//****************************************************************************
// Normalize magnetic field
// restore initial value of the magnetic field on axis
void CStconfig::restoreB0()
{
  mBnorm = 1;    // set initial scale factor
}

//****************************************************************************
//Function returns minimum value of the magnetic field |B| on axis
double CStconfig::getB0()
{
  if(!mOK) return 0;
  if(mB0==0) {
    double Bmi=1e20;
    double Bma=-1e20;
    Vector3d magCoord(0.,0.,0.);
    int ip = 200;
    double dphi = mPeriod/ip;
    for(int i=0; i<ip; i++) {
      magCoord[2]=i*dphi;
      double b = B(magCoord);
      Bmi = mmin(Bmi,b);
      Bma = mmax(Bma,b);
    }
    mB0max = Bma;
    return Bmi;
  }
  else
    return fabs(mB0*mBnorm);
}

//****************************************************************************
//Function returns maximum value of the magnetic field |B| on axis
double CStconfig::getB0max()
{
  return fabs(mB0max*mBnorm);
}

//****************************************************************************
//Function returns the magnetic field |B| value on axis at cylindrical angle ficyl
double CStconfig::getB0(double ficyl)
{
  if(!mOK) return 0;
  Vector3d magCoord(0.,0.,0.);
  magCoord[2] = this->phimag(0.,0.,ficyl); // find toroidal angle in magnetic coordinates
  return B(magCoord);  // return |B|
}

//****************************************************************************
//Function returns the value of |B_00| on the magnetic axis
double CStconfig::getB00()
{
  if(!mOK) return 0;
  return fabs(Bmn_sp(0.,0,0));

  double B0=0;
  Vector3d magCoord(0.,0.,0.);
  int ip = 600;
  double dphi = mPeriod/ip;
  for(int i=0; i<ip; i++) {
    magCoord[2]=i*dphi;
    B0 += B(magCoord);
  }
  B0/=ip;
  double B00 = fabs(Bmn_sp(0.,0,0));
  double Bav = Bavrg(0.);
  return B00;
}

//****************************************************************************
//return the sign of B-normalization factor
int CStconfig::getSignOfBnormFactor() const
{
  return mBnorm<0?-1:1;
}

//****************************************************************************
//return the direction of B with respect to the right handed cylindrical coordinate system
int CStconfig::getBdirection()
{
  if(!mOK) return 0;
  double Bfi;
  double  truncSav = truncate();
  truncate(1e-7);
  NJacobianL(Vector3d(1e-2,0.,0.));
  Bfi = getBcyl()[1];
  truncate(truncSav);
  return Bfi<0?-1:1;
}

//****************************************************************************
//set the direction of B with respect to the right handed cyl. coordinate system
void CStconfig::setBdirection(int sign)
{
  if(!mOK) return;
  if(sign==0) return;
  sign = sign<0?-1:1;
  int signcurr  = getBdirection();
  if(sign==signcurr) return;
  mBnorm *= -1;
}

//****************************************************************************
//
void CStconfig::restoreBdirection()
{
  if(!mOK) return;
  mBnorm = fabs(mBnorm);
}



//****************************************************************************
// Normalize magnetic field
// Set minimum magnetic field on axis equals  B0
void CStconfig::setB0(double B0)
{
  if(!mOK) return;
  mBnorm = setSign(B0/mB0,mBnorm);    // set new scale factor
}

//****************************************************************************
// Normalize magnetic field
// Set magnetic field on axis at cylindrical angle 'ficyl' equals to B0
void CStconfig::setB0(double B0, double ficyl)
{
  if(!mOK) return;
  Vector3d magCoord(0.,0.,0.);
  magCoord[2] = this->phimag(0.,0.,ficyl); // find magnetic toroidal angle
  double mBnormSaved = mBnorm;
  mBnorm = 1;   // in order to find the very initial value, see B()
  mBnorm = setSign(B0/B(magCoord),mBnormSaved);
}

//****************************************************************************
// Normalize magnetic field
// Set B00-component of magnetic field on axis equals B0
void CStconfig::setB00(double B0)
{
  if(!mOK) return;
  double mBnormSaved = mBnorm;
  mBnorm = 1;   // in order to find the very initial value, see Bmn_sp()
  double B_00 = fabs(Bmn_sp(0.,0,0));   // get |B_00| on axis
  mBnorm = B0/B_00;                     // set new scale factor
  mBnorm = setSign(B0/B_00,mBnormSaved);
}

//****************************************************************************
// Normalize magnetic field
// Set Ipol at LCMS equals I in Ampers
void CStconfig::setIcoil(double I)
{
  if(!mOK) return;
  double mBnormSaved = mBnorm;
  mBnorm = 1;   // in order to find the very initial value, see Ip()
  double I1 = Ip(1.);   // get Ipol at LCMS
  mBnorm = I/I1;        // set new scale factor
  mBnorm = setSign(mBnorm,mBnormSaved);
}

// get Ipol at LCMS equals I in Ampers
double CStconfig::getIcoil()
{
  if(!mOK) return 0;
  return fabs(Ip(1.));   // get Ipol at LCMS
}

double CStconfig::getBscaleFactor()
{
  return fabs(mBnorm);
}

//****************************************************************************
// Scale magnetic field
void CStconfig::scaleB(double multiplier)
{
  if(!mOK) return;
  mBnorm = setSign(multiplier,mBnorm); // set new scale factor keeping direction
}

#if 0  // arithmetic average
//****************************************************************************
// Set average toroidal magnetic field B_phi(cylindrical component) on magnetic axis.
void CStconfig::setAverageBphi0(double B0)
{
  if(!mOK) return;
  mBnorm = 1;   // in order to find the very initial value
  double Bphi0 = getAverageBphi0();   // find average Bphi on axis
  mBnorm = B0/fabs(Bphi0);            // set new scale factor
}

//****************************************************************************
//Function returns the average toroidal magnetic field B_phi(cylindrical component) on magnetic axis.
double CStconfig::getAverageBphi0()
{
  if(!mOK) return 0;

  volatile ctStack::saveState state(ct);

  double Bfi=0;
  Vector3d magCoord(0.,0.,0.);
  int ip = 200;
  double dphi = mPeriod/ip;
  for(int i=0; i<ip; i++) {
    magCoord[2]=i*dphi;
    NJacobian(magCoord);
    Bfi += getBcyl()[1];
  }
  return Bfi/ip;
}
#endif

//****************************************************************************
// Get magnetic coordinates of a field line at angle phi.
// This line starts from the point magCoord0=(s0,theta0,phi0)
//
// The method returns magCoord=(s0,theta,phi) calculating the
//  theta(phi)=theta0+iota(s0)*(phi-phi0) - [lambda(s0,theta,phi)-lambda(theta0,phi0)],
//
//  [theta(phi)+lambda(s0,theta,phi)]   -   [theta0+lambda(theta0,phi0) + iota(s0)*(phi-phi0) ] = 0
// where lambda is non zero for VMEC coordinates
//
Vector3d CStconfig::magCoordFieldLine(const Vector3d &magCoord0, double phi, double iotA)
{ 

  double     s0 = magCoord0[0];
  double theta0 = magCoord0[1];
  double   phi0 = magCoord0[2];
  iotA = iotA==0?iota(s0):iotA;    // iota(s) does not change ct
  double theta = theta0+iotA*(phi-phi0);
 
  Vector3d magCoord(s0,theta,phi);
  if( !(vmecCoordinates||fpCoordinates) ) return magCoord;

  volatile ctStack::saveState state(ct);   // save ct 

#if 1
  calcVmecLgB(magCoord0);
  theta += ct().lambda; // make theta1 = theta0+ct().lambda(theta0,phi0) + iotA*(phi-phi0);
  double theta1  = theta;
  magCoord[1] = theta;
  Vector3d vmec(magCoord);
  // find theta for VMEC coordinates as root of the expression
  //   theta+lambda(s0,theta,phi) - theta1 = 0
  bool ok=false;
  int iter = 1;
  double eps=1e-10;
  vmec[1] = theta;
  calcVmecLgB(vmec);
  double f = (theta+ct().lambda - theta1);  // residual
  for(int i=0; i<300; ++i,++iter) {
    double dtheta = -f/(1+ct().lambdat);      // correction to root
    // backtracking loop: 
    double thetaOld = theta;
    double fOldabs  = fabs(f);
    double step = 1;
    int N=5;
    while(--N) {  // backtracking loop:
      theta = thetaOld + step*dtheta; // try the Newton step
      { //this block calculates f = theta^* - theta1 (or f = theta+lambda(s0,theta,phi) - theta1)
        vmec[1] = theta;
        calcVmecLgB(vmec);
        f = (theta+ct().lambda - theta1);  // f is the new residual
      }
      if(fabs(f)<eps) return vmec;  
      if(fabs(f)>fOldabs) step *= 0.5; // shorten the step size if the residual does not decrease
      else break; // step is OK
    }
  }
std::cerr << "CStconfig::magCoordFieldLine(): # of iterations " <<iter<<std::endl;
return vmec;
#else
  calcVmecLgB(magCoord0);
  theta += ct().lambda;
  double theta1  = theta;
  magCoord[1] = theta;
  // find theta for VMEC coordinates
  // find theta as root of the expression
  //   theta+lambda(s0,theta,phi) - theta1 = 0
  bool ok=false;
  int iter = 0;
  double eps=1e-10;
  for(int i=0; i<300; ++i,++iter) {
    magCoord[1] = theta;
    calcVmecLgB(magCoord);
    double dtheta = -(theta+ct().lambda - theta1)/(1+ct().lambdat); // Newton's iterations
    if(i>5) { 
      dtheta *=0.5; //try to backtrack if convergence is too slow
      //std::cerr <<iter<<std::endl<<std::flush;
    }
    theta += dtheta;
    if(fabs(dtheta)<eps) {ok=true; break;}
  }
  if(!ok) {
    std::cerr << "CStconfig::magCoordFieldLine(): # of iterations " <<iter<<std::endl;
  }
  return magCoord;
#endif 
}


//****************************************************************************
// TODO: does not work properly
// Get magnetic coordinates of a field line at angle theta.
// This line starts from the point magCoord0=(s0,theta0,phi0)
//
// The method returns magCoord=(s0,theta,phi) calculating the
//  phi(theta)=...............
// where lambda is non zero for VMEC coordinates
//
Vector3d CStconfig::magCoordFieldLine2(const Vector3d &magCoord0, double theta, double iotA)
{ 
  volatile ctStack::saveState state(ct);

  double     s0 = magCoord0[0];
  double theta0 = magCoord0[1];
  double   phi0 = magCoord0[2];
  iotA = iotA==0?iota(s0):iotA;
  double phi1 = phi0+(theta-theta0)/iotA;
//  phi1 = modPeriod(phi1);
  
  Vector3d magCoord(s0,theta,phi1);
return magCoord;

  if( vmecCoordinates||fpCoordinates ) {
    calcVmecLgB(magCoord0);
    double phi10  = phi1 - ct().lambda/iotA;
    // find phi for VMEC coordinates
    // find phi as root of the expression
    //   phi-lambda(s0,theta,phi)/iota - phi10 = 0
    //////calcVmecLgB(magCoord);
    //////double phi=phi10+ct().lambda/iotA;
    double phi=phi1;
    bool ok=false;
    int iter = 0;
    double eps=1e-6;
    for(int i=0; i<300; ++i,++iter) {
      magCoord[2] = phi;
      calcVmecLgB(magCoord);
      double dphi = -(phi-ct().lambda/iotA - phi10)/(1-ct().lambdap/iotA); // Newton's iterations
      if(i>50) { 
        dphi *=0.5; //try to backtrack if convergence is too slow
        //std::cerr <<iter<<std::endl<<std::flush;
      }
      phi += dphi;
      if(fabs(dphi)<eps) {ok=true; break;}
    }
    if(!ok) {
      std::cerr << "CStconfig::magCoordFieldLine2(): # of iterations " <<iter<<std::endl;
    }
  }
 
  return magCoord;
}

//****************************************************************************
// Get B on field line, that passes through the point magCoord0=(s0,theta0,phi0)
// method  calculates function
//  theta(phi)=theta0+iota(s0)*(phi-phi0) - [lambda(s0,theta,phi)-lambda(theta0,phi0)]
// where lambda is non zero only for VMEC coordinates
double CStconfig::BonFieldLine(const Vector3d &magCoord0, double phi)
{
  return B( magCoordFieldLine(magCoord0, phi) );
}


//****************************************************************************
// advance mixed coordinates and return magnetic coordinates  
Vector3d CStconfig::mixCoordFieldLine(const Vector3d &mixCoord0, double fcyl)
{
  if(!mOK) return Vector3d(0.,0.,0.);
  
  if( tokamakEQDSK||vmecCoordinates||fpCoordinates ) {
    return magCoordFieldLine(mixCoord0,fcyl);
  }
 // otherwise this is Boozer coordinates 
  Vector3d magCoord0(mixCoord0);
  magCoord0[2] = phimag(mixCoord0); 
  if(mixCoord0[0]==0) magCoord0[1]=0;
  if(mixCoord0[2]==fcyl) return magCoord0;

#if 1
  // find magnetic toroidal angle phi for the given cylindrical angle fi
  //     fi(is,theta,phi) - f = 0
  double eps = 1e-9;
  bool onFildLine=true;
  double dfi_dphi;
  double phi  = fcyl;  // assign guess
  Vector3d magCoord = magCoordFieldLine(magCoord0,phi);
  double f1 = fi(magCoord,dfi_dphi,onFildLine);
  double dF = fcyl-f1;   // residual
  for(int i=1,iter=400; i<=iter; ++i) {
    double dphi = dF/dfi_dphi;  // correction to root
    // backtracking loop: 
    double phiOld = phi;
    double dFold  = fabs(dF);
    double step = 1;
    int N=5;
    while(--N) {  // backtracking loop:
      phi = phiOld + step*dphi; // try the Newton step
      magCoord = magCoordFieldLine(magCoord0,phi);
      f1 = fi(magCoord,dfi_dphi,onFildLine);
      dF = fcyl-f1;   // new residual
      if(fabs(dF)<eps) break;
      if(fabs(dF)>dFold) step *= 0.5; // shorten the step size if the residual does not decrease
      else break; // step is OK
    }
    if(fabs(dF)<eps) break;   
    if(i<iter) continue;
    std::cerr << "CStconfig::mixCoordFieldLine(): # of iterations " <<i<<std::endl;
    break;
  }
  return magCoord;  
#else
  double eps = 1e-9;
  double phi  = fcyl;  // assign guess
  bool ok=false;
  bool onFildLine=true;
  Vector3d magCoord;
  int iter = 0;
  for(int i=0; i<400; ++i,++iter) {
    double dfi_dphi;
    magCoord = magCoordFieldLine(magCoord0,phi);
    double f1 = fi(magCoord,dfi_dphi,onFildLine);
    if(fabs(fcyl-f1)<eps) {ok=true; break;}
    double dphi = (fcyl-f1)/dfi_dphi;  // Newton's iterations
    phi += dphi;
    //if(fabs(dphi)<eps) {ok=true; break;}
  }
  if(!ok) {
    std::cerr << "CStconfig::mixCoordFieldLine(): # of iterations " <<iter<<std::endl;
  }
  //magCoord = magCoordFieldLine(magCoord0,ph);
  return magCoord;
#endif
}


double CStconfig::r   (double s) const {return ma*sqrt(s);};      // effective minor radius[m], VMEC definition
double CStconfig::rdkes(double s) const {return sqrt(fabs(Flux(s)/(Bmn_sp(s,0,0)*pi)));};  // minor effective radius[m], DKES definition
double CStconfig::Flux(double s) const {return mBnorm*mFlux*s; };       // toroidal flux
//double CStconfig::iota(double s) const {return SPLINEy__(s,miota);};    // iota
double CStconfig::iota(double s) const {return miota(s,this);};    // iota
double CStconfig::pressure(double s) const {return vmecCoordinates?mpres(s,this):0;};  
double CStconfig::iotaPrime(double s) const {return miota.prime(s,this);};    // diota/ds

double CStconfig::Vprime(double s) const {return g00(s)*mNp;}           // dV/ds

double CStconfig::pp  (double s) 
{
  ////const_cast<CStconfig*>(this)->calculateCurrents();
  ////double Ipp = mBnorm*mNp*mIpol.prime(s);    //Ipol prime
  ////double Itp = mBnorm*mItor.prime(s);        //Ipol prime
  ////if(fabs(It(1.)/Ip(1.))<1e-7) Itp=0; 
  ////return -(Ipp+iota(s)*Itp)*Flux(1.)/Vprime(s);

  if(!mpprime.created) {
    mpprime.create(this);
  }
  return mBnorm*mBnorm*mpprime(s,this);
}

double CStconfig::sPol(double s) // normalized poloidal flux
{
  if(!msPol.created) { 
// calculate poloidal flux using iota; create s_poloidal-grid (normalized poloidal flux)
    double mBnormSave = mBnorm;
    mBnorm = 1;
    msPol.rw()[0]=0;
    for(int is=1; is<=mNs; is++) msPol.rw()[is]=msPol[is-1]+FluxPoloidal(ms[is-1],ms[is]);
    msPol.create(this);
    double fluxPolMax = sPol(1.);
    mPolFlux = fluxPolMax*mFlux;
    for(int i=0; i<=mNs; i++) msPol.rw()[i] /= fluxPolMax;
    msPol.create(this);
    mBnorm = mBnormSave;
  }
  return msPol(s,this);
}

double CStconfig::PoloidalFlux(double s) // poloidal flux
{
  double spol = sPol(s);
  return mBnorm*mPolFlux*spol;
}

double CStconfig::sToroidal(double sPoloidal)
{
  double spol = sPol(1.);  // create arrays and splines
  return ms.interp2(sPoloidal,msPol);
}


double CStconfig::Ip  (double s) const // poloidal current
{
  const_cast<CStconfig*>(this)->calculateCurrents();
  return mBnorm*mNp*mIpol(s,this);
} 

double CStconfig::It  (double s) const  // toroidal current
{
  const_cast<CStconfig*>(this)->calculateCurrents();
  return mBnorm*mItor(s,this);
}

double CStconfig::g00 (double s) const 
{
  const_cast<CStconfig*>(this)->calculateVprime();
  return mg00(s,this);
}

double CStconfig::Volume(double s) 
{
  this->calculateVprime();
  if(!mVolume.created) { 
    mVolume.rw()[0]=0;
    for(int is=1; is<=mNs; is++) mVolume.rw()[is]=mVolume[is-1]+VolumeInt(ms[is-1],ms[is]);
    mVolume.create(this);
  }
  return mVolume(s,this);
}



double CStconfig::S11(double s) // S11(s)  in (Flux,theta,phi)-coordinates
{
  if(!SMatrixReady) calculateSMatrix(); 
  return mBnorm*mS11.interp2(s,msK);
};
double CStconfig::S12(double s) // S12(s)
{
  if(!SMatrixReady) calculateSMatrix(); 
  return mBnorm*mS12.interp2(s, msK);
};
double CStconfig::S21(double s) // S21(s)
{
  if(!SMatrixReady) calculateSMatrix();
  return mBnorm*mS21.interp2(s, msK);
};
double CStconfig::S22(double s) // S22(s)
{
  if(!SMatrixReady) calculateSMatrix();
  return mBnorm*mS22.interp2(s, msK);
};

double CStconfig::Bmin(double s) // Bmin(s)
{
  if(!BminmaxReady) calculateBminmax();
  return fabs(mBnorm)* interp2(sqrt(s), mreK, mBmin);
};
double CStconfig::Bmax(double s) // Bmax(s)
{
  if(!BminmaxReady) calculateBminmax();
  return fabs(mBnorm)* interp2(sqrt(s), mreK, mBmax);
};

double CStconfig::Bmin() // Bmin()
{
  if(!BminmaxReady) calculateBminmax();
  return fabs(mBnorm)* mBminGlobal;
};
double CStconfig::Bmax() // Bmax()
{
  if(!BminmaxReady) calculateBminmax();
  return fabs(mBnorm)* mBmaxGlobal;
};

double CStconfig::Bavrg(double s) // <B(s)>
{
  if(!ftrappedReady) calculateFtrappedBminmax();
  return fabs(mBnorm)* interp2(sqrt(s), mreK, mBavrg);
};

double CStconfig::B2avrg(double s) // <B^2(s)>
{
  if(!ftrappedReady) calculateFtrappedBminmax();
  return mBnorm*mBnorm*interp2(sqrt(s), mreK, mB2avrg);
};

double CStconfig::B2VolAvrg()  // volume average of B^2
{
  if(!ftrappedReady) calculateFtrappedBminmax();
  return mBnorm*mBnorm*mB2VolAvrg;
}

double CStconfig::ftrapped(double s) // fraction of trapped particles
{
  if(!ftrappedReady) calculateFtrappedBminmax();
  return sqrt(interp2(sqrt(s), mreK, mFtrap));
};

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
double CStconfig::B_2avrg(double s) // <B^-2(s)>
{
  if(!Averages1Ready) calculateFluxAverages1();
  return interp2(sqrt(s), mreK, mB_2avrg)/(mBnorm*mBnorm);
};

double CStconfig::grads2B_2avrg(double s) // <grads^2*B^-2(s)>
{
  if(!Averages1Ready) calculateFluxAverages1();
  return interp2(sqrt(s), mreK, mGrads2B_2avrg)/(mBnorm*mBnorm);
};

double CStconfig::R2avrg(double s) // <R^2(s)>
{
  if(!Averages1Ready) calculateFluxAverages1();
  return interp2(sqrt(s), mreK, mR2avrg);
};

double CStconfig::R_2avrg(double s) // <R^-2(s)>
{
  if(!Averages1Ready) calculateFluxAverages1();
  return interp2(sqrt(s), mreK, mR_2avrg);
};

double CStconfig::NeoPolarizationPass(double s)
{
  if(!Averages1Ready) calculateFluxAverages1();
  return interp2(sqrt(s), mreK, mNeoPolarztnPass)/square(Bmax(s));
};

double CStconfig::NeoPolarizationTrap(double s)
{
  if(!Averages1Ready) calculateFluxAverages1();
  return interp2(sqrt(s), mreK, mNeoPolarztnTrap)/square(Bmax(s));
};

double CStconfig::NeoPolarization(double s)
{
  return NeoPolarizationPass(s)+NeoPolarizationTrap(s);
}

double CStconfig::db12Avrg(double s) // <(1-b)^0.5>
{
  if(!Averages1Ready) calculateFluxAverages1();
  return interp2(sqrt(s), mreK, mb12);
}

double CStconfig::db32Avrg(double s) // <(1-b)^1.5>
{
  if(!Averages1Ready) calculateFluxAverages1();
  return interp2(sqrt(s), mreK, mb32);
}

//****************************************************************************
// <ro_paralell^-1>
// returns f(s,x) = <b/sqrt(1-x*b)>, where 0<=x<=1, 0<=s<=1, and b=B(s,th,ph)/Bmax(s)
// linear interpolation is used
double CStconfig::roParInvAvrg(double s, double x) const
{
  if(!Averages1Ready) 
     const_cast<CStconfig*>(this)->calculateFluxAverages1();

  int ir = mreK.bsearch(sqrt(s));
  double dr = mreK[ir+1]-mreK[ir];

  double dx = mxtr[1];
  double px = x/dx;
  int jx = int(px);
  jx = mmin(jx,int(mxtr.size()-2));
  jx = mmax(jx,0);
  px -= jx;
  double f1=mpavrb_sx[ir]  [jx];
  double f2=mpavrb_sx[ir+1][jx];
  f1 += (mpavrb_sx[ir]  [jx+1]-f1)*px;
  f2 += (mpavrb_sx[ir+1][jx+1]-f2)*px;
  f1 += (f2-f1)*(sqrt(s)-mreK[ir])/dr;
  return 1./sqrt(f1);
};
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************


//****************************************************************************
// f(s,x) = <sqrt(1-x*b)>, where 0<=x<=1, 0<=s<=1, and b=B(s,th,ph)/Bmax(s)
// linear interpolation is used
double CStconfig::pavrg(double s, double x) const
{
  if(!ftrappedReady) 
    const_cast<CStconfig*>(this)->calculateFtrappedBminmax();

  int ir = mreK.bsearch(sqrt(s));
  double dr = mreK[ir+1]-mreK[ir];

  double dx = mxtr[1];
  double px = x/dx;
  int jx = int(px);
  jx = mmin(jx,int(mxtr.size()-2));
  jx = mmax(jx,0);
  px -= jx;
  double f1=mpavr_sx[ir]  [jx];
  double f2=mpavr_sx[ir+1][jx];
  f1 += (mpavr_sx[ir]  [jx+1]-f1)*px;
  f2 += (mpavr_sx[ir+1][jx+1]-f2)*px;
  f1 += (f2-f1)*(sqrt(s)-mreK[ir])/dr;
  return sqrt(f1);
};

//****************************************************************************
// @return Integral(from x to 1) dx'/f(s,x') ,
//         where f(s,x) = <sqrt(1-x*B(s,th,ph)/Bmax(s))> 0<=x<=1, 0<=s<=1
//         linear interpolation is used
double CStconfig::integralInvPavrg(double s, double x) const
{
  if(!ftrappedReady) 
    const_cast<CStconfig*>(this)->calculateFtrappedBminmax();

  int ir = mreK.bsearch(sqrt(s));
  double dr = mreK[ir+1]-mreK[ir];

  double dx = mxtr[1];
  double px = x/dx;
  int jx = int(px);
  jx = mmin(jx,int(mxtr.size()-2));
  jx = mmax(jx,0);
  px -= jx;
  double f1=mpavrI_sx[ir]  [jx];
  double f2=mpavrI_sx[ir+1][jx];
  f1 += (mpavrI_sx[ir]  [jx+1]-f1)*px;
  f2 += (mpavrI_sx[ir+1][jx+1]-f2)*px;
  f1 += (f2-f1)*(sqrt(s)-mreK[ir])/dr;
  return sqrt(f1);
};

//****************************************************************************
// function returns  B(s,theta,phi)  = sum_mn B_mn(s)*cos(m*theta-n*Np*phi)
// where  (s,theta,phi) is magnetic coordinates.
//
double CStconfig::B(const Vector3d &u) const
{
  if(!mOK) return 0;
  double modB=0;
  double s     = u[0];
  double theta = u[1];
  double phi   = u[2];

  int M,N;
  int is = getSsegMN(s, M, N); // search nearest segment in array ms;
  if(s==0) M=0;
  const double * fC = mCOEF.constArray()+mN*offsetN+ halfMesh.is0*offsetS;
  const splineCoeffs * hMesh = &halfMesh;

  double csNp = tokamakEQDSK?0:(::cos(mNp*phi)); //if(tokamakEQDSK) phi=pi/2;
  double snNp = tokamakEQDSK?1:(::sin(mNp*phi));

  CComplex expNph(csNp,snNp);               // exp(1i*mNp*phi), where 1i is the sqrt(-1)
  CComplex exp_Nph = expNph.conj();         // exp(-1i*mNp*phi)
  CComplex expth (::cos(theta), ::sin(theta));  // exp(1i*theta)
  CComplex expmth(1,0);                     // will be exp(1i*m*theta)
  for(int m=0; m<=M; m++) {
    CComplex expmthNph =expmth; // exp(1i*m*theta+1i*mNp*phi), where 1i is the sqrt(-1)
    CComplex expmth_Nph=expmth; // exp(1i*m*theta-1i*mNp*phi), where 1i is the sqrt(-1)
    const double * fCoeff = fC + m*offsetM;
    for(int n=0; n<= N; n++) {
      double iBmn, iBm_n;       // intrepolated value of B(m,n) and B(m,-n) in point s
      int idx  =  n*offsetN;
      SPLINE2(bmn,idx,iBmn,iBm_n,hMesh)
      if(n==0) {
        modB += iBmn*expmth.re;  // *cos(m*theta)
        continue;
      }
      expmth_Nph *= exp_Nph; // exp(1i*m*theta-1i*mNp*phi)
      expmthNph  *= expNph;  // exp(1i*m*theta+1i*mNp*phi)
      modB += iBmn*expmth_Nph.re + iBm_n*expmthNph.re; //iBmn*cos(m*theta-n*Np*phi)+iBm_n*cos(m*theta+n*Np*phi)
    }
    expmth *= expth;   // exp(1i*m*theta)
  }
  return fabs(modB*mBnorm);    // mBnorm  can be negativ
}

//****************************************************************************
// Bcontravariant is optional
// function returns  Bcontr(s,theta,phi)  = sum_mn B_mn^theta,phi(s)*cos(m*theta-n*Np*phi)
//    where  (s,theta,fi) is the VMEC coordinates.
// only for vmecCoordinates or fpCoordinates 
void CStconfig::calcVmecLgB(const Vector3d &u, Vector3d * Bcontrav, bool boozTrans)
{
  Vector3d Bcontra(0.,0.,0.);
  if(boozTrans) {
    ct().p       = 0;  
    ct().ps      = 0;  // d(p)/ds
    ct().pt      = 0;  
    ct().pp      = 0; 
  }
  ct().jacVmec = 0;  // Vmec Jacobian
  ct().lambda  = 0;  // Vmec lambda
  ct().lambdas = 0;  // d(lambda)/ds
  ct().lambdat = 0;  // d(lambda)/dtheta
  ct().lambdap = 0;  // d(lambda)/dphi

  if(!(vmecCoordinates||fpCoordinates)) {
    if(Bcontrav) 
      *Bcontrav = 0;
    return;
  }

#define LAMBDAS_NEEDED 1

  double s     = u[0];
  double theta = u[1];
  double phi   = u[2];

  double p(0);     
  double ps(0);     
  double pt(0);     
  double pp(0);     
  double g(0);     // Vmec Jacobian
  double lam(0);   // lambda
  double lams(0);  // d(lambda)/ds
  double lamt(0);  // d(lambda)/dtheta
  double lamp(0);  // d(lambda)/dphi

  int M,N;
  int is = getSsegMN(s, M, N); // search nearest segment in array ms;
  if(s==0) M=0;
  const double * fC = mCOEF.constArray()+mN*offsetN+ halfMesh.is0*offsetS;
  const splineCoeffs * hMesh = &halfMesh;

  CComplex expNph(::cos(mNp*phi),::sin(mNp*phi)); // exp(1i*mNp*phi), where 1i is the sqrt(-1)
  CComplex exp_Nph = expNph.conj();               // exp(-1i*mNp*phi)
  CComplex expth (::cos(theta), ::sin(theta));    // exp(1i*theta)
  CComplex expmth(1,0);                           // will be exp(1i*m*theta)
  for(int m=0; m<=M; m++) {
    CComplex expmthNph =expmth; // exp(1i*m*theta+1i*mNp*phi), where 1i is the sqrt(-1)
    CComplex expmth_Nph=expmth; // exp(1i*m*theta-1i*mNp*phi), where 1i is the sqrt(-1)
    const double * fCoeff = fC + m*offsetM;
    double iBthmn(0), iBthm_n(0);       // intrepolated value of BcontraTheta_mn in point s
    double iBfimn(0), iBfim_n(0);       // intrepolated value of BcontraFi_mn    in point s
    double igmn(0), igm_n(0);           // intrepolated value of gmn    in point s
    double ilmn(0), ilm_n(0);           // intrepolated value of lmn    in point s
    double ipmn(0), ipm_n(0);           // intrepolated value of pmn    in point s
    for(int n=0; n<= N; n++) {
      int idx  =  n*offsetN;
      double ipmns(0),ipm_ns(0);   // intrepolated value of d(pmn)/ds    in point s
      if(!fpCoordinates) {
        if(Bcontrav) {
          SPLINE(bContraTheta_mn,  idx,iBthmn ,hMesh)
          SPLINE(bContraTheta_mn, -idx ,iBthm_n,hMesh)
          SPLINE(bContraFi_mn,     idx,iBfimn ,hMesh)
          SPLINE(bContraFi_mn,    -idx ,iBfim_n,hMesh)
        }
        SPLINE(gmn,              idx,igmn,   hMesh)
        SPLINE(gmn,             -idx,igm_n,  hMesh)
        if(boozTrans) {
          SPLINED(Phimn,            idx,ipmn, ipmns,  hMesh)
          SPLINED(Phimn,           -idx,ipm_n,ipm_ns, hMesh)
        }
      }
#if LAMBDAS_NEEDED
      double ilmns,ilm_ns;   // intrepolated value of d(lmn)/ds    in point s
      SPLINED(lmn,             idx,ilmn, ilmns, hMesh)
      SPLINED(lmn,            -idx,ilm_n,ilm_ns,hMesh)
#else
      SPLINE(lmn,              idx,ilmn,   hMesh)
      SPLINE(lmn,             -idx ,ilm_n,  hMesh)
#endif
      if(n==0) {
        double csmt = expmth.re;  // cos(m*theta)
        double snmt = expmth.im;  // sin(m*theta)
        if(!fpCoordinates) {
          if(Bcontrav) {
            Bcontra[1] += iBthmn*csmt;
            Bcontra[2] += iBfimn*csmt;
          }
          g    += igmn*csmt;
          if(boozTrans) {
            ps   += ipmns*snmt;    //d(p)/ds
            p    += ipmn*snmt;     //p
            pt   += m*ipmn*csmt;   //d(p)/dtheta
          }
        }
#if LAMBDAS_NEEDED
        lams    += ilmns*snmt;      //d(lambda)/ds
#endif
        lam     +=   ilmn*snmt;     //lambda
        lamt    += m*ilmn*csmt;     //d(lambda)/dtheta
        continue;
      }
      expmth_Nph *= exp_Nph;    // exp(1i*m*theta-1i*mNp*phi)
      expmthNph  *= expNph;     // exp(1i*m*theta+1i*mNp*phi)
      double sn1 = expmth_Nph.im;  //sin(m*theta-n*Np*phi) is stored into sn1
      double cs1 = expmth_Nph.re;  //cos(m*theta-n*Np*phi)
      double sn2 = expmthNph.im;   //sin(m*theta+n*Np*phi)
      double cs2 = expmthNph.re;   //cos(m*theta+n*Np*phi)
      if(!fpCoordinates) {
        if(Bcontrav) {
          Bcontra[1] += iBthmn*cs1 + iBthm_n*cs2;
          Bcontra[2] += iBfimn*cs1 + iBfim_n*cs2;
        }
        g   +=    igmn*cs1  + igm_n*cs2;
        if(boozTrans) {
          p   +=    ipmn*sn1  + ipm_n*sn2;      //p
          ps  +=   ipmns*sn1  + ipm_ns*sn2;     //ps
          pt  += m*(ipmn*cs1  + ipm_n*cs2);     //d(p)/dtheta - poloidal derivative
          pp  -= n*(ipmn*cs1  - ipm_n*cs2);     //d(p)/dphi
        }
      }
#if LAMBDAS_NEEDED
      lams    += ilmns*sn1  + ilm_ns*sn2;     //lambdas
#endif
      lam   +=    ilmn*sn1  + ilm_n*sn2;      //lambda
      lamt  += m*(ilmn*cs1  + ilm_n*cs2);     //d(lambda)/dtheta - poloidal derivative
      lamp  -= n*(ilmn*cs1  - ilm_n*cs2);     //d(lambda)/dphi
    }
    expmth *= expth;   // exp(1i*m*theta)
  }

#if LAMBDAS_NEEDED
  #if INTERPOLATION_ORDER == 3
    lams *= halfMesh.e2rcurr;
  #endif
  ct().lambdas = lams;      // d(lambda)/ds
#endif
  if(boozTrans) {
    #if INTERPOLATION_ORDER == 3
    ps *= halfMesh.e2rcurr;
    #endif
    ct().ps      = ps; 
    ct().p       = p;  
    ct().pt      = pt;  
    ct().pp      = pp*mNp; 
  }
  ct().lambda  = lam;       // Vmec lambda
  ct().lambdat = lamt;      // d(lambda)/dtheta
  ct().lambdap = lamp*mNp;  // d(lambda)/dphi
  ct().jacVmec = g;         // Vmec Jacobian   in (s,theta,phi)
  ct().jac     = g/(mFlux*mBnorm); // Vmec Jacobian   in (flux,theta,phi)
  if(Bcontrav) 
    *Bcontrav = Bcontra*mBnorm;    // mBnorm  can be negative
#undef LAMBDAS_NEEDED
}

//****************************************************************************
//
//
Vector3d CStconfig::SFLcoordOnBlineAtTheta(const Vector3d &pest0, double theta, double iotA)
{ 

  double     s0 = pest0[0];
  double theta0 = pest0[1];
  double   phi0 = pest0[2];
  iotA = iotA==0?iota(s0):iotA;
  double phi1 = phi0+(theta-theta0)/iotA;
  
  return Vector3d(s0,theta,phi1);
}

//****************************************************************************
// VMEC coordinates vmecCoord=(s,theta,phi)
// @return PEST coordinates
Vector3d CStconfig::Vmec2Pest(const Vector3d &vmecCoord)
{
  //bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);
  //if(boozerCrd) return vmecCoord;
  if( !(vmecCoordinates||fpCoordinates) ) return vmecCoord;

  volatile CStconfig::ctStack::saveState state(ct);
  calcVmecLgB(vmecCoord,0,true);
  Vector3d pest = vmecCoord;
  pest[1] += ct().lambda; // theta + lambda
  return pest;
}

//****************************************************************************
// PEST coordinates pest=(s,theta^star,phi)
// @return VMEC coordinates
Vector3d CStconfig::Pest2Vmec(const Vector3d &pest)
{ 
  if( !(vmecCoordinates||fpCoordinates) ) return pest;

  volatile ctStack::saveState state(ct);

  double theta1 = pest[1];
  double theta  = pest[1];  // set guess to theta
  
  Vector3d vmec(pest);

  // find theta for VMEC coordinates
  // find theta as root of the expression
  //   theta+lambda(s0,theta,phi) - theta1 = 0
  bool ok=false;
  int iter = 1;
  double eps=1e-10;
  vmec[1] = theta;
  calcVmecLgB(vmec);
  double f = (theta+ct().lambda - theta1);  // residual
  for(int i=0; i<300; ++i,++iter) {
    double dtheta = -f/(1+ct().lambdat);      // correction to root
    // backtracking loop: 
    double thetaOld = theta;
    double fOldabs  = fabs(f);
    double step = 1;
    int N=5;
    while(--N) {  // backtracking loop:
      theta = thetaOld + step*dtheta; // try the Newton step
      { //this block calculates f = theta^* - theta1 (or f = theta+lambda(s0,theta,phi) - theta1)
        vmec[1] = theta;
        calcVmecLgB(vmec);
        f = (theta+ct().lambda - theta1);  // f is the new residual
      }
      if(fabs(f)<eps) return vmec;  
      if(fabs(f)>fOldabs) step *= 0.5; // shorten the step size if the residual does not decrease
      else break; // step is OK
    }
  }
  std::cerr << "CStconfig::Pest2Vmec(): # of iterations " <<iter<<std::endl; 
  return vmec;
}

//****************************************************************************
void CStconfig::Vmec2PestContraBasis(const Vector3d &vmec, Vector3d &e1, Vector3d &e2, Vector3d &e3, Vector3d &pest)
{
  volatile CStconfig::ctStack::saveState state(ct);
  NJacobian(vmec);
  e1 = ct().gradS;
  e2 = ct().gradTheta+getGradlambda();
  e3 = ct().gradPhi;
  pest = vmec;
  pest[1] += ct().lambda; // theta + lambda
}

//****************************************************************************
void CStconfig::SFLcontraBasis(const Vector3d &pest, Vector3d &e1, Vector3d &e2, Vector3d &e3, Vector3d &vmec)
{
  volatile CStconfig::ctStack::saveState state(ct);
  vmec = Pest2Vmec(pest);
  NJacobian(vmec);
  e1 = ct().gradS;
  e2 = ct().gradTheta+getGradlambda(); // grad(theta + lambda)
  e3 = ct().gradPhi;
  return;
  Vector3d pest1 = vmec;
  pest1[1] += ct().lambda; // theta + lambda
}

//****************************************************************************
// booz is the boozer coordinates
// e1,e2,e3 are boozer contravariant basis 
// mag==booz if input file is the boozer file
// mag is the corresponding VMEC coordinates  if the input file is the VMEC wout file
void CStconfig::BoozerContraBasis(const Vector3d &booz, Vector3d &e1, Vector3d &e2, Vector3d &e3, Vector3d &mag)
{
  volatile CStconfig::ctStack::saveState state(ct);
  mag = Boozer2Vmec(booz);
  NJacobian(mag);
  if(!vmecCoordinates) {
    getContraBasis(e1,e2,e3);
    return;
  }
  calcVmecLgB(mag,0,true);
// get gradient of the angle-transformation function p 
// used for transformation from VMEC to boozer coordinates
  Vector3d grp = ct().ps*ct().gradS + ct().pt*ct().gradTheta + ct().pp*ct().gradPhi;
  double p = ct().p;
  double s = booz[0];
  e1 = ct().gradS;
  e2 = ct().gradTheta+getGradlambda()+iota(s)*grp+(iotaPrime(s)*p)*e1;
  e3 = ct().gradPhi + grp;
}

//****************************************************************************
// VMEC coordinates vmecCoord=(s,theta,phi)
// @return Boozer coordinates
Vector3d CStconfig::Vmec2Boozer(const Vector3d &vmecCoord)
{
  //bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);
  //if(boozerCrd) return vmecCoord;
  if(!vmecCoordinates) return vmecCoord;

  volatile CStconfig::ctStack::saveState state(ct);
  calcVmecLgB(vmecCoord,0,true);
  Vector3d boozer = vmecCoord;
  boozer[1] += ct().lambda + iota(boozer[0])*ct().p;
  boozer[2] += ct().p;
  return boozer;
}

//****************************************************************************
// @param boozer is the Boozer coordinates (s,theta,phi)
// @return Vmec coordinates
Vector3d CStconfig::Boozer2Vmec(const Vector3d &boozer)
{
  if(!vmecCoordinates) return boozer;
  volatile CStconfig::ctStack::saveState state(ct);
  // find v as root of the expression:   boozer - Vmec2Boozer(v) = 0

  Vector3d v=boozer;
  double eps2=1e-20;
  // Vector3d bz = Vmec2Boozer(v);
  // Vector3d f = boozer - bz;
  calcVmecLgB(v,0,true);
  double iota =  SPLINEy(miota);
  Vector3d bz = v;
  bz[1] += ct().lambda + iota*ct().p;
  bz[2] += ct().p;
  Vector3d f = boozer - bz; // f is the residual
  int iter = 1;
  for(; iter<300; ++iter) {// Newton's iterations
    // Jacobian matrix
    double a11 = 1 + ct().lambdat + iota*ct().pt; // dtheta_b / dtheta
    double a12 =     ct().lambdap + iota*ct().pp; // dtheta_b / dphi
    double a21 =     ct().pt;                     // dphi_b / dtheta
    double a22 = 1 + ct().pp;                     // dphi_b / dphi
    double det = a11*a22-a12*a21;
    // inverse matrix 
    // double b11 =  a22, b12 = -a12;
    // double b21 = -a21, b22 =  a11;     
    // double dthe_v = (b11*f[1]+b12*f[2])/det;
    // double dphi_v = (b21*f[1]+b22*f[2])/det;
    double dthe_v = ( a22*f[1]-a12*f[2])/det;
    double dphi_v = (-a21*f[1]+a11*f[2])/det;
    Vector3d d(0, dthe_v, dphi_v);    // correction to root
// backtracking loop: try the full Newton step; shorten the step size 
// if the residual does not decrease
    double step = 1;
    double fOldabs2 = f.abs2();
    Vector3d vOld = v;
    int N=5;
    while(--N) {
      v = vOld + step*d; // try the full Newton step
      { //this block is calculates  bz = Vmec2Boozer(v);
        calcVmecLgB(v,0,true);
        double iota =  SPLINEy(miota);
        bz = v;
        bz[1] += ct().lambda + iota*ct().p;
        bz[2] += ct().p;
      }
      f = boozer - bz; // f is the new residual
      if(f.abs2()<eps2) return v;  
      if(f.abs2()>fOldabs2) step *= 0.5;
      else break; // step is OK
    }
  }
  std::cerr << "CStconfig::Boozer2Vmec(): # of iterations " <<iter<<"  s="<<boozer[0]<<std::endl<<std::flush;
  return v;  // return with (probably bad) root
}

//****************************************************************************
// The method is used only for transformation from vmec to boozer coordinates by interrating over vmec angles
void CStconfig::V2Btrans(const Vector3d &vmec, Vector3d &boozer, double &R, double &Z, double &p, double &B, double &jv_jb)
{
  if(!vmecCoordinates) return;
  calcVmecLgB(vmec,0,true);
  p     = ct().p;
  double iota =  SPLINEy(miota);
  boozer = vmec;
  boozer[1] += ct().lambda + iota*ct().p;
  boozer[2] += ct().p;
  jv_jb = (1+ct().lambdat)*(1+ct().pp) + ct().pt*(iota-ct().lambdap); // d(theta_booz,phi_booz)/d(theta_vmec,phi_vmec)
  Vector3d cyl = mag2cyl(vmec);
  R     = cyl[0];
  Z     = cyl[2];
  B     = CStconfig::B(vmec);
  //////{
  //////  double jvmec = ct().jacVmec;      // Jacobian in coordinates (s,theta,phi)
  //////  double iot  = SPLINEy(miota);
  //////  double Ipol = SPLINEy(mIpol)*mBnorm;  // interpolate
  //////  double Itor = SPLINEy(mItor)*mBnorm;
  //////  double Babs = CStconfig::B(vmec);
  //////  double jbooz = mSignJac_s*fabs( Flux(1.)*mu0*(Ipol*mNp+iot*Itor)/square(twopi*Babs) );// Jacobian in coordinates (s,theta,phi)
  //////  jv_jb = jvmec/jbooz;
  //////}
}

//****************************************************************************
// private: used only for transformation from vmec to boozer coordinates
void CStconfig::V2Btrans(Vector3d &boozer, double &R, double &Z, double &p, double &B)
{
  if(!vmecCoordinates) return;
  Vector3d vmec = Boozer2Vmec(boozer);    
  p     = boozer[2]-vmec[2]; //  = ct().p;
  Vector3d cyl = mag2cyl(vmec);
  R     = cyl[0];
  Z     = cyl[2];
  B     = CStconfig::B(vmec);
}

//****************************************************************************
// VMEC coordinates vmecCoord=(s,theta,phi)
// @return 
double CStconfig::Vmec2BoozerTest(const Vector3d &vmecCoord)
{
  if(!vmecCoordinates) return 0;
  volatile CStconfig::ctStack::saveState state(ct);

  Vector3d boozer, cyl;
  double R,Z,p, B, jv_jb;
  V2Btrans(vmecCoord, boozer, R,Z, p, B, jv_jb);

  Vector3d Bcontra;
  calcVmecLgB(vmecCoord,&Bcontra,true);
  double w1 = (1+ct().lambdat)*(1+ct().pp) + ct().pt*(iota(vmecCoord[0])-ct().lambdap); // d(theta_booz,phi_booz)/d(theta_vmec,phi_vmec)
  double w2 = (Bcontra[2]*(1+ct().pp) + ct().pt*Bcontra[1]) * (ct().jac*twopi);         // d(theta_booz,phi_booz)/d(theta_vmec,phi_vmec)
  
  {
    double jvmec = ct().jacVmec;      // Jacobian in coordinates (s,theta,phi)
    double iot  = SPLINEy(miota);
    double Ipol = SPLINEy(mIpol)*mBnorm;  // interpolate
    double Itor = SPLINEy(mItor)*mBnorm;
    double Babs = CStconfig::B(vmecCoord);
    double jbooz = mSignJac_s*fabs( Flux(1.)*mu0*(Ipol*mNp+iot*Itor)/square(twopi*Babs) );// Jacobian in coordinates (s,theta,phi)
    return jv_jb/(jvmec/jbooz); //must be one
  }

  return w1/jv_jb;   //must be one
  return w2/w1;      // must be one 
}

//****************************************************************************
// fi(s,theta,phi) = phi - (2*pi/Np)*sum_mn Phi_mn(s)*sin(m*theta-n*Np*phi)
// fi  is a cylindrical angle
// (s,theta,phi) is a point in magnetic coordinates.
// OUTPUT
//  dfi_dphi -- partial derivative dfi/dphi
// function returns fi -- cylindrical angle
double CStconfig::fi(double s,double theta,double phi,double &dfi_dphi,bool onFildLine) const
{
  if(tokamakEQDSK||vmecCoordinates||fpCoordinates) {
    dfi_dphi = 1;
    return phi;
  }
  return fi(s,::cos(theta),::sin(theta),::cos(mNp*phi),::sin(mNp*phi),phi,dfi_dphi,onFildLine);
}

double CStconfig::fi(const Vector3d &mCrd,double &dfi_dphi,bool onFildLine) const
{
  return fi(mCrd[0],mCrd[1],mCrd[2],dfi_dphi,onFildLine);
}

//****************************************************************************
// fi(s,theta,phi) = phi - (2*pi/Np)*sum_mn Phi_mn(s)*sin(m*theta-n*Np*phi)
// fi  is a cylindrical angle
// (s,theta,phi) is a point in magnetic coordinates.
// OUTPUT
// function returns fi -- cylindrical angle
//  dfi_dphi -- partial derivative dfi/dphi 
//  dfi_dphi along field line if onFildLine==true
// dfi/dphi_along_field is calculated asuming theta = theta0 + iota*phi 
// dfi/dphi_along_field i= phi - (twopi/Np)*sum Phi_mn(s)*(m*iota-n*Np)*cos(m*theta-n*Np*phi)
//
double CStconfig::fi(double s,double cst,double snt,double csNp,double snNp,
                     double phi,double &dfi_dphi,bool onFildLine) const
{
  if(!mOK) return 0;
  if(tokamakEQDSK||vmecCoordinates||fpCoordinates) {
    dfi_dphi = 1;
    return phi;
  }

  int M,N;
  int is = getSsegMN(s, M, N); // search nearest segment in array ms;
  if(s==0) M=0;
  const double * fC = mCOEF.constArray()+mN*offsetN+ halfMesh.is0*offsetS;
  const splineCoeffs * hMesh = &halfMesh;

  double iota = onFildLine?(SPLINEy(miota)) : 0.0;
  double iotaNp = iota/mNp;

  double fsum = 0, fpsum = 0;
  CComplex expNph(csNp,snNp);               // exp(1i*mNp*phi), where 1i is the sqrt(-1)
  CComplex exp_Nph = expNph.conj();         // exp(-1i*mNp*phi)
  CComplex expth (cst, snt);                // exp(1i*theta)
  CComplex expmth(1,0);                     // will be exp(1i*m*theta)
  for(int m=0; m<=M; m++) {
    double miotaNp = onFildLine?(m*iotaNp):0;  // m*iota/Np
    CComplex expmthNph =expmth; // exp(1i*m*theta+1i*mNp*phi), where 1i is the sqrt(-1)
    CComplex expmth_Nph=expmth; // exp(1i*m*theta-1i*mNp*phi), where 1i is the sqrt(-1)
    const double * fCoeff = fC + m*offsetM;
    for(int n=0; n<= N; n++) {
      double iPmn, iPm_n;       // intrepolated value of Pmn in point s
      int idx  =  n*offsetN;
      SPLINE(Phimn,  idx,iPmn ,hMesh)
      SPLINE(Phimn, -idx,iPm_n,hMesh)
      if(n==0) {
        fsum += iPmn*expmth.im;           // expmth.im==sin(m*theta);
      //fpsum -=       n*iPmn*expmth.re;  //dfi/dphi
        fpsum += miotaNp*iPmn*expmth.re;  //d(Phi)/dtheta
        continue;
      }
      expmth_Nph *= exp_Nph; // exp(1i*m*theta-1i*mNp*phi)
      expmthNph  *= expNph;  // exp(1i*m*theta+1i*mNp*phi)
      fsum  +=              iPmn*expmth_Nph.im+iPm_n*expmthNph.im;
      fpsum += (miotaNp-n)*(iPmn*expmth_Nph.re-iPm_n*expmthNph.re);    //dfi/dphi    /////iPmn*cos(m*theta-n*Np*phi)+iPm_n*cos(m*theta+n*Np*phi)
    }
    expmth *= expth;   // exp(1i*m*theta)
  }
  dfi_dphi = 1-twopi*fpsum; // dfi/dphi
  return phi-mPeriod*fsum;
}

//****************************************************************************
// function returns vector (x,y,z) -- Cartesian  coordinates
// of the point (s,theta,phi) in magnetic coordinates.
//
Vector3d CStconfig::mag2xyz(const Vector3d &magCoord) const
{
  return mag2xyz(magCoord[0],magCoord[1],magCoord[2]);
};

Vector3d CStconfig::mag2xyz(double s, double theta, double phi) const
{
  Vector3d c = mag2cyl(s,theta,phi);
  return Vector3d(c[0]*::cos(c[1]),c[0]*::sin(c[1]),c[2]);
}

//****************************************************************************
// return vector (R,fi,Z) for point (s,theta,phi), where
// (s,theta,phi) is a point in magnetic coordinates.
// fi -- cylindrical angle
Vector3d CStconfig::mag2cyl(const Vector3d &magCoord) const
{
  return mag2cyl(magCoord[0],magCoord[1],magCoord[2]);
}

//****************************************************************************
// function returns vector (R,fi,Z) -- cylindrical coordinates
// of point (s,theta,phi), where  (s,theta,phi) is magnetic coordinates.
//
Vector3d CStconfig::mag2cyl(double s, double theta, double phi) const
{
  if(!mOK) return Vector3d(0.,0.,0.);
  double rs=0, zs=0, fs=0;

  int M,N;
  int is = getSsegMN(s, M, N); // search nearest segment in array ms;
  if(s==0) { M=0; theta=0; }

  bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);

  const double * fCh = mCOEF.constArray()+mN*offsetN+ halfMesh.is0*offsetS;
  const double * fCf = mCOEF.constArray()+mN*offsetN+ fullMesh.is0*offsetS;
  const splineCoeffs * hMesh = &halfMesh;
  const splineCoeffs * RZmesh = vmecCoordinates? &fullMesh:&halfMesh;

  double csNp = tokamakEQDSK?0:(::cos(mNp*phi)); //if(tokamakEQDSK) phi=pi/2;
  double snNp = tokamakEQDSK?1:(::sin(mNp*phi));

  CComplex expNph(csNp,snNp);               // exp(1i*mNp*phi), where 1i is the sqrt(-1)
  CComplex exp_Nph = expNph.conj();         // exp(-1i*mNp*phi)
  CComplex expth (::cos(theta), ::sin(theta));  // exp(1i*theta)
  CComplex expmth(1,0);                     // will be exp(1i*m*theta)
  for(int m=0; m<=M; m++) {
    CComplex expmthNph =expmth; // exp(1i*m*theta+1i*mNp*phi), where 1i is the sqrt(-1)
    CComplex expmth_Nph=expmth; // exp(1i*m*theta-1i*mNp*phi), where 1i is the sqrt(-1)
    const double * fCoeff     = fCh + m*offsetM;
    const double * fCoeffFull = fCf + m*offsetM;
    const double * fCoeffPointer = vmecCoordinates?fCoeffFull:fCoeff;
    double iRmn, iZmn, iPmn(0);       // intrepolated value of Rmn in point s
    double iRm_n,iZm_n,iPm_n(0);      // intrepolated value of Rm,-n in point s
    for(int n=0; n<= N; n++) {
      int idx  =  n*offsetN;
      SPLINE(Rmn,  idx,iRmn ,RZmesh)
      SPLINE(Zmn,  idx,iZmn ,RZmesh)
      SPLINE(Rmn,  -idx  ,iRm_n,RZmesh)
      SPLINE(Zmn,  -idx  ,iZm_n,RZmesh)
      if(boozerCrd) {
        SPLINE(Phimn,idx,iPmn ,hMesh)
        SPLINE(Phimn,-idx ,iPm_n,hMesh)
      }
      if(n==0) {
        rs += iRmn*expmth.re;  // *cos(m*theta)
        zs += iZmn*expmth.im;  // *sim(m*theta)
        if(boozerCrd) 
          fs += iPmn*expmth.im;  // *sim(m*theta)
        continue;
      }
      expmth_Nph *= exp_Nph; // exp(1i*m*theta-1i*mNp*phi)
      expmthNph  *= expNph;  // exp(1i*m*theta+1i*mNp*phi)
      rs += iRmn*expmth_Nph.re + iRm_n*expmthNph.re; //iRmn*cos(m*theta-n*Np*phi)+iRm_n*cos(m*theta+n*Np*phi)
      zs += iZmn*expmth_Nph.im + iZm_n*expmthNph.im;
      if(boozerCrd) fs += iPmn*expmth_Nph.im + iPm_n*expmthNph.im;
    }
    expmth *= expth;   // exp(1i*m*theta)
  }
  fs = phi - mPeriod*fs; // cylindrical angle
  const_cast<CStconfig*>(this)->mMsum2 += M;
  const_cast<CStconfig*>(this)->mNsum2 += N;
  const_cast<CStconfig*>(this)->mNmag2cyl++;
  return Vector3d(rs,fs,zs);
}

//****************************************************************************
Vector3d CStconfig::mixcoord2xyz(const Vector3d &mixCoord,double* Phi) const
{
  Vector3d c = mixcoord2cyl(mixCoord[0],mixCoord[1],mixCoord[2],Phi);
  return Vector3d(c[0]*::cos(c[1]),c[0]*::sin(c[1]),c[2]);
}

//****************************************************************************
Vector3d CStconfig::mixcoord2cyl(const Vector3d &mixCoord,double* Phi) const
{
  return mixcoord2cyl(mixCoord[0],mixCoord[1],mixCoord[2],Phi);
}

//****************************************************************************
// Find magnetic toroidal angle phi and cylindrical coord. (R,f,Z)
// for point (s,theta,f), where theta is magnetic poloidal angle,
// f is cylindrical angle
// INPUT:
//  (s,theta,f)
// OUTPUT:
//  Phi -- magnetic toroidal angle
//  function returns vector (R,f,Z) -- cylindrical coordinates
Vector3d CStconfig::mixcoord2cyl(double s, double theta, double f,double * Phi,
                                 sumDat * sdat) const
{
  if(!mOK) return Vector3d(0,0,0);

  double phi(f);
  double csNp;
  double snNp;
  const double * fCh;
  const double * fCf;
  const splineCoeffs * RZmesh;
  int M,N;

  if(sdat) {   //&&sdat->fi==fi&&sdat->s==s) {
    M  = sdat->M;
    N  = sdat->N;
    fCh =  sdat->fCh;
    fCf =  sdat->fCf;
    RZmesh =  sdat->RZmesh;
    csNp =  sdat->csNp;
    snNp =  sdat->snNp;
  }
  else {
    getSsegMN(s, M, N); // search nearest segment in array ms;
    if(s==0) theta=0;
    fCh = mCOEF.constArray()+mN*offsetN+ halfMesh.is0*offsetS;
    fCf = mCOEF.constArray()+mN*offsetN+ fullMesh.is0*offsetS;
    RZmesh = vmecCoordinates? &fullMesh:&halfMesh;
    csNp = ::cos(mNp*phi);
    snNp = ::sin(mNp*phi);
  }

  const double cst = ::cos(theta);
  const double snt = ::sin(theta);

#if 1
  if( !(tokamakEQDSK||vmecCoordinates||fpCoordinates) ) {
    // find magnetic toroidal angle phi for the given cylindrical angle fi
    //     f - fi(is,theta,phi)  = 0
    const double eps=1e-6;
    double dfi_dphi;
    double f1 = fi(s,cst,snt,csNp,snNp,phi,dfi_dphi); //fi(is,theta,phi);
    double dF = f-f1;   // residual
    for(int i=1,iter=200; i<=iter; ++i) {
      double dphi = dF/dfi_dphi;  // correction to root
      // backtracking loop: 
      double phiOld = phi;
      double dFold  = fabs(dF);
      double step = 1;
      int N=5;
      while(--N) {  // backtracking loop:
        phi = phiOld + step*dphi; // try the Newton step
        csNp = ::cos(mNp*phi);
        snNp = ::sin(mNp*phi);
        f1 = fi(s,cst,snt,csNp,snNp,phi,dfi_dphi); //fi(is,theta,phi);
        dF = f-f1;   // new residual
        if(fabs(dF)<eps) break;
        if(fabs(dF)>dFold) step *= 0.5; // shorten the step size if the residual does not decrease
        else break; // step is OK
      }
      if(fabs(dF)<eps) break;  
      if(i<iter) continue;
      std::cerr << "CStconfig::mixcoord2cyl(): # of angle iterations " <<i;
      //if(s>1) std::cerr << "; s = "<<s<<", must be <=1";
      std::cerr<<std::endl;
      break;
    }
  }
#else  // old version
  if( !(tokamakEQDSK||vmecCoordinates||fpCoordinates) ) {
    for(int i=1,iter=200; i<=iter; ++i) {
      const double eps=1e-6;
      double dfi_dphi;
      double f1 = fi(s,cst,snt,csNp,snNp,phi,dfi_dphi); //fi(is,theta,phi);
      if(fabs(f-f1)<eps) break;
      double dphi = (f-f1)/dfi_dphi;  // Newton's iterations
      ////if(i>5) { 
      ////  dphi *=0.5; //try to backtrack if convergence is too slow
      ////  //std::cerr <<i<<std::endl<<std::flush;
      ////}
      phi += dphi;
      csNp = ::cos(mNp*phi);
      snNp = ::sin(mNp*phi);
      if(i<iter) continue;
      std::cerr << "CStconfig::mixcoord2cyl(): # of angle iterations " <<i;
      //if(s>1) std::cerr << "; s = "<<s<<", must be <=1";
      std::cerr<<std::endl;
      break;
    }
  }
#endif
//phi=phimag(s,theta,f);
//csNp = ::cos(mNp*phi);
//snNp = ::sin(mNp*phi);


  csNp = tokamakEQDSK?0:csNp; //if(tokamakEQDSK) phi=pi/2;
  snNp = tokamakEQDSK?1:snNp;
  CComplex expNph(csNp,snNp);               // exp(1i*mNp*phi), where 1i is the sqrt(-1)
  CComplex exp_Nph = expNph.conj();         // exp(-1i*mNp*phi)
  CComplex expth (cst, snt);                // exp(1i*theta)
  CComplex expmth(1,0);                     // will be exp(1i*m*theta)

  double rs=0, zs=0;
  for(int m=0; m<=M; m++)  {
    CComplex expmthNph =expmth; // exp(1i*m*theta+1i*mNp*phi), where 1i is the sqrt(-1)
    CComplex expmth_Nph=expmth; // exp(1i*m*theta-1i*mNp*phi), where 1i is the sqrt(-1)
    const double * fCoeff     = fCh + m*offsetM;
    const double * fCoeffFull = fCf + m*offsetM;
    const double * fCoeffPointer = vmecCoordinates?fCoeffFull:fCoeff;
    for(int n=0; n<= N; n++) {
      double iRmn, iZmn ;       // intrepolated value of Rmn in point s
      double iRm_n,iZm_n;      // intrepolated value of Rm,-n in point s
      int idx  =  n*offsetN;
      SPLINE(Rmn,  idx,iRmn ,RZmesh)
      SPLINE(Zmn,  idx,iZmn ,RZmesh)
      SPLINE(Rmn, -idx,iRm_n,RZmesh)
      SPLINE(Zmn, -idx,iZm_n,RZmesh)
      if(n==0) {
        rs += iRmn*expmth.re;  // *cos(m*theta)
        zs += iZmn*expmth.im;  // *sim(m*theta)
        continue;
      }
      expmth_Nph *= exp_Nph; // exp(1i*m*theta-1i*mNp*phi)
      expmthNph  *= expNph;  // exp(1i*m*theta+1i*mNp*phi)
      rs += iRmn*expmth_Nph.re + iRm_n*expmthNph.re; //iRmn*cos(m*theta-n*Np*phi)+iRm_n*cos(m*theta+n*Np*phi)
      zs += iZmn*expmth_Nph.im + iZm_n*expmthNph.im;
    }
    expmth *= expth;   // exp(1i*m*theta)
  }
  if(Phi) *Phi=phi;
  const_cast<CStconfig*>(this)->mMsum2 += M;
  const_cast<CStconfig*>(this)->mNsum2 += N;
  const_cast<CStconfig*>(this)->mNmag2cyl++;
  return Vector3d(rs,f,zs);
}

//****************************************************************************
// Find magnetic toroidal angle phi for point (s,theta,f), where
// theta is magnetic poloidal angle, f is cylindrical angle
// INPUT:
//  (s,theta,f)
// OUTPUT:
//  Phi -- magnetic toroidal angle
//  function returns magnetic toroidal angle
double CStconfig::phimag(const Vector3d &mixCoord) const {
  return phimag(mixCoord[0],mixCoord[1],mixCoord[2]);
};

double CStconfig::phimag(double s,double theta,double f) const
{
  if(!mOK) return 0;

  if( tokamakEQDSK||vmecCoordinates||fpCoordinates ) return f;
  
  if(s==0) theta=0;

  double phi = f;
  double csNp = ::cos(mNp*phi); 
  double snNp = ::sin(mNp*phi);
  double cst = ::cos(theta);
  double snt = ::sin(theta);
#if 1
  // find magnetic toroidal angle phi for the given cylindrical angle fi
  //     fi(is,theta,phi) - f = 0
  const double eps=1e-8;
  double dfi_dphi;
  double f1 = fi(s,cst,snt,csNp,snNp,phi,dfi_dphi); //fi(is,theta,phi);
  double dF = f-f1;   // residual
  for(int i=1,iter=200; i<=iter; ++i) {
    double dphi = dF/dfi_dphi;  // correction to root
    // backtracking loop: 
    double phiOld = phi;
    double dFold  = fabs(dF);
    double step = 1;
    int N=5;
    while(--N) {  // backtracking loop:
      phi = phiOld + step*dphi; // try the Newton step
      csNp = ::cos(mNp*phi);
      snNp = ::sin(mNp*phi);
      f1 = fi(s,cst,snt,csNp,snNp,phi,dfi_dphi); //fi(is,theta,phi);
      dF = f-f1;   // new residual
      if(fabs(dF)<eps) break;
      if(fabs(dF)>dFold) step *= 0.5; // shorten the step size if the residual does not decrease
      else break; // step is OK
    }
    if(fabs(dF)<eps) break;  
    if(i<iter) continue;
    std::cerr << "CStconfig::phimag(): # of poloidal angle iterations " <<i;
    //if(s>1) std::cerr << "; s = "<<s<<", must be <=1";
    std::cerr<<std::endl;
    break;
  }
  return phi;  
#else  
// find Boozer toroidal angle ph for the given cylindrical angle fi
  for(int i=1,iter=200; i<=iter; ++i) {
    const double eps=1e-5;
    double dfi_dphi;
    double f1 = fi(s,cst,snt,csNp,snNp,phi,dfi_dphi); //fi(is,theta,phi);
    double dphi = (f-f1)/dfi_dphi;  // Newton's iterations
    phi += dphi;
    csNp = ::cos(mNp*phi);
    snNp = ::sin(mNp*phi);
    if(fabs(dphi)<eps) break;
    if(i<iter) continue;
    std::cerr << "CStconfig::phimag(): # of poloidal angle iterations " <<i;
    if(s>1) std::cerr << "; s = "<<s<<", must be <=1";
    std::cerr<<std::endl;
    break;
  }
  return phi;
#endif
}

//****************************************************************************
// function returns the size of rectangle circumscribed about the cross
// section of the 's' surface, 'f' is the cylindrical angle
void CStconfig::getExtent(double &Rmin, double &Rmax, double &Zmin, double &Zmax, double s, double f) const
{
  int itmax = tokamakConfig?500:150;  // # of theta-points (poloidal)
  Rmin = Zmin =  1e15;
  Rmax = Zmax = -1e15;
  sumDat sdat(const_cast<CStconfig&>(*this),s,f);
  for(int i=0; i<itmax; i++) {
    double theta = i*twopi/itmax;
    Vector3d c = mixcoord2cyl(s,theta,f,0,&sdat);
    Rmin = mmin(Rmin,c[0]);
    Rmax = mmax(Rmax,c[0]);
    Zmin = mmin(Zmin,c[2]);
    Zmax = mmax(Zmax,c[2]);
  }
}

//****************************************************************************
// function returns the sizes of parallelepiped circumscribed about the segment f1-f2<pi
// of the 's' surface, 'f' is the cylindrical angle
void CStconfig::getExtent(Vector3d &rmin, Vector3d &rmax, double s, double f1, double f2) const
{
  rmin.set(1e15,1e15,1e15);
  rmax.set(-1e15,-1e15,-1e15);
  sortAngles(f1,f2);
  int i,nt = tokamakConfig?500:150;  // # of theta-points (poloidal)
  int j,nf = 1;
  double df=0;
  if(f2>f1) {
//    f1-=degree; // add some margin
//    f2+=degree;
    nf = 20;
    double df = (f2-f1)/nf;
    if(df>degree) df=degree;
    nf=int((f2-f1)/df);
    df=(f2-f1)/(nf-1);
  }

  nf = tokamakConfig?1:nf;

  for(j=0; j<nf; j++) {
    double f = f1+j*df;
    for(i=0; i<nt; i++) {
      double theta = i*twopi/nt;
      Vector3d b(s,theta,f);
      Vector3d c = mixcoord2xyz(b);
      rmin.x() = mmin(rmin.x(),c[0]);
      rmin.y() = mmin(rmin.y(),c[1]);
      rmin.z() = mmin(rmin.z(),c[2]);
      rmax.x() = mmax(rmax.x(),c[0]);
      rmax.y() = mmax(rmax.y(),c[1]);
      rmax.z() = mmax(rmax.z(),c[2]);
    }
  }
}

//****************************************************************************
// function sorts angles, so that f1<f2 and f2-f1<pi
void CStconfig::sortAngles(double &f1, double &f2) const
{
  f1 = mod2pi(f1);
  f2 = mod2pi(f2);
  if(f2<f1) Swap(f2,f1);
  if(f2-f1>pi) { f1+=twopi; Swap(f2,f1);}
}

//****************************************************************************
// function returns the size of rectangle circumscribed about the cross
// section of the surface 's'
void CStconfig::getExtent(double &Rmin, double &Rmax, double &Zmin, double &Zmax, double s) const
{
  int itmax = tokamakConfig?360:100;
  int ipmax = tokamakConfig?1:100;  // Tokamak ???

  Rmin = Zmin =  1e15;
  Rmax = Zmax = -1e15;
  double dphi = mPeriod/ipmax;

  if(Ntor()==0) ipmax=1; // tokamakSymmetry=true;

  for(int k=0; k<ipmax; k++) {
    double ph = k*dphi;
    for(int i=0; i<itmax; i++) {
      double theta = i*twopi/itmax;
      Vector3d c = mag2cyl(s, theta, ph);
      Rmin = mmin(Rmin,c[0]);
      Rmax = mmax(Rmax,c[0]);
      Zmin = mmin(Zmin,c[2]);
      Zmax = mmax(Zmax,c[2]);
    }
  }
}

//****************************************************************************
// function returns the sizes of rectangle circumscribed about the cross
// section of the last closed magnetic surface(LCMS)
void CStconfig::getExtentLCMS(double &Rmin, double &Rmax, double &Zmin, double &Zmax) const
{
  const_cast<CStconfig*>(this)->lcms.create(this);   // will be created if needed
  Rmin = lcms.mRmin;
  Rmax = lcms.mRmax;
  Zmin = lcms.mZmin;
  Zmax = lcms.mZmax;
}

//****************************************************************************
//   note: Jacobian in coordinates (s,theta,phi)
// Calculate Jacobian:
//    Boozer coordinates  Flux(1.)*mu0*(Ipol*Np+iota*Itor)/B^2/twopi^2
//    PEST  coordinates   Flux(1.)*R*R/(mu0*Ipol) 
//    VMEC coordinates  sum of Fourier mods
//    fpCoordinates    mixed product of covariant basis vectors
// INPUT:
//  Vector3d (s,theta,phi) is a point in magnetic coordinates,
// OUTPUT:
// Function returns the value of Jacobian in the point (s,theta,phi)
double CStconfig::Jacobian(const Vector3d &magCoord, double *Bmod) const
{
  if(vmecCoordinates||fpCoordinates) {
    if(Bmod) *Bmod = B(magCoord);
    return const_cast<CStconfig*>(this)->JacobianVmec(magCoord); //  Jacobian in coordinates (s,theta,phi)
  }

  const_cast<CStconfig*>(this)->calculateCurrents();

  double s=magCoord[0];
  int M,N;
  int is = getSsegMN(s, M, N);
 
  double Ipol = SPLINEy(mIpol)*mBnorm;  // interpolate
  if(tokamakEQDSK) {
    double R = mag2cyl(magCoord)[0];
    if(Bmod) *Bmod = B(magCoord);
    return mSignJac_s*fabs( Flux(1.)*R*R/(mu0*Ipol) );
  }
  double iota = SPLINEy(miota);
  double Itor = SPLINEy(mItor)*mBnorm;
  double Babs = B(magCoord);
  if(Bmod) *Bmod = Babs;
  return mSignJac_s*fabs( Flux(1.)*mu0*(Ipol*mNp+iota*Itor)/square(twopi*Babs) );
}

//****************************************************************************
// note: VMEC Jacobian in coordinates (s,theta,phi)
double CStconfig::JacobianVmec(const Vector3d &magCoord)
{
  volatile CStconfig::ctStack::saveState state(ct);
  calcVmecLgB(magCoord);
  if(!fpCoordinates) 
    return mSignJac_s*fabs(ct().jacVmec); //  Jacobian in coordinates (s,theta,phi)
  if(fpCoordinates) {
    NJacobian(magCoord);
    return  getJ2();
  }
  return 0;
}

//****************************************************************************
// Volume[m^3]
// Integral(0,s) V'*ds
// Volume is negative if (s,theta,phi)-coordinate system is left handed
double CStconfig::VolumeInt(double s) const
{
  double sum=0, ds = 0.002;
  if(s==0) return 0;
  int n=int(s/ds);
  if(n<8) n=8;
  n=n+n%2;  //n must be even for simpson integration
  ds=s/n;
//Simpson integration
  for(int k=0; k<=n; k++) {
    double w = (k==0||k==n)?1:(2*(1+k%2));   // weight for Simpson integration
    sum += w*Vprime(k*ds);
  }
  return (sum*ds/3);
};

//****************************************************************************
// Volume[m^3]
// Integral(s1,s2) V'*ds
// Volume is negative if (s,theta,phi)-coordinate system is left handed
double CStconfig::VolumeInt(double s1, double s2) const
{
  double sum=0, ds=0.004;
  double Ds=s2-s1;
  if(Ds==0) return 0;
  int n=int(Ds/ds);
  if(n<4) n=4;
  n=n+n%2;  //must be even for simpson integration
  ds=Ds/n;
  for(int k=0; k<=n; k++) {  //Simpson integration
    double w = (k==0||k==n)?1:(2*(1+k%2));   // weight for Simpson integration
    sum += w*Vprime(s1+k*ds);
  }
  return (sum*ds/3);
};

//****************************************************************************
//
void CStconfig::adjustNumbersForSimpson(int &Npol, int &Ntor, double &dpol, double &dtor) const
{
  int ipol = tokamakConfig?500:Npol;        // theta-points (poloidal)
  int itor = tokamakConfig?4:(Ntor/mNp);    // phi-points (toroidal)

  ipol = mmax(10,  ipol);
  ipol = mmax(Npol,ipol);

  itor = mmax(4,itor);

  Npol = ipol + ipol%2;  //must be even for simpson integration
  Ntor = itor + itor%2;

  Ntor = tokamakConfig?0:Ntor;  // Tokamak ???

  dtor = tokamakConfig?(3*twopi):mPeriod/Ntor;
  dpol = twopi/Npol;
}

//****************************************************************************
// dV/ds    [m^3]
// V' as integral from Jacobian
//  for testing only, very slow, don't use it
//    or use to tabulate mg00
//      for(is=0; is<=mNs; is++) mg00[is] = VprimeIntJ(ms[is])/mNp;
// itmax, ipmax are the number of integration points
//  itmax       theta-points (poloidal)
//  ipmax       phi-points (toroidal)
// Note: Ipol anr Itor must be known here, see also Jacobian()
double CStconfig::VprimeIntJ(double s, int ithmax, int ipmax)
{
  double sum=0;
  volatile ctStack::saveState state(ct);
  double dtheta, dphi;
  adjustNumbersForSimpson(ithmax,ipmax,dtheta,dphi);

//Simpson integration
  for(int i=0; i<=ithmax; i++) {
    double theta = i*dtheta;
    double wt = (i==0||i==ithmax)?1:(2*(1+i%2));    // weight for Simpson integration
    for(int k=0; k<=ipmax; k++) {
      double phi = k*dphi;
      double w = (k==0||k==ipmax)?1:(2*(1+k%2));   // weight for Simpson integration
      w *= wt;
      Vector3d magCoord(s,theta,phi);
      double jac = Jacobian(magCoord); 
//j    NJacobian(magCoord);
//j    jac = getJ(); // Jacobian in coordinates (s,theta,phi)
//j    jac = getJ2(); // jac=r*(Rs*(Rt^Rp)), Jacobian in coordinates (s,theta,phi)
//test  Vector3d Isum = iota*Ft+Fp;  //  (iota*dX/dtheta+dX/dphi)==B*jac*twopi
//test  sum += w*Isum.abs2();
      sum += w*jac;
    }
  }
  return sum * mNp*dtheta*dphi/9;
};

//****************************************************************************
// Susceptance matrix: S11,S12,S22 in (Flux,theta,phi)-coordinates
//
//  S_jk = 1/(4pi^2)*integral(dtheta*dphi*jac*g_jk/g)
//
// s is the flux surf. label: norm.flux;
//
// itmax, ipmax are the number of integration points
//  itmax       theta-points (poloidal)
//  ipmax       phi-points (toroidal)
//
Vector3d *CStconfig::SMatrix(double s, int ithmax, int ipmax)
{
  volatile ctStack::saveState state(ct);

  if(s<1e-4) s=1e-4;
  Vector3d SM(0.,0.,0.);
  double S21(0);

  double gradpsi=0;
  double dVdpsiG=0;
  double Isum=0;         // = Ipol+iota*Itor

  double dtheta, dphi;
  adjustNumbersForSimpson(ithmax,ipmax,dtheta,dphi);

  //Simpson integration
  for(int i=0; i<=ithmax; i++) {
    double theta = i*dtheta;
    double wt = (i==0||i==ithmax)?1:(2*(1+i%2));    // weight for Simpson integration
    for(int k=0; k<=ipmax; k++) {
      double phi = k*dphi;
      double w = (k==0||k==ipmax)?1:(2*(1+k%2));   // weight for Simpson integration
      w *= wt;
      Vector3d fluxCoord(s,theta,phi);
      NJacobian(fluxCoord);
      Vector3d F (ct().Jmatr[0]); // (R,fi,Z)
      Vector3d Fs(ct().Jmatr[1]);// d(R,fi,Z)/ds   -- partial derivative on s
      Vector3d Ft(ct().Jmatr[2]);// d(R,fi,Z)/dtheta
      Vector3d Fp(ct().Jmatr[3]);// d(R,fi,Z)/dphi
      double R  = F[0];
    // Multiplying by R we take into account the cylindrical coordinate system
      Ft[1] *=R; // R*dfi/dtheta
      Fp[1] *=R; // R*dfi/dphi
      // metric tensor in coordinates (flux,theta,phi)
      // g_ij = dX/du^i * dX/du^j , X - cartesian vector
      //  (u^0,u^1,u^2) = (flux,theta,phi)
      double g11 = Ft*Ft;
      double g12 = Ft*Fp;
      double g21 = g12;
      double g22 = Fp*Fp;
      double jac;
      if(tokamakEQDSK) {
        if(polCurrentNeeded) // be careful, currents can be unknown here
          jac = getJac2(); // jac=r/mFlux*(Rs*(Rt^Rp)), Jacobian in coordinates (Flux,theta,phi)
        else
          jac = getJac();  // jac = R*R/(mu0*Ipol), be careful, currents can be unknown here
      }
      else {
        if(torCurrentNeeded||polCurrentNeeded) // be careful, currents can be unknown here
          jac = getJac2(); // jac=r/mFlux*(Rs*(Rt^Rp)), Jacobian in coordinates (Flux,theta,phi)
        else
          jac = getJac();  // jac = I/B^2 if Boozer coordinates or sum of gmm if VMEC coordinates
      }

      S21+= w*        (g21)/jac;
      SM += w*Vector3d(g11,
                       g12*(1+ct().lambdat) - g11*ct().lambdap,
                       g22*(1+ct().lambdat) - g21*ct().lambdap)/jac;
#if 0
//////////begin:  following fragment is not needed
      Vector3d  BcylJac;         // Bcyl*jac*twopi == (iota-lambdap)*dX/dtheta+(1+lambdat)*dX/dphi
    //BcylJac = iota*Ft+Fp;      //  (iota*dX/dtheta+dX/dphi)==Bcyl*jac*twopi  
      BcylJac = ct().Bcyl*jac*twopi;        //  (iota-lambda,p)*dX/dtheta+(1+lambda,t)*dX/dphi
      dVdpsiG += w*fabs(jac);                   // will be dV/dpsi, dV/dpsi from geometry
      gradpsi +=w*getmodGradFlux()*fabs(jac);   // will be <|Grad(psi)|>
      Isum    +=w*BcylJac.abs2();         // will be Ipol+iota*Itor
//////////end   ct().BcylJac == Bcyl*jac*twopi == (iota-lambdap)*dX/dtheta+(1+lambdat)*dX/dphi
#endif
    }
  }
  Smatrix[0]=0;
  Smatrix[1]=0;
  Smatrix[2]=0;

  double h=mNp*dtheta*dphi/9;
  SM *= h/square(twopi);
  Smatrix[1][1] = SM[0];        //S_jk = 1/(4pi^2)*integral(dtheta*dphi*jac*g_jk/g)
  Smatrix[1][2] = SM[1];
  Smatrix[2][1] = S21*h/square(twopi);
  Smatrix[2][2] = SM[2];

  Smatrix[1][0]  = (iota(s)*Smatrix[1][1]+Smatrix[1][2])/mu0; // toroidal current
  Smatrix[2][0]  = (iota(s)*Smatrix[2][1]+Smatrix[2][2])/mu0; // poloidal current
//  Smatrix[2][0]  = Isum-iota(s)*Smatrix[1][0];  // poloidal current

#if 0
//////////begin:  following fragment is not needed
    dVdpsiG  *= h;                         // dVdpsi=dV/dpsi=integral(dtheta*dphi*jac)
  //V  dVds    =  dVdpsi*mFlux;            // dVds  =dV/ds  =integral(dtheta*dphi*jac)*mFlux
    double dVds    =  mSignJac_s*fabs(Vprime(s)); // dVds  =dV/ds  =integral(dtheta*dphi*jac)*mFlux
    double dVdpsi  =  dVds/Flux(1);

    gradpsi *= h/dVdpsi;                        // gradpsi=<|Grad(Psi)|>
    double grads = gradpsi/Flux(1);          // grads  =<|Grad(s)|>

    Isum    *= h/(dVdpsiG*mu0); // Isum=Ipol+iota*Itor=integral(dtheta*dphi*(iota*dX/dtheta+dX/dphi)^2)/(dVdpsi*mu0)

    double GradS    = fabs(grads);              // |<|Grad(s)|>|
    double Gradreff = GradS*0.5*ma/sqrt(s);     // |<|Grad(r_eff)|>|
    Smatrix[0][0]  = dVds;
    Smatrix[0][1]  = Gradreff;               // |<|Grad(r_eff)|>|
    Smatrix[0][2]  = GradS;                  // |<|Grad(s)|>|
//////////end
#endif
 
  return Smatrix;
};

//****************************************************************************
// Calculate flux surface average curvatures
// itmax, ipmax are the number of integration points
//  itmax       theta-points (poloidal)
//  ipmax       phi-points (toroidal)
//
// @return Vector3d k with k[0] = <|kN*R0|^2> and k[1] = <|kG*R0|^2>
//
Vector3d CStconfig::curvatureAvrg(double s, int ithmax, int ipmax)
{
  double sumkN=0,sumkG=0,sum=0;

  double dtheta, dphi;
  adjustNumbersForSimpson(ithmax,ipmax,dtheta,dphi);

  //Simpson integration
  for(int i=0; i<=ithmax; i++) {
    double theta = i*dtheta;
    double wt = (i==0||i==ithmax)?1:(2*(1+i%2));    // weight for Simpson integration
    for(int j=0; j<=ipmax; j++) {
      double phi = j*dphi;
      double w = (j==0||j==ipmax)?1:(2*(1+j%2));   // weight for Simpson integration
      w *= wt;
      Vector3d magCoord(s,theta,phi);
      Vector3d k = curvature(magCoord);
////Vector3d r;
////double B,gradR_kG;
////Vector3d k = advanceDataForEpsEff(magCoord, r, B, gradR_kG, 0);
      double kN  = square(k[0]*mR0);
      double kG  = square(k[1]*mR0);
      double jac = k[2];
      sum   += w*jac;
      sumkN += w*jac*kN;
      sumkG += w*jac*kG;
    }
  }
  return Vector3d(sumkN/sum,sumkG/sum,0);

// for test :  must be  sum = V'/V at return
  sum*=mNp*dtheta*dphi/9/Volume(1);
  return Vector3d(sum,sum,0);
};

//****************************************************************************
// @return covariant components of  grad(|B|) in flux coordinates
// dB/ds   -- partial derivative on s
// dB/dtheta
// dB/dphi
Vector3d  CStconfig::getCovaGradB(const Vector3d &u) const
{
  double s     = u[0];
  double phi   = u[2];
  double B (0);          // B
  double Bs(0);          // dB/ds   -- partial derivative on s
  double Bt(0);          // dB/dtheta
  double Bp(0);          // dB/dphi

  //double iot = iota(s); //04jun2011
  int M,N;
  int is = getSsegMN(s, M, N);

  const double * fC = mCOEF.constArray()+mN*offsetN+ halfMesh.is0*offsetS;
  const splineCoeffs * hMesh = &halfMesh;

  double csNp = tokamakEQDSK?0:(::cos(mNp*phi)); //if(tokamakEQDSK) phi=pi/2;
  double snNp = tokamakEQDSK?1:(::sin(mNp*phi));

  CComplex expNph(csNp,snNp);               // exp(1i*mNp*phi), where 1i is the sqrt(-1)
  CComplex exp_Nph = expNph.conj();         // exp(-1i*mNp*phi)
  CComplex expth (::cos(u[1]), ::sin(u[1]));  // exp(1i*theta)
  CComplex expmth(1,0);                     // will be exp(1i*m*theta)
  for(int m=0; m<=M; m++) {
    CComplex expmthNph =expmth; // exp(1i*m*theta+1i*mNp*phi), where 1i is the sqrt(-1)
    CComplex expmth_Nph=expmth; // exp(1i*m*theta-1i*mNp*phi), where 1i is the sqrt(-1)
    const double * fCoeff = fC + m*offsetM;
    for(int n=0; n<= N; n++) {
      double ibmn;       // intrepolated value of bmn in point s
      double ibm_n;      // intrepolated value of bm,-n in point s
      double ibmns;      // intrepolated value of dbmn/ds in point s
      double ibm_ns;     // intrepolated value of dbm,-n/ds in point s
      int idx  =  n*offsetN;
      //here
      // bmn(is,  m,n)==fCoeff[idx]
      // bmn(is+1,m,n)==fCoeff[idx1]
      // bmn(is,  m,n)==mCOEF[(n+mN)*offsetN + m*offsetM + is*offsetS]
      // it looks ugly, but exactly here I need speed
      SPLINED(bmn,  idx,ibmn, ibmns ,hMesh)
      SPLINED(bmn, -idx,ibm_n,ibm_ns,hMesh)
      if(n==0) {
        double csmt = expmth.re;  // cos(m*theta)
        double snmt = expmth.im;  // sin(m*theta)
        B  += ibmn*csmt;
        Bs += ibmns*csmt;
        Bt -= m*ibmn*snmt;
        continue;
      }
      expmth_Nph *= exp_Nph;    // exp(1i*m*theta-1i*mNp*phi)
      expmthNph  *= expNph;     // exp(1i*m*theta+1i*mNp*phi)
      double sn1 = expmth_Nph.im;  //sin(m*theta-n*Np*phi) is stored into sn1
      double cs1 = expmth_Nph.re;  //cos(m*theta-n*Np*phi)
      double sn2 = expmthNph.im;   //sin(m*theta+n*Np*phi)
      double cs2 = expmthNph.re;   //cos(m*theta+n*Np*phi)

////if(fabs(m*iot-n*mNp)<1e-1) continue;  //04jun2011
////if(fabs(m*iot+n*mNp)<1e-1) continue;  //04jun2011

      B  += ibmn*cs1  + ibm_n*cs2;
      Bs += ibmns*cs1 + ibm_ns*cs2;
      Bp += n*(ibmn*sn1 - ibm_n*sn2);
      Bt -= m*(ibmn*sn1 + ibm_n*sn2);
    }
    expmth *= expth;   // exp(1i*m*theta)
  }
  Bp *=mNp;
#if INTERPOLATION_ORDER == 3
//  if(s==0) std::cerr<< "NJacobian1Q: s==0"<<std::endl;
  Bs *= halfMesh.e2rcurr;
#endif
  if(tokamakEQDSK) Bp = 0;

  B  *= fabs(mBnorm);
  Bs *= fabs(mBnorm);// dB/ds   -- partial derivative on s
  Bt *= fabs(mBnorm);// dB/dtheta
  Bp *= fabs(mBnorm);// dB/dphi

  return Vector3d(Bs,Bt,Bp);
}

//****************************************************************************
// @return grad(|B|) in cylindrical coordinates
Vector3d  CStconfig::gradB(const Vector3d &u)
{
  volatile ctStack::saveState state(ct);
  Vector3d gr = getCovaGradB(u);
  NJacobian(u);
  Vector3d gradB = gr[0]*getGrads() + gr[1]*getGradtheta() + gr[2]*getGradphi();
  return gradB;
}

//****************************************************************************
// Get B and gradients in cylindrical coordinates
//  @param[in] Vector3d u(s,theta,phi) is a point in magnetic coordinates,
//  @param[out]  B  is the magnetic field  in cylindrical coordinates
//  @param[out]  gradB  is the grad(|B|) 
//  @param[out]  grads  is the grad(s)
//  @param[out]  gradTh is the grad(theta)
//  @param[out]  gradPh is the grad(phi) 
void  CStconfig::getBandGradients(const Vector3d &u, Vector3d &B, Vector3d &gradB, 
                                  Vector3d &grads, Vector3d &gradTh, Vector3d &gradPh)
{
  volatile ctStack::saveState state(ct);
  Vector3d gr = getCovaGradB(u);
  NJacobian(u);
  grads  = getGrads();
  gradTh = getGradtheta();
  gradPh = getGradphi();
  B      = getBcyl();
  gradB  = gr[0]*grads + gr[1]*gradTh + gr[2]*gradPh;
}

//****************************************************************************
// Get B and gradients in cartesian coordinates
//  @param[in] Vector3d u(s,theta,phi) is a point in magnetic coordinates,
//  @param[out]  B  is the magnetic field  in cartesian coordinates
//  @param[out]  gradB  is the grad(|B|) in cartesian coordinates
//  @param[out]  grads  is the grad(s)
//  @param[out]  gradTh is the grad(theta)
//  @param[out]  gradPh is the grad(phi) 
void  CStconfig::getBandGradientsxyz(const Vector3d &u, Vector3d &B, Vector3d &gradB, 
                                  Vector3d &grads, Vector3d &gradTh, Vector3d &gradPh)
{
  volatile ctStack::saveState state(ct);
  NJacobian(u);
  grads  = getGradsxyz();
  Vector3d g = getGradtheta();
  Vector3d cyl=getRcyl();
  double fi = cyl[1];
  double cs = ::cos(fi);
  double sn = ::sin(fi);
  // transform to cartesian coordinates  gradTh =  g.cylVector2Cartesian(cyl) 
  gradTh.x() = g[0]*cs-g[1]*sn;  
  gradTh.y() = g[0]*sn+g[1]*cs;  
  gradTh.z() = g[2];
  g = getGradphi();
  // transform to cartesian coordinates  gradPh =  g.cylVector2Cartesian(cyl) 
  gradPh.x() = g[0]*cs-g[1]*sn;  
  gradPh.y() = g[0]*sn+g[1]*cs;  
  gradPh.z() = g[2];
  B      = getBxyz();
  g = getCovaGradB(u);
  gradB  = g[0]*grads + g[1]*gradTh + g[2]*gradPh;
}

//****************************************************************************
// Get B and basis vectors in cartesian coordinates
//  @param[in]   u = (s,theta,phi) is the point magnetic-coordinates,
//  @param[out]  xyz  is the point cartesian-coordinates,
//  @param[out]  Bxyz  is the magnetic field  in cartesian coordinates
//  @param[out]  gradBxyz  is the grad(|B|) in cartesian coordinates
//  @param[out]  e1con, e2con, e3con are the contravariant-basis vectors e^i = grad(u^i) in cartesian coordinates
//  @param[out]  e1cov, e2cov, e3cov are the covariant-basis vectors e_i = dX/du^i in cartesian coordinates
void  CStconfig::getBandBasisVectorsxyz(const Vector3d &u, Vector3d &xyz, Vector3d &Bxyz, Vector3d &gradBxyz, 
                                     Vector3d &e1con, Vector3d &e2con, Vector3d &e3con,
                                     Vector3d &e1cov, Vector3d &e2cov, Vector3d &e3cov)
{
  volatile ctStack::saveState state(ct);
  NJacobian(u);
  Bxyz   = getBxyz();
  Vector3d cyl = getRcyl(); 
  xyz = cyl.toCartesian();

  Vector3d e1,e2,e3;
  getContraBasis(e1, e2, e3);
  e1con = e1.cylVector2Cartesian(cyl);  // transform to cartesian coordinates
  e2con = e2.cylVector2Cartesian(cyl);  
  e3con = e3.cylVector2Cartesian(cyl);  

  getCovaBasis(e1, e2, e3);
  e1cov = e1.cylVector2Cartesian(cyl);  
  e2cov = e2.cylVector2Cartesian(cyl);  
  e3cov = e3.cylVector2Cartesian(cyl);  
  
  Vector3d gradBcov, Bcon, Bcov;
  Bcov = getBcova();
  Bcon = getBcontra();
  gradBcov = getCovaGradB(u); //  Vector3d(dB/ds,dB/dtheta,dB/dphi)
  gradBxyz  = gradBcov[0]*e1con + gradBcov[1]*e2con + gradBcov[2]*e3con;
}

//****************************************************************************
// Calculate Jacobian matrix for Newton's iteration
//
// INPUT:
//  Vector3d (s,theta,phi) is a point in magnetic coordinates,
// OUTPUT:
// Function returns pointer to the array of 4 vectors
//    (R )(Rs )(Rt )(Rp )
// J= (fi)(fis)(fit)(fip)     in the point (s,theta,phi)
//    (Z )(Zs )(Zt )(Zp )
//
//  where Rs = dR/ds, Rt = dR/dtheta, Rp = dR/dphi
//  (s,theta,phi) is a point in magnetic coordinates,
//  (R,fi,Z) is the corresponding point in cylindrical coordinates
//
//       (R )         (Rs )
// J[0]= (fi) , J[1]= (fis)
//       (Z )         (Zs )
//
Vector3d *CStconfig::NJacobian(const Vector3d &u)
{
  NJacobian1Q(u);
  NJacobian2(INTERPOLATION_ORDER);
  ct().saveLastJacobian(mBnorm);
  return ct().Jmatr;
}

//****************************************************************************
// The same as above, but linear interpolation is used
Vector3d *CStconfig::NJacobianL(const Vector3d &u)
{
  NJacobian1L(u);
  NJacobian2(1);
  ct().saveLastJacobian(mBnorm);
  return ct().Jmatr;
}

//****************************************************************************
// 1st stage,  spline interpolation on r_eff
Vector3d *CStconfig::NJacobian1Q(const Vector3d &u)
{
  double s     = u[0];
  double phi   = u[2];
  Vector3d F (0.,0.,0.);          //  (R,fi,Z)
  Vector3d Fs(0.,0.,0.);          // d(R,fi,Z)/ds   -- partial derivative on s
  Vector3d Ft(0.,0.,0.);          // d(R,fi,Z)/dtheta
  Vector3d Fp(0.,0.,0.);          // d(R,fi,Z)/dphi

  //double iot = iota(s); //04jun2011
  bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);   //27Oct2011

  int M,N;
  int is = getSsegMN(s, M, N);
  ct().s  = s;  // save for NJacobian2
  ct().is = is; // save for NJacobian2
  const double * fCh = mCOEF.constArray()+mN*offsetN+ halfMesh.is0*offsetS;
  const double * fCf = mCOEF.constArray()+mN*offsetN+ fullMesh.is0*offsetS;
  const splineCoeffs * hMesh = &halfMesh;
  const splineCoeffs * RZmesh = vmecCoordinates? &fullMesh:&halfMesh;

  double csNp = tokamakEQDSK?0:(::cos(mNp*phi)); //if(tokamakEQDSK) phi=pi/2;
  double snNp = tokamakEQDSK?1:(::sin(mNp*phi));

  CComplex expNph(csNp,snNp);               // exp(1i*mNp*phi), where 1i is the sqrt(-1)
  CComplex exp_Nph = expNph.conj();         // exp(-1i*mNp*phi)
  CComplex expth (::cos(u[1]), ::sin(u[1]));  // exp(1i*theta)
  CComplex expmth(1,0);                     // will be exp(1i*m*theta)
  for(int m=0; m<=M; m++) {
    CComplex expmthNph =expmth; // exp(1i*m*theta+1i*mNp*phi), where 1i is the sqrt(-1)
    CComplex expmth_Nph=expmth; // exp(1i*m*theta-1i*mNp*phi), where 1i is the sqrt(-1)
    const double * fCoeff     = fCh + m*offsetM;
    const double * fCoeffFull = fCf + m*offsetM;
    const double * fCoeffPointer = vmecCoordinates?fCoeffFull:fCoeff;
    for(int n=0; n<= N; n++) {
      double iRmn, iZmn, iPmn=0;       // interpolated value of Rmn in point s
      double iRm_n,iZm_n,iPm_n=0;      // interpolated value of Rm,-n in point s
      double iRmns, iZmns, iPmns=0;    // interpolated value of dRmn/ds in point s
      double iRm_ns,iZm_ns,iPm_ns=0;   // interpolated value of dRm,-n/ds in point s
      int idx  =  n*offsetN;
      //int idx1 = idx + offsetS;
      //int id   = -n*offsetN;
      //int id1  = id  + offsetS;
      //here
      // bmn(is,  m,n)==fCoeff[idx]
      // bmn(is+1,m,n)==fCoeff[idx1]
      // bmn(is,  m,n)==mCOEF[(n+mN)*offsetN + m*offsetM + is*offsetS]
      // it looks ugly, but exactly here I need speed
      SPLINED(Rmn,  idx,iRmn, iRmns,RZmesh)
      SPLINED(Zmn,  idx,iZmn, iZmns,RZmesh)
      if(boozerCrd) {
        SPLINED(Phimn,idx,iPmn, iPmns,hMesh)
      }
      SPLINED(Rmn,  -idx ,iRm_n,iRm_ns,RZmesh)
      SPLINED(Zmn,  -idx ,iZm_n,iZm_ns,RZmesh)
      if(boozerCrd) {
        SPLINED(Phimn,-idx ,iPm_n,iPm_ns,hMesh)
      }
      if(n==0) {
        double csmt = expmth.re;  // cos(m*theta)
        double snmt = expmth.im;  // sin(m*theta)
        F [0] += iRmn*csmt;       //R
        F [1] += iPmn*snmt;       //fi
        F [2] += iZmn*snmt;       //Z
        Fs[0] += iRmns*csmt;     //Rs
        Fs[1] += iPmns*snmt;     //fis
        Fs[2] += iZmns*snmt;     //Zs
        Ft[0] -= m*iRmn*snmt;   //Rt
        Ft[1] += m*iPmn*csmt;   //fit
        Ft[2] += m*iZmn*csmt;   //Zt
        continue;
      }
      expmth_Nph *= exp_Nph;    // exp(1i*m*theta-1i*mNp*phi)
      expmthNph  *= expNph;     // exp(1i*m*theta+1i*mNp*phi)
      double sn1 = expmth_Nph.im;  //sin(m*theta-n*Np*phi) is stored into sn1
      double cs1 = expmth_Nph.re;  //cos(m*theta-n*Np*phi)
      double sn2 = expmthNph.im;   //sin(m*theta+n*Np*phi)
      double cs2 = expmthNph.re;   //cos(m*theta+n*Np*phi)

////if(fabs(m*iot-n*mNp)<1e-1) continue;  //04jun2011
////if(fabs(m*iot+n*mNp)<1e-1) continue;  //04jun2011

      F [0] += iRmn*cs1  + iRm_n*cs2;   //R
      Fs[0] += iRmns*cs1 + iRm_ns*cs2;  //Rs
      Ft[2] += m*(iZmn*cs1 + iZm_n*cs2);    //Zt
      Fp[2] -= n*(iZmn*cs1 - iZm_n*cs2);    //Zp
      if(boozerCrd) {
        Fp[1] -= n*(iPmn*cs1 - iPm_n*cs2);    //fip
        Ft[1] += m*(iPmn*cs1 + iPm_n*cs2);    //fit
        Fs[1] += iPmns*sn1 + iPm_ns*sn2;  //fis
        F [1] += iPmn*sn1  + iPm_n*sn2;   //fi
      }
      F [2] += iZmn*sn1  + iZm_n*sn2;   //Z
      Fs[2] += iZmns*sn1 + iZm_ns*sn2;  //Zs
      Fp[0] += n*(iRmn*sn1 - iRm_n*sn2);    //Rp
      Ft[0] -= m*(iRmn*sn1 + iRm_n*sn2);    //Rt
    }
    expmth *= expth;   // exp(1i*m*theta)
  }
  Fp *= mNp;    //d(R,fi,Z)/dphi
  F [1] *=-mPeriod; // cylindrical angle fi
  Fs[1] *=-mPeriod; // dfi/ds
  Ft[1] *=-mPeriod; // dfi/dtheta
  Fp[1] *=-mPeriod; // dfi/dphi
  F [1] += phi;     // cylindrical angle fi
  Fp[1] += 1;
#if INTERPOLATION_ORDER == 3
//  if(s==0) std::cerr<< "NJacobian1Q: s==0"<<std::endl;
  Fs *= halfMesh.e2rcurr;
#endif
  if(tokamakEQDSK) {
    Fp[0] = 0;
    Fp[1] = 1;
    Fp[2] = 0;
  }
// ct().Jmatr - Jacobian for Newton iteration
  ct().Jmatr[0]=F; // (R,fi,Z)
  ct().Jmatr[1]=Fs;// d(R,fi,Z)/ds   -- partial derivative on s
  ct().Jmatr[2]=Ft;// d(R,fi,Z)/dtheta
  ct().Jmatr[3]=Fp;// d(R,fi,Z)/dphi
  ct().Jmatr[4]=u;

  mMsum +=M;
  mNsum +=N;
  mNjac++;     // number of calls
  return ct().Jmatr;
}

//****************************************************************************
// 1st stage, linear interpolation
Vector3d *CStconfig::NJacobian1L(const Vector3d &u)
{
  double s     = u[0];
  double phi   = u[2];
  Vector3d F (0.,0.,0.);          //  (R,fi,Z)
  Vector3d Fs(0.,0.,0.);          // d(R,fi,Z)/ds   -- partial derivative on s
  Vector3d Ft(0.,0.,0.);          // d(R,fi,Z)/dtheta
  Vector3d Fp(0.,0.,0.);          // d(R,fi,Z)/dphi

  bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);   //27Oct2011

  int M,N;
  int is = getSsegMN(s, M, N);
  const double * fCh = mCOEF.constArray()+mN*offsetN+ halfMesh.is0*offsetS;
  const double * fCf = mCOEF.constArray()+mN*offsetN+ fullMesh.is0*offsetS;
  const splineCoeffs * hMesh = &halfMesh;
  const splineCoeffs * RZmesh = vmecCoordinates? &fullMesh:&halfMesh;

  ct().s  = s;  // save for NJacobian2
  ct().is = is; // save for NJacobian2

  double csNp = tokamakEQDSK?0:(::cos(mNp*phi)); //if(tokamakEQDSK) phi=pi/2;
  double snNp = tokamakEQDSK?1:(::sin(mNp*phi));

  CComplex expNph(csNp,snNp);               // exp(1i*mNp*phi), where 1i is the sqrt(-1)
  CComplex exp_Nph = expNph.conj();         // exp(-1i*mNp*phi)
  CComplex expth (::cos(u[1]), ::sin(u[1]));  // exp(1i*theta)
  CComplex expmth(1,0);                     // will be exp(1i*m*theta)
  for(int m=0; m<=M; m++) {
    CComplex expmthNph =expmth; // exp(1i*m*theta+1i*mNp*phi), where 1i is the sqrt(-1)
    CComplex expmth_Nph=expmth; // exp(1i*m*theta-1i*mNp*phi), where 1i is the sqrt(-1)
    const double * fCoeff     = fCh + m*offsetM;
    const double * fCoeffFull = fCf + m*offsetM;
    const double * fCoeffPointer = vmecCoordinates?fCoeffFull:fCoeff;
    for(int n=0; n<= N; n++) {
      double iRmn, iZmn, iPmn=0;       // interpolated value of Rmn in point s
      double iRm_n,iZm_n,iPm_n=0;      // interpolated value of Rm,-n in point s
      double iRmns, iZmns, iPmns=0;    // interpolated value of dRmn/ds in point s
      double iRm_ns,iZm_ns,iPm_ns=0;   // interpolated value of dRm,-n/ds in point s
      int idx  =  n*offsetN;
      //int idx1 = idx + offsetS;
      //int id   = -n*offsetN;
      //int id1  = id  + offsetS;
      //here
      // bmn(is,  m,n)==fCoeff[idx]
      // bmn(is+1,m,n)==fCoeff[idx1]
      // bmn(is,  m,n)==mCOEF[(n+mN)*offsetN + m*offsetM + is*offsetS]
      // it looks ugly, but exactly here I need speed
      LINTERPD(Rmn,  idx,iRmn, iRmns,RZmesh)
      LINTERPD(Zmn,  idx,iZmn, iZmns,RZmesh)
      LINTERPD(Rmn,  -idx,iRm_n,iRm_ns,RZmesh)
      LINTERPD(Zmn,  -idx,iZm_n,iZm_ns,RZmesh)
      if(boozerCrd) {
        LINTERPD(Phimn,idx,iPmn, iPmns,hMesh)
        LINTERPD(Phimn,-idx,iPm_n,iPm_ns,hMesh)
      }
      if(n==0) {
        double csmt = expmth.re;  // cos(m*theta)
        double snmt = expmth.im;  // sin(m*theta)
        F [0] += iRmn*csmt;       //R
        F [1] += iPmn*snmt;       //fi
        F [2] += iZmn*snmt;       //Z
        Fs[0] += iRmns*csmt;     //Rs
        Fs[1] += iPmns*snmt;     //fis
        Fs[2] += iZmns*snmt;     //Zs
        Ft[0] -= m*iRmn*snmt;   //Rt
        Ft[1] += m*iPmn*csmt;   //fit
        Ft[2] += m*iZmn*csmt;   //Zt
        continue;
      }
      expmth_Nph *= exp_Nph;    // exp(1i*m*theta-1i*mNp*phi)
      expmthNph  *= expNph;     // exp(1i*m*theta+1i*mNp*phi)
      double sn1 = expmth_Nph.im;  //sin(m*theta-n*Np*phi) is stored into sn1
      double cs1 = expmth_Nph.re;  //cos(m*theta-n*Np*phi)
      double sn2 = expmthNph.im;   //sin(m*theta+n*Np*phi)
      double cs2 = expmthNph.re;   //cos(m*theta+n*Np*phi)
      F [0] += iRmn*cs1  + iRm_n*cs2;   //R
      Fs[0] += iRmns*cs1 + iRm_ns*cs2;  //Rs
      Ft[2] += m*(iZmn*cs1 + iZm_n*cs2);    //Zt
      Fp[2] -= n*(iZmn*cs1 - iZm_n*cs2);    //Zp
      if(boozerCrd) {
        Fp[1] -= n*(iPmn*cs1 - iPm_n*cs2);    //fip
        Ft[1] += m*(iPmn*cs1 + iPm_n*cs2);    //fit  == dfi/dtheta
        Fs[1] += iPmns*sn1 + iPm_ns*sn2;  //fis
        F [1] += iPmn*sn1  + iPm_n*sn2;   //fi
      }
      F [2] += iZmn*sn1  + iZm_n*sn2;   //Z
      Fs[2] += iZmns*sn1 + iZm_ns*sn2;  //Zs
      Fp[0] += n*(iRmn*sn1 - iRm_n*sn2);    //Rp
      Ft[0] -= m*(iRmn*sn1 + iRm_n*sn2);    //Rt
    }
    expmth *= expth;   // exp(1i*m*theta)
  }
  Fp *=mNp;    //d(R,fi,Z)/dphi
  F [1] *=-mPeriod; // cylindrical angle fi
  Fs[1] *=-mPeriod; // dfi/ds
  Ft[1] *=-mPeriod; // dfi/dtheta
  Fp[1] *=-mPeriod; // dfi/dphi
  F [1] += phi;     // cylindrical angle fi
  Fp[1] += 1;
  if(tokamakEQDSK) {
    Fp[0] = 0;
    Fp[1] = 1;
    Fp[2] = 0;
  }
// ct().Jmatr - Jacobian for Newton iteration
  ct().Jmatr[0]=F; // (R,fi,Z)
  ct().Jmatr[1]=Fs;// d(R,fi,Z)/ds   -- partial derivative on s
  ct().Jmatr[2]=Ft;// d(R,fi,Z)/dtheta
  ct().Jmatr[3]=Fp;// d(R,fi,Z)/dphi
  ct().Jmatr[4]=u; // save magnetic coordinates

  mMsum +=M;
  mNsum +=N;
  mNjac++;     // number of calls
  return ct().Jmatr;
}


//****************************************************************************
// 2nd stage, calculate B, grad(s), metric tensor etc.
// linear interpolation of iota,Ipol,Itor if(typeOfInterp==1)
// spline interpolation of iota,Ipol,Itor if(typeOfInterp==3)
void CStconfig::NJacobian2(int typeOfInterp)
{
  double s=ct().s;
  int   is=ct().is;
  Vector3d F (ct().Jmatr[0]); // (R,fi,Z)
  Vector3d Fs(ct().Jmatr[1]);// d(R,fi,Z)/ds   -- partial derivative on s
  Vector3d Ft(ct().Jmatr[2]);// d(R,fi,Z)/dtheta
  Vector3d Fp(ct().Jmatr[3]);// d(R,fi,Z)/dphi

  Vector3d fluxCoord(ct().Jmatr[4]);// (s,theta,phi)

  double maxFlux = mFlux*mBnorm;
  double iota =   (typeOfInterp==1)?Linterp(miota.constArray()+is) : SPLINEy(miota);

  ct().modB    =  B(fluxCoord);

  Vector3d Bcontra(0.,0.,0.); 
  calcVmecLgB(fluxCoord, &Bcontra);

  double R  = F[0];
  if(vmecCoordinates) {
    ct().jac = ct().jacVmec/maxFlux;  // Jacobian in coordinates (Flux,theta,phi) -> ALPXP 21.08.2019 this does not work for s>1
  }
  else if(fpCoordinates) {
    ct().jac  = R*Fs*(Ft^Fp)/maxFlux; // Jacobian in coordinates (Flux,theta,phi)
  }
  else if(tokamakEQDSK) { // mIpol is defined in loadEfit.cpp
    double Ipol = ( (typeOfInterp==1)?Linterp(mIpol.constArray()+is) : SPLINEy(mIpol) ) * mBnorm;
    ct().jac = R*R/(mu0*Ipol);
  }
  else { // Jacobian for Boozer coordinates, the sign is important
    double Ipol = ( (typeOfInterp==1)?Linterp(mIpol.constArray()+is) : SPLINEy(mIpol) ) * mBnorm;
    double Itor = ( (typeOfInterp==1)?Linterp(mItor.constArray()+is) : SPLINEy(mItor) ) * mBnorm;
    ct().jac  =  mu0*(Ipol*mNp+iota*Itor)/(ct().modB*ct().modB*twopi*twopi);  // in (Flux,theta,phi) -> ALPXP 21.08.2019 this does not work for s>1
    ct().jacBoozer  = ct().jac*maxFlux;  // (s,theta,phi)
  }
// Multiplying by R we take into account the cylindrical coordinate system
  Fs[1] *=R; // R*dfi/ds
  Ft[1] *=R; // R*dfi/dtheta
  Fp[1] *=R; // R*dfi/dphi
// we must derive the correct sign of the maxFlux despite what was stored in the file
// next statesment is obsolete, 06 March 2008
// maxFlux = setSign(maxFlux,mSignJac_s*ct().jac);
  Vector3d Ff(Fs); // create a copy of Fs
  Ff /= maxFlux;   // from s to flux, Ff=d(R,fi,Z)/dflux -- partial derivative on flux
  double jacMixProd  = Ff*(Ft^Fp); // jacobian in (Flux,theta,phi)-coordinates,
  ct().jacMixProd = setSign(jacMixProd, ct().jac);  // set to jacMixProd the same sign as the sign of jac

  //double testJacobian1 = jacMixProd      - ct().jac;
  //double testJacobian2 = ct().jacMixProd - ct().jac;
  //std::cout<<testJacobian1<<" "<<testJacobian2<<std::endl;

  if(jacobianIsMixedProduct) ct().jac = ct().jacMixProd; // see comment -> ALPXP 21.08.2019 this does not work for s>1

  if((vmecCoordinates&&useVmecLambda)||fpCoordinates) {
    Bcontra.set(0., (iota-ct().lambdap)/(ct().jac*twopi), (1.+ct().lambdat)/(ct().jac*twopi) );
  }
  else {
    Bcontra.set(0., iota/(ct().jac*twopi), 1/(ct().jac*twopi) );
  }

  ct().Bcyl = Bcontra[1]*Ft + Bcontra[2]*Fp;
  ct().Bcontra = Bcontra;
  ct().Bcova.set(ct().Bcyl*Fs, ct().Bcyl*Ft, ct().Bcyl*Fp);
  ct().Bxyz.zero();

  //Bcyl = (iota*Ft   +Fp)/(jac*twopi);     // B  =(iota*dX /dtheta+dX /dphi)/J/twopi
  //Br  =  (iota*Ft[0]+Fp[0])/(jac*twopi);  // Br =(iota*dR /dtheta+dR /dphi)/J/twopi
  //Bfi =  (iota*Ft[1]+Fp[1])/(jac*twopi);  // Bfi=(iota*dfi/dtheta+dfi/dphi)*R/J/twopi
  //Bz  =  (iota*Ft[2]+Fp[2])/(jac*twopi);  // Bz =(iota*dZ /dtheta+dZ /dphi)/J/twopi
  //gradPsi_r  = R*(dfi/dtheta*dZ /dphi - dZ /dtheta*dfi/dphi)/jac
  //gradPsi_fi =   (dZ /dtheta*dR /dphi - dR /dtheta*dZ /dphi)/jac
  //gradPsi_z  = R*(dR /dtheta*dfi/dphi - dfi/dtheta*dR /dphi)/jac
  //gradPsi = (dX/dtheta^dX/dphi)/jac

  ct().gradtorFlux = (Ft^Fp)/ct().jac; // grad(Psi_tor) in cyl. coord, Psi_tor is the toroidal flux
  ct().gradS = ct().gradtorFlux;
  ct().gradS/= maxFlux; // grad(s) in cyl. coord
  ct().gradReff = ct().gradS;
  s = (s==0)?1e-8:s;
  ct().gradReff*= 0.5*ma/sqrt(s); // grad(reff) in cyl. coord

  ct().gradTheta = (Fp^Ff)/ct().jac; // grad(theta) in cyl. coord; theta is the poloidal angle
  
  if(tokamakEQDSK){
    ct().gradPhi.set(0.,1/R,0.);
  }
  else {
   ct().gradPhi = (Ff^Ft)/ct().jac; // grad(phi)   in cyl. coord; phi is the toroidal angle
  }

#if 0
  //////double jac_s = mSignJac_s*fabs(maxFlux*ct().jac);
  //////Vector3d gT = (Fp^Fs)/jac_s; // grad(theta) in cyl. coord, theta is the poloidal angle
  //////Vector3d gP = (Fs^Ft)/jac_s; // grad(phi)   in cyl. coord, phi is the toroidal angle

// metric tensor in coordinates (flux,theta,phi)
// g_ij = dX/du^i * dX/du^j , X - cartesian vector
//  (u^0,u^1,u^2) = (flux,theta,phi)
//g_00=Ff*Ff; g_01=Ff*Ft
//  Vector3d *g = ct().gTensor; in coordinates (flux,theta,phi)
  //////ct().gTensor[0].set(Ff*Ff,Ff*Ft,Ff*Fp);
  //////ct().gTensor[1].set(Ff*Ft,Ft*Ft,Ft*Fp);
  //////ct().gTensor[2].set(Ff*Fp,Ft*Fp,Fp*Fp);
  //////ct().detOfgTensor = ct().gTensor[0]*(ct().gTensor[1]^ct().gTensor[2]);  // det(g_ij)
  //////if(ct().detOfgTensor<0) std::cerr<<"why |g_ij|=-1 ";

//test  double jac3 = sqrt(fabs(ct().detOfgTensor));  // sqrt(det(g_ij)) must be equal to jacobian jacMixProd
//test  if(jac<0) jac3 = -jac3;
//test res =1e16*(jacMixProd - jac3);

//test  double modgradPsi = sqrt( fabs(ct().gTensor[1][1]*ct().gTensor[2][2]-ct().gTensor[1][2]*ct().gTensor[1][2] )/(jac*jac));
//test  modgradPsi  must be equal to |ct().gradtorFlux|
//test  res =1e16*(ct().gradtorFlux.abs() - modgradPsi);
#endif
}

//****************************************************************************
void CStconfig::zeroStatistics()
{
  mMsum=mNsum=mNjac=0;
  mNcyl2mag=mNmag2cyl=0;        // number of calls
  mMsum2=mNsum2=0;
  mNGuessCalls = 0;  
}

//****************************************************************************
double CStconfig::MAverage() const 
{
  double i = mMsum2/double(mNmag2cyl?mNmag2cyl:1);
  double j = mMsum/double(mNjac?mNjac:1);
  return mmax(i,j);
}

//****************************************************************************
double CStconfig::NAverage() const 
{
  double i = mNsum2/double(mNmag2cyl?mNmag2cyl:1);
  double j = mNsum/double(mNjac?mNjac:1);
  return mmax(i,j);
}

double CStconfig::NNjacCalls() const { return mNjac;}
double CStconfig::Ncyl2magCalls() const { return mNcyl2mag;}
double CStconfig::Nmag2cylCalls() const { return mNmag2cyl;}
double CStconfig::NGuessCalls() const { return mNGuessCalls;}

//BEGIN_methods_to retrieve_results*******************************************
//****************************************************************************

#define IS_CYL_VALID(cyl)  if( !ct().isCylValid(cyl,mBnorm) ) cyl2mag(cyl);
#define IS_XYZ_VALID(xyz)  if( !ct().isXyzValid(xyz,mBnorm) ) xyz2mag(xyz);

//****************************************************************************
double CStconfig::getmodB(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl); 
  return getmodB(); 
}

double CStconfig::getmodB() const 
{ 
  return ct().modB; 
}

//****************************************************************************
// Jac(torFlux,theta,phi) ~ mu0/(4*pi^2)*Ipol/B^2
double CStconfig::getJac(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl); 
  return getJac(); 
}

double CStconfig::getJac() const 
{ 
  return ct().jac; 
}

//****************************************************************************
// Jacobian in coordinates (s,theta,phi)
// J(s,theta,phi) ~ Flux(a)*mu0*Ipol/(4*pi^2*B^2) for Boozer coord
double CStconfig::getJ(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl);  
  return getJ(); 
}

double CStconfig::getJ() const 
{ 
  return mSignJac_s*fabs(Flux(1.)*ct().jac); 
}

//****************************************************************************
// Jacobian in coordinates (torFlux,theta,phi)
// Jac(torFlux,theta,phi) = r/mFlux*(Rs*(Rt^Rp)) -- jacobian, where ^ is the cross product sign,
//   where Rs = dR/ds - partial derivative
double CStconfig::getJac2(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl);  
  return getJac2(); 
}

double CStconfig::getJac2() const
{
  return ct().jacMixProd; 
}

//****************************************************************************
// J(s,theta,phi)=r*(Rs*(Rt^Rp)) Jacobian in coordinates (s,theta,phi)
double CStconfig::getJ2(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl);  
  return getJ2(); 
}

double CStconfig::getJ2() const 
{ 
  return mSignJac_s*fabs(Flux(1.)*ct().jacMixProd); 
}

//****************************************************************************
const Vector3d & CStconfig::getBcyl(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl); 
  return getBcyl(); 
}

const Vector3d & CStconfig::getBcyl() const 
{ 
  return ct().Bcyl; 
}

//****************************************************************************
// @return grad(|B|) in cylindrical coordinates at point cyl
Vector3d  CStconfig::getGradBcyl(const Vector3d &cyl)
{
  IS_CYL_VALID(cyl); 
  return getGradBcyl();
}

//****************************************************************************
// @return grad(|B|) in cylindrical coordinates after coordinate transformation
Vector3d  CStconfig::getGradBcyl()
{
  Vector3d u = ct().magLast;
  Vector3d gr = getCovaGradB(u);
  Vector3d gradB = gr[0]*getGrads() + gr[1]*getGradtheta() + gr[2]*getGradphi();
  return gradB;
}

//****************************************************************************
// @return grad(|B|) in cartesian coordinates at point xyz
Vector3d  CStconfig::getGradBxyz(const Vector3d &xyz)
{
  IS_XYZ_VALID(xyz);
  return getGradBxyz();
}

//****************************************************************************
// @return grad(|B|) in cartesian coordinates after coordinate transformation
Vector3d  CStconfig::getGradBxyz()
{
  Vector3d gs = getGradBcyl();
  Vector3d c  = ct().cylLast;
  return gs.cylVector2Cartesian(c);
}


//****************************************************************************
const Vector3d & CStconfig::getBcova() const 
{ 
  return ct().Bcova; 
}

//****************************************************************************
const Vector3d & CStconfig::getBcontra() const 
{ 
  return ct().Bcontra; 
}

//****************************************************************************
void CStconfig::getCovaBasis(Vector3d &e_s, Vector3d &e_t, Vector3d &e_p) const
{
  e_s = ct().Jmatr[1];
  e_t = ct().Jmatr[2];
  e_p = ct().Jmatr[3];
    // Multiplying by R we take into account the cylindrical coordinate system
  double R = ct().Jmatr[0][0];
  e_s[1] *=R; // R*dfi/ds
  e_t[1] *=R; // R*dfi/dtheta
  e_p[1] *=R; // R*dfi/dphi
}

//****************************************************************************
void CStconfig::getContraBasis(Vector3d &e1, Vector3d &e2, Vector3d &e3) const
{
  e1 = ct().gradS;
  e2 = ct().gradTheta;
  e3 = ct().gradPhi;
}

//****************************************************************************
//in cartesian coordinates
const Vector3d & CStconfig::getBxyz(const Vector3d &xyz) 
{ 
  IS_XYZ_VALID(xyz);
  return getBxyz();
}

const Vector3d & CStconfig::getBxyz() 
{
  if(ct().Bxyz==Vector3d(0.,0.,0.)) {
    double fi = ct().Jmatr[0][1];
    double cs = ::cos(fi);
    double sn = ::sin(fi);
    double Bx=ct().Bcyl[0]*cs-ct().Bcyl[1]*sn;
    double By=ct().Bcyl[0]*sn+ct().Bcyl[1]*cs;
    double Bz=ct().Bcyl[2];
    (ct().Bxyz).set(Bx,By,Bz);
  }
  return ct().Bxyz; 
}

//****************************************************************************
//in cyl. coordinates
const Vector3d & CStconfig::getGradFlux(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl);  
  return ct().gradtorFlux; 
}

const Vector3d & CStconfig::getGradFlux() const 
{ 
  return ct().gradtorFlux; 
}

//****************************************************************************
//in cyl. coord
const Vector3d & CStconfig::getGrads(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl);  
  return ct().gradS; 
}

const Vector3d & CStconfig::getGrads() const   
{ 
  return ct().gradS; 
}

//****************************************************************************
//in cyl. coord
const Vector3d & CStconfig::getGradreff(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl);  
  return ct().gradReff; 
}

const Vector3d & CStconfig::getGradreff() const 
{ 
  return ct().gradReff; 
}

//****************************************************************************
//in cyl. coord
const Vector3d & CStconfig::getGradtheta(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl);  
  return ct().gradTheta; 
}

const Vector3d & CStconfig::getGradtheta() const   
{ 
  return ct().gradTheta; 
}

//****************************************************************************
// @return grad(lambda) in cylindrical coordinates
Vector3d CStconfig::getGradlambda(const Vector3d &cyl)
{
  IS_CYL_VALID(cyl);
  Vector3d gradLambda = ct().lambdas*ct().gradS + ct().lambdat*ct().gradTheta + ct().lambdap*ct().gradPhi;
  return gradLambda; 
}

Vector3d CStconfig::getGradlambda() const
{
  Vector3d gradLambda = ct().lambdas*ct().gradS + ct().lambdat*ct().gradTheta + ct().lambdap*ct().gradPhi;
  return gradLambda; 
}

//****************************************************************************
//in cyl. coord
const Vector3d & CStconfig::getGradphi(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl);  
  return ct().gradPhi; 
}

const Vector3d & CStconfig::getGradphi() const   
{ 
  return ct().gradPhi; 
}

//****************************************************************************
// transform to Cartesian coordinates
Vector3d  CStconfig::getGradsxyz(const Vector3d &xyz) 
{ 
  IS_XYZ_VALID(xyz);
  return getGradsxyz();
}

Vector3d  CStconfig::getGradsxyz() const
{
  Vector3d gs = getGrads();
  Vector3d c  = ct().cylLast;
  return gs.cylVector2Cartesian(c);
}

//****************************************************************************
const Vector3d *CStconfig::getNJmatrix(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl); 
  return getNJmatrix(); 
}

const Vector3d *CStconfig::getNJmatrix() 
{
  return ct().Jmatr; 
}

//****************************************************************************
const Vector3d *CStconfig::getMetricTensor(const Vector3d &cyl) 
{ 
  IS_CYL_VALID(cyl); 
  return getMetricTensor(); 
}

const Vector3d *CStconfig::getMetricTensor() 
{
  Vector3d Fs(ct().Jmatr[1]);// d(R,fi,Z)/ds   -- partial derivative on s
  Vector3d Ft(ct().Jmatr[2]);// d(R,fi,Z)/dtheta
  Vector3d Fp(ct().Jmatr[3]);// d(R,fi,Z)/dphi
  // Multiplying by R we take into account the cylindrical coordinate system
  double R = ct().Jmatr[0][0]; 
  Fs[1] *=R; // R*dfi/ds
  Ft[1] *=R; // R*dfi/dtheta
  Fp[1] *=R; // R*dfi/dphi
  // metric tensor in coordinates (s,theta,phi)
  // g_ij = dX/du^i * dX/du^j , X - cartesian vector
  //  (u^0,u^1,u^2) = (s,theta,phi)
  //g_00=Fs*Fs; g_01=Ff*Ft
  //  Vector3d *g = ct().gsTensor; in coordinates (s,theta,phi)
  ct().gsTensor[0].set(Fs*Fs,Fs*Ft,Fs*Fp);
  ct().gsTensor[1].set(Fs*Ft,Ft*Ft,Ft*Fp);
  ct().gsTensor[2].set(Fs*Fp,Ft*Fp,Fp*Fp);
  return ct().gsTensor; 
}

//****************************************************************************
// get the cylindrical coordinates of the last point used in coordinate transformation
Vector3d  CStconfig::getRcyl() const 
{ 
  return ct().Jmatr[0];// ==(R,fi,Z)
} 

//****************************************************************************
// get the cartesian coordinates of the last point used in coordinate transformation
Vector3d  CStconfig::getRxyz() const   
{ 
  Vector3d r = ct().Jmatr[0]; // ct().Jmatr[0]=(R,fi,Z)
  return r.toCartesian();     // cyl-->xyz transformation
}

//END_methods_to retrieve_results*********************************************



//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void CStconfig::epsEffInvalidate() 
{
  epsilonEff.invalidate();
  epsilonEff.clear();
}

void CStconfig::epsEffSetXeffParam(double xmin, double xmax, int size) 
{
  epsilonEff.maxP = size;
  epsilonEff.xmin = xmin;
  epsilonEff.xmax = xmax;
}

void CStconfig::epsEffSetTracingParam(int turns, bool doAccuracyTest, double dphi)
{
  epsilonEff.turns= turns;
  epsilonEff.dphi = dphi;     // twopi/512 == 0.0122718
  epsilonEff.doAccuracyTest = doAccuracyTest;
}

void CStconfig::epsEffSetIotaParam(bool avoidResonances, bool useRationalIota, int iotaDenominator) 
{
  epsilonEff.useRationalIota=useRationalIota;
  epsilonEff.iotaDenominator=iotaDenominator;
  epsilonEff.avoidResonances=avoidResonances;
}

void CStconfig::epsEffSetMagnMomentParam(int numBLevels)
{
  epsilonEff.numBLevels = numBLevels;
}



void CStconfig::epsEffCreate() 
{
  if(epsilonEff.useRationalIota) {
    epsilonEff.avoidResonances=false;
    epsilonEff.doAccuracyTest=false;
  }
  epsilonEff.create(this); 
}

double CStconfig::epsEffGraz(double s)  
{
  epsilonEff.create(this); 
  return epsilonEff.Graz(s);
}

double CStconfig::epsEffDkes(double s) 
{
  epsilonEff.create(this); 
  return epsilonEff.Dkes(s);
}

double CStconfig::epsEffVmec(double s) 
{
  epsilonEff.create(this); 
  return epsilonEff.Vmec(s);
}

double CStconfig::epsEffNorm(double s) 
{
  epsilonEff.create(this); 
  return epsilonEff.Norm(s);
}

double CStconfig::epsEffHat(double s) 
{
  epsilonEff.create(this); 
  return epsilonEff.Hat(s);
}

double CStconfig::Gw(double s) 
{
  epsilonEff.create(this); 
  return epsilonEff.Gw(s);
}

double CStconfig::Gs(double s) 
{
  epsilonEff.create(this); 
  return epsilonEff.Gs(s);
}

double CStconfig::Gv(double s) 
{
  epsilonEff.create(this); 
  return epsilonEff.Gv(s);
}


//****************************************************************************
void CStconfig::calculateEpsEff()
{
  if(!mOK) return;
  int size = epsilonEff.maxP;
  Threads2<CStconfig,CStconfig> threads(this, &CStconfig::calculateEpsEffExe,this,"calculateEpsEff");  
  int low = 0, upper = low+size-1;
  threads.run(low, upper, nThreads); 
}

//****************************************************************************
// This function is called from Threads to allocate
// memory for arrays and  initialize them, for example
//    x.resize(size);
//    y.resize(size);
//    double dx = 1./(size-1);
//    for(int i=0; i<size; i++) 
//      x[i] = dx*i; 
//****************************************************************************
// This function is called from theads to fill arrays allocated 
// by the setup-function calculateEpsEffSetup();
//   i1 and i2 range of array indexes to be filled;
//   use them as follows:
//      CStconfig *mc = new CStconfig;
//      *mc = *this;  // need a new copy of CStconfig to be thread safe
//      for(int i=i1; i<i2; i++) 
//         this->y[i] = mc->epsEff(x[i]*x[i]);
void CStconfig::calculateEpsEffExe(CStconfig * mc, int low, int upper, bool setup)
{
  if(setup) {  // setup part; it allocates arrays and initialize them
    int maxP = upper-low+1;
    epsilonEff.resize(); // detach data if they are shared 
    if(epsilonEff.x.empty()) return;
    if(epsilonEff.x.size()!=maxP) return;
    for(int i=0; i<maxP; i++) {
      double dx = (epsilonEff.xmax-epsilonEff.xmin)/(maxP-1);
      double x  = epsilonEff.xmin + dx*i;  // x=sqrt(s)
      if(epsilonEff.avoidResonances) {
        if(isRationalIota(x*x)) x += dx*0.3;  // try to avoid rational iota
        if(isRationalIota(x*x)) x += dx*0.1;  // try to avoid rational iota
      }
      epsilonEff.x.rw()[i]  = x;   // x=sqrt(s)
    }
    double Bmi,Bma;
    getBminmax(1,Bmi,Bma);
    epsilonEff.Bmin = Bmi*0.9;
    epsilonEff.Bmax = Bma*1.1;
    epsilonEff.Btruncation = truncate(); // save truncation level of B
    return;
  }
  else 
  {  // exe part; it will be called from threads to fill arrays
    if(epsilonEff.x.empty()) return;
    __allEps eps;
//    int turns = mc->epsilonEff.turns;
    for(int i=low; i<=upper; i++) {
      double x = mc->epsilonEff.x[i];
      double s = x*x;
      mc->epsEff(s,eps);
      epsilonEff.epsGraz.w(i) = eps.Graz;
      epsilonEff.epsDkes.w(i) = eps.Dkes;
      epsilonEff.epsVmec.w(i) = eps.Vmec;
      epsilonEff.epsNorm.w(i) = eps.Norm;
      epsilonEff.epsHat.w(i)  = eps.Hat;
      epsilonEff.GwGraz.w(i)  = eps.GwGraz;
      epsilonEff.GsGraz.w(i)  = eps.GsGraz;
      epsilonEff.GvGraz.w(i)  = eps.GvGraz;

      ////double grdsavrg = mc->getGradsAverages(s)[2];
      ////epsilonEff.epsHat.w(i)  = grdsavrg;
    }
  }
}

//****************************************************************************
bool CStconfig::isRationalIota(double s, double eps)
{
//04jun2011
  double iot = iota(s);
  int M,N;
  int is = getSsegMN(s, M, N);
  for(int m=2; m<=M; m++) {
    for(int n=2; n<=N; n++) {
      double e = fabs(m*iot-n*mNp);
      if(e<eps) {
        TESTPRINT(std::cerr<<"abs(m*iota-n)="<<e<<"<"<<eps<<"; iota="<<iot<<"("<<n*mNp<<"/"<<m<<"); r/a="<<sqrt(s)<<std::endl<<std::flush);
        return true;  
      }
      e = fabs(m*iot+n*mNp);
      if(e<eps) {
        TESTPRINT(std::cerr<<"abs(m*iota-n)="<<e<<"<"<<eps<<"; iota="<<iot<<"("<<n*mNp<<"/"<<m<<"); r/a="<<sqrt(s)<<std::endl<<std::flush);
        return true;  
      }
    }
  }
  return false;
#if 0
  for(int m=1; m<20; ++m) {
    double im = fabs(iota(s)*m);
    double e = mantissa(im);
    if(1-e<0.005) return true;
    if(int(im)&&e<0.005) return true;
  }
  return false;
#endif
}

//****************************************************************************
// int numBlevels = 513;        // # of b levels
// int numTurns = 1000;         // # of full turns of Bfield line around torus
void CStconfig::epsEff(double s, __allEps &eps)
{
  volatile ctStack::saveState state(ct);

  const double relEps=0.01;   // if( |resultBsPrev-resultBs| < |resultBs|*relEps ) ++accuracyCnt;
                              //    if(accuracyCnt>10) break;

  bool BOOZER_COORDINATES = epsilonEff.useBoozerCoord; // if true then use properties of boozer coordinates
  bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);
  if(!boozerCrd) BOOZER_COORDINATES = false;

  double dphi = epsilonEff.dphi;
  dphi = dphi==0?(0.4*degree):dphi;   // toroidal-angle increment.

  int numBlevels = epsilonEff.numBLevels+1;
  int numTurns   = epsilonEff.turns;

  double iotaRatnl=iota(s);
  int num=epsilonEff.iotaDenominator;
  int den=epsilonEff.iotaDenominator;

  if(epsilonEff.useRationalIota) {
    //double iotaR = cfraction(iota(s), num,den);                                                                
    double iotaR = iotaRational(s,epsilonEff.iotaDenominator,num,den);
    iotaRatnl = iotaR!=0?iotaR:iotaRatnl;
    numTurns = den;
  }

  REAL *b = new REAL[6*numBlevels]; // array of b levels
  REAL *eb =  b + numBlevels;
  REAL *H  = eb + numBlevels;
  REAL *I  =  H + numBlevels;
  REAL *dgdb  = I + numBlevels;
  REAL *dIdb  = dgdb + numBlevels;
  char *flag = new char[numBlevels];

  REAL R0 = (REAL)this->R0();    // R0 is the normalization factor, take the major radius from the input file 
  //REAL R0 = (REAL)Rmn_sp(1,0,0); // R0 is the normalization factor
  REAL B00  = (REAL)Bmn_sp(s,0,0);
  REAL B00p = (REAL)Bmn_sp_prime(s,0,0); // d B_00/ds
  REAL B0 = B00;
  ////double Bmi,Bma; getBminmax(s,Bmi,Bma); Bmi *= 0.9/B0; Bma *= 1.1/B0;

  double bmi = epsilonEff.Bmin/B0;
  double bma = epsilonEff.Bmax/B0;
  double Bmaxabs;

{
Vector3d magcoordBmin(s,0.,0.), magcoordBmax(s,0.,0.);
findBminBmaxPosition(-s, magcoordBmin, magcoordBmax);
Bmaxabs = B(magcoordBmax);
}

  double db = (bma - bmi)/(numBlevels-1);

  for(int k=0; k<numBlevels; k++) {
    b[k] = REAL(bmi + k*db);   // b' in original paper
    eb[k] = 1/b[k];
    H[k] = 0;
    I[k] = 0;
    dgdb[k] = 0;
    dIdb[k] = 0;
    flag[k] = 0; // 0-undefined; 1-do H,I integrals
  }

  Vector3d mag0(s,pi,0); 
  Vector3d mag1(mag0),r1, r0, B1, gradS1;
  REAL Bmod,gradSmod,kG,dl;
  //1REAL Bprev;

  advanceDataForEpsEff02(mag1, r1, B1, gradS1, Bmod, gradSmod, kG, dl,0, 0,iotaRatnl);  //initialize
 
  r0 = r1;    // save starting position
  Vector3d rprev = r1;
  //1Bprev = Bmod;

  REAL Int(0), Norm(0), gradSAverage(0),gradS2overB2(0);
  REAL result(0), resultPrev(0);
  REAL IntGw(0);           // main integral in \Gamma_w expression
  REAL resultGw(0), resultGwPrev(0);
  REAL gradSsqrtAverage(0);

  int nsteps = int(twopi/dphi);
  dphi = twopi/nsteps;
  int accuracyCnt(0);
  int turnsCnt(0);
  //1double dphi0 = dphi;
  for(int i=0 ; i<numTurns; i++) {
    ++turnsCnt;
    for(int j=0 ; j<nsteps; j++) {
      REAL B,dl_B,gradSkG;
      double BcontraPhi;
      if(BOOZER_COORDINATES) {
        mag1 = magCoordFieldLine(mag1, mag1[2]+dphi,iotaRatnl); // advance along field line
        Vector3d gr = getCovaGradB(mag1);
        NJacobian(mag1);
        r1     = getRxyz();
        dl = (rprev-r1).abs(); 
        rprev=r1;
        Vector3d grads  = getGrads();
        Vector3d gradTh = getGradtheta();
        Vector3d gradPh = getGradphi();
        Vector3d Bvec   = getBcyl();
        B = REAL(getmodB());
        gradSmod = grads.abs();
        BcontraPhi = getBcontra()[2];
        dl_B = fabs(dphi/BcontraPhi);   //==dphi/B^phi is the same as dl/B;
        dl_B = dl/B;
       //2107 Vector3d gradB  = gr[1]*gradTh + gr[2]*gradPh + gr[0]*grads;
       //2107 REAL gradSkG = ((Bvec^gradB)*grads)/(B*B);  // the same as gradSmod*kG;
        Vector3d gradBpar  = gradTh;
        gradBpar  *= gr[1];
        gradBpar  += gr[2]*gradPh;
        gradSkG = -((Bvec^grads)*gradBpar)/(B*B);  // must be the same as gradSmod*kG;
      }
      else {
  #define  USE_KG 1    // use geodesic curvature if true
  #if USE_KG
        advanceDataForEpsEff02(mag1, r1, B1, gradS1, Bmod, gradSmod, kG, dl,0, dphi,iotaRatnl);
        //1dphi = dphi0*(1+fabs(B-Bprev)/(Bmaxabs*dphi*R0));
        B = Bmod;
        BcontraPhi = getBcontra()[2];   // == (Bvec*gradPh)
        dl_B = fabs(dphi/BcontraPhi);   //==dphi/B^phi is the same as dl/B;
        gradSkG = gradSmod*kG; 
  #else
        mag1 = magCoordFieldLine(mag1, mag1[2]+dphi,iotaRatnl); // advance along field line
        Vector3d gr = getCovaGradB(mag1);
        NJacobian(mag1);
        r1     = getRxyz();
        dl = (rprev-r1).abs(); 
        rprev=r1;
        Vector3d grads  = getGrads();
        Vector3d gradTh = getGradtheta();
        Vector3d gradPh = getGradphi();
        Vector3d Bvec   = getBcyl();
        B = REAL(getmodB());
        gradSmod = grads.abs();
        BcontraPhi = getBcontra()[2];
        dl_B = fabs(dphi/BcontraPhi);   //==dphi/B^phi is the same as dl/B;
       //2107 Vector3d gradB  = gr[1]*gradTh + gr[2]*gradPh + gr[0]*grads;
       //2107 REAL gradSkG = ((Bvec^gradB)*grads)/(B*B);  // the same as gradSmod*kG;
        Vector3d gradBpar  = gradTh;
        gradBpar  *= gr[1];
        gradBpar  += gr[2]*gradPh;
        gradSkG = -((Bvec^grads)*gradBpar)/(B*B);  // must be the same as gradSmod*kG;
  #endif
      }
#undef USE_KG
      Norm += dl_B;
      gradSAverage += dl_B*gradSmod;
      gradS2overB2 += dl_B*square(gradSmod/B);
      gradSsqrtAverage += dl_B*gradSmod*sqrt( B>=Bmaxabs?0:(1-B/Bmaxabs) );
      
      REAL B_B0 = B/B0;      
      REAL B04_B = 4/B_B0; //REAL B04_B = 4*B0/B;

      for(int k=0; k<numBlevels; k++) {
        REAL f = 1 - B_B0*eb[k]; //REAL f = 1 - B/(B0*b[k]);
        if(f>0) {
          if(flag[k]==0) continue;
          if(flag[k]==2) flag[k] = 1; 
          if(flag[k]==1) {
            REAL rf = sqrt(f);
            REAL drf = B_B0*square(eb[k])/(2*rf);
            REAL u  = B04_B-eb[k];
            REAL du = square(eb[k]);
            H[k] += dl_B*gradSkG * rf*sqrt(b[k])*u;
            I[k] += dl_B * rf;
            dgdb[k] += dl_B*gradSkG * (drf*u+rf*du);
            dIdb[k] += dl_B * drf;
          }
        }
        else {
          if(flag[k]==0) flag[k] = 2;
          if(flag[k]==2) continue;
          if(flag[k]==1) {
            flag[k] = 2;
            if(I[k]!=0)    Int  += square(H[k]*eb[k])/I[k];    // integral over b'
            if(dIdb[k]!=0) IntGw+= square(dgdb[k]/3 )/dIdb[k]; // integral over b'
            H[k] = 0;
            I[k] = 0;
            dgdb[k] = 0;
            dIdb[k] = 0;
          }
        }
      }

    }
    result = Int/Norm;
    resultGw = IntGw/Norm;

    Vector3d magerr = mag1-mag0;
    magerr[1] = mod2pi(magerr[1]+pi)-pi;
    magerr[2] = mod2pi(magerr[2]+pi)-pi;
    if(magerr.abs()<1e-4&&turnsCnt>10) {
      break;
      std::cerr.precision(3);
      std::cerr.setf(std::ios::scientific|std::ios::showpos);
      std::cerr<<"periodicity="<<magerr.abs()<<" x="<<sqrt(s);
      if(epsilonEff.useRationalIota) {
        std::cerr<<" iota="<<num<<"/"<<den
                 <<" iota_err="<< iota(s)-double(num)/double(den);
      }
      else {
        std::cerr<<" iota="<<iota(s);
      }
      std::cerr<<" turns="<<turnsCnt<<std::endl<<std::flush;
      break;
    }

    ////double dr = (r0-r1).abs();
    ////dr /= this->r(s);
    ////if(turnsCnt>20&&dr<1*degree) break;
    if(!epsilonEff.useRationalIota && epsilonEff.doAccuracyTest) {
      if(fabs(resultPrev-result) <= result*relEps) ++accuracyCnt;
      if(accuracyCnt>5&&turnsCnt>30) break; //if(accuracyCnt>=4) break;
    }
    resultPrev   = result;
    resultGwPrev = resultGw;
    //1Bprev = Bmod;
  }
  delete b;
  delete flag;
  //std::cerr<<"epsEff::s="<<s<<" turns="<<turnsCnt<<std::endl;
  gradSAverage /= Norm;  // < |grad(s)| >
  gradS2overB2 /= Norm;  // < (grad(s)/B)^2 >
  gradSsqrtAverage /= Norm;

  double epsNormFactor = 1/(square(B0)*gradS2overB2);  // 1/(B0^2 <(grad(s)/B)^2 > )
  double eps32hat  = pi*square(R0)*db*result/(8*sqrt(2.)); //  pi*square(R0)*db*Int/(8*sqrt(2.)*Norm);    
  double eps32Graz = eps32hat/square(gradSAverage);
  double eps32norm = eps32hat*epsNormFactor; // see eq.(9) in 28th EPS, http://epsppd.epfl.ch/Madeira/html/pdf/P5.055.pdf 
  double vmecFactor = square(gradSAverage)*square(r(s)/(2*s));   // (dr/dr_graz)^2 = <|grad(s)|>^2 * (r/(2s))^2, where r^2 = a^2*s, dr_graz = ds/<|grad(s)|>
  double r_dk2 = fabs(Flux(s)/(pi*B00));     
  double dkesFactor = square(gradSAverage)*r_dk2/(4*s*s);        // (dr/dr_graz)^2 = <|grad(s)|>^2 * (r/(2s))^2, where r^2 = Flux(s)/pi/B_00(s), dr_graz = ds/<|grad(s)|>
  dkesFactor *= square(1-s*B00p/B00);

  double eps32dkes = eps32Graz * dkesFactor;
  double eps32vmec = eps32Graz * vmecFactor;

  double gradFactorVmec = square(r(s)/(2*s)); // grad(r)^2 = (r/(2s))^2 * grad(s)^2, where r^2 = a^2*s
  double gradFactorDkes = r_dk2/(4*s*s);      // grad(r)^2 = (dr/ds)^2 * grad(s)^2=(r/(2s)*(1-s*B00p/B00))^2 * grad(s)^2, where r^2 = Flux(s)/pi/B_00(s)
  gradFactorDkes *= square(1-s*B00p/B00);

      //eps32Graz /= square(R0);
      //R0 = 8; //  R0 = Volume/(2(pi*a_graz)^2) for w7x-sc1
      //eps32Graz *= square(R0);

  double GwGraz =  pi*square(R0)/sqrt(2.)*db*resultGw/square(gradSAverage);  // \Gamma_w  eq(21) in Phys.Plasmas 12,112507(2005) http://link.aip.org/link/doi/10.1063/1.2131848
  double GsGraz =  pi/(2*sqrt(2.)) * gradSsqrtAverage/gradSAverage;
  double GvGraz =  sqrt(GwGraz/GsGraz);

  eps.GwGraz = GwGraz;
  eps.GsGraz = GsGraz;
  eps.GvGraz = GvGraz;

  eps.Graz = pow(eps32Graz,2./3); 
  eps.Dkes = pow(eps32dkes,2./3); 
  eps.Vmec = pow(eps32vmec,2./3);
  eps.Norm = pow(eps32norm,2./3);
      //////eps.Hat  = pow(eps32hat*gradFactorVmec,2./3);  // why i did this? for testing
      //////eps.Hat  = pow(eps32hat*gradFactorDkes,2./3);  // why i did this? for testing
  eps.Hat  = pow(eps32hat,2./3);
  eps.gradSmodAverage = gradSAverage;
  eps.vmecFactor = vmecFactor;
  eps.dkesFactor = dkesFactor;
  eps.xi = square(gradSAverage)*epsNormFactor; // xi = <|grad(s)|>^2 / (B0^2 <(grad(s)/B)^2 > ); see eq.(8) in 28th EPS, http://epsppd.epfl.ch/Madeira/html/pdf/P5.055.pdf  

  std::cerr.setf(std::ios::fixed,std::ios::showpos|std::ios::scientific);
}

//****************************************************************************
// all parameters are advanced along B-field line, except dphi and iotaRatnl
//
void CStconfig::advanceDataForEpsEff02(Vector3d &mag, Vector3d &r1, Vector3d &B1, Vector3d &n1, 
                                       REAL &Bmod, REAL &grads, REAL &kG, REAL &dl, Vector3d *curvature, 
                                       const double &dphi, const double &iotaRatnl)
{
//see epsEff()  volatile ctStack::saveState state(ct);
 // mag is the vector with the components (s, theta, phi); where theta is the poloidal angle
 // advance along field line from (s, theta, phi)  to (s, thetaNew, phi+dphi)
 // for case of Boozer coordinates: (s, theta, phi)  -->  (s, theta+iota*dphi, phi+dphi)
  mag = magCoordFieldLine(mag, mag[2]+dphi,iotaRatnl); // advance along field line
 // mag here is advanced
  NJacobian(mag);             // calculate all quantities at point mag
  Vector3d r2 = getRxyz();    // get position in cartesian coordinates
  Vector3d B2 = getBxyz();    // get B in cartesian coordinates
  Vector3d n2 = getGradsxyz(); // gradient of normalized toroidal flux
  Bmod = REAL(getmodB());      // |B|

  Vector3d b1(B1);
  Vector3d b2(B2);
  grads = REAL(n2.abs());
  Vector3d dr;
  if(dphi!=0) {
    dr = r2 - r1;
    dl = REAL(dr.abs());
  }
  else { 
    dl = 0;
  }
  r1 = r2;
  B1 = B2;
  n1 = n2;
  if(dphi==0) return;

  // Calculate the geodesic curvature, 
  // its sign is not important for eps_eff calculation  
  b1.norm();
  b2.norm();
  b2 -= b1;
  Vector3d k = b2/dl;  //  db/dl is the curvature  
  Vector3d g = n2^B2;  // ^ is the vector product
  g.norm();
  kG = REAL(k*g);    // geodesic curvature
  
// choose correct sign of dl, see page 58
// "arc length increases in the direction in which B points"  in 
// W. D. D'haeseleer, W. N. G. Hitchon, W. I. van Rij, S. P. Hirshman, 
// and J. L. Shohet, Flux Coordinates and Magnetic Field Structure.
// (Springer-Verlag, New York, 1991).
  if(curvature) {
    *curvature = k*sign(B2*dr);
  }
}

#if 0



  Vector3d b = b1+b2;
  b.norm();
  Vector3d n = n1+n2;
  n.norm();
  Vector3d dr = r2-r1;
  dl = REAL(dr.abs());

// choose correct sign of dl, see page 58
// "arc length increases in the direction in which B points"  in 
// W. D. D'haeseleer, W. N. G. Hitchon, W. I. van Rij, S. P. Hirshman, 
// and J. L. Shohet, Flux Coordinates and Magnetic Field Structure.
// (Springer-Verlag, New York, 1991).
  int sgn = sign(b*dr);
  Vector3d k = (b2*sgn)/dl;
  double kN = k*n;               // normal curvature
  double kG = k*(n^b)*sgn;       // geodesic curvature
//7  k -= kN*n1;                 // geodesic curvature vector
//7  double kG = k.abs();        // geodesic curvature
  


#endif



//****************************************************************************
//****************************************************************************
// not used
void CStconfig::advanceDataForEpsEff(Vector3d &mag, Vector3d &r, double &B, double &grads, double &kG, const double &dphi)
{
//see epsEff()  volatile ctStack::saveState state(ct);
  
  mag = magCoordFieldLine(mag, mag[2]+dphi); // advance along field line

  {
    double dphi=0.5*degree; 1./1024;        // 5mm/Rmajor, where Rmajor=5m for W7X  or  (1./1024) rad == 0.056 degree
    Vector3d r1, r2;            // positions(cartesian)
    Vector3d n1, n2;            // normals
    Vector3d b1 = magCoordFieldLine(mag,mag[2]-dphi); // move slightly along field line
    Vector3d b2 = magCoordFieldLine(mag,mag[2]+dphi); // move slightly along field line

    NJacobian(b1);
    b1 = getBxyz();
    r1 = getRxyz();
    n1 = getGradsxyz();
    B   = getmodB();

    NJacobian(b2);
    b2 = getBxyz();
    r2 = getRxyz();
    n2 = getGradsxyz();
    B   += getmodB();
    B   /= 2;
    r = (r1+r2)/2;

    Vector3d b = b1+b2;
    b1.norm();
    b2.norm();

    n1 += n2;  

    grads = n1.abs()/2;

    r2 -= r1;
    b2 -= b1;

    double dl = r2.abs();
    Vector3d k = b2/dl;
    Vector3d g = n1^b;
    g.norm();
    kG = k*g;                   // geodesic curvature
  }
}

//****************************************************************************
//****************************************************************************
void CStconfig::FbsInvalidate() 
{
  Fbs.invalidate();
  Fbs.clear();
}

void CStconfig::FbsSetXeffParam(double xmin, double xmax, int size) 
{
  Fbs.maxP = size;
  Fbs.xmin = xmin;
  Fbs.xmax = xmax;
}

void CStconfig::FbsSetSlabelParam(double smin, double smax, int size) 
{
  Fbs.maxP = size;
  Fbs.xmin = sqrt(smin);
  Fbs.xmax = sqrt(smax);
}

void CStconfig::FbsSetTracingParam(int turns, bool doAccuracyTest, double dphi)
{
  Fbs.turns= turns;
  Fbs.dphi = dphi;     // twopi/360 = 0.017453 -> 1degree;
  Fbs.doAccuracyTest = doAccuracyTest;
}

void CStconfig::FbsSetIotaParam(bool avoidResonances, bool useRationalIota, int iotaDenominator) 
{
  Fbs.useRationalIota=useRationalIota;
  Fbs.iotaDenominator=iotaDenominator;
  Fbs.avoidResonances=avoidResonances;
}

void CStconfig::FbsSetMagnMomentParam(int lambdaLevels, double atanhLambdaMax) 
{
  Fbs.numMuLevels = lambdaLevels;
  Fbs.xmumax = atanhLambdaMax;
}

void CStconfig::FbsCreate() 
{
  if(Fbs.useRationalIota) {
    Fbs.avoidResonances=false;
    Fbs.doAccuracyTest=false;
  }
  Fbs.create(this); 
}

double CStconfig::Fbsg2(double s) 
{
  Fbs.create(this); 
  return Fbs.g2(s);
}

double CStconfig::Fbsu(double s) 
{
  Fbs.create(this); 
  return Fbs.u(s);
}

double CStconfig::FbsuB2(double s) 
{
  Fbs.create(this); 
  return Fbs.uB2(s);
}

double CStconfig::Fbsu2B2(double s) 
{
  Fbs.create(this); 
  return Fbs.u2B2(s);
}
double CStconfig::FbsB2(double s) 
{
  Fbs.create(this); 
  return Fbs.B2(s);
}

double CStconfig::Fbsg2Vmec(double s) 
{
  Fbs.create(this); 
  return Fbs.g2Vmec(s);
}

double CStconfig::Fbsg2Graz(double s) 
{
  Fbs.create(this); 
  return Fbs.g2Graz(s);
}

double CStconfig::Fbsg2Dkes(double s) 
{
  Fbs.create(this); 
  return Fbs.g2Dkes(s);
}

double CStconfig::FbsVmec(double s) 
{
  Fbs.create(this); 
  return Fbs.Vmec(s);
}

double CStconfig::FbsGraz(double s)  
{
  Fbs.create(this); 
  return Fbs.Graz(s);
}

double CStconfig::FbsDkes(double s) 
{
  Fbs.create(this); 
  return Fbs.Dkes(s);
}

double CStconfig::GrazLambdab(double s) 
{
  Fbs.create(this); 
  return Fbs.GrazLambda(s);
}

double CStconfig::FbsFtrap(double s) 
{
  Fbs.create(this); 
  return Fbs.ftr(s);
}

//****************************************************************************
// linear interpolation is used
double CStconfig::Fbsg4(double s, double mu)
{
  Fbs.create(this); 

  int ir = Fbs.x.bsearch(sqrt(s));
  double dr = Fbs.x[ir+1]-Fbs.x[ir];
  
  double xm = 0.5*log((1+mu)/(1-mu));
  double dx = Fbs.xmm[1];
  double px = xm/dx;
  int jx = int(px);
  jx = mmin(jx,int(Fbs.xmm.size()-2));
  jx = mmax(jx,0);
  px -= jx;
  double f1=Fbs.g4_[ir]  [jx];
  double f2=Fbs.g4_[ir+1][jx];
  f1 += (Fbs.g4_[ir]  [jx+1]-f1)*px;
  f2 += (Fbs.g4_[ir+1][jx+1]-f2)*px;
  f1 += (f2-f1)*(sqrt(s)-Fbs.x[ir])/dr;
  return f1;
}

double CStconfig::Fbsg4Vmec(double s, double mu)
{
  return Fbsg4(s,mu)*r(s)/(2*s);
}

//****************************************************************************
void CStconfig::calculateFbs()
{
  if(!mOK) return;
  int size = Fbs.maxP;
  Threads2<CStconfig,CStconfig> threads(this, &CStconfig::calculateFbsExe,this,"calculateFbs");  
  int low = 0, upper = low+size-1;
  threads.run(low, upper, nThreads); 
}

//****************************************************************************
void CStconfig::calculateFbsExe(CStconfig * mc, int low, int upper, bool setup)
{
  if(setup) {  // setup part: it allocates arrays and initialize them
    int maxP = upper-low+1;
    Fbs.resize(); // detach data if they are shared 
    if(Fbs.x.empty()) return;
    if(Fbs.x.size()!=maxP) return;
    for(int i=0; i<maxP; i++) {
      double dx = (Fbs.xmax-Fbs.xmin)/(maxP-1);
      double x  = Fbs.xmin + dx*i; 
      if(Fbs.avoidResonances) {
        if(isRationalIota(x*x)) x += dx*0.3;  // try to avoid rational iota
        if(isRationalIota(x*x)) x += dx*0.1;  // try to avoid rational iota
      }
      Fbs.x.rw()[i]  = x;   // x=sqrt(s)
    }
    Fbs.Btruncation = truncate(); // save truncation level of B
    for(int k=0; k<Fbs.numMuLevels; k++) {
      const double dx = Fbs.xmumax/(Fbs.numMuLevels-1);
      Fbs.xmm.rw()[k] = k*dx; //  mu = tanh(xmm), where 0<=mu<=1 is the magnetic moment 
    }
    return;
  }
  else 
  {  // exe part: it will be called from threads to fill arrays
    if(mc->Fbs.x.empty()) return;
    __allFbs retvalue;
    int turns = mc->Fbs.turns;
    for(int i=low; i<=upper; i++) {
      double x = mc->Fbs.x[i];
      double s = x*x;
      mc->bscurrentFactor(s,retvalue);
      Fbs.fbsGraz.w(i) = retvalue.Graz;
      Fbs.fbsDkes.w(i) = retvalue.Dkes;
      Fbs.fbsVmec.w(i) = retvalue.Vmec;
      Fbs.ftrapped.w(i) = retvalue.ftrapped;
      Fbs.lambda_b.w(i) = retvalue.lambda_b;
      Fbs.g2vmec.w(i)   = retvalue.g2vmec;
      Fbs.g2graz.w(i)   = retvalue.g2graz;
      Fbs.g2dkes.w(i)   = retvalue.g2dkes;
      Fbs.g2_.w(i)      = retvalue.g2_;
      for(int k=0; k<mc->Fbs.numMuLevels; k++) {
        Fbs.g4_(i,k) = retvalue.g4_[k];
      }
      Fbs.u_.w(i)    = retvalue.u_   ;
      Fbs.uB2_.w(i)  = retvalue.uB2_ ;
      Fbs.u2B2_.w(i) = retvalue.u2B2_;
      Fbs.B2_.w(i)   = retvalue.B2_  ;
    }
  }
}

//******************************************************************
// Set random starting point.
//  srand()  from <stdlib.h>
// if sw==0  then random seed value
static void SeedRandom(int sw=0)
{
  if(sw)  srand(sw);
  else    srand((unsigned)time(NULL));
}
//******************************************************************
// Returns uniformly distributed random number between -0.5 and 0.5
//  rand() from <stdlib.h>
static double Random(void)
{
  static double norm = 1./(double)RAND_MAX;
  return (rand()*norm - 0.5);
}

//****************************************************************************
// dphi = 0.01227 = twopi/512
// int numMuLevels = 513;       // # of b levels
// int numTurns = 100;          // # of full turns of Bfield line around torus
void CStconfig::bscurrentFactor(double s, __allFbs &fbs)
{
  volatile ctStack::saveState state(ct);

  const double relEps=0.02;   // if( |resultBsPrev-resultBs| < |resultBs|*relEps ) ++accuracyCnt;
                              //    if(accuracyCnt>10) break;
  double xmax = Fbs.xmumax;   // mu_max = tanh(xmax)
  double dphi = Fbs.dphi;

  int numMuLevels = Fbs.numMuLevels;
  int numTurns    = Fbs.turns;
  bool BOOZER_COORDINATES = Fbs.useBoozerCoord; // if true then use properties of boozer coordinates
  bool boozerCrd = !(tokamakEQDSK||vmecCoordinates||fpCoordinates);
  if(!boozerCrd) BOOZER_COORDINATES = false;

  double iotaRatnl=iota(s);
  int num,den; 
  num=Fbs.iotaDenominator;
  den=Fbs.iotaDenominator;
  if(Fbs.useRationalIota) {
    //double iotaR = cfraction(iota(s), num,den);                                                                
    double iotaR = iotaRational(s,Fbs.iotaDenominator,num,den);
    iotaRatnl = iotaR!=0?iotaR:iotaRatnl;
    numTurns = den;
  }

  //SeedRandom();

  dphi = dphi==0?(0.4*degree):dphi;   // toroidal-angle increment.
  
  REAL *mu_dmu  = new REAL[numMuLevels];  // array of mu*dmu levels;  mu = tanh(x), where x > 10 ( tanh(10)==1-4e-9 )
  REAL *mu_Bmax  = new REAL[numMuLevels]; // mu/Bmax
  REAL *g4overXi  = new REAL[numMuLevels];
  REAL *g4Average = new REAL[numMuLevels];
  REAL *xiAverage = new REAL[numMuLevels]; 
  REAL *g5        = new REAL[numMuLevels]; 

  Vector3d magcoordBmin(s,0.,0.), magcoordBmax(s,0.,0.);

  int restartCnt=0;
restart:
  findBminBmaxPosition(-s, magcoordBmin, magcoordBmax, restartCnt>0);

#if CRAIG_BMN_TEST_FIELD
  magcoordBmax[1] = pi;
  magcoordBmax[2] = pi/10;
  std::cout<<" magcoordBmax="<<magcoordBmax/degree<<std::endl;
#endif

  REAL Bmax;
  Vector3d r0,r1,mag; 
  mag = magcoordBmax;          // starting point for integrating
  r1 = r0 = mag2xyz(mag);
  Bmax = REAL(B(magcoordBmax));//*(1+1e-8);
//////04jun2011
////NJacobian(magcoordBmax);
////Vector3d Bvec   = getBcyl();
////Bmax = Bvec.abs(); 
//////04jun2011

  const double dx = xmax/(numMuLevels-1);

  for(int k=0; k<numMuLevels; k++) {
    double x = k*dx;
    double mu = REAL(tanh(x));
    mu_dmu[k] = dx*mu*(1-square(mu));
    mu_Bmax[k] = mu/Bmax;
    xiAverage[k] = 0;
    g4overXi[k]  = 0;
    g4Average[k] = 0;
    g5[k] = 0;
  }

  REAL gradSAverage(0); // <|grad(s)|>
  REAL g2overB2(0);   
  REAL Norm(0),g2(0),g2Average(0), B2Average(0);
  REAL resultBs(0), resultBsPrev(0), resultFtr(0), resultFtrPrev(0);
  REAL u = 0 ;   
  REAL u_ = 0 ;  // =<u>
  REAL uB2_ = 0; 
  REAL u2B2_ = 0;

  gradSAverage=g2overB2=Norm=g2=g2Average=B2Average=0;
  resultBs=resultBsPrev=resultFtr=resultFtrPrev=0;

  int accuracyCnt=0;
  int turnsCnt=0;
  int nsteps = int(twopi/dphi);
  dphi = twopi/nsteps;

#define  USE_KG 0    // use geodesic curvature if true
#if USE_KG
  Vector3d B1, gradS1;
  REAL Bmod,gradSmod,kG,dl;
  advanceDataForEpsEff02(mag, r1, B1, gradS1, Bmod, gradSmod, kG, dl, 0,0,iotaRatnl);  //initialize
//  advanceDataForEpsEff02(mag, r1, B1, gradS1, Bmod, gradSmod, kG, dl, 0, iotaRatnl);  //initialize
  r0 = r1;    // save starting position 
#endif
  for(int i=0 ; i<numTurns; i++) {
    ++turnsCnt;
    if(BOOZER_COORDINATES) { //begin BOOZER_COORDINATES
      for(int j=0 ; j<nsteps; j++) { //this loop is valid only for boozer coordinates
        mag = magCoordFieldLine(mag, mag[2]+dphi,iotaRatnl); //,200./231.); // advance along field line
        Vector3d gr = getCovaGradB(mag);
        NJacobian(mag);
        Vector3d grads  = getGrads();
        REAL gradSmod = grads.abs();
        REAL B = REAL(getmodB());
        double BGsGBdl_B = mu0/Flux(1.)*(Ip(s)*gr[1]-It(s)*gr[2])*dphi;   // [Bxgrad(s)]*grad(B)*dl/B == mu0*Ipol/Flux(1) * dB/dtheta * dphi
        REAL Bp2 = B*B;
        double dl_B = dphi/Bp2;     //valid only for boozer coordinates 
        g2overB2  += BGsGBdl_B*(-2/(Bp2*B));    // I2 = integral from l_max to l
        g2 = Bp2*g2overB2;                      // g2 = B^2*I2
        g2Average += g2*dl_B;                   // = integral {g2*dl/B} from l_max to infinity
        B2Average += Bp2*dl_B;                  // = integral {B^2 *dl/B} from l_max to infinity
        Norm += dl_B;                           // = integral {dl/B} from l_max to infinity
        gradSAverage += gradSmod*dl_B;
        if(B/Bmax>1) {
          std::cerr<<"s="<<s<<" turns="<<turnsCnt<<" steps="<<j<<" B/Bmax="<<B/Bmax;
          if(++restartCnt<4) {
            std::cerr<<"  ->calculation restarted"<<std::endl<<std::flush;
            magcoordBmax = mag;
            goto restart;
          }
          else {
            std::cerr<<std::endl<<std::flush;
          }
        }
        for(int k=0; k<numMuLevels; k++) {
          double w = 1-mu_Bmax[k]*B;  // 1 - mu*B/Bmax
          if(w<0) continue;
          REAL xi = sqrt(w);
          g4overXi[k]  += BGsGBdl_B *(0.5*mu_Bmax[k]/cube(xi)); // I4 = integral from l_max to l
          REAL g4 = xi*g4overXi[k];                 // g4 = xi*I4
          g4Average[k] += g4*dl_B;                  // = integral {g4*dl/B} from l_max to infinity
          xiAverage[k] += xi*dl_B;                  // = integral {xi*dl/B} from l_max to infinity
        }
      } // one turn
    } //end BOOZER_COORDINATES
  else {
      for(int j=0 ; j<nsteps; j++) {
        REAL dl_B, B, gradSmod;
        double BGsGB,BGsGBdl_B,BcontraPhi;
  #if USE_KG
        advanceDataForEpsEff02(mag, r1, B1, gradS1, Bmod, gradSmod, kG, dl,0, dphi, iotaRatnl);
        B = Bmod;
        BcontraPhi = getBcontra()[2];
        //dl_B = fabs(dphi/BcontraPhi);   // the same as dl/B;
        dl_B = dphi/BcontraPhi;        //7Nov2011 dl/B must have sign 
        BGsGB = gradSmod*kG*(B*B);  // the same as (Bvec^gradB)*grads;
        BGsGBdl_B = BGsGB*dl_B;
  #else
        mag = magCoordFieldLine(mag, mag[2]+dphi,iotaRatnl); //,200./231.); // advance along field line
        Vector3d gr = getCovaGradB(mag);
        NJacobian(mag);
        Vector3d grads  = getGrads();
        Vector3d gradTh = getGradtheta();
        Vector3d gradPh = getGradphi();
        Vector3d Bvec   = getBcyl();
        //Vector3d gradBpar  = gr[1]*gradTh + gr[2]*gradPh;// + gr[0]*grads;
        Vector3d gradBpar  = gradTh;
        gradBpar  *= gr[1];
        gradBpar  += gr[2]*gradPh;
        //gradB   += gr[0]*grads;
        B = REAL(getmodB());
   ////B = Bvec.abs(); //04jun2011, see findBminBmaxPosition
        BGsGB = (Bvec^grads)*gradBpar;
        BcontraPhi = getBcontra()[2]; //double BcontraPh = Bvec*gradPh;
        //dl_B = fabs(dphi/BcontraPhi);        // the same as dl/B;
        dl_B = dphi/BcontraPhi;        //7Nov2011 dl/B must have sign
        BGsGBdl_B = BGsGB*dl_B;
        gradSmod = grads.abs();
    #define AVOID_SINGULARITY  0   // see N.Nakajima et al, NF, vol.29, 605(1989)
    #if AVOID_SINGULARITY
          double a = mu0*Ip(s)/(iotaRatnl*Flux(1.));
          a = 0;
          Vector3d gradB = gradBpar + gr[0]*grads;
          double BGs_a_GB = BGsGB - a*Bvec*gradB;
    #endif
  #endif
  #undef USE_KG

  #if 0
        Vector3d r2     = getRxyz();
        REAL dl = REAL((r2-r1).abs());
        r1 = r2;
        dl_B = dl/B;
  #else
        //r1 = getRxyz();
  #endif

        REAL Bp2 = B*B;
        g2overB2  += BGsGBdl_B*(-2/(Bp2*B));  // I2 = integral from l_max to l
        g2 = Bp2*g2overB2;                    // g2 = B^2*I2
#if 1
        u  -= BGsGBdl_B*(-2/(Bp2*B));         //  integral from l_max to l
        u_    += u*dl_B;                      // = integral {u*dl/B} from l_max to infinity
        uB2_  += Bp2*u*dl_B;                  // = integral {u*B^2*dl/B} from l_max to infinity
        u2B2_ += Bp2*u*u*dl_B;                // = integral {u^2*B^2*dl/B} from l_max to infinity
#endif

   //g2 = (2*s)/r(s);    //TEST
        g2Average += g2*dl_B;                 // = integral {g2*dl/B} from l_max to infinity
        B2Average += Bp2*dl_B;                // = integral {B^2 *dl/B} from l_max to infinity
        Norm += dl_B;                         // = integral {dl/B} from l_max to infinity
        gradSAverage += gradSmod*dl_B;

        if(B/Bmax>1) {
          std::cerr<<"s="<<s<<" turns="<<turnsCnt<<" steps="<<j<<" BGsGB="<<BGsGB<<" B/Bmax="<<B/Bmax;
          ////std::cerr<<std::endl<<std::flush;
          ////Vector3d c = magcoordBmax;
          ////Vector3d m = mag;
          ////c[1] = mod2pi(c[1])/degree;
          ////c[2] = modPeriod(c[2])/degree;
          ////std::cerr<<c<<std::endl<<std::flush;
          ////m[1] = mod2pi(m[1])/degree;
          ////m[2] = modPeriod(m[2])/degree;
          ////std::cerr<<m<<std::endl<<std::flush;
          if(++restartCnt<3) {
            std::cerr<<"  ->calculation restarted"<<std::endl<<std::flush;
            magcoordBmax = mag;
            goto restart;
          }
          else {
            std::cerr<<std::endl<<std::flush;
          }
        }
        for(int k=0; k<numMuLevels; k++) {
          double w = 1-mu_Bmax[k]*B;
          if(w<0) continue;
          REAL xi = sqrt(w);
          g4overXi[k]  += BGsGBdl_B *(0.5*mu_Bmax[k]/cube(xi)); // I4 = integral from l_max to l
          REAL g4 = xi*g4overXi[k];                 // g4 = xi*I4
  #if AVOID_SINGULARITY
          g5[k]        += BGs_a_GB*(0.5*mu_Bmax[k]/cube(xi))*dl_B;
          g4 = xi*g5[k];
          g4 -= a*(B==Bmax?0:(1-xi/sqrt(1-mu_Bmax[k]*Bmax)));
  #endif
   // g4 = (2*s)/r(s);   //TEST
          g4Average[k] += g4*dl_B;                  // = integral {g4*dl/B} from l_max to infinity
          xiAverage[k] += xi*dl_B;                  // = integral {xi*dl/B} from l_max to infinity
        }
      } // one turn 
    }
    REAL IntBs =0;
    REAL IntFtr=0;
    for(int k=0; k<numMuLevels; k++) {
      double w = mu_dmu[k]/xiAverage[k];
      IntBs  += w*g4Average[k];
      IntFtr += w*Norm;
      ////IntBs  += dx*g4Average[k]*mu[k]*(1-square(mu[k]))/xiAverage[k];
      ////IntFtr += dx*        Norm*mu[k]*(1-square(mu[k]))/xiAverage[k];
    }

    resultBs  = (g2Average - 0.75*B2Average*IntBs /square(Bmax) )/Norm;
    resultFtr = (Norm      - 0.75*B2Average*IntFtr/square(Bmax) )/Norm;   //the fraction of trapped particles

    Vector3d magerr = mag-magcoordBmax;
    magerr[1] = mod2pi(magerr[1]+pi)-pi;
    magerr[2] = mod2pi(magerr[2]+pi)-pi;

    //////double dr = (r0-r1).abs();
    //////dr /= this->r(s);
    //////if(turnsCnt>30&&dr<0.3*degree) {
    if(magerr.abs()<1e-4) {
     // break;
      std::cerr.precision(3);
      std::cerr.setf(std::ios::scientific|std::ios::showpos);
      std::cerr<<"periodicity="<<magerr.abs()
               <<" g2="<<g2<<" sqrt(s)="<<sqrt(s);
      if(Fbs.useRationalIota) {
        std::cerr<<" iota="<<num<<"/"<<den
                 <<" iota_err="<< iota(s)-double(num)/double(den);
      }
      else {
        std::cerr<<" iota="<<iota(s);
      }
      std::cerr<<" turns="<<turnsCnt<<std::endl<<std::flush;
      //break;
    }
    if(!Fbs.useRationalIota && Fbs.doAccuracyTest) {
      if(fabs(resultBsPrev - resultBs)  <= fabs(resultBs)*relEps) ++accuracyCnt;
      if(accuracyCnt>10&&turnsCnt>30) break;
    }
    resultBsPrev  = resultBs;
    resultFtrPrev = resultFtr;
    ////double ds = s*0.02*Random();
    ////double s1 = s+ds;
    ////mag[0] = s1<=0?s:s1;
  }   // for(int i=0 ; i<numTurns; i++) 

  //std::cerr<<"bscurrentFactor(): s="<<s<<"  turns="<<turnsCnt<<std::endl;

  gradSAverage /= Norm;  // <|grad(s)|>
  B2Average    /= Norm;  // <B^2> 
  g2Average    /= Norm;
  u_        /= Norm;
  uB2_      /= Norm;
  u2B2_     /= Norm;

  double Fbs      = resultBs;                // for s as the flux label 
  double Fbs_graz = resultBs/gradSAverage;   // for r_graz as the flux label:  dr_graz/ds = 1/<|grad(s)|>
  double Fbs_vmec = resultBs*r(s)/(2*s);     // for r_vmec as the flux label:  dr_vmec/ds = a*sqrt(s)/(2s)
  
  double g2graz   = g2Average/gradSAverage;
  double g2vmec   = g2Average*r(s)/(2*s);

  double dkesFactor = drdkes_dr(s)*r(s)/(2*s);  // dr_dkes/ds = drdkes_dr_vmec(s)*dr_vmec/ds
  double Fbs_dkes = resultBs*dkesFactor;     // for r_dkes as the flux label  
  double g2dkes   = g2Average*dkesFactor;

  //REAL B0 = (REAL)Bmn_sp(s,0,0);              // Bmn_sp(0,0,0);
  REAL B0 = Bmax;
  double lambda_b = Fbs_graz*B0*B0/B2Average; // lambda_b in Nemov,Kalyuzhnyj,Kasilov et al, PPCF vol.46.179(2004)
  
  g2graz *= B0*B0/B2Average; // lambda_b in Nemov,Kalyuzhnyj,Kasilov et al, PPCF vol.46.179(2004)


  fbs.g2_   = g2Average;
  fbs.B2_   = B2Average;
  fbs.u_    = u_;
  fbs.uB2_  = uB2_;
  fbs.u2B2_ = u2B2_;

  fbs.g4_.assign(g4Average,numMuLevels); 
  fbs.g4_.scale(1./Norm);
  //fbs.g4_.scale(r(s)/(2*s));  //for r_vmec as the flux label
  
  delete mu_dmu;
  delete mu_Bmax;
  delete g5;
  delete g4overXi;
  delete g4Average;
  delete xiAverage;

  fbs.Graz = Fbs_graz; 
  fbs.Dkes = Fbs_dkes;    
  fbs.Vmec = Fbs_vmec; 
  fbs.lambda_b = lambda_b;
  fbs.ftrapped = resultFtr;
  fbs.g2vmec = g2vmec;
  fbs.g2graz = g2graz;
  fbs.g2dkes = g2dkes;

  //for(int k=0; k<numMuLevels; k++) g4Average[k] *= (r(s)/(2*s)/Norm);
  //fbs.g2vmec = -g4Average[numMuLevels/10];

  std::cerr.setf(std::ios::fixed,std::ios::showpos|std::ios::scientific);

  ////return Vector3d(Fbs_vmec,resultFtr,lambda_b);
  
  //double lambda_b_over_Fbs_vmec = B0*B0/B2Average/gradSAverage * (2*s)/r(s);
  //return Vector3d(Fbs_vmec,lambda_b_over_Fbs_vmec,lambda_b);
  //return Vector3d(Fbs_vmec,lambda_b/Fbs_vmec,lambda_b);
  //return Vector3d(Fbs_vmec,resultFtr,g2_vmec);
}


/*
** find rational approximation to given real number
** David Eppstein / UC Irvine / 8 Aug 1993
**
** With corrections from Arno Formella, May 2008
**
** usage: a.out r d
**   r is real number to approx
**   d is the maximum denominator allowed
**
** based on the theory of continued fractions
** if x = a1 + 1/(a2 + 1/(a3 + 1/(a4 + ...)))
** then best approximation is found by truncating this series
** (with some adjustments in the last term).
**
** Note the fraction can be recovered as the first column of the matrix
**  ( a1 1 ) ( a2 1 ) ( a3 1 ) ...
**  ( 1  0 ) ( 1  0 ) ( 1  0 )
** Instead of keeping the sequence of continued fraction terms,
** we just keep the last partial product of these matrices.
*/

//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//#include <iostream>


#ifdef __GNUC__
 #define LLONG int64_t 
#else
 #define LLONG __int64 
#endif

double cfraction(double iota, int &numi, int &deni)
{
    LLONG m[2][2];
    LLONG maxden,maxden1,maxden2;
    LLONG ai;
    LLONG num(0),den(0);  // returned values
    double x, startx, err=2;

    startx = x = iota;
    maxden1 = numi;
    maxden2 = deni;

    /* initialize matrix */
    m[0][0] = m[1][1] = 1;
    m[0][1] = m[1][0] = 0;

 for (maxden=maxden1; maxden<=maxden2; maxden+=10 ) {
    /* loop finding terms until denom gets too big */
    while (m[1][0]*(ai=LLONG(x))+m[1][1] <= maxden) {
        LLONG t = m[0][0] * ai + m[0][1];
        m[0][1] = m[0][0];
        m[0][0] = t;
        t = m[1][0] * ai + m[1][1];
        m[1][1] = m[1][0];
        m[1][0] = t;
        if(x==double(ai)) 
          break;     // AF: division by zero
        x = 1/(x-double(ai));
        if(x>9e18)   // if(x>double(0x7FFFFFFFFFFFFFFF)) 
          break;  // AF: representation failure
    }
       
    /* now remaining x is between 0 and 1/ai */
    /* approx as either 0 or 1/m where m is max that will fit in maxden */
    /* first try zero */
    LLONG num1 = m[0][0];
    LLONG den1 = m[1][0];
    /* now try other possibility */
    LLONG c = (maxden - m[1][1]) / m[1][0];
    LLONG num2 = m[0][0] * c + m[0][1];
    LLONG den2 = m[1][0] * c + m[1][1];

    ////std::cerr <<num1<<"/"<<den1<<"  err="<< startx-double(num1)/double(den1)<<std::endl;
    ////std::cerr <<num2<<"/"<<den2<<"  err="<< startx-double(num2)/double(den2)<<std::endl;

    if(den2>den1) {
      double e = fabs(startx-double(num2)/double(den2));
      if(e<err) {
        err=e;
        num = num2;
        den = den2;
      }
  //    std::cerr<<"return="<<num2<<"/"<<den2<<"  err="<< startx-double(num2)/double(den2)<<std::endl;
    }
    else {
      double e = fabs(startx-double(num1)/double(den1));
      if(e<err) {
        err=e;
        num = num1;
        den = den1;
      }
   //   std::cerr<<"return="<<num1<<"/"<<den1<<"  err="<< startx-double(num1)/double(den1)<<std::endl;
    }
  
 }

 numi = int(num);
 deni = int(den);

 std::cerr<<"=final=return="<<numi<<"/"<<deni<<"  err="<< (startx-double(num)/double(den) ) <<std::endl;

 return double(num)/double(den);
}

//   estimate the best rational approach to iota 
// --------------------------------------------------------------------
//     Author:   H. Maassberg 
//     created:  Feb. 1995 
// --------------------------------------------------------------------
//     The best rational approach to the given value of the rotational 
//     transform, XIO, is defined by min{ abs(xio - MP/NT) } where the 
//     poloidal and toroidal modes, MP and NT, have no common divisor, 
//     and with the condition: NTG-INC .le. NT .le. NTG+INC. 
// --------------------------------------------------------------------
//     input parameters: 
//     XIO       value of rotational transform, iota 
//               -> the rational approximation is defined by MP/NT 
//     NTG       input guess for toroidal mode number 
//     INC       increment for estimating NT: NTG-INC.le.NT.le.NTG+INC 
//     output parameters: 
//     MP        poloidal mode number 
//     NT        toroidal mode number 
//     IER       error control parameter 
//               = 0   no error 
//               = 1   no rational approach for XIO found 
//               = 2   check XIO and NTG 
// ---------------------------------------------------------------------
//     subroutines required:    CHCMDV 
// ---------------------------------------------------------------------
static int chcmdv_(int &m, int &n, int &nd);

static int bncraio_2(double xio, int ntg, int inc, int &mp, int &nt)
{
    double xior;
    double acc;
    double aio;
    double r__1;
    int nt0;
    int ncd;
    int lio, mpl, ntl, mpu, ntu;
    int ier=2;
    nt = 1;
    mp = 0;
    aio = fabs(xio);
    if (aio == 0  || ntg <= 1) {
      return ier;
    }
    ier = 1;
    lio = 1;
    if (xio < 0) lio = -1;
    ntl = ntg - abs(inc);
    ntu = ntg + abs(inc);
    nt0 = ntl;
    acc = aio;
L1:
    mpl = (int) (aio * nt0);
    chcmdv_(mpl, nt0, ncd);
    if (ncd <= 1) {
      xior = mpl / (double) nt0;
      if ((r__1 = aio - xior, fabs(r__1)) <= acc) {
          acc = (r__1 = aio - xior, fabs(r__1));
          mp = mpl * lio;
          nt = nt0;
          ier = 0;
      }
    }
    mpu = mpl + 1;
    chcmdv_(mpu, nt0, ncd);
    if (ncd <= 1) {
      xior = mpu / (double) nt0;
      if ((r__1 = aio - xior, fabs(r__1)) <= acc) {
          acc = (r__1 = aio - xior, fabs(r__1));
          mp = mpu * lio;
          nt = nt0;
          ier = 0;
      }
    }
    if (nt0 >= ntu) return ier;
    ++nt0;
    goto L1;
} 
//   check on common divisor 
//   Author:   H. Maassberg 
//   created:  Feb. 1995 
// -------------------------------------------------------------------- 
//   The integers N and M are analyzed for a common divisor. On exit, 
//   ND contain the lowest common divisor, e.g., ND = 1 if no common 
//   divisor is found. 
static int chcmdv_(int &m, int &n, int &nd)
{
    int i0, id, ma, na, kl, iu, ku;
    nd = 0;
    na = abs(n);
    ma = abs(m);
    if (na == 0 || ma == 0) return 0;
    nd = 1;
    if (na == 1 || ma == 1) return 0;
    i0 = 0;
    id = 2;
    iu = ((na<ma)?na:ma) / 2; // iu = min(na,ma) / 2;
L1:
    if (id * (ma / id) == ma) {
      if (id * (na / id) == na) {
        nd = id;
        return 0;
      }
    }
    ++i0;
    id = (i0 << 1) + 1;
    if (id <= iu)  goto L1;
    kl = ((na<ma)?na:ma); // kl = min(na,ma);
    ku = ((na>ma)?na:ma); // ku = max(na,ma);
    if (kl * (ku / kl) == ku)  nd = kl; 
    return 0;
} 


//   estimate the best rational approach to iota 
// --------------------------------------------------------------------
//     Author:   H. Maassberg 
//     created:  Feb. 1995 
//     revision: translated to C++ by y.turkin, June 2011 
// --------------------------------------------------------------------
//     The best rational approach to the given value of the rotational 
//     transform, XIO, is defined by min{ abs(xio - MP/NT) } where the 
//     poloidal and toroidal modes, MP and NT, have no common divisor, 
//     and with the condition: NTG-INC .le. NT .le. NTG+INC. 
// --------------------------------------------------------------------
//     input parameters: 
//     XIO       value of rotational transform, iota 
//               -> the rational approximation is defined by MP/NT 
//     NTG       input guess for toroidal mode number 
//     INC       increment for estimating NT: NTG-INC.le.NT.le.NTG+INC 
//     output parameters: 
//     MP        poloidal mode number 
//     NT        toroidal mode number 
//     IER       error control parameter 
//               = 0   no error 
//               = 1   no rational approach for XIO found 
//               = 2   check XIO and NTG 
// ---------------------------------------------------------------------
//     subroutines required: gcd 
// ---------------------------------------------------------------------
static int gcd(int a, int b)
{
// greatest common divisor (GCD) by Euclidean algorithm, Knuth, page 318.
// Knuth, Donald E. (1997). The Art of Computer Programming, 
// Volume 1: Fundamental Algorithms (3rd ed.). Addison-Wesley. ISBN 0-201-89683-4.  
  if(a==0) return b;
  a=abs(a);
  b=abs(b);
  while(b!=0) {
    if(a > b) a -= b;
    else      b -= a;
  }
  return a;
}

static int bncraio(double xio, int ntg, int inc, int &mp, int &nt)
{
  int ier = 2;
  if (xio==0||ntg<=1) return ier;
  
  ier = 1;
  nt = 1;
  mp = 0;
  double aio = fabs(xio);
  int lio = xio<0?-1:1;
  int ntl = ntg - abs(inc);
  int ntu = ntg + abs(inc);
  double acc = aio;

  for(int nt0 = ntl;nt0<=ntu;nt0++) {
    int mpl = int(aio*nt0);
    int mpu = mpl+1;
    for( ;mpl<=mpu; mpl++) {
      if(gcd(mpl,nt0) <= 1) {
        double iotaEps = fabs(aio-mpl/double(nt0));
        if (iotaEps <= acc) {
          acc = iotaEps;
          mp = mpl * lio;
          nt = nt0;
          ier = 0;
        }
      }
    }
  }
  return ier;
} 

// estimate the best rational aproximation to the iota: iota = numerator/denominator
double CStconfig::iotaRational(double s,int denom_guess,int &numerator,int &denominator)
{
  int inc = 11;
  int ier = bncraio(iota(s), denom_guess,inc,numerator,denominator);
  return ier>0?0:(numerator/double(denominator));
}

//****************************************************************************
//****************************************************************************
Vector3d CStconfig::curvature(const Vector3d &magCoord)
{
  volatile ctStack::saveState state(ct);

  Vector3d b,b1(magCoord),b2(magCoord); // positions(magCoord)
  Vector3d r1, r2;            // positions(cartesian)
  Vector3d n1, n2;            // normals
  double jac;
  double dphi=1./2048;       // 5mm/Rmajor, where Rmajor=5m for W7X  or  (1./1024) rad == 0.056 degree
  //////if(isVmecCoordinates()) {  // B-tracing if VMEC
  //////  Bfield Bfld(*this);
  //////  btrace.init(Bfld,dphi);
  //////  Vector3d c0 = mag2cyl(magCoord);
  //////  Vector3d c  = c0;
  //////  btrace.advance(-dphi,c);
  //////  b1 = cyl2mag(c);
  //////  b1 = getBxyz();
  //////  r1 = getRxyz();
  //////  n1 = getGradsxyz();
  //////  jac = getJ();

  //////  c  = c0;
  //////  btrace.advance(dphi,c);
  //////  b2 = cyl2mag(c);
  //////  b2 = getBxyz();
  //////  r2 = getRxyz();
  //////  n2 = getGradsxyz();
  //////  jac += getJ();
  //////}
  //////else 
  {
    //double i = iota(magCoord[0]);
    //b1[1]-=i*dphi;              // theta-iota*dphi
    //b1[2]-=dphi;                // phi - dhi

    //b2[1]+=i*dphi;              // theta+iota*dphi
    //b2[2]+=dphi;                // phi + dhi

    b1 = magCoordFieldLine(b1,b1[2]-dphi); // move slightly along field line
    b2 = magCoordFieldLine(b2,b2[2]+dphi); // move slightly along field line

    NJacobian(b1);
    b1 = getBxyz();
    r1 = getRxyz();
    n1 = getGradsxyz();
    jac = getJ();

    NJacobian(b2);
    b2 = getBxyz();
    r2 = getRxyz();
    n2 = getGradsxyz();
    jac += getJ();
  }


  b = b1+b2;
  b.norm();
  b1.norm();
  b2.norm();

//  n1.norm();
//  n2.norm();
  n1 += n2;
  n1.norm();
  r2 -= r1;
  b2 -= b1;
// choose correct sign of dl, see page 58
// "arc length increases in the direction in which B points"  in 
// W. D. D'haeseleer, W. N. G. Hitchon, W. I. van Rij, S. P. Hirshman, 
// and J. L. Shohet, Flux Coordinates and Magnetic Field Structure.
// (Springer-Verlag, New York, 1991).
  int sgn = sign(b*r2);
  double dl = sgn*r2.abs();
  Vector3d k = b2/dl;
  double kN = k*n1;     // normal curvature
  double kG = k*(n1^b); // geodesic curvature
//18.10.2011  kG = k*(n1^b)*sgn; 
#if 0  //31may2012
 {
    Vector3d Bvec   = getBcyl();
    Vector3d grads  = getGrads();
    Vector3d gradTh = getGradtheta();
    Vector3d gradPh = getGradphi();
    Vector3d gradients = gradTh - iota(magCoord[0])*gradPh;
    Vector3d clebschB =  grads^gradients;
    sgn = sign(Bvec*clebschB);
    kG = k*(n1^b)*sgn;      // geodesic curvature
  }
#endif
//7  k -= kN*n1;                 // geodesic curvature vector
//7  double kG = k.abs();        // geodesic curvature
  
  return Vector3d(kN,kG,jac/2);
}

//****************************************************************************
// see cyl2mag(Vector3d &cyl, Vector3d &magCoordGuess)
Vector3d CStconfig::xyz2mag(const Vector3d &R,const Vector3d &magCoordGuess)
{
  Vector3d cyl = R.toCylindrical();
  return cyl2mag(cyl,magCoordGuess);
}

//****************************************************************************
//  Solve F(u) = 0 by Newton's iteration
//    where   F=(R(u),fi(u),Z(u))-(r,f,z)
//
// INPUT:
//   Vector3d magCoordGuess == (s,theta,phi) -- guess value
//   Vector3d cyl == (r,fi,z) is the point in cylindrical coordinates,
//     for which we want to know magnetic coordinates(s,theta,phi)
// OUTPUT:
//  function returns magnetic coordinates(s,theta,phi) of the point (r,fi,z)
Vector3d CStconfig::cyl2mag(const Vector3d &cyl, const Vector3d &magCoordGuess)
{
  return cyl2mag(cyl, magCoordGuess,0);
}

// Jguess is the Jacobian matrix at guess point (=NJacobian(magCoordGuess) )
Vector3d CStconfig::cyl2mag(const Vector3d &cyl, const Vector3d &magCoordGuess, Vector3d *Jguess)
{
  Vector3d u = magCoordGuess;  // guess (s,theta,phi)
  if(fabs(u[0])<1e-14) {
    u[0]=1e-14;       // magnetic axis, set s=1e-14
    u[1]=0;           // theta = 0
  }
  Vector3d *J,F,Fs,Ft,Fp,du;
  bool Ok=false;
  bool prout = false;
  double invDet,invDet0;
//  double len=sqrt(cyl[0]*cyl[0]+cyl[2]*cyl[2]); // sqrt(r^2+z^2)
//  double err2 = square(mmax(mepsA,mepsR*len));
  double err2 = square(mepsA);
  double err16= 16*err2;
  double dR2old;
  double dR2;
  double epsGuess = 1e-4;
  int nItrBeforeNewGuess = 5; 
  int newGuessCount = nItrBeforeNewGuess; 
  int i=1, imax = 1000;    // max number of iterations
  invDet0=0;
#if 0
  J = NJacobian1Q(u);        // Jacobian using spline interpolation
  dR2=J[0].diffCyl2(cyl);    // distance^2 between J[0] and cyl
  if(dR2<=err2) { Ok=true; goto Finalize; } // convergence test
#else
  J = Jguess?Jguess:NJacobian1Q(u); // use linear interpolation
  dR2 = J[0].diffCyl2(cyl);  // distance^2 between J[0] and cyl
#endif

  for(i=1;i<imax;i++) {  // iteration loop
    F  = J[0]-cyl;        // F=(R(u),Fi(u),Z(u))-(r,fi,z), we solve F(u)=0; r,fi,z are given
    Fs = J[1];          // d(R,fi,Z)/ds   -- partial derivative on s
    Ft = J[2];          // d(R,fi,Z)/dtheta
    Fp = J[3];          // d(R,fi,Z)/dphi

    if(u[0]<s1st)  {
      if(invDet0==0) {  // calculate invDet0 if it zero   
        J    = NJacobian1Q(Vector3d(s1st,u[1],u[2]));
        invDet0 = 1./(J[1]*(J[2]^J[3])); // 1/determinant
      }
      invDet = invDet0;  // use determinant at point s1st for small 's'
    }
    else {
      invDet = 1./(Fs*(Ft^Fp)); // 1/determinant, where ^ is a cross product sign
    }

    du = -(Vector3d(F*(Ft^Fp),Fs*(F^Fp),Fs*(Ft^F)) * invDet);  // solve system of equation
    if(u[0]+du[0]<0)  du *= -0.9*u[0]/du[0];  //scale 'du' to avoid s<0
    if(fabs(du[1])>1) du /= fabs(du[1]);      //scale 'du' if delta(theta) is too big
#if 1
// backtracking loop: first try the full Newton step; 
// shorten the step size if the residual does not decrease
    double step = 1;
    int N=5;
    Vector3d uOld = u;
    /////double F2old = F.abs2();
    dR2old = dR2;
    while(--N) {
      u = uOld + step*du;        // try the Newton step
      J = NJacobian1Q(u);        // Jacobian using spline interpolation
      //////F = J[0]-cyl;              // F is the new residual
      //////if(F.abs2()>F2old) step *= 0.5;  
      dR2=J[0].diffCyl2(cyl);    // distance^2 between J[0] and cyl
      if(dR2>dR2old) step *= 0.5;   // shorten the step size if the residual does not decrease
      else  break;   // exit from the backtracking loop
    }
    //////dR2=J[0].diffCyl2(cyl);    // distance^2 between J[0] and cyl
    if(dR2<=err2) { Ok=true; break; } // convergence test
#else
    u += du;
#endif
    if(i==newGuessCount && epsGuess>2e-9) {
      u = cyl2guess(cyl,epsGuess);   // try new guess
      J = NJacobian1Q(u);        // Jacobian using spline interpolation
      dR2=J[0].diffCyl2(cyl);    // distance^2 between J[0] and cyl
      if(dR2<=err2) { Ok=true; break; } // convergence test
      invDet0=0;
      newGuessCount += 10;
      epsGuess  *= 0.1;
    } 
    else if((i+1)%100==0&&err2<err16) err2 *=1.4;
    else if(i==300) prout = true;
  }

Finalize:
  NJacobian2(INTERPOLATION_ORDER); // finalize Jacobian calculation

  //u[1] = mod2pi(u[1]);
//??Danger!! Don't do: u[2] = mod2pi(u[2]);
  ct().saveLastJacobian(mBnorm);     // save for next transformation, see cyl2guess

  if(!Ok) {
    std::cerr <<"CStconfig::cyl2mag(): NO CONVERGENCE!\n";
    ct().invalidate();   // destroy the data saved for guess obtaining 
    prout = true;  //    u.set(1e-10,0,cyl[1]); //??
  }
  if(prout) {
    std::cerr <<"# of iterations="<<i<<" Det= "<<1./invDet<<
    "\n   guess="<<magCoordGuess<<"\n  magCoord="<<u<<"\n     cyl="<<cyl<<std::endl;
    if(i>imax*0.9) ct().invalidate();   // destroy tha data saved for guess obtaining
  }

  ////if(i>=nItrBeforeNewGuess) {
  ////  std::cerr <<"# of iterations="<<i<<" Det= "<<1./invDet<< " Jguess= "<<Jguess<<
  ////  "\n     guess="<<magCoordGuess<<
  ////  "\n  magCoord="<<u<<
  ////  "\n  cyl="<<cyl<<std::endl;
  ////}
  if(u[0]<1e-13) {
    u[0]=1e-13;       // magnetic axis, set s=1e-13
    u[1]=0;           // theta = 0
  }
  mNcyl2mag++;
  return u;
}

//****************************************************************************
// INPUT:
//   Vector3d xyz == (x,y,z) is the point in Cartesian coordinates,
// OUTPUT:
//  function returns magnetic coordinates(s,theta,phi) of the point xyz
//
Vector3d CStconfig::xyz2mag(const Vector3d &xyz)
{
  Vector3d cyl = xyz.toCylindrical();
  return cyl2mag(cyl);
}

//****************************************************************************
// INPUT:
//   Vector3d cyl == (r,fi,z) is the point in cylindrical coordinates,
// OUTPUT:
//  function returns magnetic coordinates(s,theta,phi) of the point 'cyl'
//
Vector3d CStconfig::cyl2mag(const Vector3d &cyl)
{
  if(!mOK) return Vector3d(0.,0.,0.);
  Vector3d guess;
  Vector3d *Jlast=0;
  double eps = 0.05;
#if 1
  if(ct().isCoordValid() && ct().magLast[0]<1 ) {
    double reff = r(ct().magLast[0]);
    double dr = (reff<0.25*ma)?(ma*0.04):(ma*0.1);  // 2cm or 5cm for W7-X
    if(ct().diffCyl2(cyl) < dr*dr) { // if distance is small then try to use previous results as guess 
      Vector3d ax = ct().axisLast;
      if(!ct().isAxisValid(cyl)) {
        ax = mixcoord2cyl(0,0,cyl[1]);  // magnetic axis
        ct().axisLast = ax;
      }
      double len2=(cyl-ax).abs2();
      if(len2 > square(0.1*ma) ) {
        guess = ct().magLast;   //  use previous result as guess and
        Jlast = ct().Jlast;     //  previous Jacobian calculated at guess
      }
      else eps = 1e-4;
    }
  }
  if(Jlast==0) {
    guess = cyl2guess(cyl,eps); //,0.06); // find guess //std::cerr <<"CStconfig::guess"<<guess<<"\n";
  }
  return cyl2mag(cyl,guess,Jlast);
#else
  if( ct().isCoordValid()
      && (ct().diffCyl2(cyl) < ma*ma*16e-4) // 2cm fo W7X
      && (ct().magLast[0]<1) ) {  // if distance is small then
    guess = ct().magLast;   //  use previous result as guess and
    Jlast = ct().Jlast;   //  previous Jacobian calculated at guess
  }
  else {
    guess = cyl2guess(cyl,0.05); //,0.06); // find guess
  }
  return cyl2mag(cyl,guess,Jlast);
#endif
}

//****************************************************************************
static double findIntersection(const Vector3d &cyl, const Vector3d &ax, const Vector3d &p1, const Vector3d &p2)
{
  if(p1==p2) {
    return (p1-ax).abs(); 
  }
  Vector3d n(p2-p1);   n[1]=0; 
  Vector3d d(cyl-ax);  d[1]=0; 
  n ^= Vector3d(0,1,0); // n is the normal to the line p2-p1. (p1-r)*n = 0 is this line
  d.norm(); // d is the direction of the line cyl-ax. r=ax + d*l is this line 
  // intersection of the line r=ax + d*l with the plane (p1-r)*n = 0
  return fabs( ((p1-ax)*n) / (d*n) ); // intersection point
}

//****************************************************************************
// approximation to atan2(y,x)
// 
inline static double fastAtan2(double y, double x)
{  //return atan2(y,x);
   const double pi2 = 3.1415926535897932384626433832795/2;
   double absy  = fabs(y)+1e-14;  // kludge to prevent 0/0 condition
   double angle = (x>=0)?(absy/(x+absy) ):(  (x+x-absy)/(x-absy)  );
   return pi2*( y<0?-angle:angle );  
}

//****************************************************************************
// Find guess using geometry, 2D-iteration in R-Z plane
// INPUT:
//   Vector3d cyl == (r,fi,z) is the point in cylindrical coordinates,
// OUTPUT:
//  function returns guess value of magnetic coordinates(s,theta,phi) of the point 'cyl'
//
// see also
//  J.A. Ford (1995), Improved Algorithms of Illinois-type for the Numerical 
//  Solution of Nonlinear Equations, Technical Report CSM-257, University of Essex, 1995
//  http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.53.8676
//
Vector3d CStconfig::cyl2guess(const Vector3d &cyl,double eps, double smax, bool *inside)
{
  const double epsMin=1e-9;
  double r  = cyl[0];
  double fi = cyl[1];
  double z  = cyl[2];
  double s(0),t(0),phi(fi);   // magnetic coordinates that will be returned
 
  sumDat sdat(*this,0,fi);

  Vector3d ax = mixcoord2cyl(s,t,fi,&phi,&sdat); // find axis coordinates on plane (r,z), fi-cut
  ct().axisLast = ax;
  double r0=ax[0],
         z0=ax[2];
  double len2=(square(r-r0)+square(z-z0));
  if(len2<square(ma*(eps>0.01?0.01:eps)) ) { //near axis??
    if(inside!=0) *inside = true;
    return Vector3d(0,0,phi); 
  }
  if(len2<square(this->r(s1st)) && eps>0.01) eps = 0.001; //near axis??
  
// Find the crosspoint p of the ray (ax-->cyl) with the countour s,
// first theta is iterated for given s, then s adusted by secant method
// using two distances from the magnetic axis to the point cyl and the point p.
// In the crosspoint p the magnetic coordinates are known,
// these coordinates will be used as the guess for root finding by Newton method 

//int method=3;
//int methodS=3;

  int dirTheta = getThetaDirection();
  Vector3d AC = cyl - ax; AC[1] = 0;
  double ac = fastAtan2(AC[2],AC[0]); //  angle on plane (r,z) w.r.t. R-axis
  double len = sqrt(len2);
  double sa=0, La=-len2;
  double sb=s=smax, Lb=0;
  int outsideCount = 0, sSide=0, itrS = 32;  // 2^-30 = 1e-9
  for(int is = 1; is<itrS; ++is) {           // iterate s
    if(is>1) s = (La*sb - Lb*sa) / (La - Lb);   // false_position_method (regula falsi) for s, comment-out this statesment to use secant method.
    sdat.set(*this,s);
  //begin of false_position_method (regula falsi) for finding the theta
        double t1=0,t2=twopi;          // poloidal magnetic angles; left and right ends of the interval
        Vector3d p,p1 = mixcoord2cyl(s,t1,fi,&phi,&sdat);  // use of mixcoord2cyl is important for QPS where fi and phi are stronly different
        Vector3d p2(p1);
        double a0 = fastAtan2(p1[2]-z0,p1[0]-r0);   // reference angle on plane (r,z) to avoid twopi - 0 interface 
        // we need to find t (poloidal angle theta) that corresponds to ap
        //   or solve equation: mod2pi(atan2(z(t)-z0,r(t)-r0)-a0)  - ap = 0
        double ap = mod2piNoWrap(ac-a0);     // angle on plane (r,z), w.r.t. a0  
        double a1=dirTheta>0?0    :twopi,   
               a2=dirTheta>0?twopi:0; 
        a1-=ap;             // values at the left and right ends of the interval
        a2-=ap;  
        int itrTheta = 34;  // twopi/2^32 = 1.5e-9
        int side=0;
        if(is==1) t = ap;   // theta = alpha is initial guess
        double thetaEps = mmin(0.1,mmax(eps/sqrt(s),2*epsMin)); // eps is the relative accuracy for length, transform to angular measure
        double sEps     = mmin(0.1,2*thetaEps);                 // eps = dr/a = ds/(2*sqrt(s)); ds/s=2eps/sqrt(s)
        double len2Eps  = 2*eps*ma*len;                         // Ls = Len^2-len^2 = (Len-len)(Len+len) < (eps*ma)*2*len
        if(ap<thetaEps||twopi-ap<thetaEps) {  // if close to root
          t = 0;   //std::cerr<<"t="<<ap<<" ";
        }
        else {
          for(int it = 1; it<itrTheta; it++) {     // iterate Theta by regula falsi, it converges as 1.442
//if(it>5) std::cerr<<" it"<<it;
            if(it>1||t==0||t==twopi) { // if (it==1) then use t from previous s-iteration, but t must not be 0 or twopi 
              t = (a1*t2-a2*t1)/(a1-a2);    
            }
//if(t<epsMin||twopi-t<epsMin) std::cerr<<"t==0,2pi"<<t;
            if(fabs(t1-t2)<thetaEps && fabs(a1-a2)<thetaEps) break; 
            p = mixcoord2cyl(s,t,fi,&phi,&sdat); 
            double at = mod2piNoWrap(fastAtan2(p[2]-z0,p[0]-r0)-a0)  - ap;
            if(fabs(at)<thetaEps) { // at is small enough
              p1=p2=p;  
              break;
            }
            else if(at*a2>0) {
              if(side==1) {
                double f0=at/a2;
                double f1=(a1!=0)?(a2/a1):-1;
                double g=1-f0/(1-f1);
                a1 *= (g>0)?g:0.5;
              }
              side = 1;
              a2 = at; t2 = t;  p2 = p;
            }
            else if(a1*at>0) {
              if(side==-1) {
                double f0=at/a1;
                double f1=(a2!=0)?(a1/a2):-1;
                double g=1-f0/(1-f1);
                a2 *= (g>0)?g:0.5;
              }
              side = -1;
              a1 = at; t1 = t;  p1 = p;
            }
          }
        }
        double Len  = findIntersection(cyl, ax, p1, p2); // distance from ax to s
        double Len2 = Len*Len;
  //end of false_position_method (regula falsi) for finding the theta

    double s1 = s;       // save
    s  = s1*len/Len;     // set new s by secant method
    if(s1>=smax && s>smax) outsideCount++;
    if(outsideCount>0) break;  //// if(outsideCount>0 && inside!=0 ) break;
    double Ls = Len2 - len2;
    if(fabs(s-s1)<s*sEps||fabs(Ls)<len2Eps) {
      break;
    }
    if(is==1) {
      Lb = Ls;
      continue; // right end is found, continue with regula falsi
    }
// false_position_method (regula falsi) for s; adjust interval
    if (Ls * Lb> 0) {
      if(side==1) {
        double f0=Ls/Lb;
        double f1=(La!=0)?(Lb/La):-1;
        double g=1-f0/(1-f1);
        La *= (g>0)?g:0.5;
      }
      sSide = 1;
      sb = s1; Lb = Ls;
    }
    else if (La * Ls > 0) {
      if(side==-1) {
        double f0=Ls/La;
        double f1=(Lb!=0)?(La/Lb):-1;
        double g=1-f0/(1-f1);
        Lb *= (g>0)?g:0.5;
      }
      sSide = -1;
      sa = s1;  La = Ls;
    }
  }
  s = mmin(smax,s);  
  if(inside!=0) *inside = (outsideCount==0);
  mNGuessCalls++;  
  return Vector3d(s,t,phi); 
}

//****************************************************************************
// create border -- 's'-contour at cylindrical angle 'ficyl'
// @param borderSize -- # of point in contour
// @param ficyl -- cylindrical angle
// @return false if errors
bool CStconfig::createBorder(double s, double ficyl,int borderSize)
{
  if(!mOK) return false;

  if(!mBorder.init(borderSize)) {
    return false;
  }
  if(!mBorderMagCoord.init(borderSize)) {
    mBorder.clear();
    return false;
  }
  mBorderMagCoord.init(borderSize);
  Vector3d * ptrBorder = mBorder.array();
  Vector3d * ptrBorderMag = mBorderMagCoord.array();

  double phi, dtheta = twopi/(borderSize-1);
  int ka;
  mBorderRmin = mBorderZmin =  1e15;
  mBorderRmax = mBorderZmax = -1e15;
  sumDat sdat(*this,s,ficyl);
  for(ka=0; ka<borderSize-1; ka++) {
    ptrBorder[ka] = mixcoord2cyl(s,ka*dtheta,ficyl,&phi,&sdat);
    mBorderRmin = mmin(mBorderRmin,ptrBorder[ka][0]); // sizes of s-surface
    mBorderRmax = mmax(mBorderRmax,ptrBorder[ka][0]);
    mBorderZmin = mmin(mBorderZmin,ptrBorder[ka][2]);
    mBorderZmax = mmax(mBorderZmax,ptrBorder[ka][2]);
    ptrBorderMag[ka].set(s,ka*dtheta,phi); // magnetic coordinates
  }
  ptrBorder[borderSize-1]     = ptrBorder[0]; // periodicity
  ptrBorderMag[borderSize-1] = ptrBorderMag[0];
  return true;
}

//****************************************************************************
// 'CStconfig::createBorder' must be callled first!
// INPUT:
//   Vector3d cyl == (r,fi,z) is the point in cylindrical coordinates,
// RETURN:
//  function returns true if point lies inside contour 'mBorder'
//
bool CStconfig::cylIsInsideBorder(const Vector3d &cyl,double * distance, Vector3d * point, Vector3d * magCoord) const
{
  if(distance) *distance = 1.e30;
  if(!mOK) return false;
  if(mBorder.empty()) return false;

  if(distance==NULL)
    if(cyl[0]<mBorderRmin||cyl[0]>mBorderRmax||
       cyl[2]<mBorderZmin||cyl[2]>mBorderZmax ) return false;

  const Vector3d * ptrBorder = mBorder.constArray();
  const Vector3d * ptrBorderMag = mBorderMagCoord.constArray();

  int ka,nc=0,sh,nsh;
  double ua=ptrBorder[0][0]-cyl[0];
  double va=ptrBorder[0][2]-cyl[2];

  Vector3d b, a(ua,0,va);
  double blen,alen=ua*ua+va*va;
  double distMin=alen;
  Vector3d pnt = a;  // nearest point on  the border
  Vector3d bza = ptrBorderMag[0];  // nearest point on the border
  Vector3d bz  = bza;
  double dtheta = ptrBorderMag[1][1] - ptrBorderMag[0][1];

  sh = (va<0)?-1:1;
  for(ka=0; ka<int(mBorder.size())-1; ka++) {
    int kb = ka+1;
    double ub=ptrBorder[kb][0]-cyl[0];
    double vb=ptrBorder[kb][2]-cyl[2];
    Vector3d bzb = ptrBorderMag[kb];
    nsh = (vb<0)?-1:1;
    if(sh!=nsh) {
      sh = nsh;
      if(ua>0&&ub>0) nc++;       // if true then there is an intersection with border
      else if(!(ua<0&&ub<0)) if((va*ub-vb*ua)*(va-vb)>0) nc++; // compute intersection
//    else if(!(ua<0&&ub<0)) if((va*ub-vb*ua)/(va-vb)>0) nc++; // compute intersection
    }
// do we need the distance from cyl to Border?
    if(distance) {
      // triangle with vertixes (0, a, b)
      // here we find distance form 0-vertex to opposite side of the triangle
      b.set(ub,0.,vb);
      blen = ub*ub+vb*vb;
      if(blen<distMin) {
        distMin = blen;
        pnt = b;
        bz = bzb;
      }
      Vector3d c = b-a;
      if(a*c<0&&b*c>0) {             // all triangle angles are less then pi/2
        Vector3d n(-c.z(),0,c.x());  // normal to c
        n.norm();
        double h = square(a*n);      // square of height
        if(h<distMin) {
          distMin = h;
          if(point||magCoord) {
            int dir = (a^b)[1]<0?-1:1;
            pnt = sqrt(h)*n*dir;
            double p = (a-pnt).abs2()/c.abs2();
            bz[1]=bza[1] + dtheta*sqrt(p); // linear interpolation on theta
            double dphi = bzb[2]-bza[2];
            if(fabs(dphi)<pi)  bz[2] = bza[2]+dphi*sqrt(p); // linear interpolation on phi;
          }
        }
      }
      a  = b;
      bza = bzb;
      alen = blen;
    }
    va = vb;
    ua = ub;
  }

  nc %= 2;
  if(point) {
    *point = pnt;
    *point += cyl;
  }
  if(magCoord) *magCoord = bz;
  if(distance) *distance = (nc==1)?-sqrt(distMin):sqrt(distMin);
  return nc==1?true:false; // if nc==1 then point lies inside border
}

//****************************************************************************
// The method returns the pointer to the 2d-array of vertexes of LCMS.
// @return the pointer to the 2d-array of vertexes of LCMS, the vertexes are in cartesian coordinates.
const CArray2d<Vector3d> * CStconfig::LCMS() 
{
  lcms.create(this);  // will be created if needed
  return lcms.vertex.isOK() ? &lcms.vertex : 0;
} 

//****************************************************************************
// The method returns the pointer to the 2d-array of the magnetic field value at the LCMS.
const CArray2d<Vector3d> * CStconfig::BLCMS()  
{
  lcms.create(this);  // will be created if needed
  return lcms.Bfield.isOK() ? &lcms.Bfield : 0;
}

//****************************************************************************
// Create Last Closed Magnetic Surface
bool CStconfig::createLCMS()
{
  if(!mOK) return false;
  clock_t tStart = clock();  //timing
  TESTPRINT( std::cerr <<"CStconfig: Calculating LCMS..." ) 

  volatile ctStack::saveState state(ct);

  double  truncSav = truncate();
  if(truncate()<5.e-6&&!tokamakConfig) truncate(5.e-6);

  double s   = tokamakConfig?1.0:1.02;
  int nTheta = tokamakConfig?200:100;
  int nPhi   = 401; 

  lcms.mshift = 0.005/pi;
  lcms.mdPhi = 0;
  lcms.mRmin = lcms.mZmin =  1e15;
  lcms.mRmax = lcms.mZmax = -1e15;
  double dtheta = twopi/(nTheta-1);
  double dPhi =  twopi;
  int nPhiOnPeriod=0; // number of points on period
  if(mNp>1) {
    nPhiOnPeriod = nPhi/mNp; // mNp is # of periods
    if(nPhiOnPeriod<20) nPhiOnPeriod=20;
    nPhi = nPhiOnPeriod*mNp+1;
  }
  if(!lcms.vertex.resize(nTheta,nPhi)) {
   if(truncSav != truncate()) truncate(truncSav); // restore truncation level
   return false;
  }
  if(!lcms.Bfield.resize(nTheta,nPhi)) {
    lcms.vertex.clear(); 
  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
    return false;
  }
  dPhi /= (nPhi-1);
  int i,j,k;
  if(mNp>1) {  // here we use periodicity
    for(i=0; i<nTheta-1; i++) {
      double theta = i*dtheta + lcms.mshift;
      for(j=0; j<nPhiOnPeriod; j++) {
        double phi = j*dPhi + lcms.mshift;
        double phib; //toroiodal angle in magnetic coordinate
        Vector3d c = mixcoord2cyl(s,theta,phi,&phib); // point (r,fi,z) on surface 's'
        lcms.mRmin = mmin(lcms.mRmin,c[0]); // sizes of LCMS
        lcms.mRmax = mmax(lcms.mRmax,c[0]);
        lcms.mZmin = mmin(lcms.mZmin,c[2]);
        lcms.mZmax = mmax(lcms.mZmax,c[2]);
        lcms.vertex(i,j)=c;
        Vector3d mag(s,theta,phib);
        //NJacobianL(mag);
        //lcms.Bfield(i,j) = getmodB() ;//getBcyl();
        lcms.Bfield(i,j) = B(mag);
        for(k=1; k<mNp; k++) {
          lcms.vertex(i,j+k*nPhiOnPeriod).set(c[0],c[1]+k*mPeriod,c[2]);
          lcms.Bfield(i,j+k*nPhiOnPeriod) = lcms.Bfield(i,j);
        }
      }
    }
    for(j=0; j<nPhi-1; j++) {
      lcms.vertex(nTheta-1,j) = lcms.vertex(0,j);
      lcms.Bfield(nTheta-1,j) = lcms.Bfield(0,j);
    }
    for(i=0; i<nTheta; i++) { // cyl-->xyz transformation
      lcms.vertex(i,nPhi-1) = lcms.vertex(i,0);
      lcms.Bfield(i,nPhi-1) = lcms.Bfield(i,0);
      for(j=0; j<nPhi; j++) {
        double r = lcms.vertex(i,j)[0];
        double f = lcms.vertex(i,j)[1];
        lcms.vertex(i,j)[0] = r*::cos(f);
        lcms.vertex(i,j)[1] = r*::sin(f);
      }
    }
  }
  else {
    for(i=0; i<nTheta-1; i++) {
      double theta = i*dtheta + lcms.mshift;
      double phib; //toroidal angle in magnetic coordinates
      Vector3d c;  // cylindrical coordinates
      double Bmod(0);
      if(tokamakConfig) {
        c = mixcoord2cyl(s,theta,0,&phib); // point (r,fi,z) on surface 's'
        Vector3d mag(s,theta,0);
        Bmod = B(mag);
      }
      for(j=0; j<nPhi-1; j++) {
        double phi = j*dPhi+lcms.mshift;
        if(!tokamakConfig) {
          c = mixcoord2cyl(s,theta,phi,&phib); // point (r,fi,z) on surface 's'
          Vector3d mag(s,theta,phib);
  //        NJacobianL(mag);
  //        lcms.Bfield(i,j) = getmodB(); //getBcyl();
          lcms.Bfield(i,j) = B(mag);
        }
        else {
          c[1] = phi;
          lcms.Bfield(i,j) = Bmod;
        }
        lcms.mRmin = mmin(lcms.mRmin,c[0]);
        lcms.mRmax = mmax(lcms.mRmax,c[0]);
        lcms.mZmin = mmin(lcms.mZmin,c[2]);
        lcms.mZmax = mmax(lcms.mZmax,c[2]);
        lcms.vertex(i,j).set(c[0]*::cos(c[1]),c[0]*::sin(c[1]),c[2]);
      }
    }
    for(j=0; j<nPhi-1; j++) {
      lcms.vertex(nTheta-1,j) = lcms.vertex(0,j);
      lcms.Bfield(nTheta-1,j) = lcms.Bfield(0,j);
    }
    for(i=0; i<nTheta; i++) {
      lcms.vertex(i,nPhi-1)   = lcms.vertex(i,0);
      lcms.Bfield(i,nPhi-1)   = lcms.Bfield(i,0);
    }
  }
  lcms.mdPhi = dPhi;
  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
  TESTPRINT( std::cerr <<"time="<<double(clock()-tStart)/CLOCKS_PER_SEC<<std::endl );
  return true;
}

//****************************************************************************
// Create Surface s
bool CStconfig::createSurface(double s,int nTheta,int nPhi,CArray2d<Vector3d> &vertex, CArray2d<Vector3d> *Bfld, bool cylPhi)
{
  if(!mOK) return false;
  clock_t tStart = clock();  //timing
  TESTPRINT( std::cerr <<"CStconfig: createSurface..." ) 

  volatile ctStack::saveState state(ct);
  double  truncSav = truncate();
  if(truncate()<5.e-6&&!tokamakConfig) truncate(5.e-6);

  double shift = 0.005/pi;
  double dtheta = twopi/(nTheta-1);
  double dPhi =  twopi;
  int nPhiOnPeriod=0; // number of points on period
  if(mNp>1) {
    nPhiOnPeriod = nPhi/mNp; // mNp is # of periods
    if(nPhiOnPeriod<20) nPhiOnPeriod=20;
    nPhi = nPhiOnPeriod*mNp+1;
  }
  if(!vertex.resize(nTheta,nPhi)) {
    return false;
  }
  if(Bfld)
    if(!(*Bfld).resize(nTheta,nPhi)) {
      vertex.clear(); 
      return false;
    }
  dPhi /= (nPhi-1);
  int i,j,k;
  if(mNp>1) {  // here we use periodicity
    for(i=0; i<nTheta-1; i++) {
      double theta = i*dtheta + shift;
      for(j=0; j<nPhiOnPeriod; j++) {
        double phi = j*dPhi + shift;
        double phib = phi; //toroidal angle in magnetic coordinates
        Vector3d c;
        if(cylPhi)
          c = mixcoord2cyl(s,theta,phi,&phib); // point (r,fi,z) on surface 's'
        else
          c = mag2cyl(s,theta,phi); // point (r,fi,z) on surface 's'
        vertex(i,j)=c;
        if(Bfld) {
          Vector3d mag(s,theta,phib);
          ////NJacobianL(mag);
          (*Bfld)(i,j) = B(mag) ;//getBcyl();
        }
        for(k=1; k<mNp; k++) {
          vertex(i,j+k*nPhiOnPeriod).set(c[0],c[1]+k*mPeriod,c[2]);
          if(Bfld) (*Bfld)(i,j+k*nPhiOnPeriod) = (*Bfld)(i,j);
        }
      }
    }
    for(j=0; j<nPhi-1; j++) {
      vertex(nTheta-1,j) = vertex(0,j);
      if(Bfld) (*Bfld)(nTheta-1,j) = (*Bfld)(0,j);
    }
    for(i=0; i<nTheta; i++) { // cyl-->xyz transformation
      vertex(i,nPhi-1) = vertex(i,0);
      if(Bfld) (*Bfld)(i,nPhi-1) = (*Bfld)(i,0);
      for(j=0; j<nPhi; j++) {
        double r = vertex(i,j)[0];
        double f = vertex(i,j)[1];
        vertex(i,j)[0] = r*::cos(f);
        vertex(i,j)[1] = r*::sin(f);
      }
    }
  }
  else {
    for(i=0; i<nTheta-1; i++) {
      double theta = i*dtheta + shift;
      double phib=0; //toroidal angle in magnetic coordinates
      Vector3d c;  // cylindrical angle
      if(tokamakConfig)
        c = mixcoord2cyl(s,theta,0,&phib); // point (r,fi,z) on surface 's'
      for(j=0; j<nPhi-1; j++) {
        double phi = j*dPhi+shift;
          phib = phi;
        if(!tokamakConfig) {
          if(cylPhi)
            c = mixcoord2cyl(s,theta,phi,&phib); // point (r,fi,z) on surface 's'
          else
            c = mag2cyl(s,theta,phi); // point (r,fi,z) on surface 's'
        }
        else {
          c[1] = phi;
          phib = phi;
        }
        vertex(i,j).set(c[0]*::cos(c[1]),c[0]*::sin(c[1]),c[2]);
        if(Bfld) {
          Vector3d mag(s,theta,phib);
          ////NJacobianL(mag);
          (*Bfld)(i,j) = B(mag) ;//getBcyl();
        }
      }
    }
    for(j=0; j<nPhi-1; j++) {
      vertex(nTheta-1,j) = vertex(0,j);
      if(Bfld) (*Bfld)(nTheta-1,j) = (*Bfld)(0,j);
    }
    for(i=0; i<nTheta; i++) {
      vertex(i,nPhi-1)   = vertex(i,0);
      if(Bfld) (*Bfld)(i,nPhi-1)   = (*Bfld)(i,0);
    }
  }
  if(truncSav != truncate()) truncate(truncSav); // restore truncation level
  TESTPRINT( std::cerr <<"time="<<double(clock()-tStart)/CLOCKS_PER_SEC<<std::endl );
  return true;
}

//****************************************************************************
// The method creates array of normals to the surface Vertexes.
// @param Vertexes -- array of vertexes of a surface.
// @param sign -- direction of normals
// @param normals -- array of normals
// @return normals -- array of normals
void CStconfig::createNormals(const CArray2d<Vector3d> &Vertex, CArray2d<Vector3d> &normal, int sign) const
{
#if 1
  {  // find correct sign of normals
    CStconfig* This = const_cast<CStconfig*>(this);
    volatile ctStack::saveState state(This->ct);
    int i=2,j=2;
    Vector3d p1 = Vertex(i,j);
    Vector3d p0 = This->xyz2mag(p1);    // find magnetic coordinates
    p0 = mag2xyz(0,p0[1],p0[2]);  // p0 is the point on the magnetic axis at the same toroidal angle as point p1
    Vector3d g = p1-p0;           // vector from the mag.axis to the Vertex(i,j)
    Vector3d w = (Vertex(i,  j+1)-Vertex(i,  j-1))^
                 (Vertex(i-1,j  )-Vertex(i+1,j  ));
	  sign *= (g*w<0?-1:1) ;
    sign /= abs(sign);
  }
#else
//find direction of theta - clockwise or counterclockwise?
//we need this to find correct direction of normals
  Vector3d w = Vertex(0,0)^Vertex(1,0); //find direction of theta
  sign *= w[1]<0?-1:1;
  sign /= abs(sign);
#endif
//Create vertex normals, averaging all adjacent faces normals
  int nTheta=Vertex.size1();
  int nPhi  =Vertex.size2();
  normal.resize(nTheta,nPhi);
  for(int i=0;i<nTheta;i++)
//???    for(int j=1;j<nPhi-1;j++) {
    for(int j=0;j<nPhi;j++) {
      int ip1=(i+1)         %(nTheta-1);
      int im1=(i-1+nTheta-1)%(nTheta-1);
      int jp1=(j+1)         %(nPhi-1);
      int jm1=(j-1+nPhi-1)  %(nPhi-1);
      normal(i,j) = (Vertex(i,  jp1)-Vertex(i,  jm1))^
                    (Vertex(im1,j  )-Vertex(ip1,j  ));
      normal(i,j)*= sign;
      normal(i,j).norm();  // normalize
    }
}

//*********************************************************************
// INPUT:
//   Vector3d xyz == (x,y,z) is the point in Cartesian coordinates,
// OUTPUT:
//  function returns 'true' if the point lies inside LCMS
// note: tracing the LCMS
bool CStconfig::xyzIsInside1(Vector3d &xyz) const
{
  double rc = sqrt(xyz.x()*xyz.x()+xyz.y()*xyz.y());
  double fi = mod2pi(atan2(xyz.y(),xyz.x()));
  double zc = xyz.z();
  Vector3d cyl(rc,fi,zc);
  return CStconfig::cylIsInside1(cyl);
}

//*********************************************************************
// INPUT:
//   Vector3d cyl is a point in cylindrical coordinates,
// OUTPUT:
//  function returns 'true' if the point lies inside LCMS
// note: tracing the LCMS
bool CStconfig::cylIsInside1(Vector3d &cyl) const
{  
  const_cast<CStconfig*>(this)->lcms.create(this);  // will be created if needed

  const CArray2d<Vector3d> & vertex = lcms.vertex;
  if(!vertex.isOK()) return false;

  double rc = cyl[0];
  double fi = mod2pi(cyl[1]);
  double zc = cyl[2];

  if( rc<lcms.mRmin||rc>lcms.mRmax||zc<lcms.mZmin||zc>lcms.mZmax ) return false;

  Vector3d xyz(rc*::cos(fi),rc*::sin(fi),zc);

// find segment
  int i;
  int nfi=vertex.size2()-1;
  int j=int((fi-lcms.mshift)/lcms.mdPhi+nfi+0.5);
  j  %= nfi;
  int jm1=j-1;
  int jp1=j+1;

  int hitCount = 0;
  for(i=0; i<vertex.size1()-1; i++)
    for(j=jm1; j<jp1; j++) {
      int k  = (j+nfi  )%nfi;
      int k1 = (j+nfi+1)%nfi;
      Vector3d r[3];                       // vertexes of a triangle
      r[0] = vertex(i,  k  );
      r[1] = vertex(i+1,k  );
      r[2] = vertex(i+1,k1);
      hitCount += isTriangleHitted(r,xyz); // check for intersection with triangle (r0,r1,r2)
      r[1] = vertex(i,  k1);
      hitCount += isTriangleHitted(r,xyz);
    }
   hitCount %= 2;
   return (hitCount==1);
}

//*********************************************************************
// INPUT:
//   Vector3d xyz == (x,y,z) is the point in Cartesian coordinates,
// OUTPUT:
//  function returns 'true' if the point lies inside LCMS
bool CStconfig::xyzIsInside(const Vector3d &xyz,double * distance) const
{
  double rc = xyz.x()*xyz.x()+xyz.y()*xyz.y();
  double zc = xyz.z();
  double fi = mod2pi(atan2(xyz.y(),xyz.x()));
  Vector3d cyl(sqrt(rc),fi,zc);
  return CStconfig::cylIsInside(cyl,distance);
}

//*********************************************************************
// INPUT:
//   Vector3d cyl is a point in cylindrical coordinates,
// OUTPUT:
//  function returns 'true' if the point lies inside LCMS
// note: counting of intersections with the countor of LCMS
//       this function works faster then cylIsInside1
//
bool CStconfig::cylIsInside(const Vector3d &cyl, double * distance) const
{
  const_cast<CStconfig*>(this)->lcms.create(this);   // will be created if needed

  if(distance) *distance = 1.e30;

  double rc = cyl[0];
  double zc = cyl[2];

  if(distance==NULL)
    if( rc>lcms.mRmax||rc<lcms.mRmin||zc<lcms.mZmin||zc>lcms.mZmax ) return false;

  double fi = modPeriod(cyl[1]-lcms.mshift);
  Vector3d cylWrk(rc,fi,zc);

  return cylIsInsideEx(cylWrk,lcms.vertex,lcms.mdPhi,distance);
}

//*********************************************************************
//
//  fi=cyl[1] must be given on a period, see cylIsInside()
//
bool CStconfig::cylIsInsideEx(const Vector3d &cyl,const CArray2d<Vector3d> &vertex, double dPhi,double * distance) const
{
  if(!vertex.isOK()) return false;

  double rc = cyl[0];
  double fi = cyl[1];
  double zc = cyl[2];

  // find segment in toroidal direction
  int nfi=vertex.size2()-1;
  double p = fi/dPhi;
  int j=int(p);
  int jp1=j+1;
  p -= j;

  int bSize=vertex.size1()-1;

// linear interpolation
  const Vector3d * vrtx = vertex[0]; // 1st index is the poloidal direction
  Vector3d v = vrtx[j];
  double r1 = sqrt(v.x()*v.x()+v.y()*v.y());
  double z1 = v.z();
  v = vrtx[jp1];
  double r2 = sqrt(v.x()*v.x()+v.y()*v.y());
  double z2 = v.z();
  r1 += (r2-r1)*p;   // interpolate on fi, toroidal direction
  z1 += (z2-z1)*p;
  double r0 = r1;
  double z0 = z1;
//
  int nc=0,sh,nsh;
  double ub,ua=r1-rc;
  double vb,va=z1-zc;

  Vector3d b, a(ua,va,0);
  double blen,alen=ua*ua+va*va;
  double distMin=alen;

  sh = (va<0)?-1:1;

  for(int i=1; i<bSize; i++) {
    vrtx = vertex[i];
    v  = vrtx[j];
    r1 = sqrt(v.x()*v.x()+v.y()*v.y());
    z1 = v.z();
    v  = vrtx[jp1];
    r2 = sqrt(v.x()*v.x()+v.y()*v.y());
    z2 = v.z();
    r1 += (r2-r1)*p;
    z1 += (z2-z1)*p;
    ub=r1-rc;
    vb=z1-zc;
    nsh = (vb<0)?-1:1;
    if(sh!=nsh) {
      sh = nsh;
      if(ua>0&&ub>0) nc++;       // if true then there is an intersection with border
      else if(!(ua<0&&ub<0)) if((va*ub-vb*ua)*(va-vb)>0) nc++; // compute intersection
    }
// do we need the distance from cyl to LCMS?
    if(distance) {
      // triangle with vertixes (0, a, b)
      // here we find distance form 0-vertex to opposite side of the triangle
      b.set(ub,vb,0.);
      blen = ub*ub+vb*vb;
      distMin = mmin(distMin,blen);
      Vector3d c = b-a;
      if(a*c<0&&b*c>0) {           // all triangle angles are less then pi/2
        Vector3d n(-c.y(),c.x(),0.);  // normal to c
        n.norm();
        double h = square(a*n);    // square of height
        distMin = mmin(distMin,h);
      }
      a  = b;
      alen = blen;
    }
    va = vb;
    ua = ub;
  }
// test last segment
  ub=r0-rc;
  vb=z0-zc;
  nsh = (vb<0)?-1:1;
  if(sh!=nsh) {
    sh = nsh;
    if(ua>0&&ub>0) nc++;       // if true then there is an intersection with border
    else if(!(ua<0&&ub<0)) if((va*ub-vb*ua)*(va-vb)>0) nc++; // compute intersection
  }
  nc %= 2;

  if(distance) {
    // triangle with vertixes (0, a, b)
    // here we find distance form 0-vertex to opposite side of the triangle
    b.set(ub,vb,0.);
    blen = ub*ub+vb*vb;
    distMin = mmin(distMin,blen);
    Vector3d c = b-a;
    if(a*c<0&&b*c>0) {           // all triangle angles are less then pi/2
      Vector3d n(-c.y(),c.x(),0.);  // normal to c
      n.norm();
      double h = square(a*n);    // square of the height
      distMin = mmin(distMin,h);
    }
    *distance = (nc==1)?-sqrt(distMin):sqrt(distMin);
  }
  return (nc==1); // if nc==1 then point lies inside border
}

//*********************************************************************
// Test intersection of the ray with the triangle
// The ray is R(t)=rayR0+(0,0,1)*t  with t>0
// INPUT
//  r -- vertices of the triangle
//  rayR0 -- origin of the ray
// OUTPUT
//  function returns 1 if the ray  hits the triangle
int CStconfig::isTriangleHitted(const Vector3d *r,const Vector3d &rayR0) const
{
  Vector3d n = (r[1]-r[0])^(r[2]-r[0]);     // normal to the triangle
  Vector3d rayRd(0,0,1); // only z-direction
  Vector3d ri;
  double vd = n*rayRd;           // dot product
  if(fabs(vd)<1e-13) return 0;   // the ray parallels the triangle
  double t = (n*(r[0]-rayR0))/vd;
  if(t<0) return 0;
  ri = rayR0+rayRd*t; // there is an intersection with the plane of the triangle
  return isInsideTriangle(r,n,ri);
}

//*********************************************************************
//  is intersection point within triangle
// INPUT
//  r -- vertices of triangle
//  n -- normal to the plane of the triangle
//  ri -- intersection
// OUTPUT
//  function returns 'true' if intersection point lies inside triangle
int CStconfig::isInsideTriangle(const Vector3d *r,const Vector3d &n, Vector3d &ri) const
{
//n = (r[1]-r[0])^(r[2]-r[0]);     // normal to the triangle (r0,r1,r2)
  int k;
  Vector3d R[4];
  for(k=0; k<3; k++) {R[k]=r[k]; if(r[k]==ri) return 1;} //return if intersection point coincides with vertex
  int m = n.maxComponent();  // find dominant coordinate in normal vector
  for(k=0; k<3; k++) {
    R[k] -= ri;
    R[k].delComponent(m);             // throw away the dominant coordinate
    if(R[k].y()==0.) R[k].y()+=1e-13; // if vertex lies on test ray
  }
  R[3]=R[0];
  int nc=0,sh,nsh;
  double ub,ua=R[0].x();
  double vb,va=R[0].y();
  sh = (va<0)?-1:1;
  for(k=1; k<=3; k++) {
    ub=R[k].x();
    vb=R[k].y();
    nsh = (vb<0)?-1:1;
    if(sh!=nsh) {
      sh = nsh;
      if(ua>0&&ub>0) nc++;       // if true then there is an intersection
      else if(!(ua<0&&ub<0)) if((va*ub-vb*ua)*(va-vb)>0) nc++; // compute intersection
//    else if(!(ua<0&&ub<0)) if((va*ub-vb*ua)/(va-vb)>0) nc++; // compute intersection
    }
    va = vb;
    ua = ub;
  }
  nc %= 2;
  return nc; // intersection point lies inside if 1
}

//****************************************************************************
// find direction of theta(poloidal angle)
//  signof(grads
int CStconfig::getThetaDirection()
{
  if(mDirTheta==0) signOfJacobianAndAngles();
  return mDirTheta;
    //////Vector3d ax = mixcoord2cyl(0,   0,   0); // axis position in cyl. coord.
    //////Vector3d p  = mixcoord2cyl(0.1, 0,   0); // p  = (s=0.1,th=0,   ficyl=0)
    //////Vector3d p1 = mixcoord2cyl(0.1, 0.03,0); // p1 = (s=0.1,th=0.03,ficyl=0)
    //////p  -= ax;
    //////p1 -= ax;
    //////p [1] = 0;
    //////p1[1] = 0;
    //////p = p^p1; // == p(th=0)^p(th=0.03)
    //////mDirTheta = p[1]<0?1:-1; //direction of theta - clockwise or counterclockwise?
}

//****************************************************************************
// find direction of phi(toroidal angle) w.r.t cylindrical angle
int CStconfig::getPhiDirection()
{
  if(mDirPhi==0) signOfJacobianAndAngles();
  return mDirPhi;
    //////double dfi_dphi;
    //////fi(0.,0.,0.,dfi_dphi);     // get dfi/dphi in point (s=0,theta=0,phi=0)
    //////mDirPhi = dfi_dphi>0?1:-1; //direction of phi
}

//****************************************************************************
// @return the sign of jacobian in (s,theta,phi)-coordinates
int CStconfig::getJacobianSign()
{
  if(mSignJac_s==0) signOfJacobianAndAngles();
  return mSignJac_s;
}

//****************************************************************************
//find the sign of jacobian
void CStconfig::signOfJacobianAndAngles()
{
  if(!mOK) return;
  volatile CStconfig::ctStack::saveState state(ct);
  Vector3d p(0.25,0,0);       //(s,theta,phi)
  NJacobian1L(p);
  Vector3d F (ct().Jmatr[0]); // (R,fi,Z)
  Vector3d Fs(ct().Jmatr[1]); // d(R,fi,Z)/ds   -- partial derivative on s
  Vector3d Ft(ct().Jmatr[2]); // d(R,fi,Z)/dtheta
  Vector3d Fp(ct().Jmatr[3]); // d(R,fi,Z)/dphi
  double R   = F[0];       // take into account cylindrical coordinate system
  Vector3d gS = R*(Ft^Fp); // gS is proportional to gradS, where ^ is a cross product sign

  double jac = Fs*gS;      // jacobian
  ////////disabled 06Nov2011
  ////////gS /= jac;               // gS is true gradS now
  ////////// dalfa_dtheta is proportional to dalfa/dtheta, where alfa = atan((Z-Zaxis)/(R-Raxis))
  ////////double dalfa_dtheta = Ft[2]*gS[0] - Ft[0]*gS[2]; // Ft = d(R,fi,Z)/dtheta
  ////////mDirTheta  = sign(dalfa_dtheta); // direction of theta
  mDirTheta  = sign(-Fp[1]/jac); // direction of theta
  mDirPhi    = sign(Fp[1]);        // Fp[1]=dfi/dphi gives direction of phi w.r.t cyl angle
  mSignJac_s = sign(jac);          // sign of jacobian in (s,theta,phi)-coordinates
  return;
}

//****************************************************************************
//****************************************************************************
//****************************************************************************
// Quadratic interpolation in u. x and y are input arrays
// on return y(u)
double CStconfig::interp2(double u, const CArray1d<double> &x, const CArray1d<double> &y) const
{
  int N = x.size();
  int L = x.bsearch(u);
  if(L+2>=N) L=N-3;
  double y0 = y[L  ];
  double y1 = y[L+1];
  double y2 = y[L+2];
  double x0 = x[L  ];
  double x1 = x[L+1];
  double x2 = x[L+2];
  double c1=(y1-y0)/(x1-x0);
  double c2=(y2-y1)/(x2-x1);
  c2=(c2-c1)/(x2-x0);
  return  y0 + (c1 + c2*(u-x1))*(u-x0);
}

//****************************************************************************
// Linear interpolation in u. x and y are input arrays
// on return y(u)
double CStconfig::interp1(double u, const CArray1d<double> &x, const CArray1d<double> &y) const
{
  int L = x.bsearch(u);
  int R = L+1;
  double yp= (y[R]-y[L])/(x[R]-x[L]); // derivative in point u
  return(y[L] + yp*(u-x[L]));
}


#if 0
//****************************************************************************
//****************************************************************************
//****************************************************************************
// perform a binary search of the segment in which u lies. x input array
int CStconfig::bsearch(double u, const CArray1d<double> &xa) const
{
  int k;
  int L = 0;    // Left bracket
  int R = xa.size()-1;  // Right bracket
  const double * x = xa.constArray();

  if(u<=x[0]) return L;
  else if(u>=x[R]) L=R-1;
  else
    while(L+1 < R) {  // bisection
      k = (L + R)/2;
      if(u < x[k]) R = k;
      else L = k;
    }
  return L;
}

//****************************************************************************
// perform a binary search of the segment in which u lies. x input array
static int bsearch(double u, const double *x, int N)
{
  int k;
  int L = 0;    // Left bracket
  int R = N-1;  // Right bracket

  if(u<=x[0]) return L;
  else if(u>=x[R]) L=R-1;
  else
    while(L+1 < R) {  // bisection
      k = (L + R)/2;
      if(u < x[k]) R = k;
      else L = k;
    }
  return L;
}

//****************************************************************************
// Quadratic interpolation in u. x and y are input arrays
// on return y(u)
static double interp2(double u, const double *x, const double *y, int N)
{
  int L = bsearch(u,x,N);
  if(L+2>=N) L=N-3;
  double y0 = y[L  ];
  double y1 = y[L+1];
  double y2 = y[L+2];
  double x0 = x[L  ];
  double x1 = x[L+1];
  double x2 = x[L+2];
  double c1=(y1-y0)/(x1-x0);
  double c2=(y2-y1)/(x2-x1);
  c2=(c2-c1)/(x2-x0);
  return  y0 + (c1 + c2*(u-x1))*(u-x0);
}

//****************************************************************************
// Quadratic interpolation in u. x and y are input arrays with sizes N>2
// on return y(u) and yp==y'(u)
static double interp2(double u, double &yp, double *x, double *y, int N)
{
  int L = bsearch(u,x,N);
  if(L+2>=N) L=N-3;
  double y0 = y[L  ];
  double y1 = y[L+1];
  double y2 = y[L+2];
  double x0 = x[L  ];
  double x1 = x[L+1];
  double x2 = x[L+2];
  double c1=(y1-y0)/(x1-x0);
  double c2=(y2-y1)/(x2-x1);
  c2=(c2-c1)/(x2-x0);
  yp = c1 + c2*(u-x1)+c2*(u-x0); // derivative
  return  y0 + (c1 + c2*(u-x1))*(u-x0);
}

//****************************************************************************
// Quadratic interpolation in u. x and y are input arrays with size N>2
// on return y(u) and yp==y'(u)
static double interp2(double u, double &yp, double *x, double *y)
{
  double c1=(y[1]-y[0])/(x[1]-x[0]);
  double c2=(y[2]-y[1])/(x[2]-x[1]);
  c2=(c2-c1)/(x[2]-x[0]);
  yp = c1 + c2*(u-x[1])+c2*(u-x[0]); // derivative
  return  y[0] + (c1 + c2*(u-x[1]))*(u-x[0]);
}

//****************************************************************************
// Quadratic interpolation in u. x and y are input arrays
// on return y(u)
static double interp2(double u, double *x, double *y)
{
  double c1=(y[1]-y[0])/(x[1]-x[0]);
  double c2=(y[2]-y[1])/(x[2]-x[1]);
  c2=(c2-c1)/(x[2]-x[0]);
  return  y[0] + (c1 + c2*(u-x[1]))*(u-x[0]);
}

//****************************************************************************
// Quadratic interpolation in u. x and y are input arrays
// on return y(u)
static double interp2(double u, double *x, double *y, int i0, int i1, int i2)
{
  double c1=(y[i1]-y[i0])/(x[i1]-x[i0]);
  double c2=(y[i2]-y[i1])/(x[i2]-x[i1]);
  c2=(c2-c1)/(x[i2]-x[i0]);
  return  y[i0] + (c1 + c2*(u-x[i1]))*(u-x[i0]);
}

//****************************************************************************
// Linear interpolation in u. x and y are input arrays
// on return y(u) and yp==y'(u)
static double interp(double u, double &yp, const double *x, const double *y, int N)
{
  int L = bsearch(u,x,N);
  int R = L+1;
  yp= (y[R]-y[L])/(x[R]-x[L]); // derivative in point u
  return(y[L] + yp*(u-x[L]));
}

//****************************************************************************
// Linear interpolation in u. x and y are input arrays
// on return y(u)
static double interp(double u, const double *x, const double *y, int N)
{
  int L = bsearch(u,x,N);
  int R = L+1;
  double yp= (y[R]-y[L])/(x[R]-x[L]); // derivative in point u
  return(y[L] + yp*(u-x[L]));
}

//****************************************************************************
// Linear interpolation/extrapolation in u. x and y are input arrays with sizes == 2
// on return y(u) and yp==y'(u)
static inline double interp(double u, double &yp, double *x, double *y)
{
  yp= (y[1]-y[0])/(x[1]-x[0]); // derivative in point u
  return(y[0] + yp*(u-x[0]));
}

#endif


//****************************************************************************
// FluxPoloidal
// Integral(s1,s2) iota(s)*ds
double CStconfig::FluxPoloidal(double s1, double s2) const
{
  double sum=0;
  int n=10;
  n=n+n%2;         //must be even for simpson integration
  double dp=(s2-s1)/n;
  if(dp==0) return 0;
  for(int k=0; k<=n; k++) {  //Simpson integration
    double w = (k==0||k==n)?1:(2*(1+k%2));   // weight for Simpson integration
    sum += w*iota(s1+k*dp);
  }
  return fabs(sum*dp/3);
};

//****************************************************************************
// function returns the average area of the cross section of the surface 's'
double CStconfig::averageCrossSectionArea(double s)
{
  if(s==1.&&mCrossArea>0) return mCrossArea;

  double  truncSav = truncate();
  if(truncate()<1.e-6&&!tokamakConfig) truncate(1.e-6);

  int itor = tokamakConfig?1:100;   // Tokamak ???
  double dfi = mPeriod/itor;        // toroidal angle increment
  double area=0;
  for(int k=0; k<itor; k++)
    area += crossSectionArea(s, k*dfi);
  area/=itor;

  if(truncSav != truncate()) truncate(truncSav); // restore truncation level

  if(s==1.) mCrossArea = area;

  //double a2=0;
  //for(int m=1; m<=mM; m++) {
  //  int lowN = tokamakEQDSK?0:-mN;
  //  for(int n=lowN; n<= mN; n++) 
  //    a2 +=  m*Rmn_sp(s,m,n)*Zmn_sp(s,m,n);
  //}
  // a2 for vmec coordinates does not coincide with that in boozer
  // a2 for vmec coordinates coincides with my calculation
  //a2 = pi*fabs(a2);

  return area;
}

//****************************************************************************
// function returns the area of the cross section of the surface 's'
// at the cylindrical angle 'fi'
double CStconfig::crossSectionArea(double s, double fi)
{
  if(!mOK) return 0;
  int ipol = tokamakConfig?300:200; // number of poloidal points
  double dtheta = twopi/(ipol-1);
  sumDat sdat(*this,0,fi);
  Vector3d ax   = mixcoord2cyl(0,0,fi,0,&sdat);  // magnetic axis
  sdat.set(*this,s);
  Vector3d brd1 = mixcoord2cyl(s,0,fi,0,&sdat);  //1st border point
  double area=0;
  for(int i=1; i<ipol; i++) {
    Vector3d brd2 = mixcoord2cyl(s,i*dtheta,fi,0,&sdat);  //2nd border point
    Vector3d a = (brd1-ax)^(brd2-ax);
    area += a.abs();
    brd1 = brd2;
  }
  return area/2;
}

//***************************************************************
double CStconfig::getdBGradsxyz(const Vector3d &xyz,Vector3d &B,Vector3d &dBdx,Vector3d &dBdy,Vector3d &dBdz,Vector3d &grads )
{

  B = 0;
  dBdx = 0; // dB/dx
  dBdy = 0; // dB/dy
  dBdz = 0; // dB/dz
  grads = 0;

  Vector3d cyl(sqrt(xyz.x()*xyz.x()+xyz.y()*xyz.y()), mod2pi(atan2(xyz.y(),xyz.x())), xyz.z());
  Vector3d b,dbdr,dbdfir,dbdz,gs;           // in cyl. coordinates
  double s = getdBGradscyl(cyl,b,dbdr,dbdfir,dbdz,gs);

  //if(s>1) {
  //  return s;  // return zeros (B and derivatives) if point lies outside
  //}

// transform to Cartesian coordinates
  double r  = cyl[0];
  double fi = cyl[1];
  double cs = ::cos(fi);
  double sn = ::sin(fi);
  double Bx = b[0]*cs-b[1]*sn;  // Bx=Br*cos(fi)-Bfi*sin(fi);
  double By = b[0]*sn+b[1]*cs;  // By=Br*sin(fi)+Bfi*cos(fi);
  double Bz = b[2];
  B.set(Bx,By,Bz);
// grad(s)
  grads.x() = gs[0]*cs-gs[1]*sn;  // gx=gr*cos(fi)-gfi*sin(fi);
  grads.y() = gs[0]*sn+gs[1]*cs;  // gy=gr*sin(fi)+gfi*cos(fi);
  grads.z() = gs[2];

  double dBxdr = dbdr[0]*cs-dbdr[1]*sn;            // dBxdr = dBr/dr*cs-dBfi/dr*sn;
  double dBydr = dbdr[0]*sn+dbdr[1]*cs;            // dBydr = dBr/dr*sn+dBfi/dr*cs;

  double dBxdfir = dbdfir[0]*cs-dbdfir[1]*sn-By/r;  // dBx/dfi/r = dBr/dfi/r*cs-dBfi/dfi/r*sn-By/r;
  double dBydfir = dbdfir[0]*sn+dbdfir[1]*cs+Bx/r;  // dBy/dfi/r = dBr/dfi/r*sn+dBfi/dfi/r*cs+Bx/r;

  dBdx.x() = dBxdr*cs-dBxdfir*sn;      // dBx/dx=dBx/dr*cs-dBx/dfi/r*sn
  dBdy.x() = dBxdr*sn+dBxdfir*cs;      // dBx/dy=dBx/dr*sn+dBx/dfi/r*cs
  dBdz.x() = dbdz[0]*cs-dbdz[1]*sn;    // dBx/dz=dBr/dz*cs-dBfi/dz*sn;

  dBdx.y() = dBydr*cs-dBydfir*sn;      // dBy/dx=dBy/dr*cs-dBy/dfi/r*sn
  dBdy.y() = dBydr*sn+dBydfir*cs;      // dBy/dy=dBy/dr*sn+dBy/dfi/r*cs
  dBdz.y() = dbdz[0]*sn+dbdz[1]*cs;    // dBy/dz=dBr/dz*sn+dBfi/dz*cs;

  dBdx.z() = dbdr[2]*cs-dbdfir[2]*sn;  // dBz/dx=dBz/dr*cos(fi)-dBz/dfi/r*sin(fi);
  dBdy.z() = dbdr[2]*sn+dbdfir[2]*cs;  // dBz/dy=dBz/dr*sin(fi)+dBz/dfi/r*cos(fi);
  dBdz.z() = dbdz[2];                  // dBz/dz

  return s;
}

//***************************************************************
double CStconfig::getdBGradscyl(const Vector3d &cyl,Vector3d &B,Vector3d &dBdr,Vector3d &dBdfir,Vector3d &dBdz,Vector3d &grads )
{
  B     = 0;
  dBdr  = 0; // dB/dr
  dBdfir= 0; // dB/dfi/r
  dBdz  = 0; // dB/dz
  grads = 0;

  volatile ctStack::saveState state(ct);

  Vector3d b  = cyl2mag(cyl);
  double s  = b[0];
  //if(s>1) {
  //  return s;
  //}

  B  = getBcyl();
  grads = getGrads();
  { // central differencing is more accurate than forward differencing
    double d =  ma/512;  // 1./1024; //1mm for W7-X
    double df = d/mR0;   // 1./2048;
    Vector3d dr(d,0,0),dz(0,0,d),dfi(0,df,0);

    Vector3d B1,B2, c1,c2;
    c1 = cyl - dr;
    c2 = cyl + dr;
    B1 = getBcyl(c1);
    B2 = getBcyl(c2);

    dBdr = (B2-B1)/(2*d);

    c1 = cyl - dz;
    c2 = cyl + dz;
    B1 = getBcyl(c1);
    B2 = getBcyl(c2);

    dBdz = (B2-B1)/(2*d);

    c1 = cyl - dfi;
    c2 = cyl + dfi;
    B1 = getBcyl(c1);
    B2 = getBcyl(c2);

    dBdfir = (B2-B1)/(2*df)/cyl[0];
  }
  return s;
}

//***************************************************************
// @return cylindrical coordinates
// 
Vector3d CStconfig::getBcurlBGradscyl(const Vector3d &magCrd,Vector3d &B,Vector3d &curlB,Vector3d &grads)
{
  Vector3d cyl(0.,0.,0.);
  B     = 0.;
  grads = 0.;


  volatile ctStack::saveState state(ct);

  double s  = magCrd[0];
  //if(s>1) {
  //  return cyl;
  //}

  NJacobian(magCrd);

  cyl = getRcyl();    //  (R,fi,Z)
  B   = getBcyl();
  grads = getGrads();
  double jac = getJ();   // Jacobian in coordinates (s,theta,phi)

  Vector3d dBds,dBdt,dBdp;
  Vector3d Fs,Ft,Fp;    // covariant-basis vectors
  getCovaBasis(Fs,Ft,Fp);

  {
    //volatile ctStack::saveState state(ct);
   // for central differencing
    double d_s =  1./512;        // =1mm    for W7-X; a=510mm
    double d_t =  d_s;           // =1mm/a  for W7-X
    double d_p =  d_s*ma/cyl[0]; // =1mm/r  for W7-X
    Vector3d ds(d_s,0,0), dt(0,d_t,0), dp(0,0,d_p);

    Vector3d B1,B2,m1,m2;
    m1 = magCrd - ds;
    m2 = magCrd + ds;
    NJacobian(m1);
    B1 = getBcova();
    NJacobian(m2);
    B2 = getBcova();

    dBds = (B2-B1)/(2*d_s);

    m1 = magCrd - dt;
    m2 = magCrd + dt;
    NJacobian(m1);
    B1 = getBcova();
    NJacobian(m2);
    B2 = getBcova();

    dBdt = (B2-B1)/(2*d_t);

    m1 = magCrd - dp;
    m2 = magCrd + dp;
    NJacobian(m1);
    B1 = getBcova();
    NJacobian(m2);
    B2 = getBcova();

    dBdp = (B2-B1)/(2*d_p);
  }

  // contravariant components of the curl(B)*jacobian
  double jt = dBdp[0] - dBds[2];
  double jp = dBds[1] - dBdt[0];
  double js = dBdt[2] - dBdp[1];  // must be zero, but who knows .....

  curlB = (jt*Ft + jp*Fp)/jac; //+js*Fs/jac;

  return cyl;
}

//***************************************************************
// @return double local shear
// @return cyl, B, grads in cylindrical coordinates
// 
double CStconfig::getLocalShear(const Vector3d &magCrd,Vector3d &B,Vector3d &cyl,Vector3d &grads)
{
  cyl   = 0.;
  B     = 0.;
  grads = 0.;

  volatile ctStack::saveState state(ct);

  double s  = magCrd[0];
  //if(s>1) {
  //  return 0;
  //}

  NJacobian(magCrd);
  Vector3d Fs,Ft,Fp;    // covariant-basis vectors

  cyl = getRcyl();    //  (R,fi,Z)
  B   = getBcyl();
  grads = getGrads();
  double jac = getJ();   // Jacobian in coordinates (s,theta,phi)
  getCovaBasis(Fs,Ft,Fp);
  Vector3d b  = B;
  Vector3d gs = grads;
  b.norm();
  gs.norm();
  Vector3d H = gs^b;

  Vector3d dHds,dHdt,dHdp;
  {
   // for central differencing
    double d_s =  1./512;        // =1mm    for W7-X; a=510mm
    double d_t =  d_s;           // =1mm/a  for W7-X
    double d_p =  d_s*ma/cyl[0]; // =1mm/r  for W7-X
    Vector3d ds(d_s,0,0), dt(0,d_t,0), dp(0,0,d_p);
    Vector3d Fs,Ft,Fp;    // covariant-basis vectors

    Vector3d b,gs,h,H1,H2,m1,m2;
    m1 = magCrd - ds;
    m2 = magCrd + ds;
    NJacobian(m1);
    b  = getBcyl();
    gs = getGrads();
    b.norm();
    gs.norm();
    h = gs^b;
    getCovaBasis(Fs,Ft,Fp);
    H1.set(h*Fs,h*Ft,h*Fp); 
    NJacobian(m2);
    b  = getBcyl();
    gs = getGrads();
    b.norm();
    gs.norm();
    h = gs^b;
    getCovaBasis(Fs,Ft,Fp);
    H2.set(h*Fs,h*Ft,h*Fp); 

    dHds = (H2-H1)/(2*d_s);

    m1 = magCrd - dt;
    m2 = magCrd + dt;
    NJacobian(m1);
    b  = getBcyl();
    gs = getGrads();
    b.norm();
    gs.norm();
    h = gs^b;
    getCovaBasis(Fs,Ft,Fp);
    H1.set(h*Fs,h*Ft,h*Fp); 
    NJacobian(m2);
    b  = getBcyl();
    gs = getGrads();
    b.norm();
    gs.norm();
    h = gs^b;
    getCovaBasis(Fs,Ft,Fp);
    H2.set(h*Fs,h*Ft,h*Fp); 

    dHdt = (H2-H1)/(2*d_t);

    m1 = magCrd - dp;
    m2 = magCrd + dp;
    NJacobian(m1);
    b  = getBcyl();
    gs = getGrads();
    b.norm();
    gs.norm();
    h = gs^b;
    getCovaBasis(Fs,Ft,Fp);
    H1.set(h*Fs,h*Ft,h*Fp); 
    NJacobian(m2);
    b  = getBcyl();
    gs = getGrads();
    b.norm();
    gs.norm();
    h = gs^b;
    getCovaBasis(Fs,Ft,Fp);
    H2.set(h*Fs,h*Ft,h*Fp); 

    dHdp = (H2-H1)/(2*d_p);
  }

  // contravariant components of the curl(H)*jacobian
  double jp = dHds[1] - dHdt[0];
  double js = dHdt[2] - dHdp[1];
  double jt = dHdp[0] - dHds[2];

  Vector3d curlH = (jt*Ft+jp*Fp)/jac + js*Fs/jac;

  double localShear = -(H*curlH);  //  twopi*(H*curlH)

  return localShear;
}

//***************************************************************
double CStconfig::getLocalShear(const Vector3d &magCrd)
{
  Vector3d H  = 0; 
  Vector3d dHdr  = 0; // dB/dr
  Vector3d dHdfir= 0; // dB/dfi/r
  Vector3d dHdz  = 0; // dB/dz

  volatile ctStack::saveState state(ct);
  Vector3d cyl = mag2cyl(magCrd);
  Vector3d b  = getBcyl(cyl);
  Vector3d gs = getGrads();
  b.norm();
  gs.norm();
  H = gs^b;
  { // central differencing is more accurate than forward differencing
    double d =  ma/512;  // 1./1024; //1mm for W7-X
    double df = d/mR0;   // 1./2048;
    Vector3d dr(d,0,0),dz(0,0,d),dfi(0,df,0);

    Vector3d H1,H2, c1,c2, gs,b,h;
    c1 = cyl - dr;
    c2 = cyl + dr;
    b = getBcyl(c1);
    gs = getGrads();
    b.norm();
    gs.norm();
    H1 = gs^b;

    b = getBcyl(c2);
    gs = getGrads();
    b.norm();
    gs.norm();
    H2 = gs^b;

    dHdr = (H2-H1)/(2*d);

    c1 = cyl - dz;
    c2 = cyl + dz;
    b = getBcyl(c1);
    gs = getGrads();
    b.norm();
    gs.norm();
    H1 = gs^b;

    b = getBcyl(c2);
    gs = getGrads();
    b.norm();
    gs.norm();
    H2 = gs^b;

    dHdz = (H2-H1)/(2*d);

    c1 = cyl - dfi;
    c2 = cyl + dfi;
    b = getBcyl(c1);
    gs = getGrads();
    b.norm();
    gs.norm();
    H1 = gs^b;

    b = getBcyl(c2);
    gs = getGrads();
    b.norm();
    gs.norm();
    H2 = gs^b;

    dHdfir = (H2-H1)/(2*df)/cyl[0];
  }

  double curl_r  = dHdfir[2] - dHdz[1];
  double curl_fi = dHdz[0]   - dHdr[2];
  double curl_z  = H[1]/cyl[0] + dHdr[1] - dHdfir[0];

  //double localShr = twopi*( H*Vector3d(curl_r,curl_fi,curl_z) );
  double localShr = -( H*Vector3d(curl_r,curl_fi,curl_z) );

  return localShr;
}





//*********************************************************************
// Find intersection of the Ray with the last surface
// The ray is R(t)=r0+rd*t  with t>0
// It is assumed, that the origin of ray lies outside of the last surface
//
// INPUT:
//   r0, rd  -- origin and direction of the ray in cartesian coordinates
// OUTPUT:
//  entryPoint
//  @return  if the method succeeds, the return value is true; false otherwise.
bool CStconfig::getRayEntryPoint(const Vector3d &r0,const Vector3d &rd,Vector3d &entryPoint)
{
  if(!mOK) return false;
  volatile ctStack::saveState state(ct);
  lcms.create(this);   // will be created if needed
  if(!lcms.vertex.isOK()) return false;

  if(xyzIsInside(r0)) return false;  // return if the origin lies inside LCMS

  double dl = rminor()*0.04;     // step 2cm for W7X
  Vector3d dr = dl*rd/rd.abs();
  double Rmax = lcms.mRmax + 20*dl;
  int N = int(100*Rmax/dl);      // max iterations
  Vector3d r = r0;
  double Rmax2 = square(Rmax);
  double r2 = r.abs2();
  do {                      // cycle while not inside
    double rprev2 = r2;
    r += dr;                // move forward along ray
    r2 = r.abs2();
    if(r2>rprev2&&r2>Rmax2) return false; // return if point moves far outwards
    if(--N < 1) return false; // entry point not found
  } while(!xyzIsInside(r));   // cycle while not inside
 
  bool entryFound = getRayIntersection(r,dr,entryPoint);
  return entryFound;
}

//*********************************************************************
// Find intersection of the Ray with the last surface
// The ray is R(t)=r0+rd*t  with t>0
// It is assumed, that the origin of ray lies outside of the last surface
//
// INPUT:
//   r0, rd  -- origin and direction of the ray in cartesian coordinates
// OUTPUT:
//  entryPoint
//  exitPoint
//  @return  if the method succeeds, the return value is true; false otherwise.
bool CStconfig::getRayIntersectionPoints(const Vector3d &r0,const Vector3d &rd,Vector3d &entryPoint,Vector3d &exitPoint)
{
  if(!mOK) return false;
  volatile ctStack::saveState state(ct);
  lcms.create(this);   // will be created if needed
  if(!lcms.vertex.isOK()) return false;

  if(xyzIsInside(r0)) return false;  // return if the origin lies inside LCMS

  double dl = rminor()*0.04;     // step 2cm for W7X
  Vector3d dr = dl*rd/rd.abs();
  double Rmax = lcms.mRmax + 20*dl;
  int N = int(100*Rmax/dl);      // max iterations
  Vector3d r = r0;
  double Rmax2 = square(Rmax);
  double r2 = r.abs2();
  do {                      // cycle while not inside
    double rprev2 = r2;
    r += dr;                // move forward along ray
    r2 = r.abs2();
    if(r2>rprev2&&r2>Rmax2) return false; // return if point moves far outwards
    if(--N < 1) return false; // entry point not found
  } while(!xyzIsInside(r));   // cycle while not inside
 
  bool entryFound = getRayIntersection(r,dr,entryPoint);

  N = int(100*Rmax/dl);     // max iterations
  r2 = r.abs2();
  do {                      // cycle while inside
    double rprev2 = r2;
    r += dr;                // move forward along ray
    r2 = r.abs2();
    if(r2>rprev2&&r2>Rmax2) return false; // return if point moves far outwards
    if(--N < 1) return false; // exit point not found
  } while(xyzIsInside(r));    // cycle while inside

  bool exitFound = getRayIntersection(r,-dr,exitPoint);

  return entryFound&exitFound;
}

//*********************************************************************
// Find intersection of the Ray with the last surface
bool CStconfig::getRayIntersection(Vector3d r, Vector3d dr, Vector3d &intersectionPoint)
{
  const double accuracy = 1e-4;
// we are near LCMS, bracket the LCMS between two point rOutside and rInside
  Vector3d rOutside = r; // 'nearest' outside point
  Vector3d rInside;
  int N  = 64;               // max number of iterations
  if(xyz2s(r)<=1) {   // if current position lies inside
    do {                 // cycle while inside LCMS
      r -= dr;           // move backward along the ray
      if(--N < 1) return false;
    } while(xyz2s(r) < 1); // cycle while inside LCMS
    rOutside = r;           // 'nearest' outside point
    rInside = r+dr;        // move inside and save position
  }
  else  {                     // if current position lies slightly outside
    do {                      // cycle while outside LCMS
      r += dr;                // move forward along ray
      if(--N < 1) return false;
    } while(xyz2s(r) > 1); // cycle while outside LCMS
    rInside = r;           // inside, save position
  }
  N = 100;
// Find intersection of the Ray with the last surface by bisection
  while(--N) {
    r = (rInside+rOutside)/2;
    double ds = 1 - xyz2s(r);
    if(ds< 0.)
      rOutside = r;
    else {
      rInside = r;          // if r lies inside
      if(ds<=accuracy) break;  // Ok, done
    }
  }
  intersectionPoint=rInside;
  return true;
}

#if 1
//*********************************************************************
// @return cartesian vector with coordinates at which B equals Bres
Vector3d CStconfig::xyzAtB(double Bres,const Vector3d &r0,const Vector3d &rd, const Vector3d &entryPoint)
{
  const double epsF=1.e-4;
  const double epsA=rminor()*0.002;   // 1mm for W7X 
  const double epsR=1.e-4;
  const double dR = rminor()*0.02;   // 1cm for W7X

  setAccuracy(epsA);

  Vector3d entry(entryPoint);
  if(0==entry) {
  if(!getRayEntryPoint(r0,rd,entry))  
    return Vector3d(0.,0.,0.);  
  }

  Vector3d a, b(entry);
  Vector3d dr = dR*rd/rd.abs();

  double Fb = getBxyz(b).abs()- Bres;
  int maxItr = int(4*rminor()/dR+2);
  while(--maxItr) {      // move along the ray
    a = b + dr;
    if(!xyzIsInside(a)) return Vector3d(0.,0.,0.); 
    double Fa = getBxyz(a).abs()- Bres;
    if(Fa*Fb<0.) break;  // bracketed
    b = a;
    Fb = Fa;
  }
  
  maxItr = 100;
  while(--maxItr) {
    Vector3d xyz = (a+b)/2;
    double F = getBxyz(xyz).abs()- Bres;
    if( fabs(F)<=Bres*epsF||(b-a).abs()<=epsA+epsR*xyz.abs() ) return xyz;  
    if(F*Fb<0.) a = xyz;
    else      { b = xyz; Fb = F; }
  }
  return Vector3d(0.,0.,0.);     // exceeded maxIter
}
#endif


// find Bmin and Bmax position
void CStconfig::findBminBmaxPosition(double s1, Vector3d &magcrdBmin, Vector3d &magcrdBmax, bool useInputAsGuess)
{  
  double Bmin  =  1e20;
  double Bmax  = -1e20;
  double Bmin1 =  1e20;
  double Bmax1 = -1e20;
  double dtheta = twopi/256;   //must be less than the smallest ripple size
  double dphi   = mPeriod/64;  //must be less than the smallest ripple size
  double s = s1<0?-s1:s1;

  const double amoebaTol = 1e-15;

  Vector3d magcrdBmin1, magcrdBmax1;
  
  if(!useInputAsGuess) {
    while(1) { // find the guess for magcrdBmin , magcrdBmax
      //int nth=150, nph=50;
      int nth=int(twopi/dtheta), nph=int(mPeriod/dphi);
      int itmax = tokamakConfig?500:nth;  // theta-points (poloidal)
      int ipmax = tokamakConfig?1:nph;    // phi-points (toroidal)

      itmax = mmax(nth,itmax);
      itmax = mmax(10, itmax);

      dtheta = twopi/itmax;
      dphi   = tokamakConfig?(1./4):(0.5*mPeriod/(ipmax-1)); // stellarator symmetry is assumed
      int i;
      for(i=0; i<itmax; i++) {
        double theta = i*dtheta;
        for(int k=0; k<ipmax; k++) {
          Vector3d magCoord(s,theta,k*dphi);
          double b = B(magCoord);
          if(b<Bmin) {
            Bmin = b;
            magcrdBmin = magCoord;
          } 
          else if(b>Bmax) {
            Bmax = b;
            magcrdBmax = magCoord;
          }
        }
      }
    
      for(i=0; i<itmax; i++) {
        double theta = i*dtheta;
        for(int k=0; k<ipmax; k++) {
          Vector3d magCoord(s, theta, mPeriod - k*dphi);
          double b = B(magCoord);
          if(b<Bmin1) {
            Bmin1 = b;
            magcrdBmin1 = magCoord;
          } 
          else if(b>Bmax1) {
            Bmax1 = b;
            magcrdBmax1 = magCoord;
          }
        }
      }
      break;
    }
  }
  else {
    magcrdBmax1=magcrdBmax;
    magcrdBmin1=magcrdBmin;
    Bmax = B(magcrdBmax);
    Bmin = B(magcrdBmin);
  }

  const int ndim=2;
  double  y[ndim+1],v[(ndim+1)*ndim],*p[ndim+1];
  for (int i=0;i<ndim+1;i++) p[i] = v+i*ndim;

// find exact position of Bmax
  {
    amoebaFunc func(this,s,-1);
    p[0][0] = magcrdBmax[1];    // P0.x0  =theta0
    p[0][1] = magcrdBmax[2];    // P0.x1  =phi0
    p[1][0] = p[0][0]+dtheta;   // P1.x0
    p[1][1] = p[0][1];          // P1.x1
    p[2][0] = p[0][0];          // P2.x0
    p[2][1] = p[0][1]+dphi;     // P3.x1
    y[0] =-Bmax;                // initialize
    ////y[0] = func(p[0]); //04jun2011
    y[1] = func(p[1]);
    y[2] = func(p[2]);
    amoeba(p,y,ndim,amoebaTol,func);
    if(func.ncalls()>1000)
      std::cerr<<"CStconfig::amoeba().magcrdBmax: #func="<<func.ncalls()<<std::endl;
    magcrdBmax = Vector3d(s,p[2][0],p[2][1]);  // Bmax position
  }
if(0) {  //repeat
  Bmax  = B(magcrdBmax);
  amoebaFunc func(this,s,-1);
  p[0][0] = magcrdBmax[1];    // P0.x0  =theta0
  p[0][1] = magcrdBmax[2];    // P0.x1  =phi0
  p[1][0] = p[0][0]+dtheta;   // P1.x0
  p[1][1] = p[0][1];          // P1.x1
  p[2][0] = p[0][0];          // P2.x0
  p[2][1] = p[0][1]+dphi;     // P3.x1
  y[0] =-Bmax;                // initialize
  y[1] = func(p[1]);
  y[2] = func(p[2]);
  amoeba(p,y,ndim,amoebaTol,func);
  if(func.ncalls()>1000)
    std::cerr<<"CStconfig::amoeba().magcrdBmax: #func="<<func.ncalls()<<std::endl;
  magcrdBmax = Vector3d(s,p[2][0],p[2][1]);  // Bmax position
}
  
// find exact position of Bmax
  if(!useInputAsGuess) 
  {
    amoebaFunc func(this,s,-1);
    p[0][0] = magcrdBmax1[1];    // P0.x0  =theta0
    p[0][1] = magcrdBmax1[2];    // P0.x1  =phi0
    p[1][0] = p[0][0]-dtheta;    // P1.x0
    p[1][1] = p[0][1];           // P1.x1
    p[2][0] = p[0][0];           // P2.x0
    p[2][1] = p[0][1]-dphi;      // P3.x1
    y[0] =-Bmax1;                // initialize
    ////y[0] = func(p[0]); //04jun2011
    y[1] = func(p[1]);
    y[2] = func(p[2]);
    amoeba(p,y,ndim,amoebaTol,func);
    if(func.ncalls()>1000)
      std::cerr<<"CStconfig::amoeba().magcrdBmax: #func="<<func.ncalls()<<std::endl;
    magcrdBmax1 = Vector3d(s,p[2][0],p[2][1]);  // Bmax position
    double b  = B(magcrdBmax);
    double b1 = B(magcrdBmax1);
    if(b1>b) magcrdBmax = magcrdBmax1;
  }
if(0) {  //repeat
  Bmax1  = B(magcrdBmax1);
  amoebaFunc func(this,s,-1);
  p[0][0] = magcrdBmax[1];    // P0.x0  =theta0
  p[0][1] = magcrdBmax[2];    // P0.x1  =phi0
  p[1][0] = p[0][0]-dtheta;   // P1.x0
  p[1][1] = p[0][1];          // P1.x1
  p[2][0] = p[0][0];          // P2.x0
  p[2][1] = p[0][1]-dphi;     // P3.x1
  y[0] =-Bmax1;                // initialize
  y[1] = func(p[1]);
  y[2] = func(p[2]);
  amoeba(p,y,ndim,amoebaTol,func);
  if(func.ncalls()>1000)
    std::cerr<<"CStconfig::amoeba().magcrdBmax: #func="<<func.ncalls()<<std::endl;
  magcrdBmax1 = Vector3d(s,p[2][0],p[2][1]);  // Bmax position
  double b  = B(magcrdBmax);
  double b1 = B(magcrdBmax1);
  if(b1>b) magcrdBmax = magcrdBmax1;
}
  
  if(s1<0) return; // return if only Bmax position is needed

// find exact position of Bmin
  {
    amoebaFunc func(this,s,1);
    p[0][0] = magcrdBmin[1];    // P0.x0  =theta0
    p[0][1] = magcrdBmin[2];    // P0.x1  =phi0
    p[1][0] = p[0][0]-dtheta;   // P1.x0
    p[1][1] = p[0][1];          // P1.x1
    p[2][0] = p[0][0];          // P2.x0
    p[2][1] = p[0][1]-dphi;     // P3.x1
    y[0] = Bmin;                // initialize
    ////y[0] = func(p[0]); //04jun2011
    y[1] = func(p[1]);
    y[2] = func(p[2]);
    amoeba(p,y,ndim,amoebaTol,func);
    if(func.ncalls()>1000)
      std::cerr<<"CStconfig::amoeba().magcrdBmin: #func="<<func.ncalls()<<std::endl;
    magcrdBmin = Vector3d(s,p[2][0],p[2][1]);  // Bmin position
  }

return;
}

//****************************************************************************
//***Find minimum*************************************************************
//***Downhill Simplex Method in Multidimensions ******************************
//****************************************************************************


#define NMAX  5000
#define ALPHA 1.0
#define BETA  0.5
#define GAMMA 2.0
#define NDIM  10

#define GET_PSUM for(j=0;j<ndim;j++)  \
                   for(i=0,psum[j]=0;i<mpts;i++) psum[j] += p[i][j];

void CStconfig::amoeba(double **p,double *y,int ndim,double ftol,amoebaFunc &funk)
{
  int i,j,ilo,ihi,inhi,mpts=ndim+1;
  double psum_[NDIM], ptry_[NDIM];
  double *psum=ndim>NDIM?(new double[ndim]):psum_;
  double *ptry=ndim>NDIM?(new double[ndim]):ptry_;
  GET_PSUM;
  for(;;) {
    ilo = 0;
    ihi = y[0]>y[1]?(inhi=1,0):(inhi=0,1);
    for(i=0;i<mpts;i++) {
      if(y[i] < y[ilo]) ilo=i;
      if(y[i] > y[ihi]) {
        inhi=ihi;
        ihi=i;
      } else if(y[i] > y[inhi])
        if(i!=ihi) inhi=i;
    }
    double rtol=2*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
    if(rtol < ftol) break;
    if(funk.ncalls() >= NMAX) {
      std::cerr<<"CStconfig::amoeba(): too many iterations"<<std::endl;
      break;
    }
    double ytry=amotry(p,y,psum,ndim,ihi,-ALPHA,ptry,funk);
    if (ytry <= y[ilo])
      ytry=amotry(p,y,psum,ndim,ihi,GAMMA,ptry,funk);
    else if(ytry >= y[inhi]) {
      double ysave=y[ihi];
      ytry=amotry(p,y,psum,ndim,ihi,BETA,ptry,funk);
      if (ytry >= ysave) {
        for(i=0;i<mpts;i++) {
          if(i!=ilo) {
            for(j=0;j<ndim;j++) {
              psum[j]=(p[i][j]+p[ilo][j])/2;
              p[i][j]=psum[j];
            }
            y[i]=funk(psum);
          }
        }
        GET_PSUM;
      }
    }
  }
  if(ndim>NDIM) {
    delete[] psum;
    delete[] ptry;
  }
}

double CStconfig::amotry(double **p,double *y,double *psum,int ndim,int ihi,double fac,double *ptry,amoebaFunc &funk)
{
  double fac1=(1-fac)/ndim;
  double fac2=fac1-fac;
  for(int j=0;j<ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
  double ytry=funk(ptry);
  if(ytry < y[ihi]) {
    y[ihi]=ytry;
    for(int j=0;j<ndim;j++) {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j]=ptry[j];
    }
  }
  return ytry;
}
#undef ALPHA
#undef BETA
#undef GAMMA
#undef NMAX
#undef NDIM
#undef GET_PSUM


//****************************************************************************
//
double CStconfig::getOmega_d(Vector3d & mag1, double Bma, Vector3d & curvature1,Vector3d &dBdu1 )
{
  double  Bmod1 = getmodB();         
  double  vPar2 = 1 - Bmod1/Bma;
  double s = mag1[0];
  double iotaOriginal = iota(s);
  Vector3d cyl    = getRcyl();
  Vector3d Bvec   = getBcyl();
  Vector3d grads  = getGrads();
  Vector3d gradTh = getGradtheta();
  Vector3d gradPh = getGradphi();
  Vector3d gr     = getCovaGradB(mag1); // (dB/ds,dB/dtheta,dB/dphi)
  Vector3d gradB  = gr[0]*grads + gr[1]*gradTh + gr[2]*gradPh;
  Vector3d curv = curvature1.Cartesian2cylVector(cyl); // transform to cylindrical coordinates
  Vector3d vtmp =  (vPar2*curv + gradB/Bma/2);
  Vector3d vDrift = (Bvec^vtmp)/Bvec.abs(); 
  Vector3d gradlambda = getGradlambda();
  Vector3d gradients = gradTh+gradlambda - iotaOriginal*gradPh;  //grads/grads.abs();
      
  dBdu1   = gr;
  double omega_d = vDrift*gradients;
#if 0 //31may2012
      Vector3d clebschB =  grads^gradients;
      if(firstStep) Bsign = sign(Bvec*clebschB);
      omega_d *= Bsign;
#endif 
  return omega_d;
}

//****************************************************************************
//
Vector3d CStconfig::omega_d_BounceAveraged(const Vector3d &mag, double dphi, int iotaDenominator)
{
  std::vector<double> res = Jinv_omega_d_BounceAveraged(mag,dphi,iotaDenominator);
  Vector3d retValue(res[5],res[3],res[4]);
  return retValue;
}

//****************************************************************************
//
std::vector<double> CStconfig::Jinv_omega_d_BounceAveraged(const Vector3d &mag, double dphi, int iotaDenominator)
{
  volatile ctStack::saveState state(ct);

  int iteration=dphi<0?1:0;

  dphi = dphi==0?(0.4*degree):dphi;   // toroidal-angle increment.

  int numTurns = 200;

  double s = mag[0];

  double iotaOriginal = iota(s);
  double iotaRatnl = iotaOriginal;
  int num=iotaDenominator;
  int den=iotaDenominator;

  if(iotaDenominator) {
    double iotaR = iotaRational(s,iotaDenominator,num,den);
    iotaRatnl = iotaR!=0?iotaR:iotaRatnl;
    numTurns = den;
  }

  REAL Bma = B(mag);

  Vector3d mag0(mag), mag1(mag),r1,r0,B0,B1,gradS1;
  Vector3d dBdu0,dBdu1;  // covariant components of grad(B)==(dB/ds,dB/dtheta,dB/dphi)
  REAL omega_d0,omega_d1;
  REAL Bmod0, Bmod1,gradSmod,kG,dl;
  Vector3d curvature0(0.,0.,0.);
  Vector3d curvature1(0.,0.,0.);

  int nsteps = int(fabs(twopi/dphi));
  dphi = (twopi+1*degree)/nsteps;

  Vector3d mg = magCoordFieldLine(mag0, mag0[2]+dphi,iotaRatnl); //advance along field line
  if(B(mg)>Bma) dphi = -dphi; // d|B|/dl > 0

  
  mg = magCoordFieldLine(mag0, mag0[2]+dphi,iotaRatnl); //advance along field line
  if(B(mg)>Bma) 
    return std::vector<double>(6,0.); // d|B|/dl > 0

  advanceDataForEpsEff02(mag1, r1, B1, gradS1, Bmod1, gradSmod, kG, dl,&curvature1,0,0);  //initialize
  r0 = r1;    // save starting position
  Bma = Bmod1;
  omega_d1 = getOmega_d(mag1, Bma, curvature1, dBdu1 );

////Vector3d Bvec   = getBcyl();
////Vector3d gradTh = getGradtheta();
////Vector3d gradPh = getGradphi();
////Vector3d gr = getCovaGradB(mag1);
////Bvec *= getBdirection();
////double dBdl  = Bvec*(gr[1]*gradTh + gr[2]*gradPh);

  REAL omegaDaverage(0), tbounce(0), Jinv(0);
  REAL dJds(0),dJda(0),dJda1(0),dJda2(0);
  double dphi0 = dphi;
  int turnsCnt=0;
  int numberOfPoints = 0;

  bool lastStep=false; 
  int Bsign = 1;
  for(int i=0 ; i<numTurns; i++) {
    ++turnsCnt;
    for(int j=0 ; j<nsteps; j++) {
      mag0 = mag1;
      r0   = r1;
      B0 = B1;
      Bmod0 = Bmod1;
      omega_d0 = omega_d1;
      dBdu0 = dBdu1;
      advanceDataForEpsEff02(mag1, r1, B1, gradS1, Bmod1, gradSmod, kG, dl, &curvature1, dphi,iotaRatnl);
      REAL  vPar2 = 1 - Bmod1/Bma;

      if(Bmod1>Bma) { //do bisection to find reflection point
        double left  = mag0[2];
        double right = mag1[2];
        int N = 100;
        while(--N) {  // find turning point by bisection
          double fi = (left+right)/2;
          dphi = fi-mag0[2];          
          Vector3d mg = magCoordFieldLine(mag0, mag0[2]+dphi,iotaRatnl); //advance along field line
          REAL Bbisec = B(mg);
          if(Bbisec>Bma) right = mg[2];
          else           left = mg[2];
          if(fabs(right-left)<1e-10) {  //dphi is found
            dphi = left-mag0[2];
            mg = magCoordFieldLine(mag0, mag0[2]+dphi,iotaRatnl); //advance along field line
            break;   
          }
        }
        r1 = r0;
        mag1 = mag0;
        advanceDataForEpsEff02(mag1, r1, B1, gradS1, Bmod1, gradSmod, kG, dl, &curvature1, dphi,iotaRatnl);
        Bmod1 = Bma;
        lastStep = true;
      } else if(Bmod1==Bma) {
        lastStep = true;
      }
      ////REAL BcontraPhi = getBcontra()[2];       // == (Bvec*gradPh)
      ////REAL dl_B = fabs(dphi/BcontraPhi);       //==dphi/B^phi is the same as dl/B;
      ////dl = dl_B*B1;
      omega_d1 = getOmega_d(mag1, Bma, curvature1, dBdu1 );

      //integral_l0^l0+dl ( dx/sqrt(1-B(x)/Bma) )
      double i1 = 2*dl*sqrt(Bma)/(Bmod1-Bmod0)*(sqrt(Bma-Bmod0)-sqrt(Bma-Bmod1));
      tbounce += i1;
      omegaDaverage += i1*(omega_d0+omega_d1)/2;
      dJds  += i1*(dBdu0[0]+dBdu1[0])/2;
      dJda1 += i1*(dBdu0[1]+dBdu1[1])/2;
      dJda2 += i1*(dBdu0[2]+dBdu1[2])/2;

      //integral_l0^l0+dl ( dx*sqrt(1-B(x)/Bma) )
      Jinv += 2./3.*dl/((Bmod1-Bmod0)*sqrt(Bma))*(pow(Bma-Bmod0,1.5)- pow(Bma-Bmod1,1.5));
      if(lastStep) goto ext;

      numberOfPoints++;
    }
  }
ext:
  dJds = -dJds/(2*Bma);
  dJda = -(dJda1*iota(s)-dJda2)/(2*Bma);
  omegaDaverage = tbounce==0?0:(omegaDaverage/tbounce);
  double velocity = 1; //1.876e7;  // 1keV electron velocity 
  //  2*tbounce/velocity == \oint\frac{dl}{v\sqrt{1-B/B_m}}
  std::vector<double> retValue; 
  retValue.push_back(Jinv);
  retValue.push_back(dJds);
  retValue.push_back(dJda);
  retValue.push_back(2*tbounce);
  retValue.push_back(Bma);
  retValue.push_back(omegaDaverage);
  if(iteration==1&&numberOfPoints<2) {
    retValue=std::vector<double>(6,0.);
  }
  else if(iteration==0&&numberOfPoints<4) {
    retValue=Jinv_omega_d_BounceAveraged(mag, -fabs(dphi0/2), iotaDenominator);
  }
  return retValue;
}

//*********Jinvariant*************************************************************
//*********Jinvariant*************************************************************
//*********Jinvariant*************************************************************

#define JTEST 0

#if JTEST
  double Btest(const Vector3d &mag) {
    // Btest(phi) = (bmix-bmin)*|phi-a|/a + bmin 
    const double b=2,c=1,a=3.14159265358979323846/5;
    return (b-c)*fabs(mag[2]-a)/a + c;
  }

  Vector3d getCovaGradBtest(const Vector3d &mag) {
    // Btest(phi) = (bmix-bmin)*|phi-a|/a + bmin 
    const double b=2,c=1,a=3.14159265358979323846/5;
    Vector3d retvalue(0);
    retvalue[2]=(b-c)/a*((mag[2]-a)<0?-1:1);
    return retvalue;
  }

  double integralVpardphi(double d) {
    // integral dphi*sqrt(1-Btest(phi)/d) 
    const double b=2,c=1,a=3.14159265358979323846/5;
    return 4*a*d*pow((d-c)/d,1.5)/(3*(b-c));
  }
  #define B(mag) Btest(mag)
  #define getCovaGradB(mag) getCovaGradBtest(mag)
 
#endif

std::vector<double> CStconfig::Jinvariant(const Vector3d &mag, double dphi, int iotaDenominator)
{
  volatile ctStack::saveState state(ct);

  int iteration=dphi<0?1:0;

  dphi = dphi==0?(0.4*degree):dphi;   // toroidal-angle increment.

  int numTurns = 200;

  double s = mag[0];

  double iotaOriginal = iota(s);
  double iotaRatnl = iotaOriginal;
  int num=iotaDenominator;
  int den=iotaDenominator;

  if(iotaDenominator) {
    double iotaR = iotaRational(s,iotaDenominator,num,den);
    iotaRatnl = iotaR!=0?iotaR:iotaRatnl;
    numTurns = den;
  }


  REAL Bma = B(mag);  // Bma is turning point

  Vector3d mag0(mag), mag1(mag),r1,r0;
  Vector3d dBdu0,dBdu1;  // covariant components of grad(B)==(dB/ds,dB/dtheta,dB/dphi)
  REAL Bmod0,Bmod1,dl;

  int nsteps = int(fabs(twopi/dphi));
  dphi = (twopi+1*degree)/nsteps;

  Vector3d mg = magCoordFieldLine(mag0, mag0[2]+dphi,iotaRatnl); //advance along field line
  if(B(mg)>Bma) dphi = -dphi; // d|B|/dl > 0

  
  mg = magCoordFieldLine(mag0, mag0[2]+dphi,iotaRatnl); //advance along field line
  if(B(mg)>Bma) 
    return std::vector<double>(5,0.); // d|B|/dl > 0

  advanceDataForJinvariant(mag1, r1, Bmod1,dBdu1, dl,0,0);  //initialize
  r0 = r1;    // save starting position
  Bma = Bmod1;

  REAL tbounce(0), Jinv(0);
  REAL dJds(0),dJda(0),dJda1(0),dJda2(0); 
  double dphi0 = dphi;
  int turnsCnt=0;
  int numberOfPoints = 0;

  bool lastStep=false; 
  for(int i=0 ; i<numTurns; i++) {
    ++turnsCnt;
    for(int j=0 ; j<nsteps; j++) {
      mag0 = mag1;
      r0   = r1;
      Bmod0 = Bmod1;
      dBdu0 = dBdu1;
      advanceDataForJinvariant(mag1, r1, Bmod1,dBdu1,dl, dphi,iotaRatnl);
      //double vPar2 = 1 - Bmod1/Bma;

      if(Bmod1>Bma) { //do bisection to find reflection point
        double left  = mag0[2];
        double right = mag1[2];
        int N = 100;
        while(--N) {  // find turning point by bisection
          double fi = (left+right)/2;
          dphi = fi-mag0[2];          
          Vector3d mg = magCoordFieldLine(mag0, mag0[2]+dphi,iotaRatnl); //advance along field line
          REAL Bbisec = B(mg);
          if(Bbisec>Bma) right = mg[2];
          else           left = mg[2];
          if(fabs(right-left)<1e-10) {  //dphi is found
            dphi = left-mag0[2];
            mg = magCoordFieldLine(mag0, mag0[2]+dphi,iotaRatnl); //advance along field line
            break;   
          }
        }
        r1 = r0;
        mag1 = mag0;
        advanceDataForJinvariant(mag1, r1, Bmod1,dBdu1,dl, dphi,iotaRatnl);
        //vPar2 = 1 - Bmod1/Bma;
        //if(vPar2<0) std::cout<< dl<<" " <<Bmod1<<" "<<vPar2<<" "<<std::endl;
        Bmod1 = Bma;
        lastStep = true;
      } else if(Bmod1==Bma) {
        lastStep = true;
      }

      //integral_l0^l0+dl ( dx*sqrt(1-B(x)/Bma) )
      Jinv += 2./3.*dl/((Bmod1-Bmod0)*sqrt(Bma))*(pow(Bma-Bmod0,1.5)- pow(Bma-Bmod1,1.5));
      //integral_l0^l0+dl ( dx/sqrt(1-B(x)/Bma) )
      double i1 = 2*dl*sqrt(Bma)/(Bmod1-Bmod0)*(sqrt(Bma-Bmod0)-sqrt(Bma-Bmod1));
      tbounce += i1;
      dJds  += i1*(dBdu0[0]+dBdu1[0])/2;
      dJda1 += i1*(dBdu0[1]+dBdu1[1])/2;
      dJda2 += i1*(dBdu0[2]+dBdu1[2])/2;
      if(lastStep) goto ext;
 // std::cout<<"mag1="<<mag1<<" dBdu="<<dBdu1<<" Bt="<<Bmod1<<" dl/sqrt="<<dl_vPar<<std::endl;
      numberOfPoints++;
    }
  }
ext:
  dJds = -dJds/(2*Bma);
  dJda = -(dJda1*iota(s)-dJda2)/(2*Bma);
  double velocity = 1; //1.876e7;  // 1keV electron velocity 
  //  2*tbounce/velocity == \oint\frac{dl}{v\sqrt{1-B/B_m}}
#if JTEST
  Jinv -= integralVpardphi(Bma); // Jinv-Jinvtest must be zero
#endif
  std::vector<double> retValue; 
  retValue.push_back(Jinv);
  retValue.push_back(dJds);
  retValue.push_back(dJda);
  retValue.push_back(2*tbounce);
  retValue.push_back(Bma);
  if(iteration==1&&numberOfPoints<2) {
    retValue=std::vector<double>(5,0.);
  }
  else if(iteration==0&&numberOfPoints<4) {
    retValue=Jinvariant(mag, -fabs(dphi0/2), iotaDenominator);
  }
  return retValue;

}

//****************************************************************************
void CStconfig::advanceDataForJinvariant(Vector3d &mag, Vector3d &r1,REAL &Bmod,Vector3d &covaGradB, REAL &dl, 
                                       const double &dphi, const double &iotaRatnl)
{
//see Jinvariant()  volatile ctStack::saveState state(ct);
 // mag is the vector with the components (s, theta, phi); where theta is the poloidal angle
 // advance along field line from (s, theta, phi)  to (s, thetaNew, phi+dphi)
 // for case of Boozer coordinates: (s, theta, phi)  -->  (s, theta+iota*dphi, phi+dphi)
  mag = magCoordFieldLine(mag, mag[2]+dphi,iotaRatnl); // advance along field line
#if 0
  NJacobian(mag);             // calculate all quantities at point mag
  Vector3d r2 = getRxyz();    // get position in cartesian coordinates
  Bmod = REAL(getmodB());     // |B|
  covaGradB = getCovaGradB(mag);
#else
  Vector3d r2 = mag2xyz(mag);    // get position in cartesian coordinates
  Bmod = B(mag);                 // |B|
  covaGradB = getCovaGradB(mag); // (dB/ds,dB/dtheta,dB/dphi)
#endif
  if(dphi!=0) {
    Vector3d dr = r2 - r1;
    dl = REAL(dr.abs());
  }
  else { 
    dl = 0;
  }
#if JTEST
  dl=fabs(dphi);
#endif
  r1 = r2;
}

#undef JTEST


//****************************************************************************
// The method returns mag. coord on a field line at which B(mag) == Bt
//  Bmin is seached within one machine period: 0<=phi<2*pi/Np
// @param mag0 is the initial point on the field line
Vector3d CStconfig::getMagCoordAtB(const Vector3d & mag0, double Bt) {
  // alpha = theta/iota - phi
  Vector3d mag1=mag0, mag2, magmin;

//find position of the minimum value of the magnetic field |B| on axis
  double Bmi=1e20;
  double Bma=-1e20;
  mag1 = magCoordFieldLine(mag0, 0,0); //advance along field line to phi=0
  int ip = 500;
  for(int i=0; i<ip; i++) {
    double dphi = mPeriod/ip;
    mag2 = magCoordFieldLine(mag1, i*dphi,0); //advance along field line
    double b = B(mag2);
    if(b<Bmi) {
      Bmi=b;
      magmin = mag2; // magnetic coordinate at Bmin for considered field line
    }
    Bma = mmax(Bma,b);
  }

  if(Bt>=Bma||Bt<=Bmi) return Vector3d(-1,Bmi,Bma);

  double dphi = 0.5*degree;
  double B1,B2;
  mag1 = magmin;
  B1 = B(mag1);
  while(1) {
    mag2 = magCoordFieldLine(mag1, mag1[2]+dphi,0); //advance along field line
    B2 = B(mag2);
    if(B2>Bt) break;
    B1 = B2;
    mag1 = mag2;
  }
 //do bisection to find mag. coord at which B(mag) == Bt
  double left  = mag1[2];
  double right = mag2[2];
  double Bbisec;
  Vector3d mg;
  int N = 100;
  while(--N) { 
    double phi = (left+right)/2;
    mg = magCoordFieldLine(mag1,phi,0); //advance along field line
    Bbisec = B(mg);
    if(Bbisec>Bt) right = mg[2];
    else           left = mg[2];
    if(fabs(right-left)<1e-10) {
      break;   
    }
  }
  return mg;
}

//****************************************************************************
// The method returns mag. coord on a field line at which B(mag) == Bt
std::vector<Vector3d> CStconfig::getAllMagCoordAtB(const Vector3d & mag0, double Bt) 
{
  using namespace std;
  Vector3d mag1=mag0, mag2;
  vector< vector<Vector3d> >maglist;
  vector<Vector3d> magvalues;  //Coordinates for Bt left of well
  magvalues.push_back(mag1);
  magvalues.push_back(mag1);

  mag1 = magCoordFieldLine(mag0,-pi,0);  //advance along field line to phi=-pi
  int ip = 400;
  for(int i=0; i<ip; i++) {
    double dphi = twopi/ip;
    mag2 = magCoordFieldLine(mag1, i*dphi-pi,0); //advance along field line
    if(B(mag1)>Bt && B(mag2)<Bt) {           //For Bt at the left side of the well
      magvalues[0]=mag1;
      magvalues[1]=mag2;
      maglist.push_back(magvalues);
    }
    mag1=mag2;
  }

  vector<Vector3d> Btcoordinates;
  for(unsigned int i=0;i<maglist.size();i++){
    double left  = maglist[i][0][2];
    double right = maglist[i][1][2];
    double Bbisec;
    Vector3d mg;
    int N = 64;
    while(--N) {
      double phi = (left+right)/2;
      mg = magCoordFieldLine(maglist[i][0],phi,0);
      Bbisec = B(mg);
      if(Bbisec<Bt) right = mg[2];
      else           left = mg[2];
      if(fabs(right-left)<1e-8) {
        break;
      }
    }
    Btcoordinates.push_back(mg);
  }
  return Btcoordinates;
}

//****************************************************************************
// The method returns mag. coord on a field line at which B(mag) == Bt
Vector3d CStconfig::getNextMagCoordAtB(const Vector3d & mag0, double Bt) 
{
  using namespace std;
  Vector3d mag1=mag0, mag2;
  vector<Vector3d> magvalues;  
  magvalues.push_back(mag1);
  magvalues.push_back(mag1);

  int ip = 400;
  for(int i=0; i<ip; i++) {
    double dphi = twopi/ip;
    mag2 = magCoordFieldLine(mag1, mag1[2]+dphi,0); //advance along field line
    if(B(mag1)<Bt && B(mag2)>=Bt){   // ->
      magvalues[0]=mag1;
      magvalues[1]=mag2;
      break;
    }
    mag1=mag2;
  }
  
  double left =magvalues[0][2];
  double right=magvalues[1][2];
  double Bbisec;
  Vector3d mg;
  int N = 64;
  while(--N) {
    double phi = (left+right)/2;
    mg = magCoordFieldLine(magvalues[0],phi,0); //advance along field line
    Bbisec = B(mg);
    if(Bbisec>Bt) right = mg[2];
    else           left = mg[2];
    if(fabs(right-left)<1e-8) {
      break;
    }
  }
  Vector3d Btcoordinate=mg;
  return Btcoordinate;
}


#if 0 
// TODO: desabled; need more testing!!
//****************************************************************************
// The method returns mag. coord on a field line at which B(mag) == Bt
// @param mag0 is the initial point on the field line
Vector3d CStconfig::getMagCoordAtB(const Vector3d & mag0, double Bt) {
  // alpha = theta/iota - phi
  double dphi = 0.5*degree;
  Vector3d mag1=mag0, mag2;

  double B1 = B(mag1);
  mag2 = magCoordFieldLine(mag1, mag1[2]+dphi,0); //advance along field line
  double B2 = B(mag2);
  dphi = sign(Bt-B1)*sign(B2-B1)*dphi;

  if(B1>Bt) 
      while(1) {
        mag2 = magCoordFieldLine(mag1, mag1[2]+dphi,0); //advance along field line
        B2 = B(mag2);
        if(B2<Bt) break;
        B1 = B2;
        mag1 = mag2;
      }
  else
      while(1) {
        mag2 = magCoordFieldLine(mag1, mag1[2]+dphi,0); //advance along field line
        B2 = B(mag2);
        if(B2>Bt) break;
        B1 = B2;
        mag1 = mag2;
      }
 //do bisection to find mag. coord at which B(mag) == Bt
  double left  = dphi<0?mag2[2]:mag1[2];
  double right = dphi<0?mag1[2]:mag2[2];
  double Bbisec;
  Vector3d mg;
  int N = 27;
  while(--N) {  // find turning point by bisection
    double phi = (left+right)/2;
    mg = magCoordFieldLine(mag1, phi,0); //advance along field line
    Bbisec = B(mg);
    if((Bbisec-Bt)*sign(B2-B1)>0) right = mg[2];
    else           left = mg[2];
    if(fabs(right-left)<1e-8) {
      break;   
    }
  }
  
  std::cout<<Bbisec<<" "<<Bt<<" "<<N<<std::endl;
  return mg;
}
#endif 

///****************************************************************************
/// The method returns the coefficients needed for 
/// the poloidal flux equation in Astra transport code
void CStconfig::getCoeffForAstraCode(double sqrts, double &reff, double &gradr2Avr, double &J, double &G2, double &hVprime, double &B0, double &R0, double &h)
{
  double s=sqrts*sqrts;
  double a = r(1.);
  reff = r(s);
  h = getJacobianSign();
  J = Ip(s)/Ip(1.);
  G2 = s==0?0:(h*S11(s)/Flux(1.)*a/(2*sqrts));
  hVprime = h*Vprime(s)*2*reff/(a*a);
  B0 = h*Flux(1.)/(pi*a*a);
  R0 = mu0*Ip(1.)/(2*pi*B0);
  if(s==0) s=1e-8;
  double dr_ds  = r(s)/(2*s); // dr/ds 
  gradr2Avr = fabs(Grads2Avrg(s)*dr_ds*dr_ds);

}


};   //namespace MConf

