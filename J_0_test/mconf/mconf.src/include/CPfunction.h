#ifndef MC_CPfunction_
#define MC_CPfunction_
#include <math.h>
#include "ngarray.h"

namespace MConf {

//*********************************************************************
/// This class provides the function for creating a plasma profile.
/// f(x) = g-h+(1-g+h)(1-x<sup>p</sup>)<sup>q</sup> + 
///          h(1-e<sup>-x<sup>2</sup>/w<sup>2</sup></sup>),
///    where 0<=x<=1, \par
///  properties: g=PF(1)/PF(0), PF(0) = 1, PF(x)>0,
/// \note don't use it directly, use CProfile instead
class CPfunction {
  double c;
  public:
  CPfunction() {c=0;}
  virtual ~CPfunction() {;}
  /// PF(x) = g-h+(1-g+h)(1-x<sup>p</sup>)<sup>q</sup> + 
  ///          h(1-e<sup>-x<sup>2</sup>/w<sup>2</sup></sup>),
  ///    where 0<=x<=1.
  /// @param x is the abscissa, 0<=x<=1, can be r<sub>eff</sub>/a.
  /// @param parm is the array[5] of coefficients: g,p,q,h,w.
  /// @return function's value.
  double PF(double x,const double *parm) const  // x is the r/a
  {
    if(x<0) x=-x;
    if(x>1) x=1;
    double a=parm[0], p=parm[1], q=parm[2], h=parm[3],hole=0;
    if(fabs(h)>=1e-8) { // if hollow
      a -= h;
      double w = fabs(parm[4])<1e-3?1e-3:parm[4];
      double e = fabs((x-c)/w);
      hole = (e>6)?h:(h*(1-::exp(-e*e)));
    }
    double f = a + hole;
    if(q==0) f=1+hole;
    else if(x<1) f += (1-a)*::pow(1-pow(x,p),q);
    return (f>0)?f:1e-15; // n,T or Zeff profiles must be greater then zero
  }
//28.09.15: need testing 
  void setHoleCenter(double c) {
    this->c=c;   // center of hole
  }
};

/// This class provides the plasma profile which can be set by set of points 
/// or parameterized by CPfunction().
/// F(x) =n0[g-h+(1-g+h)(1-x<sup>p</sup>)<sup>q</sup> + 
///          h(1-e<sup>-x<sup>2</sup>/w<sup>2</sup></sup>)],
///    where 0<=x<=1, \par
///  properties: g=F(1)/F(0), F(0) = 1, F(x)>0, \par
///  \anchor usageOfCProfile
///  usage:
///
/// \code
///  CProfile ne;
///  // n0 is the value in the center
///  // parm is the array[5] with coefficients: g,p,q,h,w.
///  ne.set(n0,parm);
///  double n = ne(x);
/// \endcode
/// \image html CProfile.png
class  CProfile {
  CPfunction Pf;
  double a0;      // value in the center
  double parm[5];
  double ymin;
  double ymax;
  double xma;
  pBase::ngArray<double> x_;
  pBase::ngArray<double> y_;
  bool useAnalitic;
  bool doMinMax;
public:
  CProfile() {a0=0;doMinMax=useAnalitic=false;xma=1;}
  /// Constructor sets parameters for analitic representation of 
  /// the profile function.
  /// @param v0 is the value in the center
  /// @param param is the array[5] of coefficients: g,p,q,h,w.
  /// See detailed description of CProfile::set()
  CProfile(double v0, double *param) { set(v0,param);}
  /// Constructor sets parameters for analitic representation of 
  /// the profile function.
  /// @param v0 is the value in the center
  /// @param g,p,q,h,w are the coefficients.
  /// See detailed description of CProfile::set()
  CProfile(double v0, double g, double p,double q, double h=0, double w=0, double c=0) {
    parm[0] = g; 
    parm[1] = p; 
    parm[2] = q; 
    parm[3] = h; 
    parm[4] = w; 
    set(v0,parm);
    Pf.setHoleCenter(c); 
  }

  virtual ~CProfile() {;}
  /// The method sets parameters for analitic representation of 
  /// the profile function.
  /// F(x) = c[g-h+(1-g+h)(1-x<sup>p</sup>)<sup>q</sup> + 
  ///          h(1-e<sup>-x<sup>2</sup>/w<sup>2</sup></sup>)],
  ///    where 0<=x<=1, \par
  /// @param v0 is the value in the center
  /// @param param is the array[5] of coefficients: g,p,q,h,w.
  /// @param c is the center of the hole.
  void set(double v0, double *param,double c=0) {
    clear();
    doMinMax=useAnalitic=true;
    xma=1;
    this->a0=v0;
    for(int i=0; i<5; i++) this->parm[i] = param[i];
    Pf.setHoleCenter(c); 
  }
  /// The method sets xmax for analitic representation of 
  /// the profile function.
  void setXmax(double xmax) {
    xma=xmax;
  }
  /// The method sets points of the profile function from arrays.
  /// For evaluating the function in intermediate points
  /// linear interpolation is used, see interpolate()
  /// \param x is the abscissa points, array of size n,
  ///   0\<=x\<=1 (dimentionless flux label)
  /// \param p is the ordinate set, array of size n.
  /// \param n is the size of arrays.
  void set(const double *x,const double *p,int n) {
    set(x,p,1,n);
  }
  /// The method sets points of the profile function from arrays.
  /// For evaluating the function in intermediate points
  /// linear interpolation is used, see interpolate()
  /// \param r is the abscissa points, array of size n, 0<=r<=rmax
  /// \param p is the ordinate set, array of size n.
  /// \param rmax is the max value of abscissa, can be larger than last r,
  ///    internally \c r is normalized by \c rmax
  /// \param n is the size of arrays.
  void set(const double *r,const double *p,double rmax,int n) {
    a0=0;
    useAnalitic=doMinMax=false;
    x_.init(n);
    y_.init(n);
    ymin=ymax=p[0];
    xma=1;
    for(int i=0; i<n; i++) {
      x_[i] = r[i]/rmax;
      y_[i] = p[i];
      ymin = ymin>p[i]?p[i]:ymin;
      ymax = ymax<p[i]?p[i]:ymax;
      xma = xma<x_[i]?x_[i]:xma;
    }
  }
  /// The method clears the profile, it becomes empty.
  void clear() {
    a0=0;
    x_.clear();
    y_.clear();
    useAnalitic=doMinMax=false;
    Pf.setHoleCenter(0); 
  }
 /// The method returns true if profile was not set.
  bool empty() const  {
    return useAnalitic?(a0==0):(x_.empty());
  }
  /// The operator returns function value at point x
  /// \param x is the abscissa
  /// see also \ref usageOfCProfile "usage of CProfile"
  double operator() (const double x) const {
    return useAnalitic?value(x):interpolate(x); // 0<=x<=1, x = r/a
  }
 /// The method returns minimal value of function.
  double mmin() {
    return empty()?0:(minmax(200),ymin);
  }
 /// The method returns maximal value of function.
 double mmax() {
    return empty()?0:(minmax(200),ymax);
  }
 double xmax() {
    return xma;
  }
  /// The method scales the function.
  /// \param factor is the scale-factor
 void scale(double factor) {
    a0 *=factor;
    ymin *=factor;
    ymax *=factor;
    if(useAnalitic) return;
    if(empty()) return;
    for(int i=0; i<(int)y_.size(); i++) y_[i] *=factor;
  }
private:
  // find y(x) -- linear interpolation
  double interpolate(double u) const  {
    if(x_.empty()) return 0;
    int L = 0;                // Left bracket
    int R = int(x_.size()-1); // Right bracket
    // perform a binary search of the segment in which u lies
    if(u<=x_[0]) R=1;
    else if(u>=x_[R]) return y_[R];  //L=R-1;
    else
      while(L+1 < R) {       // bisection
        int k = (L + R)/2;
        if(u < x_[k]) R = k;
        else L = k;
      }
    return (y_[L] + (y_[R]-y_[L])*(u-x_[L])/(x_[R]-x_[L]));
  }

  double value(double x) const  {
    x /= xma;
    return a0*Pf.PF(x,parm);
  }

  void minmax(int N) {
    if(!doMinMax) return;
    if(!useAnalitic) return;
    ymin=ymax=value(0);
    for(int i=1; i<=N; i++) {
      double x = double(i)/N;
      double p = value(x);
      ymin = ymin>p?p:ymin;
      ymax = ymax<p?p:ymax;
    }
    doMinMax=false;
  }

};

};   //namespace MConf

#endif
