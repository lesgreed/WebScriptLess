#include "../include/CStconfig.h"

namespace MConf {

//*************************************************************
// Evaluate the coefficients for a cubic interpolating spline
// !!Attention
//   remove lines bracketed by comments //CStConfig
//   if you want to use this function in other applications
//http://www.mech.uq.edu.au/staff/jacobs/nm_lib/cmathsrc/spline.c
//
int CStconfig::spline (int n, const double x[], const double y[],
            double b[], double c[], double d[],
            int end1, int end2,
            double slope1, double slope2)
/* Evaluate the coefficients b[i], c[i], d[i], i = 0, 1, .. n-1 for
   a cubic interpolating spline

   S(xx) = Y[i] + b[i] * w + c[i] * w**2 + d[i] * w**3
   where w = xx - x[i] and  x[i] <= xx <= x[i+1]

   The n supplied data points are x[i], y[i], i = 0 ... n-1.

Input :
   n       : The number of data points or knots (n >= 2)
   end1,
   end2    : = 1 to specify the slopes at the end points
             = 0 to obtain the default conditions
   slope1,
   slope2  : the slopes at the end points x[0] and x[n-1]
             respectively
   x[]     : the abscissas of the knots in increasing order
   y[]     : the ordinates of the knots

 Output :
   b, c, d : arrays of spline coefficients as defined above
             (See note 2 for a definition.)
 Return:
            = 0 normal return
            = 1 less than two data points; cannot interpolate
            = 2 x[] are not in ascending order

   This C code written by ...  Peter & Nigel,
   ----------------------      Design Software,
                               42 Gubberley St,
                               Kenmore, 4069,
                               Australia.

   Version ... 1.1, 30 September 1987
   -------     2.0, 6 April 1989    (start with zero subscript)
                                     remove ndim from parameter list
               2.1, 28 April 1989   (check on x[])
               2.2, 10 Oct   1989   change number order of matrix

   Notes ...
   -----
   (1) The accompanying function seval() may be used to evaluate the
       spline while deriv will provide the first derivative.
   (2) Using p to denote differentiation
       y[i] = S(X[i])
       b[i] = Sp(X[i])
       c[i] = Spp(X[i])/2
       d[i] = Sppp(X[i])/6  ( Derivative from the right )
   (3) Since the zero elements of the arrays ARE NOW used here,
       all arrays to be passed from the main program should be
       dimensioned at least [n].  These routines will use elements
       [0 .. n-1].
   (4) Adapted from the text
       Forsythe, G.E., Malcolm, M.A. and Moler, C.B. (1977)
       "Computer Methods for Mathematical Computations"
       Prentice Hall
   (5) Note that although there are only n-1 polynomial segments,
       n elements are requird in b, c, d.  The elements b[n-1],
       c[n-1] and d[n-1] are set to continue the last segment
       past x[n-1].
*----------------------------------------------------------------*/
{
  int i, nm1=n-1;
  if(n<2) return 1;     // no possible interpolation
//CStConfig  for(i=1; i<n; ++i) if(x[i] <= x[i-1])  return 2;  //x[] are not in ascending order
  if(n==2) { /* linear segment only  */
    b[1]=b[0]=(y[1]-y[0])/(x[1]-x[0]);
    c[1]=c[0]=0;
    d[1]=d[0]=0;
    return 0;
  }

//CStConfig-beg
  end1=0;
  end2=0;
   slope1 = (y[1]-y[0])/(x[1]-x[0]);
   slope2 = (y[n-1]-y[n-2])/(x[n-1]-x[n-2]);
//CStConfig-end


//    Set up the symmetric tri-diagonal system
//  b = diagonal d = offdiagonal c = right-hand-side
  d[0] = x[1] - x[0];
  c[1] = (y[1] - y[0]) / d[0];
  for(i = 1; i < nm1; ++i) {
    d[i]   = x[i+1] - x[i];
    b[i]   = 2*(d[i-1]+d[i]);
    c[i+1] = (y[i+1]-y[i])/d[i];
    c[i]   = c[i+1]-c[i];
  }
  /* ---- Default End conditions
      Third derivatives at x[0] and x[n-1] obtained
      from divided differences  */
  b[0]   = -d[0];
  b[nm1] = -d[n-2];
  c[0]   = 0;
  c[nm1] = 0;
  if(n != 3) {
    c[0]   = c[2]/(x[3]-x[1]) - c[1]/(x[2]-x[0]);
    c[nm1] = c[n-2]/(x[nm1]-x[n-3]) - c[n-3]/(x[n-2]-x[n-4]);
    c[0]   = c[0]*d[0]*d[0]/(x[3]-x[0]);
    c[nm1] = -c[nm1]*d[n-2]*d[n-2] / (x[nm1]-x[n-4]);
  }
  /* Alternative end conditions -- known slopes */
  if(end1) {
    b[0] = 2*(x[1]-x[0]);
    c[0] = (y[1]-y[0])/(x[1]-x[0]) - slope1;
  }
  if(end2) {
    b[nm1] = 2*(x[nm1]-x[n-2]);
    c[nm1] = slope2 - (y[nm1]-y[n-2])/(x[nm1]-x[n-2]);
  }
  /* Forward elimination */
  for(i = 1; i < n; ++i) {
    double t = d[i-1]/b[i-1];
    b[i] = b[i]-t*d[i-1];
    c[i] = c[i]-t*c[i-1];
  }
  /* Back substitution */
  c[nm1] = c[nm1]/b[nm1];
  for(i=n-2; i>=0; i--)
    c[i] = (c[i]-d[i]*c[i+1])/b[i];

  /* c[i] is now the sigma[i] of the text */
  /* Compute the polynomial coefficients */
  b[nm1] = (y[nm1]-y[n-2])/d[n-2]+d[n-2]*(c[n-2]+2*c[nm1]);
  for(i=0; i<nm1; ++i) {
    b[i] = (y[i+1]-y[i])/d[i] - d[i]*(c[i+1]+2*c[i]);
    d[i] = (c[i+1]-c[i])/d[i];
    c[i] = 3*c[i];
  }
  c[nm1] = 3*c[nm1];
  d[nm1] = d[n-2];
  return 0;
#if 0 
//CStConfig-beg
    i=n-2;
    b[i] = slope2;
    d[i] = c[i] = 0;
    i=n-1;
    b[i] = slope2;
    d[i] = c[i] = 0;
//CStConfig-end
  return 0;
#endif
}

//***************************************************************************
//
double CStconfig::seval (int n, double u,
              const double x[], const double y[],
              const double b[], const double c[], const double d[],
              double &yp)

/*
  Evaluate the cubic spline function
  and the derivative of the cubic spline function

  S(u)  = y[i]+b[i]*w+c[i]*w^2+d[i]*w^3
  S'(u) = b[i]+2*c[i]*w + 3*d[i]*w**2
  where w = u - x[i]
  and   x[i] <= u <= x[i+1]

  If u < x[0]   then i = 0 is used.
  If u > x[n-1] then i = n-1 is used.

  Input :
  n       : The number of data points or knots (n >= 2)
  u       : the abscissa at which the spline is to be evaluated
  x[]     : the abscissas of the knots in increasing order
  y[]     : the ordinates of the knots
  b, c, d : arrays of spline coefficients computed by spline().

  Return:
  seval   : the value of the spline function at u
  yp      : the derivative of the cubic spline function
*/
{
  int i=0;        //Left bracket
  int j=n;        //Right bracket
  if(u<=x[0]) j=1;
  else if(u>=x[j-1]) i=j-1;
  else
    while(i+1 < j) {            //bisection
      int k = (i+j)/2;
      if(u < x[k]) j = k;
      else i = k;
    }

  double dx = u - x[i];
  double s  = y[i] + dx*(b[i]+dx*(c[i]+dx*d[i])); // evaluate the spline
  yp = b[i] + dx*(2*c[i]+dx*3*d[i]);              // evaluate the derivative
  return (s);
}


//***************************************************************************
//  Given arrays x[0..n-1] and y[0..n-1] containing a tabulated function, i.e., yi = f(xi), with
//  x1 < x2 < .. . < xN, and given values yp1 and ypn for the first derivative of the interpolating
//  function at points 1 and n, respectively, this routine returns an array y2[0..n-1] that contains
//  the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
//  ypn are equal to 1e30 or larger, the routine is signaled to set the corresponding boundary
//  condition for a natural spline, with zero second derivative on that boundary.
//  Numerical Recipes in C

int CStconfig::splineNRC(int n, const double *x, const double *y, double *y2, double yp1, double ypn)
{
  double qn(0),un(0);
  if(n<3) return 1;     // no possible interpolation
  double *u=new double[n];

  if (yp1 > 0.9e30) {
    y2[0]=u[0]=0;  // natural boundary condition
  }
  else {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  if (ypn > 0.9e30) {
    qn=un=0;       // natural boundary condition
  }
  else {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }

  for (int i=1;i<=n-2;i++) {   // decomposition loop
    double sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    double p=sig*y2[i-1]+2;
    y2[i]=(sig-1)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1);   // backsubstitution loop
  for (int k=n-2;k>=0;k--) y2[k]=y2[k]*y2[k+1]+u[k];

  delete[] u;
  return 0;
}

//***************************************************************************
//  Given the arrays xa[0..n-1] and ya[0..n-1], which tabulate a function (with the xa’s in order),
//  and given the array y2a[0..n-1], which is the output from spline above, and given a value of
//  x, this routine returns a cubic-spline interpolated value y and first derivative yp.
//  Numerical Recipes in C
double CStconfig::sevalNRC(int n,const double *xa,const double *ya,const double *y2a, double x, double &yp)
{
  int L=0,R=n-1;

  if(x<=xa[0]) R=1;
  else if(x>=xa[R]) L=R-1;
  else
    while (R-L > 1) {  // bisection
      int k=(R+L) >> 1;
      if (xa[k] > x) R=k;
      else L=k;
    }

  double h=xa[R]-xa[L];
  double b=(x-xa[L])/h;
  double a=1-b;       // a=(xa[R]-x)/h;
  yp=(ya[R]-ya[L])/h-( (3*a*a-1)*y2a[L]-(3*b*b-1)*y2a[R] )*h/6;
  double y = ya[L]+b*(ya[R]-ya[L])+( (a*a*a-a)*y2a[L]+(b*b*b-b)*y2a[R] )*h*h/6;
  return y;
}

};   //namespace MConf
