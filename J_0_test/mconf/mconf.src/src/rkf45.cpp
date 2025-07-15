#include "../include/rkf45.h"

namespace MConf {

//     test to see if rkf45 is being severely impacted by too many
//     output points
static const int KOP_MAX=1000000000;
//     the expense is controlled by restricting the number
//     of function evaluations to be approximately MAX_NFE.
static const int MAX_NFE=1000000000;

void rkf45::rkf_45(RightHandSide &rhs,int neqn,double* y,
            double& t,double& tout,double& relerr,double abserr,int& iflag,
            double* work,int* iwork)
{
/*	

     fehlberg fourth-fifth order runge-kutta method

     written by h.a.watts and l.f.shampine
                   sandia laboratories
                  albuquerque,new mexico

    rkf45 is primarily designed to solve non-stiff and mildly stiff
    differential equations when derivative evaluations are inexpensive.
    rkf45 should generally not be used when the user is demanding
    high accuracy.

 abstract

    subroutine  rkf45  integrates a system of neqn first order
    ordinary differential equations of the form
             dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
              where the y(i) are given at t .
    typically the subroutine is used to integrate from t to tout but it
    can be used as a one-step integrator to advance the solution a
    single step in the direction of tout.  on return the parameters in
    the call list are set for continuing the integration. the user has
    only to call rkf45 again (and perhaps define a new value for tout).
    actually, rkf45 is an interfacing routine which calls subroutine
    rkfs for the solution.  rkfs in turn calls subroutine  fehl which
    computes an approximate solution over one step.

    rkf45  uses the runge-kutta-fehlberg (4,5)  method described
    in the reference
    e.fehlberg , low-order classical runge-kutta formulas with stepsize
                 control , nasa tr r-315

    the performance of rkf45 is illustrated in the reference
    l.f.shampine,h.a.watts,s.davenport, solving non-stiff ordinary
                 differential equations-the state of the art ,
                 sandia laboratories report sand75-0182 ,
                 to appear in siam review.

    the parameters represent-
      f -- subroutine f(t,y,yp) to evaluate derivatives yp(i)=dy(i)/dt
      neqn -- number of equations to be integrated
      y(*) -- solution vector at t
      t -- independent variable
      tout -- output point at which solution is desired
      relerr,abserr -- relative and absolute error tolerances for local
            error test. at each step the code requires that
                 abs(local error) <= relerr*abs(y) + abserr
            for each component of the local error and solution vectors
      iflag -- indicator for status of integration
      work(*) -- array to hold information internal to rkf45 which is
            necessary for subsequent calls. must be dimensioned
            at least  3+6*neqn
      iwork(*) -- integer array used to hold information internal to
            rkf45 which is necessary for subsequent calls. must be
            dimensioned at least  5

  first call to rkf45

    the user must provide storage in his calling program for the arrays
    in the call list  -      y(neqn) , work(3+6*neqn) , iwork(5)  ,
    declare f in an external statement, supply subroutine f(t,y,yp) and
    initialize the following parameters-

      neqn -- number of equations to be integrated.  (neqn >= 1)
      y(*) -- vector of initial conditions
      t -- starting point of integration , must be a variable
      tout -- output point at which solution is desired.
            t=tout is allowed on the first call only, in which case
            rkf45 returns with iflag=2 if continuation is possible.
      relerr,abserr -- relative and absolute local error tolerances
            which must be non-negative. relerr must be a variable while
            abserr may be a constant. the code should normally not be
            used with relative error control smaller than about 1.e-8 .
            to avoid limiting precision difficulties the code requires
            relerr to be larger than an internally computed relative
            error parameter which is machine dependent. in particular,
            pure absolute error is not permitted. if a smaller than
            allowable value of relerr is attempted, rkf45 increases
            relerr appropriately and returns control to the user before
            continuing the integration.
      iflag -- +1,-1  indicator to initialize the code for each new
            problem. normal input is +1. the user should set iflag=-1
            only when one-step integrator control is essential. in this
            case, rkf45 attempts to advance the solution a single step
            in the direction of tout each time it is called. since this
            mode of operation results in extra computing overhead, it
            should be avoided unless needed.

  output from rkf45

      y(*) -- solution at t
      t -- last point reached in integration.
      iflag = 2 -- integration reached tout. indicates successful retur
                   and is the normal mode for continuing integration.
            =-2 -- a single successful step in the direction of tout
                   has been taken. normal mode for continuing
                   integration one step at a time.
            = 3 -- integration was not completed because relative error
                   tolerance was too small. relerr has been increased
                   appropriately for continuing.
            = 4 -- integration was not completed because more than
                   3000 derivative evaluations were needed. this
                   is approximately 500 steps.
            = 5 -- integration was not completed because solution
                   vanished making a pure relative error test
                   impossible. must use non-zero abserr to continue.
                   using the one-step integration mode for one step
                   is a good way to proceed.
            = 6 -- integration was not completed because requested
                   accuracy could not be achieved using smallest
                   allowable stepsize. user must increase the error
                   tolerance before continued integration can be
                   attempted.
            = 7 -- it is likely that rkf45 is inefficient for solving
                   this problem. too much output is restricting the
                   natural stepsize choice. use the one-step integrator
                   mode.
            = 8 -- invalid input parameters
                   this indicator occurs if any of the following is
                   satisfied -   neqn <= 0
                                 t=tout  and  iflag .ne. +1 or -1
                                 relerr or abserr < 0.
                                 iflag == 0  or  < -2  or  > 8
      work(*),iwork(*) -- information which is usually of no interest
                   to the user but necessary for subsequent calls.
                   work(1),...,work(neqn) contain the first derivatives
                   of the solution vector y at t. work(neqn+1) contains
                   the stepsize h to be attempted on the next step.
                   iwork(1) contains the derivative evaluation counter.

  subsequent calls to rkf45

    subroutine rkf45 returns with all information needed to continue
    the integration. if the integration reached tout, the user need onl
    define a new tout and call rkf45 again. in the one-step integrator
    mode (iflag=-2) the user must keep in mind that each step taken is
    in the direction of the current tout. upon reaching tout (indicated
    by changing iflag to 2),the user must then define a new tout and
    reset iflag to -2 to continue in the one-step integrator mode.

    if the integration was not completed but the user still wants to
    continue (iflag=3,4 cases), he just calls rkf45 again. with iflag=3
    the relerr parameter has been adjusted appropriately for continuing
    the integration. in the case of iflag=4 the function counter will
    be reset to 0 and another 3000 function evaluations are allowed.

    however,in the case iflag=5, the user must first alter the error
    criterion to use a positive value of abserr before integration can
    proceed. if he does not,execution is terminated.

    also,in the case iflag=6, it is necessary for the user to reset
    iflag to 2 (or -2 when the one-step integration mode is being used)
    as well as increasing either abserr,relerr or both before the
    integration can be continued. if this is not done, execution will
    be terminated. the occurrence of iflag=6 indicates a trouble spot
    (solution is changing rapidly,singularity may be present) and it
    often is inadvisable to continue.

    if iflag=7 is encountered, the user should use the one-step
    integration mode with the stepsize determined by the code or
    consider switching to the adams codes de/step,intrp. if the user
    insists upon continuing the integration with rkf45, he must reset
    iflag to 2 before calling rkf45 again. otherwise,execution will be
    terminated.

    if iflag=8 is obtained, integration can not be continued unless
    the invalid input parameters are corrected.

    it should be noted that the arrays work,iwork contain information
    required for subsequent integration. accordingly, work and iwork
    should not be altered.
*/
  int k1,k2,k3,k4,k5,k6,k1m;
//     compute indices for the splitting of the work array
      k1m=neqn+0;
      k1=k1m+1;
      k2=k1+neqn;
      k3=k2+neqn;
      k4=k3+neqn;
      k5=k4+neqn;
      k6=k5+neqn;
//     this interfacing routine merely relieves the user of a long
//     calling list via the splitting apart of two working storage
//     arrays. if this is not compatible with the users compiler,
//     he must use rkfs directly.
  rkfs(rhs,neqn,y,t,tout,relerr,abserr,iflag,&work[0],&work[k1m],
    &work[k1],&work[k2],&work[k3],&work[k4],&work[k5],&work[k6],
    &work[k6+1],&iwork[0],&iwork[1],&iwork[2],&iwork[3],&iwork[4]);

  return;
}

void rkf45::fehl(RightHandSide &rhs,int neqn,
  double* y,double &t,double h,double* yp,double* f1,
  double* f2,double* f3,double* f4,double* f5,double* s)
{
//     fehlberg fourth-fifth order runge-kutta method
//
//    fehl integrates a system of neqn first order
//    ordinary differential equations of the form
//             dy(i)/dt=f(t,y(1),---,y(neqn))
//    where the initial values y(i) and the initial derivatives
//    yp(i) are specified at the starting point t. fehl advances
//    the solution over the fixed step h and returns
//    the fifth order (sixth order accurate locally) solution
//    approximation at t+h in array s(i).
//    f1,---,f5 are arrays of dimension neqn which are needed
//    for internal storage.
//    the formulas have been grouped to control loss of significance.
//    fehl should be called with an h not smaller than 13 units of
//    roundoff in t so that the various independent arguments can be
//    distinguished.
  double ch = h/4.0;
  int k;
  for (k=0; k<neqn; k++) f5[k]=y[k]+ch*yp[k];
  rhs(t+ch,f5,f1);
  ch=3.0*h/32.0;
  for (k=0; k<neqn; k++) f5[k]=y[k]+ch*(yp[k]+3.0*f1[k]);
  rhs(t+3.0*h/8.0,f5,f2);
  ch=h/2197.0;
  for (k=0; k<neqn; k++) f5[k]=y[k]+ch*(1932.0*yp[k]+(7296.0*f2[k]-7200.0*f1[k]));
  rhs(t+12.0*h/13.0,f5,f3);
  ch=h/4104.0;
  for (k=0; k<neqn; k++) f5[k]=y[k]+ch*((8341.0*yp[k]-845.0*f3[k])+(29440.0*f2[k]-32832.0*f1[k]));
  rhs(t+h,f5,f4);
  ch=h/20520.0;
  for (k=0; k<neqn; k++) f1[k]=y[k]+ch*((-6080.0*yp[k]+(9295.0*f3[k]-5643.0*f4[k]))+(41040.0*f1[k]-28352.0*f2[k]));
  rhs(t+h/2.0,f1,f5);
  // compute approximate solution at t+h
  ch=h/7618050.0;
  for (k=0; k<neqn; k++) s[k]=y[k]+ch*((902880.0*yp[k]+(3855735.0*f3[k]-1371249.0*f4[k]))+(3953664.0*f2[k]+277020.0*f5[k]));
}

static double d1mach()
{
  static volatile double eps(1), epsp1(1), one(1);
  if(eps==1) {
    do {
      eps /= 2;
      epsp1 = one + eps;
    } while(epsp1 != one);
    eps *= 2;
  }
  return eps;
}

static int sign(int a, int b) { return (b<0)?-abs(a):abs(a); }
static double sign(double a, double b) { return (b<0)?-fabs(a):fabs(a); }
static double max(double a, double b) {return (a>b)?a:b;}
static double min(double a, double b) {return (a<b)?a:b;}

void rkf45::rkfs(RightHandSide &rhs,int neqn,double* y,
  double& t,double& tout,double& relerr,double abserr,int& iflag,
  double* yp,double* h,double* f1,double* f2,double* f3,
  double* f4,double* f5,double* savre,double* savae,int* nfe,
  int* kop,int* init,int* jflag,int* kflag)
{
//     fehlberg fourth-fifth order runge-kutta method
//
//
//     rkfs integrates a system of first order ordinary differential
//     equations as described in the comments for rkf45 .
//     the arrays yp,f1,f2,f3,f4,and f5 (of dimension at least neqn) and
//     the variables h,savre,savae,nfe,kop,init,jflag,and kflag are used
//     internally by the code and appear in the call list to eliminate
//     local retention of variables between calls. accordingly, they
//     should not be altered. items of possible interest are
//         yp - derivative of solution vector at t
//         h  - an appropriate stepsize to be used for the next step
//         nfe- counter on the number of derivative function evaluations
  int hfaild,output;
  double ae,dt,esttol,hmin,remin,s,scale,toln,twoeps,u26;

  int k,mflag;
//  remin is the minimum acceptable value of relerr.  attempts
//  to obtain higher accuracy with this subroutine are usually
//  very expensive and often unsuccessful.
  remin=1.0e-12;
//   here two constants emboding the machine epsilon is present
//   twoesp is set to twice the machine epsilon while u26 is set
//   to 26 times the machine epsilon
//
//     data twoeps, u26/4.4d-16, 5.72d-15/
  twoeps = 2*d1mach();
  u26 = 13*twoeps;
  mflag=abs(iflag);
//     check input parameters
  if(neqn<1||relerr<0||abserr<0||mflag<1||mflag>8) {
    iflag=8;  //   invalid input
    return;
  }

  if(mflag == 1) goto label_50;   // is this the first call

  if(t==tout&& *kflag!=3) {   // check continuation possibilities
    iflag=8;  //   invalid input
    return;
  }
  if(mflag != 2) goto label_25;
//		iflag = +2 or -2
  if(*kflag == 3) goto label_45;
  if(*init == 0) goto label_45;
  if(*kflag == 4) goto label_40;
  if((*kflag == 5)&&(abserr == 0.0)) goto label_30;
  if((*kflag == 6)&&(relerr <= *savre)&&(abserr <= *savae)) goto label_30;
  goto label_50;
//     iflag = 3,4,5,6,7 or 8
  label_25:
  if (iflag == 3) goto label_45;
  if (iflag == 4) goto label_40;
  if ((iflag == 5) && (abserr > 0.0)) goto label_45;
// integration cannot be continued since user did not respond to the instructions pertaining to iflag=5,6,7 or 8
  label_30:
  return; // exit(-1);

  label_40:
  *nfe=0;  // reset function evaluation counter
  if (mflag == 2) goto label_50;
// reset flag value from previous call
  label_45:
  iflag = *jflag;
  if (*kflag == 3) mflag=abs(iflag);
// save input iflag and set continuation flag value for subsequent input checking
  label_50:
  *jflag = iflag;
  *kflag = 0;
// save relerr and abserr for checking input on subsequent calls
  *savre = relerr;
  *savae = abserr;
//  restrict relative error tolerance to be at least as large as 2*eps+remin 
//  to avoid limiting precision difficulties arising from impossible accuracy requests
  if(relerr < twoeps+remin) {
    relerr = twoeps+remin;
    iflag = *kflag = 3;//   relative error tolerance too small
    return;
  }
  dt=tout-t;

  if (mflag == 1) goto label_60;
  if (*init == 0) goto label_65;
  goto label_80;
// initialization --
//    set initialization completion indicator, init; set indicator for too many output points, kop
//    evaluate initial derivatives set counter for function evaluations,nfe
//    estimate starting stepsize
  label_60:
  *init=0;
  *kop=0;

  rhs(t,y,yp);
  *nfe = 1;
  if (t == tout) { 
    iflag = 2;
    return;
  }

  label_65:
  *init = 1;
  *h = fabs(dt);
  toln = 0;
  for(k=0; k<neqn; k++) {
    double tol=relerr*fabs(y[k])+abserr;
    if(tol > 0) {
      toln=tol;
      double ypk=fabs(yp[k]);
      if (ypk*(*h)*(*h)*(*h)*(*h)*(*h) > tol) *h = pow(tol/ypk,0.2);
    }
  }
  if(toln <= 0) *h = 0.0;
  *h = max(*h,u26*max(fabs(t),fabs(dt)));
  *jflag = sign(2,iflag);
//     set stepsize for integration in the direction from t to tout
  label_80:
  *h = sign(*h,dt);
//  test to see if rkf45 is being severely impacted by too many output points
  if(fabs(*h)>=2.0*fabs(dt)) (*kop)++;
  if((*kop) == KOP_MAX) {
    *kop = 0;
    iflag = 7;
    return;  // unnecessary frequency of output
  }
//   set smallest allowable stepsize
  hmin=u26*fabs(t);
  if (fabs(dt) <= hmin) { // if too close to output point,extrapolate and return
    for (k=0; k<neqn; k++) y[k]+=dt*yp[k];
    t=tout;
    rhs(t,y,yp);
    (*nfe)++;
    iflag=2;
    return;
  }
//   initialize output point indicator
  output= 0;
//   to avoid premature underflow in the error tolerance function, scale the error tolerances
  scale=2.0/relerr;
  ae=scale*abserr;
//   step by step integration
  label_100:
  hfaild= 0;
//   adjust stepsize if necessary to hit the output point.
//   look ahead two steps to avoid drastic changes in the stepsize and
//   thus lessen the impact of output points on the code.
  dt=tout-t;
  if(fabs(dt)<2*fabs(*h)) {
    if(fabs(dt)<=fabs(*h)) {
      output = 1; // the next successful step will complete the integration to the output point
      *h = dt;
    } 
    else {
      *h = 0.5*dt;
    }
  }
//     core integrator for taking a single step
//
//     the tolerances have been scaled to avoid premature underflow in
//     computing the error tolerance function et.
//     to avoid problems with zero crossings,relative error is measured
//     using the average of the magnitudes of the solution at the
//     beginning and end of a step.
//     the error estimate formula has been grouped to control loss of
//     significance.
//     to distinguish the various arguments, h is not permitted
//     to become smaller than 26 units of roundoff in t.
//     practical limits on the change in the stepsize are enforced to
//     smooth the stepsize selection process and to avoid excessive
//     chattering on problems having discontinuities.
//     to prevent unnecessary failures, the code uses 9/10 the stepsize
//     it estimates will succeed.
//     after a step failure, the stepsize is not allowed to increase for
//     the next attempted step. this makes the code more efficient on
//     problems having discontinuities and more effective in general
//     since local extrapolation is being used and extra caution seems
//     warranted.
//
//     test number of derivative function evaluations.
//     if okay,try to advance the integration from t to t+h

  for(;;) {  
    if (*nfe > MAX_NFE) {
      iflag = *kflag = 4;  //  too much work
      return;  
    }
//     advance an approximate solution over one step of length h
    fehl(rhs,neqn,y,t,*h,yp,f1,f2,f3,f4,f5,f1);
    *nfe += 5;
//     compute and test allowable tolerances versus local error estimates
//     and remove scaling of tolerances. note that relative error is
//     measured with respect to the average of the magnitudes of the
//     solution at the beginning and end of the step.
    double eeoet=0;
    for (k=0; k<neqn; k++) {
      double et=fabs(y[k])+fabs(f1[k])+ae;
      if(et<=0) {
        iflag=5;  // inappropriate error tolerance
        return;  
      }
      double ee=fabs((-2090.0*yp[k]+(21970.0*f3[k]-15048.0*f4[k]))+(22528.0*f2[k]-27360.0*f5[k]));
      eeoet=max(eeoet,ee/et);
    }
    esttol=fabs(*h)*eeoet*scale/752400.0;
    if(esttol <= 1.0) break;
// unsuccessful step; reduce the stepsize and try again; the decrease is limited to a factor of 1/10
    hfaild= 1;
    output= 0;
    s=0.1;
    if (esttol < 59049.0) s=0.9/pow(esttol,0.2);
    *h *= s;
    if(fabs(*h) > hmin) continue;
    iflag = *kflag = 6; // requested error unattainable at smallest allowable stepsize
    return;
  }
// successful step; store solution at t+h and evaluate derivatives there
  t += (*h);
  for (k=0; k<neqn; k++) y[k]=f1[k];
  rhs(t,y,yp);
  (*nfe)++;
// choose next stepsize; the increase is limited to a factor of 5 
// if step failure has just occurred, next stepsize is not allowed to increase
  s=5.0;
  if (esttol>1.889568e-4) s=0.9/pow(esttol,0.2);
  if (hfaild) s=min(s,1.0);
  *h = sign(max(s*fabs(*h),hmin),*h);
//     end of core integrator
//     should we take another step
  if (output) { 
    t=tout;
    iflag=2;  //  interval mode
    return;
  }
  if (iflag > 0) goto label_100;
//     integration successfully completed
//     one-step mode
  iflag=-2;
}

};   //namespace MConf

