#include <math.h>
#include <stdlib.h>
#include <vector>
#include <iostream>

static const double ATOL=1e-8;
static const double RTOL=1e-9;

namespace MConf {

//****************************************************************************
// class RightHandSide must provide the operator() (double t, double *y, double *f)
//  describing the system of ODE dy/dt = f(t,y)
class RightHandSide {
public:
  virtual void operator() (double fi,const double *y,double *yp) = 0; 
};

//****************************************************************************
///  Fehlberg fourth-fifth order Runge-Kutta method
class rkf45
{
public:
  rkf45() { init(0,1); }
  rkf45(RightHandSide *rhs, int neqn,double atol=ATOL,double rtol=RTOL) {
    init(rhs,neqn,atol,rtol);
  }
  void init(RightHandSide *rhs,int neqn, double atol=ATOL,double rtol=RTOL) {
    this->a_tol=atol;
    this->r_tol=rtol;
    this->neqn=neqn;
    this->rhs=rhs;
    iflag = 1;
    y.reserve(neqn);
    y.resize(neqn);
    work.reserve(3+6*neqn);
    work.resize(3+6*neqn);
  }
  void init() {
    iflag = 1;
  }
  virtual ~rkf45() {;}

// advance u from t to tout
  void advance(double *u, double &t, double tout)
  {
    double * wrk = &(work[0]);
    rkf_45 (*rhs,neqn,u,t,tout,r_tol,a_tol,iflag,wrk,iwork);
    if(iflag>2) {  
      std::cerr<<"rk45_flag="<<iflag<<std::endl;
      iflag = 2;
    }    
  }

private:
// set start point
  void setStart(double &t,double *u) 
  {
    t; // to avoid warning: unreferenced parameter
    for(size_t i=0; i<y.size(); ++i) y[i] = u[i];
    iflag = 1;
  }
// advance u from t to tout
  void advance1(double *u, double &t, double tout)
  {
    if(iflag == 1) {
      for(size_t i=0; i<y.size(); ++i) y[i] = u[i];
    }
    rkf_45 (*rhs,neqn,&(y[0]),t,tout,r_tol,a_tol,iflag,&(work[0]),iwork);
    if(iflag>2) { 
      std::cerr<<"rk45_flag="<<iflag<<std::endl;
      iflag = 2;
    }    
    for(size_t i=0; i<y.size(); ++i) u[i] = y[i];
  }

private:
// class RightHandSide must provide the operator() (double t, double *y, double *f)
//  describing the system of ODE dy/dt = f(t,y)
  RightHandSide *rhs; 
  int neqn,iflag;
  double a_tol;
  double r_tol;
  std::vector<double> y;        // array[neqn]
  std::vector<double> work;     // array[3+6*neqn]
  int iwork[5]; 

  void fehl(RightHandSide &rhs,int neqn,
            double* y,double &t,double h,double* yp,double* f1,
            double* f2,double* f3,double* f4,double* f5,double* s);

  void rkfs(RightHandSide &rhs,int neqn,double* y,
            double& t,double& tout,double& relerr,double abserr,int& iflag,
            double* yp,double* h,double* f1,double* f2,double* f3,
            double* f4,double* f5,double* savre,double* savae,int* nfe,
            int* kop,int* init,int* jflag,int* kflag);

  void rkf_45(RightHandSide &rhs,int neqn,double* y,
	          double& t,double& tout,double& relerr,double abserr,int& iflag,
	          double* work,int* iwork);

};

};   //namespace MConf

