#include "../include/CStconfig.h"
#include "../include/CEfit.h"
#include "../include/threads.h"

#include <time.h>


#define PRINT_POINCARE 0

#define PRINT_TIME 0

static const bool timePrintFlag = PRINT_TIME != 0;

namespace MConf {


#define  NTHREADS  8

//****************************************************************************
// Read file in EFIT-format.
//
bool CStconfig::loadEfit_p(const char * fname, double scaleBpol, double scaleBtor,int signBpol,int signBtor, 
                            int signQ, int psiOverTwopi, const CStconfig *mconf, long fPostn)
{
#if PRINT_TIME
  mathf::StopWatch watch;
  watch.setPrint(timePrintFlag);
  watch.start();
#endif
  
  int m = 128;
  int ns= 100; 
  int nLastToRemove = 1;

  CEfit efit(m,ns,nLastToRemove);
  if(mconf) {
    if(!efit.load(mconf)) return(mOK=false);
    mfname = mconf->fname();
    mfname += "-as_circular_tokamak";
  }
  else {
    if(!efit.load(fname,scaleBpol,scaleBtor,signBpol,signBtor,signQ,psiOverTwopi,fPostn)) return(mOK=false);
    mfname = fname;
  }
 
#if PRINT_TIME
  watch.printTime("loadEfit::loadEfit_p ",true);
#endif

// now we import equilibria from efit-object
  ns = efit.getNs();
  m  = efit.getM();

  mNp=1;
  if(!resize(ns,m,1)) return(mOK=false); // mNs0 will be assigned  
#if PRINT_POINCARE
  std::ofstream out("trace.prn");
  std::ofstream *outp=&out;
#else
  std::ofstream *outp=NULL;
#endif

#if  NTHREADS > 0
  calcEfitFourierMods(&efit);
#else
  for(int i=0; i<mNs0; i++) {
    int is = m1stSindex+i;
    double dVds, Itor;
    efit.createSurface(i,dVds,Itor,outp);
    ms.rw()[is]      = efit.getS(i);
    miota.rw()[is]   = efit.getIota(i);
    mpprime.rw()[is] = -efit.getPp(i);

    mIpol.rw()[is]   = efit.getIpol(i);
    mItor.rw()[is]   = Itor;
    mg00.rw()[is]    = dVds;

    for(int m=0; m<=mM; m++) {
      Rmn(is,m,0)   = efit.getRm0(m);
      Zmn(is,m,0)   = efit.getZm0(m);
      Phimn(is,m,0) = 0;
      bmn(is,m,0)   = efit.getBm0(m);

      Rmn(is,m,1)   = efit.getRm1(m);
      Zmn(is,m,1)   = efit.getZm1(m);
      Phimn(is,m,1) = 0;
      bmn(is,m,1)   = efit.getBm1(m);
    }
  }
#endif


#if PRINT_TIME
  watch.printTime("loadEfit::calcEfitFourierMods ");
#endif

  mFlux = fabs(efit.getFlux());
  mFlux *= twopi;
  mR0 = efit.getRmajor();
  ma  = sqrt(fabs(mFlux)/bmn(mLCMSindex,0,0)/pi);

  if(0) {
    mR0=Rmn(m1stSindex,0,0);
  }

#if PRINT_POINCARE
 {
    std::ofstream ooo("B_efit.prn");
    double r1=4.25;    //ITER
    double r2=8.18;
    double z=0.681;
//    double r1=0.65;  //TCV
//    double r2=1.1;
//    double z=0;
    for(int i=0;i< 100; i++) {
      double r = r1+i*(r2-r1)/99;
     Vector3d c(r,0,z);
     ooo<<c<<"  "<<efit.M3DgetBcyl(c)<<std::endl;
    }
  }
#endif
  Vector3d axis = efit.getMAxis();

 //// vprimeNeeded = true;     //for test only: signal to load-function: we need vprime(mg00)
  ////torCurrentNeeded = true; //for test only: signal to load-function:  need toroidal current
  mNewFormat = true;
  tokamakSymmetryFluxCoordinates = true;
  tokamakEQDSK   = true;
  return(mOK=true);
}


//****************************************************************************
void CStconfig::calcEfitFourierMods(void * efit)
{
  Threads2<CStconfig,CEfit> threads(this, &CStconfig::calcEfitFourierModsExe, (CEfit*)efit,"calcEfitFourierMods");  
  int low = 0, upper = mNs0-1;
  threads.run(low, upper, NTHREADS); 
}

//template <class efit> void CStconfig::calcEfitFourierModsExe(efit* ef,int low, int upper, bool setup)
void CStconfig::calcEfitFourierModsExe(CEfit* ef,int low, int upper, bool setup)
{
  if(setup) {  // setup part; it allocates arrays and initialize them
    int size = upper-low+1;
  }
  else
  {  // exe part; it will be called from threads to fill arrays
    for(int i=low; i<=upper; i++) {
      double dVds, Itor;
      ef->createSurface(i,dVds,Itor, 0);
      int is = m1stSindex+i;
      ms.w()[is]      = ef->getS(i);
      miota.w()[is]   = ef->getIota(i);
      mpprime.w()[is] = -ef->getPp(i);

      mIpol.w()[is]   = ef->getIpol(i);
      mItor.w()[is]   = Itor;
      mg00.w()[is]    = dVds;

      for(int m=0; m<=mM; m++) {
        Rmn(is,m,0)   = ef->getRm0(m);
        Zmn(is,m,0)   = ef->getZm0(m);
        Phimn(is,m,0) = 0;
        bmn(is,m,0)   = ef->getBm0(m);
        Rmn(is,m,1)   = ef->getRm1(m);
        Zmn(is,m,1)   = ef->getZm1(m);
        Phimn(is,m,1) = 0;
        bmn(is,m,1)   = ef->getBm1(m);
      }
    }
  }
}



};   //namespace MConf
