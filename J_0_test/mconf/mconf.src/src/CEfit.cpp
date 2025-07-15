#include "../include/CStconfig.h"
#include "../include/CEfit.h"
#include "../include/threads.h"
#include <time.h>
#include <string.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <sstream>

#define FIND_PSI_EDGE 0

#define DENORMALIZE_JET_PSIRZ 1  

#define PRINT_TIME  0

#define PRINT_ARRAYS 0


static const bool timePrintFlag = PRINT_TIME != 0;

namespace MConf {

#define  NTHREADS  8

static const double RKSTEP=1./1024;     // 0.056 degree
//static const double RKSTEP=1./2048;   // 0.028 degree

//****************************************************************************
CEfit::CEfit(int m, int n, int nLastToSkip)
{
  Mpol = m;
  NSurfaces = n;
  nLastToRemove = nLastToSkip;
  meshOK=false;
  pi    = 3.14159265358979323846;
  twopi = 2*pi;
  Rmin=Rmax=Zmin=Zmax=0;
  doPositionCheck=true;
  sMax   = 1.1;
}

//****************************************************************************
CEfit::CEfit()
{
  Mpol = 128;
  NSurfaces = 103;
  nLastToRemove = 2;
  meshOK=false;
  pi    = 3.14159265358979323846;
  twopi = 2*pi;
  Rmin=Rmax=Zmin=Zmax=0;
  doPositionCheck=true;
  sMax   = 1.1;
}

//****************************************************************************
void CEfit::clear()
{
  psi_.clear();
  Fp_.clear(); 
  FFp_.clear();
  Pp_.clear(); 
  q_.clear();  
  //s_.clear();      
  R.clear();         
  Z.clear();         
  psirz.clear();     
  mNormFlux.clear(); 
  LCMS.clear();      
  RBZm_.clear();     
  RBZm1_.clear();    
  RBZt_.clear();     
  RBZt_smooth.clear();
  mR_.clear();      
  mpsi_.clear();    
  ms_.clear();      
  //miota_.clear(); 
  mq_.clear();
  meshOK=false;
}

//****************************************************************************
bool CEfit::resize()
{
  clear();
  meshOK=true;
  R.init(nR);
  Z.init(nZ);
//  s_.init(nR);
  psi_.init(nR);
  Fp_.init(nR);
  FFp_.init(nR);

  q_.init(nR);
  Pp_.init(nR);
  meshOK  = psirz.resize(nR,nZ);
  meshOK &= mNormFlux.resize(nR,nZ);
  if(!meshOK) clear();
  return meshOK;
}

//****************************************************************************
// psiSign==1 if poloidal flux has minimum in the center
//  -1 otherwise
int CEfit::getPsirzSign() 
{
// find sign/direction of poloidal flux; we need this for amoebaFunc 
  dR   = xdim/(nR-1);
  dZ   = zdim/(nZ-1);
  int iz = int((Zaxis - Zmin)/dZ);             // grid-index close to magnetic axis
  int ir = int((Raxis - Rmin)/dR);             // grid-index close to magnetic axis
  double psiAx = psirz(ir,iz);  // poloidal flux at magnetic axis
  psiSign = 0;
  for(int m=1; m<12; m++) {         // we assume at least 11 points till LCMS
    double psiEd = psirz(ir-m, iz); // get poloidal flux while moving to the plasma edge
    psiSign+=(psiAx<psiEd)?1:-1;     
  }
  psiSign=(psiSign>0)?1:-1; // psiSign==1 if poloidal flux has minimum in the center
  return psiSign;
}

//****************************************************************************
// find psi at LCMS
Vector3d CEfit::findPsiEdge(double &qEdge, double &psi_edge ,double &r_edge, double &R1, double &R2)    
{
  double DR = R2-Raxis;
  R1 = Raxis+0.007*DR;
  R2 = R2-0.00001*DR;

  double qGuess = qEdge==0?q_[nR-1]:qEdge;
// find LCMS
//2015_1124  if(traceForQ(Vector3d(R2,0,Zaxis),qGuess)==0) 
// if R2 outside LCMS then find correct R2 by bisection
  double rEdge    = R2;        // Edge   point
  double rInside  = R2-0.1*DR; // inner  point
  rInside  = R1; // inner  point
  int itr=200;
  for(; itr; --itr) {
    double r = (rInside+rEdge)/2;
    double qE = traceForQ(Vector3d(r,0,Zaxis),qGuess);
    qEdge = qE;
    if(qE>1e30)
      rEdge   = r;
    else {
      rInside = r;      // if r lies inside
      if(fabs(rEdge-r)<=Raxis*1e-7) break;  // Ok, done
    }
  }
  R2    = rInside;
  r_edge = rInside;
  Vector3d edgeCoord(r_edge,0,Zaxis);
  psi_edge = cyl2psi(edgeCoord);
  return getBcyl(edgeCoord);
}

//****************************************************************************
//
bool CEfit::set(int ns, int m, CEfit::fileFormat fileFmt)
{
  mathf::StopWatch watch;
  watch.setPrint(timePrintFlag);
  watch.start();
  M = m;
  Ns = ns;

  Ntrace = 16*M;
  RBZm_.init(M+1);
  RBZm1_.init(M+1);
  RBZt_.init(Ntrace);
//  RBZt_smooth.init(Ntrace);

  mR_.init(ns);
  ms_.init(ns);
  mq_.init(ns);
//  miota_.init(ns);
  mpsi_.init(ns);

// Redge was found in loadx by bisection
  double R1=Raxis,R2=Redge;
#if FIND_PSI_EDGE
  double qEdge = q_[nR-1], psi_edge, r_edge;
  findPsiEdge(qEdge, psi_edge, r_edge, R1, R2);
  Redge = r_edge;
  watch.printTime("CEfit::set--findPsiEdge ",true);
#endif
  double dr = (R2-R1);
  R1 = Raxis + dr*0.03;
  if(fileFmt==P_EQDSK   ) R2 -=dr*0.002;  // because of inaccurate format
  if(fileFmt==P_EQDSKold) R2 -=dr*0.005;  // because of inaccurate format

  int i;
  for(int i=0; i<Ns; i++) {  // create new grid
    double x = double(i)/(Ns-1);
    mR_[i]   = R1+(R2-R1)*mapFunction(x);    //more fine grid at  x=0 and x=1
    mpsi_[i] = cyl2psi(Vector3d(mR_[i],0,Zaxis));
  }

  double R2x = mR_[Ns-1] ;
  mR_[Ns-1]   = R2;
  mpsi_[Ns-1] = cyl2psi(Vector3d(mR_[Ns-1],0,Zaxis));

  if(fileFmt==C_TOKAMAK) {
 //   for(i=0; i<Ns; i++) miota_[i] = 1/q0;
    for(int i=0; i<Ns; i++) mq_[i] = q_[0];
  }
  else 
  {
#if  NTHREADS == 0
    for(i=0; i<Ns; i++) {// calculate q by B-field line tracing
      double qGuess = q(mpsi_[i]);
      mq_[i] = traceForQ(Vector3d(mR_[i],0,Zaxis),qGuess);
      if(mq_[i]>1e30) break;
    }
#else
    calculateQ();
#endif
    int NsCopy = Ns;
    for(int i=NsCopy/2; i<NsCopy; i++) { // skip points lying outside LCMS
      if(mq_[i]>1e30) break;
      Ns = i;
    }
  }

  if(!(fileFmt==C_TOKAMAK||fileFmt==C_TOKAMAK_EQUIVALENT)) Ns -= nLastToRemove;  //for safety if bad behaviour at the edge

  watch.printTime("CEfit::set--calculateQ() ",true);


// calculate toroidal flux using q calculated in previous loop
  ms_[0]=mFluxTor(psiAxis,mpsi_[0]);
  for(i=1; i<Ns; i++) ms_[i]=ms_[i-1]+mFluxTor(mpsi_[i-1],mpsi_[i]);
  fluxTorMax = ms_[Ns-1]+mFluxTor(mpsi_[Ns-1],psiEdge);
  fluxTorMax = ms_[Ns-1];
//  std::cerr<<"edge flux skiped: "<<(psiEdge-mpsi_[Ns-1])/(psiEdge-psiAxis)<<std::endl;

  for(i=0;i< Ns; i++) ms_[i] /= fluxTorMax;

#if PRINT_ARRAYS
  std::ofstream out("MY_spol_s_q_mr.prn");
  for(i=0;i<Ns;i++) out<<mpsi_[i]<<"  "<<ms_[i]<<"  "<<mq_[i]<<"   "<<mR_[i]<<std::endl;
#endif

return true;
#if 0 
    std::ofstream ooo("traceIterEFIT202");
    ooo<<"r    x      Br,   Bfi,   Bz     |B|  1/q"<<std::endl;
    ooo.precision(7);
    ooo.setf(std::ios::scientific);
    for(i=0;i< 200; i++) {
      double r1=4.34;
      double r2=8.12;
      double r = r1+i*(r2-r1)/199;
      double z = 0.;
      Vector3d c(r,0,z);
    // ms_[]-mesh must be ready before calling cyl2s()
      double x = sqrt(cyl2s(c));
      double psi;
      Vector3d B = getBcyl(c, &psi);
      double iot = 1/q(psi);
      ooo<<r<<"  "<<x<<"  "<<getBcyl(c)<<"  "<<getBcyl(c).abs()<<"  "<<iot<<std::endl;
    }

//return true;

    double torFlux=FluxTor(psiAxis,mpsi_[0]);
    for(i=1; i<Ns; i++) torFlux +=FluxTor(mpsi_[i-1],mpsi_[i]);
    torFlux += FluxTor(mpsi_[Ns-1],psiEdge);


  std::ofstream out("sriota202.prn");
  out.precision(7);
  out.setf(std::ios::scientific);
  for(i=0;i< Ns; i++)
    out<<"  "<<ms_[i]<<"  "<<mR_[i]<<"  "<<1/mq_[i]<<"  "<<1/q(mpsi_[i])<<"   "<<(mpsi_[i]-psiAxis)/(psiEdge-psiAxis)<<std::endl;
return true;

  std::cerr<<cyl2psi(Vector3d(Redge,0,Zaxis     ))<<std::endl;
  std::cerr<<cyl2psi(Vector3d(R2,0,Zaxis        ))<<std::endl;
  std::cerr<<cyl2psi(Vector3d(R2-0.01*DR,0,Zaxis))<<std::endl;
  std::cerr<<cyl2psi(Vector3d(R2-0.02*DR,0,Zaxis))<<std::endl;
  return true;
#endif 
}

//****************************************************************************
void CEfit::calculateQ()
{
  Threads2 < CEfit,CEfit > threads(this, &CEfit::calculateQExe,this,"calculateQ");  
  int low = 0, upper = Ns-1;
  threads.run(low, upper, NTHREADS); 
}

//template <class efit> void CEfit::calculateQExe(efit * ef,int low, int upper, bool setup)
void CEfit::calculateQExe(CEfit * ef,int low, int upper, bool setup)
{
  if(setup) {  // setup part; it allocates arrays if needed
    int size = upper-low+1;
  }
  else
  {  // exe part; it will be called from threads to fill arrays
    for(int i=low; i<=upper; i++) {
      double p = ef->mpsi_[i];
      double qGuess = ef->q(p);
      double R = ef->mR_[i];
      Vector3d r(R,0,Zaxis);
      double q = ef->traceForQ(r,qGuess);
      mq_.wd(i) = q;
      if(q>1e30) break;
    }
  }
}

//****************************************************************************
// Ascertain(find out) the file format
//
CEfit::fileFormat CEfit::getfileformat(const char * fullname,long fPostn) const
{
  // we use this table to find out the format of a file
  static struct {
    CEfit::fileFormat fmt;  const char * str;
  } pattern[] = {
    {CEfit::C_TOKAMAK,"circular tokamak"}, // B0, R0, a, q  are given
    {CEfit::P_EQDSK,"EQDSK by A. Portone"},  // "free format": five numbers with blanks between them;
                                             // http://efdasql.ipp.mpg.de/EUHandCD/ECRH/equilibria_index.htm
    {CEfit::J_EQDSK,"JET"},   // JET format uses normalized psirz(i,j) with 0 on axis and 1 at the separatrix
    {CEfit::G_EQDSK,"EQDSK"}, // see G EQDSK file format(5(e16.9)); five numbers without delimeters between them
    {CEfit::G_EQDSK,"LIUQE"}, 
    {CEfit::G_EQDSK,"eqdsk"},
    {CEfit::D_EQDSK,"EFITD"}, 
    {CEfit::G_EQDSK,"EFIT"},
    {CEfit::G_EQDSK,"CHEASE"},
    {CEfit::G_EQDSK,"micdu"}
   };

  fileFormat fmt=UNDEFINED;
  FILE * fp = fopen(fullname,"rb");
  if(fp==NULL) return fmt;

  fseek(fp,fPostn,SEEK_SET); //if(fPostn!=0)  fsetpos(fp,&fPostn);    // set position to read

  char buff[256];
  int i, sz=sizeof(pattern)/sizeof(pattern[0]);

// read character strings and compare with the table records
  int k=50;
  while(k--) {
    if(fgets(buff,255,fp)==NULL) goto ret;
    if(strlen(buff)<4) continue;
    for(i=0; i<sz; i++)
      if(strstr(buff,pattern[i].str)) {
        fmt=pattern[i].fmt;
        goto ret;  // OK, we got the format
      }
  }
ret:
  if(fmt==UNDEFINED) fmt=G_EQDSK; // try  G_EQDSK format by default
  if(fmt==P_EQDSK) {
    std::istringstream str(buff);
    std::string name,date;
    str>>name>>date;
    //if(date<"Mar-2005") fmt==P_EQDSKold
  }
  fclose(fp);
  return fmt;
}

//****************************************************************************
bool CEfit::load(const char * fname, double scaleBpol, double scaleBtor, int signBpol,int signBtor, 
                  int signQ, int psiOverTwopi, long fPostn)
{
  mathf::StopWatch watch;
  watch.setPrint(timePrintFlag);
  watch.start();
  CEfit::fileFormat fileFmt = getfileformat(fname,fPostn);
  if(loadx(fname,fileFmt,scaleBpol,scaleBtor,signBpol,signBtor,signQ,psiOverTwopi,0,fPostn)) {
    if(fileFmt==C_TOKAMAK) NSurfaces = 64;
    if(fileFmt==C_TOKAMAK) Mpol  = 25;
    watch.printTime("CEfit::loadx--total ",true);
    bool result=set(NSurfaces,Mpol,fileFmt);
    watch.printTime("CEfit::set--total ");
    return result;
  }
  else
    return false;
}

//****************************************************************************
// Create circular tokamak using a,R0, B00, iota(0) from given configuration mconf
bool CEfit::load(const CStconfig *mconf)
{
  NSurfaces=64;
  Mpol=25;
  double scaleBpol=1, scaleBtor=1;
  mathf::StopWatch watch;
  watch.setPrint(timePrintFlag);
  watch.start();
  CEfit::fileFormat fileFmt = C_TOKAMAK_EQUIVALENT;
  if(loadx("C_TOKAMAK_EQUIVALENT",fileFmt,0,0,0,0,0,0,mconf)) {
    watch.printTime("CEfit::loadx--total ",true);
    bool result=set(NSurfaces,Mpol,fileFmt);
    watch.printTime("CEfit::set--total");
    return result;
  }
  else
    return false;
}

static double getNumber(std::string &buf, const char * numName, double defaultVal=0) 
{
  double num=defaultVal;  //return defaultVal if numName not found
  size_t j = buf.find(numName);
  if(j!=std::string::npos) {
    std::string bu(buf.c_str()+j+strlen(numName));
    if(strcmp("COCOS",numName)==0) bu.resize(2);  // to handle case COCOSxx, where xx are two digits 
    std::istringstream str(bu);
    str >>num;  
  }
  return num;
}

//****************************************************************************
// @return  -1 if "no"; 1 if "yes"; 0 if not found
static int getBool(std::string &buf, const char * valName, int defaultVal=0) 
{
  std::string value;
  size_t j = buf.find(valName);
  if(j!=std::string::npos) {
    std::istringstream str(buf.c_str()+j+strlen(valName));
    str >>value;  
  }
  if(value.empty()) return defaultVal; //return defaultVal if valName is not found
  else if(value=="yes") return 1;
  else if(value=="no") return -1;
  else return defaultVal;
}

//****************************************************************************
// this file is called if (fileFmt==C_TOKAMAK||fileFmt==C_TOKAMAK_EQUIVALENT)
bool CEfit::loadCircularTokamak(const char * fname, CEfit::fileFormat fileFmt, const CStconfig *mconf, long fPostn)
{
  //if(!(fileFmt==C_TOKAMAK||fileFmt==C_TOKAMAK_EQUIVALENT)) return meshOK=false;
  double B0(0), R0(0), a0(0), q0(0);  // for circular tokamak

  if(fileFmt==C_TOKAMAK) {
    std::ifstream in(fname,std::ios::in|std::ios::binary); //std::ifstream::binary
    if(!in.is_open())  return(meshOK=false);
    if(fPostn>=0) in.seekg(fPostn, std::ios::beg); 

    std::string str;
    std::getline(in,str);   
    while(str.compare(0, 8, "circular")!=0) {
      std::getline(in,str);
    } //  in.ignore(1000,'\n');  // ignore the next line
   // from this file we need only B0, R0, a0, q0
    in>>B0>>R0>>a0>>q0;
    a0 = mmin(a0, 0.9*R0);	
  } else if(fileFmt==C_TOKAMAK_EQUIVALENT) {
    if(!mconf->isOK()) return meshOK=false;
   // from mconf we need only B0, R0, a0, q0
    B0 = mconf->Bmn_sp(0,0,0);
    R0 = mconf->R0();
    a0 = mconf->r(1.);
    a0 = mmin(a0, 0.9*R0);
    q0 = 1/fabs(mconf->iota(0));
  }
  else
    return meshOK=false;
 
  nR = nZ = 513;
  xdim = zdim = 2*a0/0.9;  
  Rmajor  = R0;
  Raxis   = R0;   
  Zaxis   = 0;
  zmid    = 0;
  Rmin    = R0 - xdim/2;
  psiAxis = 0;
  psiEdge = -B0*R0*(R0-sqrt(R0*R0-a0*a0))/q0;     // ==Psi/twopi  [Weber/rad], see EFIT doc.

  if(!resize()) return meshOK=false;

  for(int i=0;i<nR;i++) {
    q_[i]  = q0;
    Fp_[i] = R0*B0; 
    FFp_[i] = 0; 
    Pp_[i] = 0; 
  }
  // create rectangular mesh
  dR   = xdim/(nR-1);
  dZ   = zdim/(nZ-1);
  edR = 1/dR;
  edZ = 1/dZ;
  for(int i=0;i< nR; i++) R[i] = Rmin+i*dR; Rmax=R[nR-1];
  for(int k=0;k< nZ; k++) Z[k] = zmid - zdim/2.0 + k*dZ; Zmin = Z[0]; Zmax=Z[nZ-1];

  if(fileFmt==C_TOKAMAK) {
    for(int k=0;k<nZ;k++)
      for(int i=0;i<nR;i++) {
        double r2 = square(R[i]-Raxis)+square(Z[k]-Zaxis);
        psirz(i,k) = -B0*R0*(R0-sqrt(R0*R0-r2))/q0;
      }
  }
  if(fileFmt==C_TOKAMAK_EQUIVALENT) {
  // create poloidal flux psi(r^2)
    int Nr = 4097;
    pBase::ngArray<double> r2(Nr);
    pBase::ngArray<double> psi(Nr);
    pBase::ngArray<double> iotaExt(Nr);
    double dr2 = Zmax*Zmax*2.1/(Nr-1);
    r2[0] = 0;
    psi[0] = 0;
    iotaExt[0] = fabs(mconf->iota(0));
    for(int k=1; k<Nr; k++) {  
      r2[k] = k*dr2;
      double s = r2[k]/(a0*a0);
      iotaExt[k] = fabs(mconf->iota(s<1.1?s:1.1));
      double dtorflux = B0*R0*( sqrt(R0*R0-r2[k-1])-sqrt(R0*R0-r2[k]) ); // =dTorFlux
      double dpsi = iotaExt[k]*dtorflux; // integral iota(r) dTorFlux 
      //double dpsi = dr2/2*B0*R0*iotaExt[k]/sqrt(R0*R0-r2[k]); // integral iota(r) dTorFlux 
      psi[k] = psi[k-1] - dpsi;
    }
    psiAxis = 0;
    psiEdge = interp2(a0*a0,psi,r2,Nr);

    for(int i=0;i<nR;i++) {  //PSI/twopi on the uniform poloidal flux grid
      double dpsi = (psiEdge-psiAxis)/(nR-1);   //for(int i=0;i<nR;i++) psi_[i]=psiAxis+i*dpsi;   //uniform poloidal flux grid
      double r2_psi = interp2(psiAxis+i*dpsi,r2,psi,Nr); // find r^2 which correspons to the  (psiAxis+i*dpsi)
      q_[i]  = 1/interp2(r2_psi,iotaExt,r2,Nr);
    }
    for(int j=0;j<nZ;j++)
      for(int i=0;i<nR;i++) {
        double r2p = square(R[i]-Raxis)+square(Z[j]-Zaxis);
         psirz(i,j) = interp2(r2p,psi,r2,Nr);
      }
  }
  Bcenter = Fp_[0]/Raxis;
  psiSign = getPsirzSign();
// uniform poloidal flux grid psi_ will be created later
  return meshOK=true;
}

//****************************************************************************
bool CEfit::loadx(const char * fname, CEfit::fileFormat fileFmt, double scale_Bpol,double scale_Btor,int sign_Bpol,int sign_Btor, 
                  int sign_Q, int psiOverTwopi, const CStconfig *mconf, long fPostn)
{
  std::ifstream in;
  std::string name,date;
  int i,k;
  int psiOver2pi = 1;
  double COCOS = 1;
  double scalBpol = 1;
  double scalBtor = 1;
  double signBpol = 0;
  double signBtor = 0;
  double signQ = 0;

  mathf::StopWatch watch;
  watch.setPrint(timePrintFlag);
  watch.start();

  if(fileFmt==C_TOKAMAK||fileFmt==C_TOKAMAK_EQUIVALENT) { 
    loadCircularTokamak(fname, fileFmt, mconf, fPostn);
  }
  else  {  // EQDSK file; start reading
    in.open(fname,std::ios::in|std::ios::binary);
    if(!in.is_open())  return(meshOK=false);
    if(fPostn>=0) in.seekg(fPostn, std::ios::beg); 
    
    std::string buf;
    std::getline(in,buf);

    scalBpol = getNumber(buf,"BpolScale=",1);
    scalBtor = getNumber(buf,"BtorScale=",1);
    signBpol = getNumber(buf,"BpolSign=",0);
    signBtor = getNumber(buf,"BtorSign=",0);
    signQ    = getNumber(buf,"qSign=",0);
    //COCOS    = getNumber(buf,"COCOS",0);
    psiOver2pi = getBool(buf,"psi/2pi=",0);      // -1 means no, i.e psi is not devided over twopi, this means the same as COCOS>=11
    //if(psiOver2pi==0&&COCOS>=10) psiOver2pi=-1;   // if("psi/2pi=" is not found and COCOS>=10) then take COCOS convention that psi is not devided over 2pi 
    if(psiOver2pi==0) psiOver2pi=1;              // otherwise it is assumed that psirz is psi/2pi

    if(scalBpol==0) scalBpol=1;
    if(scalBtor==0) scalBtor=1;
    if(scale_Bpol!=0) scalBpol=scale_Bpol;       // scale_Bpol from function parameters has priority
    if(scale_Btor!=0) scalBtor=scale_Btor;       // scale_Btor from function parameters has priority
    if(sign_Bpol!=0) signBpol=sign_Bpol;         // sign_Bpol from function parameters has priority
    if(sign_Btor!=0) signBtor=sign_Btor;         // sign_Btor from function parameters has priority
    if(sign_Q!=0)    signQ=sign_Q;               // sign_Q from function parameters has priority
    if(psiOverTwopi!=0) psiOver2pi=psiOverTwopi; // psiOverTwopi from function parameters has priority

    scalBpol = abs(scalBpol);
    scalBtor = abs(scalBtor);

    if(signBpol!=0) signBpol = signBpol>0?1:-1;
    if(signBtor!=0) signBtor = signBtor>0?1:-1;
    if(signQ!=0) signQ = signQ>0?1:-1;
    if(psiOver2pi!=0) psiOver2pi>0?1:-1;

    if(fileFmt==G_EQDSK||fileFmt==D_EQDSK) {     // read(5,'(1a48,3i4)') data,idum,nR,nZ
      std::istringstream str(buf.c_str()+48); 
      str >>i>>nR>>nZ;           
    }
    else if(fileFmt==P_EQDSK||fileFmt==P_EQDSKold) {  // TODO: check what is the position of i ,nR, nZ
      std::istringstream str(buf.c_str());
      str>>name>>date>>i>>nR>>nZ;
    }
    else if(fileFmt==J_EQDSK) {  // TODO: check what is the position of i ,nR, nZ
      std::istringstream str(buf.c_str());
      str>>nR>>nZ;
    }
    else return(meshOK=false);

    if(fileFmt==G_EQDSK||fileFmt==J_EQDSK) in.width(16); // see G EQDSK file format(5(e16.9)): five numbers without delimeters between them

    double  dummy,cpasma; // not used

    in>>xdim>>zdim>>Rmajor>>Rmin>>zmid;         //read(5,'(5(e16.9))') rdim,zdim,rcentr,rleft,zmid
    in>>Raxis>>Zaxis>>psiAxis>>psiEdge>>Bcenter; //read(5,'(5(e16.9))') rmaxis,zmaxis,psiAxis,psiEdge,bcentr
    //in>>cpasma>>psiAxis>>dummy>>Raxis>>dummy;   //read(5,'(5(e16.9))') cpasma,psiAxis,dummy,rmaxis,dummy
    //in>>Zaxis>>dummy>>psiEdge>>dummy>>dummy;    //read(5,'(5(e16.9))') zmaxis,dummy,psiEdge,dummy,dummy
    in>>cpasma>>dummy>>dummy>>dummy>>dummy;   //read(5,'(5(e16.9))') cpasma,psiAxis,dummy,rmaxis,dummy
    in>>dummy>>dummy>>dummy>>dummy>>dummy;    //read(5,'(5(e16.9))') zmaxis,dummy,psiEdge,dummy,dummy

    if(!resize()) return meshOK=false;

    for(i=0;i<nR;i++) in>>Fp_[i];       //read(5,'(5(e16.9))') (fpol(i),i=1,nR) fpol is defined on psi_ and fpol==R*Btor
    for(i=0;i<nR;i++) in>>dummy;        //read(5,'(5(e16.9))') (pressure(i),i=1,nR)
    for(i=0;i<nR;i++) in>>FFp_[i];      //read(5,'(5(e16.9))') (ffprim(i),i=1,nR)
    for(i=0;i<nR;i++) in>>Pp_[i];       //read(5,'(5(e16.9))') (pprime(i),i=1,nR)
    for(k=0;k<nZ;k++)
      for(i=0;i<nR;i++) in>>psirz(i,k); //read(5,'(5(e16.9))') ((psirz(i,j),i=1,nR),j=1,nZ)
    for(i=0;i<nR;i++) in>>q_[i];        //read(5,'(5(e16.9))') (qpsi(i),i=1,nR)
    //***End reading
    watch.printTime("CEfit::loadx--end reading ",true);
    // create rectangular mesh
    dR   = xdim/(nR-1);
    dZ   = zdim/(nZ-1);
    edR = 1/dR;
    edZ = 1/dZ;
    for(i=0;i< nR; i++) R[i] = Rmin+i*dR; Rmax=R[nR-1];
    for(k=0;k< nZ; k++) Z[k] = zmid - zdim/2.0 + k*dZ; Zmin = Z[0]; Zmax=Z[nZ-1];
#if 1
    // scale poloidal and toroidal fields
    for(i=0;i<nR;i++) { 
      Fp_[i]  *= scalBtor;
      //FFp_[i] /= scalBpol;  // FFp_ is not used  
      //Pp_[i]  /= scalBpol;  // Pp_ is not used 
       q_[i]  *= scalBtor/scalBpol; // q ~ r^2/R^2 * Ipol/Itor ~rB*tor/R*Bpol
    }
    if(psiOver2pi==-1) scalBpol /= twopi; 
    psiEdge *= scalBpol;
    psiAxis *= scalBpol;
    if(fileFmt!=J_EQDSK||!(fileFmt==P_EQDSK&&psiSign==1&&psiAxis>psiEdge&&name=="ITER-FEAT")) { // strange format; its looks like 0<=psirz(i,k)<=1
      for(k=0;k<nZ;k++)    // (0 on axis and 1 at the separatrix)
        for(i=0;i<nR;i++) 
          psirz(i,k) *= scalBpol;
    }
#endif
    Bcenter = Fp_[0]/Raxis;
    psiSign = getPsirzSign();
    double qEdge = q_[nR-1],R1=Raxis,R2=Rmin; // Raxis - (Raxis-Rmin)/3; 
    R2=(Rmin+Raxis)/2;
    qEdge = q_[nR/2];
    psi_.clear(); // disable psi_
    // find out the original signes of Bpol, Btor, and q stored in the file; see also scalBpol, scalBtor scaling above 
    double psi_edge, r_edge;  // edge valus from psirz
    Vector3d BhfsEdge = findPsiEdge(qEdge, psi_edge, r_edge, R1, R2); // BhfsEdge is the B at high field side in (x,z) plane 
    double BpolSignOrig = (BhfsEdge[2]>0)?1:-1;     // what is the sign of B derived from EQDSK file ?
    double BtorSignOrig = (BhfsEdge[1]>0)?1:-1;
    double qSignOrig = (q_[nR/3]>0)?1:-1;

    // change poloidal or toroidal fields direction if needed
    if(signBtor!=0) for(int i=0;i<nR;i++) Fp_[i] *= (BtorSignOrig*signBtor); 
    if(signQ!=0)    for(int i=0;i<nR;i++) q_[i] *= (qSignOrig*signQ); 

    Bcenter = Fp_[0]/Raxis;
    if(signBpol!=0) {
      psiEdge *= (BpolSignOrig*signBpol);
      psiAxis *= (BpolSignOrig*signBpol);
      for(k=0;k<nZ;k++)
        for(i=0;i<nR;i++) 
          psirz(i,k) *= (BpolSignOrig*signBpol);
    }
    
    deltaPsi = psiEdge - psiAxis;

#if DENORMALIZE_JET_PSIRZ
    {  // denormalize psirz in from some file. TODO: must be revised or returned with error.
      if(fileFmt==P_EQDSK&&psiSign==1&&psiAxis>psiEdge&&name=="ITER-FEAT") {
        fileFmt=P_EQDSKold;
        for(k=0;k<nZ;k++) // denormalize this strange format; its looks like 0<=psirz(i,k)<=1
          for(i=0;i<nR;i++) psirz(i,k) = psirz(i,k)*deltaPsi + psiAxis;
      }
      if(fileFmt==J_EQDSK) {
        for(k=0;k<nZ;k++) //// denormalize this strange format JET format; its looks like 0<=psirz(i,k)<=1 
          for(i=0;i<nR;i++) psirz(i,k) = psirz(i,k)*deltaPsi + psiAxis;
      }
    }
#endif

    watch.printTime("CEfit::loadx--end signes ",true);
    psiSign = getPsirzSign();  // renew psiSign
    if(!(fileFmt==C_TOKAMAK||fileFmt==C_TOKAMAK_EQUIVALENT)) {  // find magnetic axis from psirz
      amoebaFunc func(this,psiSign);
      int ndim(2);
      double  y[3],v[6],*p[3];
	    for (int i=0;i<ndim+1;i++) p[i] = v+i*2;
      double R0 = (Rmin+Rmax)/2;
      double Z0 = (Zmin+Zmax)/2;
      p[0][0] = R0;    p[0][1] = Z0;
      p[1][0] = R0+dR; p[1][1] = Z0;
      p[2][0] = R0;    p[2][1] = Z0+dZ;
      y[0] = func(p[0]);
      y[1] = func(p[1]);
      y[2] = func(p[2]);
      amoeba(p,y,ndim,1e-7,func);
      if(func.ncalls()>1000) std::cerr<<"CStconfig::amoeba().magcrdBmin: #func="<<func.ncalls()<<std::endl;
      Vector3d mAxis = Vector3d(p[2][0],0,p[2][1]);  // magnetic axis in cyl. coord.
      Raxis = mAxis[0];
      Zaxis = mAxis[2];
      double psiAxis1 = psiAxis;
      psiAxis = cyl2psi(mAxis);
      double psiAxis2 = psiAxis;
    }
    watch.printTime("CEfit::loadx--find magnetic axis ",true);
    // end reading of EQDSK file 
  }  
// psi_[], the uniform poloidal flux grid is not defined here, it will be created later


  psi_.clear(); // disable psi_
  if(!(fileFmt==C_TOKAMAK||fileFmt==C_TOKAMAK_EQUIVALENT)) {
#if FIND_PSI_EDGE
    double R1=Raxis,R2=Rmin;
    double r_edge, qEdge = q_[nR-1];
    findPsiEdge(qEdge, psiEdge, r_edge,R1, R2); // find psiEdge in order to create array psi_  
    watch.printTime("CEfit::loadx--findPsiEdge ",true);
#else
    double psiEdgeOriginal = psiEdge;
    psiEdge = psiAxis + deltaPsi;
#endif
  }

  psi_.init(nR);
  dpsi = (psiEdge-psiAxis)/(nR-1);
  for(i=0;i<nR;i++) psi_[i]=psiAxis+i*dpsi;   //initial(EFIT's) uniform poloidal flux grid

  {// Find Redge by bisection
    Vector3d rOutside = Vector3d(Rmax, 0,Zaxis); // outside point
    //Vector3d rOutside = Vector3d(Rmin, 0,Zaxis); // outside point ; Rmin as outside point also works fine
    Vector3d rInside  = Vector3d(Raxis,0,Zaxis); // inside point
    for(int itr=200; itr; --itr) {
      Vector3d r = (rInside+rOutside)/2;
      double dPsi = (psiEdge - cyl2psi(r))/(psiEdge-psiAxis); //normalized poloidal flux
      if(dPsi< 0.)
        rOutside   = r;
      else {
        rInside = r;      // if r lies inside
        if(dPsi<=1e-9) break;  // Ok, done
      }
    }
    Redge = rInside[0];
  }
    watch.printTime("CEfit::loadx--Find Redge by bisection ");

#if PRINT_ARRAYS  // read LCMS
  if(fileFmt==G_EQDSK) {
    pBase::ngArray<Vector3d> LCMS; ///< last closed magnetic surface
    int nbdry,limitr;
    if(fileFmt==G_EQDSK) in.width(0);
    in>>nbdry>>limitr;                  //read(5,2022)  nbdry,limitr

    if(fileFmt==G_EQDSK) in.width(16);
    LCMS.clear();
    if(nbdry) {
      LCMS.init(nbdry);
      //if(fmt==P_EQDSKold) {
      //  for(i=0;i<nbdry;i++) in>>LCMS[i].x();
      //  for(i=0;i<nbdry;i++) in>>LCMS[i].z();
      //}
      //else
      for(i=0;i<nbdry;i++) in>>LCMS[i].x()>>LCMS[i].z();  //read(5,2020) (rbdry(i),zbdry(i),i=1,nbdry)

     // find edge segment of LCMS
      int iS=0;
      for(i=0;i<nbdry-1;i++) {
        if(LCMS[i].x()>Raxis&&LCMS[i].z()==Zaxis) {iS=i; break;}
        if(LCMS[i].x()>Raxis&&(LCMS[i].z()-Zaxis)*(LCMS[i+1].z()-Zaxis)<0)  {iS=i; break;}
      }
      double Redge = LCMS[iS].x() +(Zaxis - LCMS[iS].z())*(LCMS[iS].x()-LCMS[iS+1].x())/(LCMS[iS].z()-LCMS[iS+1].z());
      FILE *fb = fopen("EFIT_lcms.prn","w");
      if(fb) for(i=0;i<nbdry;i++) fprintf(fb,"%g %g\n", LCMS[i].x(),LCMS[i].z());
      fclose(fb);
    }
    else
      std::cerr<<" CEfit::load() : EQDSK file: plasma boundary not found."<<std::endl;
  }
//// calculate toroidal flux using EFIT's q
//  for(i=1, s_[0]=0; i<nR; i++) s_[i]=s_[i-1]+FluxTor(psi_[i-1],psi_[i]);
//// s-grid from EFIT's q
//  fluxTorMax = s_[nR-1];
//  for(i=0;i< nR; i++) s_[i] /= fluxTorMax;
//  std::ofstream out("EFIT's_spol_s_q.prn");
//  for(i=0;i<nR;i++) out<<(psi_[i]-psiAxis)/(psiEdge-psiAxis)<<"  "<<s_[i]<<"  "<<q_[i]<<std::endl;

    FILE *ff = fopen("psirz.prn","w");
    for(k=0;k<nZ;k++) {
      for(i=0;i<nR;i++) fprintf(ff,"%g ", psirz(i,k));
      fprintf(ff,"\n");
    }
    fclose(ff);
#endif


//TODO: ????? 2021_0916
  ////new_nR=nR;
  ////for(int i=0; i<nR; i++) {  //added 2015_1124
  ////  if(q_[i]<0) { 
  ////    new_nR = i; 
  ////    break;
  ////  }
  ////}
  return meshOK=true;
}

//****************************************************************************
// calculate toroidal flux using EFIT's q
// Integral(psi1,psi2) q(psi)*dpsi
double CEfit::FluxTor(double psi1, double psi2) const
{
  double sum=0;
  int n=20;
  n=n+n%2;         //must be even for simpson integration
  double dp=(psi2-psi1)/n;
  if(dp==0) return 0;
  for(int k=0; k<=n; k++) {  //Simpson integration
    double w = (k==0||k==n)?1:(2*(1+k%2));   // weight for Simpson integration
    sum += w*q(psi1+k*dp);
  }
  return fabs(sum*dp/3);
};

//****************************************************************************
// calculate toroidal flux
// Integral(psi1,psi2) qTraced(psi)*dpsi
// where qTraced has been calculated by field line tracing
double CEfit::mFluxTor(double p1, double p2) const
{
  double sum=0;
  int n=20;
  n=n+n%2;         //must be even for simpson integration
  double dp=(p2-p1)/n;
  if(dp==0) return 0;
  for(int k=0; k<=n; k++) {  //Simpson integration
    double w = (k==0||k==n)?1:(2*(1+k%2));   // weight for Simpson integration
    sum += w*qTraced(p1+k*dp);
  }
  return fabs(sum*dp/3);
};

//***************************************************************
// Calculate coefficients and derivatives for the polynomial interpolation
inline static void LagrangeCoeffs(double &p, double *c, double *cd)
{
  static const double i6 = 1./6.;
  c[0] =-(p-1)*(p-2)*(p-3)/6;
  c[1] = p*(p-2)*(p-3)/2;
  c[2] =-p*(p-1)*(p-3)/2;
  c[3] = p*(p-1)*(p-2)/6;
  cd[0]=-( (p-2)*(p-3) + (p-1)*(p-3) + (p-1)*(p-2) )/6;  // derivative d c[0] / dp
  cd[1]= ( (p-2)*(p-3) +     p*(p-3) +     p*(p-2) )/2;
  cd[2]=-( (p-1)*(p-3) +     p*(p-3) +     p*(p-1) )/2;
  cd[3]= ( (p-1)*(p-2) +     p*(p-2) +     p*(p-1) )/6;
}

//***************************************************************
// function finds the integer coordinates on the mesh
//  and calculates coefficients for interpolation
bool CEfit::FindMeshCoordinates(const Vector3d &cyl)
{
  if(!meshOK) return false;

  if( cyl[0]<Rmin||cyl[0]>Rmax||cyl[2]<Zmin||cyl[2]>Zmax ) return false; // return if outside

  double pr = (cyl[0]-Rmin)*edR;
  double pz = (cyl[2]-Zmin)*edZ;

// set the 1st point for the Lagrange interpolation
  j0 = int(pr) -1; // int R
  k0 = int(pz) -1; // int Z

  j0=mmax(j0,psirz.lowbound1());
  j0=mmin(j0,psirz.upperbound1()-3);
  k0=mmax(k0,psirz.lowbound2());
  k0=mmin(k0,psirz.upperbound2()-3);

  pr=(cyl[0]-R[j0])*edR;
  pz=(cyl[2]-Z[k0])*edZ;
  LagrangeCoeffs(pr, cr, crd );
  LagrangeCoeffs(pz, cz, czd );

  if(!doPositionCheck) return true;

  return true;
}

#define FIND_COORD   (const_cast<CEfit*>(this))->FindMeshCoordinates

//***************************************************************
// Find flux label for the point with cylindrical coordinates 'cyl'
//  function returns  s - normalized tor. flux
double CEfit::cyl2s(const Vector3d &cyl) const
{
  double s(0);
  if(!FIND_COORD(cyl)) return 1001;
  s = CEfit::sNew(cyl2psi(cyl));
  if(s<0) s=0;
  return (s<=sMax)?s:1001;
}

//***************************************************************
// Find poloidal flux at the point with cylindrical coordinates 'cyl'
//  function returns  psi - poloidal flux(which was tabulated on R-Z grid)
double CEfit::cyl2psi(const Vector3d &cyl) const
{
  double psi(0);
  if(!FIND_COORD(cyl)) return 1001;

  int j,k;
  for(j=0;j<4;j++) { // loop over cube (j0,k0)..(j0+4,k0+4)
    const double * pflux = psirz[j+j0];
    for(k=0;k<4;k++)
      psi  += pflux[k+k0]*(cr[j]*cz[k]);
  }
  return psi;
}

//***************************************************************
// cylindrical coordinates are used
//   function retuns B in  cylindrical coordinates
// return 0 if outside Rmin,Rmax,Zmin,Zmax
Vector3d CEfit::getBcyl(const Vector3d &cyl, double *Psi) const
{
  Vector3d b(0.,0.,0);
  if(!FIND_COORD(cyl)) return b; // return zero if point lies outside
  double psi=0;
  int j,k;
  for(j=0;j<4;j++) {
    const double * pflux = psirz[j+j0];
    for(k=0;k<4;k++) {
      psi  += pflux[k+k0]*(cr [j]*cz[k]);
      b[0] -= pflux[k+k0]*(cr [j]*czd[k]);  // -dpsi/d(z/dZ)
      b[2] += pflux[k+k0]*(crd[j]*cz [k]);  //  dpsi/d(r/dR)
    }
  }
  double R = cyl[0];
  b[0] *= edZ;
  b[1] = psi_.empty()?Fp_[0]:RBtor(psi);  //Fp_[0] ~ Bcenter*R
  b[2] *= edR;
  b /= R; // -dpsi/dz/R, RBtor/R, dpsi/dr/R
  if(Psi) *Psi = psi;
  return b;
}

//****************************************************************************
// Quadratic interpolation in psi using initial(EFIT's) uniform poloidal flux grid
// x and y are input arrays
// on return y(u)
double CEfit::interp2psi_(double u, const pBase::ngArray<double> &y) const 
{
  const pBase::ngArray<double> &x = psi_; // initial mesh
  int N = (int)x.size();
  ////N = new_nR;
  int L = int((u-x[0])/dpsi);
  if(L+2>=N) L=N-3;
  else if(L<0) L=0;

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
// perform a binary search of the segment in which u lies. x input array
// ! size <= x.size()
int CEfit::bsearch(double u, const pBase::ngArray<double> &xa,int size) const
{
  int k;
  int Lb = 0;    // Left bracket
  int Rb = size-1;  // Right bracket
  const double * x = xa.constArray();
  int sign = x[Rb]>x[Lb]?1:-1;
  if(sign*(u-x[0])<=0) return Lb;
  else if(sign*(u-x[Rb])>=0) Lb=Rb-1;
  else
    while(Lb+1 < Rb) {  // bisection
      k = (Lb + Rb)/2;
      if(sign*(u-x[k])<0) Rb = k;
      else Lb = k;
    }
  return Lb;
}

//****************************************************************************
// Quadratic interpolation in u. x and y are input arrays
// on return y(u)
// ! size must be < x.size()
double CEfit::interp2(double u, const pBase::ngArray<double> &y, const pBase::ngArray<double> &x,int size,double *yp)  const
{
  int N = size;
  int L = bsearch(u,x,size);
  if(L+2>=N) L=N-3;
  else if(L<0) L=0;

  double y0 = y[L  ];
  double y1 = y[L+1];
  double y2 = y[L+2];
  double x0 = x[L  ];
  double x1 = x[L+1];
  double x2 = x[L+2];
  double c1=(y1-y0)/(x1-x0);
  double c2=(y2-y1)/(x2-x1);
  c2=(c2-c1)/(x2-x0);
  if(yp)
    *yp = c1 + c2*(u-x1)+c2*(u-x0); // derivative
  return  y0 + (c1 + c2*(u-x1))*(u-x0);
}

inline double CEfit::traceStep(double q)
{
  double iota = q==0?100:(1/fabs(q));
  double factor = mmax(iota, 1e-2);
  if(factor>1) factor=1;
  return RKSTEP*factor;
}

//****************************************************************************
// calculate q by B-field line tracing.
// after full poloidal turn the q=increase_of_torAngle/2pi
// @param cyl is the cylindrical coordinates
double CEfit::traceForQ(const Vector3d &cyl, double q_guess)
{
  if(q_guess==0) q_guess = 10; 
  double step = traceStep(q_guess);
  Bfield Bfld(*this);
  btrace.init(Bfld,step);

  double dToroidalAngle = RKSTEP;
  Vector3d c(cyl);
  //double fi0 = c[1];
  double z0  = c[2];
  double qFactor = 2e30;  // if q>= 1e30 or the caller routine is signaled that the field-line goes outside plasma

  double t = 0;           // poloidal angle
  double fi0 = t*q_guess;
  double dtheta = 1/1024.; // 1/128 ~ 0.44*degree is the increment for poloidal angle

  int i=0;
  int Ntr = 10000*int(twopi/dToroidalAngle);
  btrace.init();
  for(int j=0; j<Ntr; j++) {
    Vector3d cprev = c;
    double fiprev = t*q_guess;   // toroidal angle
    btrace.advancePoloidally(dtheta,q_guess,c);
    t += dtheta;
    if(c[0]<Rmin||c[0]>Rmax||c[2]<Zmin||c[2]>Zmax) break; // return 0 if outside
    if(j>0&&(cprev.z()-z0)*(c.z()-z0)<=0) {
      i++;
      if(cprev.z()==z0) i--;
    }
    if(i==2) {   // full poloidal turn if i==2
      double fi = t*q_guess;   // fi is the toroidal angle
      fi = fiprev +(z0-cprev.z())*(fiprev-fi)/(cprev.z()-c.z());
      fi -= fi0;
      qFactor = fi/twopi;
      break;
    }
  }

  return qFactor;
}

//****************************************************************************
// B field line tracing
// RBZt_[] is filled
double CEfit::traceForSurface(int is, double &dVds, double &Itor)
{
  double q = mq_[is];
  double step = traceStep(q);
  Bfield Bfld(*this);
  btrace.init(Bfld,step);

  double dPoloidalAngle = twopi/Ntrace;

  Vector3d c(mR_[is],0,Zaxis);
  double B = getBcyl(c).abs();

  RBZt_[0]    = c;
  RBZt_[0][1] = B;                  // Vector3d RBZt_[0] holds (R,B,Z)
  dVds = 0;
  Itor = 0;

  btrace.init();
  for(int i=1; i<Ntrace; i++) {
    btrace.advancePoloidally(dPoloidalAngle,q,c);
    Vector3d b = getBcyl(c);
    B = b.abs();
    RBZt_[i]    = c;
    RBZt_[i][1] = B;                // Vector3d RBZt_[i] holds (R,B,Z)
    double jacobian = fluxTorMax*c[0]/b[1];  //Jacobian in (s,theta,phi) fluxTorMax*R/Btor ; fluxTorMax is flux_inWeber/twopi
    dVds += jacobian;
 // Itor = 1/mu0*Integral(dtheta*B_theta); where B_theta = B*e_theta, e_theta = dX/dtheta; e_theta is the covariant basis vector
    Vector3d e_theta = RBZt_[i]-RBZt_[i-1];
    double Bcova = b*e_theta;
    Itor += Bcova;
  }

  dVds *= dPoloidalAngle*twopi;
  const double mu0 = 4e-7*pi; // [N/A^2 = henry/m]
  Itor /= mu0;
  return dPoloidalAngle;
}

static void gsmooth(Vector3d *u,const Vector3d *y,int N,int b);

//****************************************************************************
void CEfit::createSurface(int is, double &dVds, double &Itor, std::ofstream *out)
{
  if(is>=Ns) is=Ns-1;

  double dPoloidalAngle = traceForSurface(is,dVds,Itor);

//  gsmooth(RBZt_smooth.array(),RBZt_.constArray(),RBZt_.size(),20);
// see also in CEfit::set() //  RBZt_smooth.init(Ntrace);

//TODO: ? poloidal_angle_direction:
//  direction is not needed; tested using TCV, ITER configurations
  Vector3d r1(RBZt_[0]); r1[1]=0;
  Vector3d r2(RBZt_[4]); r2[1]=0;
  Vector3d r3 = r1^r2;
  int sign = r3[1]>0?1:-1;
//Discrete Fourier Transform
  for(int m=0; m<=M; m++) {
    RBZm_[m]=0;
    RBZm1_[m]=0;
    double argm = m*dPoloidalAngle;
    double cm = ::cos(argm);
    double sm = ::sin(argm);
    double cs1, cs(1), sn(0);
    for(int i=0; i<Ntrace; i++) {
      RBZm_[m][0] += RBZt_[i][0]*cs;     
      RBZm_[m][1] += RBZt_[i][1]*cs;    
      RBZm_[m][2] += RBZt_[i][2]*sn;    

      RBZm1_[m][0] += RBZt_[i][0]*sn;
      RBZm1_[m][1] += RBZt_[i][1]*sn;
      RBZm1_[m][2] -= RBZt_[i][2]*cs;
      cs1= cs*cm - sn*sm;
      sn = sn*cm + cs*sm;
      cs = cs1;
    }
    RBZm_[m]  *= 2./Ntrace;
    RBZm1_[m] *= 2./Ntrace;
  }
  RBZm_[0]  /= 2;
  RBZm1_[0] /= 2;

  if(out) for(int i=0; i<Ntrace; i++) *out<<RBZt_[i]<<std::endl;
}


//************************************************************
// EXPm(x) == exp(-x)
static double EXPm(double x)
{
  const double WGauss=20;
  static int max = 512;
  static double *dat=NULL;
  static double dx,edx;
  int i;
  if(dat==NULL) {
    dat = new double[max+1];
    dx = WGauss/max; edx = max/WGauss;
    for(i=0;i<=max;i++) dat[i] = ::exp(-i*dx);
  }
  double p = x*edx;
  i = int(p);
  if(i>=max) return(dat[max]);
  p = x*edx - i;
  return( (1.-p)*dat[i] + p*dat[i+1] );
}

//****************************************************************************
// Returns 'u' N-element vector created by smoothing
// using a Gaussian kernel to return weighted averages of the elements in y,
//  the elements of x must be in ascending order
//  b is the bandwidth of the smoothing window, it should be set to
//  a few times the spacing between x data points, try b=5
static void gsmooth(Vector3d *u,const Vector3d *y,int N,int b)
{
  int i,j;
  double twopi = 2*3.14159265358979323846;
  double c=3.6523/(b*b);  // c=1/(b*b*0.2738);
  for (i=0; i<N; i++) {
    Vector3d s1=y[i];
    double s2=1;
    for (j=i+1; ; j++) {
      double a = c*(j-i)*(j-i);
      if(a>15) break;
      double K = EXPm(a);  // K = exp(-a);
      s1 += K*y[j%N];
      s2 += K;
    }
    for (j=i-1; ; j--) {
      double a = c*(j-i)*(j-i);
      if(a>15) break;
      double K = EXPm(a);  // K = exp(-a);
      int j1=(j<0)?(j+N):j;
      s1 += K*y[j1];
      s2 += K;
    }
    u[i] = s1/s2;
  }
}

//****************************************************************************
//***Find minimum*************************************************************
//***Downhill Simplex Method in Multidimensions ******************************
//****************************************************************************
#define NMAX  19000
#define ALPHA 1.0
#define BETA  0.5
#define GAMMA 2.0
#define NDIM  10

#define GET_PSUM for(j=0;j<ndim;j++)  \
                   for(i=0,psum[j]=0;i<mpts;i++) psum[j] += p[i][j];

void CEfit::amoeba(double **p,double *y,int ndim,double ftol,amoebaFunc &funk)
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
    if(funk.ncalls() >= NMAX) std::cerr<<"CStconfig::amoeba(): too many iterations"<<std::endl;
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

double CEfit::amotry(double **p,double *y,double *psum,int ndim,int ihi,double fac,double *ptry,amoebaFunc &funk)
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

};   //namespace MConf
