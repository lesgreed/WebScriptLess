#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "../include/C3dMesh.h"
#include "../include/CPfunction.h"

#include <iostream>
#include <sstream>
#include "mconf_matlab.h"

using MConf::C3dMesh;
using MConf::CStconfig;
using MConf::Vector3d;
using MConf::CProfile;
using MConf::CPfunction;

#undef True
#undef False

const static int True(1);
const static int False(0);

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
//******************************************************************************
// Checks to see if a directory exists. Note: This method only checks the
// existence of the full path AND if path leaf is a dir.
// 
// @return  >0 if dir exists AND is a dir,
//           0 if dir does not exist OR exists but not a dir,
//          <0 if an error occurred (errno is also set)
static int isDir(const char* path)
{
  struct stat info;
  int statRC = stat( path, &info );
  if(statRC != 0 ){
    if(errno == ENOENT)  { return 0; } // something along the path does not exist
    if(errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
    return -1;
  }
  int k = (info.st_mode&S_IFDIR)?1:0;
  return k;
}

//******************************************************************************
// Create temp folder for files
// @return path to the temp dir with a trailing slash
static std::string getTempDir()
{
#ifdef WIN32
  char *key[] = {"TMP", "TEMP", "USERPROFILE", NULL};
#else
  char *key[] = {"TMPDIR", "TMP", "TEMP", "TEMPDIR","USERPROFILE", NULL};
#endif
  char *val;
  for(int i=0; key[i]!=NULL; i++) {
    val = getenv(key[i]);
    if(val) break;
  }
#ifdef WIN32
  std::string path = (val==NULL)?std::string(""):std::string(val);
  path += "\\MConf";
  std::string mkdir("mkdir "); 
  mkdir += path;
  int k = isDir(path.c_str());
  path += "\\";
#else
  std::string path = (val==NULL)?std::string("/tmp"):std::string(val);
  path += "/MConf";
  std::string mkdir("mkdir -p "); 
  mkdir += path;
  int k = isDir(path.c_str());
  path += "/";
#endif
  if(k==0) system(mkdir.c_str());
  return path;
}

// create unique name for temp file
static std::string createFileName(const std::string &name)
{
  std::string tmpname = name;
  for(int n=1;;++n) {
    FILE* fp = ::fopen(tmpname.c_str(),"r");   // try to open
    if(fp==0) break; // break if not found
    ::fclose(fp);
    std::ostringstream os;
    os<<name<<"_"<<n;
    tmpname=os.str();
  }
  return tmpname;
}



#ifdef __cplusplus
extern "C" {
#endif

//*********************************************************************
// constructor
static MC_HANDLE mc_load(const char * fname, CStconfig::fileFormat fmt)
{
  C3dMesh * mc = new C3dMesh;  
  MC_HANDLE mConf = (MC_HANDLE)mc;
  bool ok = mc->load (fname,fmt);
  if(!ok) {
    std::cout << "MCLOAD: Loading error, file:'"<<fname<<"'\n";
    delete (mc);
    mConf = 0;    
  }
  return mConf;
}

WINDLLEXPORT MC_HANDLE MCload(char * fname)
{
  return mc_load(fname, CStconfig::UNDEFINED);
}

WINDLLEXPORT MC_HANDLE MCloadVMEC(char * fname)
{
  return mc_load(fname, CStconfig::VMEC);
}

WINDLLEXPORT MC_HANDLE MCloadLHD(char * fname)
{
  return mc_load(fname, CStconfig::LHD);
}

WINDLLEXPORT MC_HANDLE MCloadEQDSK(char * fname, double scaleBpol, double scaleBtor, int signBpol,int signBtor,int signQ,int psiOverTwopi)
{
  C3dMesh * mc = new C3dMesh;  
  MC_HANDLE mConf = (MC_HANDLE)mc;
  bool ok = mc->loadEfitAndScale(fname, scaleBpol, scaleBtor,signBpol,signBtor,signQ,psiOverTwopi);
  if(!ok) {
    std::cout << "MCLOAD: Loading error, file:'"<<fname<<"'\n";
    delete (mc);
    mConf = 0;    
  }
  return mConf;
}

//*********************************************************************
// copy constructor
WINDLLEXPORT MC_HANDLE MCcopy(MC_HANDLE mConfSrc)
{
  C3dMesh * dest = new C3dMesh;
  C3dMesh &src  = *((C3dMesh*)mConfSrc);
  if(dest) *dest = src;
  return (MC_HANDLE)dest;
}

//*********************************************************************
// destructor
WINDLLEXPORT void MCfree(MC_HANDLE mConf)
{
  if (mConf) delete ((C3dMesh*)mConf);
}

WINDLLEXPORT double MCgetB00(MC_HANDLE mConf)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.getB00();
}

WINDLLEXPORT void MCsetB00(MC_HANDLE mConf, double B00)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.setB00(B00);
}

WINDLLEXPORT void MCsetsmax(MC_HANDLE mConf, double smax)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.setsMax(smax);
}

WINDLLEXPORT double MCgetB0(MC_HANDLE mConf,double ficyl)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.getB0(ficyl);
}

WINDLLEXPORT void MCsetB0(MC_HANDLE mConf,double B0, double ficyl)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.setB0(B0,ficyl);
}

WINDLLEXPORT void MCsetIpLCMS(MC_HANDLE mConf,double Ip)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.setIcoil(Ip);
}


WINDLLEXPORT void MCtruncate(MC_HANDLE mConf, double epsTrunc)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.truncate(epsTrunc);
}

WINDLLEXPORT void MCsetAccuracy(MC_HANDLE mConf, double epsA)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.setAccuracy(epsA);
}

WINDLLEXPORT int MCgetRayIntersectionPoints(MC_HANDLE mConf,const double *r0,const double *rd,double *entryPoint,double *exitPoint)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R0(r0),Rd(rd);
  Vector3d entry;
  Vector3d exit;
  bool b = mc.getRayIntersectionPoints(R0,Rd,entry,exit);
  entry.copyTo(entryPoint);
  exit.copyTo(exitPoint);
  return b?True:False;
}

WINDLLEXPORT double MCxyz2s(MC_HANDLE mConf,const double *xyz)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(xyz);
  return mc.xyz2s(R);
}

WINDLLEXPORT void MCxyz2mag(MC_HANDLE mConf,const double *xyz,double *mag)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(xyz);
  Vector3d mg = mc.xyz2mag(R);
  mg.copyTo(mag);
}

WINDLLEXPORT double MCgetBxyz(MC_HANDLE mConf,const double *xyz,double *B)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(xyz);

  double s = mc.xyz2s(R);
  Vector3d b = mc.getBxyz(R);
  b.copyTo(B);
  return s;
}

WINDLLEXPORT double MCgetGradBxyz(MC_HANDLE mConf,const double *xyz,double *gradB)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(xyz);

  double s = mc.xyz2s(R);
  Vector3d gb = mc.getGradBxyz(R);
  gb.copyTo(gradB);
  return s;
}


WINDLLEXPORT double MCreff (MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.r(s);
}

WINDLLEXPORT double MCR0(MC_HANDLE mConf)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.R0();
}

WINDLLEXPORT double MCtorFlux2polFlux(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.sPol(s);
}

WINDLLEXPORT double MCPoloidalFlux(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.PoloidalFlux(s);
}

WINDLLEXPORT double MCsToroidal(MC_HANDLE mConf,double spol)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.sToroidal(spol);
}

WINDLLEXPORT void MCmag2xyz(MC_HANDLE mConf,const double *magCoord, double *xyz)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d mag(magCoord);
  Vector3d c = mc.mag2xyz(mag);
  c.copyTo(xyz);            
}

WINDLLEXPORT void MCmix2xyz(MC_HANDLE mConf,const double *mixCoord, double *xyz)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d mix(mixCoord);
  Vector3d c = mc.mixcoord2xyz(mix);
  c.copyTo(xyz);            
}

WINDLLEXPORT double MCVprime(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.Vprime(s);
}

WINDLLEXPORT double MCVolume(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.Volume(s);
}

WINDLLEXPORT double MCiota(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.iota(s);
}

WINDLLEXPORT double MCpressure(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.pressure(s);
}

WINDLLEXPORT double MCiotaPrime(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.iotaPrime(s);
}

WINDLLEXPORT double MCIp(MC_HANDLE mConf,double s) 
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.Ip(s);
}

WINDLLEXPORT double MCIt(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.It(s);
}

//WINDLLEXPORT double MCBmax(MC_HANDLE mConf,double s)
//{
//  C3dMesh &mc = *((C3dMesh*)mConf);
//  return mc.Bmax(s);
//}

WINDLLEXPORT double MCFlux(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.Flux(s);
}

WINDLLEXPORT double MCgradsAvrg(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.GradsAvrg(s);
}


WINDLLEXPORT double MCgrads2Avrg(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.Grads2Avrg(s);
}

WINDLLEXPORT double MCgradRhoAvrg(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  double r   = mc.r(s);
  double a   = mc.r(1);
  double dr_ds  = r/(2*s); // dr_vmec / ds 
  return mc.GradsAvrg(s)*dr_ds/a;
}


WINDLLEXPORT double MCgradRho2Avrg(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  double r   = mc.r(s);
  double a   = mc.r(1);
  double dr_ds  = r/(2*s); // dr_vmec / ds 
  return mc.Grads2Avrg(s)*(dr_ds*dr_ds)/(a*a);
}

WINDLLEXPORT double MCftrapped(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.ftrapped(s);
}

WINDLLEXPORT double MCFbs(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.FbsVmec(s)*2*s/mc.r(s);
}

WINDLLEXPORT double MCFbsVmec(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.FbsVmec(s);
}

WINDLLEXPORT double MCFbsg2(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.Fbsg2(s);
}

WINDLLEXPORT double MCFbsu(MC_HANDLE mConf,double s){
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.Fbsu(s);
}
 
WINDLLEXPORT double MCFbsuB2(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.FbsuB2(s);
}
 
WINDLLEXPORT double MCFbsu2B2(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.Fbsu2B2(s);
}
 
WINDLLEXPORT double MCFbsB2(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.FbsB2(s);
}
 

WINDLLEXPORT double MCFbsg4(MC_HANDLE mConf,double s, double mu)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.Fbsg4(s,mu);
}

WINDLLEXPORT void MCFbsSetSlabelParam(MC_HANDLE mConf,double smin, double smax, int size)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.FbsSetSlabelParam(smin, smax, size);
}

WINDLLEXPORT void MCFbsSetXeffParam(MC_HANDLE mConf,double xmin, double xmax, int size)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.FbsSetXeffParam(xmin, xmax, size);
}

WINDLLEXPORT void MCFbsSetTracingParam(MC_HANDLE mConf,int turns, bool doAccuracyTest, double dphi)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.FbsSetTracingParam(turns, doAccuracyTest, dphi);
}

WINDLLEXPORT void MCFbsSetIotaParam   (MC_HANDLE mConf,bool avoidResonances, bool useRationalIota, int iotaDenominator)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.FbsSetIotaParam(avoidResonances, useRationalIota, iotaDenominator);
}

WINDLLEXPORT void MCFbsSetMagnMomentParam (MC_HANDLE mConf,int muLevels, double yMax)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.FbsSetMagnMomentParam (muLevels,yMax);
}

//Mesh prefix is M3D

WINDLLEXPORT int MCcreateMesh (MC_HANDLE mConf,double dr,double dz,double dfi)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  bool b = mc.createMesh(dr,dz,dfi);
  return b?True:False;
}

WINDLLEXPORT int MCcreateMeshUsingSymmetry(MC_HANDLE mConf,double dr,double dz,double dfi)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  bool b = mc.createMeshUsingSymmetry(dr,dz,dfi);
  return b?True:False;
}

WINDLLEXPORT int M3DgetRayEntryPoint(MC_HANDLE mConf,const double *r0,const double *rd,double *entryPoint)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R0(r0),Rd(rd);
  Vector3d entry;
  bool b = mc.M3DgetRayEntryPoint(R0,Rd,entry);
  entry.copyTo(entryPoint);
  return b?True:False;
}

WINDLLEXPORT double M3DgetBxyz     (MC_HANDLE mConf,const double *xyz,double *B)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(xyz);

  double s = mc.M3Dxyz2s(R);
  (mc.M3DgetBxyz(R)).copyTo(B);
  return s;
}

WINDLLEXPORT double M3DgetGradBxyz(MC_HANDLE mConf,const double *xyz,double *gradB)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(xyz);

  double s = mc.M3Dxyz2s(R);
  Vector3d gb = mc.M3DgetGradBxyz(R);
  gb.copyTo(gradB);
  return s;
}

WINDLLEXPORT double M3DgetBandGradBxyz(MC_HANDLE mConf,const double *xyz,double *B,double *gradB,double *gradS)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(xyz);
  Vector3d b, gb, gs;

  double s = mc.M3DgetBandGradBxyz(R, b, gb, gs);
  gb.copyTo(gradB);
  gs.copyTo(gradS);
  b.copyTo(B);
  return s;
}

WINDLLEXPORT double M3DgetBandGradB(MC_HANDLE mConf,const double *cyl,double *B,double *gradB,double *gradS)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(cyl);
  Vector3d b, gb, gs;

  double s = mc.M3DgetBandGradB(R, b, gb, gs);
  gb.copyTo(gradB);
  gs.copyTo(gradS);
  b.copyTo(B);
  return s;
}

WINDLLEXPORT double M3Dxyz2s       (MC_HANDLE mConf,const double *xyz)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(xyz);
  return mc.M3Dxyz2s(R);
}

WINDLLEXPORT double M3DgetBGradsxyz(MC_HANDLE mConf,const double *xyz,double *B, double *grads)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(xyz);

  Vector3d b  = mc.M3DgetBxyz(R);
  Vector3d gs = mc.M3DgetGradsxyz(R);

  b.copyTo(B);
  gs.copyTo(grads);
  return mc.M3Dxyz2s(R);
}

WINDLLEXPORT double MCgetdB_Gradsxyz(MC_HANDLE mConf,const double *xyz,double *B,double *dBdx,double *dBdy,double *dBdz,double *gradS)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(xyz);
  Vector3d b, gs;
  Vector3d dbdx,dbdy,dbdz;

  double s = mc.getdBGradsxyz (R,b,dbdx,dbdy,dbdz,gs);
  dbdx.copyTo(dBdx);
  dbdy.copyTo(dBdy);
  dbdz.copyTo(dBdz);
  gs.copyTo(gradS);
  b.copyTo(B);
  return s;
}


WINDLLEXPORT double M3DgetdB_Gradsxyz(MC_HANDLE mConf,const double *xyz,double *B,double *dBdx,double *dBdy,double *dBdz,double *gradS)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(xyz);
  Vector3d b, gs;
  Vector3d dbdx,dbdy,dbdz;

  double s = mc.M3DgetdBGradsxyz(R, b, dbdx,dbdy,dbdz, gs);
  dbdx.copyTo(dBdx);
  dbdy.copyTo(dBdy);
  dbdz.copyTo(dBdz);
  gs.copyTo(gradS);
  b.copyTo(B);
  return s;
}

WINDLLEXPORT double M3DgetdGradsxyz(MC_HANDLE mConf,const double *xyz,double *gradS,double *dGdx,double *dGdy,double *dGdz)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d R(xyz);
  Vector3d dgdx,dgdy,dgdz;

  Vector3d gs =  mc.M3Dgetdgradsxyz(R,dgdx,dgdy,dgdz);
  dgdx.copyTo(dGdx);
  dgdy.copyTo(dGdy);
  dgdz.copyTo(dGdz);
  gs.copyTo(gradS);
  return mc.M3Dxyz2s(R);
}


WINDLLEXPORT double MCepsEff(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.epsEffVmec(s);
}
 

WINDLLEXPORT void MCgetBandGradientsxyz(MC_HANDLE mConf,const double *magCoord, double *B, 
                               double *gradB, double *gradS, double *gradTh, double *gradPh)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d mag(magCoord);
  Vector3d b, gradb, grads, gradth, gradph;

  mc.getBandGradientsxyz(mag, b, gradb, grads, gradth, gradph);
  b.copyTo(B);
  gradb.copyTo(gradB);
  grads.copyTo(gradS);
  gradth.copyTo(gradTh); 
  gradph.copyTo(gradPh);
}

WINDLLEXPORT void  MCgetBandBasisVectorsxyz(MC_HANDLE mConf, const double *u, double *xyz, double *Bxyz, double *gradBxyz, 
                           double *e1con, double *e2con, double *e3con,
                           double *e1cov, double *e2cov, double *e3cov)

{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d mag(u);
  Vector3d b, gradb, xy;
  Vector3d e1_con, e2_con, e3_con;
  Vector3d e1_cov, e2_cov, e3_cov;

  mc.getBandBasisVectorsxyz(mag, xy, b, gradb,e1_con, e2_con, e3_con, e1_cov, e2_cov, e3_cov);
  b.copyTo(Bxyz);
  xy.copyTo(xyz);
  gradb.copyTo(gradBxyz);
  e1_con.copyTo(e1con);
  e2_con.copyTo(e2con);
  e3_con.copyTo(e3con);
  e1_cov.copyTo(e1cov);
  e2_cov.copyTo(e2cov);
  e3_cov.copyTo(e3cov);
}


WINDLLEXPORT double MCgetLocalShear(MC_HANDLE mConf,const double *magCoord)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d mag(magCoord);
  Vector3d B, cyl, grads;
  double localshear = mc.getLocalShear(mag,B,cyl,grads);
  return localshear;
}

WINDLLEXPORT int MCsetLCMS(MC_HANDLE mConf, double newLastSurface) //, const char * fname)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  int nsmax=0;
  double newVolume=0;
  bool varNumberMods=true;
  bool magAxis=false;
  bool circular=false;
  bool append=false;
 // bool ok = mc.writeasciiReduced(fname,nsmax,newVolume,newLastSurface,varNumberMods,magAxis,circular,append);

  std::string path = getTempDir();
  std::string name = "mconfig";
  path += name;
  std::string pname = createFileName(path);
  bool ok = mc.writeasciiReduced(pname.c_str(),nsmax,newVolume,newLastSurface,varNumberMods,magAxis,circular,append);
  bool ok2 = mc.load (pname.c_str(), CStconfig::UNDEFINED);

  ::remove(pname.c_str());

  return ok&ok2?1:0;
}


//*********************************************************************
// save magnetic coordinate file
// *.bin4; *.bin8 are supported
WINDLLEXPORT int MCwrite(MC_HANDLE mConf, const char * fname)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  bool ok= mc.write(fname);
  return ok?1:0;
}


//*********************************************************************
WINDLLEXPORT void MCgetCoeffForAstraCode(MC_HANDLE mConf, double sqrts, double *reff, double *gradr2Avr, double *J, double *G2, double *hVprime, double *B0, double *R0, double *h)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
 double reff_, gradr2Avr_, J_, G2_, hVprime_, B0_, R0_, h_;
 mc.getCoeffForAstraCode(sqrts, reff_, gradr2Avr_, J_, G2_, hVprime_, B0_, R0_, h_);
  *reff = reff_;
  *gradr2Avr = gradr2Avr_;
  *J  = J_; 
  *G2 = G2_; 
  *hVprime = hVprime_; 
  *B0 = B0_;
  *R0 = R0_; 
  *h  = h_;
}

/// Susceptance matrix: S11 in (Flux,theta,phi)-coordinates.
WINDLLEXPORT double MCS11(MC_HANDLE mConf,double s) 
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.S11(s);  
}
/// Susceptance matrix: S12 in (Flux,theta,phi)-coordinates.
WINDLLEXPORT double MCS12(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.S12(s);  
}
/// Susceptance matrix: S21 in (Flux,theta,phi)-coordinates.
WINDLLEXPORT double MCS21(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.S21(s);  
}
/// Susceptance matrix: S22 in (Flux,theta,phi)-coordinates.
WINDLLEXPORT double MCS22(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.S22(s);  
}

WINDLLEXPORT double MCBoozerBmn(MC_HANDLE mConf,double s,int m,int n)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.BoozerBmn(s,m,n);  
}


WINDLLEXPORT double MCelongation(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  double B00 = mc.BoozerBmn(s,0,0);
  double B10 = mc.BoozerBmn(s,1,0);
  double et  = mc.r(s)/mc.R0();
  B10 /= B00;  
  double elongation = pow(et/B10,2.0); 
  return elongation;
}

/// <B^2(s)>
WINDLLEXPORT double MCB2avrg(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.B2avrg(s);  
}
    
/// <B(s)>
WINDLLEXPORT double MCBavrg(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.Bavrg(s);  
}
    
/// Bmin
WINDLLEXPORT double MCBmin(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.Bmin(s);  
}

/// Bmax
WINDLLEXPORT double MCBmax(MC_HANDLE mConf,double s)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.Bmax(s);  
}
WINDLLEXPORT double MCJacobian(MC_HANDLE mConf,const double *magCoord)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  Vector3d mag(magCoord);
  return mc.Jacobian(mag);            
}

WINDLLEXPORT void  MCsetBdirection(MC_HANDLE mConf, int dir)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.setBdirection(dir);
}

WINDLLEXPORT int MCgetBdirection(MC_HANDLE mConf)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.getBdirection();
}

WINDLLEXPORT void MCrestoreBdirection(MC_HANDLE mConf)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.restoreBdirection();
}

WINDLLEXPORT void MCuseMixedProductForJacobian(MC_HANDLE mConf,int yesNo)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.useMixedProductForJacobian(yesNo>0?true:false);
}

#ifdef __cplusplus
}

#endif

