#include "../include/C3dMesh.h"
#include "../include/threads.h"
#include <iostream>
#include <time.h>

namespace MConf {

static const double sLCMS(1.0001);
static const bool NoDerivatives(true);


// undefine this if you don't want debug print
#undef _TESTPRINT

#if (defined( _DEBUG ) || defined(_TESTPRINT) )
 #define  TESTPRINT(x) {(x);}
#else 
 #define  TESTPRINT(x) 
#endif


//****************************************************************************
C3dMesh::C3dMesh()
{
  meshOK=fullMesh=false;
  useSymmetry=false;
  tokamakSymmetry=false;
  Rmin=Rmax=Zmin=Zmax=Fmin=Fmax=0;
  B0new=B0mesh=1;
  sMax = 1.3;
  doPositionCheck=true;
  doDivergenceCorrection=false;
}

//****************************************************************************
void C3dMesh::freeMesh()
{
  mNormFlux.clear();
  mBfield.clear();
  mgrads.clear();
  mgradB.clear();
  mWhere.clear();
  R.clear();
  Z.clear();
  Fi.clear();
  sMaxSurface.clear();

  meshOK=false;
}

//****************************************************************************
bool C3dMesh::resize(bool array3d)
{
  freeMesh();
  meshOK=true;
// R [nR], Z [nZ], mFi[nFi+2]
  R.init(nR);
  Z.init(nZ);
  Fi.init(-1,nFi,double(0));
  if(array3d) {
    if(!mBfield.resize  (-1,nFi)) meshOK=false;
    if(!mgrads.resize   (-1,nFi)) meshOK=false;
    if(!mgradB.resize   (-1,nFi)) meshOK=false;
    if(!mNormFlux.resize(-1,nFi)) meshOK=false;
    if(!mWhere.resize   (-1,nFi)) meshOK=false;
  }
  if(!meshOK) freeMesh();
  return meshOK;
}

//****************************************************************************
bool C3dMesh::createMeshUsingSymmetry(double dr,double dz,double df, int numThreads)
{
  useSymmetry=true;
  return createMeshEx(dr,dz,df,0,0,numThreads);
}

//****************************************************************************
bool C3dMesh::createMesh(double dr,double dz,double df,int numThreads)
{
  useSymmetry=false;
  return createMeshEx(dr,dz,df,0,0,numThreads);
}

//****************************************************************************
bool C3dMesh::createMeshEx(double dr,double dz,double df,double fi1,double fi2,int nThreads)
{
  if(!isOK()) return false;

  tokamakSymmetry = tokamakSymmetryFluxCoordinates;
  hasGradB = true;

  volatile CStconfig::ctStack::saveState state(ct);

  getExtentLCMS(Rmin, Rmax, Zmin, Zmax);  // run this method to create LCMS
  getExtent(Rmin, Rmax, Zmin, Zmax,sMax);

  if(tokamakSymmetry) {
    double Rmi, Rma, Zmi, Zma;
    getExtent(Rmi , Rma , Zmi , Zma ,1.);
    Rmin=mmax(Rmin,Rmi-5*dr);
    Rmax=mmin(Rmax,Rma+5*dr);
    Zmin=mmax(Zmin,Zmi-5*dz);
    Zmax=mmin(Zmax,Zma+5*dz);
  }

  Rmin-=5*dr;  // add margin
  Rmax+=5*dr;
  Zmin-=5*dz;
  Zmax+=5*dz;

  dR =dr;
  dZ =dz;
  nR = int((Rmax-Rmin)/dR+1);
  nZ = int((Zmax-Zmin)/dZ+1);

  Fmin = mmin(fi1,fi2);
  Fmax = mmax(fi1,fi2);

  double period = useSymmetry?mPeriod/2:mPeriod; // mPeriod is the period size

  fullMesh = (Fmax-Fmin>0.9*period||Fmin==Fmax)?true:false; // if true then create mesh,which can be
                                                            // used for the whole machine

  if(fullMesh==false) useSymmetry=false;  // for safety: useSymmetry only with fullMesh

  if(tokamakSymmetry) {
    fullMesh   =true;
    useSymmetry=false;
    period=1;
    df    =1;
  }

  dFi =  fullMesh?period:(Fmax-Fmin);
  Fmin=  fullMesh?0:Fmin;

  nFi = int(dFi/df+1.5);
  if(nFi<=1) nFi=2;
  dFi = dFi/(nFi-1);   // if(fullMesh) then (nFi-1)*dFi is the same as 0 (periodicity)

  if(!resize()) return false;

  int i,j,k;
  for(i=-1;i<=nFi;i++) Fi.rw()[i]= Fmin+i*dFi;Fmax=Fi[nFi];
  for(j= 0;j< nR; j++) R.rw()[j] = Rmin+j*dR; Rmax=R[nR-1];
  for(k= 0;k< nZ; k++) Z.rw()[k] = Zmin+k*dZ; Zmax=Z[nZ-1];

  edFi = 1/dFi;
  edR  = 1/dR;
  edZ  = 1/dZ;

  B0mesh = CStconfig::getB0();  // |B|
  B0new  = B0mesh;

  truncationLevel=truncate();  // save actual truncation level
  epsALevel=mepsA;             // save actual epsA

  // if(fullMesh&!useSymmetry==true) is for tokamak or stellarator using periodicity
  // otherwise  stellarator using symmetry(fullMesh&useSymmetry==true) or toroidal segment
  int ifirst = (fullMesh&&!useSymmetry)? 0     :-1;
  int ilast  = (fullMesh&&!useSymmetry)?(nFi-2):nFi;

  TESTPRINT( std::cerr<<"C3dMesh: creating mesh planes for cyl. angles from "<<Fi[ifirst]/degree<<" to "<<Fi[ilast]/degree<<" degree"<<std::endl<<std::flush );

#ifdef _WIN32
  #ifndef _MT
  nThreads = 1;  //Singlethread application
  #endif
#endif
  nThreads = nThreads>0?nThreads:1;      // number of threads
  int nCuts = (ilast-ifirst)/nThreads;   // number of cuts per thread
  nCuts    = (nCuts<=0)?1:nCuts;
  nThreads = (nCuts==1)?1:nThreads;      // number of threads
  bool runThreads=(nThreads>1)?true:false;
  meshMutex = ui_mutex_init();
  std::vector<uiThread> aThreads;
  for(int i1=ifirst,i2=0,n=0; ; n++) {
    i2 = (i1+nCuts>ilast)?ilast:(i1+nCuts);
    uiThread th = createMeshPlanes(i1, i2, runThreads,*this);
    if(th!=0 ) { TESTPRINT( std::cerr<<"C3dMesh: thread "<<th<<" is started\n"<<std::flush ); }
    aThreads.push_back(th);  // add to threads array
    i1 = i2+1;
    if(i1>ilast) break;
  }
  // Wait until threads terminate
  for(i=0; i<(int)aThreads.size(); i++) ui_waitthread(aThreads[i]);
  TESTPRINT( std::cerr<<std::endl );

  ui_mutex_destroy(meshMutex);
  aThreads.clear();

  meshOK = true;
  for(i=ifirst; i<=ilast; i++) {  // do we succeed?
    meshOK &= mNormFlux[i].isOK();
    meshOK &= mBfield  [i].isOK();
    meshOK &= mgrads   [i].isOK();
    meshOK &= mgradB   [i].isOK();
    meshOK &= mWhere   [i].isOK();
  }
  if(!meshOK) return meshOK;      // return false if no memory
  copyCuts();
  CreateSmaxSurface();
  return meshOK;
}

//*********************************************************************
//  Routine that begins execution of new thread
//    this function is used as argument in  ui_beginthreadex
extern "C" uiThreadFunc C3dMesh::solveInThread(void * A)
{
  struct Args { C3dMesh *mesh; int i1, i2; };
  Args *arg  =(Args*)A;
  C3dMesh *rootMesh = arg->mesh;
  int i1 = arg->i1;
  int i2 = arg->i2;
  delete arg;
  ui_mutex_lock(rootMesh->meshMutex);
  CStconfig *mc = rootMesh->CStconfig::clone(); //new CStconfig;  *mc = *rootMesh;  // C3dMesh needs a new copy of CStconfig for multithreading to work properly
  ui_mutex_unlock(rootMesh->meshMutex);
  rootMesh->createMeshPlanes(i1, i2, false,*mc);
  delete mc;
  ui_endthreadex();
  return 0;
}

//***************************************************************
// create RZ-mesh-plane for cyl. angles from Fi[i1] to Fi[i2]
// if(startThread==true) run this function in new thread
// @return thread handle or 0 if thread not started or if(startThread==false)
uiThread C3dMesh::createMeshPlanes(int i1, int i2, bool startThread, CStconfig &mc)
{
  if(startThread) {
    struct Args { C3dMesh *mesh; int i1, i2; };
    Args *arg = new Args;
    arg->mesh = this;
    arg->i1 = i1;
    arg->i2 = i2;
    uiThread tId = ui_beginthreadex(solveInThread,(void*)arg);
    if(tId!=uiThread(0)) return tId; // return if thread is started successfully 
    delete arg;
  }
  for(int i=i1; i<=i2; i++) {
    if(!createCut(i, mc)) return uiThread(0);
    if(i%2==0) { TESTPRINT( std::cerr<<"*"<<std::flush ); }
  }
  return uiThread(0);
}

//***************************************************************
// create mesh in RZ-plane for cyl. angle Fi[iCyl]
// @return false if no memory
bool C3dMesh::createCut(int iCyl, CStconfig &mc)
{
  double fiCyl=Fi[iCyl];

  mc.createBorder(sMax, fiCyl,tokamakSymmetry?360:200);

  int r0 = int((mc.getBorderRmin()-Rmin)*edR)-2;
  int r1 = int((mc.getBorderRmax()-Rmin)*edR)+3;
  int z0 = int((mc.getBorderZmin()-Zmin)*edZ)-2;
  int z1 = int((mc.getBorderZmax()-Zmin)*edZ)+3;

  r0=mmax(r0,0);
  r1=mmin(r1,nR-1);
  z0=mmax(z0,0);
  z1=mmin(z1,nZ-1);

  if(!mNormFlux[iCyl].resize(r0,r1,z0,z1)) return false;   // return false if no memory
  if(!mBfield  [iCyl].resize(r0,r1,z0,z1)) return false;   // return false if no memory
  if(!mgrads   [iCyl].resize(r0,r1,z0,z1)) return false;   // return false if no memory
  if(!mgradB   [iCyl].resize(r0,r1,z0,z1)) return false;   // return false if no memory
  if(!mWhere   [iCyl].resize(r0,r1,z0,z1)) return false;   // return false if no memory

  int signBnorm = getSignOfBnormFactor(); //20-08-06

  int sizeZ = z1-z0+1;   // number of elements in z-direction

  Vector3d cylprev(1e7,0.,0.);
  Vector3d magcoord(0.,0.,0.);
  int j,k,k_0,incr,cnt;
//  for(j=r1; j>=r0;j--) {
  for(j=r0; j<=r1;j++) {
    if(j%2) { k_0=z0; incr= 1;}
    else    { k_0=z1; incr=-1;}
    for(k=k_0,cnt=0; cnt<sizeZ; k+=incr,cnt++) { // for(k=z0; k<=z1; k++) {
      Vector3d cyl(R[j],fiCyl,Z[k]);
      double distance;
      Vector3d r;
      bool inside = mc.cylIsInsideBorder(cyl,&distance, &r,&magcoord);
      // for ITER the border can go inside if sMax==1.1 and truncation==2e-5
      if(!inside)  // test against LCMS
        inside = mc.cylIsInside(cyl);  // is inside LCMS?
      if(inside) {                     // if inside border
        magcoord = mc.cyl2mag(cyl);
        mNormFlux[iCyl](j,k) = magcoord[0];
        mBfield  [iCyl](j,k) = mc.getBcyl()*signBnorm;
        mgrads   [iCyl](j,k) = mc.getGrads();
        mgradB   [iCyl](j,k) = mc.getGradBcyl();
        mWhere   [iCyl](j,k) = magcoord[0]<sLCMS?3:2; // is inside LCMS
      }
      else {  // here magcoord is the magnetic coordinates of the nearest point on the border
        mc.NJacobianL(magcoord);
        mNormFlux[iCyl](j,k) = sMax+0.1*distance*edR; // magcoord[0];
        mBfield  [iCyl](j,k) = mc.getBcyl()*signBnorm;
        mgrads   [iCyl](j,k) = mc.getGrads();
        mgradB   [iCyl](j,k) = mc.getGradBcyl();
        mWhere   [iCyl](j,k) = 1;  // outside sMax
      }
    }
  }

  return true;

// find nodes which are outside 'sMax-contour' and has Neighbour inside 'sMax'
  for(j=r0; j<=r1;j++)
    for(k=z0; k<=z1; k++)
      if(hasNeighbourInside(iCyl,j,k)) {
        mWhere   [iCyl](j,k) = 1; // set flag that this node lies outside 'sMax' and has Neighbour inside 'sMax'
        Vector3d cyl(R[j],fiCyl,Z[k]);
        mNormFlux[iCyl](j,k) = mc.cyl2mag(cyl)[0];
        mBfield  [iCyl](j,k) = mc.getBcyl();
        mBfield  [iCyl](j,k) *= signBnorm;   //20-08-06
        mgrads   [iCyl](j,k) = mc.getGrads();
        mgradB   [iCyl](j,k) = mc.getGradBcyl();
      }
  return true;
};

//****************************************************************************
// return: true if node (iCyl,j,k) is outside of sMax and has the neighbour inside sMax
//         false othewise
bool C3dMesh::hasNeighbourInside(int iCyl,int j,int k) const
{
  const CArray2d<char> &node = mWhere[iCyl];
  if(node(j,k)<=1) {  // if outside
    if(node.isIndexOK(j-1,k-1)) if(node(j-1,k-1)>1) return true;
    if(node.isIndexOK(j-1,k  )) if(node(j-1,k  )>1) return true;
    if(node.isIndexOK(j-1,k+1)) if(node(j-1,k+1)>1) return true;
    if(node.isIndexOK(j+1,k-1)) if(node(j+1,k-1)>1) return true;
    if(node.isIndexOK(j+1,k  )) if(node(j+1,k  )>1) return true;
    if(node.isIndexOK(j+1,k+1)) if(node(j+1,k+1)>1) return true;
    if(node.isIndexOK(j  ,k+1)) if(node(j  ,k+1)>1) return true;
    if(node.isIndexOK(j  ,k-1)) if(node(j  ,k-1)>1) return true;
  }
  return false;
}

//****************************************************************************
// Copy cuts, use periodicity
// this function works only for fullMesh==true && useSymmetry=false
void C3dMesh::copyCuts()
{
  if(!meshOK)   return;
  if(!fullMesh) return;  // full mesh means that the mesh can be used for arbitrary cylindrical angles
  if(useSymmetry) return;

// (fullMesh&&!useSymmetry)  here are the tokamak or stellarator using periodicity(not Symmetry)

  mNormFlux[nFi-1] = mNormFlux[0];     meshOK &= mNormFlux[nFi-1].isOK();
  mNormFlux[nFi  ] = mNormFlux[1];     meshOK &= mNormFlux[nFi  ].isOK();
  mNormFlux[-1   ] = mNormFlux[nFi-2]; meshOK &= mNormFlux[   -1].isOK();

  mBfield  [nFi-1] = mBfield  [0];     meshOK &= mBfield[nFi-1].isOK();
  mBfield  [nFi  ] = mBfield  [1];     meshOK &= mBfield[nFi  ].isOK();
  mBfield  [-1   ] = mBfield  [nFi-2]; meshOK &= mBfield[   -1].isOK();

  mgrads   [nFi-1] = mgrads   [0];     meshOK &= mgrads[nFi-1].isOK();
  mgrads   [nFi  ] = mgrads   [1];     meshOK &= mgrads[nFi  ].isOK();
  mgrads   [-1   ] = mgrads   [nFi-2]; meshOK &= mgrads[   -1].isOK();

  mgradB   [nFi-1] = mgradB   [0];     meshOK &= mgradB[nFi-1].isOK();
  mgradB   [nFi  ] = mgradB   [1];     meshOK &= mgradB[nFi  ].isOK();
  mgradB   [-1   ] = mgradB   [nFi-2]; meshOK &= mgradB[   -1].isOK();

  mWhere   [nFi-1] = mWhere   [0];     meshOK &= mWhere[nFi-1].isOK();
  mWhere   [nFi  ] = mWhere   [1];     meshOK &= mWhere[nFi  ].isOK();
  mWhere   [-1   ] = mWhere   [nFi-2]; meshOK &= mWhere[   -1].isOK();
}

//****************************************************************************
bool C3dMesh::writeMesh(const char * fname)  const
{
  if(!meshOK) return false;
  const_cast<C3dMesh*>(this)->calculateVprime();
  const_cast<C3dMesh*>(this)->calculateCurrents();  
  const char * ext = fname_ext(fname);
  int i = (int)strlen(ext);
  if(i==0)  return writeMeshbin (fname,static_cast<float>(0));
  if(     strncmp(ext,".bin4",mmin(i,5))==0) return writeMeshbin (fname,static_cast<float>(0));
  else if(strncmp(ext,".bin8",mmin(i,5))==0) return writeMeshbin (fname,static_cast<double>(0));
  else return writeMeshbin (fname,static_cast<float>(0));
}

//****************************************************************************
// write the mesh and the configuration file
//   see   CStconfig::writebin (fname,szType);
template <class T> bool C3dMesh::writeMeshbin(const char * fname, T szType) const
{
  FILE * fp = fopen(fname,"wb");
  if(fp==NULL) return(false);

  std::string str;

  str = (sizeof(T) == sizeof(double))?"<8":"<4";
  str += "-byte mesh file derived from the file ";
  str += fname_name(mfname.c_str());
  str += (sizeof(T) == sizeof(double))?">double\n":">float \n";

  fwrite(str.c_str(), str.size(),1,fp);

  char buf[512];
  sprintf(buf,"dR=%g,dZ=%g,dFi=%gdegree,sMax=%g,trunc=%g,eps=%g%s",
     dR,dZ,dFi/degree,sMax,truncationLevel,epsALevel,useSymmetry?",symmetry used\n":"\n");

  fwrite(buf, strlen(buf),1,fp);

  size_t  ws = sizeof(szType);
  fwrite(&ws, sizeof(int),1,fp);  // save the size of real number

  double wd[14];
  wd[0] = Rmin;
  wd[1] = Rmax;
  wd[2] = Zmin;
  wd[3] = Zmax;
  wd[4] = Fmin;
  wd[5] = Fmax;
  wd[6] = dR;
  wd[7] = dZ;
  wd[8] = dFi;
  wd[9] = sMax;
  wd[10] = B0mesh;
  wd[11] = truncationLevel;
  wd[12] = epsALevel;
  int flags = (fullMesh?1:0) + (useSymmetry?2:0) + (hasGradB?4:0);
  wd[13] = flags;

  fwrite(wd, sizeof(wd[0]),sizeof(wd)/sizeof(wd[0]),fp);
  fwrite(&nR, sizeof(int),1,fp);
  fwrite(&nZ, sizeof(int),1,fp);
  fwrite(&nFi,sizeof(int),1,fp);

  bool writeOK=true;
  char charType=0;
  int i1=0,i2=0;                        //####
  i2 = fullMesh&&!useSymmetry?(nFi-2):0; //####
  if(0==mBfield.write  (szType,  fp,i1,i2)) writeOK=false;
  if(0==mgrads.write   (szType,  fp,i1,i2)) writeOK=false;
  if(hasGradB) {
    if(0==mgradB.write   (szType,  fp,i1,i2)) writeOK=false;
  }
  if(0==mNormFlux.write(szType,  fp,i1,i2)) writeOK=false;
  if(0==mWhere.write   (charType,fp,i1,i2)) writeOK=false;

  if(writeOK) {
    double  truncSav = truncate();
    (const_cast<C3dMesh*>(this))->truncate(truncationLevel);
    if(sizeof(T)==sizeof(double))
      writeOK=writebin8(fname,fp);
    else
      writeOK=writebin4(fname,fp);
    (const_cast<C3dMesh*>(this))->truncate(truncSav);
  }

  fclose(fp);
  return writeOK;
}

//****************************************************************************
// write the mesh and the configuration file
//   see   CStconfig::writebin (fname,szType);
bool C3dMesh::writeMeshAscii(const char * fname) const
{
  if(!meshOK) return false;

  FILE * fp = fopen(fname,"w");
  if(fp==NULL) return(false);

  std::string str = "# mesh file derived from the file '";
  str += fname_name(mfname.c_str());
  str += "'\n";

  fputs(str.c_str(), fp);
  fprintf(fp,"%13d    ! number of periods\n", nPeriods());
  fprintf(fp,"%13d    ! stellarator symmetry is used if 1\n", useSymmetry?1:0);

  char buf[1024];
  sprintf(buf,"dR=%g,dZ=%g,dFi=%gdegree,sMax=%g,trunc=%g,eps=%g%s",
     dR,dZ,dFi/degree,sMax,truncationLevel,epsALevel,useSymmetry?",symmetry used\n":"\n");
  //fputs(buf, fp);

  int i1=0,i2;                           
  i2 = fullMesh&&!useSymmetry?(nFi-2):0;
  if(fullMesh&&useSymmetry) {
    i1=Fi.lowbound();
    i2=Fi.upperbound();
  }

  fprintf(fp,"# Phi-grid: Phimin+i*dphi, where i=%d..Nphi-1\n",i1);
  fprintf(fp,"%13.5E    ! Phimin,degree\n", Fmin/degree);
  fprintf(fp,"%13.5E    ! dphi,degree\n",   dFi/degree);
  fprintf(fp,"%13d    ! Nphi\n", i2+1);
  fprintf(fp,"# R-grid: Rmin+j*dR, where j=0..Nr-1\n");
  fprintf(fp,"%13.5E    ! Rmin\n", Rmin);
  fprintf(fp,"%13.5E    ! dR\n",   dR);
  fprintf(fp,"%13d    ! Nr\n", nR);
  fprintf(fp,"# Z-grid: Zmin+k*dZ, where k=0..Nz-1\n");
  fprintf(fp,"%13.5E    ! Zmin\n", Zmin);
  fprintf(fp,"%13.5E    ! dZ\n",   dZ);
  fprintf(fp,"%13d    ! Nz\n", nZ);
  fprintf(fp,"#----------------------------\n");


  for(int i=i1; i<=i2; i++) {
    int dim1 = mBfield[i].size1();
    int dim2 = mBfield[i].size2();
    int low1 = mBfield[i].lowbound1();
    int low2 = mBfield[i].lowbound2();
    int upp1 = mBfield[i].upperbound1();
    int upp2 = mBfield[i].upperbound2();

    fprintf(fp,"%13.5E    ! toroidal cut at cylindrical angle (degree)\n",(Fmin+i*dFi)/degree);
    fprintf(fp,"%6d %6d    ! R-grid indexes %6d<=j<=%-6d \n", low1,upp1,low1,upp1);
    fprintf(fp,"%6d %6d    ! Z-grid indexes %6d<=k<=%-6d \n", low2,upp2,low2,upp2);
    fprintf(fp,"%13d    ! number of nodes\n", dim1*dim2);
    fprintf(fp,"#%3s %4s %8s %13s %13s %13s\n", "j","k","s","Br","Bphi","Bz");
//    fprintf(fp,"#%3s %4s %8s %8s %13s %13s %13s\n", "j","k","flag","s","Br","Bphi","Bz");

    for(int j=low1;j<=upp1;j++) {
      for(int k=low2;k<=upp2;k++) {
        fprintf(fp,"%4d %4d ", j,k);
//        fprintf(fp,"%4d ",    mWhere[i](j,k));
        fprintf(fp,"%13.5E ", mNormFlux[i](j,k));
        fprintf(fp,"%13.5E ", mBfield[i](j,k)[0]);
        fprintf(fp,"%13.5E ", mBfield[i](j,k)[1]);
        fprintf(fp,"%13.5E ", mBfield[i](j,k)[2]);
        fprintf(fp,"\n");
     }
    }
  }
  fclose(fp);
  return true;
}

//****************************************************************************
// we use this table to find out the format of a file
static struct {
  CStconfig::fileFormat fmt;  const char * str;
} pattern[] = {
  {CStconfig::W7X_BIN4,"<4-byte mesh file"},
  {CStconfig::W7X_BIN8,"<8-byte mesh file"}
};

//****************************************************************************
// Ascertain(find out) the file format
//
CStconfig::fileFormat C3dMesh::getfileformat(const char * fullname,const char *nam,const char *ext) const
{
  fileFormat fmt=UNDEFINED;
  FILE * fp = fopen(fullname,"r");
  if(fp==NULL) return fmt;

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
  fclose(fp);
  return fmt;
}

//****************************************************************************
// clear all arrays
void C3dMesh::clear()
{
  C3dMesh::freeMesh();
  CRayTrace::freeTrace();
  CStconfig::freeConf();
}

//****************************************************************************
// clear mesh arrays
void C3dMesh::clearMesh()
{
  C3dMesh::freeMesh();
}

//****************************************************************************
// Load mesh file or bc-file or binary bc-file
//
bool C3dMesh::load(const char * fname, CStconfig::fileFormat frmt)
{
  const char * ext = fname_ext (fname);
  const char * nam = fname_name(fname);
  
  doDivergenceCorrection = false;

  this->clear(); // clear all

  fileFormat fmt = getfileformat(fname,nam,ext);  //find out the format of the file
  if     (fmt==W7X_BIN4) loadMeshbin(fname,static_cast<float>(0),fmt);
  else if(fmt==W7X_BIN8) loadMeshbin(fname,static_cast<double>(0),fmt);
  else if(frmt==UNDEFINED) CStconfig::loadfile(fname,fmt);  // not a mesh-file, transfer control to the base class
  else   CStconfig::loadfile(fname,frmt);  // not a mesh-file, transfer control to the base class

  tokamakSymmetry = false;
  if(!CStconfig::mOK) {
     this->clear(); 
  }
  else {
    tokamakSymmetry = tokamakSymmetryFluxCoordinates;
    //tokamakSymmetry |= (nPeriods()==1&&Ntor()==0)?true:false;
    //14Aug2007tokamakSymmetry=(nPeriods()==1&&(Ntor()==1||Ntor()==0))?true:false;
    //?TODO:  tokamakSymmetry=(Ntor()==0)?true:false;
  }

  return CStconfig::mOK;
}

//****************************************************************************
template <class T> bool C3dMesh::loadMeshbin(const char * fname, T szType,fileFormat fmt)
{
  FILE * fp = fopen(fname,"rb");
  if(fp==NULL) return(meshOK=false);

  //CStconfig::freeConf(); // not needed, see this->clear in C3dMesh::load 

  char str[512];
  fgets(str,511,fp);
  fgets(str,511,fp);
  int ws;
  fread(&ws, sizeof(int),1,fp);

  double wd[14];
  fread(wd, sizeof(wd[0]),sizeof(wd)/sizeof(wd[0]),fp);
  Rmin = wd[0];
  Rmax = wd[1];
  Zmin = wd[2];
  Zmax = wd[3];
  Fmin = wd[4];
  Fmax = wd[5];
  dR   = wd[6];
  dZ   = wd[7];
  dFi  = wd[8];
  sMax = wd[9];
  B0mesh = wd[10];
  truncationLevel = wd[11];
  epsALevel = wd[12];
  int flags = int(wd[13]+0.5);
  fullMesh    = (flags&1)>0;
  useSymmetry = (flags&2)>0;
  hasGradB    = (flags&4)>0;

  if(fullMesh==false) useSymmetry=false;  // for safety

  fread(&nR, sizeof(int),1,fp);
  fread(&nZ, sizeof(int),1,fp);
  fread(&nFi,sizeof(int),1,fp);

  if(!resize(false)) { // allocate memory only for R,Fi,Z
    fclose(fp);
    return meshOK=false;
  }

  int i,j,k;
  for(i=-1;i<=nFi;i++) Fi.rw()[i]= Fmin+i*dFi;Fmax=Fi[nFi];
  for(j= 0;j< nR; j++) R.rw()[j] = Rmin+j*dR; Rmax=R[nR-1];
  for(k= 0;k< nZ; k++) Z.rw()[k] = Zmin+k*dZ; Zmax=Z[nZ-1];

  edFi = 1/dFi;
  edR  = 1/dR;
  edZ  = 1/dZ;

  meshOK=true;
  char charType=0;
  if(0==mBfield.read(szType,fp))   meshOK=false;
  if(0==mgrads.read(szType,fp))    meshOK=false;
  if(hasGradB) {
    if(0==mgradB.read(szType,fp))    meshOK=false;
  }
  if(0==mNormFlux.read(szType,fp)) meshOK=false;
  if(0==mWhere.read(charType,fp))  meshOK=false;

  copyCuts(); // ####

  if(meshOK) {
    if(sizeof(T)==sizeof(double))
      meshOK=loadbin8(fname,fp);
    else
      meshOK=loadbin4(fname,fp);
  }
  else
    CStconfig::mOK=false;  // if mesh loading is failed then set also bad flag for configuration

  fclose(fp);

  if(initAfterLoad(fmt)) B0new=CStconfig::getB0();
  if(CStconfig::mOK) setAccuracy(epsALevel);
  if(CStconfig::mOK) CreateSmaxSurface();

  BnormMesh = B0new/B0mesh;  // set scale factor for magn.field
  BnormMesh *= CStconfig::getSignOfBnormFactor();  //20-08-06

  return CStconfig::mOK;
}

// ****************************************************************************
bool C3dMesh::CreateSmaxSurface()
{

  int nTheta=120;

  shift = 1.e-4/pi;

  if(!sMaxSurface.resize(0,nTheta-1,-1,nFi)) return meshOK=false;

  double dtheta = twopi/(nTheta-1);
  int i,j;
  for(i=0; i<nTheta-1; i++) {          // poloidal direction
    double theta = i*dtheta + shift;
    for(j=-1; j<=nFi; j++) {           // toroidal direction
      double phi = j*dFi + shift;
      Vector3d mixcrd(sMax,theta,phi);
      sMaxSurface(i,j)=mixcoord2xyz(mixcrd); // point (x,y,z) on surface 'sMax'
    }
  }

  for(j=-1; j<=nFi; j++) sMaxSurface(nTheta-1,j) = sMaxSurface(0,j);

  return meshOK=true;
}

//*********************************************************************
// INPUT:
//   Vector3d xyz == (x,y,z) is the point in Cartesian coordinates,
// OUTPUT:
//  function returns 'true' if the point lies inside LCMS
bool C3dMesh::M3DxyzIsInsideMesh(const Vector3d &xyz,double * distance) const
{
  double rc = sqrt(xyz.x()*xyz.x()+xyz.y()*xyz.y());
  double fi = mod2pi(atan2(xyz.y(),xyz.x()));
  double zc = xyz.z();
  Vector3d cyl(rc,fi,zc);
  return M3DcylIsInsideMesh(cyl,distance);
}

//*********************************************************************
// INPUT:
//   Vector3d cylP is a point in cylindrical coordinates,
// OUTPUT:
//  function returns 'true' if the point lies inside sMax
//
bool C3dMesh::M3DcylIsInsideMesh(const Vector3d &cylP, double * distance) const
{
  if(!meshOK) return false;
  if(!sMaxSurface.isOK()) return false;

  Vector3d cyl(cylP);
  double f  = fullMesh?modPeriod(cyl[1]-shift) : (cyl[1]-shift);

  if(tokamakSymmetry) {
    f=0;
  }
  else
    if(useSymmetry&&f>mPeriod/2) {  // use stellarator symmetry property
      f = mPeriod-f;
      cyl[2]=-cyl[2];
    }

  cyl[1]=f;

  if(distance==NULL)
    if( cyl[0]<Rmin||cyl[0]>Rmax||cyl[2]<Zmin||cyl[2]>Zmax ) return false; // return if outside

  return cylIsInsideEx(cyl,sMaxSurface,dFi,distance);
}

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
};

inline static void LagrangeCoeffs(double &p, double *c)
{
  c[0] =-(p-1)*(p-2)*(p-3)/6;
  c[1] = p*(p-2)*(p-3)/2;
  c[2] =-p*(p-1)*(p-3)/2;
  c[3] = p*(p-1)*(p-2)/6;
}

//***************************************************************
// function finds the integer coordinates on the mesh
//  and calculates coefficients for interpolation
bool C3dMesh::FindMeshCoordinates(const Vector3d &cylcoord, bool noDerivatives)
{
  if(!meshOK) return false;
  Vector3d cyl(cylcoord);

  double f  = fullMesh?modPeriod(cyl[1]) : cyl[1];

  symmetryActivated=false;  //!$!
  if(tokamakSymmetry) {
    f=0;
  }
  else
    if(useSymmetry&&f>mPeriod/2) {  // use stellarator symmetry property
      symmetryActivated=true;
      f = mPeriod-f;
      cyl[2]=-cyl[2];
    }                       //!$!

  cyl[1]=f;

  if( cyl[0]<Rmin||cyl[0]>Rmax||cyl[2]<Zmin||cyl[2]>Zmax ) return false; // return if outside

  i0 = int((f     -Fmin)*edFi); // int fi
  j0 = int((cyl[0]-Rmin)*edR);  // int R
  k0 = int((cyl[2]-Zmin)*edZ);  // int Z

  if(!fullMesh) {
    if(i0==nFi-1) i0--;
    if(i0<0||i0>nFi-2) {
      std::cerr <<"C3dMesh: Cylindrical angle is out of definition range\n"; // 0<=i0<=nFi-2 must be always
      return false;   // return if outside
    }
  }
  else
    if(i0==nFi-1) i0--;  //it's needed for case useSymmetry&&cyl[1]==mPeriod/2

// set the 1st point for the Lagrange interpolation
  i0--;
  j0--;
  k0--;

  int i;
  for(i=i0;i<i0+4;i++) {
    const CArray2d<char> &w= mWhere[i];
    j0=mmax(j0,w.lowbound1());
    j0=mmin(j0,w.upperbound1()-3);
    k0=mmax(k0,w.lowbound2());
    k0=mmin(k0,w.upperbound2()-3);
  }

  double pf=(f     -Fi[i0])*edFi;
  double pr=(cyl[0]-R [j0])*edR;
  double pz=(cyl[2]-Z [k0])*edZ;
  if(noDerivatives) {
    LagrangeCoeffs(pf, cfi);
    LagrangeCoeffs(pr, cr );
    LagrangeCoeffs(pz, cz );
  }
  else {
    LagrangeCoeffs(pf, cfi,cfid);
    LagrangeCoeffs(pr, cr, crd );
    LagrangeCoeffs(pz, cz, czd );
  }
  BnormMesh = B0new/B0mesh;  // set scale factor for magn.field
  BnormMesh *= CStconfig::getSignOfBnormFactor();  //20-08-06

  if(!doPositionCheck) return true;

// Scan interpolation cube, look for nodes which were not visited
  int cubeHasNodesInsideSmax = 0;
  int cubeHasNodesVisited = 0;
  int j,k;
  // we must have at least one node inside sMax-surface
  for(i=i0;i<i0+4;i++) {    // loop over cube from (i0,j0,k0) to (i0+4,j0+4,k0+4)
    const CArray2d<char> &where = mWhere[i];
    for(j=j0;j<j0+4; j++) {
      const char * w = where[j];
      for(k=k0;k<k0+4; k++) {
        if(w[k]>1) cubeHasNodesInsideSmax++;  // if node inside sMax-surface
        if(w[k]>0) cubeHasNodesVisited++;     // if node was visited
      }
    }
  }
  if(cubeHasNodesInsideSmax<64) {               // some nodes are not inside sMax if less then 4*4*4
    if(cubeHasNodesInsideSmax==0) return false; // can not be interpolated, too far from sMax-surface
    if(doPositionCheck) {
      double distance=0;
      cylIsInsideEx(cyl,sMaxSurface,dFi,&distance);
      if(distance>dR*0.5)  return false; // return if cyl outside sMax-surface
    }
    if(cubeHasNodesVisited<64) {
      double epsSav   = mepsA;
      double B0saved  = CStconfig::getB0();  // save current B0
      int signBnorm = CStconfig::getSignOfBnormFactor(); //20-08-06
      B0saved  *= signBnorm;                 // save current B0 with sign
      CStconfig::setB0(B0mesh);              // back to mesh initial B
      signBnorm = CStconfig::getSignOfBnormFactor(); //20-08-06
      setAccuracy(epsALevel);
    // scan interpolation cube, look for nodes which were not visited
      bool jacobianDataWereSaved=false;
      volatile ctStack::saveState state(CStconfig::ct);
      for(i=i0;i<i0+4;i++) {    // loop over cube from (i0,j0,k0) to (i0+4,j0+4,k0+4)
        for(j=j0;j<j0+4; j++)
          for(k=k0;k<k0+4; k++)
            if(mWhere[i](j,k)==0) { // if node[i](j,k) was not visited then
              if(!jacobianDataWereSaved) { // if first time here
                //19-Apr-2008  saveNJacData();
                jacobianDataWereSaved=true;
              }
              Vector3d c(R[j],Fi[i],Z[k]);
              mNormFlux[i](j,k) = CStconfig::cyl2mag(c)[0];
              mBfield  [i](j,k) = CStconfig::getBcyl();
              mBfield  [i](j,k) *= signBnorm;   //20-08-06
              mgrads   [i](j,k) = CStconfig::getGrads();
              mgradB   [i](j,k) = CStconfig::getGradBcyl();
              mWhere   [i](j,k) = 1;   // set flag that this node lies outside 'sMax' and has Neighbour inside 'sMax'
            }
      }
      //19-Apr-2008  if(jacobianDataWereSaved) restoreNJacData();
      setAccuracy(epsSav);             // restore
      CStconfig::setB0(B0saved);       // restore B0
    }
  }

  return true;
}

#define FIND_COORD   (const_cast<C3dMesh*>(this))->FindMeshCoordinates

//***************************************************************
// Find flux label for the point with cylindrical coordinates 'cyl'
//  function returns  s - normalized tor. flux
double C3dMesh::M3Dcyl2s(const Vector3d &cyl) const
{
  double s(0);
  if(!FIND_COORD(cyl,NoDerivatives)) return 1001;

  int i,j,k;
  for(i=0;i<4;i++) {    // loop over cube (i0,j0,k0)..(i0+4,j0+4,k0+4)
    const CArray2d<double> &S1 = mNormFlux[i+i0];
    for(j=0;j<4;j++) {
      const double   *S2 = S1[j+j0];
      double w=cfi[i]*cr[j];
      for(k=0;k<4;k++)
        s += S2[k+k0]*w*cz[k];
    }
  }
  if(s<0) {
    double epsSav   = mepsA;
    const_cast<C3dMesh*>(this)->setAccuracy(1e-9);
    s=const_cast<C3dMesh*>(this)->cyl2s(cyl);
    const_cast<C3dMesh*>(this)->setAccuracy(epsSav);
  }
  return (s<=sMax)?s:1001;
}

//***************************************************************
// Find flux label for the point with cartesian coordinates 'xyz'
//  function returns  s - normalized tor. flux
double C3dMesh::M3Dxyz2s(const Vector3d &xyz) const
{
  Vector3d cyl(sqrt(xyz.x()*xyz.x()+xyz.y()*xyz.y()),
               mod2pi(atan2(xyz.y(),xyz.x())),
               xyz.z());
  return M3Dcyl2s(cyl);
}

//***************************************************************
//   function retuns grad(|B|) = Vector3d(dB/dr, dB/dfi/r, dB/dz)
Vector3d C3dMesh::M3DgetGradB(const Vector3d &cyl) const
{
  Vector3d gB(0.,0.,0.);
  if(!FIND_COORD(cyl,NoDerivatives)) return gB; // return zero if point lies outside
  if(!hasGradB) return gB;

  int i,j,k;
  for(i=0;i<4;i++) {    // loop over cube (i0,j0,k0)..(i0+4,j0+4,k0+4)
    const CArray2d<Vector3d> &G1 = mgradB[i+i0];
    for(j=0;j<4;j++) {
      const Vector3d *G2 = G1[j+j0];
      for(k=0;k<4;k++)
        gB += G2[k+k0]*(cfi[i]*cr[j]*cz[k]);
    }
  }
  gB *= BnormMesh;
  if(symmetryActivated) {  //!$!
    gB[1]=-gB[1];
    gB[2]=-gB[2];
  };
  return gB;
}

//***************************************************************
//   function retuns grad(|B|) = Vector3d(dB/dx, dB/dy, dB/dz)
Vector3d C3dMesh::M3DgetGradBxyz(const Vector3d &xyz) const
{
  Vector3d cyl = xyz.toCylindrical();
  Vector3d g = M3DgetGradB(cyl);
// transform to Cartesian coordinates and return
  return g.cylVector2Cartesian(cyl);
}

//***************************************************************
// cylindrical coordinates are used
//   function retuns grad(s) = Vector3d(ds/dr, ds/dfi/r, ds/dz)
// this function is more accurate then M3DgetGrads1
Vector3d C3dMesh::M3DgetGrads(const Vector3d &cyl) const
{
  Vector3d gs(0.,0.,0.);
  if(!FIND_COORD(cyl,NoDerivatives)) return gs; // return zero if point lies outside

  int i,j,k;
  for(i=0;i<4;i++) {    // loop over cube (i0,j0,k0)..(i0+4,j0+4,k0+4)
    const CArray2d<Vector3d> & G1 = mgrads   [i+i0];
    for(j=0;j<4;j++) {
      const Vector3d * G2 = G1[j+j0];
      for(k=0;k<4;k++)
        gs += G2[k+k0]*(cfi[i]*cr[j]*cz[k]);
    }
  }

  if(symmetryActivated) {  //!$!
    gs[1]=-gs[1];
    gs[2]=-gs[2];
  };
  return gs;
}

//***************************************************************
Vector3d C3dMesh::M3DgetGradsxyz(const Vector3d &xyz) const
{
  Vector3d cyl(sqrt(xyz.x()*xyz.x()+xyz.y()*xyz.y()),
               mod2pi(atan2(xyz.y(),xyz.x())),
               xyz.z());
  Vector3d g = M3DgetGrads(cyl);
// transform to Cartesian coordinates
  double fi = cyl[1];
  double cs = ::cos(fi);
  double sn = ::sin(fi);
  double gx = g[0]*cs-g[1]*sn;  // gr*cos(fi)-gfi*sin(fi);
  double gy = g[0]*sn+g[1]*cs;  // gr*sin(fi)+gfi*cos(fi);

  return Vector3d(gx,gy,g[2]);
}

//***************************************************************
// cylindrical coordinates are used
//   function retuns B in  cylindrical coordinates
// return 0 if outside sMax-surface
Vector3d C3dMesh::M3DgetBcyl(const Vector3d &cyl) const
{
  Vector3d b(0.,0.,0.);
  if(!FIND_COORD(cyl,NoDerivatives)) return b; // return zero if point lies outside
  int i,j,k;
  for(i=0;i<4;i++) {    // loop over cube (i0,j0,k0)..(i0+4,j0+4,k0+4)
    const CArray2d<Vector3d> & B1 = mBfield[i+i0];
    for(j=0;j<4;j++) {
      const Vector3d * B2 = B1[j+j0];
      for(k=0;k<4;k++)
        b += B2[k+k0]*(cfi[i]*cr[j]*cz[k]);
    }
  }
  b *= BnormMesh;
  if(symmetryActivated)
    b[0] =-b[0];  // Br = - Br
  return b;
}

//***************************************************************
// cylindrical coordinates are used
//   function retuns B in  cylindrical coordinates
// return 0 if outside sMax-surface
Vector3d C3dMesh::M3DgetBcyl(const Vector3d &cyl, double &s) const
{
  Vector3d b(0.,0.,0.);
  s = 1001;
  if(!FIND_COORD(cyl,NoDerivatives)) return b; // return zero if point lies outside
  s = 0;
  int i,j,k;
  for(i=0;i<4;i++) {    // loop over cube (i0,j0,k0)..(i0+4,j0+4,k0+4)
    const CArray2d<Vector3d> & B1 = mBfield[i+i0];
    const CArray2d<double>  &S1 = mNormFlux[i+i0];
    for(j=0;j<4;j++) {
      const Vector3d * B2 = B1[j+j0];
      const double   *S2 = S1[j+j0];
      double w=cfi[i]*cr[j];
      for(k=0;k<4;k++) {
        double w1 = w*cz[k];
        s += S2[k+k0]*w1;
        b += B2[k+k0]*w1;
      }
    }
  }
  if(s<0) {
    double epsSav   = mepsA;
    const_cast<C3dMesh*>(this)->setAccuracy(1e-9);
    s=const_cast<C3dMesh*>(this)->cyl2s(cyl);
    const_cast<C3dMesh*>(this)->setAccuracy(epsSav);
  }
  s = (s<=sMax)?s:1001;
  b *= BnormMesh;
  if(symmetryActivated)
    b[0] =-b[0];  // Br = - Br
  return b;
}

//***************************************************************
// cartesian coordinates are used
//   function retuns B in  cartesian coordinates
Vector3d C3dMesh::M3DgetBxyz(const Vector3d &xyz) const
{
  Vector3d cyl(sqrt(xyz.x()*xyz.x()+xyz.y()*xyz.y()),
               mod2pi(atan2(xyz.y(),xyz.x())),
               xyz.z());
  Vector3d b = M3DgetBcyl(cyl);
// transform to Cartesian coordinates
  double fi = cyl[1];
  double cs = ::cos(fi);
  double sn = ::sin(fi);
  double Bx = b[0]*cs-b[1]*sn;  // Br*cos(fi)-Bfi*sin(fi);
  double By = b[0]*sn+b[1]*cs;  // Br*sin(fi)+Bfi*cos(fi);

  return Vector3d(Bx,By,b[2]);
}

//***************************************************************
// cartesian coordinates are used
//   function retuns B in  cartesian coordinates
Vector3d C3dMesh::M3DgetBxyz(const Vector3d &xyz, double &s) const
{
  Vector3d cyl(sqrt(xyz.x()*xyz.x()+xyz.y()*xyz.y()),
               mod2pi(atan2(xyz.y(),xyz.x())),
               xyz.z());
  Vector3d b = M3DgetBcyl(cyl,s);
// transform to Cartesian coordinates
  double fi = cyl[1];
  double cs = ::cos(fi);
  double sn = ::sin(fi);
  double Bx = b[0]*cs-b[1]*sn;  // Br*cos(fi)-Bfi*sin(fi);
  double By = b[0]*sn+b[1]*cs;  // Br*sin(fi)+Bfi*cos(fi);

  return Vector3d(Bx,By,b[2]);
}

void C3dMesh::M3DdivergenceFree(const Vector3d &cyl,Vector3d &B,Vector3d &dBdr,Vector3d &dBdfir,Vector3d &dBdz) const
{
  if(tokamakSymmetry) return;
  const int r=0,fi=1,z=2;  // indexes
  dBdfir[fi] = -dBdr[r] - dBdz[z] - B[r]/cyl[r];  //dBf/df/R = -(dBr/dr +dBz/dz) -Br/R

  //curlB = 0;
  //dBdr[fi] =  dBdfir[r]-B[fi]/cyl[r];             //dBf/dr =  dBr/df/R-Bf/R
  //dBdz[fi] =  dBdfir[z];                          //dBf/dz =  dBz/df/R
  //dBdr[z ] =  dBdz[r];                            //dBz/dr =  dBr/dz
}

//***************************************************************
// get B and derivatives dB/dr, dB/dfi/r, dB/dz in point cyl
//  Input:
//    Vector3d cyl -- cylindrical coordinates
//  Output:
//    Vector3d B -- cylindrical coordinates of B
//    Vector3d gradB   -- cylindrical coordinates of grad(|B|)
//    Vector3d grads   -- cylindrical coordinates of grad(s)
//  Return:   s
//
double C3dMesh::M3DgetBandGradB(const Vector3d &cyl,Vector3d &B, Vector3d &gradB, Vector3d &grads) const
{
  double s(0);
  B = 0;
  gradB=0;
  grads=0;
  Vector3d &b  = B; 
  Vector3d &gB = gradB; 
  Vector3d &gs = grads; 
  if(!FIND_COORD(cyl,(hasGradB?NoDerivatives:false))) return 1001; // return zero if point lies outside

  int i,j,k;
  for(i=0;i<4;i++) {    // loop over cube (i0,j0,k0)..(i0+4,j0+4,k0+4)
    const CArray2d<Vector3d> &g1 = mgradB [i+i0];
    const CArray2d<Vector3d> &G1 = mgrads [i+i0];
    const CArray2d<Vector3d> &B1 = mBfield[i+i0];
    const CArray2d<double>   &S1 = mNormFlux[i+i0];
    for(j=0;j<4;j++) {
      const Vector3d *g2 = g1[j+j0];
      const Vector3d *G2 = G1[j+j0];
      const Vector3d *B2 = B1[j+j0];
      const double   *S2 = S1[j+j0];
      double w=cfi[i]*cr[j];
      for(k=0;k<4;k++) {
        if(!hasGradB) {
          double B = B2[k+k0].abs();
          gB[0] += B*(cfi [i]*crd[j]*cz [k]);  // dB/dr
          gB[1] += B*(cfid[i]*cr [j]*cz [k]);  // dB/dfi
          gB[2] += B*w*czd[k];
        }
        else {
          gB += g2[k+k0]*w*cz[k];
        }
        b     += B2[k+k0]*w*cz[k];
        gs    += G2[k+k0]*w*cz[k];
        s     += S2[k+k0]*w*cz[k];
     }
    }
  }
  gB *= BnormMesh;
  b  *= BnormMesh;
  if(!hasGradB) {
    gB[0] *= edR;         // dB/dr
    gB[1] *= edFi/cyl[0]; // dB/dfi/r
    gB[2] *= edZ;         // dB/dz
  }
  if(symmetryActivated) {  //!$!
    b[0] =-b[0];
    gs[1]=-gs[1];
    gs[2]=-gs[2];
    gB[1]=-gB[1];
    gB[2]=-gB[2];
  }

  if(s<0) {
    double epsSav   = mepsA;
    const_cast<C3dMesh*>(this)->setAccuracy(1e-9);
    s=const_cast<C3dMesh*>(this)->cyl2s(cyl);
    const_cast<C3dMesh*>(this)->setAccuracy(epsSav);
  }
  return (s<=sMax)?s:1001;
}

//***************************************************************
//  Input:
//    Vector3d xyz -- cartesian coordinates
//  Output:
//    Vector3d B -- cartesian coordinates of B
//    Vector3d gradB   -- cartesian coordinates of grad(|B|)
//    Vector3d grads   -- cartesian coordinates of grad(s)
//  Return:   s
//
double C3dMesh::M3DgetBandGradBxyz(const Vector3d &xyz,Vector3d &B, Vector3d &gradB, Vector3d &grads) const
{
  Vector3d cyl = xyz.toCylindrical();
  Vector3d g = M3DgetGradB(cyl);
  double s = M3DgetBandGradB(cyl,B, gradB, grads);
// transform to Cartesian coordinates and return
  B     = B.cylVector2Cartesian(cyl);
  gradB = gradB.cylVector2Cartesian(cyl);
  grads = grads.cylVector2Cartesian(cyl);
  return s;
}

//***************************************************************
// get B and derivatives dB/dr, dB/dfi/r, dB/dz in point cyl
//  Input:
//    Vector3d cyl -- cylindrical coordinates
//  Output:
//    Vector3d dBdr   -- dB/dr, where dBdr[0]=dBr/dr, dBdr[1]=dBfi/dr, dBdr[2]=dBz/dr
//    Vector3d dBdfir -- dB/dfi/r
//    Vector3d dBdz   -- dB/dz
//  Return:
//    Vector3d B -- cylindrical coordinates of B
//
Vector3d C3dMesh::M3DgetdBcyl(const Vector3d &cyl,Vector3d &dBdr,Vector3d &dBdfir,Vector3d &dBdz) const
{
  Vector3d b(0.,0.,0.);
  dBdr  = 0; // dB/dr
  dBdfir= 0; // dB/dfi
  dBdz  = 0; // dB/dz
  if(!FIND_COORD(cyl)) return b; // return zero if point lies outside

  int i,j,k;
  for(i=0;i<4;i++) {    // loop over cube (i0,j0,k0)..(i0+4,j0+4,k0+4)
    const CArray2d<Vector3d> & B1 = mBfield[i+i0];
    for(j=0;j<4;j++) {
      const Vector3d * B2 = B1[j+j0];
      for(k=0;k<4;k++) {
        b     += B2[k+k0]*(cfi [i]*cr [j]*cz [k]);
        dBdr  += B2[k+k0]*(cfi [i]*crd[j]*cz [k]);
        dBdfir+= B2[k+k0]*(cfid[i]*cr [j]*cz [k]);
        dBdz  += B2[k+k0]*(cfi [i]*cr [j]*czd[k]);
      }
    }
  }
  dBdr  *= edR;         // dB/dr
  dBdfir*= edFi/cyl[0]; // dB/dfi/r
  dBdz  *= edZ;         // dB/dz
  dBdr  *= BnormMesh;
  dBdfir*= BnormMesh;
  dBdz  *= BnormMesh;
  b     *= BnormMesh;
  if(symmetryActivated) {  //!$!
    b[0]     =-b[0];
    dBdr[0]  =-dBdr[0];
    dBdfir[1]=-dBdfir[1];
    dBdfir[2]=-dBdfir[2];
    dBdz  [1]=-dBdz  [1];
    dBdz  [2]=-dBdz  [2];
  };
  if(doDivergenceCorrection) M3DdivergenceFree(cyl,b,dBdr,dBdfir,dBdz);
  return b;
}

//***************************************************************
// get B and derivatives dB/dx,dB/dy,dB/dz in point xyz
//  Input:
//    Vector3d xyz -- cartesian coordinates
//  Output:
//    Vector3d dBdx -- dB/dx, where dBdx(1)=dBx/dx, dBdx(2)=dBy/dx, dBdx(3)=dBz/dx
//    Vector3d dBdy -- dB/dy
//    Vector3d dBdz -- dB/dz
//  Return:
//    Vector3d B -- cartesian coordinates of B
//
Vector3d C3dMesh::M3DgetdBxyz(const Vector3d &xyz,Vector3d &dBdx,Vector3d &dBdy,Vector3d &dBdz) const
{
  Vector3d cyl(sqrt(xyz.x()*xyz.x()+xyz.y()*xyz.y()),
               mod2pi(atan2(xyz.y(),xyz.x())),
               xyz.z());
  Vector3d B,b,dbdr,dbdfir,dbdz;

  B = 0;
  dBdx = 0; // dB/dx
  dBdy = 0; // dB/dy
  dBdz = 0; // dB/dz

  b=M3DgetdBcyl(cyl,dbdr,dbdfir,dbdz);

  if(b==B) return b; // return zeros if point lies outside

// transform to Cartesian coordinates
  double r  = cyl[0];
  double fi = cyl[1];
  double cs = ::cos(fi);
  double sn = ::sin(fi);
  double Bx = b[0]*cs-b[1]*sn;  // Bx=Br*cos(fi)-Bfi*sin(fi);
  double By = b[0]*sn+b[1]*cs;  // By=Br*sin(fi)+Bfi*cos(fi);
  double Bz = b[2];

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

  return Vector3d(Bx,By,Bz);
}

//***************************************************************
// get B and derivatives dB/dx,dB/dy,dB/dz in point xyz
//  Input:
//    Vector3d xyz -- cartesian coordinates
//  Output:
//    Vector3d dBdx -- dB/dx, where dBdx(1)=dBx/dx, dBdx(2)=dBy/dx, dBdx(3)=dBx/dx
//    Vector3d dBdy -- dB/dy
//    Vector3d dBdz -- dB/dz
//    Vector3d B    -- cartesian coordinates of B
//    Vector3d grads -- cartesian coordinates of grad(s)
//  Return:   s
//
double C3dMesh::M3DgetdBGradsxyz(const Vector3d &xyz,Vector3d &B,Vector3d &dBdx,Vector3d &dBdy,Vector3d &dBdz, Vector3d &grads ) const
{
  Vector3d cyl(sqrt(xyz.x()*xyz.x()+xyz.y()*xyz.y()),
               mod2pi(atan2(xyz.y(),xyz.x())),
               xyz.z());
  Vector3d b,dbdr,dbdfir,dbdz,gs;
  B    = 0;
  dBdx = 0; // dB/dx
  dBdy = 0; // dB/dy
  dBdz = 0; // dB/dz
  grads= 0;

  double s = M3DgetdBGradscyl(cyl,b,dbdr,dbdfir,dbdz,gs);
  if(s>sMax) return 1001;

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

  return (s<=sMax)?s:1001;
}


//***************************************************************
double C3dMesh::M3DgetdBGradscyl(const Vector3d &cyl,Vector3d &B,Vector3d &dBdr,Vector3d &dBdfir,Vector3d &dBdz,Vector3d &grads) const
{
  double s(0);
  B     = 0;
  dBdr  = 0; // dB/dr
  dBdfir= 0; // dB/dfi/r
  dBdz  = 0; // dB/dz
  grads = 0;
  if(!FIND_COORD(cyl)) return 1001; // return zeros if point lies outside

  int i,j,k;
  for(i=0;i<4;i++) {    // loop over cube (i0,j0,k0)..(i0+4,j0+4,k0+4)
    const CArray2d<Vector3d> & B1 = mBfield  [i+i0];
    const CArray2d<Vector3d> & G1 = mgrads   [i+i0];
    const CArray2d<double>   & S1 = mNormFlux[i+i0];
    for(j=0;j<4;j++) {
      const Vector3d * B2 = B1[j+j0];
      const Vector3d * G2 = G1[j+j0];
      const double   * S2 = S1[j+j0];
      for(k=0;k<4;k++) {
        B     += B2[k+k0]*(cfi [i]*cr [j]*cz [k]);
        dBdr  += B2[k+k0]*(cfi [i]*crd[j]*cz [k]);
        dBdfir+= B2[k+k0]*(cfid[i]*cr [j]*cz [k]);
        dBdz  += B2[k+k0]*(cfi [i]*cr [j]*czd[k]);
        grads += G2[k+k0]*(cfi [i]*cr [j]*cz [k]);
        s     += S2[k+k0]*(cfi [i]*cr [j]*cz [k]);
      }
    }
  }
  dBdr   *= edR;         // dB/dr
  dBdfir *= edFi/cyl[0]; // dB/dfi/r
  dBdz   *= edZ;         // dB/dz
  dBdr   *= BnormMesh;
  dBdfir *= BnormMesh;
  dBdz   *= BnormMesh;
  B      *= BnormMesh;
  if(symmetryActivated) {  //!$!
    B[0]     =-B[0];
    dBdr[0]  =-dBdr[0];
    dBdfir[1]=-dBdfir[1];
    dBdfir[2]=-dBdfir[2];
    dBdz  [1]=-dBdz  [1];
    dBdz  [2]=-dBdz  [2];
    grads [1]=-grads [1];
    grads [2]=-grads [2];
  };

  if(doDivergenceCorrection) M3DdivergenceFree(cyl,B,dBdr,dBdfir,dBdz);

  if(s<0) {
    double epsSav   = mepsA;
    const_cast<C3dMesh*>(this)->setAccuracy(1e-9);
    s=const_cast<C3dMesh*>(this)->cyl2s(cyl);
    const_cast<C3dMesh*>(this)->setAccuracy(epsSav);
  }
  return (s<=sMax)?s:1001;
}

//****************************************************************************
// Set normalization factor for value of the magnetic field
// Set minimum magnetic field on axis equals  B
void C3dMesh::setB0(double B)
{
  CStconfig::setB0(B);
  B0new = CStconfig::getB0();
}

//****************************************************************************
// Set normalization factor the magnetic field
// restore initial value of the magnetic field on axis
void C3dMesh::restoreB0()
{
  CStconfig::restoreB0();
  B0new = CStconfig::getB0();
}

//****************************************************************************
// Set normalization factor for value of the magnetic field
// INPUT:
// B -- value of magnetic field on axis at angle ficyl
void C3dMesh::setB0(double B,double ficyl)
{
  CStconfig::setB0(B,ficyl);
  B0new = CStconfig::getB0();
}

//****************************************************************************
// Set the value of B_00 on the axis
// INPUT:
// B00 -- B_00 value of magnetic field on axis
void C3dMesh::setB00(double B00)
{
  CStconfig::setB00(B00);
  B0new = CStconfig::getB0();
}

//****************************************************************************
// Normalize magnetic field
// Set Ipol at LCMS equals I in Ampers
void C3dMesh::setIcoil(double I)
{
  CStconfig::setIcoil(I);
  B0new = CStconfig::getB0();
}

double C3dMesh::getBscaleFactor()
{
  return CStconfig::getBscaleFactor();
}

//****************************************************************************
// Set normalization factor for value of the magnetic field
void C3dMesh::scaleB(double multiplier)
{
  CStconfig::scaleB(multiplier);
  B0new = CStconfig::getB0();
}

#if 0
//****************************************************************************
// Set normalization factor for the value of the magnetic field
// INPUT:
// Bphi0 -- average value of Bphi on axis
void C3dMesh::setAverageBphi0(double Bphi0)
{
  CStconfig::setAverageBphi0(Bphi0);
  B0new = CStconfig::getB0();
}
#endif

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
bool C3dMesh::M3DgetRayEntryPoint(const Vector3d &r0,const Vector3d &rd,Vector3d &entryPoint) const
{
  const double accuracy = 1e-3;
  if(!meshOK) return false;  // return if no mesh
  if(xyzIsInside(r0)) return false;  // return if the origin lies inside LCMS

  double dl = rminor()*0.04;     // step 2cm for W7X
  Vector3d dr = dl*rd/rd.abs();

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

// we are near LCMS, bracket the LCMS between two point rOutside and rInside
  bool seachResult = false;
  bool posCheck = doPositionCheck;   // save flag
  (const_cast<C3dMesh*>(this))->doPositionCheck = false;
  Vector3d rOutside = r; // 'nearest' outside point
  Vector3d rInside;
  N  = 32;               // max number of iterations
  if(M3Dxyz2s(r)<=1) {   // if current position lies inside
    do {                 // cycle while inside LCMS
      r -= dr;           // move backward along the ray
      if(--N < 1) goto ret;
    } while(M3Dxyz2s(r) < 1); // cycle while inside LCMS
    rOutside = r;             // 'nearest' outside point
    rInside = r+dr;           // move inside and save position
  }
  else  {                     // if current position lies slightly outside
    do {                      // cycle while outside LCMS
      r += dr;                // move forward along ray
      if(--N < 1) goto ret;
    } while(M3Dxyz2s(r) > 1); // cycle while outside LCMS
    rInside = r;              // inside, save position
  }
  N = 100;
// Find intersection of the Ray with the last surface by bisection
  while(--N) {
    r = (rInside+rOutside)/2;
    double ds = 1 - M3Dxyz2s(r);
    if(ds< 0.)
      rOutside   = r;
    else {
      rInside = r;             // if r lies inside
      if(ds<=accuracy) break;  // Ok, done
    }
  }
  seachResult=true;
  entryPoint=rInside;
ret:
  (const_cast<C3dMesh*>(this))->doPositionCheck = posCheck;   // restore flag
  return seachResult;
}



//***************************************************************
//***************************************************************
//*************OBSOLETE******************************************
//***************************************************************
//***************************************************************


//***************************************************************
// INPUT:
//  Vector3d cyl -- point coord. in cylindrical coordinates
//  function returns 'true' if the point lies inside LCMS
//
bool C3dMesh::M3DcylIsInside(const Vector3d &cyl)
{
  if(!meshOK) return CStconfig::cylIsInside(cyl);
  bool save = doPositionCheck;
  doPositionCheck = true;
  bool isInside = (M3Dcyl2s(cyl)<sLCMS)?true:false;
  doPositionCheck = save;  // restore
  return  isInside;
}

//***************************************************************
// INPUT:
//  Vector3d xyz -- point coord. in cylindrical coordinates
//  function returns 'true' if the point lies inside LCMS
//
bool C3dMesh::M3DxyzIsInside(const Vector3d &xyz)
{
  if(!meshOK) return CStconfig::xyzIsInside(xyz);
  bool save = doPositionCheck;
  doPositionCheck = true;
  bool isInside = (M3Dxyz2s(xyz)<sLCMS)?true:false;
  doPositionCheck = save;  // restore
  return  isInside;
}


//****obsolete***********************************************************
//   function retuns grad(s) = Vector3d(ds/dr, ds/dfi/r, ds/dz)
//
// not used; obsolete: see M3DgetGrads,
//   it's better to save grad(s) on a mesh, then interpolate.
Vector3d C3dMesh::M3DgetGrads1(const Vector3d &cyl) const
{
  Vector3d gs(0.,0.,0.);
  if(!FIND_COORD(cyl)) return gs; // return zero if point lies outside

  int i,j,k;
  for(i=i0;i<i0+4;i++)     // loop over cube (i0,j0,k0)..(i0+4,j0+4,k0+4)
    for(j=j0;j<j0+4; j++)
      for(k=k0;k<k0+4; k++) {
        gs[0] += mNormFlux[i](j,k)*(cfi [i-i0]*crd[j-j0]*cz [k-k0]); // ds/dr
        gs[1] += mNormFlux[i](j,k)*(cfid[i-i0]*cr [j-j0]*cz [k-k0]); // ds/dfi
        gs[2] += mNormFlux[i](j,k)*(cfi [i-i0]*cr [j-j0]*czd[k-k0]); // ds/dz
      }
  gs[0] *= edR;         // grad(s)_r  == ds/dr
  gs[1] *= edFi/cyl[0]; // grad(s)_fi == ds/dfi/r
  gs[2] *= edZ;         // grad(s)_z  == ds/dz
  if(symmetryActivated) {  //!$!
    gs[1]=-gs[1];
    gs[2]=-gs[2];
  };
  return gs;
}

#if 1

//***************************************************************
// get grads and derivatives dgrads/dr, dgrads/dfi/r, dgrads/dz in point cyl
//  Input:
//    Vector3d cyl -- cylindrical coordinates
//  Output:
//    Vector3d dgradsdr   -- dgrads/dr, where dgradsdr[0]=dgradsr/dr, dgradsdr[1]=dgradsfi/dr, dgradsdr[2]=dgradsz/dr
//    Vector3d dgradsdfir -- dgrads/dfi/r
//    Vector3d dgradsdz   -- dgrads/dz
//  Return:
//    Vector3d grads -- cylindrical coordinates of grads
//
Vector3d C3dMesh::M3Dgetdgradscyl(const Vector3d &cyl,Vector3d &dgradsdr,Vector3d &dgradsdfir,Vector3d &dgradsdz) const
{
  Vector3d grads(0.,0.,0.);
  dgradsdr  = 0; // dgrads/dr
  dgradsdfir= 0; // dgrads/dfi
  dgradsdz  = 0; // dgrads/dz
  if(!FIND_COORD(cyl)) return grads; // return zero if point lies outside

  int i,j,k;
  for(i=0;i<4;i++) {    // loop over cube (i0,j0,k0)..(i0+4,j0+4,k0+4)
    const CArray2d<Vector3d> &G1 = mgrads [i+i0];
    for(j=0;j<4;j++) {
      const Vector3d * G2 = G1[j+j0];
      for(k=0;k<4;k++) {
        grads     += G2[k+k0]*(cfi [i]*cr [j]*cz [k]);
        dgradsdr  += G2[k+k0]*(cfi [i]*crd[j]*cz [k]);
        dgradsdfir+= G2[k+k0]*(cfid[i]*cr [j]*cz [k]);
        dgradsdz  += G2[k+k0]*(cfi [i]*cr [j]*czd[k]);
      }
    }
  }
  dgradsdr  *= edR;         // dgrads/dr
  dgradsdfir*= edFi/cyl[0]; // dgrads/dfi/r
  dgradsdz  *= edZ;         // dgrads/dz
  if(symmetryActivated) {   //!$!
    grads[1] = -grads[1];
    grads[2] = -grads[2];

    dgradsdr[1]   = -dgradsdr[1];
    dgradsdr[2]   = -dgradsdr[2];

    dgradsdfir[0] = -dgradsdfir[0];
    dgradsdz  [0] = -dgradsdz  [0];
  }
  return grads;
}

//***************************************************************
// get gradS and partial derivatives at point xyz
//  Input:
//    Vector3d xyz -- cartesian coordinates
//  Output:
//    Vector3d dGdx -- dgrads/dx, where dgradsdx(1)=dgradsx/dx, dgradsdx(2)=dgradsy/dx, dgradsdx(3)=dgradsz/dx
//    Vector3d dGdy -- dgrads/dy
//    Vector3d dGdz -- dgrads/dz
//  Return:
//    Vector3d grads -- cartesian coordinates of grads
//
Vector3d C3dMesh::M3Dgetdgradsxyz(const Vector3d &xyz,Vector3d &dGdx,Vector3d &dGdy,Vector3d &dGdz) const
{
  Vector3d cyl(sqrt(xyz.x()*xyz.x()+xyz.y()*xyz.y()),
               mod2pi(atan2(xyz.y(),xyz.x())),
               xyz.z());
  Vector3d g, dgcdr,dgcdfir,dgcdz; // grads and derivatives in cyl. coordinates, see M3Dgetdgradscyl 

  dGdx = 0; // dG/dx
  dGdy = 0; // dG/dy
  dGdz = 0; // dG/dz

  g=M3Dgetdgradscyl(cyl,dgcdr,dgcdfir,dgcdz);

  if(g==0) return g; // return zeros if point lies outside

// transform to Cartesian coordinates
  double r  = cyl[0];
  double fi = cyl[1];
  double cs = ::cos(fi);
  double sn = ::sin(fi);
  double Gx = g[0]*cs-g[1]*sn;  // Gx=Gr*cos(fi)-Gfi*sin(fi);
  double Gy = g[0]*sn+g[1]*cs;  // Gy=Gr*sin(fi)+Gfi*cos(fi);
  double Gz = g[2];

  double dGxdr = dgcdr[0]*cs-dgcdr[1]*sn;            // dGxdr = dGr/dr*cs-dGfi/dr*sn;
  double dGydr = dgcdr[0]*sn+dgcdr[1]*cs;            // dGydr = dGr/dr*sn+dGfi/dr*cs;

  double dGxdfir = dgcdfir[0]*cs-dgcdfir[1]*sn-Gy/r;  // dGx/dfi/r = dGr/dfi/r*cs-dGfi/dfi/r*sn-Gy/r;
  double dGydfir = dgcdfir[0]*sn+dgcdfir[1]*cs+Gx/r;  // dGy/dfi/r = dGr/dfi/r*sn+dGfi/dfi/r*cs+Gx/r;

  dGdx.x() = dGxdr*cs-dGxdfir*sn;      // dGx/dx=dGx/dr*cs-dGx/dfi/r*sn
  dGdy.x() = dGxdr*sn+dGxdfir*cs;      // dGx/dy=dGx/dr*sn+dGx/dfi/r*cs
  dGdz.x() = dgcdz[0]*cs-dgcdz[1]*sn;    // dGx/dz=dGr/dz*cs-dGfi/dz*sn;

  dGdx.y() = dGydr*cs-dGydfir*sn;      // dGy/dx=dGy/dr*cs-dGy/dfi/r*sn
  dGdy.y() = dGydr*sn+dGydfir*cs;      // dGy/dy=dGy/dr*sn+dGy/dfi/r*cs
  dGdz.y() = dgcdz[0]*sn+dgcdz[1]*cs;    // dGy/dz=dGr/dz*sn+dGfi/dz*cs;

  dGdx.z() = dgcdr[2]*cs-dgcdfir[2]*sn;  // dGz/dx=dGz/dr*cos(fi)-dGz/dfi/r*sin(fi);
  dGdy.z() = dgcdr[2]*sn+dgcdfir[2]*cs;  // dGz/dy=dGz/dr*sin(fi)+dGz/dfi/r*cos(fi);
  dGdz.z() = dgcdz[2];                  // dGz/dz

  return Vector3d(Gx,Gy,Gz);
}


#endif

};   //namespace MConf
