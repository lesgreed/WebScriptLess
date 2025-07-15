#include "../include/CRayTrace.h"
#include <iostream>
#include <time.h>

namespace MConf {

//****************************************************************************
CRayTrace::CRayTrace()
{
  sizeTrace = 0;
  maxTrace  = 0;
}

//****************************************************************************
// clear trace arrays
void CRayTrace::clearTrace()
{
  CRayTrace::freeTrace();

}

//****************************************************************************
// clear all arrays
void CRayTrace::clear()
{
  CRayTrace::freeTrace();
  CStconfig::freeConf();
}

//****************************************************************************
void CRayTrace::freeTrace()
{
  svertex.clear();

  len.clear();
  reffc.clear();
  vol.clear();
  Rxyz.clear();
  snumber.clear();
  maxTrace = 0;
  sizeTrace = 0;
}

//****************************************************************************
bool CRayTrace::resizeTraceList()
{
  maxTrace = mNsurf*4+20;
  sizeTrace = 0;

  Rxyz.init(maxTrace);
  len.init(maxTrace);
  reffc.init(maxTrace);
  vol.init(maxTrace);
  snumber.init(maxTrace);

  return true;
}

//*********************************************************************
void CRayTrace::addInList(int is, double s, double length, Vector3d &ri)
{
  if(sizeTrace>=maxTrace) return;       // there is no space
//?  for(int i=sizeTrace-1; i>=0; i--)
//?    if(fabs(len[i]-length)<1.e-6) return;  // return if in list
  len.rw() [sizeTrace]   = length;           // otherwise add
  Rxyz.rw()[sizeTrace]   = ri;
  reffc.rw()[sizeTrace]  = r(s);
  vol.rw() [sizeTrace]   = s;   // not a volume! see sort
  snumber.rw()[sizeTrace] = is;
  sizeTrace++;
}

//****************************************************************************
bool CRayTrace::buildSurface(double s, int nTheta, int nPhi, double Phi1, double Phi2)
{
  if(s==0) return false;
  static double shift = 0.01/pi;
  bool wrapPhi = ((Phi1==0.)&&(Phi2==0.))?true:false;
  double dPhi =  wrapPhi?twopi:(Phi2-Phi1);
  int nPhiOnPeriod=0; // number of points on period
  if(wrapPhi&&mNp>1) {
    nPhiOnPeriod = nPhi/mNp; // mNp is # of periods
    if(nPhiOnPeriod<20) nPhiOnPeriod=20;
    nPhi = nPhiOnPeriod*mNp+1;
  }
  if(!svertex.resize(nTheta,nPhi)) return false;
  dPhi /= (nPhi-1);
  double dtheta = twopi/(nTheta-1);
  int i,j,k;
  if(wrapPhi) {  // here we use periodicity
    for(i=0; i<nTheta-1; i++) {          // poloidal direction
      double theta = i*dtheta + shift;
      for(j=0; j<nPhiOnPeriod; j++) {    // toroidal direction
        double phi = j*dPhi + shift;
        svertex(i,j)=mag2cyl(s,theta,phi); // point (x,y,z) on surface 's'
        for(k=1; k<mNp; k++) {
          svertex(i,j+k*nPhiOnPeriod)[0]=svertex(i,j)[0];
          svertex(i,j+k*nPhiOnPeriod)[1]=svertex(i,j)[1]+k*mPeriod;
          svertex(i,j+k*nPhiOnPeriod)[2]=svertex(i,j)[2];
        }
      }
    }
    for(j=0; j<nPhi-1; j++) svertex(nTheta-1,j) = svertex(0,j);
    for(i=0; i<nTheta; i++) { // cyl-->xyz transformation
      svertex(i,nPhi-1) = svertex(i,0);
      for(j=0; j<nPhi; j++) {
        double r = svertex(i,j)[0];
        double f = svertex(i,j)[1];
        svertex(i,j)[0] = r*::cos(f);
        svertex(i,j)[1] = r*::sin(f);
      }
    }
  }
  else {
    for(i=0; i<nTheta-1; i++) {
      double theta = i*dtheta + shift;
      for(j=0; j<nPhi; j++) {
        double phi = j*dPhi+Phi1;
        svertex(i,j)=mag2xyz(s,theta,phi); // point (x,y,z) on surface 's'
      }
    }
    for(j=0; j<nPhi; j++) svertex(nTheta-1,j) = svertex(0,j);
  }
  return true;
}

//*********************************************************************
// Find intersections with the surface #is
//
void CRayTrace::traceSurface(int is, double s, const CArray2d<Vector3d> &vertex, int step)
{
  if(s==0) return;
  if(!vertex.isOK()) return;
  if(maxTrace==0) resizeTraceList();
  if(maxTrace==0) return;               // return if no memory
  int i,i1, ni=vertex.size1()-1;
  int j,j1, nj=vertex.size2()-1;
  for(i=0; i<ni; i+=step) {
    i1=i+step; if(i1>ni) i1=ni;
    for(j=0; j<nj; j+=step) {
      Vector3d r[4],ri;  // vertexes of a rectangle, intersection point
      j1=j+step; if(j1>nj) j1=nj;
      r[0] = vertex(i, j );
      r[1] = vertex(i1,j );
      r[2] = vertex(i1,j1);
      r[3] = vertex(i, j1);
      // triangle 012
      double length = isIntersection(r,ri); // check for intersection with triangle (r0,r1,r2)
      if(length>=0) addInList(is,s,length,ri);  // if true then there is an intersection with the triangle
// printf("123 i=%02d j=%02d s=%g  ri=(%g %g %g)\n",i,j,s,ri.x(),ri.y(),ri.z());
      r[1] = r[3]; // triangle 032
      length = isIntersection(r,ri);
      if(length>=0) addInList(is,s,length,ri);
    }
  }
}

//*********************************************************************
// Test intersection of the ray with the triangle
// INPUT
//  r -- vertices of triangle
//  the ray is R(t)=rayR0+rayRd*length  with length>0, see setRay(...
// OUTPUT
//  ri -- intersection coordinates (cartesian)
//  function returns distance from origin of the ray to the intersection
double CRayTrace::isIntersection(const Vector3d *r, Vector3d &ri) const
{
  Vector3d n = (r[1]-r[0])^(r[2]-r[0]);     // normal to the triangle
  ri = 0;
  double vd = n*rayRd;             // dot product
  if(fabs(vd)<1e-13) return -1.;   // the ray parallels the triangle
  double length = (n*(r[0]-rayR0))/vd;  //
  if(length<0) return -1.;
  ri = rayR0+rayRd*length; // there is an intersection with the plane of the triangle
  return isInsideTriangle(r,n,ri)>0?length:-1.;
}

//*******************************************************************
// Routine to sort intersections into increasing order of distances.
// For simplicity we use a bubble sort - ok for sizeTrace small.
void CRayTrace::sort()
{
  bool swp=true;
  if(sizeTrace==0) return;
  int i,i1,Top;
  for(Top=sizeTrace-1; Top>0&&swp; Top--) {
    swp=false;
    for(i=0,i1=1; i<Top; i++,i1++)
      if(len[i] > len[i1]) {
        Swap (len.rw() [i], len.rw() [i1]);
        Swap (reffc.rw()[i],reffc.rw()[i1]);
        Swap (vol.rw() [i], vol.rw() [i1]);
        Swap (Rxyz.rw()[i], Rxyz.rw()[i1]);
        swp = true;
      }
  }
  for(i=1; i<sizeTrace; i++) len.rw()[i] -= len[0];
  len.rw()[0]=0; // we measure distances starting from the first intersection

  for(i=1; i<sizeTrace; i++)  // remove 2nd intersections with plasma
    if(snumber[i]==snumber[0]) {
      sizeTrace = i+1;
      break;
    }
  for(i=0; i<sizeTrace; i++) { //set volume (see addInList)
    double s = vol[i];
    vol.rw()[i]   = Volume(s);
  }
return;

//@
  if(sizeTrace%2==0) {
    int sz2=sizeTrace/2;  // sizeTrace must be even number!!
    for(i=0; i<sz2; i++) {
      snumber.rw()[i    ]=sz2-i;
      snumber.rw()[i+sz2]=i+1;
    }
  }
  else {   // sizeTrace is odd, it is possible when ray only touch surface
    int sz2=(sizeTrace-1)/2;
    for(i=0; i<sz2; i++) {
      snumber.rw()[i    ]=sz2-i;
      snumber.rw()[i+sz2+1]=i+1;
    }
    snumber.rw()[sz2]=0;
  }
}

//*********************************************************************
// Find intersections of the ray with all surfaces
// The ray is R(t)=r0+rd*t  with t>0
// It is assumed, that the origin of ray lies outside of the last surface
//
// INPUT:
//   r0, rd  -- origin and direction of the ray in cartesian coordinates
//   dtheta, dphi  -- poloidal and totoidal angle mesh steps
//   nSurf  -- number of surfaces to trace
// OUTPUT:
//   use the following methods to retrieve trace information:
//     int    getTraceSize()
//     Vector3d *getTraceRxyz()
//     double *getTraceLen()
//     double *getTracereff()
//     double *getTraceVolume()
//     int    *getTraceSurfNum()
//
void CRayTrace::traceAll(const Vector3d &r0,const Vector3d &rd, double dtheta, double dphi,int nSurf)
{
  CRayTrace::freeTrace();
  if(xyzIsInside(r0)) return; // return if inside plasma
  setRay(r0,rd);  // the ray is R(t)=R0+Rd*t  with t>0
  int nTheta, nPhi,i;
  double phi0, phi1, dr;
  if(!findSegment(dtheta,dphi,nTheta,nPhi,phi0,phi1)) return;

  mNsurf = nSurf; // number of surfaces to trace, see resizeTraceList()
  dr=1.0/(nSurf-1);
  for(i=0; i<nSurf; i++) {
    double s = (i*dr)*(i*dr);
    if(s==0) s=dr*dr/9;
    if(i==nSurf-1) s=1;
    buildSurface(s, nTheta, nPhi,phi0,phi1);
    traceSurface(i,s,svertex);
  }
  sort();
return;

  dr=1.0/nSurf;
  for(i=1; i<=nSurf; i++) {
    double s = (i*dr)*(i*dr);
    if(i==nSurf) s=1;
    buildSurface(s, nTheta, nPhi,phi0,phi1);
    traceSurface(i,s,svertex);
  }
  sort();
}

//*********************************************************************
int CRayTrace::traceLCMS(const Vector3d &r1,const Vector3d &r2, Vector3d &entrypoint, Vector3d &exitpoint, double dpol,double dtor)
{
//1  clock_t tStart = clock();  //timing
  entrypoint=0;
  exitpoint =0;
  if(xyzIsInside(r1)) return -2; // if true then point lies inside LCMS

  Vector3d rd = r2-r1;         // direction of ray
  traceLast(r1, rd, dpol*degree, dtor*degree);
  if(sizeTrace==0) return -1;  // ray does not enter in plasma

  entrypoint = Rxyz[0];      // point where ray enters into plasma
  exitpoint  = Rxyz[1];      // point where ray exits  plasma
  double r21 = (r2-r1).abs2();
  double e1 = (entrypoint-r1).abs2();
  double x1 = (exitpoint-r1).abs2();
//  CRayTrace::freeTrace(); // ?
//1  std::cerr<< double(clock() - tStart)/CLOCKS_PER_SEC<<" traceLCMS"<<std::endl;
  if(e1> r21) return 0;  // segment r1--r2 is outside plasma
  if(x1> r21) return 1;  // only entry point lies on the segment r1--r2, i.e. r2 is inside plasma
  if(x1<=r21) return 2;  // entry and exit points lie on the segment r1--r2
  return 0;

}

//*********************************************************************
// Find intersection of the Ray with the last surface
// The ray is R(t)=r0+rd*t  with t>0
// It is assumed, that the origin of ray lies outside of the last surface
//
// INPUT:
//   r0, rd  -- origin and direction of the ray in cartesian coordinates
//   dtheta, dphi  -- poloidal and totoidal angle mesh steps
// OUTPUT:
//   use the following methods to retrieve trace information:
//     int    getTraceSize()
//     Vector3d *getTraceRxyz()
//     double *getTraceLen()
//     double *getTracereff()
//     double *getTraceVolume()
//     int    *getTraceSurfNum()
//
void CRayTrace::traceLast(const Vector3d &r0,const Vector3d &rd, double dtheta, double dphi)
{
  CRayTrace::freeTrace();
  if(xyzIsInside(r0)) return; // return if inside plasma
  setRay(r0,rd);  // the ray is R(t)=R0+Rd*t  with t>0
  int nTheta, nPhi;
  double phi0, phi1;
  if(!findSegment(dtheta,dphi,nTheta,nPhi,phi0,phi1)) return;

  double s = 1;     // last surface
  mNsurf = 1; // number of surfaces to trace, see resizeTraceList()
  buildSurface(s, nTheta, nPhi,phi0,phi1); //fill  svertex //svertex.print();
  traceSurface(1,s,svertex);
  sort();
}

//*********************************************************************
// Find toroidal angles of intersections with the last surface
bool CRayTrace::findSegment(double dtheta, double dphi, int &nTheta, int &nPhi, double &phi0, double &phi1)
{
  double s=1;
  mNsurf = 1; // number of surfaces to trace, see resizeTraceList()
//1  clock_t tStart = clock();  //timing
  const CArray2d<Vector3d> & vertex = *LCMS();

  traceSurface(1,s,vertex,2);
  sort();
  if(sizeTrace==0) return false;

//2  phi0 = mod2pi( atan2(Rxyz[0].y(), Rxyz[0].x()) );
//2  phi1 = mod2pi( atan2(Rxyz[1].y(), Rxyz[1].x()) );

  Vector3d magCrd0 = xyz2mag(Rxyz[0]);  // get magnetic coordinates
  Vector3d magCrd1 = xyz2mag(Rxyz[1]);
  phi0 = magCrd0[2];
  phi1 = magCrd1[2];
// sort angles
///  phi0 = mod2pi(magCrd0[2]);     // mod2pi( added 09.04.05
///  phi1 = mod2pi(magCrd1[2]);     // mod2pi( added 09.04.05
///  if(phi1<phi0) Swap(phi1,phi0);
///  if(phi1-phi0>pi) {phi0+=twopi;Swap(phi1,phi0);}
  sortAngles(phi0,phi1);      //added 09.04.05
  phi0-=0.02; // add some margin 0,0175rad = 1degree
  phi1+=0.02;
  dtheta = (dtheta==0.)?degree:dtheta;
  dphi   = (dphi==0.)  ?degree:dphi;
  nTheta = int(twopi/dtheta+0.5);
  nPhi   = int(twopi/dphi+0.5);
  dphi = twopi/nPhi;
  nPhi = int(fabs(phi1-phi0)/dphi)+1;
  CRayTrace::freeTrace();
//1  std::cerr<< double(clock() - tStart)/CLOCKS_PER_SEC<<" findSegment"<<std::endl;
  return true;
}

};   //namespace MConf
