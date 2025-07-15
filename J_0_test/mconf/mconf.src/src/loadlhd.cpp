#include "../include/CStconfig.h"
#include <fstream>
#include <iostream>
#include <time.h>

namespace MConf {

#ifdef _WIN32
  #define  Strnicmp   _strnicmp
#else
  #define  Strnicmp   strncasecmp
#endif

//****************************************************************************
// replace D with E in scientific format notation: 1.2D-02 --> 1.2E-02
// return false if not a number
static bool translate(char *p)
{
  if(0!=strchr(p,'N')) return false;  // search for a 'NaN', 'INF'
  while(p) {
    if((p=strchr(p,'D'))!=0)  { *p='E'; p++;}
    else break;
  }
  return true;
}

//****************************************************************************
// Read bc-file of LHD-format (produced by NEWBOZ)
//   here fscanf is used because it works faster than C++ streams >>
bool CStconfig::loadlhd(const char * fname)
{
  FILE * fp = fopen(fname,"r");
  if(fp==NULL) return(mOK=false);
  mfname = fname;

  char buff[512];

  int nC=0;
  while(1) {
    if(fgets(buff,sizeof(buff)-1,fp)==NULL) {fclose(fp); return(mOK=false);}
    int len = (int)strlen(buff);
    if(strncmp(buff,"CC",mmin(len,2))!=0&&strncmp(buff,"cc",mmin(len,2))!=0) break;
    mComments.push_back(buff);
    nC++;
  }
  nC++;

  int i,ih,ns,m,n, numOfMod;

  translate(buff);
//bool needR0 = (6!=sscanf(buff,"%d %d %d %d %lf %lf",&numOfMod,&ns,&i,&mNp,&ma,&mR0));
//bool needR0 = (6!=sscanf(buff,"%d %d %d %d %lf %lf",&numOfMod,&ns,&mNp,&i,&ma,&mR0));

// try to find out the format of a file:
//   qas2-fb.bc contains six numbers in the first line:  578  60  16  2 0.4133 1.471
//   lhd_boz10.r360q100b004a8020.FMT  contains only three:  830  60  10
  int nRead = sscanf(buff,"%d %d %d %d %lf %lf",&numOfMod,&ns,&i,&mNp,&ma,&mR0);
  if(nRead == 3) mNp = i;

  bool needR0 = (6!=nRead);

{ // find numbers of mods and resize arrays
  fgets(buff,sizeof(buff)-1,fp);
  if(!translate(buff)) {
    std::cerr <<"Bad number in the file "<<fname<<std::endl;
    fclose(fp); 
    return(mOK=false);
  } 
  double s,iot,itor,ipol;
  sscanf(buff,"%lf %lf %lf %lf",&s,&iot,&itor,&ipol);
  if(s==0) m1stSindex = 0;  // file has magnetic axis if true
  for(i=1; i<ns; i++) fgets(buff,sizeof(buff)-1,fp); 
  m=n=0;
  for(ih=0; ih<numOfMod; ih++) {
    int mm,nn;
    fscanf(fp,"%d %d",&mm,&nn);
    fgets(buff,sizeof(buff)-1,fp); // ignore rest of line
    nn = abs(nn)/mNp;
    m=mmax(mm,m);
    n=mmax(nn,n);
    for(i=0; i<ns; i++) fgets(buff,sizeof(buff)-1,fp); // for(i=0; i<ns; i++) fscanf(fp,"%lf %lf %lf %lf",&dum,&dum,&dum,&dum);
  }
  if(!resize(ns,m,n)) {fclose(fp); return(mOK=false);}
  rewind(fp);
  while(nC--) fgets(buff,sizeof(buff)-1,fp); // skip header lines
}


  for(i=0; i<mNs0; i++) {
    int is = m1stSindex+i;
    fgets(buff,sizeof(buff)-1,fp);
    if(!translate(buff)) {
      std::cerr <<"Bad number in the file "<<fname<<std::endl;
      fclose(fp); 
      return(mOK=false);
    } 
    double s,iot,itor,ipol;
    sscanf(buff,"%lf %lf %lf %lf",&s,&iot,&itor,&ipol);
    ms.rw()[is] =s;
    miota.rw()[is] = iot;
    mItor.rw()[is] = itor*twopi/mu0;
    mIpol.rw()[is] = ipol*twopi/(mNp*mu0); //convert to w7x-format, poloidal current on a period
    mpprime.rw()[is]=0;               //undefined
    mg00.rw()[is]=1;                  //will be calculated later
  }

  mFlux = ms[mLCMSindex];
  for(i=0; i<mNs0; i++)   
    ms.rw()[m1stSindex+i]/=mFlux;
  mFlux *= twopi;


// LHD NEWBOZ produces boozer representation:
//  B = sum Bmn*cos(n*phi - m*theta)
//  R = sum Rmn*cos(n*phi - m*theta)
//  Z = sum Zmn*sin(n*phi - m*theta)
//  fi = phi + sum Phimn*sin(n*phi - m*theta)
//   where fi is the cylindrical angle,
//         phi is the boozer toroidal angle,
//         theta is the boozer poloidal angle.
// see in Plasma and Fusion Research - Vol.5, 016 (2010).
//
// W7-X format is 
//  R = sum Rmn*cos(m*theta-n*Np*phi)
//  Z = sum Zmn*sin(m*theta-n*Np*phi)
//  fi = phi - 2pi/Np*sum Phimn*sin(m*theta-n*Np*phi)
  for(ih=0; ih<numOfMod; ih++) {
    fscanf(fp,"%d %d",&m,&n);
    fgets(buff,sizeof(buff)-1,fp); // ignore rest of line
    n /= mNp;
    for(i=0; i<mNs0; i++) {
      int is = m1stSindex+i;
      fgets(buff,sizeof(buff)-1,fp);
      if(!translate(buff)) {
        std::cerr <<"Bad number in the file "<<fname<<std::endl;
        fclose(fp); 
        return(mOK=false);
      } 
      double b,r,z,p; 
      sscanf(buff,"%lf %lf %lf %lf",&b,&r,&z,&p);
      bmn(is,m,n) = b;
      Rmn(is,m,n) = r;
      Zmn(is,m,n) = z;
      Phimn(is,m,n) = p;
      Zmn  (is,m,n) = -Zmn(is,m,n); //convert to w7x-format
      Phimn(is,m,n)*= mNp/twopi;    //convert to w7x-format
    }
  }

  fclose(fp);
  if(needR0) {
    ma  = 0;
    mR0 = 0;
  }
  vprimeNeeded = true; // signal to load-function: we need vprime(mg00)
  mNewFormat = true;
  return(mOK=true);
}

};   //namespace MConf
