#include "../include/CStconfig.h"
#include <fstream>
#include <iostream>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <math.h>

#ifndef _WIN32
#include <unistd.h>
#include <sys/wait.h>
#else
#include <process.h>
#endif

#include <sys/stat.h>
#include <errno.h>

namespace MConf {

using std::ifstream;
using std::ofstream;
using std::string;
using std::endl;
using std::getline;

#ifdef _WIN32
  #ifndef M_PI
  #define M_PI    3.14159265358979323846
  #endif
#endif

bool CStconfig::loadwout(const char * fname) {return false;}

#define TRANSFORM_TO_BOOZER 1
#define SkipLine1(fp) fgets(sdummy,buff.size(),fp);
#define SkipLineX(fp) fscanf(fp,"%*[^\n]\n",NULL)  // need testing

//********************************************************************
//*********************************************************
// load VMEC output file (wout-file):
//  read the fourier coefficients from the wout-file:
//          Rmn  (j, m, n) = rmnc;    //  on full mesh
//          Zmn  (j, m, n) = zmns;    //  on full mesh
//          Phimn(j, m, n) = 0;                // on half-mesh
//          bmn  (j, m, n) = bmn1;             // |B|   on half-mesh
//          gmn  (j, m, n) = gmnc;             // gmnc   on half-mesh
//          lmn  (j, m, n) = lmns;             // lmns   on half-mesh
//          bContraTheta_mn(j, m, n) = bsupu;  // contravariant poloidal component of B on half-mesh
//          bContraFi_mn   (j, m, n) = bsupv;  // contravariant toroidal component of B on half-mesh
// Note:
//     fscanf is 2 times faster than >>
bool CStconfig::loadvmec(const char * fname, long fPostn)
{
  mOK = false;
  FILE * fp = fopen(fname,"rb");
  if(fp==NULL) return(mOK=false);
  mfname = fname;

  fseek(fp,fPostn,SEEK_SET); //if(fPostn!=0)  fsetpos(fp,&fPostn);           // set position to read
 
  CArray1d<char> buff(8192);
  char *sdummy=buff.array();
  
  int iversion;
  {
    if(fPostn) {
      //for(int c=fgetc(fp);(c!='\n'&&c!=EOF);c=fgetc(fp)) {;} 
      fgets(sdummy,buff.size(),fp);
      if(strncmp(sdummy,"VMEC VERSION",12)!=0) {
        while(1) {
          if(fgets(sdummy,buff.size(),fp)==NULL) return(mOK=false);
          int len = (int)strlen(sdummy);
          if(strncmp(sdummy,"CC",mmin(len,2))!=0&&strncmp(sdummy,"cc",mmin(len,2))!=0) break;
        }
      }
    } else {
      fgets(sdummy,buff.size(),fp);
    }
    if(strncmp(sdummy,"VMEC VERSION",12)!=0) return (mOK=false);
    double vers = atof(sdummy+15);
    iversion = int(vers*10000+0.1);
    if(iversion < 62000) return (mOK=false);
  }
  double  wb, wp, gamma, pfac, rmax_surf, rmin_surf, zmax_surf;
  fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&wb,&wp,&gamma,&pfac,&rmax_surf,&rmin_surf,&zmax_surf);

  int nfp, ns, mpol, ntor,mnmax,mnmax_nyq(0), itfsq, niter, iasym, ireconstruct, ier_flag;
  fscanf(fp,"%d %d %d %d %d ",&nfp,&ns,&mpol,&ntor,&mnmax);
  if(iversion > 80000) fscanf(fp,"%d", &mnmax_nyq);
  fscanf(fp,"%d %d %d %d %d",&itfsq,&niter,&iasym,&ireconstruct,&ier_flag);

  vmecErrorCode = ier_flag;
  if(iasym) { 
    std::cerr<< "VMEC import: asymmetry mode"<<std::endl;
    vmecErrorCode = -3000;
    std::cerr<< "VMEC import: at the moment MConf can not handle wout-file in asymmetry mode"<<std::endl;
    return false;    // at the moment can not handle asymmetry
  }
  
  if(ier_flag!=0) {  // check VMEC's error flag
    switch(ier_flag) {
      case 1:
        std::cerr<< "VMEC file '"<<fname<<"' contains error flag: "<<ier_flag<<"(bad initial jacobian)"<<std::endl; 
        break;
      case 4:
//        std::cerr<< "VMEC error flag: "<<ier_flag<<".  more iterations to reach force limit is needed"<<std::endl; 
        break;
      case 6:
        std::cerr<< "VMEC error flag: "<<ier_flag<<".  bad jacobian"<<std::endl; 
        break;
      default:
        std::cerr<< "VMEC error flag: "<<ier_flag<<".  unknown error"<<std::endl; 
        break;
    }
    if(!(ier_flag==4))  return false; 
  }
  
  vmecCoordinates = true; // resize needs this flag
  vmecAsym = iasym != 0;  // resize needs this flag
  // m is the poloidal mode number, 0<=m<=mpol-1.
  // n is the toroidal mode number. -ntor<=n<=ntor.
  // if(m==0)    0<=n<=ntor.
  // The magnetic axis is included in VMEC output. - 
  // However, we will discard it here (m1stSindex must be >= 1), see (ns-1)
  if (!resize(ns-1 , mpol-1, ntor))  return vmecCoordinates=false;
  mNewFormat = true;
  mNp=nfp;

  int imse2_minus_1,itse, nbsets, nobd, nextcur, nstore_seq;
  fscanf(fp,"%d %d %d %d %d %d",&imse2_minus_1,&itse,&nbsets,&nobd,&nextcur,&nstore_seq);
  if(nbsets>0) {
    double ddummy;
    for (int i = 0; i < nbsets; i++) fscanf(fp,"%lf", &ddummy);
  }

  SkipLine1(fp);  // ignore rest of line
  SkipLine1(fp);  // ignore mgrid_file name

  int *mnum = new int [mnmax];
  int *nnum = new int [mnmax];
  int *mnum_nyq = new int [mnmax_nyq];
  int *nnum_nyq = new int [mnmax_nyq];
  
  if(iversion > 80000) {
      for(int i = 0; i < ns; i++) {
        int j = i-1+ m1stSindex;
        if(iversion >= 90000) SkipLine1(fp);  // ignore line   'JS: ' n
        for(int mn = 0; mn < mnmax; mn++) { {
        ////for(m = 0; m < mpol; m++) {
        ////  int nmin = (m==0)?0:-ntor;
        ////  for(n = nmin; n <=ntor; n++) {
            if (i==0)   fscanf(fp,"%d %d",&mnum[mn],&nnum[mn]); // read  m, n*nfp into arrays
            int m = mnum[mn];
            int n = nnum[mn]/nfp;
            double rmnc, zmns;    // are defined on full-mesh
            double lmns;          // are defined on half-mesh
            double rmns,zmnc,lmn_c;
            fscanf(fp,"%lf %lf %lf",&rmnc,&zmns,&lmns);
            if(iasym) { // additional fourier coefficients in asymmetrical mode 
              fscanf(fp,"%lf %lf %lf",&rmns,&zmnc,&lmn_c);
            }
            if(i==0) continue; // skip magnetic axis if i==0
            Rmn(j, m, n) = rmnc;    //  on full mesh
            Zmn(j, m, n) = zmns;    //  on full mesh
            lmn(j, m, n) = lmns;    //  on half-mesh, sin(mu+nv) stellarator symmetric part
            if(iasym) { // additional fourier coefficients in asymmetrical mode 
              Rmns(j, m, n) = rmns;    //  on full mesh, asymmetric part
              Zmnc(j, m, n) = zmnc;    //  on full mesh, asymmetric part
              lmnc(j, m, n) = lmn_c;   //  on half-mesh, cos(mu+nv) asymmetric part
            }
          }
        }

        for(int mn = 0; mn < mnmax_nyq; mn++) { {
        ////for (m = 0; m < mpol; m++) { // could be a problem if vmec2000 was run with lnyquist==.true. 
        ////  int nmin = (m==0)?0:-ntor;
        ////  for (n = nmin; n <=ntor; n++) {
            if (i==0) fscanf(fp,"%d %d",&mnum_nyq[mn],&nnum_nyq[mn]); // read  m, n*nfp into arrays
            double bmnc,gmnc,bsubu,bsubv,bsubs,bsupu,bsupv;  // are defined on half-mesh
            double bmn_s,gmn_s,bsubumns,bsubvmns,bsubsmnc,bsupumns,bsupvmns;
            fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&bmnc,&gmnc,&bsubu,&bsubv,&bsubs,&bsupu,&bsupv);
            if(iasym) { // additional fourier coefficients in asymmetrical mode 
              fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&bmn_s,&gmn_s,&bsubumns,&bsubvmns,&bsubsmnc,&bsupumns,&bsupvmns);
            }
            if(i==0) continue;     // skip magnetic axis if i==0
            // m is the poloidal mode number, 0<=m<=mpol-1.
            // n is the toroidal mode number. -ntor<=n<=ntor.
            // if(m==0)    0<=n<=ntor.
            int m = mnum_nyq[mn];
            int n = nnum_nyq[mn]/nfp;
            if(m>mpol-1) continue; // skip extra modes 
            if(n>ntor)  continue;  // skip extra modes 
            if(n<-ntor) continue;  // skip extra modes
            {  // symmetrical mode
              bmn  (j, m, n) = bmnc;             // |B|   on half-mesh
              gmn  (j, m, n) = gmnc;             // gmnc   on half-mesh
              bContraTheta_mn(j, m, n) = bsupu;  // B^u(mn)c contravariant poloidal component of B on half-mesh, cos(mu+nv) stellarator symmetric part
              bContraFi_mn   (j, m, n) = bsupv;  // B^v(mn)c contravariant toroidal component of B on half-mesh,  cos(mu+nv) stellarator symmetric part
              //bCovaTheta_mn(j, m, n) = bsubu;  // B_u(mn)c covariant poloidal component of B on half-mesh, cos(mu+nv) stellarator symmetric part
              //bCovaFi_mn   (j, m, n) = bsubv;  // B_v(mn)c covariant toroidal component of B on half-mesh, cos(mu+nv) stellarator symmetric part
              Phimn(j, m, n) = 0;
#if TRANSFORM_TO_BOOZER  
              double wrk;
              if(m!=0) wrk = bsubu/m;
              else if(n!=0) wrk = bsubv/(-n*mNp);
              else wrk = 0;
              Phimn(j, m, n) = wrk;
#endif
            }
            if(iasym) { // additional fourier coefficients in asymmetrical mode 
              bmns  (j, m, n) = bmn_s;   // half mesh: bmns sin(mu+nv)            asymmetric part
              gmns  (j, m, n) = gmn_s;   // gmns   on half-mesh
              bContraTheta_mns(j, m, n) = bsupumns;  // B^u(mn)s contravariant poloidal component of B on half-mesh, sin(mu+nv) asymmetric part
              bContraFi_mns   (j, m, n) = bsupvmns;  // B^v(mn)s contravariant toroidal component of B on half-mesh, sin(mu+nv) asymmetric part
              Phimns(j, m, n) = 0;   //sin(), asymmetric part
#if TRANSFORM_TO_BOOZER  
              double wrk;
              if(m!=0) wrk = bsubumns/m;
              else if(n!=0) wrk = bsubvmns/(-n*mNp);
              else wrk = 0;
              Phimns(j, m, n) = wrk;
#endif              
            }
          }
        }
        if(iversion >= 90000) SkipLine1(fp);  // skip rest of line
      }  // end loop ns
  }

  else  {  //for VMECversion <= 8.00

      for(int i = 0; i < ns; i++) {
        int j = i-1+ m1stSindex;
        for(int m = 0; m < mpol; m++) {
          int nmin = (m==0)?0:-ntor;
          for(int n = nmin; n <=ntor; n++) {
            int idummy;
            if (i==0)   fscanf(fp,"%d %d",&idummy,&idummy); // read  m, n*nfp
            double rmnc, zmns,currv;                              // are defined on full-mesh
            double lmns,bmnc,gmnc,bsubu,bsubv,bsubs,bsupu,bsupv;  // are defined on half-mesh
            fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                      &rmnc,&zmns,&lmns,&bmnc,&gmnc,&bsubu,&bsubv,&bsubs,&bsupu,&bsupv,&currv);
            if(i>0) {  // skip magnetic axis if i==0
              Rmn  (j, m, n) = rmnc;    //  on full mesh
              Zmn  (j, m, n) = zmns;    //  on full mesh
              bmn  (j, m, n) = bmnc;             // |B|   on half-mesh
              gmn  (j, m, n) = gmnc;             // gmnc   on half-mesh
              lmn  (j, m, n) = lmns;             // lmns   on half-mesh, sin(mu+nv) stellarator symmetric part
              bContraTheta_mn(j, m, n) = bsupu;  // B^u(mn)c contravariant poloidal component of B on half-mesh, cos(mu+nv) stellarator symmetric part
              bContraFi_mn   (j, m, n) = bsupv;  // B^v(mn)c contravariant toroidal component of B on half-mesh,  cos(mu+nv) stellarator symmetric part
              //bCovaTheta_mn(j, m, n) = bsubu;  // B_u(mn)c covariant poloidal component of B on half-mesh, cos(mu+nv) stellarator symmetric part
              //bCovaFi_mn   (j, m, n) = bsubv;  // B_v(mn)c covariant toroidal component of B on half-mesh, cos(mu+nv) stellarator symmetric part
              Phimn(j, m, n) = 0;
#if TRANSFORM_TO_BOOZER  
              double wrk;
              if(m!=0) wrk = bsubu/m;
              else if(n!=0) wrk = bsubv/(-n*mNp);
              else wrk = 0;
              Phimn(j, m, n) = wrk;
              //////if(m!=0&&n!=0) { //test: the following must be zero
              //////  wrk = bsubu/m - bsubv/(-n*mNp);
              //////  wrk = bsubu/m - bsubv/(-n*mNp);
              //////}
#endif
            }
            if(iasym) { // additional fourier coefficients in asymmetrical mode 
              double rmns,zmnc;
              double lmn_c, bmn_s,gmn_s,bsubumns,bsubvmns,bsubsmnc,bsupumns,bsupvmns;
              fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &rmns,&zmnc,&lmn_c,&bmn_s,&gmn_s,&bsubumns,&bsubvmns,&bsubsmnc,&bsupumns,&bsupvmns);
              fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&bmn_s,&gmn_s,&bsubumns,&bsubvmns,&bsubsmnc,&bsupumns,&bsupvmns);
              if(i>0) {  // skip magnetic axis if i==0
                Rmns(j, m, n) = rmns;    //  on full mesh, asymmetric part
                Zmnc(j, m, n) = zmnc;    //  on full mesh, asymmetric part
                bmns  (j, m, n) = bmn_s;   // half mesh: bmns sin(mu+nv)            asymmetric part
                gmns  (j, m, n) = gmn_s;   // gmns   on half-mesh
                lmnc(j, m, n) = lmn_c;    //  on half-mesh, cos(mu+nv) asymmetric part
                bContraTheta_mns(j, m, n) = bsupumns;  // B^u(mn)s contravariant poloidal component of B on half-mesh, sin(mu+nv) asymmetric part
                bContraFi_mns   (j, m, n) = bsupvmns;  // B^v(mn)s contravariant toroidal component of B on half-mesh, sin(mu+nv) asymmetric part
                Phimns(j, m, n) = 0;   //sin(), asymmetric part
#if TRANSFORM_TO_BOOZER  
                double wrk;
                if(m!=0) wrk = bsubumns/m;
                else if(n!=0) wrk = bsubvmns/(-n*mNp);
                else wrk = 0;
                Phimns(j, m, n) = wrk;
#endif
              } // end if(i>0)
            }  
          }
        }
      }
  }

  delete mnum_nyq;
  delete nnum_nyq;
  delete mnum;
  delete nnum;

  bool standardVMEC = true;      // try to read standard format of wout-file
//bool standardVMEC = false;     // try to read Zarnstorff's format of wout-file
//if(iversion> 80000) standardVMEC = true;  //08.11.2011 we accept only standard format of wout-file
  fpos_t filePosition;
  int iwrk = fgetpos(fp,&filePosition); // save position in file
  for(int k = 0; k < ns; k++) {
    int j = k + m1stSindex;
    msFull.rw()[j] = double(k+1) / double(ns-1);  // VMEC full mesh
    ms.rw()    [j] = (0.5 + k)/double(ns-1);      // VMEC half mesh
  }
tryToReadAnotherFormat:                    // it is possible we return here trying to read Zarnstorff's format of wout-file
  double presPrev=0;
  if(standardVMEC) {
 
//!     HALF-MESH QUANTITIES (except phi, jcuru, jcurv which are FULL MESH - computed in eqfor)
//!     NOTE: jcuru, jcurv are local current densities, NOT integrated in s and normed to twopi
//!     NOTE: In version <= 6.00, mass, press are written out in INTERNAL units
//!     and should be multiplied by 1/mu0 to transform to pascals. In version > 6.00,
//!     the pressure, mass are in correct (physical) units
//!
//!     NOTE: phipf has a twopi * signgs factor compared with phips...
    //WRITE (nwout2, *) (iotaf(js), presf(js)/mu0,twopi*signgs*phipf(js),phi(js), jcuru(js)/mu0, jcurv(js)/mu0, js=1,ns)
    //WRITE (nwout2, *) (iotas(js), mass(js)/mu0, pres(js)/mu0,beta_vol(js), phip(js), buco(js), bvco(js), vp(js),overr(js), specw(js),js=2,ns)

    for(int i = 0; i < ns; i++) { // read the radial profiles
      int j = i + m1stSindex;
      double  iotaf, presf, phipf, phi, jcuru, jcurv;
      fscanf(fp,"%lf %lf %lf %lf %lf %lf", &iotaf,&presf,&phipf,&phi,&jcuru,&jcurv);
      mFlux = phi;
      if(i) mpprime.rw()[j-1] = (presf - presPrev) /(msFull[j]-msFull[j-1]);
      if(i>=ns-2)
        mpprime.rw()[j] = extrapolateLCMS(msFull[j],msFull.constArray(),mpprime.constArray(),j);
      presPrev = presf;
    }
  }
  /* read the radial profiles */
  for(int i = 0; i < (ns-1); i++) {
    int j = i + m1stSindex;
//half-mesh quantities (except phi, jcuru, jcurv which are full mesh)
    double iotas, mass, pressure, beta_vol, phip, buco, bvco;
    fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&iotas,&mass,&pressure,&beta_vol,&phip,&buco,&bvco);    
    miota.rw()[j] = iotas;
    mpres.rw()[j] = pressure;
    mIpol.rw()[j] = bvco*twopi/(mNp*mu0); //convert to w7x-format, poloidal current on a period
    mItor.rw()[j] = buco*twopi/mu0;
    double phi,vp,overr, jcuru, jcurv, specw;
    if(standardVMEC)
      fscanf(fp,"%lf %lf %lf",&vp,&overr,&specw);
    else {
      fscanf(fp,"%lf %lf %lf %lf %lf %lf",&phi,&vp,&overr,&jcuru,&jcurv,&specw);
      mFlux = phi;
      // calculate dp/ds, ?? must be on half mesh
      if(i) mpprime.rw()[j-1] = (pressure - presPrev) /(ms[j]-ms[j-1]);
      if(i==ns-2) mpprime.rw()[j] = extrapolateLCMS(ms[j],ms.constArray(),mpprime.constArray(),j);
      presPrev = pressure;
    }
    mg00.rw()[j] = vp * 4.0 * M_PI * M_PI / mNp; // mg00 is the dV/ds/Np
  }

  double aspect, betatot, betapol, betator, betaxis;
#if 0
  double b0;
  // before 10.04.2015
  fscanf(fp,"%lf %lf %lf %lf %lf %lf",&aspect,&betatot,&betapol,&betator,&betaxis,&b0);
#else
  // changed to handle wout file from USA
  fscanf(fp,"%lf %lf %lf %lf %lf",&aspect,&betatot,&betapol,&betator,&betaxis);
  SkipLine1(fp);    // skip rest of line
#endif
  if (iversion >= 62000) {
    double ddummy;
    fscanf(fp,"%lf",&ddummy); // must be +1 or -1 meaning the sign of Jacobian
    SkipLine1(fp);             // skip rest of line
    fgets(sdummy,buff.size(),fp);    // sdummy is the aaaaaa from a filename in form wout_aaaaaa.txt
#if  0    //10.04.2015    
    bool mustBeTrue =  (abs(ddummy)==1)&&(sdummy[0]!=' ');   // sdummy[0] must not be blank
#else
    bool mustBeTrue =  (abs(ddummy)==1);
#endif

    if(mustBeTrue!=true) {
      if(standardVMEC) {             // if it was standard VMEC try to read another format
        standardVMEC = false;
        fsetpos(fp,&filePosition);   //   restore position
        goto tryToReadAnotherFormat; //   try to read another format
      }
      else {   
        fclose(fp);
        return (mOK=false);
      }
    }
    double IonLarmor, VolAvgB, rbtor0, rbtor, ctor, Aminor_p, Rmajor_p, volume_plasma;
    fscanf(fp,"%lf %lf %lf %lf",&IonLarmor, &VolAvgB, &rbtor0, &rbtor);
    fscanf(fp,"%lf %lf %lf %lf",&ctor,&Aminor_p,&Rmajor_p,&volume_plasma);
    ma  = Aminor_p;
    mR0 = Rmajor_p;
  }

#if TRANSFORM_TO_BOOZER
  {
    for(int i = 0; i < (ns-1); i++) {
      int j = i + m1stSindex;
      for(int m = 0; m < mpol; m++) {
        int nmin = (m==0)?0:-ntor;
        for(int n = nmin; n <=ntor; n++) {
          double wrk = Phimn(j,m,n);
          double iotas = miota[j];
          double bvco  = mIpol[j]/(twopi/(mNp*mu0));
          double buco  = mItor[j]/(twopi/mu0);
          double pmn = (wrk - buco*lmn(j,m,n))/(bvco+iotas*buco);
          Phimn(j,m,n) = pmn;
  //      Phimn(j,m,n) = pmn*mNp/twopi; 
        }
      }  
    }
  }
#endif

  fclose(fp);
  return mOK = true;
}

#ifdef NETCDF
#include <netcdf.h>
#ifdef  NDEBUG
  #undef NDEBUG
  #define NDEBUG_IS_DISABLE
#endif
#include <assert.h>

////////////////////////////////////////////////////////////////////////////////
// Read a scalar integer variable from a CDF file by name (or an array, with  //
//  no bounds checking).                                                      //
////////////////////////////////////////////////////////////////////////////////
void cdf_read_int(int fid, const char *vname, int *vbuf)
{
  int vid;
  assert(vbuf);
  assert(nc_inq_varid(fid, vname, &vid) == NC_NOERR);
  assert(nc_get_var_int(fid, vid, vbuf) == NC_NOERR);
}

////////////////////////////////////////////////////////////////////////////////
// Read a 1D integer array from a CDF file by name, with bounds checking.     //
////////////////////////////////////////////////////////////////////////////////
void cdf_read_int_vec(int fid, const char *vname, int ld, int *vbuf)
{
  size_t lengthp;
  int    dimid, vid, ndims;

  assert(vbuf);
  assert(nc_inq_varid(fid, vname, &vid) == NC_NOERR);
  assert(nc_inq_varndims(fid, vid, &ndims) == NC_NOERR);
  assert(ndims == 1);
  assert(nc_inq_vardimid(fid, vid, &dimid) == NC_NOERR);
  assert(nc_inq_dimlen(fid, dimid, &lengthp) == NC_NOERR);
  assert(static_cast<int>(lengthp) == ld);
  assert(nc_get_var_int(fid, vid, vbuf) == NC_NOERR);
}

////////////////////////////////////////////////////////////////////////////////
// Read a 1D double array from a CDF file by name, with bounds checking.     //
////////////////////////////////////////////////////////////////////////////////
void cdf_read_double_vec(int fid, const char *vname, int ld, double *vbuf)
{
  size_t lengthp;
  int    dimid, vid, ndims;

  assert(vbuf);
  assert(nc_inq_varid(fid, vname, &vid) == NC_NOERR);
  assert(nc_inq_varndims(fid, vid, &ndims) == NC_NOERR);
  assert(ndims == 1);
  assert(nc_inq_vardimid(fid, vid, &dimid) == NC_NOERR);
  assert(nc_inq_dimlen(fid, dimid, &lengthp) == NC_NOERR);
  assert(static_cast<int>(lengthp) == ld);
  assert(nc_get_var_double(fid, vid, vbuf) == NC_NOERR);
}

////////////////////////////////////////////////////////////////////////////////
// Read a scalar double variable from a CDF file by name (or an array, with   //
//  no bounds checking).                                                      //
////////////////////////////////////////////////////////////////////////////////
void cdf_read_double(int fid, const char *vname, double *vbuf)
{
  int vid;
  assert(vbuf);
  assert(nc_inq_varid(fid, vname, &vid) == NC_NOERR);
  assert(nc_get_var_double(fid, vid, vbuf) == NC_NOERR);
}

////////////////////////////////////////////////////////////////////////////////
// Read the final 1D vector from a 2D array of doubles by name.               //
////////////////////////////////////////////////////////////////////////////////
void cdf_read_last_double_subarr(int fid, const char *vname, int ld,
				 double *vbuf)
{
  size_t start[2], count[2];
  int    dimids[2], vid, ndims;

  assert(vbuf);
  assert(nc_inq_varid(fid, vname, &vid) == NC_NOERR);
  assert(nc_inq_varndims(fid, vid, &ndims) == NC_NOERR);
  assert(ndims == 2);
  assert(nc_inq_vardimid(fid, vid, dimids) == NC_NOERR);
  assert(nc_inq_dimlen(fid, dimids[0], &count[0]) == NC_NOERR);
  assert(static_cast<int>(count[0]) > 0);
  assert(nc_inq_dimlen(fid, dimids[1], &count[1]) == NC_NOERR);
  assert(static_cast<int>(count[1]) == ld);

  start[0] = count[0]-1;  start[1] = 0;  count[0] = 1;
  assert(nc_get_vara_double(fid, vid, start, count, vbuf) == NC_NOERR);
}

////////////////////////////////////////////////////////////////////////////////
// Read the final 1D vector from a 2D array of doubles by name.               //
////////////////////////////////////////////////////////////////////////////////
void cdf_read_double_subarr(int fid, const char *vname, int ld, int lt,
				 double *vbuf)
{
  size_t start[2], count[2];
  int    dimids[2], vid, ndims;

  assert(vbuf);
  assert(nc_inq_varid(fid, vname, &vid) == NC_NOERR);
  assert(nc_inq_varndims(fid, vid, &ndims) == NC_NOERR);
  assert(ndims == 2);
  assert(nc_inq_vardimid(fid, vid, dimids) == NC_NOERR);
  assert(nc_inq_dimlen(fid, dimids[0], &count[0]) == NC_NOERR);
  assert(static_cast<int>(count[0]) > 0);
  assert(nc_inq_dimlen(fid, dimids[1], &count[1]) == NC_NOERR);
  assert(static_cast<int>(count[1]) == ld);

//  start[0] = count[0]-1;  start[1] = 0;  count[0] = 1;
  start[0] = lt;  start[1] = 0;  count[0] = 1;
  assert(nc_get_vara_double(fid, vid, start, count, vbuf) == NC_NOERR);
}

#ifdef  NDEBUG_IS_DISABLE
  #define NDEBUG
  #undef NDEBUG_IS_DISABLE
#endif

bool CStconfig::loadvmecnetcdf(const char * fname, long fPostn)
{
  int ncid;
  mOK = false;
// Open the NetCDF file, if present
  if (nc_open(fname, NC_NOWRITE, &ncid)) {
     fprintf(stderr, "Could not open VMEC NetCDF file %s for reading.\n",fname);
     return(mOK=false);
  }
  mfname = fname;

  int iasym, nfp, nmodes,m,n,ns,mpol,ntor, ier_flag, nmodes_nyq;
  double vers;

  cdf_read_double(ncid, "version_", &vers);
  int iversion = int(vers*10000+0.1);
  cdf_read_int(ncid, "lasym__logical__", &iasym);
  cdf_read_int(ncid, "nfp", &nfp);
  cdf_read_int(ncid, "ns", &ns);
  cdf_read_int(ncid, "mpol", &mpol);
  cdf_read_int(ncid, "ntor", &ntor);
  cdf_read_int(ncid, "mnmax", &nmodes);
  cdf_read_int(ncid, "mnmax_nyq", &nmodes_nyq);
  cdf_read_int(ncid, "ier_flag", &ier_flag);
 
  vmecErrorCode = ier_flag;
  if(iasym) { 
    std::cerr<< "VMEC import: asymmetry mode"<<std::endl;
    vmecErrorCode = -3000;
    std::cerr<< "VMEC import: at the moment MConf can not handle wout-file in asymmetry mode"<<std::endl;
    nc_close(ncid);
    return false;    // at the moment can not handle asymmetry
  }
  
  if(ier_flag!=0) {  // check VMEC's error flag
    switch(ier_flag) {
      case 1:
        std::cerr<< "VMEC file '"<<fname<<"' contains error flag: "<<ier_flag<<"(bad initial jacobian)"<<std::endl; 
        break;
      case 4:
//        std::cerr<< "VMEC error flag: "<<ier_flag<<".  more iterations to reach force limit is needed"<<std::endl; 
        break;
      case 6:
        std::cerr<< "VMEC error flag: "<<ier_flag<<".  bad jacobian"<<std::endl; 
        break;
      default:
        std::cerr<< "VMEC error flag: "<<ier_flag<<".  unknown error"<<std::endl; 
        break;
    }
    if(!(ier_flag==4)) {
      nc_close(ncid);
      return false; 
    }
  }

  vmecCoordinates = true; // resize needs this flag
  vmecAsym = iasym != 0;  // resize needs this flag
  // m is the poloidal mode number, 0<=m<=mpol-1.
  // n is the toroidal mode number. -ntor<=n<=ntor.
  // if(m==0)    0<=n<=ntor.
  // The magnetic axis is included in VMEC output,
  // however we will discard it here (m1stSindex must be >= 1), see (ns-1)
  if (!resize(ns-1 , mpol-1, ntor)) {
    nc_close(ncid);
    return vmecCoordinates=false;
  }
  mNewFormat = true;
  mNp=nfp;

  

  double *mnum = new double [nmodes]; assert(mnum);
  double *nnum = new double [nmodes]; assert(nnum);
  double *rmnc = new double [nmodes]; assert(rmnc);
  double *zmns = new double [nmodes]; assert(zmns);
  double *lmns = new double [nmodes]; assert(lmns);
  double *rmns = new double [nmodes]; assert(rmns);
  double *zmnc = new double [nmodes]; assert(zmnc);
  double *lmn_c = new double [nmodes]; assert(lmn_c);
  double *mnum_nyq = new double [nmodes_nyq]; assert(mnum_nyq);
  double *nnum_nyq = new double [nmodes_nyq]; assert(nnum_nyq);
  double *bmnc = new double [nmodes_nyq]; assert(bmnc);
  double *gmnc = new double [nmodes_nyq]; assert(gmnc);
  double *bsupu = new double [nmodes_nyq]; assert(bsupu);
  double *bsupv = new double [nmodes_nyq]; assert(bsupv);
  double *bsubu = new double [nmodes_nyq]; assert(bsubu);
  double *bsubv = new double [nmodes_nyq]; assert(bsubv);
  double *bmn_s = new double [nmodes_nyq]; assert(bmn_s);
  double *gmn_s = new double [nmodes_nyq]; assert(gmn_s);
  double *bsupus = new double [nmodes_nyq]; assert(bsupus);
  double *bsupvs = new double [nmodes_nyq]; assert(bsupvs);
  double *bsubus = new double [nmodes_nyq]; assert(bsubus);
  double *bsubvs = new double [nmodes_nyq]; assert(bsubvs);

  // Read mode numbers
  cdf_read_double_vec(ncid, "xm", nmodes, mnum); // poloidal
  cdf_read_double_vec(ncid, "xn", nmodes, nnum); // toroidal
  for(int i = 1; i < ns; i++)  // i==0 is skiped : mag.axis is skiped 
  {
     int l = i-1+ m1stSindex;
     cdf_read_double_subarr(ncid, "rmnc", nmodes, i, rmnc); // R cos
     cdf_read_double_subarr(ncid, "zmns", nmodes, i, zmns); // Z sin
     cdf_read_double_subarr(ncid, "lmns", nmodes, i, lmns); // L sin
     if (iasym)
     {
        cdf_read_double_subarr(ncid, "rmns", nmodes, i, rmns); // R sin
        cdf_read_double_subarr(ncid, "zmnc", nmodes, i, zmnc); // Z cos
        cdf_read_double_subarr(ncid, "lmnc", nmodes, i, lmn_c); // L cos
     }
     for (int j = 0; j < nmodes; j++)
     {
        m = int(mnum[j]);
        n = int(nnum[j]/nfp);
        Rmn  (l, m, n) = rmnc[j];
        Zmn  (l, m, n) = zmns[j];
        lmn  (l, m, n) = lmns[j];
        if (iasym)
        {
           Rmns  (l, m, n) = rmns[j];
           Zmnc  (l, m, n) = zmnc[j];
           lmnc  (l, m, n) = lmn_c[j];
        }
     }
  }
  // Read mode numbers NYquist
  cdf_read_double_vec(ncid, "xm_nyq", nmodes_nyq, mnum_nyq); // poloidal
  cdf_read_double_vec(ncid, "xn_nyq", nmodes_nyq, nnum_nyq); // toroidal
  for(int i = 1; i < ns; i++)
  {
     int l = i-1+ m1stSindex;
     cdf_read_double_subarr(ncid, "bmnc", nmodes_nyq, i, bmnc); // B cos
     cdf_read_double_subarr(ncid, "gmnc", nmodes_nyq, i, gmnc); // G cos
     cdf_read_double_subarr(ncid, "bsupumnc", nmodes_nyq, i, bsupu); // B^u cos
     cdf_read_double_subarr(ncid, "bsupvmnc", nmodes_nyq, i, bsupv); // B^v cos
     cdf_read_double_subarr(ncid, "bsubumnc", nmodes_nyq, i, bsubu); // B_u cos
     cdf_read_double_subarr(ncid, "bsubvmnc", nmodes_nyq, i, bsubv); // B_v cos
     if (iasym)
     {
        cdf_read_double_subarr(ncid, "bmns", nmodes_nyq, i, bmn_s); // B sin
        cdf_read_double_subarr(ncid, "gmns", nmodes_nyq, i, gmn_s); // G sin
        cdf_read_double_subarr(ncid, "bsupumns", nmodes_nyq, i, bsupus); // B^u sin
        cdf_read_double_subarr(ncid, "bsupvmns", nmodes_nyq, i, bsupvs); // B^v sin
        cdf_read_double_subarr(ncid, "bsubumns", nmodes_nyq, i, bsubus); // B_u sin
        cdf_read_double_subarr(ncid, "bsubvmns", nmodes_nyq, i, bsubvs); // B_v sin
     }
     for (int j = 0; j < nmodes_nyq; j++)
     {
        // m is the poloidal mode number, 0<=m<=mpol-1.
        // n is the toroidal mode number. -ntor<=n<=ntor.
        // if(m==0)    0<=n<=ntor.
        m = int(mnum_nyq[j]);
        n = int(nnum_nyq[j]/nfp);
        if(m>mpol-1) continue; // skip extra modes
        if(n>ntor)  continue;  // skip extra modes
        if(n<-ntor) continue;  // skip extra modes
        bmn  (l, m, n) = bmnc[j];
        gmn  (l, m, n) = gmnc[j];
        bContraTheta_mn(l, m, n) = bsupu[j];
        bContraFi_mn(l,m,n) = bsupv[j];
        Phimn(l, m, n) = 0;
#if TRANSFORM_TO_BOOZER  
        double wrk;
        if(m!=0) wrk = bsubu[j]/m;
        else if(n!=0) wrk = bsubv[j]/(-n*mNp);
        else wrk = 0;
        Phimn(l, m, n) = wrk;
#endif
        if (iasym)
        {
           bmns            (l, m, n) = bmn_s[j];
           gmns            (l, m, n) = gmn_s[j];
           bContraTheta_mns(l, m, n) = bsupus[j];
           bContraFi_mns   (l, m, n) = bsupvs[j];
           Phimns          (l, m, n) = 0;
#if TRANSFORM_TO_BOOZER  
           double wrk;
           if(m!=0) wrk = bsubus[j]/m;
           else if(n!=0) wrk = bsubvs[j]/(-n*mNp);
           else wrk = 0;
           Phimns(l, m, n) = wrk;
#endif
        }
     }
  }
  for(int k = 0; k < ns; k++) {
    int j = k + m1stSindex;
    msFull.rw()[j] = double(k+1) / double(ns-1);  // VMEC full mesh
    ms.rw()    [j] = (0.5 + k)/double(ns-1);      // VMEC half mesh
  }
  // Read Profiles
  double *presf = new  double [ns]; assert(presf);
  double *phi = new  double [ns]; assert(phi);
  double *iotas = new  double [ns]; assert(iotas);
  double *pres = new  double [ns]; assert(pres);
  double *buco = new  double [ns]; assert(buco);
  double *bvco = new  double [ns]; assert(bvco);
  double *vp = new  double [ns]; assert(vp);

  cdf_read_double_vec(ncid, "presf", ns, presf);
  cdf_read_double_vec(ncid, "phi", ns, phi);
  cdf_read_double_vec(ncid, "iotas", ns, iotas);
  cdf_read_double_vec(ncid, "pres", ns, pres);
  cdf_read_double_vec(ncid, "buco", ns, buco);
  cdf_read_double_vec(ncid, "bvco", ns, bvco);
  cdf_read_double_vec(ncid, "vp", ns, vp);

  double presPrev=0;
  for(int i = 0; i < ns; i++) { // read the radial profiles
      int j = i + m1stSindex;
      //double  iotaf, presf, phipf, phi, jcuru, jcurv;
      //fscanf(fp,"%lf %lf %lf %lf %lf %lf", &iotaf,&presf,&phipf,&phi,&jcuru,&jcurv);
      mFlux = phi[i];
      if(i) mpprime.rw()[j-1] = (presf[i] - presPrev) /(msFull[j]-msFull[j-1]);
      if(i>=ns-2)
        mpprime.rw()[j] = extrapolateLCMS(msFull[j],msFull.constArray(),mpprime.constArray(),j);
      presPrev = presf[i];
  }
  for(int i = 0; i < (ns-1); i++) {
    int j = i + m1stSindex;
    miota.rw()[j] = iotas[i+1];
    mpres.rw()[j] = pres[i+1];
    mIpol.rw()[j] = bvco[i+1]*twopi/(mNp*mu0); //convert to w7x-format, poloidal current on a period
    mItor.rw()[j] = buco[i+1]*twopi/mu0;
    mg00.rw()[j] = vp[i+1] * 4.0 * M_PI * M_PI / mNp; // mg00 is the dV/ds/Np
  }

  double Aminor_p, Rmajor_p;
  cdf_read_double(ncid, "Rmajor_p", &Rmajor_p);
  cdf_read_double(ncid, "Aminor_p", &Aminor_p);

  ma  = Aminor_p;
  mR0 = Rmajor_p;
#if TRANSFORM_TO_BOOZER
  {
    for(int i = 0; i < (ns-1); i++) {
      int j = i + m1stSindex;
      for(int m = 0; m < mpol; m++) {
        int nmin = (m==0)?0:-ntor;
        for(int n = nmin; n <=ntor; n++) {
          double wrk = Phimn(j,m,n);
          double iotas = miota[j];
          double bvco  = mIpol[j]/(twopi/(mNp*mu0));
          double buco  = mItor[j]/(twopi/mu0);
          double pmn = (wrk - buco*lmn(j,m,n))/(bvco+iotas*buco);
          Phimn(j,m,n) = pmn;
  //      Phimn(j,m,n) = pmn*mNp/twopi; 
        }
      }  
    }
  }
#endif

  delete mnum;
  delete nnum;
  delete rmnc;
  delete zmns;
  delete lmns;
  delete rmns;
  delete zmnc;
  delete lmn_c;
  delete mnum_nyq;
  delete nnum_nyq;
  delete bmnc;
  delete gmnc;
  delete bsupu; 
  delete bsupv; 
  delete bsubu; 
  delete bsubv; 
  delete bmn_s; 
  delete gmn_s; 
  delete bsupus;
  delete bsupvs;
  delete bsubus;
  delete bsubvs;

  delete presf;
  delete pres ;
  delete phi  ;
  delete iotas;
  delete buco ;
  delete bvco ;
  delete vp   ;

  nc_close(ncid);
  return mOK = true; 
}
#endif



//*********************************************************
bool CStconfig::loadnemec(const char * fname, long fPostn)
{
  mOK = false;
    std::cerr<< "NEMEC import is not ready"<<std::endl;
  return mOK;
//=================================================
  FILE * fp = fopen(fname,"rb");
  if(fp==NULL) return(mOK=false);
  mfname = fname;

  fseek(fp,fPostn,SEEK_SET); //if(fPostn!=0)  fsetpos(fp,&fPostn);           // set position to read
 
  double gam,enfp,enrho,empol,entor,empnt,eiasym,phiedge;
  fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf",&gam,&enfp,&enrho,&empol,&entor,&empnt,&eiasym,&phiedge);
  int nfp, ns, nrho, mpol, ntor, mpnt, iasym;
  nfp    = int(enfp);
  ns=nrho= int(enrho);
  mpol   = int(empol);
  ntor   = int(entor);
  mpnt   = int(empnt);
  iasym  = int(eiasym);
  int nsin   = nrho-1;
  int mpol1  = mpol-1;

  if(iasym) { 
    std::cerr<< "NEMEC import: asymmetry mode"<<std::endl;
  }

  m1stSindex = 0;  // The magnetic axis is included
  vmecCoordinates = true; // resize needs this flag
  vmecAsym = iasym != 0;  // resize needs this flag
  if (!resize(ns , mpol-1, ntor))  return vmecCoordinates=false;
  mNewFormat = true;
  mNp=nfp;


  // VMEC full mesh and half mesh
  double ds = 1./double(ns-1);
  for(int i = 0; i <ns; i++) {
    msFull.rw()[i] = double(i)/double(ns-1);  // VMEC full mesh
    ms.rw()    [i] =   (i-0.5)/double(ns-1);  // half mesh
  }
  ms.rw()[0] = 0;
  msFull.rw()[ns-1] = 1;

  

  for(int j = 0; j < ns; j++) {
    for(int m = 0; m < mpol; m++) {
      int nmin = (m==0)?0:-ntor;
      for(int n = nmin; n <=ntor; n++) {
        double rmnc, zmns, rmns,zmnc;        // are defined on full-mesh
        double bsupuc, bsupvc,bsupus,bsupvs;   
        double lmns, lmn_c;
        fscanf(fp,"%lf %lf %lf %lf",&rmnc,&zmns,&rmns,&zmnc);
        fscanf(fp,"%lf %lf %lf %lf",&bsupuc,&bsupvc,&bsupus,&bsupvs);
        fscanf(fp,"%lf %lf",&lmns,&lmn_c);
        Rmn (j, m, n) = rmnc;    //  on full mesh
        Zmn (j, m, n) = zmns;    //  on full mesh
        Rmns(j, m, n) = rmns;    //  on full mesh, asymmetric part
        Zmnc(j, m, n) = zmnc;    //  on full mesh, asymmetric part
        bContraTheta_mn (j, m, n) = bsupuc; // B^u(mn)c contravariant poloidal component of B on half-mesh, cos(mu+nv) stellarator symmetric part
        bContraFi_mn    (j, m, n) = bsupvc; // B^v(mn)c contravariant toroidal component of B on half-mesh, cos(mu+nv) stellarator symmetric part
        bContraTheta_mns(j, m, n) = bsupus; // B^u(mn)s contravariant poloidal component of B on half-mesh, sin(mu+nv) asymmetric part
        bContraFi_mns   (j, m, n) = bsupvs; // B^v(mn)s contravariant toroidal component of B on half-mesh, sin(mu+nv) asymmetric part
        lmn  (j, m, n) = lmns;    // lmns   on half-mesh, sin(mu+nv) stellarator symmetric part
        lmnc (j, m, n) = lmn_c;   //  on half-mesh, cos(mu+nv) asymmetric part
        double bsubuc,bsubvc, bsubss;  
        double bsubus,bsubvs, bsubsc;  
        fscanf(fp,"%lf %lf %lf",&bsubuc,&bsubvc,&bsubss);
        fscanf(fp,"%lf %lf %lf",&bsubus,&bsubvs,&bsubsc);
        //bCovaTheta_mn(j, m, n) = bsubuc;  // B_u(mn)c covariant poloidal component of B on half-mesh, cos(mu+nv) stellarator symmetric part
        //bCovaFi_mn   (j, m, n) = bsubvc;  // B_v(mn)c covariant toroidal component of B on half-mesh, cos(mu+nv) stellarator symmetric part
        //bCovaS_mn    (j, m, n) = bsubss;  // B_s(mn)s covariant toroidal component of B on half-mesh, sin(mu+nv) stellarator symmetric part
#if ____TRANSFORM_TO_BOOZER  
          double wrk;
          if(m!=0) wrk = bsubuc/m;
          else if(n!=0) wrk = bsubvc/(-n*mNp);
          else wrk = 0;
          Phimn(j, m, n) = wrk;
          //////if(m!=0&&n!=0) { //test: the following must be zero
          //////  wrk = bsubu/m - bsubv/(-n*mNp);
          //////  wrk = bsubu/m - bsubv/(-n*mNp);
          //////}
#endif
        }
      }
    }

  double presPrev=0;
  for(int js = 1; js < ns; js++) { // read the radial profiles
    double hiota,hmass,hpres,hphip,hbuco,  
           hbvco,hphi,hvp,hoverr,fjcuru,fjcurv,hspecw;
    fscanf(fp,"%lf %lf %lf %lf %lf",&hiota,&hmass,&hpres,&hphip,&hbuco);    
    fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&hbvco,&hphi,&hvp,&hoverr,&fjcuru,&fjcurv,&hspecw);    
    miota.rw()[js] = hiota;
    mIpol.rw()[js] = -hbvco*twopi/(mNp*mu0); //convert to w7x-format, poloidal current on a period
    mItor.rw()[js] = -hbuco*twopi/mu0;
    mg00.rw() [js] = hvp * 4.0 * M_PI * M_PI; // mg00 is the dV/ds/Np
    // calculate dp/ds, ?? must be on half mesh
    if(js>2) mpprime.rw()[js-1] = (hpres - presPrev)/mu0/(ms[js]-ms[js-1]);
    if(js==ns-1) mpprime.rw()[js] = extrapolateLCMS(ms[js],ms.constArray(),mpprime.constArray(),js);
    presPrev = hpres;
  }
  miota.rw()[0] = miota[1];
  mIpol.rw()[0] = mIpol[1];
  mItor.rw()[0] = 0;
  mg00.rw() [0] = mg00[1];
  mpprime.rw()[0] = mpprime[1]; 

  mFlux=phiedge;  
  ma  = 0;
  mR0 = 0;

  fclose(fp);
#if ____TRANSFORM_TO_BOOZER
  {
    for(int i = 0; i < (ns-1); i++) {
      int j = i + m1stSindex;
      for(int m = 0; m < mpol; m++) {
        int nmin = (m==0)?0:-ntor;
        for(int n = nmin; n <=ntor; n++) {
          double wrk = Phimn(j,m,n);
          double iotas = miota[j];
          double bvco  = mIpol[j]/(twopi/(mNp*mu0));
          double buco  = mItor[j]/(twopi/mu0);
          double pmn = (wrk - buco*lmn(j,m,n))/(bvco+iotas*buco);
          Phimn(j,m,n) = pmn;
  //      Phimn(j,m,n) = pmn*mNp/twopi; 
        }
      }  
    }
  }
#endif
  return mOK = true;
}

};   //namespace MConf
