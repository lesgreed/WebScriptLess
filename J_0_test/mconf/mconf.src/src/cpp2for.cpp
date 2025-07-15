//*********************************************************************
// Fortran Interface to C++
/*
!USAGE:
  implicit none
  character (len=80) fname
  real(8) cyl(3),boozer(3),Bfield(3),epsTrunc/1d-6/,epsA/1d-4/,ficyl/0/
  integer(8) mcObj/0/         ! this variable holds the address of C++ object
  integer(4) False/0/,mcload

  EXTERNAL mcload,mcsetB0,mctruncate,mcsetAccuracy,mcsetB0,mccyl2mag
! Here EXTERNAL operator is used, but for true FORTRAN-90 program you have to
! use FORTRAN-statement INTERFACE for defining prototypes of subroutine called:
!
!   INTERFACE
!    INTEGER(4) FUNCTION MCLOAD(mcObj,fname)
!     INTEGER(8) mcObj
!     CHARACTER *(*) fname
!    END FUNCTION MCLOAD
!   END INTERFACE

  fname = 'w7x-sc1.bc'
  if  (mcload (mcObj,fname)==False ) stop     ! stop if error

  call mctruncate   (mcObj,epsTrunc)          ! truncate spectrum
  call mcsetAccuracy(mcObj,epsA)              ! set accuracy of coordinate transformation
  call mcsetB0      (mcObj,2.5d0,ficyl)       ! scale B, set value of B on magn. axis

  cyl = (/6d0,1d-2,5d-1/)

  call mccyl2mag (mcObj,cyl,boozer,Bfield) ! coordinate transformation cyl -> boozer

  if(boozer(1) > 1)  print '(/a)', 'Point outside plasma!'
  print '(/a,3f14.6 )','cyl   = (r, fi, z)  =',cyl
  print '( a,3f14.6 )','boozer=(s,theta,phi)=',boozer
  print '( a,3f14.6/)','Bfield=(Br, Bfi, Bz)=',Bfield

  call mcfree(mcObj)                          ! destroy C++ object and free memory
  end
*/
//
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "../include/C3dMesh.h"
#include "../include/CPfunction.h"

#include <iostream>

using MConf::C3dMesh;
using MConf::CStconfig;
using MConf::Vector3d;
using MConf::CProfile;
using MConf::CPfunction;

using namespace std;

const static double   pi = 3.14159265358979323846;
const static double   degree = pi/180;

#undef True
#undef False

const static int True(1);
const static int False(0);

// see option -nomixed_str_len_arg for Intel Fortran 
// or /iface:nomixed_str_len_arg  for Visual Fortran 6
#define NOMIXED_STR_LEN_ARG

//_WIN32 is defined for Win32 and Win64 applications. Always defined.
#if ( defined(_WIN32) && !( defined(__MINGW32__) || defined(__MINGW64__) ) )
    #ifdef _MSC_VER
 // for Visual Studio 2012 and Intel Fortran 2013 with default calling convention
   //   #if _MSC_VER > 1600      // Visual Studio 2012 or higher
      #if _MSC_VER > 1400      // Visual Studio 2008 or higher
        #define STDCALL 
      #elif _MSC_VER == 1200
        // for Visual Fortran 6 with default calling convention
        // and /iface:nomixed_str_len_arg
        #define STDCALL __stdcall 
      #else
        #define STDCALL __stdcall 
      #endif
    #else
      #define STDCALL  
    #endif
#endif

#if ( defined(_WIN32) && !( defined(__MINGW32__) || defined(__MINGW64__) ) )
    #define WINDLLEXPORT  
    #ifdef _MT
      #pragma message( "===> Multithread library" )
      const static char* libDescription = "Multithread MConf library";
    #else
      #pragma message( "===> Singlethread library")
     const static char* libDescription = "Singlethread MConf library";
    #endif
#endif


// Note: for Intel Fortran 2013
// if Fortran program under Windows uses BIND(C), like in the the following declaration 
//    INTEGER(4) FUNCTION mcload(mconf,name) BIND(C)
// then all function names must be lowcase names

#if (defined( __unix__ ) || defined(linux) || defined(__sun) || defined(__APPLE__) || defined(__MINGW32__) || defined(__MINGW64__) )
  #define WINDLLEXPORT  
  #define STDCALL
  #define PROFILEFUNCTION  profilefunction_
  #define MCLOAD        mcload_
  #define MCISVMEC      mcisvmec_
  #define MCISTOKPEST   mcistokpest_
  #define MCLOADVMEC    mcloadvmec_
  #define MCLOADVMECNETCDF    mcloadvmecnetcdf_
  #define MCLOADLHD     mcloadlhd_
  #define MCLOADEFIT    mcloadefit_
  #define MCISOK        mcisok_
  #define MCISMESHOK    mcismeshok_
  #define MCFREE        mcfree_
  #define MCCLONE       mcclone_
  #define MCTRUNCATE    mctruncate_
  #define MCSETACCURACY mcsetaccuracy_
  #define MCSETMINB0    mcsetminb0_
  #define MCSETB0       mcsetb0_
  #define MCGETB0       mcgetb0_
  #define MCSETB00      mcsetb00_
  #define MCGETB00      mcgetb00_
  #define MCSCALEB      mcscaleb_
  #define MCSETICOIL    mcseticoil_
  #define MCSETBDIRECTION mcsetbdirection_
  #define MCGETBDIRECTION mcgetbdirection_
  #define MCRESTOREBDIRECTION mcrestorebdirection_

  #define MCGETCOVAGRADB  mcgetcovagradb_
  #define MCGETGRADB      mcgetgradb_
  #define MCGETBANDGRADIENTS  mcgetbandgradients_
  #define MCGETBCURLBGRADFCYL     mcgetbcurlbgradfcyl_
  #define MCCYL2MAG     mccyl2mag_
  #define MCCYL2MAG2    mccyl2mag2_
  #define MCXYZ2MAG2    mcxyz2mag2_
  #define MCXYZ2S       mcxyz2s_
  #define MCCYL2S       mccyl2s_
  #define MCREFF        mcreff_
  #define MCTORFLUX     mctorflux_
  #define MCVPRIME      mcvprime_
  #define MCIOTA        mciota_
  #define MCIOTAPRIME   mciotaprime_
  #define MCVOLUME      mcvolume_
  #define MCFTRAPPED    mcftrapped_
  #define MCBAVRG       mcbavrg_
  #define MCB2AVRG      mcb2avrg_
  #define MCR2AVRG      mcr2avrg_
  #define MCB_2AVRG     mcb_2avrg_
  #define MCGRADS2B_2AVRG     mcgrads2b_2avrg_
  #define MCR_2AVRG     mcr_2avrg_
  #define MCDB12AVRG    mcdb12avrg_
  #define MCDB32AVRG    mcdb32avrg_
  #define MCBDB12AVRG   mcbdb12avrg_
  #define MCNEOPOLARIZATION   mcneopolarization_
  #define MCNEOPOLARIZATIONPASS   mcneopolarizationpass_
  #define MCNEOPOLARIZATIONTRAP   mcneopolarizationtrap_
  #define MCRHPOTENCIAL   mcrhpotencial_

  #define MCIPOL         mcipol_
  #define MCITOR         mcitor_
  #define MCB            mcb_
  #define MCJACOBIAN     mcjacobian_
  #define MCBANDJACOBIAN     mcbandjacobian_
  #define MCROPARINVAVRG mcroparinvavrg_

  #define MCBPHI        mcbphi_
  #define MCBMIN        mcbmin_
  #define MCBMAX        mcbmax_

  #define MCPAVRG       mcpavrg_
  #define MCINTINVPAVRG mcintinvpavrg_
  #define MCTORFLUX2POLFLUX  mctorflux2polflux_

  #define MCJACOBIANB2  mcjacobianb2_

  #define MCRMN          mcrmn_
  #define MCZMN          mczmn_
  #define MCFMN          mcfmn_
  #define MCBMN          mcbmn_
  #define MCBOOZERBMN    mcboozerbmn_
     

  #define MCAVERAGECROSSSECTIONAREA  mcaveragecrosssectionarea_
  #define MCR0          mcr0_
  #define MCRMINOR      mcrminor_
  #define MCNPERIODS    mcnperiods_
  #define MCCURRENTSIGNTOINCREASEIOTA mccurrentsigntoincreaseiota_

  #define MCTRACELCMS     mctracelcms_
  #define MCCREATEBORDER  mccreateborder_
  #define MCISINSIDEBORDER mcisinsideborder_
  #define MCWRITE          mcwrite_

// mesh
  #define MCM3DGETRAYENTRYPOINT  mcm3dgetrayentrypoint_
  #define MCCREATEMESH    mccreatemesh_
  #define MCCREATEMESHUSINGSYMMETRY    mccreatemeshusingsymmetry_
  #define MCWRITEMESH     mcwritemesh_
  #define MCM3DCYL2S      mcm3dcyl2s_
  #define MCM3DXYZ2S      mcm3dxyz2s_
  #define MCM3DGETDBCYL   mcm3dgetdbcyl_
  #define MCM3DGETDBXYZ   mcm3dgetdbxyz_
  #define MCM3DGETGRADS   mcm3dgetgrads_
  #define MCM3DGETGRADSXYZ mcm3dgetgradsxyz_
  #define MCGETGRADSXYZ    mcgetgradsxyz_

  #define MCXYZISINSIDE   mcxyzisinside_
  #define MCCYLISINSIDE   mccylisinside_

  #define MCGETALLXYZ     mcgetallxyz_
  #define MCM3DGETALLXYZ  mcm3dgetallxyz_
  #define MCM3DGETALLCYL  mcm3dgetallcyl_

  #define MCM3DGETSBXYZ   mcm3dgetsbxyz_
  #define MCM3DGETSBCYL   mcm3dgetsbcyl_
  #define MCGETSBCYL      mcgetsbcyl_
  #define MCGETSBXYZ      mcgetsbxyz_
  #define MCGETBCYL       mcgetbcyl_
  #define MCGETDBCYL      mcgetdbcyl_

  #define MCM3DSETSMAX    mcm3dsetsmax_
  #define MCM3DGETSMAX    mcm3dgetsmax_
  #define MCGETBMINMAX    mcgetbminmax_
  #define MCGETBMINMAXGLOBAL mcgetbminmaxglobal_
  #define MCRAXIS         mcraxis_
  #define MCGETEXTENT       mcgetextent_
  #define MCGETEXTENTLCMS   mcgetextentlcms_
  #define MCMIXCOORD2CYL    mcmixcoord2cyl_
  #define MCMIXCOORD2XYZ    mcmixcoord2xyz_
  #define MCMAG2CYL         mcmag2cyl_
  #define MCMAG2XYZ         mcmag2xyz_
  #define MCPEST2VMEC       mcpest2vmec_
  #define MCVMEC2PEST       mcvmec2pest_

  #define MCBOOZER2VMEC     mcboozer2vmec_
  #define MCVMEC2BOOZER     mcvmec2boozer_

  #define MCCOORDFIELDLINE  mccoordfieldline_
  #define MCSAVEREDUCED     mcsavereduced_
  #define MCGETGFACTORS     mcgetgfactors_
  #define MCGETEPSEFF       mcgetepseff_
  #define MCEPSEFFPARAMS    mcepseffparams_

  #define MCSFLCOORDONBLINEATTHETA mcsflcoordonblineattheta_
  #define MCSFLCONTRABASIS    mcsflcontrabasis_
  #define MCBOOZERCONTRABASIS mcboozercontrabasis_
  #define MCPRESSUREPRIME mcpressureprime_
  #define MCPRESSURE      mcpressure_
  #define SETARGV SETARGV_
  #define MCGETMAGCOORDATB mcgetmagcoordatb_
  #define MCJINVARIANT     mcjinvariant_
  #define MCUSEMIXEDPRODUCTFORJACOBIAN mcusemixedproductforjacobian_
#endif


template <class T> static inline T mmin (T x, T y) { return (x<y)?x:y; }
template <class T> static inline T mmax (T x, T y) { return (x>y)?x:y; }

static void strtrim(char *f)  // remove initial and trailing blanks
{
  char *p,*p2;
  int i=strlen(f);
  for(p=f; p!=f+i-1; p++) if(*p!=' ') break;
  i=strlen(p);
  for(p2=p+i-1; p2!=p; p2--) if(*p2!=' ') {*(p2+1)=0; break; }
  memmove(f,p,p2-p+2);
};

// remove initial and trailing blanks from FORTRAN string: character(i) fname
// size of dest must be large than fname
static void strtrim(char *dest, char *fname, int i)
{
  char *p;
  // remove initial blanks
  for(p=fname; p!=fname+i-1; p++) if(*p!=' ') break;
  i -= p-fname;
  strncpy(dest,p,i);
  dest[i]=0;
  // remove trailing blanks
  for(p=dest+i-1; p!=dest; p--) if(*p!=' ') {*(p+1)=0; break; }
}


//*********************************************************************
extern "C" WINDLLEXPORT double STDCALL PROFILEFUNCTION(double *x,double *g)
{
  CPfunction profile;
  return profile.PF(*x,g);
}

//*********************************************************************
// 'destructor'
extern "C" WINDLLEXPORT void STDCALL MCFREE(C3dMesh **mConf)
{
  if (*mConf) delete (*mConf);
  *mConf=0;
}

//*********************************************************************
// clone object  C3dMesh
extern "C" WINDLLEXPORT C3dMesh * STDCALL MCCLONE(C3dMesh **mConf)
{
  C3dMesh &mc = **mConf;
  return mc.clone();
}

//*********************************************************************
// create object  C3dMesh
// and load magnetic coordinate data file: *.bc; *.bc_old; *.bin4; *.bin8
// This function must be called before others.
// The first call must be made with *mConf==0, the following calls
//  will reload configuration into the same object
//
static int MC_load(C3dMesh **mConf, char * fname, int i, CStconfig::fileFormat fmt)
{
  char *fn = new char[i+2];
  strtrim(fn,fname,i);
  if (*mConf==0) *mConf = new C3dMesh;   // create object C3dMesh and return to caller the address
  bool ok = (*mConf)->load(fn,fmt);
  if(!ok) {
    cerr << "MCLOAD: Loading error, file:'"<<fn<<"'\n";
    delete (*mConf);
    *mConf = 0;      // return NULL if errors
  }
  delete fn;
//  (*mConf)->M3DpositionTestOff();
  return ok?True:False;
}

//*********************************************************************
// create object  C3dMesh
// and load magnetic coordinate data file: *.bc; *.bc_old; *.bin4; *.bin8
// This function must be called before others.
// The first call must be made with *mConf==0, the following calls
//  will reload configuration into the same object
//
extern "C" WINDLLEXPORT int STDCALL MCLOAD(C3dMesh **mConf, char * fname, int i)
{
  return MC_load(mConf, fname, i, CStconfig::UNDEFINED);
}

extern "C" WINDLLEXPORT int STDCALL MCLOADVMEC(C3dMesh **mConf, char * fname, int i)
{
  return MC_load(mConf, fname, i, CStconfig::VMEC);
}
#ifdef NETCDF
extern "C" WINDLLEXPORT int STDCALL MCLOADVMECNETCDF(C3dMesh **mConf, char * fname, int i)
{
  return MC_load(mConf, fname, i, CStconfig::VMECNETCDF);
}
#endif
extern "C" WINDLLEXPORT int STDCALL MCLOADLHD(C3dMesh **mConf, char * fname, int i)
{
  return MC_load(mConf, fname, i, CStconfig::LHD);
}

extern "C" WINDLLEXPORT int STDCALL MCLOADEFIT(C3dMesh **mConf, char * fname, int i)
{
  return MC_load(mConf, fname, i, CStconfig::EFIT);
}

//*********************************************************************
extern "C" WINDLLEXPORT int STDCALL MCISMESHOK(C3dMesh **mConf)
{
  C3dMesh &mc = **mConf;
  return mc.isMeshOK()?True:False;
}

//*********************************************************************
extern "C" WINDLLEXPORT int STDCALL MCISOK(C3dMesh **mConf)
{
  C3dMesh &mc = **mConf;
  return mc.isOK()?True:False;
}

//*********************************************************************
// set truncation level
extern "C" WINDLLEXPORT void STDCALL MCTRUNCATE(C3dMesh **mConf, double *epsTrunc)
{
  C3dMesh &mc = **mConf;
  mc.truncate(*epsTrunc);
}

//*********************************************************************
// set accuracy of coordinate transformation
extern "C" WINDLLEXPORT void STDCALL MCSETACCURACY(C3dMesh **mConf, double *epsA)
{
  C3dMesh &mc = **mConf;
  mc.setAccuracy(*epsA);
}

//*********************************************************************
// set value of B on magn. axis at cylindrical angle ficyl
extern "C" WINDLLEXPORT void STDCALL MCSETB0(C3dMesh **mConf, double *b, double *ficyl)
{
  C3dMesh &mc = **mConf;
  mc.setB0(*b, *ficyl);
}

//*********************************************************************
// set minimum value of B on magn. axis
extern "C" WINDLLEXPORT void STDCALL MCSETMINB0(C3dMesh **mConf, double *b)
{
  C3dMesh &mc = **mConf;
  mc.setB0(*b);
}

//*********************************************************************
// Get value of B on magn. axis at cylindrical angle ficyl
extern "C" WINDLLEXPORT double STDCALL MCGETB0(C3dMesh **mConf,double *ficyl)
{
  C3dMesh &mc = **mConf;
  return mc.getB0(*ficyl);
}

//*********************************************************************
// set average B on magn. axis
extern "C" WINDLLEXPORT void STDCALL MCSETB00(C3dMesh **mConf, double *b)
{
  C3dMesh &mc = **mConf;
  mc.setB00(*b);
}

//*********************************************************************
// 
extern "C" WINDLLEXPORT void STDCALL MCSETICOIL(C3dMesh **mConf, double *I)
{
  C3dMesh &mc = **mConf;
  mc.setIcoil(*I);
}

//*********************************************************************
// 
extern "C" WINDLLEXPORT void STDCALL MCSETBDIRECTION(C3dMesh **mConf, int *dir)
{
  C3dMesh &mc = **mConf;
  mc.setBdirection(*dir);
}

//*********************************************************************
// 
extern "C" WINDLLEXPORT int STDCALL MCGETBDIRECTION(C3dMesh **mConf)
{
  C3dMesh &mc = **mConf;
  return mc.getBdirection();
}

//*********************************************************************
// 
extern "C" WINDLLEXPORT void STDCALL MCRESTOREBDIRECTION(C3dMesh **mConf)
{
  C3dMesh &mc = **mConf;
  mc.restoreBdirection();
}

//*********************************************************************
// Get average B on magn. axis
extern "C" WINDLLEXPORT double STDCALL MCGETB00(C3dMesh **mConf)
{
  C3dMesh &mc = **mConf;
  return mc.getB00();
}

//*********************************************************************
// get PEST-coordinates 'pest' of the point lying on the B-line at poloidal angle 'thetaStar';  
// the B-line is going through point pest0
//
extern "C" WINDLLEXPORT void STDCALL MCSFLCOORDONBLINEATTHETA(C3dMesh **mConf, const double *pest0, 
                        const double *thetaStar, double *pest)
{
  C3dMesh &mc = **mConf;
  Vector3d mgc(pest0),b,grB,grs,grT,grP;
  Vector3d pst  = mc.SFLcoordOnBlineAtTheta(mgc, *thetaStar);
  pst.copyTo(pest);   
}

//*********************************************************************
// Get PEST gradients in cylindrical coordinates
// @param[in]  pest is the Boozer or PEST coordinate
// @param[out] mag  is the Boozer coordinates (if input file is Boozer file) or 
//  VMEC coordinates (if input file is VMEC-wout file)
// @param[out] grad... are gradients at point pest
// mConf must be originated from Boozer or VMEC file
extern "C" WINDLLEXPORT void STDCALL MCSFLCONTRABASIS(C3dMesh **mConf, const double *pest, 
                        double *gradS,double *gradThetaStar,double *gradPhi,double *mag)
{
  C3dMesh &mc = **mConf;
  Vector3d pst(pest),grS, grT, grP, mg;
  
  mc.SFLcontraBasis(pst,grS,grT, grP,mg);
 
  grS.copyTo(gradS);       // copy to array grads
  grT.copyTo(gradThetaStar);  
  grP.copyTo(gradPhi);    
  mg.copyTo(mag);     
}

//*********************************************************************
// Get Boozer gradients in cylindrical coordinates
// u    is the Boozer coordinate
// mag  is the Boozer or VMEC coordinate
// grad... are gradients at point u
// mConf must be originated from Boozer or VMEC file
extern "C" WINDLLEXPORT void STDCALL MCBOOZERCONTRABASIS(C3dMesh **mConf, const double *u, 
                        double *gradS,double *gradThetaStar,double *gradPhi,double *mag)
{
  C3dMesh &mc = **mConf;
  Vector3d u1(u),grS, grT, grP, mg;
  
  mc.BoozerContraBasis(u1,grS,grT, grP,mg);
 
  grS.copyTo(gradS);       // copy to array grads
  grT.copyTo(gradThetaStar);  
  grP.copyTo(gradPhi);    
  mg.copyTo(mag);     
}


//*********************************************************************
// Get gradients in cylindrical coordinates
//  gradB  = grad(|B|) 
//  gradF  = grad(torFlux)
//  gradTh = grad(theta)
//  gradPh = grad(phi) 
extern "C" WINDLLEXPORT void STDCALL MCGETBANDGRADIENTS(C3dMesh **mConf, double *magCoord, 
                     double *B,double *gradB,double *gradF,double *gradTh,double *gradPh)
{
  C3dMesh &mc = **mConf;
  Vector3d mgc(magCoord),b,grB,grs,grT,grP;
  
  mc.getBandGradients(mgc,b,grB,grs,grT,grP);
  double Fluxmax = mc.Flux(1.);
  Vector3d grF = grs*Fluxmax; // grad(torFlux)

  b.copyTo(B);             // B=(Br,Bfi,Bz) in cyl. coordinates
  grB.copyTo(gradB);       // copy to array gradB
  grF.copyTo(gradF);       // copy to array gradF
  grT.copyTo(gradTh);      // grad(theta)  in cyl. coord, theta is the poloidal angle
  grP.copyTo(gradPh);      // grad(phi)    in cyl. coord, phi is the toroidal angle
}

/// The method returns B and curl(B) in cylindrical coordinates
/// @param[in]  mConf is the pointer to memory at which the address of C3dMesh-object is stored.
///  @param[in]  u is the point (s,theta,phi) in magnetic coordinates.
///  @param[out]  B  is the magnetic field  in cylindrical coordinates.
///  @param[out]  curlB  is the curl(B). 
///  @param[out]  gradF  is the grad(torFlux)
extern "C" WINDLLEXPORT void STDCALL MCGETBCURLBGRADFCYL(C3dMesh **mConf, double *u, 
                     double *B,double *curlB,double *gradF)
{
  C3dMesh &mc = **mConf;
  Vector3d mgc(u),b,curl,grs;
  
  mc.getBcurlBGradscyl(mgc,b,curl,grs);
  double Fluxmax = mc.Flux(1.);
  Vector3d grF = grs*Fluxmax; // grad(torFlux)

  b.copyTo(B);             // B=(Br,Bfi,Bz) in cyl. coordinates
  curl.copyTo(curlB);      
  grF.copyTo(gradF);       // copy to array gradF
}



//*********************************************************************
//  @return grad(|B|) in cylindrical coordinates
extern "C" WINDLLEXPORT void STDCALL MCGETGRADB(C3dMesh **mConf, double *magCoord, double *gradB)
{
  C3dMesh &mc = **mConf;
  Vector3d bz,gr;
  bz.setFrom(magCoord);   // copy magCoord to vector bz
  gr = mc.gradB(bz);
  gr.copyTo(gradB);       // copy to array gradb
}

//*********************************************************************
// @return covariant components of  grad(|B|) in magnetic coordinates
// dB/ds, dB/dtheta, dB/dphi
extern "C" WINDLLEXPORT void STDCALL MCGETCOVAGRADB(C3dMesh **mConf, double *magCoord, double *gradB)
{
  C3dMesh &mc = **mConf;
  Vector3d bz,gr;
  bz.setFrom(magCoord);       // copy booz to vector bz
  gr = mc.getCovaGradB(bz);
  gr.copyTo(gradB);       // copy to array gradb
}

//*********************************************************************
// set average B on magn. axis
extern "C" WINDLLEXPORT void STDCALL MCSCALEB(C3dMesh **mConf, double *factor)
{
  C3dMesh &mc = **mConf;
  mc.scaleB(*factor);
}


//*********************************************************************
// Coordinate transformation
// Input :
//  cyl=(r,fi,z) - point in cyl. coordinates
// Output:
//  magCoord=(s,theta,phi)
//  Bfield=(Br,Bfi,Bz) in cyl. coordinates
extern "C" WINDLLEXPORT void STDCALL MCCYL2MAG(C3dMesh **mConf,double *cyl,double *magCoord,double *Bfield)
{
  Vector3d c,b;
  c.setFrom(cyl);            // copy cyl to vector c
  C3dMesh &mc = **mConf;
  if(mc.cylIsInside(c) ) {   // function cylIsInside returns 'TRUE' if the point lies inside LCMS
    b = mc.cyl2mag(c);       // b will be Vector3d(s,theta,phi)
    b.copyTo(magCoord);      // copy to fortran array
    b = mc.getBcyl();
    b.copyTo(Bfield);        // Bfield=(Br,Bfi,Bz) in cyl. coordinates
  }
  else {
    b = 0;
    b.copyTo(Bfield);
    b[0]=1e4;              // point lies outside plasma, set huge value of 's' as signal
    b.copyTo(magCoord);    // copy to array
  }
}

//*********************************************************************
// Faster version of MCCYL2MAG, without testing whether point lies inside LCMS
// Coordinate transformation
// Input :
//  cyl=(r,fi,z) - point in cyl. coordinates
// Output:
//  magCoord=(s,theta,phi)
//  Bfield=(Br,Bfi,Bz) in cyl. coordinates
extern "C" WINDLLEXPORT void STDCALL MCCYL2MAG2(C3dMesh **mConf,double *cyl,double *magCoord,double *Bfield)
{
  C3dMesh &mc = **mConf;
  Vector3d c(cyl),b;
  b = mc.cyl2mag(c);    // b is a Vector3d(s,theta,phi)
  b.copyTo(magCoord);      // copy to array
  b = mc.getBcyl();
  b.copyTo(Bfield);        // Bfield=(Br,Bfi,Bz) in cyl. coordinates
}

//*********************************************************************
// Coordinate transformation
// Input :
//  xyz -- point in cartesian coordinates
// Output:
//  magCoord=(s,theta,phi)
//  Bfield  in cartesian coordinates
extern "C" WINDLLEXPORT void STDCALL MCXYZ2MAG2(C3dMesh **mConf,double *xyz,double *magCoord,double *Bfield)
{
  C3dMesh &mc = **mConf;
  Vector3d c(xyz),b;
  b = mc.xyz2mag(c);    // b is a Vector3d(s,theta,phi)
  b.copyTo(magCoord);      // copy to array
  b = mc.getBxyz();
  b.copyTo(Bfield);        // Bfield=(Bx,By,Bz) in cart. coordinates
}


//*********************************************************************
// Coordinate transformation,
// Find label of flux surface
// Input :
//  xyz -- point in cartesian coordinates
// Output:
// function returns label of flux surface -- normalized flux 's'
//
extern "C" WINDLLEXPORT double STDCALL MCXYZ2S(C3dMesh **mConf,double *xyz)
{
  C3dMesh &mc = **mConf;
  Vector3d c(xyz);
  return mc.xyz2s(c);
}

//*********************************************************************
// Coordinate transformation,
// Find label of flux surface
// Input :
//  cyl -- point in cartesian coordinates
// Output:
// function returns label of flux surface -- normalized flux 's'
//
extern "C" WINDLLEXPORT double STDCALL MCCYL2S(C3dMesh **mConf,double *cyl)
{
  C3dMesh &mc = **mConf;
  Vector3d c(cyl);
  return mc.cyl2s(c);
}

//*********************************************************************
extern "C" WINDLLEXPORT void STDCALL MCGETSBXYZ(C3dMesh **mConf,double *xyz,double *S,double *B)
{
  C3dMesh &mc = **mConf;
  Vector3d c(xyz),b;
  *S=mc.xyz2s(c);
  b = mc.getBxyz();
  b.copyTo(B);               // Bfield=(Bx,By,Bz) in cart. coordinates
}

//*********************************************************************
extern "C" WINDLLEXPORT void STDCALL MCGETSBCYL(C3dMesh **mConf,double *cyl,double *S,double *B)
{
  C3dMesh &mc = **mConf;
  Vector3d c(cyl),b;
  *S=mc.cyl2s(c);
  b = mc.getBcyl();
  b.copyTo(B);               // Bfield=(Br,Bfi,Bz) in cyl. coordinates
}


//*********************************************************************
// Input :
//  magCoord=(s,theta,phi) -- point in magnetic coordinates
// Output:
//  B  -- magnetic field  in cylindrical coordinates
extern "C" WINDLLEXPORT void STDCALL MCGETBCYL(C3dMesh **mConf,double *magCoord,double *B)
{
  C3dMesh &mc = **mConf;
  Vector3d b,mag(magCoord);
  mc.NJacobian(mag);
  b = mc.getBcyl();
  b.copyTo(B);               // Bfield=(Br,Bfi,Bz) in cyl. coordinates
}

//*********************************************************************
// s to toroidal Flux
extern "C" WINDLLEXPORT double STDCALL MCTORFLUX(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.Flux(*s);
}

//*********************************************************************
// The method returns the direction of the toroidal current, which increases 
// the absolute value of the iota. 
// @return 1 or -1 depending on the sign of the current 
//         in the right-handed cylindrical coordinates system.
extern "C" WINDLLEXPORT int STDCALL MCCURRENTSIGNTOINCREASEIOTA(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.getSignOfCurrentThatIncreasesAbsIota();
}

//*********************************************************************
// s to r_eff
extern "C" WINDLLEXPORT double STDCALL MCREFF(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.r(*s);
}

//*********************************************************************
// Vprime(s) == dVolume/ds
// Vprime is negative if (s,theta,phi)-coordinate system is left handed
extern "C" WINDLLEXPORT double STDCALL MCVPRIME(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.Vprime(*s);
}

//*********************************************************************
// Volume(s)
// volume is negative if (s,theta,phi)-coordinate system is left handed
extern "C" WINDLLEXPORT double STDCALL MCVOLUME(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.Volume(*s);
}

//*********************************************************************
// Fraction of trapped particles on surface s
extern "C" WINDLLEXPORT double STDCALL MCFTRAPPED(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.ftrapped(*s);
}

//*********************************************************************
// <B> on surface s
extern "C" WINDLLEXPORT double STDCALL MCBAVRG(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.Bavrg(*s);
}

//*********************************************************************
// <B^2> on surface s
extern "C" WINDLLEXPORT double STDCALL MCB2AVRG(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.B2avrg(*s);
}

//*********************************************************************
// <B^-2> on surface s
extern "C" WINDLLEXPORT double STDCALL MCB_2AVRG(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.B_2avrg(*s);
}

//*********************************************************************
// <grads^2*B^-2> on surface s
extern "C" WINDLLEXPORT double STDCALL MCGRADS2B_2AVRG(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.grads2B_2avrg(*s);
}

//*********************************************************************
// <R^2> on surface s
extern "C" WINDLLEXPORT double STDCALL MCR2AVRG(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.R2avrg(*s);
}

//*********************************************************************
// <R^-2> on surface s
extern "C" WINDLLEXPORT double STDCALL MCR_2AVRG(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.R_2avrg(*s);
}

//*********************************************************************
//
extern "C" WINDLLEXPORT double STDCALL MCNEOPOLARIZATION(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.NeoPolarization(*s);
}

//*********************************************************************
//
extern "C" WINDLLEXPORT double STDCALL MCNEOPOLARIZATIONPASS(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.NeoPolarizationPass(*s);
}

//*********************************************************************
//
extern "C" WINDLLEXPORT double STDCALL MCNEOPOLARIZATIONTRAP(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.NeoPolarizationTrap(*s);
}

//*********************************************************************
// toroidal covariant component of B in Tesla*m
// Tokamak's people name it as poloidal current function
extern "C" WINDLLEXPORT double STDCALL MCBPHI(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  double pi    = 3.14159265358979323846;
  double twopi = 2*pi;
  double mu0   = 4e-7*pi; // [N/A^2 = henry/m]
  double Bphi  = mc.Ip(*s)*mu0/twopi;
  return Bphi;
}

//*********************************************************************
// poloidal current in Amp
extern "C" WINDLLEXPORT double STDCALL MCIPOL(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.Ip(*s);
}

//*********************************************************************
// toroidal current in Amp
extern "C" WINDLLEXPORT double STDCALL MCITOR(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.It(*s);
}

//*********************************************************************
// |B|
extern "C" WINDLLEXPORT double STDCALL MCB(C3dMesh **mConf,double *s,double *polAngle,double *torAngle)
{
  Vector3d magCoord(*s,*polAngle,*torAngle);
  C3dMesh &mc = **mConf;
  return mc.B(magCoord);
}

//*********************************************************************
// Calculate Jacobian in coordinates (s,theta,phi)
extern "C" WINDLLEXPORT double STDCALL MCJACOBIAN(C3dMesh **mConf,double *s,double *polAngle,double *torAngle)
{
  Vector3d magCoord(*s,*polAngle,*torAngle);
  C3dMesh &mc = **mConf;
  return mc.Jacobian(magCoord);
}


//*********************************************************************
// Calculate Jacobian in coordinates (s,theta,phi)
extern "C" WINDLLEXPORT void STDCALL MCBANDJACOBIAN(C3dMesh **mConf,double *s,double *polAngle,
                                                    double *torAngle,double *B,double *jac)
{
  Vector3d magCoord(*s,*polAngle,*torAngle);
  C3dMesh &mc = **mConf;
  *jac = mc.Jacobian(magCoord,B);
}


//*********************************************************************
// for A.Mishchenko
// Rosenbluth-Hinton residual zonal flow for tokamak geometry
extern "C" WINDLLEXPORT double STDCALL MCRHPOTENCIAL(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  double pi    = 3.14159265358979323846;
  double twopi = 2*pi;
  double mu0   = 4e-7*pi; // [N/A^2 = henry/m]
  double I  = mc.Ip(*s)*mu0/twopi;
  double L  = mc.NeoPolarization(*s);
  double R2 = mc.R2avrg(*s);
  double B_2= mc.B_2avrg(*s);
  double wrk = 1+I*I*L/( R2 - I*I*B_2 );
  return 1/wrk;
}

//*********************************************************************
// <(1-b)^0.5>; b = B(s,th,phi)/Bmax(s)
extern "C" WINDLLEXPORT double STDCALL MCDB12AVRG(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.db12Avrg(*s);
}

//*********************************************************************
// <(1-b)^1.5>
extern "C" WINDLLEXPORT double STDCALL MCDB32AVRG(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.db32Avrg(*s);
}

//*********************************************************************
// <b(1-b)^0.5>
extern "C" WINDLLEXPORT double STDCALL MCBDB12AVRG(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.db12Avrg(*s)-mc.db32Avrg(*s);;
}

//*********************************************************************
// flux average of (b/sqrt(1-x*b)) on the surface
// f(s,x) = <b/sqrt(1-x*b)>, where 0<=x<=1, 0<=s<=1, and b=B(s,th,ph)/Bmax(s)
extern "C" WINDLLEXPORT double STDCALL MCROPARINVAVRG(C3dMesh **mConf,double *s,double *x)
{
  C3dMesh &mc = **mConf;
  return mc.roParInvAvrg(*s,*x);
}

//*********************************************************************
// <p> -- flux average of pitch angle on the surface
//   p_avrg(s,x) = <sqrt(1-x*B(s,th,ph)/Bmax(s))>, where 0<=x<=1, 0<=s<=1,
// and definition of x is the following x=(1-p^2)/b, b=B/Bmax
extern "C" WINDLLEXPORT double STDCALL MCPAVRG(C3dMesh **mConf,double *s,double *x)
{
  C3dMesh &mc = **mConf;
  return mc.pavrg(*s,*x);
}

//****************************************************************************
// @return Integral(from x to 1) dx'/p_avrg(s,x') ,
//         where p_avrg(s,x) = <sqrt(1-x*B(s,th,ph)/Bmax(s))> 0<=x<=1, 0<=s<=1
extern "C" WINDLLEXPORT double STDCALL MCINTINVPAVRG(C3dMesh **mConf,double *s,double *x)
{
  C3dMesh &mc = **mConf;
  return mc.integralInvPavrg(*s,*x);
}

#if 0 
//*********************************************************************
// fabs(Jacobian*B^2) on surface s for (s,theta,phi)-coordinates
extern "C" WINDLLEXPORT double STDCALL MCJACOBIANB2(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.JacobianB2(*s);
}
#endif

//*********************************************************************
// iota(s)
extern "C" WINDLLEXPORT double STDCALL MCIOTA(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.iota(*s);
}

//*********************************************************************
// iotaPrime(s)== diota/ds
extern "C" WINDLLEXPORT double STDCALL MCIOTAPRIME(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.iotaPrime(*s);
}

//*********************************************************************
// PressurePrime(s)== dP/ds
extern "C" WINDLLEXPORT double STDCALL MCPRESSUREPRIME(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.pp(*s);
}

//*********************************************************************
// Pressure(s) in Pa
extern "C" WINDLLEXPORT double STDCALL MCPRESSURE(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.pressure(*s);
}

//*********************************************************************
// @return normalized poloidal flux
extern "C" WINDLLEXPORT double STDCALL MCTORFLUX2POLFLUX(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.sPol(*s);
}

//*********************************************************************
// R0(s)
extern "C" WINDLLEXPORT double STDCALL MCR0(C3dMesh **mConf)
{
  C3dMesh &mc = **mConf;
  return mc.R0();
}

//*********************************************************************
// Rmn(s)
extern "C" WINDLLEXPORT double STDCALL MCRMN(C3dMesh **mConf, int *m, int *n, double *s)
{
  C3dMesh &mc = **mConf;
  return mc.Rmn_sp(*s,*m,*n);
}

//*********************************************************************
// Zmn(s)
extern "C" WINDLLEXPORT double STDCALL MCZMN(C3dMesh **mConf, int *m, int *n, double *s)
{
  C3dMesh &mc = **mConf;
  return mc.Zmn_sp(*s,*m,*n);
}

//*********************************************************************
// Phimn(s)
extern "C" WINDLLEXPORT double STDCALL MCFMN(C3dMesh **mConf, int *m, int *n, double *s)
{
  C3dMesh &mc = **mConf;
  return mc.Phimn_sp(*s,*m,*n);
}

//*********************************************************************
// Bmn(s)
extern "C" WINDLLEXPORT double STDCALL MCBMN(C3dMesh **mConf, int *m, int *n, double *s)
{
  C3dMesh &mc = **mConf;
  return mc.Bmn_sp(*s,*m,*n);
}

//*********************************************************************
// BoozerBmn(s)
extern "C" WINDLLEXPORT double STDCALL MCBOOZERBMN(C3dMesh **mConf, int *m, int *n, double *s)
{
  C3dMesh &mc = **mConf;
  return mc.BoozerBmn(*s,*m,*n);
}

//*********************************************************************
// The method returns the average area of the cross section of the surface \e s
// @param s is the normalized toroidal flux.
extern "C" WINDLLEXPORT double STDCALL MCAVERAGECROSSSECTIONAREA(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  return mc.averageCrossSectionArea(*s);
}

//*********************************************************************
// rminor(s)
extern "C" WINDLLEXPORT double STDCALL MCRMINOR(C3dMesh **mConf)
{
  C3dMesh &mc = **mConf;
  return mc.rminor();
}

//*********************************************************************
// rminor(s)
extern "C" WINDLLEXPORT int STDCALL MCNPERIODS(C3dMesh **mConf)
{
  C3dMesh &mc = **mConf;
  return mc.nPeriods ();
}

//*********************************************************************
// Find intersection of the Ray with the last surface.
// The ray is R(t)=r0+rd*t  with t>0
// The origin of ray must be outside of the last surface
//
// INPUT:
//   Real *8 R1(3)  -- the first(origin) and
//   Real *8 R2(3)  --   the second point of the ray (cartesian coord.)
// OUTPUT:
//   Real *8 Entry(3) -- entry point of the ray into plasma (cartesian coord.)
//   Real *8 Exit(3)  -- exit  point of the ray from plasma (cartesian coord.)
//
/// The method traces the Last Closed Magnetic Surface (LCMS).
/// The method finds intersections of the ray with the LCMS, where the
/// ray is R(t)=r1+(r2-r1)/|r2-r1|*t,  and t>=0.
/// @param mConf is the address of the variable which contains the address of C3dMesh-object
/// @param R1 is the cartesian coordinates of the first point of the ray, must be outside with respect to the LCMS.
/// @param R2 is the cartesian coordinates of the second point of the ray
/// @param Entry this vector at return contains the coordinates of the entry point of the ray into plasma
/// @param Exit  this vector at return contains the coordinates of the exit point
/// @return
/// \li -2 -- the origin r1 of the ray lies inside LCMS;
/// \li -1 -- the ray does not enters into plasma;
/// \li  0 -- segment r1--r2 is outside plasma;
/// \li  1 -- only entry point lies on the segment r1--r2, i.e. r2 is inside plasma;
/// \li  2 -- entry and exit points lie on the segment r1--r2, i.e. r1 and r2 are outside plasma;
///
///  if returned value >= 0 then entry and exit points of the ray are in \e entrypoint and \e exitpoint
extern "C" WINDLLEXPORT int STDCALL MCTRACELCMS(C3dMesh **mConf,double *R1,double *R2,
                                                double *Entry,double *Exit)
{
  C3dMesh &mc = **mConf;
  Vector3d r1(R1),r2(R2); // 1st and 2nd points of the ray
  Vector3d r0(R1),rd(R2); // origin and direction of the ray
  rd  = r2-r0;            // direction of the ray
  rd /= rd.abs();         // normalize; the ray is r(t)=r0+rd*t  with t>0
  Vector3d ent(0.,0.,0.), ex(0.,0.,0.); // zero output
  ent.copyTo(Entry);
  ex.copyTo(Exit);

  if(!mc.isOK()) return -10;
  if(mc.xyzIsInside(r0)) return -2; // if TRUE then point r0 lies inside LCMS

  double dphi = 1*degree; // change this if you need more accurate calculations

  mc.traceLast(r0,rd,dphi,dphi);

  int size     = mc.getTraceSize();
  const Vector3d *R  = mc.getTraceRxyz();
  if(size==0) return -1; // ray does not enter in plasma

  ent = R[0];  // point where ray enters into plasma
  ex  = R[1];  // point where ray exits  plasma
  ent.copyTo(Entry);
  ex.copyTo(Exit);

  double r21 = (r2-r1).abs2();
  double e1 = (ent-r1).abs2();
  double x1 = (ex -r1).abs2();
//  CRayTrace::freeTrace(); // ?
  if(e1> r21) return 0;  // segment r1--r2 is outside plasma
  if(x1> r21) return 1;  // only entry point lies on the segment r1--r2, i.e. r2 is inside plasma
  if(x1<=r21) return 2;  // entry and exit points lie on the segment r1--r2
  return 0;
}

//*********************************************************************
// Find entry point of the Ray into the plasma(last closed surface)
// The ray is R(t)=r0+rd*t  with t>0
// It is assumed, that the origin of ray lies outside of the last surface
//
// INPUT:
//   r0, rd  -- origin and direction of the ray in cartesian coordinates
// OUTPUT:
//  entryPoint
//  @return  if the method succeeds, the return value is true; false otherwise.
extern "C" WINDLLEXPORT int STDCALL MCM3DGETRAYENTRYPOINT(C3dMesh **mConf, 
                                                          double *r0,double *rd, double *entryPoint)
{
  C3dMesh &mc = **mConf;
  Vector3d R0(r0),Rd(rd); // origin and direction of the ray
  Vector3d entry;
  bool b = mc.M3DgetRayEntryPoint(R0,Rd,entry);
  entry.copyTo(entryPoint);
  return b?True:False;
}

//****************************************************************************
// create border 'slast' of last contour at cylindrical angle 'ficyl'
// borderSize -- # of point in contour
// ficyl -- cylindrical angle
//  function returns 0 if no memory
extern "C" WINDLLEXPORT int STDCALL MCCREATEBORDER(C3dMesh **mConf,double *slast,double *ficyl,int *BrdSize)
{
  C3dMesh &mc = **mConf;
  bool b = mc.createBorder(*slast,*ficyl,*BrdSize);
  return b?True:False;
}

//****************************************************************************
// MCCREATEBORDER must be callled first!
// INPUT:
//  REAL*8 XYZ(3) -- point in cartesian coordinates
// OUTPUT:
//  function returns 'TRUE' if the point lies inside Border
//
extern "C" WINDLLEXPORT int STDCALL MCISINSIDEBORDER(C3dMesh **mConf,double *XYZ)
{
  C3dMesh &mc = **mConf;
  Vector3d xyz(XYZ);
  double rc = sqrt(xyz.x()*xyz.x()+xyz.y()*xyz.y());
  double fi = mc.mod2pi(atan2(xyz.y(),xyz.x()));
  double zc = xyz.z();
  Vector3d cyl(rc,fi,zc);
  bool b = mc.cylIsInsideBorder(cyl);
  return b?True:False;
}

//*********************************************************************
// Get B and derivatives using
// Program returns B and derivatives dB/dr, dB/dfi/r, dB/dz in point cyl
//  Input:
//    Real *8 cyl(3) -- cylindrical coordinates
//  Output:
//    Real *8 B(3) -- cylindrical coordinates of B
//    Real *8 dBdr(3)   -- dB/dr, where dBdr(1)=dBr/dr, dBdr(2)=dBfi/dr, dBdr(3)=dBz/dr
//    Real *8 dBdfir(3) -- dB/dfi/r
//    Real *8 dBdz(3)   -- dB/dz
//
extern "C" WINDLLEXPORT void STDCALL MCGETDBCYL(C3dMesh **mConf,double *cyl,
                                                double *B,double *dBdr,double *dBdfir,double *dBdz)
{
  C3dMesh &mc = **mConf;
  Vector3d c(cyl);
  Vector3d b,dbdr, dbdfir, dbdz, grads;
  double s = mc.getdBGradscyl(c,b,dbdr,dbdfir,dbdz,grads);
  b.copyTo(B);
  dbdr.copyTo(dBdr);
  dbdfir.copyTo(dBdfir);
  dbdz.copyTo(dBdz);
}

/// The method returns magnetic coordinates of the point of a field line, 
/// that passes through the given point (s0,theta0,ficyl0). 
/// The method  calculates function
///  \f[\theta(\varphi)=\theta_{0}+\iota(s_{0})(\varphi-\varphi_{0}) - [\lambda(s_{0},\theta,\varphi)-\lambda(s_{0},\theta_{0},\varphi_{0})]\f]
/// where \f$\lambda\f$ is non zero for VMEC coordinates.
/// @param mConf is the address of the variable which contains the address of C3dMesh-object
/// @param[in] mixCoord0 is the point of B-line in mixed coordinates,
///     where s and theta coincide with correspondent magnetic coordinates, 
///     but toroidal angle is cylindricyl angle
/// @param[in]  fcyl is the cylindricyl angle
/// @param[out] magCoord is the  magCoord(s0,theta(ficyl),phi(ficyl)) 
extern "C" WINDLLEXPORT void STDCALL MCCOORDFIELDLINE(C3dMesh **mConf,double *mixCoord0,
                                                      double *fcyl,double *magCoord)
{
  C3dMesh &mc = **mConf;
  Vector3d mxCoord(mixCoord0);
  Vector3d mag = mc.mixCoordFieldLine(mxCoord, *fcyl);
  mag.copyTo(magCoord);
}

//****************************************************************************
//****************************************************************************
//****************************************************************************
//The function creates 3d-mesh.
// The function tabulates the magnetic field, the flux surface label, and grad(s) in cylindrical coordinates.
extern "C" WINDLLEXPORT int STDCALL MCCREATEMESH(C3dMesh **mConf,double *dr,double *dz,double *dfi)
{
  C3dMesh &mc = **mConf;
  bool b = mc.createMesh(*dr,*dz,*dfi);
  return b?True:False;
}

//The function creates 3d-mesh.
// The function tabulates the magnetic field, the flux surface label, and grad(s) in cylindrical coordinates.
// Stellarator symmetry is used
extern "C" WINDLLEXPORT int STDCALL MCCREATEMESHUSINGSYMMETRY(C3dMesh **mConf,
                                                              double *dr,double *dz,double *dfi)
{
  C3dMesh &mc = **mConf;
  bool b = mc.createMeshUsingSymmetry(*dr,*dz,*dfi);
  return b?True:False;
}

//*********************************************************************
// get B and derivatives dB/dx,dB/dy,dB/dz in point xyz
//  Input:
//    Real *8 xyz(3) -- cartesian coordinates
//  Output:
//    Real *8 B(3) -- cartesian coordinates of B
//    Real *8 dBdx(3) -- dB/dx, where dBdx(1)=dBx/dx, dBdx(2)=dBy/dx, dBdx(3)=dBx/dx
//    Real *8 dBdy(3) -- dB/dy
//    Real *8 dBdz(3) -- dB/dz
//
extern "C" WINDLLEXPORT void STDCALL MCM3DGETDBXYZ(C3dMesh **mConf,double *xyz,
                                                   double *B,double *dBdx,double *dBdy,double *dBdz)
{
  C3dMesh &mc = **mConf;
  Vector3d R(xyz);
  Vector3d dbdx, dbdy, dbdz;

  Vector3d b=mc.M3DgetdBxyz(R,dbdx,dbdy,dbdz);
  b.copyTo(B);
  dbdx.copyTo(dBdx);
  dbdy.copyTo(dBdy);
  dbdz.copyTo(dBdz);
}

//*********************************************************************
// Get B and derivatives using 3d mesh interpolation
// Program returns B and derivatives dB/dr, dB/dfi/r, dB/dz in point cyl
//  Input:
//    Real *8 cyl(3) -- cylindrical coordinates
//  Output:
//    Real *8 B(3) -- cylindrical coordinates of B
//    Real *8 dBdr(3)   -- dB/dr, where dBdr(1)=dBr/dr, dBdr(2)=dBfi/dr, dBdr(3)=dBz/dr
//    Real *8 dBdfir(3) -- dB/dfi/r
//    Real *8 dBdz(3)   -- dB/dz
//
extern "C" WINDLLEXPORT void STDCALL MCM3DGETDBCYL(C3dMesh **mConf,double *cyl,double *B,
                                                   double *dBdr,double *dBdfir,double *dBdz)
{
  C3dMesh &mc = **mConf;
  Vector3d R(cyl);
  Vector3d dbdr, dbdfir, dbdz;

  Vector3d b=mc.M3DgetdBcyl(R,dbdr,dbdfir,dbdz);
  b.copyTo(B);
  dbdr.copyTo(dBdr);
  dbdfir.copyTo(dBdfir);
  dbdz.copyTo(dBdz);
}


//*********************************************************************
// Coordinate transformation using 3d mesh interpolation
// Find label of flux surface
// Input :
//  xyz -- point in cartesian coordinates
// Output:
// function returns label of flux surface -- normalized flux 's'
//
extern "C" WINDLLEXPORT double STDCALL MCM3DXYZ2S(C3dMesh **mConf,double *xyz)
{
  C3dMesh &mc = **mConf;
  Vector3d c(xyz);
  return mc.M3Dxyz2s(c);
}

//*********************************************************************
// Coordinate transformation using 3d mesh interpolation
// Find label of flux surface
// Input :
//  cyl -- point in cartesian coordinates
// Output:
// function returns label of flux surface -- normalized flux 's'
extern "C" WINDLLEXPORT double STDCALL MCM3DCYL2S(C3dMesh **mConf,double *cyl)
{
  Vector3d c(cyl);
  C3dMesh &mc = **mConf;
  return mc.M3Dcyl2s(c);
}

//*********************************************************************
//   function retuns grad(s) = Vector3d(ds/dr, ds/dfi/r, ds/dz)
//  Input:
//    Real *8 cyl(3) -- cylindrical coordinates
//  Output:
//    Real *8 grads(3) -- Vector3d(ds/dr, ds/dfi/r, ds/dz)
extern "C" WINDLLEXPORT void STDCALL MCM3DGETGRADS(C3dMesh **mConf,double *cyl,double *grads)
{
  C3dMesh &mc = **mConf;
  Vector3d R(cyl);
  Vector3d g=mc.M3DgetGrads(R);
  g.copyTo(grads);
}

//*********************************************************************
//   function retuns grad(s) in cartesian coordinates
//  Input:
//    Real *8 xyz(3) -- cartesian coordinates
//  Output:
//    Real *8 grads(3) --  grads in cartesian coordinates
extern "C" WINDLLEXPORT void STDCALL MCM3DGETGRADSXYZ(C3dMesh **mConf,double *xyz,double *grads)
{
  C3dMesh &mc = **mConf;
  Vector3d R(xyz);
  Vector3d g=mc.M3DgetGradsxyz(R);
  g.copyTo(grads);
}

//*********************************************************************
//   function retuns grad(s) in cartesian coordinates
//  Input:
//    Real *8 xyz(3) -- cartesian coordinates
//  Output:
//    Real *8 grads(3) --  grads in cartesian coordinates
extern "C" WINDLLEXPORT void STDCALL MCGETGRADSXYZ(C3dMesh **mConf,double *xyz,double *grads)
{
  C3dMesh &mc = **mConf;
  Vector3d R(xyz);
  double s=mc.xyz2s(R);
  Vector3d g = mc.getGradsxyz();
  g.copyTo(grads);
}
//****************************************************************************
// INPUT:
//  REAL*8 XYZ(3) -- point in cartesian coordinates
// OUTPUT:
//  function returns 'TRUE' if the point lies inside LCMS
extern "C" WINDLLEXPORT int STDCALL MCXYZISINSIDE(C3dMesh **mConf,double *XYZ)
{
  C3dMesh &mc = **mConf;
  Vector3d xyz(XYZ);
  bool b = mc.xyzIsInside(xyz);
  return b?True:False;
}

//****************************************************************************
// INPUT:
//  REAL*8 cyl(3) -- point in cartesian coordinates
// OUTPUT:
//  function returns 'TRUE' if the point lies inside LCMS
extern "C" WINDLLEXPORT int STDCALL MCCYLISINSIDE(C3dMesh **mConf,double *cyl)
{
  C3dMesh &mc = **mConf;
  Vector3d c(cyl);
  bool b = mc.cylIsInside(c);
  return b?True:False;
}

//*********************************************************************
// save mesh file
// *.bin4; *.bin8 are supported
extern "C" WINDLLEXPORT int STDCALL MCWRITEMESH(C3dMesh **mConf,char * fname,int i)
{
  C3dMesh &mc = **mConf;
  char *fn = new char[i+2];
  strtrim(fn,fname,i);
  bool ok= mc.writeMesh(fn);
  if(!ok) cerr << "mcWriteMesh: writing error\n";
  delete fn;
  return ok?True:False;
}


//*********************************************************************
// save magnetic coordinate file
// *.bin4; *.bin8 are supported
extern "C" WINDLLEXPORT int STDCALL MCWRITE(C3dMesh **mConf,char * fname,int i)
{
  C3dMesh &mc = **mConf;
  char *fn = new char[i+2];
  strtrim(fn,fname,i);
  bool ok= mc.write(fn);
  if(!ok) cerr << "mcWrite: writing error\n";
  delete fn;
  return ok?True:False;
}

//*********************************************************************
// get B and derivatives dB/dx,dB/dy,dB/dz in point xyz
//  Input:
//    Real *8 xyz(3) -- cartesian coordinates
//  Output:
//    Real *8 s (normalized toroidal flux), which normally satisfies \b 0<=s<=1;
//            s>1 means that point lies outside from LCMS; 
//    Real *8 B(3) -- cartesian coordinates of B
//    Real *8 dB(9)
//
//    dB(1:3) -- dB/dx, where dB(1)=dBx/dx, dB(2)=dBy/dx, dB(3)=dBz/dx
//    dB(4:6) -- dB/dy, where dB(4)=dBx/dy, dB(5)=dBy/dy, dB(6)=dBz/dy
//    dB(7:9) -- dB/dz, where dB(7)=dBx/dz, dB(8)=dBy/dz, dB(9)=dBz/dz
// 
extern "C" WINDLLEXPORT void STDCALL MCGETALLXYZ(C3dMesh **mConf,double *xyz,double *s,
                                                    double *B,double *gradS,double *dB)
{
  C3dMesh &mc = **mConf;
  Vector3d R(xyz);
  Vector3d b,dbdx,dbdy,dbdz,gs;

  *s = mc.getdBGradsxyz(R,b,dbdx,dbdy,dbdz,gs);
  b.copyTo(B);
  gs.copyTo(gradS);
//#define DB(row,col)  dB[(col)*3 + (row)] // Fortran77 style addressing
  dbdx.copyTo(dB);  // dB(i,j) = dB_i/dx_j
  dbdy.copyTo(dB+3);
  dbdz.copyTo(dB+6);
}

//*********************************************************************
// get B and derivatives dB/dx,dB/dy,dB/dz in point xyz
//  Input:
//    Real *8 xyz(3) -- cartesian coordinates
//  Output:
//    Real *8 B(3) -- cartesian coordinates of B
//    Real *8 dB(9)
//    dB(1:3) -- dB/dx, where dB(1)=dBx/dx, dB(2)=dBy/dx, dB(3)=dBz/dx
//    dB(4:6) -- dB/dy, where dB(4)=dBx/dy, dB(5)=dBy/dy, dB(6)=dBz/dy
//    dB(7:9) -- dB/dz, where dB(7)=dBx/dz, dB(8)=dBy/dz, dB(9)=dBz/dz
//
/// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
/// <b>1<s<1.5</b> means that point lies outside LCMS in the extrapolated region; \n
/// \b s>1000 means that point lies far outside from LCMS; \n
///
extern "C" WINDLLEXPORT void STDCALL MCM3DGETALLXYZ(C3dMesh **mConf,double *xyz,double *S,
                                                    double *B,double *gradS,double *dB)
{
  C3dMesh &mc = **mConf;
  Vector3d R(xyz);
  Vector3d b,dbdx,dbdy,dbdz,gs;

  *S = mc.M3DgetdBGradsxyz(R,b,dbdx,dbdy,dbdz,gs);
  b.copyTo(B);
  gs.copyTo(gradS);
//#define DB(row,col)  dB[(col)*3 + (row)] // Fortran77 style addressing
  dbdx.copyTo(dB);  // dB(i,j) = dB_i/dx_j
  dbdy.copyTo(dB+3);
  dbdz.copyTo(dB+6);
}

//*********************************************************************
// Coordinate transformation using 3d mesh interpolation
// Find label of flux surface and B
// Input :
//  xyz -- point in cartesian coordinates
// Output:
// function returns label of flux surface -- normalized flux 's'
//
extern "C" WINDLLEXPORT void STDCALL MCM3DGETSBXYZ(C3dMesh **mConf,double *xyz,double *S,double *B)
{
  C3dMesh &mc = **mConf;
  Vector3d R(xyz);

  *S = mc.M3Dxyz2s(R);
  (mc.M3DgetBxyz(R)).copyTo(B);
}

//*********************************************************************
/// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
/// <b>1<s<1.5</b> means that point lies outside LCMS in the extrapolated region; \n
/// \b s>1000 means that point lies far outside from LCMS; \n
///
extern "C" WINDLLEXPORT void STDCALL MCM3DGETALLCYL(C3dMesh **mConf,double *cyl,double *S,
                                                    double *B,double *gradS,double *gradB)
{
  C3dMesh &mc = **mConf;
  Vector3d R(cyl);
  Vector3d b,dbdr,dbdfir,dbdz,gs;

  *S = mc.M3DgetdBGradscyl(R,b,dbdr,dbdfir,dbdz,gs);
  b.copyTo(B);
  gs.copyTo(gradS);
  dbdr.copyTo(gradB);  // dB(i,j) = dB_i/dx_j
  dbdfir.copyTo(gradB+3);
  dbdz.copyTo(gradB+6);
}

//*********************************************************************
// Coordinate transformation using 3d mesh interpolation
// Find label of flux surface and B
// Input :
//  xyz -- point in cartesian coordinates
// Output:
// function returns label of flux surface -- normalized flux 's'
//
extern "C" WINDLLEXPORT void STDCALL MCM3DGETSBCYL(C3dMesh **mConf,double *cyl,double *S,double *B)
{
  C3dMesh &mc = **mConf;
  Vector3d R(cyl);
  *S = mc.M3Dcyl2s(R);
  if(*S<0) *S = mc.cyl2s(R);
  (mc.M3DgetBcyl(R)).copyTo(B);
}


//*********************************************************************
//
extern "C" WINDLLEXPORT void STDCALL MCGETBMINMAX(C3dMesh **mConf,double *s,double *Bmin,double *Bmax )
{
  C3dMesh &mc = **mConf;
  if(*s>=0) {
    *Bmin=mc.Bmin(*s);
    *Bmax=mc.Bmax(*s);
  }
  else {
    *Bmin=mc.Bmin();  // global minimum
    *Bmax=mc.Bmax();
  }
}

//*********************************************************************
//
extern "C" WINDLLEXPORT void STDCALL MCGETBMINMAXGLOBAL(C3dMesh **mConf,double *Bmin,double *Bmax )
{
  C3dMesh &mc = **mConf;
  *Bmin=mc.Bmin();  // global minimum
  *Bmax=mc.Bmax();
}

//*********************************************************************
//
extern "C" WINDLLEXPORT double STDCALL MCBMIN(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  if(*s>=0) return mc.Bmin(*s);
  else      return mc.Bmin();
}

//*********************************************************************
//
extern "C" WINDLLEXPORT double STDCALL MCBMAX(C3dMesh **mConf,double *s)
{
  C3dMesh &mc = **mConf;
  if(*s>=0) return mc.Bmax(*s);
  else      return mc.Bmax();
}


//*********************************************************************
//
extern "C" WINDLLEXPORT double STDCALL MCRAXIS(C3dMesh **mConf,double *ficyl)
{
  C3dMesh &mc = **mConf;
  return (mc.mixcoord2cyl(0,0,*ficyl))[0];
}

//*********************************************************************
//
extern "C" WINDLLEXPORT void STDCALL MCMAG2CYL(C3dMesh **mConf,double *magCoord,double *cyl)
{
  C3dMesh &mc = **mConf;
  Vector3d mag(magCoord);
  Vector3d c = mc.mag2cyl(mag);
  c.copyTo(cyl);               // copy to array cyl
}

//*********************************************************************
//
extern "C" WINDLLEXPORT void STDCALL MCMAG2XYZ(C3dMesh **mConf,double *magCoord,double *cyl)
{
  C3dMesh &mc = **mConf;
  Vector3d mag(magCoord);
  Vector3d c = mc.mag2xyz(mag);
  c.copyTo(cyl);               // copy to array cyl
}

//*********************************************************************
// input:
//    mixCoord = (s,theta,ficyl) coincides with magnetic for VMEC and tokamak symmetry flux coordinates
//    mixCoord = (s,theta,ficyl) s,theta coincide with boozer coordinates for bc-files
//     ficyl is the cylindrical angle
// output:
//    cyl  cylindrical coordinates
extern "C" WINDLLEXPORT void STDCALL MCMIXCOORD2CYL(C3dMesh **mConf,double *mixCoord,double *cyl)
{
  C3dMesh &mc = **mConf;
  Vector3d mCoord(mixCoord);
  Vector3d c = mc.mixcoord2cyl(mCoord);
  c.copyTo(cyl);               // copy to array cyl
}

//*********************************************************************
// input:
//    mixCoord = (s,theta,ficyl) coincides with magnetic for VMEC and tokamak symmetry flux coordinates
//    mixCoord = (s,theta,ficyl) s,theta coincide with boozer coordinates for bc-files
// output:
//    xyz  cartesian coordinates
extern "C" WINDLLEXPORT void STDCALL MCMIXCOORD2XYZ(C3dMesh **mConf,double *mixCoord,double *xyz)
{
  C3dMesh &mc = **mConf;
  Vector3d mCoord(mixCoord);
  Vector3d r = mc.mixcoord2xyz(mCoord);
  r.copyTo(xyz);               // copy to array xyz
}

//*********************************************************************
//  Return the sizes of a rectangle circumscribed about the cross section of the surface s at cylindrical angle fi.
extern "C" WINDLLEXPORT void STDCALL MCGETEXTENT(C3dMesh **mConf,double *s,double *fi,
                                                 double *Rmin,double *Rmax,double *Zmin,double *Zmax)
{
  C3dMesh &mc = **mConf;
  mc.getExtent(*Rmin,*Rmax,*Zmin,*Zmax,*s,*fi);
}

//*********************************************************************
//  Return the sizes of a rectangle circumscribed about the cross section of the LCMS.
extern "C" WINDLLEXPORT void STDCALL MCGETEXTENTLCMS(C3dMesh **mConf,double *Rmin,double *Rmax,
                                                     double *Zmin,double *Zmax)
{
  C3dMesh &mc = **mConf;
  mc.getExtentLCMS(*Rmin,*Rmax,*Zmin,*Zmax);
}

//*********************************************************************
//
extern "C" WINDLLEXPORT void STDCALL MCM3DSETSMAX(C3dMesh **mConf,double *sMax)
{
  C3dMesh &mc = **mConf;
  mc.setsMax(*sMax);
}

//*********************************************************************
//
extern "C" WINDLLEXPORT double STDCALL MCM3DGETSMAX(C3dMesh **mConf)
{
  C3dMesh &mc = **mConf;
  return mc.getsMax();
}

//*********************************************************************
/// The function transforms the current magnetic configuration to new volume or to new mesh on s.
/// @param[in]  mConf is the pointer to memory at which the address of C3dMesh-object is stored,
/// @param[in] *nsmax is the number of surfaces to store, if zero then the original flux surfaces are stored,
///   otherwise store the surfaces using a new mesh which is equidistant on sqrt(s).
/// @param[in] *eps is the truncation level to reduce the number of harmonics used in summation,
///  the only harmonics with m\<M and |n|\<N satisfying B<SUB>mn</SUB>/B<SUB>00</SUB> > eps 
///  are used in summation, eps<=1e-6 is a good choice.
/// @param[in] *appendWrite if true then the magnetic configuration is appended to a file 
///            set by next parameter,
/// @param[in] fname is the filename under which the configuration is saved.
/// @param[in] i is the strlen(fname), skip this parameter if this function is called from fortran.
/// @return  if the method succeeds, the return value is  1;  0 otherwise.

extern "C" WINDLLEXPORT int STDCALL MCSAVEREDUCED(C3dMesh **mConf, int *nsmax, double *eps, int *appendWrite, 
                                                  char * fname,int i)
{
  C3dMesh &mc = **mConf;
  char *fn = new char[i+2];
  strtrim(fn,fname,i);
   
  double  truncSav =mc.truncate();
  mc.truncate(*eps);

  bool append = appendWrite!=0;

  bool ok= mc.writeasciiReduced(fn, *nsmax,0,1,true,false,false,append);

  if(!ok) cerr << "mcSaveReduced: writing error\n";
  mc.truncate(truncSav);
  delete fn;
  return ok?True:False;
}

//*********************************************************************
/// The methods return the bounce-averade gradB drift velocity of trapped particles averaged over flux surface.
/// See eq(21,26,28) for Gw,Gs,Gv in Physics of Plasmas 12,112507(2005), http://link.aip.org/link/doi/10.1063/1.2131848
extern "C" WINDLLEXPORT void STDCALL MCGETGFACTORS(C3dMesh **mConf,double *s,double *Gw,double *Gs,double *Gv)
{
  C3dMesh &mc = **mConf;
  *Gw = mc.Gw(*s);
  *Gs = mc.Gs(*s);
  *Gv = mc.Gv(*s);
}

//*********************************************************************
/// The method returns the effective helical ripple for \f$1/\nu\f$ transport.  
/// \sa V.V.Nemov, S.V.Kasilov, W. Kernbichler, and M.F.Heyn, Physics of Plasmas, vol. 6, 4622(1999),
///   http://link.aip.org/link/PHPAEN/v6/i12/p4622/s1 
/// \sa 28th EPS, http://epsppd.epfl.ch/Madeira/html/pdf/P5.055.pdf , EGA Vol.25A (2001) 1985-1988
/// @param[in]  mConf is the pointer to memory at which the address of C3dMesh-object is stored,
/// @param[in]  s is the normalized toroidal flux
/// @param[out] epsEffGraz -- using Graz definition of the minor radius dr=ds/\<|grad(s)|\>, 
/// @param[out] epsEffDkes -- using DKES definition of minor radius r=sqrt(Flux(s)/(pi*B_00(s))), 
/// @param[out] epsEffVmec -- using VMEC definition of minor radius r=a*sqrt(s), 
extern "C" WINDLLEXPORT void STDCALL MCGETEPSEFF(C3dMesh **mConf,double *s,double *epsEffGraz,
                                                      double *epsEffDkes,double *epsEffVmec)
{
  C3dMesh &mc = **mConf;
  *epsEffGraz = mc.epsEffGraz(*s);
  *epsEffDkes = mc.epsEffDkes(*s);
  *epsEffVmec = mc.epsEffVmec(*s);
}


//*********************************************************************
/// The method sets parameters for calculating eps_eff and Gw, Gs, Gv.
/// @param xmin is the minimal value of sqrt(s), where s is the normalized toroidal flux,
/// @param xmax is the maximal value of sqrt(s),
/// @param maxPoints is the number of points for creating look-up table
extern "C" WINDLLEXPORT void STDCALL MCEPSEFFPARAMS(C3dMesh **mConf,double *xmin,double *xmax, int *maxPoints)
{
  C3dMesh &mc = **mConf;
  mc.epsEffSetXeffParam(*xmin, *xmax, *maxPoints);
}

//*********************************************************************
extern "C" WINDLLEXPORT int STDCALL MCISVMEC(C3dMesh **mConf)
{
  C3dMesh &mc = **mConf;
  return mc.isVmecCoordinates()?True:False;
}

//*********************************************************************
extern "C" WINDLLEXPORT int STDCALL MCISTOKPEST(C3dMesh **mConf)
{
  C3dMesh &mc = **mConf;
  return mc.isTokamakSymmetryFluxCoordinates()?True:False;
}

//*********************************************************************
//
extern "C" WINDLLEXPORT void STDCALL MCPEST2VMEC(C3dMesh **mConf,double *pest,double *vmec)
{
  C3dMesh &mc = **mConf;
  Vector3d p(pest);
  Vector3d v = mc.Pest2Vmec(p);
  v.copyTo(vmec);              
}

//****************************************************************************
// @param[in]  vmec is the  VMEC coordinates =(s,theta,phi)
// @param[out] pest  is the PEST coordinates
extern "C" WINDLLEXPORT void STDCALL MCVMEC2PEST(C3dMesh **mConf,double *vmec,double *pest)
{
  C3dMesh &mc = **mConf;
  Vector3d v(vmec);
  Vector3d p = mc.Vmec2Pest(v);
  p.copyTo(pest);              
}

//*********************************************************************
//
extern "C" WINDLLEXPORT void STDCALL MCBOOZER2VMEC(C3dMesh **mConf,double *boozer,double *vmec)
{
  C3dMesh &mc = **mConf;
  Vector3d b(boozer);
  Vector3d v = mc.Boozer2Vmec(b);
  v.copyTo(vmec);               // copy to array vmec
}

//****************************************************************************
// @param[in]  vmec is the  VMEC coordinates =(s,theta,phi)
// @param[out] boozer is the  Boozer coordinates
extern "C" WINDLLEXPORT void STDCALL MCVMEC2BOOZER(C3dMesh **mConf,double *vmec,double *boozer)
{
  C3dMesh &mc = **mConf;
  Vector3d v(vmec);
  Vector3d b = mc.Vmec2Boozer(v);
  b.copyTo(boozer);               // copy to array boozer
}

//****************************************************************************
extern "C" WINDLLEXPORT void STDCALL SETARGV(int n,char * opt, int i)
{
  char *op = new char[i+2];
  strtrim(op,opt,i);
  static std::vector<std::string> argv;
  std::string str=op;

  if(n==0) {
    argv.clear();
    argv.push_back(str);
  }
  else if(n>0) {
    argv.push_back(str);
  }
  else {
    const char **arg = new const char*[argv.size()];
    for(unsigned int i=0; i<argv.size(); i++) {
      arg[i] = argv[i].c_str();
    }
    // call c_main(argv.size(),arg);
  }
}


/// The method returns the parallel adiabatic invariant with the integral taken between two turning points according formula
/// \f$\int \sqrt{1-B/B_t} dl\f$ 
/// @param mag is the magnetic coordinates of the point 
/// from which integrating along field line is started and where \f$B_t=B(mag)\f$
/// @param dphi is the toroidal step in radians for fild line integration (default value is 0.4degree if dphi==0)
/// @param iotaDenominator if not zero then this number is
///        the denominator in iota aproximation by rational number.
/// @param[out] Jinv is the array with the parallel adiabatic invariant in the first element, dJ/ds in the second element, dJ/dalpha in the third element
/// @param[out] Tbounce is the  bounce time \f$\oint\frac{dl}{\sqrt{1-B/B_t}}\f$ , and 
/// @param[out] Bt  is the  \f$B_t=B(mag)\f$ 
/// \sa  getMagCoordAtB()
extern "C" WINDLLEXPORT void STDCALL MCJINVARIANT(C3dMesh **mConf,double *mag,double *dphi, int *iotaDenominator, 
                                                      double *Jinv, double *Tbounce, double *Bt)
{
  C3dMesh &mc = **mConf;
  Vector3d v(mag);
  std::vector<double> j = mc.Jinvariant(v,*dphi,*iotaDenominator);
  Jinv[0] = j[0];
  Jinv[1] = j[1];
  Jinv[2] = j[2];
  *Tbounce = j[3];
  *Bt      = j[4];
}

/// The method returns the magnetic coordinates of the point (that belongs to the given field line)
/// at which the magnetic field is equal to the given value.
/// @param mag defines the field line; 
/// @param Bt is the magnetic field value, must be within the interval Bmin..Bmax for the considered field line;  
/// @@param[out] magBt array with the magnetic coordinates at which the magnetic field value is Bt or
///            array with elements (-1,Bmin,Bmax) if Bt is not within the interval Bmin..Bmax;                 
/// \note the method seaches the point of interest starting from the point where B=Bmin;
  ///        Bmin is seached within one machine period (0<=phi<2*pi/Np) along given field line.
/// \sa Jinvariant(), Bmin(double s), Bmax(double s)
extern "C" WINDLLEXPORT void STDCALL MCGETMAGCOORDATB(C3dMesh **mConf,double *mag,double *Bt,double *magBt)
{
  C3dMesh &mc = **mConf;
  Vector3d v(mag);
  Vector3d m = mc.getMagCoordAtB(v, *Bt);
  m.copyTo(magBt);               // copy to array magBt
}

/// If yesNo>0 then the method forces to use Jacobian=dX/dFlux*(dX/dtheta^dX/dphi) for calculating B field.
extern "C" WINDLLEXPORT void STDCALL MCUSEMIXEDPRODUCTFORJACOBIAN(C3dMesh **mConf,int *yesNo)
{
  C3dMesh &mc = **mConf;
  mc.useMixedProductForJacobian(*yesNo>0?true:false);
}
