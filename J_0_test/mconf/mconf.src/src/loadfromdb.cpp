#include "../include/CStconfig.h"

#ifdef NO_MCDB

  namespace MConf {

  bool CStconfig::writeBoozerCoord(const void *mdb, int magneticgeometryrunlogid, int Nsmax,bool magAxis) const {return false;};
  bool CStconfig::loadFromDB  (const void *mcdb, int equilibrium_id) {return false;};
  bool CStconfig::deleteFromDB(const void *mcdb, int equilibrium_id) const {return false;};
  bool CStconfig::writeIntoDB (const void *mcdb, const char* confname,const char* machine_name,int Nsmax,bool magAxis) const {return false;};

  };   //namespace MConf

#else

#include "../../MCDB/MCAPI.h"
#include "../include/CStconfig.h"
#include <fstream>
#include <iostream>
#include <time.h>

namespace MConf {

#define SEARCH_FLUXSURFACE_INDEX_IN_FLUXSURFACETABLE     1

//****************************************************************************
static int sorts(const void *f1,const void *f2) {
  return (((FLUXSURFACE *)f1)->s < ((FLUXSURFACE *)f2)->s)?-1:1;
}

//****************************************************************************
static int sortid(const void *f1,const void *f2) {
  return ((FLUXSURFACE *)f1)->id - ((FLUXSURFACE *)f2)->id;
}

//****************************************************************************
static int is_bsearch(FLUXSURFACE * fxs, int nrec, int fx_id) {
  FLUXSURFACE fkey;
  fkey.id=fx_id;
  FLUXSURFACE *f=(FLUXSURFACE *)bsearch((void *)&fkey,(void *)fxs,
                     (size_t)nrec,sizeof(FLUXSURFACE),sortid);
  return int(f->s);
}

//****************************************************************************
// Read bc-file  from data base.
//
bool CStconfig::loadFromDB(const void *mdb, int equilibrium_id)
{
  MCDB *mcdb = (MCDB *)mdb;
  if(!mcdb->connected) return mOK=false;
  std::string confname;

  init();

  freeConf();

  int i,runlog_id;
  char condition[64];
  clock_t tStart = clock();  //timing

  {
    EQUILIBRIUM      *eq    = NULL;
    EQUILIBRIUMINPUT *eqinp = NULL;
    EXPERIMENTSETUP  *ex    = NULL;
    MACHINESETUP     *mch   = NULL;
    MAGNETICGEOMETRY *mg    = NULL;
    GLOBALQUANTITIES *gb    = NULL;
    mOK=false;
    sprintf(condition, "%s=%d",FN_EQ_ID,equilibrium_id);
    if(load_equilibrium(mcdb,&eq,condition)>0) {
      sprintf(condition, "%s=%d",FN_EI_ID,eq->equilibriuminput_id);
      if(load_equilibriuminput(mcdb,&eqinp,condition)>0) {
        sprintf(condition, "%s=%d",FN_ES_ID,eqinp->experimentsetup_id);
        if(load_experimentsetup(mcdb,&ex,condition)>0) {
          sprintf(condition, "%s=%d",FN_MS_ID,ex->machinesetup_id);
          if(load_machinesetup(mcdb,&mch,condition)>0) {
            mNp = mch->fieldperiods;
            sprintf(condition, "%s=%d",FN_MG_EQUILIBRIUM_ID, equilibrium_id);
            if(load_magneticgeometry(mcdb,&mg,condition)>0) {
              runlog_id = mg->runlog_id;
              sprintf(condition, "%s=%d",FN_GQ_RUNLOG_ID,runlog_id);
              if(load_globalquantities(mcdb,&gb,condition)>0) {
                mR0   = gb->majorradius;
                ma    = gb->minorradius;
                mFlux = gb->flux;
                mOK=true;
              }
            }
          }
        }
      }
    }
    if(mOK) confname = eq->configurationname;
    free(eq);
    free(eqinp);
    free(ex);
    free(mch);
    free(mg);
    free(gb);
    if(mOK==false)  return mOK;
  }

  FLUXSURFACEGEOMETRY *fxsg=NULL;
  FLUXSURFACE         *fxs =NULL;
  int mpol,ntor,nosurfaces;
  int fluxsurfacegeometry_id;
  mOK=false;
  sprintf(condition, "%s=%d",FN_FG_MAGNETICGEOMETRY_RUNLOG_ID, runlog_id);
  if(load_fluxsurfacegeometry(mcdb,&fxsg,condition)>0) {
    mpol  = fxsg->nopoloidalmodes;
    ntor  = fxsg->notoroidalmodes;
    nosurfaces = fxsg->nosurfaces;     // ????
    fluxsurfacegeometry_id = fxsg->id; // save for future use
    sprintf(condition, "%s=%d ORDER BY S",FN_FS_FLUXSURFACEGEOMETRY_ID,fluxsurfacegeometry_id);
    nosurfaces=load_fluxsurface(mcdb,&fxs,condition);
    if(nosurfaces>0) mOK=true; else free(fxs);
  }

  free(fxsg);
  if(mOK==false) return mOK=false;

  if(!resize(nosurfaces,mpol,ntor)) {free(fxs);return mOK=false;}

  for(i=0; i<nosurfaces; i++) { // nosurfaces must be equal to mNs0
    int is = m1stSindex+i;
    ms.rw()   [is]  = fxs[i].s;
    miota.rw()[is]  = fxs[i].iota;
    mIpol.rw()[is]  = fxs[i].currpol;
    mItor.rw()[is]  = fxs[i].currtor;
    mpprime.rw()[is]= fxs[i].pprime;
    mg00.rw() [is]  = fxs[i].sqrtg00;
    fxs[i].s   = is;        // we will need this later
  }
  if(ms.rw()[m1stSindex]==0.) ms.rw()[m1stSindex]=1e-8;  // TODO: fix this

#if SEARCH_FLUXSURFACE_INDEX_IN_FLUXSURFACETABLE != 0
  qsort((void *)fxs,(size_t)nosurfaces,sizeof(FLUXSURFACE),sortid); // sort on fxs[i].id
#endif

  FOURIERBLOCKS *fourierblocks=NULL;
  mOK=false;
  sprintf(condition, "%s=%d",FN_FB_FLUXSURFACEGEOMETRY_ID,fluxsurfacegeometry_id);
  int nrec=load_fourierblocks(mcdb,&fourierblocks,condition);
  if(nrec>0) {
    mOK=true;
    int i,j,k;
    for (i=0; i<nrec; i++) {
      k=fourierblocks[i].data_size/sizeof(FMODE);
      std::cerr <<"fourierblocks_size="<<fourierblocks[i].data_size <<"\n";
      FMODE *fmode = (FMODE*)fourierblocks[i].data;
      int is,fx_id=-1;
      for(j=0;j<k;j++) {
#if SEARCH_FLUXSURFACE_INDEX_IN_FLUXSURFACETABLE != 0
        if(fx_id!=fmode[j].fluxsurface_id) {    // if new surface
          fx_id = fmode[j].fluxsurface_id;      //  then
          is = is_bsearch(fxs,nosurfaces,fx_id);//    find its number
        }
#else
        is           = fmode[j].is + m1stSindex;
#endif
        int m        = fmode[j].m;
        int n        = fmode[j].n;
        Rmn(is,m,n)  = fmode[j].r_mn,
        Zmn(is,m,n)  = fmode[j].z_mn,
        Phimn(is,m,n)= fmode[j].phi_lam_mn,
        bmn(is,m,n)  = fmode[j].b_mn;
      }
    }
  }

  free(fxs);
  free_fourierblocks(fourierblocks,nrec);

  std::cerr <<"loadFromDB time="<<double(clock()-tStart)/CLOCKS_PER_SEC<<"\n";
  mOK = initAfterLoad(W7X_DB);
  mfname = confname;
  if(!mOK) {
    mfname = "no data loaded";
    freeConf();
  }
  return mOK;
}



//*************************************************************************
//
//
bool CStconfig::writeIntoDB(const void *mdb, const char* confname,const char* machine_name,int Nsmax,bool magAxis) const
{  /* complex table insert operation (with dummy content) test in dependency order */

  if(vmecCoordinates)  return false;

  MCDB *mcdb = (MCDB*)mdb;

  { //test whether confname exist
    EQUILIBRIUM *eq;
    char condition[512];
    sprintf(condition, "%s='%s'",FN_EQ_CONFIGURATIONNAME,confname);
    if(load_equilibrium(mcdb,&eq,condition)>0) {
      printf("'%s' already exist in DB\n",confname);
      free(eq);
      return false;
    }
  }

  int   periods = mNp;
  double majorradius = mR0;
  double minorradius = ma;

  int i,codeid,runlogid,profileid,
      machinesetupid,experimentsetupid,
      makexgridparameterid,vacuumfieldrunlogid,
      currentmappingid,equilibriuminputid,
      equilibriumid,magneticgeometryrunlogid,
      gridprofileid;

  CODE code;
  RUNLOG runlog;
  PROFILE profile;
  GRIDPROFILE gridprofile;
  MACHINESETUP machinesetup;
  EXPERIMENTSETUP experimentsetup;
  EQUILIBRIUMINPUT equilibriuminput;
  EQUILIBRIUM equilibrium;
  MAGNETICGEOMETRY magneticgeometry;
  MAKEXGRIDPARAMETER makexgridparameter;
  MGRID_VACUUMFIELD vacuumfield;
  CURRENTMAPPING currentmapping;
  GLOBALQUANTITIES globalquantities; //added by turkin

  // preparing/inserting code
 clear_code(&code);
 strcpy(code.codename,"mcviewer");
 strcpy(code.version,"1.1");
 codeid=insert_code(mcdb,&code);
 printf("%s %s=%d\n",TN_CO,FN_CO_ID,codeid);
 if(codeid<=0) return false;

 // preparing/inserting runlog
 time_t aclock;
 time( &aclock );
 struct tm *ltm = localtime( &aclock );
 clear_runlog(&runlog);
 runlog.code_id=codeid;
 strcpy(runlog.author,"Turkin");
 strcpy(runlog.compiler,"Visual C++ 2005");
 strcpy(runlog.host,"P4 3.06GHz");
 runlog.rundate.day=ltm->tm_mday;
 runlog.rundate.month=ltm->tm_mon;
 runlog.rundate.year=ltm->tm_year+1900;
 runlogid=insert_runlog(mcdb,&runlog);
 printf("%s %s=%d\n",TN_RL,FN_RL_ID,runlogid);
 if(runlogid<=0) return false;

 // preparing/inserting profile
 clear_profile(&profile);
 sprintf(profile.description, "dummydata%d",runlogid);
 sprintf(profile.physquantity,"dummydata%d",runlogid);
 profile.profiletype=1;
 profile.profilestoragetype=1;
 profileid=insert_profile(mcdb,&profile);
 printf("%s %s=%d\n",TN_PR,FN_PR_ID,profileid);
 if(profileid<=0) return false;

 // preparing/inserting grid profile
 for (i=0; i<1; i++) {
   clear_gridprofile(&gridprofile);
   gridprofile.profile_id=profileid;
   gridprofile.valuex=i;
   gridprofile.valuey=4.5 + i;
   gridprofileid=insert_gridprofile(mcdb,&gridprofile);
   printf("%s %s=%d\n",TN_GP,FN_GP_ID,gridprofileid);
   if(gridprofileid<=0) return false;
 }

 // preparing/inserting machine setup
 clear_machinesetup(&machinesetup);
 strcpy(machinesetup.machinename, machine_name);    //  strcpy(machinesetup.machinename,"W7X");
 sprintf(machinesetup.machineversion,"V%d",runlogid);
 machinesetup.fieldperiods=periods;
 machinesetup.radius = majorradius;
 machinesetupid=insert_machinesetup(mcdb,&machinesetup);
 printf("%s %s=%d\n",TN_MS,FN_MS_ID,machinesetupid);
 if(machinesetupid<=0) return false;

 // preparing/inserting experiment setup
 clear_experimentsetup(&experimentsetup);
 experimentsetup.machinesetup_id=machinesetupid;
 sprintf(experimentsetup.description, "dummydata%d",runlogid);
 experimentsetup.expsetupid=1;
 experimentsetupid=insert_experimentsetup(mcdb,&experimentsetup);
 printf("%s %s=%d\n",TN_ES,FN_ES_ID,experimentsetupid);
 if(experimentsetupid<=0) return false;

 // preparing/inserting makexgridparameter
 clear_makexgridparameter(&makexgridparameter);
 makexgridparameterid=insert_makexgridparameter(mcdb,&makexgridparameter);
 printf("%s %s=%d\n",TN_MX,FN_MX_ID,makexgridparameterid);
 if(makexgridparameterid<=0) return false;

 // preparing/inserting currentmapping
 clear_currentmapping(&currentmapping);
 currentmapping.machinesetup_id=machinesetupid;
 currentmappingid=insert_currentmapping(mcdb,&currentmapping);
 printf("%s %s=%d\n",TN_CM,FN_CM_ID,currentmappingid);
 if(currentmappingid<=0) return false;

 // preparing/inserting vacuumfield
 clear_mgrid_vacuumfield(&vacuumfield);
 vacuumfield.currentmapping_id=currentmappingid;
 vacuumfield.makexgridparameter_id=makexgridparameterid;
 vacuumfield.runlog_id=runlogid;
 vacuumfield.mgrid_size=16;
 vacuumfield.mgrid=(BLOBPTR)calloc(1,vacuumfield.mgrid_size);
 vacuumfieldrunlogid=insert_mgrid_vacuumfield(mcdb,&vacuumfield);
 printf("%s %s=%d\n",TN_MV,FN_MV_RUNLOG_ID,vacuumfieldrunlogid);
 free(vacuumfield.mgrid);
 if(vacuumfieldrunlogid<=0) return false;

 // preparing/inserting equilibrium input
 clear_equilibriuminput(&equilibriuminput);
 equilibriuminput.experimentsetup_id=experimentsetupid;
 equilibriuminput.vacuumfield_runlog_id=vacuumfieldrunlogid;
 equilibriuminput.inputprofile1_id=profileid;
 equilibriuminput.inputprofile2_id=profileid;
 equilibriuminputid=insert_equilibriuminput(mcdb,&equilibriuminput);
 printf("%s %s=%d\n",TN_EI,FN_EI_ID,equilibriuminputid);
 if(equilibriuminputid<=0) return false;

 // preparing/inserting equilibrium
 clear_equilibrium(&equilibrium);
 strcpy(equilibrium.configurationname,confname);
 sprintf(equilibrium.description, "Only Boozer coordinates are stored. Remeshed and saved by mcviewer. Truncation level Bmn/B0=%g\n",mmodBtol);
 equilibrium.equilibriuminput_id=equilibriuminputid;
 equilibriumid=insert_equilibrium(mcdb,&equilibrium);
 printf("%s %s=%d\n",TN_EQ,FN_EQ_ID,equilibriumid);
 if(equilibriumid<=0) return false;

 // preparing/inserting magnetic geometry
 clear_magneticgeometry(&magneticgeometry);
 magneticgeometry.equilibrium_id=equilibriumid;
 magneticgeometry.runlog_id=runlogid;
 magneticgeometryrunlogid=insert_magneticgeometry(mcdb,&magneticgeometry);
 printf("%s %s=%d\n",TN_MG,FN_MG_RUNLOG_ID,magneticgeometryrunlogid);
 if(magneticgeometryrunlogid<=0) return false;

 // preparing/inserting globalquantities      //added by turkin, see also Boozer.cpp
 clear_globalquantities(&globalquantities);
 globalquantities.runlog_id=runlogid;
 globalquantities.majorradius=majorradius;
 globalquantities.minorradius=minorradius;
 globalquantities.flux=1;
 i=insert_globalquantities(mcdb,&globalquantities);
 printf("%s %s=%d\n",TN_GQ,FN_GQ_RUNLOG_ID,i);
 if(i<=0) return false;

 return writeBoozerCoord(mdb,magneticgeometryrunlogid,Nsmax,magAxis);
}

#define USEFOURIERBLOCKS 1

//****************************************************************************
// Remesh and write Boozer Coord. with surfaces equidistantly distributed on r_eff.
// Nsmax is the number of surfaces.
// Don't remesh if Nsmax==0, use original s-mesh
//
bool CStconfig::writeBoozerCoord(const void *mdb, int magneticgeometryrunlogid, int Nsmax,bool magAxis) const
{
  MCDB *mcdb = (MCDB*)mdb;

  bool remesh = (Nsmax!=0)?true:false;

  if(!remesh) Nsmax=mNs0;
  Nsmax=abs(Nsmax);

  bool mNewFormatCC=true;

  int Mmax = mmin(1,mM);
  int Nmax = mmin(1,mN);
  int i,is,m,n;
  for(i=0; i<=mLCMSindex; i++) {
    Mmax = mmax(Mmax,mpM[i]);
    Nmax = mmax(Nmax,mpN[i]);
  }
  int Mmax0 = Mmax;
  int Nmax0 = Nmax;

  double Fluxmax = Flux(1.);
  double rmax = r(1.);

//1  fprintf(fp,"%4d %4d %4d %4d ",Mmax0,Nmax0,Nsmax,mNp);
//1  fprintf(fp,"%13.6E %8.5f %8.5f\n",Fluxmax,rmax,mR0);
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// insert now a FluxSurfaceGeometry record
  FLUXSURFACEGEOMETRY fluxsurfacegeometry;
  clear_fluxsurfacegeometry(&fluxsurfacegeometry);
  fluxsurfacegeometry.magneticgeometry_runlog_id=magneticgeometryrunlogid;
  fluxsurfacegeometry.nopoloidalmodes=Mmax0; // max# of poloidal harmonic  0 <= m <= M
  fluxsurfacegeometry.notoroidalmodes=Nmax0; // max# of toroidal harmonic -N <= n <= N
  fluxsurfacegeometry.nosurfaces=Nsmax;      //do we need this field? // # of flux surface

  int fluxsurfacegeometryid=insert_fluxsurfacegeometry(mcdb,&fluxsurfacegeometry);
  printf(" %s %s %d\n",TN_FG,FN_FG_ID,fluxsurfacegeometryid);

// insert now a globalquantities record

  GLOBALQUANTITIES *gb;
  char condition[64];
  sprintf(condition, "%s=%d",FN_GQ_RUNLOG_ID,magneticgeometryrunlogid);
  if(load_globalquantities(mcdb,&gb,condition)>0) {
    (*gb).majorradius=mR0;          // major radius [m]
    (*gb).minorradius=rmax;         // minor radius [m]
    (*gb).flux=Fluxmax;             // max Flux [Tm^2]
    int i=update_globalquantities(mcdb,gb,condition);
    if(i==1) printf("%s %s %d\n",TN_GQ,FN_GQ_RUNLOG_ID,magneticgeometryrunlogid);
    free(gb);
  }
  else {
    GLOBALQUANTITIES globalquantities;
    clear_globalquantities(&globalquantities);
    globalquantities.runlog_id=magneticgeometryrunlogid;
    globalquantities.majorradius=mR0;       // major radius [m]
    globalquantities.minorradius=rmax;      // minor radius [m]
    globalquantities.flux=Fluxmax;          // max Flux [Tm^2]
    int globalquantitiesid=insert_globalquantities(mcdb,&globalquantities);
    printf("%s %s %d\n",TN_GQ,FN_GQ_RUNLOG_ID,globalquantitiesid);
  }

  FOURIERBLOCKS fourierblocks;
  clear_fourierblocks(&fourierblocks);
  fourierblocks.fluxsurfacegeometry_id=fluxsurfacegeometryid;
  if(USEFOURIERBLOCKS) fourierblocks.data = (BLOBPTR) new FMODE[( (Nmax0+1)+Mmax0*(2*Nmax0+1) )*Nsmax];
  FMODE *FMode = (FMODE*)fourierblocks.data;
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  Mmax = mmin(1,mM);
  Nmax = mmin(1,mN);
  double s1st  = CStconfig::s1st;
  double sLast = CStconfig::sLast;
  if(magAxis) s1st = 0;

  double dr = (sqrt(sLast)-sqrt(s1st))/(Nsmax-1);
  for(i=0; i<Nsmax; i++) {  // loop over surfaces
    double r = sqrt(s1st) + i*dr;
    if(i==0&&magAxis) r = dr/1024;
    double s = (i==Nsmax-1)?sLast:(r*r);
    double iota_,Ip_,It_,pp_,g00_;
    if(remesh) {
      is = SearchSindx(s);    // perform a binary search of the segment in which s lies
      iota_=iota(s);
      Ip_  =Ip(s)/mNp;
      It_  =It(s);
      pp_  =const_cast<CStconfig*>(this)->pp(s);
      g00_ =g00(s);
    }
    else {
      is = m1stSindex+i;
      s    =ms[is];
      iota_=miota[is];
      Ip_  =mIpol[is]*mBnorm;
      It_  =mItor[is]*mBnorm;
      pp_  =mpprime[is];
      g00_ =mg00[is];
    }
//2    fprintf(fp,"%17.10E %16.8E %16.8E %16.8E %16.8E %16.8E\n", s,iota_,Ip_,It_,pp_,g00_);
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// insert now a FluxSurface record
    FLUXSURFACE fluxsurface;
    clear_fluxsurface(&fluxsurface);
    fluxsurface.fluxsurfacegeometry_id=fluxsurfacegeometryid;
    fluxsurface.s      =s;
    fluxsurface.iota   =iota_;
    fluxsurface.currpol=Ip_;
    fluxsurface.currtor=It_;
    fluxsurface.pprime =pp_;
    fluxsurface.sqrtg00=g00_;
    int fluxsurfaceid=insert_fluxsurface(mcdb,&fluxsurface);
    printf("surface %d",i);
    printf(" %s %s %d\n",TN_FS,FN_FS_ID,fluxsurfaceid);
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    int Mma=mpM[is];
    int Nma=mpN[is];
    Mmax = mmax(Mmax,Mma);    // # of harmonics must be monotonically increasing
    Nmax = mmax(Nmax,Nma);    //  function of surface number;
    Mmax = mmin(Mmax,Mmax0);  //  we need this to avoid ill-behaved spline;
    Nmax = mmin(Nmax,Nmax0);  //  see also calculation of spline coef. in addSurfaces()

    for(m=0; m<=Mmax; m++) {
      int lowN = (mNewFormatCC&&m==0)?0:-Nmax;
      for(n=lowN; n<= Nmax; n++) {
        double r,z,p,b;
        if(remesh) {
          r=Rmn_sp(s,m,n);
          z=Zmn_sp(s,m,n);
          p=Phimn_sp(s,m,n);
          b=Bmn_sp(s,m,n);
        }
        else {
          r=Rmn(is,m,n);
          z=Zmn(is,m,n);
          p=Phimn(is,m,n);
          b=Bmn(is,m,n);
        }
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if(r==0&&z==0&&p==0&&b==0) continue; // don't insert zero-harmonics into the DB
//3     fprintf(fp,"% 5d %4d ", m, n);
//3     fprintf(fp,"%16.8E %16.8E %16.8E %16.8E\n",r,z,p,b);
// insert now a FourierModes record
        if(USEFOURIERBLOCKS) {
          FMode->fluxsurface_id=fluxsurfaceid;
          FMode->m=m;   // poloidal harmonic number, 0 <= m <= M
          FMode->n=n;   // toroidal harmonic number,-N <= n <= N
          FMode->r_mn=r;
          FMode->z_mn=z;
          FMode->phi_lam_mn=p;
          FMode->b_mn=b;
          FMode++;
          fourierblocks.data_size += sizeof(FMODE);
        }
        else {
          FOURIERMODES fouriermodes;
          clear_fouriermodes(&fouriermodes);
          fouriermodes.fluxsurface_id=fluxsurfaceid;
          fouriermodes.m=m;
          fouriermodes.n=n;
          fouriermodes.r_mn=r;
          fouriermodes.z_mn=z;
          fouriermodes.phi_lam_mn=p;
          fouriermodes.b_mn=b;
          int fouriermodesid=insert_fouriermodes(mcdb,&fouriermodes);
          printf(" %s %s %d n",TN_FM,FN_FM_ID,fouriermodesid);
          }
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      }
    }
  }

  if(USEFOURIERBLOCKS) {
    int fourierblocksid=insert_fourierblocks(mcdb,&fourierblocks);
    printf(" %s %s %d\n",TN_FB,FN_FB_ID,fourierblocksid);
    delete [] fourierblocks.data;
  }

  return(true);
}


//****************************************************************************
// delete 'equilibrium_id' from MCDB
//
bool CStconfig::deleteFromDB(const void *mdb, int equilibrium_id) const
{
  MCDB *mcdb = (MCDB *)mdb;
  if(!mcdb->connected) return false;

  EQUILIBRIUM         *equilibrium        = NULL;
  EQUILIBRIUMINPUT    *equilibriuminput   = NULL;
  EXPERIMENTSETUP     *experimentsetup    = NULL;
  MAGNETICGEOMETRY    *magneticgeometry   = NULL;
  FLUXSURFACEGEOMETRY *fluxsurfacegeometry= NULL;
  FLUXSURFACE         *fluxsurface = NULL;
  RUNLOG              *runlog      = NULL;
  MGRID_VACUUMFIELD   *vacuumfield = NULL;

  int nr,nr1,runlog_id;
  char condition[64];

  sprintf(condition, "%s=%d",FN_EQ_ID,equilibrium_id);
  if(0==load_equilibrium(mcdb,&equilibrium,condition)) return false;

  sprintf(condition, "%s=%d",FN_MG_EQUILIBRIUM_ID, equilibrium_id);
  load_magneticgeometry(mcdb,&magneticgeometry,condition);
  runlog_id = magneticgeometry->runlog_id;

  sprintf(condition, "%s=%d",FN_FG_MAGNETICGEOMETRY_RUNLOG_ID, runlog_id);
  load_fluxsurfacegeometry(mcdb,&fluxsurfacegeometry,condition);
//-----------------------

  sprintf(condition, "%s=%d",FN_FB_FLUXSURFACEGEOMETRY_ID,fluxsurfacegeometry->id);
  if((nr=delete_fourierblocks(mcdb,condition))>0)
    printf("delete %s->%s\n",TN_FB,condition);

  if(nr<0) return false;

//-----------------------
  sprintf(condition, "%s=%d",FN_FS_FLUXSURFACEGEOMETRY_ID,fluxsurfacegeometry->id);
  nr=load_fluxsurface(mcdb,&fluxsurface,condition);
  while(nr--) {
    sprintf(condition, "%s=%d",FN_FM_FLUXSURFACE_ID,fluxsurface[nr].id);
    if(delete_fouriermodes(mcdb,condition)>0)
      printf("delete %s->%s\n",TN_FM,condition);
  }
  sprintf(condition, "%s=%d",FN_FS_FLUXSURFACEGEOMETRY_ID,fluxsurfacegeometry->id);
  if(delete_fluxsurface(mcdb,condition)>0)
    printf("delete %s->%s\n",TN_FS,condition);
//-----------------------
  sprintf(condition, "%s=%d",FN_FG_MAGNETICGEOMETRY_RUNLOG_ID, runlog_id);
  if(delete_fluxsurfacegeometry(mcdb,condition)>0)
    printf("delete %s->%s\n",TN_FG,condition);
//-----------------------
  sprintf(condition, "%s=%d",FN_GQ_RUNLOG_ID,runlog_id);
  if(delete_globalquantities(mcdb,condition)>0)
    printf("delete %s->%s\n",TN_GQ,condition);
//-----------------------
  sprintf(condition, "%s=%d",FN_MG_EQUILIBRIUM_ID, equilibrium_id);
  if(delete_magneticgeometry(mcdb,condition)>0)
    printf("delete %s->%s\n",TN_MG,condition);
//-----------------------
  sprintf(condition, "%s=%d",FN_EQ_ID,equilibrium_id);
  nr=load_equilibrium(mcdb,&equilibrium,condition);
  if(delete_equilibrium(mcdb,condition)>0)
    printf("delete %s->%s\n",TN_EQ,condition);
//-----------------------
  sprintf(condition, "%s=%d",FN_EI_ID,equilibrium->equilibriuminput_id);
  nr = nr1= load_equilibriuminput(mcdb,&equilibriuminput,condition);
  if(delete_equilibriuminput(mcdb,condition)>0)
    printf("delete %s->%s\n",TN_EI,condition);

  while(nr--) {
    PROFILE *profile=NULL;
    sprintf(condition, "%s=%d",FN_PR_ID,equilibriuminput->inputprofile1_id);
    int i = load_profile(mcdb,&profile,condition);
    while(i--) {
      sprintf(condition, "%s=%d",FN_GP_PROFILE_ID,profile[i].id);
      if(delete_gridprofile(mcdb,condition)>0)
        printf("delete %s->%s\n",TN_GP,condition);
    }
    free(profile); profile=NULL;
    sprintf(condition, "%s=%d",FN_PR_ID,equilibriuminput->inputprofile1_id);
    if(delete_profile(mcdb,condition)>0)
      printf("delete %s->%s\n",TN_PR,condition);
  }

  while(nr1--) {
    PROFILE *profile=NULL;
    sprintf(condition, "%s=%d",FN_PR_ID,equilibriuminput->inputprofile2_id);
    int i = load_profile(mcdb,&profile,condition);
    while(i--) {
      sprintf(condition, "%s=%d",FN_GP_PROFILE_ID,profile[i].id);
      if(delete_gridprofile(mcdb,condition)>0)
        printf("delete %s->%s\n",TN_GP,condition);
    }
    free(profile); profile=NULL;
    sprintf(condition, "%s=%d",FN_PR_ID,equilibriuminput->inputprofile1_id);
    if(delete_profile(mcdb,condition)>0)
      printf("delete %s->%s\n",TN_PR,condition);
  }

//=======================================================
  sprintf(condition, "%s=%d",FN_MV_RUNLOG_ID,runlog_id);
  nr=nr1=load_mgrid_vacuumfield(mcdb,&vacuumfield,condition);
  if(delete_mgrid_vacuumfield(mcdb,condition)>0)
    printf("delete %s->%s\n",TN_MV,condition);

  while(nr--) {
    sprintf(condition, "%s=%d",FN_CM_ID,vacuumfield[nr].currentmapping_id);
    if(delete_currentmapping(mcdb,condition)>0)
      printf("delete %s->%s\n",TN_CM,condition);
  }

  while(nr1--) {
    sprintf(condition, "%s=%d",FN_MX_ID,vacuumfield[nr1].makexgridparameter_id);
    if(delete_makexgridparameter(mcdb,condition)>0)
      printf("delete %s->%s\n",TN_MX,condition);
  }
//-----------------------
  sprintf(condition, "%s=%d",FN_ES_ID,equilibriuminput->experimentsetup_id);
  nr=load_experimentsetup(mcdb,&experimentsetup,condition);
  if(delete_experimentsetup(mcdb,condition)>0)
    printf("delete %s->%s\n",TN_ES,condition);

  while(nr--) {
    sprintf(condition, "%s=%d",FN_MS_ID,experimentsetup[nr].machinesetup_id);
    if(delete_machinesetup(mcdb,condition)>0)
      printf("delete %s->%s\n",TN_MS,condition);
  }
//------------------------
  sprintf(condition, "%s=%d",FN_RL_ID,runlog_id);
  nr = load_runlog(mcdb,&runlog,condition);
  if(delete_runlog(mcdb,condition)>0)
    printf("delete %s->%s\n",TN_RL,condition);

  sprintf(condition, "%s=%d",FN_CO_ID,runlog->code_id);
  if(delete_code(mcdb,condition)>0)
    printf("delete %s->%s\n",TN_CO,condition);

  free(runlog);
  free(equilibrium);
  free(equilibriuminput);
  free(experimentsetup);
  free(magneticgeometry);
  free(fluxsurfacegeometry);
  free(fluxsurface);
  free(vacuumfield);

  return true;
}

};   //namespace MConf

#endif //NO_MCDB

