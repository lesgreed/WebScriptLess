#ifndef MC_CEFIT_
#define MC_CEFIT_

#include "CArray2d.h"
#include "CBtrace.h"
#include "threadtypes.h"
#include <fstream>

namespace MConf {

class Bfield;

//****************************************************************************
/// Helper class that reads EFIT GEQDSK file and transforms it to
///  the Tokamak Symmetry flux coordinates.
/// In this coordinates the toroidal angle coinside with the cylindrical
/// angle, field lines are straight,
/// and Jacobian is proportional to R*R/(mu0*Ipol),
/// class is used by CStconfig::loadEfit()
class CEfit {
/// file formats, for internal use
  enum fileFormat  {
    UNDEFINED,     ///< unknown file format
    P_EQDSKold,    ///<  old EQDSK file format by A. Portone (before march 2005)
    P_EQDSK,       ///<  new EQDSK file format by A. Portone
    G_EQDSK,       ///<  old G EQDSK file format
    D_EQDSK,       ///<  EFITD in header 
    J_EQDSK,       ///<  jet  EQDSK file format
    C_TOKAMAK,     ///<  circular tokamak
    C_TOKAMAK_EQUIVALENT
  };

  /// Function object returns vector of magnetic field for given point in cylindrical coordinates.
  class Bfield {
    CEfit *mC;
  public:
    Bfield(CEfit &efit) { mC = &efit; }
    Vector3d operator() ( const Vector3d &cyl ) const { 
      return mC->getBcyl(cyl); }
  };
  CBtrace<Bfield> btrace;
  bool meshOK;           ///< result of loading: true if Okey
  double pi;             ///< == 3.14159265358979323846
  double twopi;
  double cr [4],cz [4];  ///< Lagrange coefficients for polynomial interpolation.
  double crd[4],czd [4]; ///< Derivatives of Lagrange coefficients.
  double Rmajor;
  double xdim,zdim,zmid;
  double Raxis, Redge, Zaxis, Bcenter;
  double Rmin, Rmax, Zmin, Zmax;   ///< mesh size in meters
  double dR;                  ///< step of mesh
  double dZ;                  ///< step of mesh
  double edR, edZ;            ///< Inverse mesh steps.
  double sMax;                ///< mesh covers surfaces from 0 to smax
  int j0;                     ///< the 1st point for the Lagrange interpolation in R-direction
  int k0;                     ///< the 1st point for the Lagrange interpolation in Z-direction
/// \name EFIT arrays and others
//@{
  int nR;                     ///< number of mesh points in R-direction
  int nZ;                     ///< number of mesh points in Z-direction
 
  ////int new_nR;

  double psiAxis;  ///< poloidal flux at magnetic axis in Weber/rad
  double psiEdge;  ///< poloidal flux at the plasma boundary in Weber/rad
  double deltaPsi; // =  psiEdge - psiAxis 
  double dpsi;     ///< = (psiEdge-psiAxis)/(nR-1);

  pBase::ngArray<double> psi_; ///< Uniform poloidal flux grid, of nR points from psiAxis to psiEdge
  pBase::ngArray<double> Fp_;  ///< FPOL[nR] - Poloidal current function in m*T, F = R*Btor on flux grid
  pBase::ngArray<double> FFp_; ///< FF'(psi) in (mT)2 / (Weber/rad) on uniform flux grid
  pBase::ngArray<double> Pp_;  ///< P'(psi) in (nt/m2) / (Weber/rad) on uniform flux grid
  pBase::ngArray<double> q_;   ///< q values on uniform flux grid from axis to boundary
//  pBase::ngArray<double> s_;   ///< norm. tor. flux  on uniform pol. flux grid, see FluxTor

  pBase::ngArray<double> R;   ///< mesh R [nR]
  pBase::ngArray<double> Z;   ///< mesh Z [nZ]

  CArray2d<double> psirz;     ///< Poloidal flux in Weber/rad on the rectancular mesh, psirz(nR,nZ)
  CArray2d<double> mNormFlux; ///< mesh array double mNormToroidalFlux(nR,nZ)

  pBase::ngArray<Vector3d> LCMS; ///< last closed magnetic surface
  //@}

/// \name arrays from B-tracing
//@{
  int Ns;       ///< number of surfaces
  int M;        ///< number of poloidal modes
  int Ntrace;
  pBase::ngArray<Vector3d> RBZm_;  ///< RBZm_[M+1] Fourier mods
  pBase::ngArray<Vector3d> RBZm1_; ///< RBZm1_[M+1]
  pBase::ngArray<Vector3d> RBZt_;  ///< RBZt_[Ntrace]
  pBase::ngArray<Vector3d> RBZt_smooth;  ///< RBZt_[Ntrace]
  pBase::ngArray<double> mR_;      ///< mR_[Ns]
  pBase::ngArray<double> mpsi_;    ///< mpsi_[Ns]
  pBase::ngArray<double> ms_;      ///< ms_[Ns]
  //pBase::ngArray<double> miota_;   ///< miota_[Ns] is the 1/q values recieved by tracing
  pBase::ngArray<double> mq_;      ///< mq_[Ns]  is the q values recieved by tracing 
  double fluxTorMax;
  //@}

  uiMutex efitMutex;
  
  int NSurfaces;
  int nLastToRemove;
  int Mpol; 

public:
  /// The method creates a clone of CEfit using allocation by new.
  /// The caller is responsible for deleting the clone
  CEfit * clone() {
    CEfit *ef = new CEfit;
    *ef = *this; 
    return ef;
  }
  /// The method creates a clone of CEfit using allocation by new.
  //std::auto_ptr<CEfit>  aclone() {
  //  CEfit *ef = new CEfit;
  //  *ef = *this;
  //  std::auto_ptr<CEfit> a(ef);
  //  return a;
  //}

  /// The constructor creates an empty instance of CEfit.
  CEfit();
  /// @param[in] M  is the number of poloidal mode for Fourier expantion of surfaces.
  /// @param[in] Ns is the number of flux surfaces to build
  CEfit(int M, int Ns, int nLastToSkip);
  /// The destructor is a virtual method that frees all memory used.
  virtual ~CEfit() {clear();}
  /// Free all memory used by CEfit.
  /// Usually no needs to use this method. Destructor does this.
  void clear();
  /// The method loads magnetic configuration stored in the EQDSK file \c fname.
  /// The method understands format of a file by analyzing its contents.
  /// \sa the detailed description of CStconfig::load().
  /// The method import equilibrium from the EQDSK File into straight-field-line coordinates 
  /// system (s,theta,phi) with the cylindrical angle phi as the toroidal angle. 
  /// The right handed cylindrical coodinate system (r,phi,z) is used. The system (s,theta,phi) can be left or right handed; 
  ///  the direction of theta and correcpondingly the Jacobian sign depends on Bpol, Btor directions and the sign of safety factor, see below.
  ///
  ///  \n The value of the following parameters BpolScale, BtorScale, BpolSign, BtorSign, qSign, psiOverTwopi 
  ///  can be put in the first line of the EQDSK file, for example 
  /// \code
  /// RC-ITER...             06-DEC-02                  3 129 129  EQDSK   BpolScale=0 BtorScale=0 BpolSign=1 BtorSign=1 qSign=1 psi/2pi=yes    
  /// \endcode
  /// @param[in] fname is the name of G EQDSK file
  /// @param[in] BpolScale is the scale factor for the poloidal field, the factor sign is ignored, 
  ///   \li      BpolScale=0 means that the scale factor is taken from the the first line of the EQDSK file 
  ///            if it found, otherwise no scaling is performed.
  /// @param[in] BtorScale is the scale factor for the toroidal field, the same logic as for poloidal field scaling.

  /// @param[in] BpolSign sets the sign of the poloidal field (by changing sign of the poloidal flux psi):
  ///    \li     +1 means positive Bpol, positive Bpol is created by the positive currrent density j_phi flowing in direction of the angle phi 
  ///             in the right handed cylindrical coodinate system (r,phi,z);
  ///    \li     -1 means negative Bpol created by the negative currrent density flowing in counter direction of the cylindrical angle phi; 
  ///    \li      0 means that Btor direction is taken from the BpolSign=n (where n is one from =-1,0,1) if it stored in the first line of the EQDSK file,
  ///    \li     if BpolSign not found in the file or BpolSign=0 then Bpol direction is defined by the psirz array stored in the EQDSK file.
  ///    \li       B = (-dpsirz/dz, Fpol, dpsirz/dr)/r 

  /// @param[in] BtorSign sets the sign of the toroidal field:
  ///    \li     +1 means positive Btor in direction of the angle phi in the right handed cylindrical coodinate system (r,phi,z),
  ///    \li     -1 means negative Btor in counter direction of the cylindrical angle phi, 
  ///    \li      0 means that Btor direction is taken from the BtorSign=n (where n is one from =-1,0,1) if it stored in the first line of the EQDSK file,
  ///    \li     if BtorSign not found in the file or BtorSign=0 then Btor direction is defined by the poloidal current function stored in the EQDSK file.
  ///    \li       B = (-dpsirz/dz, Fpol, dpsirz/dr)/r 

  /// @param[in] qSign sets the sign of the safety factor, qSign is taken from the EQDSK file if this parameter is zero, 
  ///    \li     if qSign not found in the file or qSign=0 in it then the safety factor is defined by the q values stored in the EQDSK file.

  /// @param[in] psiOverTwopi defines whether the poloidal flux psi is divided by 2pi:
  ///    \li      -1 means no, i.e psi is not devided over 2pi, this corresponds to COCOS>=11.
  ///    \li      0,1 assume that the flux psi stored in the file is already divided by 2pi.
  ///    \li      in case psiOverTwopi is 0 the psiOverTwopi is derived from the string psi/2pi=val stored in the EQDSK file hader, where val is yes or no.

  /// @param[in] fPostn is the position in the file fname from where to start reading EQDSK staff. 
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  virtual bool load(const char * fname, double BpolScale=0, double BtorScale=0,int BpolSign=0,int BtorSign=0, 
                  int qSign=0, int psiOverTwopi=0, long fPostn=0);
  /// Create circular tokamak using a,R0, B00, iota(0) from given configuration mconf
  virtual bool load(const CStconfig *mconf);

  /// The method returns \c true if no error occurs during last operation.
  bool isOK() const {return meshOK;}
// Methods to retrieve data from the mesh
  /// \name Methods to retrieve data from the EFIT-mesh
  //@{
  /// The method returns cylindrical B-vector.
  /// 4-points Lagrange interpolation is used.
  /// @param[in] cyl is the cylindrical coordinates
  /// @param[out] psi if not NULL then this address at return contains the poloidal flux in the point cyl.
  /// @return cylindrical B-vector or zero-vector if \c cyl lies outside mesh \c (Rmin,Rmax,Zmin,Zmax)
  Vector3d getBcyl (const Vector3d &cyl,double *psi=0) const;  // B-vector in cylindrical coord. in point cyl
  /// The method finds the flux label \c s (normalized toroidal flux)
  ///  which corresponds to the point \c cyl in space.
  /// @param cyl is the input parameter with the cylindrical coordinates of the point of interest.
  /// @return \c s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// \c s\>1 means that the point lies outside the LCMS, retuned parameters are not reliable here; \n
  double   cyl2s   (const Vector3d &cyl) const;     // map cylindrical coordinates to normalized flux
  //@}

  int getNs() const {return Ns;}
  int getM() const  {return M;}

  /// 0\<=is\<=getNs()
  double getS(int is) const    {return ms_[is];    }
  double getIota(int is) const {return 1./mq_[is]; }  // iota from B-field line tracing
  ////2015_1129double getIotaPrime(int is) const {
  ////  double ip;
  ////  interp2(miota_[is], miota_, ms_,Ns, &ip);// iota calculated by B-field line tracing
  ////  return ip;
  ////}

  double getPp(int is) const   {return Pp(mpsi_[is]);}
  double getIpol(int is) const { // in Ampers
    return RBtor(mpsi_[is])*0.5e7; //*twopi/mu0 = 0.5e7;
  }
  ////2015_1129double getJtor(int is) const {  // do i need this method?
  ////  return mR_[is]*interp2psi_(mpsi_[is], Pp_)+interp2psi_(mpsi_[is],FFp_)/mR_[is];
  ////}

  double getFlux() const {return fluxTorMax;}
  double getRmajor() const {return Rmajor;}
  Vector3d getMAxis() const {return Vector3d(Raxis,0,Zaxis);}

  /// \name Methods to create Fourier decomposition of flux surfaces
  //@{
  /// 0\<=is\<=getNs()
  void createSurface(int is, double &dVds, double &Itor, std::ofstream *out=NULL);
  /// 0\<=m\<=getM()
  double getRm0(int m) const {return RBZm_[m][0];}   // R_m  stellarator symmetric part cos(m*theta)
  double getBm0(int m) const {return RBZm_[m][1];}   // B_m  stellarator symmetric part cos(m*theta)
  double getZm0(int m) const {return RBZm_[m][2];}   // Z_m  stellarator symmetric part sin(m*theta)
  double getRm1(int m) const {return RBZm1_[m][0];}  // R_m  asymmetric part sin(m*theta)
  double getBm1(int m) const {return RBZm1_[m][1];}  // B_m  asymmetric part sin(m*theta)
  double getZm1(int m) const {return RBZm1_[m][2];}  // Z_m  asymmetric part cos(m*theta)
  //@}

private:
  int getPsirzSign();
  Vector3d findPsiEdge(double &qEdge,double &psi_edge ,double &r_edge, double &R1, double &R2);    // find LCMS

  double mapFunction(double x) {
    //return (tanh(4*(x-0.5))/tanh(4*(1.-0.5))+1)/2;
    return (1-pow(1-x,2.3));
  }

  bool loadCircularTokamak(const char * fname, CEfit::fileFormat fileFmt, const CStconfig *mconf, long fPostn);

  bool loadx(const char * fname, CEfit::fileFormat fileFmt, double scale_Bpol=0,double scale_Btor=0,int sign_Bpol=0,int sign_Btor=0, 
                  int sign_Q=0, int psiOverTwopi=0, const CStconfig *mconf=0, long fPostn=0);

  void calculateQ();
 // template <class efit> void calculateQExe(efit * ef,int low, int upper, bool setup);
  void calculateQExe(CEfit * ef,int low, int upper, bool setup);

  bool set(int Ns, int M, fileFormat fileFmt);
  CEfit::fileFormat getfileformat(const char * fname, long fPostn) const;
  /// Get memory
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool resize();
  /// The function finds the integer coordinates on the mesh and calculates
  /// interpolation coefficients.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool FindMeshCoordinates(const Vector3d &cyl);

  template <class T> T square(const T x) const   { return( x*x );  }
  template <class T> T mmin ( T x,  T y) const { return (x<y)?x:y; }
  template <class T> T mmax ( T x,  T y) const { return (x>y)?x:y; }

  double iota    (const double psi) const { return 1/interp2  (psi, mq_, mpsi_,Ns);  }  // iota calculated by B-field line tracing
  double qTraced (const double psi) const { return interp2    (psi, mq_, mpsi_,Ns);  }  // q using mq_ calculated by B-field line tracing


  double q    (const double psi) const { return interp2psi_(psi, q_ ); }
  double Pp   (const double psi) const { return interp2psi_(psi, Pp_)/qTraced(psi)*fluxTorMax; } ///< dP/ds in (nt/m2)
  double RBtor(const double psi) const { return interp2psi_(psi, Fp_); } ///< R*Btor - poloidal current function in m*T

  double sNew (const double psi) const { return interp2    (psi, ms_, mpsi_,Ns);  }

  double cyl2psi (const Vector3d &cyl) const;     // map cylindrical coordinates to poloidal flux

  int bsearch(double u, const pBase::ngArray<double> &x,int size) const;
  // interpolation on new mesh
  double interp2(double u, const pBase::ngArray<double> &y, const pBase::ngArray<double> &x,int size, double *yp=0) const;
  // interpolation on old mesh
  double interp2psi_(double fluxPol, const pBase::ngArray<double> &y) const;
  double FluxTor(double psi1, double psi2) const;
  double mFluxTor(double R1, double R2)  const;

  double traceStep(double qGuess);
  double traceForQ(const Vector3d &cyl, double qguess);
  double traceForSurface(int is, double &dVds, double &Itor);

private: // downhill simplex method, see NR
  /// \name Downhill simplex method
  //@{
  int psiSign;
  //double funk(double *x) { return 1+psiSign*cyl2psi(Vector3d(x[0],0,x[1])); }
/// Function object returns the the poloidal flux.
  class amoebaFunc;
  friend class amoebaFunc;

  class amoebaFunc {
    int sign_;
    int ncall;
    CEfit *ef;
  public:
    amoebaFunc(CEfit *Ef, int psiSign) { ef = Ef; sign_=psiSign;ncall=0;}
/// The operator returns the poloidal flux
    double operator() (const double *x) {
      ncall++;
      return 1+sign_*ef->cyl2psi(Vector3d(x[0],0,x[1]));
    }
    int ncalls() {return ncall;}
  };
  void   amoeba(double **p,double *y,int ndim,double ftol,amoebaFunc &func);
  double amotry(double **p,double *y,double *psum,int ndim,int ihi,double fac,double *ptry,amoebaFunc &func);
//@}

private:
  /// \name Obsolete
  //@{
  /// The method enables the test of a point position before data retrieving.
  /// The methods cyl2s(), getBcyl(), and others take as input
  /// parameter cylindrical/cartesian coordinates of the point and
  /// return values of <b>s, B,...</b> in this point. With the position test enabled
  /// the all function, which retrieve data from the mesh, at first do test
  /// whether point lies inside sMax-surface.
  /// Use this feature if you don't know where the point is.
  void positionTestOn()  {doPositionCheck=true;}
  /// The method disables the test of a point position.
  /// If you know that the point, in which you need s and \b B lies inside sMax-surface
  /// then use this method to disable test of a position with respect to the sMax.
  /// This speeds up calculation. \sa positionTestOn(), cyl2s().
  void positionTestOff() {doPositionCheck=false;}
  /// The method returns \b true if position test is enabled.
  /// \sa positionTestOn().
  bool positionTest() const { return doPositionCheck;}
  bool doPositionCheck;     ///< do test whether point lies inside sMax
  double mod2pi   (double x) const { return (x>=0)?fmod(x,twopi):(fmod(x,twopi)+twopi);} // remove period
  template <class T> int sign(T x) const { return x>=0?1:-1; }
  //@}
};

};   //namespace MConf

#endif
