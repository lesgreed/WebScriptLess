#ifndef MC_C3DMESH_
#define MC_C3DMESH_

#include "threadtypes.h"
#include "CRayTrace.h"
#include "CArray1d.h"
#include "CArray3d.h"

namespace MConf {

//****************************************************************************
//                          2004-2005,   Author: yuriy.turkin@ipp.mpg.de
// Magnetic configuration of a stellarator
// C3dMesh is derived from CStconfig and CRayTrace classes.
//  C3dMesh is used for
// Goal:
// 1. Tabulate magnetic field, flux label, grad(s) on 3d-mesh
//    in cylindrical coordinates.
// 2. Then use sequential cubic interpolation to find the magnetic field
//    in arbitrary space point.
//
/// This class tabulates magnetic field, flux surface label, grad(s) on a 3d-mesh
/// in cylindrical coordinates and provides functions for interpolation.
/// 4 point Lagrange interpolation in each direction is used.
/// \par Usage:
/// \code
/// #include "C3dMesh.h"
/// #include <iostream>
/// int main() {
///   C3dMesh magConf;
///   if(!magConf.load ("w7x-sc1.bin4")) exit(1);   // exit if not OK
///   if(!magConf.isMeshOK()) // create mesh if needed
///      magConf.createMesh(0.02,0.02,pi/180); // mesh=2cm X 2cm X 1degree
///   if(!magConf.isMeshOK()) exit(1);    // exit if not OK
///
///   Vector3d cyl(6,0,0); double s;
///   Vector3d B, dBdr,dBdfir,dBdz, grads, curlB;
///   s = magConf.M3Dcyl2s(cyl);
///   B = magConf.M3DgetdBcyl(cyl,dBdr,dBdfir,dBdz);
///   std::cout <<dBdr<<std::endl;
///   if(magConf.cylIsInside(cyl)) { ;/*do something*/ }
///   if(s<1.001) { ;/*do something*/ }
///
///   return 0;
/// }
/// \endcode
/// \sa class CStconfig, class CRayTrace.
class C3dMesh : public CRayTrace {
  bool meshOK;        ///< result of the last operation
  double BnormMesh;   ///< normalizing factor for B
  double B0new;       ///< current value of B0
  double B0mesh;      ///< value of B0 used when the mesh was created
  double cfi [4],cr [4],cz [4];  ///< Lagrange coefficients for polynomial interpolation.
  double cfid[4],crd[4],czd [4]; ///< Derivatives of Lagrange coefficients.
  double Rmin, Rmax, Zmin, Zmax, Fmin, Fmax;   ///< mesh size in meters
  double dR;                  ///< step of mesh
  double dZ;                  ///< step of mesh
  double dFi;                 ///< step of mesh
  double edR, edZ, edFi;      ///< Inverse mesh steps.
  double sMax;                ///< mesh covers surfaces from 0 to smax
  double truncationLevel;     ///< truncation level of Bmn used when the mesh was created
  double epsALevel;           ///< accuracy of coordinate transformation used when the mesh was created
  int i0;                     ///< the 1st point for the Lagrange interpolation in fi-direction
  int j0;                     ///< the 1st point for the Lagrange interpolation in R-direction
  int k0;                     ///< the 1st point for the Lagrange interpolation in Z-direction
  int nR;                     ///< number of mesh points in R-direction
  int nZ;                     ///< number of mesh points in Z-direction
  int nFi;                    ///< number of mesh points in fi-direction
  CArray1d<double> R;   ///< mesh R [nR]
  CArray1d<double> Z;   ///< mesh Z [nZ]
  CArray1d<double> Fi;  ///< mesh Fi[-1..nFi]

  CArray3d<Vector3d> mgrads;     ///< mesh array Vector3d mgrads (-1:nFi,nRc,nZc),2nd and 3rd dimentions depend on 1st
  CArray3d<Vector3d> mgradB;     ///< mesh array Vector3d mgradB (-1:nFi,nRc,nZc),2nd and 3rd dimentions depend on 1st
  CArray3d<Vector3d> mBfield;    ///< mesh array Vector3d mBfield(-1:nFi,nRc,nZc),2nd and 3rd dimentions depend on 1st
  CArray3d<double>   mNormFlux;  ///< mesh array double mNormFlux(-1:nFi,nRc,nZc),2nd and 3rd dimentions depend on 1st
  CArray3d<char>     mWhere;     ///< mesh array char mWhere(-1:nFi,nRc,nZc),2nd and 3rd dimentions depend on 1st
//CArray3d<double>   mJacobian;  ///< mesh array double mJacobian(-1:nFi,nRc,nZc),2nd and 3rd dimentions depend on 1st
  CArray2d<Vector3d> sMaxSurface;///< 2d-Array of Vector3d to hold vertixes of the sMax-surface
  double shift;

  bool   hasGradB;            ///< if true then the mesh contains grad(|B|)
  bool   fullMesh;            ///< if true then the mesh was created for whole machine, using periodicity condition
  bool   useSymmetry;         ///< if useSymmetry==true&&fullMesh==true then the mesh was created for whole stellarator using symmetry conditions
  bool   tokamakSymmetry;     ///< if true then no dependencies on toroidal angle
  bool   symmetryActivated;   ///< if true then symmetry is used
  bool   doPositionCheck;     ///< do test whether point lies inside sMax
  bool   doDivergenceCorrection;
  uiMutex meshMutex;

public:
  /// The method creates a clone of C3dMesh using operator new.
  /// The caller is responsible for deleting the clone.
  C3dMesh * clone() {
    C3dMesh *mc = new C3dMesh;
    *mc = *this; 
    return mc;
  }
  /// The constructor creates an empty instance of C3dMesh.
  C3dMesh();
  /// The destructor is a virtual method that frees all memory used.
  virtual ~C3dMesh() {
    freeMesh();};
  /// The method loads the 3d-mesh and/or magnetic configuration stored in the file \c fname.
  /// The method understands format of a file by analyzing its contents.
  /// \sa the detailed description of CStconfig::load()
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  virtual bool load(const char * fname) { return load(fname,CStconfig::UNDEFINED );}  
  virtual bool load(const char * fname, CStconfig::fileFormat fmt);
  /// The method returns \b true if no error occurs during last operation.
  bool isMeshOK() const {return meshOK;};
  /// This method frees internal mesh arrays.
  void clearMesh();
  /// This method frees all internal arrays (including base classes) and makes this object empty.
  virtual void clear();
  /// Restore initial value of the magnetic field
  virtual void restoreB0();              // restore initial value of the magnetic field on axis
  /// Set minimum magnetic field on magnetic axis equals B
  /// @param B is the value[Tesla] of magnetic field to set
  virtual void setB0(double B);              // Set minimum magnetic field on axis equals B
  /// Set the value of magnetic field on magnetic axis at given cylindrical angle.
  /// @param B is the value[Tesla] of magnetic field to set
  /// @param ficyl is the value[radians] of cylindrical angle
  virtual void setB0(double B,double ficyl); // Set magnetic field on axis equals B at the cyl angle ficyl
  /// Set B00-component of magnetic field on magnetic axis.
  /// @param B00 is the value[Tesla] of magnetic field to set
  virtual void setB00(double B00);
  /// Set scaling factor for value of the magnetic field
  /// @param multiplier is the scaling factor
  virtual void scaleB(double multiplier);
  /// Set the I<SUB>pol</SUB> value of the coil currents at LCMS.
  /// @param I is the value[Ampers] of the current to set
  virtual void setIcoil(double I);
  /// Get the scaling factor for the magnetic field value.
  virtual double getBscaleFactor();

#if 0
  /// Set average toroidal magnetic field B_phi(cylindrical component) on magnetic axis.
  /// @param Bphi0 is the value[Tesla] of magnetic field to set
  virtual void setAverageBphi0(double Bphi0);
#endif
  /// Set the last surface for the mesh creation.
  /// @param smax is the label value of the outermost surface used for the mesh creation,
  /// recommended value is 1.3, default value is 1.3
  /// This surface is needed for creating the extrapolated region, see also M3Dcyl2s(),getsMax().
  void setsMax(double smax) { sMax=smax<1.01?1.3:smax;};
  /// Get the value of normalized toroidal flux on the outermost surface used for the mesh creation.
  /// @return sMax the label value of the outermost surface used for the mesh creation.
  /// This surface is needed for creating the extrapolated region, see also M3Dcyl2s(),setsMax().
  double getsMax() const { return sMax;};
  /// The method creates 3d-mesh in cylindrical coordinates.
  /// The method tabulates the magnetic field, the flux surface label, and grad(s);
  /// the periodicity of the device is taken into account.
  /// @param dr is the mesh step in r-direction.
  /// @param dz is the mesh step in z-direction.
  /// @param dfi is the mesh step in fi-direction, where \e fi is the cylindrical angle in radians.
  /// @param numThreads is the number of execution threads used to speed up calculation
  ///        on a multiprocessor system.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// \note Recomended value of mesh parameters for W7-X are as follows:\n
  ///  dr = 0.02; dz = 0.02; dfi = pi/180;  (1degree) \n
  ///  please keep in mind that the function adjusts \b dfi parameter
  ///  in accordance with the periodicity.
  /// \attention The current truncation level of spectrum, value of magnetic field, and accuracy
  /// of coordinate transformatiom are used,
  /// see CStconfig::truncate(double level), CStconfig::setB0(), CStconfig::setAccuracy()
  /// \par Performance
  ///  In the following program the createMesh() needs 32sec to build the 3d-mesh
  ///  on Intel Core 2 Duo CPU P8400 at 2.26GHz using two execution threads; 
  ///  the code is compiled by Visual C++ 2005; the size of saved mesh file is about 13MB.
  /// \code
  ///  C3dMesh mc;
  ///  mc.load ("w7x-sc1.bc");
  ///  mc.truncate(1e-6);                // truncate spectrum
  ///  mc.setAccuracy(1e-3);             // set accuracy of coordinate transformation
  ///  mc.createMesh(0.02,0.02,pi/180,2);// create mesh  2cm X 2cm X 1degree
  ///  mc.writeMesh("w7x-sc1.mesh.bin4");// save mesh
  /// \endcode
  bool createMesh(double dr,double dz,double dfi,int numThreads=8);
  /// The method creates 3d-mesh in cylindrical coordinates.
  /// The method tabulates the magnetic field, the flux surface label, and grad(s);
  /// the \b periodicity (<b>B</b>(r,fi,z)=<b>B</b>(r,fi+period,z)) and
  ///  \b symmetry (B<sub>r</sub>(r,fi,z)=-B<sub>r</sub>(r,-fi,-z)) of the device are used.
  /// @param dr is the mesh step in r-direction.
  /// @param dz is the mesh step in z-direction.
  /// @param dfi is the mesh step in fi-direction, where \e fi is the cylindrical angle in radians.
  /// @param numThreads is the number of execution threads used to speed up calculation
  ///        on a multiprocessor system.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// \note Recomended value of mesh parameters for W7-X are as follows:\n
  ///  dr = 0.02; dz = 0.02; dfi = pi/180;  (1degree) \n
  ///  this method needs two times less memory then createMesh() and works two times faster.
  /// \par Performance
  ///  In the following program the createMeshUsingSymmetry() needs 18sec to build the 3d-mesh
  ///  on Intel Core 2 Duo CPU P8400 at 2.26GHz using two execution threads; 
  ///  the code is compiled by Visual C++ 2005 ; the size of saved mesh file is about 8MB.
  /// \code
  ///  C3dMesh mc;
  ///  mc.load ("w7x-sc1.bc");
  ///  mc.truncate(1e-6);                // truncate spectrum
  ///  mc.setAccuracy(1e-3);             // set accuracy of coordinate transformation
  ///  mc.createMeshUsingSymmetry(0.02,0.02,pi/180,2);  // create mesh  2cm X 2cm X 1degree
  ///  mc.writeMesh("w7x-sc1.mesh.bin4");// save mesh
  /// \endcode
  bool createMeshUsingSymmetry(double dr,double dz,double dfi,int numThreads=4);
  /// The method writes the mesh and current magnetic configuration into file.
  /// The format of the file is derived from the file extension.
  /// The following file extention are possible:
  ///  \li bin4  -- writeMesh("w.bin4"); //save into binary file 'w.bin4' using float numbers
  ///  \li bin8  -- writeMesh("w.bin8"); //save into binary file 'w.bin8' using double numbers
  /// @param fname -- filename
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// \attention The current truncation level of spectrum, value of magnetic field, and accuracy
  /// of coordinate transformatiom are used.
  bool writeMesh(const char * fname) const;

  ///
  bool writeMeshAscii(const char * fname) const;

// Methods to retrieve data from the mesh
  /// \name Methods to retrieve data from the mesh
  //@{
  /// The method tests whether the point cyl lies inside sMax
  /// @param cyl -- cylindrical coordinates of the point.
  /// @param distance if not NULL then this address at return
  ///        will contain the distance between cyl and the sMax,
  ///        the distance is negativ if the point lies inside sMax.
  /// @return  the return value is \b true if the point lies inside sMax; \b false otherwise.
  /// \sa setsMax()
  bool M3DcylIsInsideMesh(const Vector3d &cyl, double * distance=NULL) const;
  /// The method tests whether the point xyz lies inside sMax
  /// @param xyz -- cartesian coordinates of the point.
  /// @param distance if not NULL then this address at return
  ///        will contain the distance between xyz and the sMax,
  ///        the distance is negativ if the point lies inside sMax.
  /// @return  the return value is \b true if the point lies inside sMax; \b false otherwise.
  /// \sa setsMax()
  bool M3DxyzIsInsideMesh(const Vector3d &xyz, double * distance=NULL) const;
  /// The method enables the test of a point position before data retrieving.
  /// The methods M3Dcyl2s(), M3DgetBcyl(), and others take as input
  /// parameter cylindrical/cartesian coordinates of the point and
  /// return values of <b>s, B,...</b> in this point. With the position test enabled
  /// the all function, which retrieve data from the mesh, at first do test
  /// whether point lies inside sMax-surface.
  /// Use this method if you don't know where the point is.
  void M3DpositionTestOn()  {doPositionCheck=true;};
  /// The method disables the test of a point position.
  /// If you know that the point, in which you need s and \b B lies inside sMax-surface
  /// then use this method to disable test of a position with respect to the sMax.
  /// This speeds up calculation. \sa M3DpositionTestOn(), M3Dcyl2s().
  void M3DpositionTestOff() {doPositionCheck=false;};
  /// The method returns \b true if position test is enabled.
  /// \sa M3DpositionTestOn().
  bool M3DpositionTest() const { return doPositionCheck;};
  /// The method traces the Last Closed Magnetic Surface (LCMS).
  /// The method finds intersections of the ray with the LCMS, where the
  /// ray is R(t)=r0+rd*t,  and t>=0.
  /// @param r0 is the origin of the ray, must be outside with respect to the LCMS.
  /// @param rd is the direction of the ray
  /// @param entryPoint this vector at return contains the coordinates of the ray entry point into the plasma
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool M3DgetRayEntryPoint(const Vector3d &r0,const Vector3d &rd,Vector3d &entryPoint) const;
  /// The method finds the flux label \b s which corresponds to the point in space.
  /// @param cyl is the input parameter with the cylindrical coordinates of the point of interest.
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  /// \b s>1000 means that the point lies outside the \b sMax. \n
  ///  The extrapolated region is needed in order to do interpolation near LCMS, see also setsMax(),getsMax().
  double   M3Dcyl2s   (const Vector3d &cyl) const;     // map cylindrical coordinates to normalized flux
  /// The method finds the flux label \b s which corresponds to the point in space.
  /// @param xyz is the input parameter with the cartesian coordinates of the point of interest.
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  /// \b s>1000 means that the point lies outside the \b sMax. \n
  ///  The extrapolated region is needed in order to do interpolation near LCMS, see also setsMax().
  double   M3Dxyz2s   (const Vector3d &xyz) const;     // map cartesian  coordinates  to normalized flux
  /// The method returns cylindrical B-vector.
  /// 4-points Lagrange interpolation is used.
  /// @param cyl is the cylindrical coordinates
  /// @return cylindrical B-vector or zero-vector if \b cyl lies outside the sMax-surface, see M3Dcyl2s()
  Vector3d M3DgetBcyl (const Vector3d &cyl) const;     // B-vector in cylindrical coord. in point cyl
  /// The method returns cylindrical B-vector.
  /// 4-points Lagrange interpolation is used.
  /// @param cyl is the cylindrical coordinates
  /// @param s is the reference to the output parameter, at return it contains the flux surface label
  /// @return cylindrical B-vector or zero-vector if \b cyl lies outside the sMax-surface, see M3Dcyl2s()
  Vector3d M3DgetBcyl (const Vector3d &cyl,double &s) const;     // B-vector in cylindrical coord. in point cyl
  /// The method returns cartesian B-vector.
  /// 4-points Lagrange interpolation is used.
  /// @param xyz is the cartesian coordinates
  /// @return cartesian B-vector or zero-vector if \b xyz lies outside the sMax-surface, see M3Dcyl2s()
  Vector3d M3DgetBxyz (const Vector3d &xyz) const;     // B-vector in cartesian coord.  in point xyz(cartesian)
  /// The method returns cartesian B-vector.
  /// 4-points Lagrange interpolation is used.
  /// @param xyz is the cartesian coordinates
  /// @param s is the reference to the output parameter, at return it contains the flux surface label
  /// @return cartesian B-vector or zero-vector if \b xyz lies outside the sMax-surface, see M3Dcyl2s()
  Vector3d M3DgetBxyz (const Vector3d &xyz,double &s) const;
  /// The method returns cylindrical vectors B and partial derivatives of B.
  /// 4-points Lagrange interpolation is used.
  /// @param cyl is the input parameter with the cylindrical coordinates of the point of interest
  /// @param dBdr,dBdfir,dBdz are the references to the output parameters,
  ///        where \b B is a vector in cylindrical coordinates and
  ///        dBdr=d<b>B</b>/dr, dBdfir=d<b>B</b>/dfi/r, dBdz=d<b>B</b>/dz
  /// @return cylindrical B-vector or zero-vector if \b cyl lies outside the sMax-surface, see M3Dcyl2s()
  Vector3d M3DgetdBcyl(const Vector3d &cyl,Vector3d &dBdr,Vector3d &dBdfir,Vector3d &dBdz) const;
  /// The method returns cartesian vectors B and partial derivatives of B.
  /// 4-points Lagrange interpolation is used.
  /// @param xyz is the input parameter with the cartesian coordinates of the point of interest
  /// @param dBdx,dBdy,dBdz are the references to the output parameters,
  ///        where \b B is a vector in cartesian coordinates and
  ///        dBdx=d<b>B</b>/dx, dBy=d<b>B</b>/dy, dBdz=d<b>B</b>/dz
  /// @return cartesian B-vector or zero-vector if \b xyz lies outside the sMax-surface, see M3Dcyl2s()
  Vector3d M3DgetdBxyz(const Vector3d &xyz,Vector3d &dBdx,Vector3d &dBdy,Vector3d &dBdz) const;
  /// The method returns cartesian vectors B, partial derivatives of B, and grad(s).
  /// 4-points Lagrange interpolation is used.
  /// @param xyz is the input parameter with the cartesian coordinates of the point of interest
  /// @param B,dBdx,dBdy,dBdz are the references to the output parameters,
  ///        where \b B is a vector in cartesian coordinates and
  ///        dBdx=d<b>B</b>/dx, dBy=d<b>B</b>/dy, dBdz=d<b>B</b>/dz.
  /// @param grads is the references to the storage for grad(s).
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  /// \b s>1000 means that the point lies outside the \b sMax. \n
  ///  The extrapolated region is needed in order to do interpolation near LCMS, see also setsMax().
  double M3DgetdBGradsxyz(const Vector3d &xyz,Vector3d &B,Vector3d &dBdx,Vector3d &dBdy, Vector3d &dBdz, Vector3d &grads) const;
  /// The method returns cylindrical vectors B, partial derivatives of B, and grad(s).
  /// 4-points Lagrange interpolation is used.
  /// @param cyl is the input parameter with the cylindrical coordinates of the point of interest
  /// @param B,dBdr,dBdfir,dBdz are the references to the output parameters,
  ///        where \b B is a vector in cylindrical coordinates and
  ///        dBdr=d<b>B</b>/dr, dBdfir=d<b>B</b>/dfi/r, dBdz=d<b>B</b>/dz
  /// @param grads is the references to the storage for grad(s).
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  /// \b s>1000 means that the point lies outside the \b sMax-surface. \n
  ///  The extrapolated region is needed in order to do interpolation near LCMS, see also setsMax().
  double M3DgetdBGradscyl(const Vector3d &cyl,Vector3d &B,Vector3d &dBdr,Vector3d &dBdfir,Vector3d &dBdz,Vector3d &grads) const;
  /// The method returns cylindrical vector grad(s).
  /// @param cyl is the input parameter with the cylindrical coordinates of the point of interest.
  /// @return grad(s) or zero-vector if \b cyl lies outside the sMax-surface, see M3Dcyl2s()
  Vector3d M3DgetGrads(const Vector3d &cyl) const;
  /// The method returns cartesian vector grad(s).
  /// @param xyz is the input parameter with the cartesian coordinates of the point of interest.
  /// @return cartesian coordinates of grad(s) or zero-vector if \b xyz lies outside the sMax-surface, see M3Dcyl2s()
  Vector3d M3DgetGradsxyz(const Vector3d &xyz) const;
  /// The method returns cylindrical vector grad(B).
  /// @param cyl is the input parameter with the cylindrical coordinates of the point of interest.
  /// @return grad(|B|) or zero-vector if \b cyl lies outside the sMax-surface, see M3Dcyl2s()
  Vector3d M3DgetGradB(const Vector3d &cyl) const; 
  /// The method returns cartesian vector grad(B).
  /// @param xyz is the input parameter with the cartesian coordinates of the point of interest.
  /// @return grad(|B|) or zero-vector if \b xyz lies outside the sMax-surface, see M3Dcyl2s()
  Vector3d M3DgetGradBxyz(const Vector3d &xyz) const; 
  /// The method returns cylindrical vectors B, grad(|B|), and grad(s).
  /// @param cyl is the input parameter with the cylindrical coordinates of the point of interest
  /// @param B is the references to the B.
  /// @param gradB is the references to the storage for grad(|B|).
  /// @param grads is the references to the storage for grad(s).
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  /// \b s>1000 means that the point lies outside the \b sMax-surface. \n
  ///  The extrapolated region is needed in order to do interpolation near LCMS, see also setsMax().
  double M3DgetBandGradB(const Vector3d &cyl,Vector3d &B, Vector3d &gradB, Vector3d &grads) const;
  /// The method returns Cartesian vectors B, grad(|B|), and grad(s).
  /// @param xyz is the input parameter with the Cartesian coordinates of the point of interest
  /// @param B is the references to the B.
  /// @param gradB is the references to the storage for grad(|B|).
  /// @param grads is the references to the storage for grad(s).
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  /// \b s>1000 means that the point lies outside the \b sMax-surface. \n
  ///  The extrapolated region is needed in order to do interpolation near LCMS, see also setsMax().
  double M3DgetBandGradBxyz(const Vector3d &xyz,Vector3d &B, Vector3d &gradB, Vector3d &grads) const;
  Vector3d M3Dgetdgradscyl(const Vector3d &cyl,Vector3d &dgradsdr,Vector3d &dgradsdfir,Vector3d &dgradsdz) const;
  Vector3d M3Dgetdgradsxyz(const Vector3d &xyz,Vector3d &dgradsdx,Vector3d &dgradsdy,Vector3d &dgradsdz) const;
  
  
  //@}

private:
  /// Free all memory used by C3dMesh; only mesh related things are released.
  void freeMesh();
  /// Get memory
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool resize(bool array3d=true);
  /// create mesh in RZ-plane for cylindrical angle Fi[iCyl]
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool createCut(int iCyl, CStconfig &mc);
  /// create RZ-mesh-planes for cyl. angles from Fi[i1] to Fi[i2]
  /// @return thread ID
  uiThread createMeshPlanes(int i1, int i2, bool thread, CStconfig &mc);
  /// copy mesh planes using periodicity
  void copyCuts();
  static uiThreadFunc solveInThread (void *A);
  /// The function finds the integer coordinates on the mesh and calculates interpolation coefficients.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool FindMeshCoordinates(const Vector3d &cyl, bool noDerivatives=false);
  /// Ascertain the file format
  /// @return  fileFormat, see CStconfig::fileFormat.
  fileFormat getfileformat(const char * fullname,const char *nam,const char *ext) const;
  /// helper function for C3dMesh::load(), see also CStconfig::loadbin();
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  template <class T> bool loadMeshbin(const char * fname, T szType, fileFormat fmt);
  /// helper function for C3dMesh::writeMesh(), see also CStconfig::writebin();
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  template <class T> bool writeMeshbin(const char * fname, T szType) const;
  ///
  bool hasNeighbourInside(int iCyl,int j,int k) const;
  /// helper function for createMeshUsingSymmetry(), createMesh()
  bool createMeshEx(double dr,double dz,double df,double fi1,double fi2,int numThreads);

  bool CreateSmaxSurface();


private:  // disabled methods
  // not used; obsolete 
  /// \name Methods not used; obsolete
  //@{
  Vector3d M3DgetGrads1(const Vector3d &cyl) const;   // grad(s)
  /// The method tests whether the point lies inside LCMS
  /// @param cyl -- cylindrical coordinates of the point.
  /// @return  the return value is \b true if the point lies inside LCMS; \b false otherwise.
  bool M3DcylIsInside(const Vector3d &cyl); // check whether point lies inside LCMS, cyl (cylindrical coord.)
  /// The method tests whether the point lies inside LCMS
  /// @param xyz -- cartesian coordinates of the point.
  /// @return  the return value is \b true if the point lies inside LCMS; \b false otherwise.
  bool M3DxyzIsInside(const Vector3d &xyz); // check whether point lies inside LCMS, xyz (cartesian coord.)
  /// The method returns \b true if the mesh was created for whole machine; \b false otherwise.
  bool isFullMesh() const {return meshOK?fullMesh:false;};
  /// The method returns minimum value of mesh cylindrical angle.
  double getFimin() const {return fullMesh?0:Fmin;};
  /// The method returns maximum value of mesh cylindrical angle.
  double getFimax() const {return fullMesh?twopi:Fmax;};
  /// The method enables correction of the magnetic field derivatives.
  /// The following derivatives are expressed through others using
  ///  divB = 0:
  /// dBf/df/R = -(dBr/dr +dBz/dz) -Br/R
#if 0
  /// dBf/dr =  dBr/df/R-Bf/R
  /// dBf/dz =  dBz/df/R
  /// dBz/dr =  dBr/dz
  ///  curlB = 0
  /// \note  curlB = 0 can be used only in case of stellarator with Itor << Ipol,
  /// the corrections are disabled in case of tokamak.
#endif
  /// The method disabled. don't use it!
  void M3DBdivergenceCorrectionOn()  {doDivergenceCorrection=true;};
  /// The method disables correction of the magnetic field derivatives.
  /// \sa M3DBdivergenceCorectionOn().
  void M3DBdivergenceCorrectionOff() {doDivergenceCorrection=false;};
  void M3DdivergenceFree(const Vector3d &cyl,Vector3d &B,Vector3d &dBdr,Vector3d &dBdfir,Vector3d &dBdz) const;
  //@}

/** \example 3dmeshtest.cpp
 * This is an example of how to use the C3dMesh class.
 * It can be used to study error distribution after mesh creation, especially near 1<s<sMax, see setsMax()
 */

};

};   //namespace MConf

#endif
