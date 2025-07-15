#ifndef MC_CRAYTRACE_
#define MC_CRAYTRACE_

#include "CStconfig.h"
#include "CArray1d.h"
#include "CArray2d.h"

namespace MConf {

//****************************************************************************
//                          2004-2005,   Author: yuriy turkin at ipp.mpg.de
// Magnetic configuration of a stellarator
/*! \class CRayTrace
   \brief Find intersections of straight ray with the flux surfaces, this class is obsolete.

   CRayTrace is derived from CStconfig and is used for
   straight line tracing through flux surfaces.
   Usually you need only one method from this class,namely
   CRayTrace::traceLCMS(), to find the first intersection with LCMS. 
   Then you can use CStconfig::xyz2mag() to get magnetic coordinates along the ray,
   see \ref raytrace2.cpp.
   
   \note This class is obsolete; \n 
   use the method CStconfig::getRayEntryPoint() and CStconfig::xyz2mag(); \n
   or C3dMesh::M3DgetRayEntryPoint() and C3dMesh::M3Dxyz2s()
   \sa class CStconfig, class C3dMesh, and \ref raytrace2.cpp
 */
class CRayTrace : public CStconfig {
protected:
  Vector3d rayR0;           ///< origin of ray=rayR0+rayRd*len
  Vector3d rayRd;           ///< direction of ray=rayR0+rayRd*len
  CArray2d<Vector3d> svertex;  ///< 2d-array to hold vertixes of the surface
  CArray1d<Vector3d> Rxyz;  ///< array of intersection points of the ray with the surfaces
  CArray1d<double> len;     ///< array of distances from Ray Origin to intersection points
  CArray1d<double> reffc;   ///< radius reff, which corresponds intersection points
  CArray1d<double> vol;     ///< volume inside a surface
  CArray1d<int> snumber;    ///< array of surface numbers
  int maxTrace;                ///< max number of allocated elements in above arrays
  int sizeTrace;               ///< number of actual intersections
  int mNsurf;                  ///< number of surfaces to trace
  /// Free all memory used by CRayTrace.
  void   freeTrace();

public:
  /// The constructor creates an empty instance of CRayTrace.
  CRayTrace();
  /// The destructor is a virtual method that frees all memory used.
  virtual ~CRayTrace() {
    freeTrace(); }
  /// This method frees trace arrays.
  void clearTrace();
  /// This method frees all internal arrays (including base class) and makes this object empty.
  virtual void clear();

  /// The method traces the Last Closed Magnetic Surface (LCMS).
  /// The method finds intersections of the ray with the LCMS, where the
  /// ray is R(t)=r1+(r2-r1)/|r2-r1|*t,  and t>=0.
  /// @param r1 is the cartesian coordinates of the first point of the ray,
  ///        must be outside with respect to the LCMS.
  /// @param r2 is the cartesian coordinates of the second point of the ray
  /// @param entrypoint this vector at return contains the coordinates of the
  ///         entry point of the ray into plasma
  /// @param exitpoint this vector at return contains the coordinates of the exit point
  /// @param dpol is the angle step in degree along poloidal direction
  /// @param dtor is the angle step in degree along toroidal direction
  /// @return
  /// \li -2 -- the origin r1 of the ray lies inside LCMS;
  /// \li -1 -- the ray does not enters into plasma;
  /// \li  0 -- segment r1--r2 is outside plasma;
  /// \li  1 -- only entry point lies on the segment r1--r2, i.e. r2 is inside plasma;
  /// \li  2 -- entry and exit points lie on the segment r1--r2, i.e. r1 and r2 are outside plasma;
  ///
  ///  if returned value >= 0 then entry and exit points of the ray are
  ///   saved in \e entrypoint and \e exitpoint
  int traceLCMS(const Vector3d &r1,const Vector3d &r2, Vector3d &entrypoint,
         Vector3d &exitpoint, double dpol=1.,double dtor=1.);
  /// The method traces the Last Closed Magnetic Surface (LCMS); obsolete method.
  /// The method finds intersections of the ray with the LCMS.
  /// The ray is R(t)=r0+rd*t/|rd| with t>=0.
  /// @param r0 is the origin of the ray, must be outside with respect to the LCMS.
  /// @param rd is the direction of the ray
  /// @param dtheta is the poloidal step for the triangulization of the LCMS
  /// @param dphi is the toroidal step for the triangulization of the flux suface
  ///  \note  Obsolete function; use traceLCMS() in stead of this method.
  void   traceLast(const Vector3d &r0, const Vector3d &rd, double dtheta=0., double dphi=0.);
  /// The method traces magnetic surfaces; obsolete method.
  /// The method finds intersections of the ray with the flux surfaces.
  /// The ray is R(t)=r0+rd*t/|rd| with t>=0.
  /// @param r0 is the origin of the ray, must be outside with respect to the LCMS.
  /// @param rd is the direction of the ray
  /// @param dtheta is the poloidal step for the triangulization of the flux sufaces
  /// @param dphi is the toroidal step for the triangulization of the flux sufaces
  /// @param nSurf is the number of sufaces to trace, the surfaces are distributed
  ///         equidistantly on r<sub>eff</sub>.
  ///
  /// Use getTraceSize(),getTraceRxyz(),getTraceLen() methods to get trace information
  /// \note Obsolete function; it is better(and faster) to use
  ///  traceLCMS() and then xyz2mag(), see \ref raytrace1.cpp
  void   traceAll(const Vector3d &r0,const Vector3d &rd,double dtheta=0.,double dphi=0.,int nSurf=50);
  /// The method returns numbers of intersections of the ray with the magnetic surface or surfaces.
  /// This method is used after traceLast() or traceAll() to get numbers of intersections.
  /// @return numbers of intersections or 0 if no intersections were found.
  int    getTraceSize() const    {return sizeTrace;}
  /// The method returns pointer to the array of intersections.
  /// The array contains cartesian coordinates of intersections of the ray with the magnetic surface/s.
  /// This method is used after traceAll() or traceLast().
  /// The array size is getTraceSize()
  const Vector3d *getTraceRxyz() const {return Rxyz.constArray(); }
  /// The method returns pointer to the array of lengths.
  /// The array contains lengths between the first and successive intersections of the ray
  /// with the magnetic surface/s.
  /// This method is used after traceAll() or traceLast().
  /// The array size is getTraceSize()
  const double *getTraceLen() const    {return len.constArray();  }
  /// The method returns pointer to the array of effective radii.
  /// The array contains effective radii of surfaces, which the ray intersects.
  /// This method is used after traceAll() or traceLast().
  /// The array size is getTraceSize()
  const double *getTracereff() const   {return reffc.constArray(); }
  /// The method returns pointer to the array of volumes.
  /// The array contains volumes of surfaces, which the ray intersects.
  /// This method is used after traceAll() or traceLast().
  /// The array size is getTraceSize()
  const double *getTraceVolume() const {return vol.constArray();  }
  /// The method returns pointer to the array of surface numbers.
  /// The array contains surface numbers, which the ray intersects.
  /// This method is used after traceAll() or traceLast().
  /// The array size is getTraceSize()
  const int    *getTraceSurfNum() const {return snumber.constArray();}

private:
  /// Get memory
  bool   resizeTraceList();
  /// Save origin and direction of the ray, the ray is R(t)=r0+rd*t/|rd| with t>=0.
  void   setRay(const Vector3d &r0,const Vector3d &rd) {rayR0=r0; rayRd=rd/rd.abs();}
  /// The method creates flux surface \b s.
  /// nTheta, nPhi are the number of points in poloidal and toroidal directions
  bool   buildSurface(double s, int nTheta, int nPhi, double Phi1, double Phi2);
  /// Find intersections with the surface
  void   traceSurface (int is, double s, const CArray2d<Vector3d> &vertex, int step=1);
  /// Helper function for tracing functions
  bool   findSegment(double dtheta,double dphi,int &nTheta,int &nPhi,double &phi0,double &phi1);
  /// Find intersection of the ray with the triangle, the ray is R(t)=rayR0+rayRd*t  with t>=0.
  /// @param  r are the vertices of the triangle
  /// @param  ri this vector at return contains intersection coordinates (cartesian)
  /// @return distance (>=0) from the origin of the ray to the intersection point or
  ///                  -1 if no intersection was found.
  double isIntersection(const Vector3d *r, Vector3d &ri) const;
  /// Add intersection point in the list
  void   addInList(int is, double s, double ln, Vector3d &ri);
  /// Sort the list of intersections into increasing order of length along the ray.
  void   sort();

/** \example raytrace1.cpp
   This is an example of how to use the CRayTrace::traceLCMS().
   \note This example is obsolete; \n 
   use the method MConf::CStconfig::getRayEntryPoint() and MConf::CStconfig::xyz2mag(); \n
   or MConf::C3dMesh::M3DgetRayEntryPoint() and MConf::C3dMesh::M3Dxyz2s()
   \sa \ref raytrace2.cpp "raytrace2.cpp"
 */
/** \example raytrace2.cpp
   This is an example of how to get plasma and magnetic configuration parameters
   along line-of-sight (ray tracing through plasma)
 */
};

};   //namespace MConf

#endif // C_CRAYTRACE_
