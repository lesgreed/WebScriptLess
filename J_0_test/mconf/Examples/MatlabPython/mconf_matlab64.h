#ifndef __CPP2MATLAB_
#define __CPP2MATLAB_

#if defined( _WIN32 )
  #define WINDLLEXPORT  __declspec( dllexport )
#else  //#if (defined( __unix__ ) || (defined(linux)) || defined(__sun) )
  #define WINDLLEXPORT 
#endif

#ifdef __GNUC__
 #define __int64 int64_t
  #if defined( _LP64 )
    #define _M_X64  1 
  #endif
#endif

#if BUILD_LIB  
  /// A platform-specific type that is used to represent a pointer.
  #if defined( _M_X64 )
  //typedef __int64 MC_HANDLE;  // for 64-bit architecture
  typedef unsigned long long int   MC_HANDLE;  // don't remove spaces between int and MC_HANDLE
  #else
  typedef unsigned int  MC_HANDLE;  // for 32-bit architecture
  #endif
#else
typedef unsigned long long int MC_HANDLE;  // for 64-bit architecture
#endif

#ifdef __cplusplus
extern "C" {
#endif

  /// The function loads the magnetic configuration stored in the file \b fname.
  /// The function recognizes the format of the file by analyzing its contents;
  /// bc-, bc-binary-, LHD- and VMEC(wout, version>=6.20)-format are supported.
  /// For bc-format see \ref W7X_Format "W7-X format".
  /// On exit the function returns the handle of the object, that holds magnetic configuration.
  /// @param[in] fname is the name of the file to load.
  /// @return  if the function succeeds, the return value is \b non-zero, \b zero otherwise.
  /// \sa the detailed description of CStconfig::load()
  WINDLLEXPORT MC_HANDLE MCload      (char * fname);
  /// The function loads the \ref VmecCoordinates "VMEC wout-file" stored in the file \b fname.
  /// No transformation to Boozer representation is performed.
  /// On exit the function returns the handle of the object, that holds magnetic configuration.
  /// @param[in] fname is the name of the file to load.
  /// @return  if the function succeeds, the return value is \b non-zero, \b zero otherwise.
  WINDLLEXPORT MC_HANDLE MCloadVMEC  (char * fname);
  /// The function loads the magnetic configuration with Boozer coordinates from the file in LHD-format. 
  /// On exit the function returns the handle of the object, that holds magnetic configuration.
  /// @param[in] fname is the name of the file to load.
  /// @return  if the function succeeds, the return value is \b non-zero, \b zero otherwise.
  WINDLLEXPORT MC_HANDLE MCloadLHD   (char * fname);
  /// The function loads the \ref GEQDSK "G EQDSK file" 
  /// and transforms it to the \ref TokamakSymmetryFluxCoordinates "tokamak-symmetry flux coordinates".
  /// On exit the function returns the handle of the object, that holds magnetic configuration.
  /// @param[in] fname is the name of the file to load.
  /// @return  if the function succeeds, the return value is \b non-zero, \b zero otherwise.
  WINDLLEXPORT MC_HANDLE MCloadEFIT  (char * fname);
  /// The function creates the copy of the object.
  /// The function works as a C++ copy constructor.
  /// On exit the function returns the handle of the object, that holds magnetic configuration.
  /// @param[in] mConfSrc is the handle of the magnetic configuration to be copied.
  /// @return  if the function succeeds, the return value is \b non-zero, \b zero otherwise.
  /// \sa \ref ObjectCanBeCopied "Objects can be copied"
  /// \note use MCfree to free memory
  WINDLLEXPORT MC_HANDLE MCcopy  (MC_HANDLE mConfSrc);
  /// Free all memory used.
  /// The function works as a C++ destructor.
  /// @param[in] mConf is the handle of the magnetic configuration to be released.
  WINDLLEXPORT void   MCfree  (MC_HANDLE mConf);
  /// Get the B<SUB>00</SUB> value of the magnetic field on axis.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @return the B<SUB>00</SUB> value of the magnetic field on magnetic axis.
  WINDLLEXPORT double MCgetB00(MC_HANDLE mConf);
  /// Get the value of the magnetic field on axis at cylindrical angle.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] ficyl is the value[radians] of cylindrical angle.
  /// @return the value of the magnetic field on magnetic axis at cylindrical angle \e ficyl.
  WINDLLEXPORT double MCgetB0 (MC_HANDLE mConf,double ficyl);
  /// Set the B<SUB>00</SUB> value of magnetic field on magnetic axis.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] B00 is the value[Tesla] of magnetic field to set.
  WINDLLEXPORT void   MCsetB00(MC_HANDLE mConf,double B00);
  /// Set the value of smax.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] smax is the maximum value of s for extrapolation.
  WINDLLEXPORT void   MCsetsmax(MC_HANDLE mConf,double smax);
  /// Set the value of magnetic field on magnetic axis at given cylindrical angle.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] B0 is the value[Tesla] of magnetic field to set.
  /// @param[in] ficyl is the value[radians] of cylindrical angle.
  WINDLLEXPORT void   MCsetB0 (MC_HANDLE mConf,double B0,double ficyl);
  /// Set the value of magnetic field by setting the poloidal current Ip at LCMS (s=1).
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] Ip is the value[A] of the poloidal current to set.
  WINDLLEXPORT void   MCsetIpLCMS (MC_HANDLE mConf,double Ip);
  /// Set accuracy of coordinate transformation.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] epsA is the absolute accuracy[m] of the transformation from cylindrycal to magnetic coordinates.
  /// \note the default accuracy \e epsA is 10<sup>-5</sup>m if this function is not used.
  WINDLLEXPORT void   MCsetAccuracy(MC_HANDLE mConf,double epsA);
  /// Set truncation level to reduce number of harmonics used in summation.
  /// For each flux surface this function calculates M and N
  /// such that B<SUB>mn</SUB>/B<SUB>00</SUB>>epsTrunc for m<M and |n|<N.
  /// Then only harmonics with m<M and |n|<N are used in summation.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] epsTrunc is the truncation level.
  WINDLLEXPORT void   MCtruncate   (MC_HANDLE mConf,double epsTrunc);
  /// The function returns effective radius in [m].
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] s is the normalized toroidal flux.
  WINDLLEXPORT double MCreff (MC_HANDLE mConf,double s);
  /// The function returns V'(s) = dVolume/ds in [m^3].
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] s is the normalized toroidal flux.
  /// \note V' is negative if (s,theta,phi)-coordinate system is left handed
  WINDLLEXPORT double MCVprime(MC_HANDLE mConf,double s);
  /// The function returns plasma volume V(s) in [m^3].
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] s is the normalized toroidal flux.
  /// \note Volume is negative if (s,theta,phi)-coordinate system is left handed
  WINDLLEXPORT double MCVolume(MC_HANDLE mConf,double s);
  /// The function returns the cartesian coordinates of the point.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] magCoord is the array(s,theta,phi) in magnetic coordinates
  /// @param[out] xyz is the cartesian coordinates.
  WINDLLEXPORT void   MCmag2xyz(MC_HANDLE mConf,const double *magCoord, double *xyz);
  /// The method returns cartesian coordinates which correspond to the mixed coordinates.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] mixCoord are the mixed coordinates, array[3] =(s,theta,fi), where
  ///  s is the flux surface label,
  ///  theta is the magnetic poloidal angle,
  ///  fi is the cylindrical angle.
  /// @param[out] xyz is the cartesian coordinates.
  WINDLLEXPORT void   MCmix2xyz(MC_HANDLE mConf,const double *mixCoord, double *xyz);
  /// The function finds the flux label \b s which corresponds to the point in space.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] xyz is the input parameter with the cartesian coordinates of the point of interest.
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  WINDLLEXPORT double MCxyz2s(MC_HANDLE mConf,const double *xyz);
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] xyz is the input parameter with the cartesian coordinates of the point of interest.
  WINDLLEXPORT void MCxyz2mag(MC_HANDLE mConf,const double *xyz,double *mag);
  /// The function returns cartesian B-vector.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] xyz is the cartesian coordinates
  /// @param[out] B is the cartesian B-field vector
  /// @return \b s (normalized toroidal flux), which normally must satisfy \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here\n
  WINDLLEXPORT double MCgetBxyz(MC_HANDLE mConf,const double *xyz,double *B);
  /// The function returns cartesian grad(|B|) vector.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] xyz is the cartesian coordinates
  /// @param[out] gradB is the cartesian grad(|B|) vector
  /// @return \b s (normalized toroidal flux), which normally must satisfy \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here\n
  WINDLLEXPORT double MCgetGradBxyz(MC_HANDLE mConf,const double *xyz,double *gradB);
  /// The method traces the Last Closed Magnetic Surface (LCMS).
  /// The method finds intersections of the ray with the LCMS, 
  /// where the ray is R(t)=r0+rd*t,  and t>=0.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param r0 is the Cartesian coordinates of the ray origin, must be outside with respect to the LCMS,
  /// @param rd is the direction of the ray in Cartesian coordinates,
  /// @param entryPoint this vector at return contains the Cartesian coordinates of the ray entry point into the plasma,
  /// @param exitPoint this vector at return contains the Cartesian coordinates of the ray exit point from the plasma,
  /// @return  if the function succeeds, the return value is \b non-zero; \b zero otherwise.
  WINDLLEXPORT int MCgetRayIntersectionPoints(MC_HANDLE mConf,const double *r0,const double *rd,double *entryPoint,double *exitPoint); 
  
  /// The method returns the toroidal flux in Wb.
  /// @param s is the normalized toroidal flux.
  WINDLLEXPORT double MCFlux  (MC_HANDLE mConf,double s);  

  /// The function returns the normalized poloidal flux, where s is the normalized toroidal flux.
  WINDLLEXPORT double MCtorFlux2polFlux(MC_HANDLE mConf,double s);
 
  /// The method returns the poloidal flux in Wb.
  /// @param s is the normalized toroidal flux.
  WINDLLEXPORT double MCPoloidalFlux(MC_HANDLE mConf,double s);

  /// The method returns the normalized toroidal flux, where s is the normalized poloidal flux.
  WINDLLEXPORT double MCsToroidal(MC_HANDLE mConf,double spoloidal);

  /// The function returns the iota (=1/q) value , where s is the normalized toroidal flux.
  WINDLLEXPORT double MCiota(MC_HANDLE mConf,double s);

  /// The function returns the diota/ds value, where s is the normalized toroidal flux.
  WINDLLEXPORT double MCiotaPrime(MC_HANDLE mConf,double s);

  /// The method returns the poloidal current in [A].
  /// @param s is the normalized toroidal flux.
  WINDLLEXPORT double MCIp    (MC_HANDLE mConf,double s);   // poloidal current [A]
  /// The method returns the toroidal current in [A].
  /// @param s is the normalized toroidal flux.
  WINDLLEXPORT double MCIt    (MC_HANDLE mConf,double s);   // toroidal current [A]
  /// The method returns the maximum value of magnetic field on the surface \b s.
  WINDLLEXPORT double MCBmax  (MC_HANDLE mConf,double s);

  /// The method returns  \<|grad(s)|\>  \f$\left\langle\left|\nabla s\right|\right\rangle\f$
  WINDLLEXPORT double MCgradsAvrg(MC_HANDLE mConf,double s);
  /// The method returns  \<(grad(s))^2\> \f$\left\langle\left|\nabla s\right|^2\right\rangle\f$
  WINDLLEXPORT double MCgrads2Avrg(MC_HANDLE mConf,double s);
  /// The method returns  \<|grad(rho)|\>  \f$\left\langle\left|\nabla \rho\right|\right\rangle\f$, where rho=sqrt(s).
  WINDLLEXPORT double MCgradRhoAvrg(MC_HANDLE mConf,double s);
  /// The method returns  \<(grad(rho))^2\> \f$\left\langle\left|\nabla \rho\right|^2\right\rangle\f$, where rho=sqrt(s).
  WINDLLEXPORT double MCgradRho2Avrg(MC_HANDLE mConf,double s);

  /// The method returns the fraction of trapped particles on the surface s. 
  ///  The averaged fraction of trapped particles on a magnetic flux
  ///  surface is only defined by the magnetic field strength:
  ///
  ///    FT = 1 - 3/4 * <B^2/Bmax^2> * int (0->1) xdx/<g(x)>
  ///
  ///  where <g(x)> is the flux surface average of sqrt(1-x*B/Bmax)
  ///  with Bmax being the maximum of the magnetic field strength and
  ///  <A> = int int A*jacobian*dth*dph  /  int int jacobian*dth*dph
  ///  is the definition of the flux surface average in magnetic coordinates.
  WINDLLEXPORT double MCftrapped(MC_HANDLE mConf,double s);
  
  /// The method returns the bootstrap current geometric factor for \f$1/\nu\f$ transport, 
  /// using s as a flux surface label, where s is the normalized toroidal flux.  
  /// g2 is defined as B*grad(g2/B^2) = B x grad(s) * grad(1/B^2),
  /// g4 is defined as B*grad(g4/xi) = B x grad(s) * grad(1/xi).
  /// @param s is the normalized toroidal flux
  /// \sa P.Helander et al,Phys.Plasmas 18,092505(2011), http://link.aip.org/link/doi/10.1063/1.3633940
  /// \note  At the first call the method creates a look-up table using multithreading technique to spread 
  /// calculation among several CPU cores; this table is used for interpolating over the flux surface label.
  /// \n Use the following methods to set parameters:
  /// FbsSetXeffParam(), FbsSetTracingParam(), FbsSetIotaParam(), FbsSetMagnMomentParam()
  WINDLLEXPORT double MCFbs(MC_HANDLE mConf,double s);

  /// The method returns the bootstrap current geometric factor for \f$1/\nu\f$ transport, 
  /// using VMEC definition of minor radius r=a*sqrt(s), where s is the normalized toroidal flux.  
  /// g2 is defined as B*grad(g2/B^2) = B x grad(r) * grad(1/B^2),
  /// g4 is defined as B*grad(g4/xi) = B x grad(r) * grad(1/xi).
  /// @param s is the normalized toroidal flux
  /// \sa P.Helander et al,Phys.Plasmas 18,092505(2011), http://link.aip.org/link/doi/10.1063/1.3633940
  /// \note  At the first call the method creates a look-up table using multithreading technique to spread 
  /// calculation among several CPU cores; this table is used for interpolating over the flux surface label.
  /// \n Use the following methods to set parameters:
  /// FbsSetXeffParam(), FbsSetTracingParam(), FbsSetIotaParam(), FbsSetMagnMomentParam().
  /// \note  MCFbsVmec = MCFbs*a*sqrt(s)/(2*s)
  WINDLLEXPORT double MCFbsVmec(MC_HANDLE mConf,double s);

  /// The method returns flux surface average of g2, 
  /// where g2 is defined as B*grad(g2/B^2) = B x grad(s) * grad(1/B^2),
  /// s is the normalized toroidal flux.  
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] s is the normalized toroidal flux.
  /// \sa P.Helander et al,Phys.Plasmas 18,092505(2011), http://link.aip.org/link/doi/10.1063/1.3633940
  /// \note  At first call the method creates a look-up table using multithreading technique to spread 
  /// calculation among several CPU cores; this table is used for interpolating over the flux surface label.
  /// \n Use the following methods to set parameters:
  /// MCFbsSetSlabelParam(),MCFbsSetXeffParam(), MCFbsSetTracingParam(), MCFbsSetIotaParam(), MCFbsSetMagnMomentParam()
  WINDLLEXPORT double MCFbsg2(MC_HANDLE mConf,double s);

  /// The method returns flux surface average of u, 
  /// where u is defined as B*grad(u) = -B x grad(s) * grad(1/B^2),
  /// s is the normalized toroidal flux.  
  /// @param s is the normalized toroidal flux
  WINDLLEXPORT double MCFbsu(MC_HANDLE mConf,double s); 
  /// The method returns flux surface average of u*B^2, 
  /// where u is defined as B*grad(u) = -B x grad(s) * grad(1/B^2),
  /// s is the normalized toroidal flux.  
  /// @param s is the normalized toroidal flux
  WINDLLEXPORT double MCFbsuB2(MC_HANDLE mConf,double s); 
  /// The method returns flux surface average of u^2*B^2, 
  /// where u is defined as B*grad(u) = -B x grad(s) * grad(1/B^2),
  /// s is the normalized toroidal flux.  
  /// @param s is the normalized toroidal flux
  WINDLLEXPORT double MCFbsu2B2(MC_HANDLE mConf,double s); 
  /// The method returns flux surface average of B^2 
  /// @param s is the normalized toroidal flux
  WINDLLEXPORT double MCFbsB2(MC_HANDLE mConf,double s); 


  /// The method returns flux surface average of g4, 
  /// where g4 is defined as B*grad(g4/xi) = B x grad(s) * grad(1/xi) and
  /// xi = sqrt(1 - mu*B/Bmax),
  /// s is the normalized toroidal flux.  
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] s is the normalized toroidal flux.
  /// \sa P.Helander et al,Phys.Plasmas 18,092505(2011), http://link.aip.org/link/doi/10.1063/1.3633940
  /// \note  At first call the method creates a look-up table using multithreading technique to spread 
  /// calculation among several CPU cores; this table is used for interpolating over the flux surface label.
  /// \n Use the following methods to set parameters:
  /// MCFbsSetSlabelParam(),MCFbsSetXeffParam(), MCFbsSetTracingParam(), MCFbsSetIotaParam(), MCFbsSetMagnMomentParam()
  /// Default parameters are the following:
  /// \code
  ///      dphi = 0.017453;     // twopi/360 = 0.017453 -> 1degree; 
  ///      doAccuracyTest =true;
  ///      avoidResonances=false;
  ///      useRationalIota=false;
  ///      iotaDenominator=150;
  ///      turns=100;
  ///      yMax=15;
  ///      numMuLevels=257;
  ///      size = 101;
  ///      smin = .0025;
  ///      smax = .9;
  /// \endcode

  WINDLLEXPORT double MCFbsg4(MC_HANDLE mConf,double s, double mu);

  /// The method sets parameters for calculating the bootstrap current geometric factor.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] smin is the minimal value of s, where s is the normalized toroidal flux,
  /// @param[in] smax is the maximal value of s,
  /// @param[in] size is the number of points for creating look-up table
  WINDLLEXPORT void MCFbsSetSlabelParam(MC_HANDLE mConf,double smin, double smax, int size); 
  /// The method sets parameters for calculating the bootstrap current geometric factor.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] xmin is the minimal value of sqrt(s), where s is the normalized toroidal flux,
  /// @param[in] xmax is the maximal value of sqrt(s),
  /// @param[in] size is the number of points for creating look-up table
  WINDLLEXPORT void MCFbsSetXeffParam   (MC_HANDLE mConf,double xmin, double xmax, int size);   
  /// The method sets parameters for calculating the bootstrap current geometric factor.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] turns is the maximal number of toroidal turns for a magnetic field line following,
  /// @param[in] doAccuracyTest -- if true then integrating along magnetic field line is stopped 
  ///                          when relative accuracy is less then 2% 
  /// @param dphi is the toroidal angle step in radian for the magnetic field line following.
  WINDLLEXPORT void MCFbsSetTracingParam(MC_HANDLE mConf,int turns, bool doAccuracyTest, double dphi);  // twopi/360 == 0.017453
  /// The method sets parameters for calculating the bootstrap current geometric factor.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] avoidResonances -- if true then change of flux surface label value to avoid 
  ///                           iota resonances is attempted,
  /// @param[in] useRationalIota -- if true then rational iota aproximation is used, 
  /// @param[in] iotaDenominator is the denominator guess for the rational iota aproximation.
  WINDLLEXPORT void MCFbsSetIotaParam   (MC_HANDLE mConf,bool avoidResonances, bool useRationalIota, int iotaDenominator); 
  /// The method sets parameters for calculating the bootstrap current geometric factor.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] muLevels is the number of integration points over magnetic moment,
  /// @param[in] yMax is the new integration limit after variable change from the magnetic moment
  ///    mu to y with mu = tanh(y),
  ///    see N.Nakajima et al, Nuclear Fusion,29,605(1989), 
  ///      http://dx.doi.org/10.1088/0029-5515/29/4/006
  WINDLLEXPORT void MCFbsSetMagnMomentParam (MC_HANDLE mConf,int muLevels, double yMax); 


  /// The function creates 3d-mesh in cylindrical coordinates.
  /// The function tabulates the magnetic field, the flux surface label, and grad(s);
  /// the \b periodicity (<b>B</b>(r,fi,z)=<b>B</b>(r,fi+period,z)) and
  ///  \b symmetry (B<sub>r</sub>(r,fi,z)=-B<sub>r</sub>(r,-fi,-z)) of the device are used.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] dr is the mesh step in r-direction.
  /// @param[in] dz is the mesh step in z-direction.
  /// @param[in] dfi is the mesh step in fi-direction, where \e fi is the cylindrical angle in radians.
  /// @return  if the function succeeds, the return value is \b non-zero; \b zero otherwise.
  /// \sa the detailed description of createMeshUsingSymmetry()
  WINDLLEXPORT int    MCcreateMeshUsingSymmetry(MC_HANDLE mConf,double dr,double dz,double dfi);
  /// The function creates 3d-mesh in cylindrical coordinates.
  /// The function tabulates the magnetic field, the flux surface label, and grad(s);
  /// the periodicity of the device is taken into account.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] dr is the mesh step in r-direction.
  /// @param[in] dz is the mesh step in z-direction.
  /// @param[in] dfi is the mesh step in fi-direction, where \e fi is the cylindrical angle in radians.
  /// @return  if the function succeeds, the return value is \b non-zero; \b zero otherwise.
  /// \sa the detailed description of createMesh()
  WINDLLEXPORT int    MCcreateMesh             (MC_HANDLE mConf,double dr,double dz,double dfi);
  /// The function traces the Last Closed Magnetic Surface (LCMS).
  /// The function finds intersections of the ray with the LCMS, where the
  /// ray is R(t)=r0+rd*t,  and t>=0.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] r0 is the origin of the ray, must be outside with respect to the LCMS.
  /// @param[in] rd is the direction of the ray
  /// @param[in] entryPoint this vector at return contains the coordinates of the ray entry point into the plasma.
  /// @return  if the function succeeds, the return value is \b non-zero; \b zero otherwise.
  WINDLLEXPORT int    M3DgetRayEntryPoint(MC_HANDLE mConf,const double *r0,const double *rd,double *entryPoint);
  /// The function returns cartesian B-vector.
  /// 4-points Lagrange interpolation is used.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] xyz is the cartesian coordinates
  /// @param[out] B is the cartesian B-field vector
  /// @param[out] grads is the cartesian grads vector
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  /// \b s>1000 means that the point lies outside the \b sMax. \n
  ///  The extrapolated region is needed in order to do interpolation near LCMS, see also setsMax(),getsMax().
  WINDLLEXPORT double M3DgetBGradsxyz(MC_HANDLE mConf,const double *xyz,double *B, double *grads);
  /// The function returns cartesian B-vector.
  /// 4-points Lagrange interpolation is used.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] xyz is the cartesian coordinates
  /// @param[out] B is the cartesian B-field vector
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  /// \b s>1000 means that the point lies outside the \b sMax. \n
  ///  The extrapolated region is needed in order to do interpolation near LCMS, see also setsMax(),getsMax().
  WINDLLEXPORT double M3DgetBxyz     (MC_HANDLE mConf,const double *xyz,double *B);
  /// The function returns cartesian grad(|B|) vector.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] xyz is the cartesian coordinates
  /// @param[out] gradB is the cartesian grad(|B|) vector
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  /// \b s>1000 means that the point lies outside the \b sMax. \n
  ///  The extrapolated region is needed in order to do interpolation near LCMS, see also setsMax(),getsMax().
  WINDLLEXPORT double M3DgetGradBxyz (MC_HANDLE mConf,const double *xyz,double *gradB);
  /// The method returns Cartesian vectors B, grad(|B|), and grad(s).
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] xyz is the cartesian coordinates
  /// @param[out] B is the cartesian B-field vector
  /// @param[out] gradB is the cartesian grad(|B|) vector
  /// @param[out] gradS is the cartesian grad(s) vector
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  /// \b s>1000 means that the point lies outside the \b sMax. \n
  ///  The extrapolated region is needed in order to do interpolation near LCMS, see also setsMax(),getsMax().
  WINDLLEXPORT double M3DgetBandGradBxyz(MC_HANDLE mConf,const double *xyz,double *B,double *gradB,double *gradS);
  /// The method returns cylindrical vectors B, grad(|B|)=(dB/dr, dB/dfi/r, dB/dz), and grad(s).
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] cyl is the cylindrical coordinates
  /// @param[out] B is the cylindrical B-field vector
  /// @param[out] gradB is the cylindrical grad(|B|) vector
  /// @param[out] gradS is the cylindrical grad(s) vector
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  /// \b s>1000 means that the point lies outside the \b sMax. \n
  ///  The extrapolated region is needed in order to do interpolation near LCMS, see also setsMax(),getsMax().
  WINDLLEXPORT double M3DgetBandGradB(MC_HANDLE mConf,const double *cyl,double *B,double *gradB,double *gradS);
  /// The function finds the flux label \b s which corresponds to the point in space.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] xyz is the input parameter with the cartesian coordinates of the point of interest.
  /// @return \b s (normalized toroidal flux), which normally satisfies \b 0<=s<=1; \n
  /// <b>1<s<sMax</b> means that the point lies in the extrapolated region,
  ///    retuned parameters are not reliable here; \n
  /// \b s>1000 means that the point lies outside the \b sMax. \n
  ///  The extrapolated region is needed in order to do interpolation near LCMS, see also setsMax().
  WINDLLEXPORT double M3Dxyz2s       (MC_HANDLE mConf,const double *xyz);
  
  //***************************************************************
  // get B and derivatives dB/dx,dB/dy,dB/dz at the point xyz
  //  Input:
  //    Vector3d xyz -- cartesian coordinates
  //  Output:
  //    Vector3d dBdx -- dB/dx, where dBdx(1)=dBx/dx, dBdx(2)=dBy/dx, dBdx(3)=dBz/dx
  //    Vector3d dBdy -- dB/dy
  //    Vector3d dBdz -- dB/dz
  //    Vector3d B    -- cartesian coordinates of B
  //    Vector3d grads -- cartesian coordinates of grad(s)
  //  Return:   s
  //
  WINDLLEXPORT double M3DgetdB_Gradsxyz(MC_HANDLE mConf,const double *xyz,double *B,double *dBdx,double *dBdy,double *dBdz,double *gradS);
  
  //***************************************************************
  // get B and derivatives dB/dx,dB/dy,dB/dz at the point xyz
  //  Input:
  //    Vector3d xyz -- cartesian coordinates
  //  Output:
  //    Vector3d dBdx -- dB/dx, where dBdx(1)=dBx/dx, dBdx(2)=dBy/dx, dBdx(3)=dBz/dx
  //    Vector3d dBdy -- dB/dy
  //    Vector3d dBdz -- dB/dz
  //    Vector3d B    -- cartesian coordinates of B
  //    Vector3d grads -- cartesian coordinates of grad(s)
  //  Return:   s
  //
  WINDLLEXPORT double MCgetdB_Gradsxyz(MC_HANDLE mConf,const double *xyz,double *B,double *dBdx,double *dBdy,double *dBdz,double *gradS);


  //***************************************************************
  // get gradS and partial derivatives at point xyz, all parameters in cartesian coordinates 
  //  Input:
  //    Vector3d xyz -- cartesian coordinates
  //  Output:
  //    Vector3d gradS  -- gradient(s) 
  //    Vector3d dGdx -- dgrads/dx, where dgradsdx(1)=dgradsx/dx, dgradsdx(2)=dgradsy/dx, dgradsdx(3)=dgradsz/dx
  //    Vector3d dGdy -- dgrads/dy
  //    Vector3d dGdz -- dgrads/dz
  //  Return:
  //    s -- normalized toroidal flux
  WINDLLEXPORT double M3DgetdGradsxyz(MC_HANDLE mConf,const double *xyz,double *gradS,double *dGdx,double *dGdy,double *dGdz);

  //*********************************************************************
  /// The method returns the effective helical ripple for \f$1/\nu\f$ transport.  
  /// \sa V.V.Nemov, S.V.Kasilov, W. Kernbichler, and M.F.Heyn, Physics of Plasmas, vol. 6, 4622(1999),
  ///   http://link.aip.org/link/PHPAEN/v6/i12/p4622/s1 
  /// \sa 28th EPS, http://epsppd.epfl.ch/Madeira/html/pdf/P5.055.pdf , EGA Vol.25A (2001) 1985-1988
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in]  s is the normalized toroidal flux
  /// @return  eps_eff  using VMEC definition of minor radius r=a*sqrt(s),
  WINDLLEXPORT double MCepsEff(MC_HANDLE mConf,double s);


  /// Get B and gradients in cartesian coordinates
  ///  @param[in] Vector3d magCoord(s,theta,phi) is a point in magnetic coordinates,
  ///  @param[out]  B  is the magnetic field  in cartesian coordinates
  ///  @param[out]  gradB  is the grad(|B|) in cartesian coordinates
  ///  @param[out]  gradS  is the grad(s)
  ///  @param[out]  gradTh is the grad(theta)
  ///  @param[out]  gradPh is the grad(phi) 
  WINDLLEXPORT void MCgetBandGradientsxyz(MC_HANDLE mConf,const double *magCoord, double *B, 
                               double *gradB, double *gradS, double *gradTh, double *gradPh);




  /// The method returns the local shear.
  /// @param[in] Vector3d magCoord(s,theta,phi) is a point in magnetic coordinates,
  /// @return  local_shear = -H*curl(H), where H = grad(s)/abs(grad(s) ^ B/abs(B), ^ is the cross product
  WINDLLEXPORT double MCgetLocalShear(MC_HANDLE mConf,const double *magCoord);

  /// The methods transform the current magnetic configuration to new mesh on s.
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] newLastSurface is the value of new last closed magnetic surface 
  ///   in inits of current s (normalized toroidal flux),
  ///   the old last surface will be at s_new = 1/sqrt(newLastSurface).
  /// @param[in] fname is the filename under which the configuration is saved.
  /// @return  if the function succeeds, the return value is \b non-zero; \b zero otherwise.
  /// \attention Current truncation level and value of magnetic field are used;
  WINDLLEXPORT int MCsetLCMS(MC_HANDLE mConf, double newLastSurface); //, const char * fname); 

  /// The methods writes the current magnetic configuration into the file fname.
  WINDLLEXPORT int MCwrite(MC_HANDLE mConf, const char * fname);

  ///****************************************************************************
  /// The method returns the coefficients needed for 
  /// the poloidal flux equation in Astra transport code
  /// @param[in] mConf is the handle of the magnetic configuration.
  /// @param[in] sqrts is the sqrt(s), where s is the normalized toroidal flux. 
  /// @param[out] r = a*sqrt(s) is the effective plasma radius, where a is the minor radius.
  /// @param[out] gradr2Avr is the flux surface average of (grad(r))^2, \<(grad(r))^2\> \f$\left\langle\left|\nabla r\right|^2\right\rangle\f$
  /// @param[out] J = Ipol/Ipol(a)
  /// @param[out] G2 = h*S11 > 0 is the G2 coefficient in the poloidal flux equation,
  ///       h is +1 for the right handed (s,theta,phi) coordinate system and is -1 for the left handed one.
  /// @param[out] hVprime = h*dV/dr > 0
  /// @param[out] B0 = h*TorFlux(a) / (pi*a^2)
  /// @param[out] R0 = mu0*Ipol(a) / (2pi*B0)
  /// @param[out] h  is +1 for the right handed (s,theta,phi) coordinate system and is -1 for the left handed one.
  /// \note 
  /// Python user.
  /// declare double as c_double() and use byref(): 
  /// mconf.MCgetCoeffForAstraCode(mc, sqrts, byref(r), byref(gradr2Avr), byref(J), byref(G2), byref(hVprime), byref(B0), byref(R0), byref(h))
  /// https://docs.python.org/3.3/library/ctypes.html#passing-pointers-or-passing-parameters-by-reference
  /// https://docs.python.org/3/library/ctypes.html#module-ctypes
  WINDLLEXPORT void MCgetCoeffForAstraCode(MC_HANDLE mConf, double sqrts, double *r, double *gradr2Avr, double *J, double *G2, double *hVprime, double *B0, double *R0, double *h);

  /// The method calculates the Susceptance matrix in (Flux,theta,phi)-coordinates.
  /// See definition of the Susceptance matrix in P.I.Strand,W.A.Houlberg, Physics of Plasmas, \b 8, 2782(2001)
  /// The Susceptance matrix in (s,theta,phi)-coordinates can be expressed 
  /// through relation S_ij(Flux)=S_ij(s)* dFlux/ds
  /// @param s is the normalized toroidal flux, s = Flux/Flux_max
  /// @return S<sub>11</sub> element of susceptance matrix.
  ///  See also  S12(), S21(), S22().
  WINDLLEXPORT double MCS11(MC_HANDLE mConf,double s);
  /// Susceptance matrix: S12 in (Flux,theta,phi)-coordinates.
  WINDLLEXPORT double MCS12(MC_HANDLE mConf,double s);
  /// Susceptance matrix: S21 in (Flux,theta,phi)-coordinates.
  WINDLLEXPORT double MCS21(MC_HANDLE mConf,double s);
  /// Susceptance matrix: S22 in (Flux,theta,phi)-coordinates.
  WINDLLEXPORT double MCS22(MC_HANDLE mConf,double s);
  /// The method returns B_mn(s) for Boozer coordinates.  
  WINDLLEXPORT double MCBoozerBmn(MC_HANDLE mConf,double s,int m,int n);
  /// The method returns the major radius.
  WINDLLEXPORT double MCR0(MC_HANDLE mConf);
  /// The method returns the elongation:
  /// e = (r(s)/R0)^2/(B10(s)/B00(s))^2, where Bmn are in Boozer coordinates, see MCBoozerBmn()
  WINDLLEXPORT double MCelongation(MC_HANDLE mConf,double s);
 
  ///The method forces to use Jacobian=dX/dFlux*(dX/dtheta^dX/dphi) for calculating B field if yesNo>0 
  WINDLLEXPORT void MCuseMixedProductForJacobian(MC_HANDLE mConf,int yesNo);
  /// The method sets the direction of B with respect
  /// to cylindrical angle of the right handed coordinate system.
  /// @param sign is the positive or negative integer number.
  WINDLLEXPORT void  MCsetBdirection(MC_HANDLE mConf, int sign);
  /// The method returns the direction of B with respect
  /// to cylindrical angle of the right handed cylindrical coordinate system.
  WINDLLEXPORT int MCgetBdirection(MC_HANDLE mConf);
  /// The method sets the original direction of B with respect
  /// to cylindrical angle of the right handed coordinate system.
  /// The original direction is the direction of B stored in the equilibrium file. 
  WINDLLEXPORT void MCrestoreBdirection(MC_HANDLE mConf);

  
#ifdef __cplusplus
}
#endif


#endif
