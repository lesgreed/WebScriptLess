// this is doxygen style html-page

/*

useful things

<div>
<table border="0">
<tr> <td width="100">under construction</td> <td> \image html undercon.gif " "  </td> </tr>
</table>
</div>

<!-- [top][bootom][back] -->
<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

<!-- [back] -->
<td>[ <a href="index.html" onClick="history.back();return false;" class="el"><b>back</b></a> ] </td>

*/
//************************************************************************************

/* \example w7xFormat.f90
 * This example illustrates format of bc-file.
 */

/** @mainpage notitle

<h1>MConf Library</h1>

- \ref WhatItDoes
- \ref WhatItDoesNot
- \ref FluxCoordinates
  - \ref VmecCoordinates
    - \ref PESTCoordinates
  - \ref BoozerCoordinates
  - \ref TokamakSymmetryFluxCoordinates
- \ref CoordinateTransformation
  - \ref VMEC2Boozer
  - \ref Cyl2Mag 
    - \ref NewtonMethod

<br>

- \ref Gettingstarted  
  - \ref MinimalExample
  - \ref ObjectCanBeCopied
  - <a class="el" href="examples.html">Other examples</a>
- \ref RelatedProjects
  - <a class="el" href="../docs-mcviewer/index.html">MCviewer - Magnetic Configuration viewer</a>
  - \ref Travis
  - \ref NBIviewer
- \ref References
- \ref SoftwareLicense "Software License"

<br>

- \ref MixedLanguage "Mixed-language programming"
  - \ref cpp2for.cpp "Fortran calls MConf (cpp2for.cpp)"
- \ref MatlabInterface "Matlab (or Python) Interface"
  - \ref mconf_matlab.h "Functions"

<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

\b MConf stands for Magnetic Configuration. The goal of this library is to
provide fast and convenient tools for coordinate transformations between
flux coordinates and real space coordinates. Here you will find the
description of CStconfig, CRayTrace, and C3dMesh classes, which comprise the
package MConf. An C++ <a class="el" href="examples.html">examples</a> 
and example cpp2for.cpp of how to write
interface C-functions (see also \ref MixedLanguage) in order to call C++ methods from within FORTRAN programs
are available. 

//************************************************************************************
\section  WhatItDoes  What it does

The input file for the MConf is the Boozer magnetic configuration file which
is the result of VMEC code and JMC/xbooz_xform codes. However, MConf
understands VMEC wout-file, in this case MConf uses VMEC coordinates. Magnetic
configuration can be loaded also from \ref MCDB "Magnetic Configuration Database".

MConf makes available numerous information about magnetic configuration of a
stellarator or tokamak: magnetic field, Jacobian, iota, trapped particle
fraction, minimum and maximum magnetic field on a flux surface, volume inside
flux surface, and so on. The package contains programs for straight line
tracing through flux surfaces. MConf can also create 3d-mesh in cylindrical
coordinates on which it tabulates magnetic field, flux surface label \a s,
grad(s) and provides functions for interpolation. The mesh and magnetic
configuration can be stored in a file for future loading.

Due to object oriented approach this package easy to use in other application;
for example, the visualization program <a class="el" href="../docs-mcviewer/index.html">MCviewer</a>,
which displays various aspects of a magnetic configuration, has
been created using the library MConf. The library is written in C++, interface
routines for using MConf from within FORTRAN are provided. The library can be
used in Multithreaded Applications, see 
<a class="el" href="examples.html">examples</a>. \n

In short, MConf features are:

\addindex Features

- load function, supported formats:
    - \ref W7X_Format "W7-X" by J. Geiger
    -  LHD \ref References "[6]" 
    -  VMEC wout text file \ref References "[4]"
    -  \ref GEQDSK "EFIT G EQDSK"
    -  proprietary binary format (bin4, bin8)
- save function, supported formats:
    - \ref W7X_Format "W7-X"
    -  proprietary binary format (bin4, bin8)
- data base operation 
    - load/delete/insert magnetic configuration from/to 
    \ref MCDB "Magnetic Configuration Data Base"
- functions for:  
    - coordinate transformation from real to flux coordinates: forth and back
    - coordinate transformation from VMEC to Boozer coordinates
    - inquire position (inside or outside LCMS or relative to the given contour)
    - surface quantity: \b B, B<sub>min</sub>, B<sub>max</sub>, Jacobian, iota,
        V, V',I<sub>tor</sub> , I<sub>pol</sub>, grad(r<sub>eff</sub>), Fourier
        coefficients
    - remeshing to decrease number of flux surfaces and truncate spectrum
    - straight line ray tracing to find intersection of a ray with the flux
        surface
    - 3d-mesh generator which tabulates magnetic field, flux surface
        label \a s, grad(r<sub>eff</sub>) on a 3d-mesh in cylindrical coordinates
        and provides functions for interpolation. The mesh and its functions can be
        used time consuming computing (NBI modeling(including Monte-Carlo),
        ECE and ECRH Ray tracing)
    - mesh generator can be started in several threads thus utilizing extra cores
        in modern processors, for example, Intel Core&tm; 2 Duo or Intel
        Core&tm; 2 Quad, see MConf::C3dMesh::createMesh()
    - function for saving and loading 3d-mesh or magnetic configuration file
        in binary format for faster loading in future use
    - calculating 
       - effective helical ripple MConf::C3dMesh::epsEffVmec()
       - gradB drift velocity of trapped particles MConf::C3dMesh::Gv()
       - bootstrap current geometric factor MConf::C3dMesh::FbsVmec()


<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

//************************************************************************************
\section  WhatItDoesNot  What it doesn't do

MConf doesn't reconstruct or calculate equilibrium,
but uses VMEC(VMEC based) or EFIT EQDSK equilibrium.\n

MConf can't work outside the last closed magnetic surface (LCMS);
LCMS is the surface with the flux suface label \a s = 1,
where \a s is the normalized toroidal flux.

<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

//************************************************************************************
\section  FluxCoordinates Flux Coordinates
Flux coordinates concept is very fundamental and convenient tools for describing complex 
geometry of plasma in tokamaks or stellarators \ref References "[1]". These coordinates are aligned with 
the magnetic field structure that represents the set of nested magnetic (or flux) surfaces 
of a toroidal type. The surfaces are enumerated by the flux surface label or effective 
plasma radius; usually that is a quantity derived from the poloidal or toroidal magnetic 
flux value enclosed by the flux surface. To define position on a flux surfaces the 
two angle-like variables \f$(\theta,\varphi)\f$  are used. The flux coordinates in which these angles are 
chosen in such a way that the magnetic field lines are straight lines on 
the  \f$(\theta,\varphi)\f$  plane are called magnetic coordinates. 
//************************************************************************************
\subsection  VmecCoordinates Vmec Coordinates

\par

MConf can interpret the VMEC \ref References "[4]" wout-file version 6.20 and higher.
Present version of the library handle VMEC equilibria in stellarator symmetric mode only. 
VMEC toroidal angle coincides with cylindrical angle and additionally VMEC provides 
the stream function \f$\lambda\f$ that helps to build straight field line 
coordinates system    

\f[
\mathbf{B}=\frac{\psi_{lcms}}{2\pi}\left(\nabla s\times\nabla\theta^{*}+
           \iota(s)\nabla\varphi\times\nabla s\right)
\f]

with \f$\theta^{*}=\theta+\lambda(s,\theta,\varphi)\f$

The flux coordinates \f$(s,\theta,\varphi)\f$ 
are defined by giving the cylindrical coordinates as 
functions of the flux coordinates

\f{eqnarray*}
     r(s,\theta,\varphi) &=& \sum_{m=0}^M\sum_{n=-N}^{N}R_{mn}(s)\cos\left(m\theta-nN_{p}\varphi\right)\\
  \phi(s,\theta,\varphi) &=& \varphi\\
	   z(s,\theta,\varphi) &=& \sum_{m=0}^{M}\sum_{n=-N}^{N}Z_{mn}(s)\sin\left(m\theta-nN_{p}\varphi\right)\\
	\lambda(s,\theta,\varphi) &=& \sum_{m=0}^{M}\sum_{n=-N}^{N}\lambda_{mn}(s)\sin\left(m\theta-nN_{p}\varphi\right)
\f}
 
The Jacobian and magnetic field are calculated using formulas

\f{eqnarray*}
	J(s,\theta,\varphi) &=& \sum_{m=0}^{M}\sum_{n=-N}^{N}g_{mn}(s)\cos\left(m\theta-nN_{p}\varphi\right)\\	
	\mathbf{B}(s,\theta,\varphi) &=&\frac{\psi_{lcms}}{2\pi J}\left[\frac{\partial\mathbf{X}}{\partial\varphi}
  \left(1+\frac{\partial\lambda}{\partial\theta}\right)+
  \frac{\partial\mathbf{X}}{\partial\theta}\left(\iota-
  \frac{\partial\lambda}{\partial\varphi}\right)\right]
	
\f}

B-field line \f$(s_{0},\theta,\varphi)\f$ going through the point 
\f$(s_{0},\theta_{0},\varphi_{0})\f$ is defined by expression

\f[
\theta(\varphi)=\theta_{0}+\iota(s_{0})(\varphi-\varphi_{0}) - 
[\lambda(s_{0},\theta,\varphi)-\lambda(s_{0},\theta_{0},\varphi_{0})]
\f]

Both methods can be applied to follow field line:

- MConf::C3dMesh::magCoordFieldLine()  uses cylindrical angle,
- MConf::C3dMesh::mixCoordFieldLine()  uses cylindrical angle.

\note 

Recently J. Geiger and I have discovered that even VMEC 2000 produces slightly 
different output.

Reading the wout-file created by VMEC in use by J. Geiger.  

\code
! full-mesh quantities
  read(iunit,*)(iotaf(j),presf(j),phipf(j),phi(j),jcuru(j),jcurv(j),j=1,ns)
! half-mesh quantities 
  read(iunit,*)(iotas(j),mass(j),pres(j),beta_vol(j),phip(j),buco(j),bvco(j), &
    vp(j),overr(j),specw(j),j=2,ns)

\endcode

Reading the file created by VMEC in use by PPPL.

\code

! half-mesh quantities (except phi, jcuru, jcurv which are on full-mesh)
  read(iwout,*)(iotas(j),mass(j),pres(j),beta_vol(j),phip(j),buco(j),bvco(j), phi(j), &
     vp(j),overr(j),jcuru(j),jcurv(j),specw(j),j=2,ns)
\endcode

<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

//************************************************************************************
\subsubsection  PESTCoordinates PEST coordinates

The straight-field-line coordinates system that uses the cylindrical angle \f$\phi\f$ as the 
toroidal angle in tokamak literature is called PEST \ref References "[7]" coordinates. 

In case of VMEC input file these coordinates, \f$(s,\theta^{*},\phi)\f$ , can be easily  
build from \ref VmecCoordinates and its stream function \f$\lambda\f$ using transformation 

\f$\theta^{*}=\theta+\lambda(s,\theta,\phi)\f$

The following methods perform transformation between the PEST and VMEC coordinates:

MConf::CStconfig::Vmec2Pest() <br>
MConf::CStconfig::Pest2Vmec() <br>

The method returns contravariant-basis vectors \f$(\nabla s, \nabla \theta^{*}, \nabla \phi)\f$ in cylindrical coordinates:

MConf::CStconfig::SFLcontraBasis() <br>

<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

//************************************************************************************
\subsection  BoozerCoordinates Boozer Coordinates


In Boozer coordinates, the contravariant and covariant components of the
magnetic field are defined through the expressions \ref References "[1,2]":

\f[
\mathbf{B}=\frac{1}{2\pi}\left(\nabla\psi\times\nabla\theta+
           \iota(\psi)\nabla\varphi\times\nabla\psi\right)
\f]


\f[
\mathbf{B}=\frac{\mu_{0}}{2\pi}\left(I_{pol}(\psi)\nabla\varphi+
           I_{tor}(\psi)\nabla\theta+\beta_{*}\nabla\psi\right)
\f]  
where \f$(\psi,\theta,\varphi)\f$ are the toroidal flux, 
poloidal, and toroidal angles; \f$0\leq\psi\leq\psi_{lcms}\f$ (\f$\psi=0\f$ at the magnetic axis),
\f$0\leq\theta\leq2\pi\f$, \f$0\leq\varphi\leq2\pi\f$.

The magnetic coordinate system \f$(\psi,\theta,\varphi)\f$ is defined by giving 
the cylindrical coordinates as functions of magnetic coordinates, i.e. the 
cylindrical coordinates of a point on a magnetic surface
and the magnetic field value are:
 
\f{eqnarray*}
  r(\psi,\theta,\varphi) &=& \sum_{m=0}^M\sum_{n=-N}^{N}R_{mn}(\psi)\cos\left(m\theta-nN_{p}\varphi\right)\\
  \phi(\psi,\theta,\varphi) &=& \varphi-\frac{2\pi}{N_{p}}\sum_{m=0}^{M}\sum_{n=-N}^{N}\Phi_{mn}(\psi)\sin\left(m\theta-nN_{p}\varphi\right)\\
	z(\psi,\theta,\varphi) &=& \sum_{m=0}^{M}\sum_{n=-N}^{N}Z_{mn}(\psi)\sin\left(m\theta-nN_{p}\varphi\right)\\
	B(\psi,\theta,\varphi) &=& \sum_{m=0}^{M}\sum_{n=-N}^{N}B_{mn}(\psi)\cos\left(m\theta-nN_{p}\varphi\right)
\f}

where \f$N_p\f$ is the number of field periods.

Poloidal and toroidal currents \f$I_{pol},\:I_{tor}\f$ and Fourier coefficients 
\f$R_{mn},\:\Phi_{mn},\:Z_{mn},\:B_{mn}\f$ are stored in a \ref W7X_Format "Boozer file".
This information is enough to calculate the Jacobian and magnetic field using formulas

\f{eqnarray*}
	J(\psi,\theta,\varphi) &=& \frac{\mu_{0}}{4\pi^{2}}\frac{I_{pol}+\iota I_{tor}}{B^{2}}\\
	\mathbf{B}(\psi,\theta,\varphi) &=& \frac{1}{2\pi J}\left(\frac{\partial\mathbf{X}}{\partial\varphi}+
	\iota\frac{\partial\mathbf{X}}{\partial\theta}\right)
\f}

where \b X is the spatial position given by

\f[
\mathbf{X}=r(\psi,\theta,\varphi)\cos(\phi)\mathbf{\hat{x}}+
           r(\psi,\theta,\varphi)\sin(\phi)\mathbf{\hat{y}}+
           z(\psi,\theta,\varphi)\hat{\mathbf{z}}
\f]

MConf library uses normalized toroidal flux \f$s=\psi/\psi_{lcms}\f$ as a 
flux surface label, which is zero on a magnetic axes and one at the LCMS. 
The flux surface label \a s is the same as in original VMEC file.
Then expressions for the Jacobian and magnetic field are as follows 

\f{eqnarray*}
	J_{s}(s,\theta,\varphi) &=& \frac{\mu_{0} \psi_{lcms}}{4\pi^{2}}\frac{I_{pol}+\iota I_{tor}}{B^{2}}\\
	\mathbf{B}(s,\theta,\varphi) &=& \frac{\psi_{lcms}}{2\pi J_{s}}\left(\frac{\partial\mathbf{X}}{\partial\varphi}+
	\iota\frac{\partial\mathbf{X}}{\partial\theta}\right)
\f}

Boozer coordinate system represent the magnetic coordinates, so the B-field line going through the point \f$(s_{0},\theta_{0},\varphi_{0})\f$ is defined by linear expression

\f[
\theta(\varphi)=\theta_{0}+\iota(s_{0})(\varphi-\varphi_{0})
\f]

There are two library methods that can be applied to follow field line:

- MConf::C3dMesh::magCoordFieldLine()  uses Boozer toroidal angle,
- MConf::C3dMesh::mixCoordFieldLine()  uses cylindrical angle.

<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

//************************************************************************************
\subsection  TokamakSymmetryFluxCoordinates Tokamak Symmetry Flux Coordinates

MConf can import equilibrium from \ref GEQDSK "G EQDSK File" into straight-field-line coordinates system
with the symmetry angle \f$\phi\f$ as the toroidal angle \f$\varphi\f$. 
This coordinate system is called Tokamak Symmetry Flux Coordinates in  \ref References "[1]".  
Other name of this system is the PEST coordinates \ref References "[7]". 


\f[
\mathbf{B}=\frac{1}{2\pi}\left(\nabla\psi\times\nabla\theta+\iota(\psi)\nabla\varphi\times\nabla\psi\right)
\f]
or in the form used in the tokamak literature
\f[
\mathbf{B}=\mathbf{B_{\mathit{tor}}}+\mathbf{B_{\mathrm{\mathit{pol}}}}=
\frac{1}{2\pi}\left(\mu_{0}I_{pol}(\psi)\nabla\varphi+
\iota(\psi)\nabla\varphi\times\nabla\psi\right)
\f]

In MConf, the magnetic coordinates 
\f$(s,\theta,\varphi)\f$ are defined as a Fourier decomposition by giving 
the cylindrical coordinates as functions of magnetic coordinates
 
\f{eqnarray*}
  r(s,\theta,\varphi) &=& \sum_{m=0}^M\sum_{n=0}^{1}R_{mn}(s)\cos\left(m\theta-n\frac{\pi}{2}\right)\\
  \phi(s,\theta,\varphi) &=& \varphi\\
	z(s,\theta,\varphi) &=& \sum_{m=0}^{M}\sum_{n=0}^{1}Z_{mn}(s)\sin\left(m\theta-n\frac{\pi}{2}\right)\\
\f}

With this representation (summation over \em n from 0 to 1) the tokamak 
equilibrium from \ref GEQDSK "G EQDSK File" can be saved in \ref W7X_Format "W7-X format", where 
\f$n \pi/2\f$ term provides sin and cos terms in expansions. 
The coefficients \f$R_{mn},\:Z_{mn}\f$ are calculated by B-field 
line tracing while creating straight-field-line coordinates system. 

Fourier decomposition of the tokamak equilibrium allows us 
to use all MConf-library functions without need to be rewritten;  
the Jacobian is only different. The import is provided by class MConf::CEfit.

The Jacobian and magnetic field are calculated using formulas

\f{eqnarray*}
	J(s,\theta,\varphi) &=& \frac{\psi_{lcms} r^2}{\mu_{0}I_{pol}}\\
	\mathbf{B}(s,\theta,\varphi) &=& \frac{\psi_{lcms}}{2\pi J}\left(\frac{\partial\mathbf{X}}{\partial\varphi}+
	\iota\frac{\partial\mathbf{X}}{\partial\theta} \right)
\f}

Poloidal current \f$I_{pol}=2\pi rB_{tor}/\mu_{0}\f$ is provided by G EQDSK File.
  
B-field line \f$(s_{0},\theta,\varphi)\f$ going through the point 
\f$(s_{0},\theta_{0},\varphi_{0})\f$ is defined by expression

\f[
\theta(\varphi)=\theta_{0}+\iota(s_{0})(\varphi-\varphi_{0})
\f]

Both methods can be applied to follow field line:

- MConf::C3dMesh::magCoordFieldLine()  uses cylindrical angle,
- MConf::C3dMesh::mixCoordFieldLine()  uses cylindrical angle.

\note

To distinguish this coordinates from "true" \ref W7X_Format 
"W7-X format" the following comment string
 
\code CC Tokamak Symmetry Flux Coordinates \endcode

is written into files exported by 
MConf::CStconfig::write,  MConf::CStconfig::writeasciiReduced
  
<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

//************************************************************************************
\section  CoordinateTransformation Coordinate Transformation

\subsection  VMEC2Boozer Transformation from VMEC to Boozer Coordinates

Substituting the Boozer angles \f$(\theta_{B}, \varphi_{B})\f$  expressed through 
the VMEC angles \f$(\theta_{V}, \varphi_{V})\f$ 

\f{eqnarray*}
\theta_{B} & = & \theta_{V}+\lambda(s,\theta_{V},\varphi_{V})+\iota h(s,\theta_{V},\varphi_{V})\\
\varphi_{B} & = & \varphi_{V}+h(s,\theta_{V},\varphi_{V})
\f}
 
into the covariant representation of the magnetic field in \ref BoozerCoordinates we obtain
the equation for determining the double-periodic transformation function \f$h\f$

\f[
\mathbf{B}=\frac{\mu_{0}}{2\pi}\left[I_{pol}\nabla\varphi_{V}+I_{tor}\left(\nabla\theta_{V}+\nabla\lambda\right)+
\left(I_{pol}+\iota I_{tor}\right)\nabla h+\left(\beta_{*}+I_{tor}h\frac{d\iota}{d\psi}\right)\nabla\psi\right]
\f]

The above expression yields the following covariant components of the magnetic field in VMEC coordinates:  

\f[
B_{\theta_{V}}=\frac{\mu_{0}}{2\pi}\left[I_{tor}+I_{tor}\frac{\partial\lambda}{\partial\theta_{V}}+\left(I_{pol}+
\iota I_{tor}\right)\frac{\partial h}{\partial\theta_{V}}\right]
\f]

\f[
B_{\varphi_{V}}=\frac{\mu_{0}}{2\pi}\left[I_{pol}+I_{tor}\frac{\partial\lambda}{\partial\varphi_{V}}+
\left(I_{pol}+\iota I_{tor}\right)\frac{\partial h}{\partial\varphi_{V}}\right]
\f]

from which the Fourier coefficients \f$h_{mn}\f$  of the angle transformation function \f$h\f$ 
are expressed through the known VMEC quantities:

\f{eqnarray*}
h_{00} & = & 0\\
I_{tor}\lambda_{mn}+\left(I_{pol}+\iota I_{tor}\right)h_{mn} & = & \frac{2\pi}{\mu_{0}} \frac{1}{m} \left(B_{\theta_{V}}\right)_{mn}\,\,\, if\, m\neq0\\
 & = & -\frac{2\pi}{\mu_{0}} \frac{1}{nN_p} \left(B_{\varphi_{V}}\right)_{mn}\,\,\, if\, n\neq0
\f}

The Boozer spectra of quantities of interest are calculated by integrating over the Boozer angles

\f{eqnarray*}
    B_{mn}^{B} & = & \frac{1}{2\pi^2}\int^{2\pi}_0 \int^{2\pi}_0
    d\theta_B d\varphi_B\, cos(m\theta_B - nN_p\varphi_B)\,B(\theta_B,\varphi_B)\\
    R_{mn}^{B} & = & \frac{1}{2\pi^2}\int^{2\pi}_0 \int^{2\pi}_0
    d\theta_B d\varphi_B\, cos(m\theta_B - nN_p\varphi_B)\,r(\theta_B,\varphi_B)\\
    Z_{mn}^{B} & = & \frac{1}{2\pi^2}\int^{2\pi}_0 \int^{2\pi}_0
    d\theta_B d\varphi_B\,sin(m\theta_B - nN_p\varphi_B)\,z(\theta_B,\varphi_B)\\
    \frac{2\pi}{nN_p}\Phi_{mn}^{B} & = & \frac{1}{2\pi^2}\int^{2\pi}_0 \int^{2\pi}_0
    d\theta_B d\varphi_B\,sin(m\theta_B - nN_p\varphi_B)\,h(\theta_B,\varphi_B)\\
\f}

 or VMEC angles

\f{eqnarray*}
    B_{mn}^{B} & = & \frac{1}{2\pi^2}\int^{2\pi}_0 \int^{2\pi}_0
    d\theta_V d\varphi_V\, \frac{\partial(\theta_B,\varphi_B)}{\partial(\theta_V,\varphi_V)}\,
    cos(m\theta_B - nN_p\varphi_B)\,B(\theta_V,\varphi_V)\\
    R_{mn}^{B} & = & \frac{1}{2\pi^2}\int^{2\pi}_0 \int^{2\pi}_0
    d\theta_V d\varphi_V\, \frac{\partial(\theta_B,\varphi_B)}{\partial(\theta_V,\varphi_V)}\,
    cos(m\theta_B - nN_p\varphi_B)\,r(\theta_V,\varphi_V)\\
    Z_{mn}^{B} & = & \frac{1}{2\pi^2}\int^{2\pi}_0 \int^{2\pi}_0
    d\theta_V d\varphi_V\,\frac{\partial(\theta_B,\varphi_B)}{\partial(\theta_V,\varphi_V)}\,
    sin(m\theta_B - nN_p\varphi_B)\,z(\theta_V,\varphi_V)\\
    \frac{2\pi}{nN_p}\Phi_{mn}^{B} & = & \frac{1}{2\pi^2}\int^{2\pi}_0 \int^{2\pi}_0
    d\theta_V d\varphi_V\,\frac{\partial(\theta_B,\varphi_B)}{\partial(\theta_V,\varphi_V)}\,
    sin(m\theta_B - nN_p\varphi_B)\,h(\theta_V,\varphi_V)\\
\f}


where \f$m=0..M,\, n=-N..N,\f$ \f$R_{00}^B=R_{mn}^B/2,\f$  \f$B_{00}^B=B_{mn}^B/2,\f$
 and 

\f[
\frac{\partial(\theta_B,\varphi_B)}{\partial(\theta_V,\varphi_V)} =
\left(1+\frac{\partial\lambda}{\partial\theta_{V}}  \right) 
\left(1+\frac{\partial h}{\partial\varphi_{V}}\right)+
\frac{\partial h}{\partial\theta_{V}}
\left(\iota-\frac{\partial\lambda}{\partial\varphi_{V}}\right)
\f]

<br>
The following methods perform transformation between the Boozer and VMEC coordinates:<br><br>
MConf::CStconfig::Vmec2Boozer() <br>
MConf::CStconfig::Boozer2Vmec() <br>
MConf::CStconfig::writeVmec2Boozer() <br>

The example of how to transform the VMEC wout-file to
Boozer-coordinate data file in \ref W7X_Format "W7X format" is as follows
  \code
  #include "CStconfig.h"
  int main() {
    MConf::CStconfig mConf;
    if(!mConf.load("wout.w7x-sc1.txt")) exit(1);
    mConf.writeVmec2Boozer("w7x-sc1.bc");
  }
  \endcode

see also \ref Gettingstarted and \ref MinimalExample   

<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>


\subsection  Cyl2Mag Transformation from cylindrical to flux Coordinates

The formulas in section \ref FluxCoordinates give prescription of how to find
spatial position and magnetic field for point given in flux coordinates.
However in most cases one need to know the flux surface label and
the magnetic field at point given in cylindrical coordinates. Let
us consider some examle. In pencil-beam approach the attenuation of
the neutral beam injected into plasma can be described by the following
equations 

\f{eqnarray*}
\frac{dI(l)}{dl} & = & -N_{e}(\mathbf{X})\sigma_{eff}I\\
\frac{d\mathbf{X}}{dl} & = & \mathbf{\hat{v}}_{b}
\f}

where \f$I\f$ is the beam current, \f$l\f$ is the length along beam, \f$N_{e}\f$is
the plasma density, \f$\sigma_{eff}\f$is the effective cross section,
\f$\mathbf{X}\f$is the spatial position, \f$\mathbf{\hat{v}}_{b}\f$is the
unit vector of beam direction. The second equation gives the beam
trajectory in real space, though the plasma density or plasma profile
\f$N_{e}(s)\f$ is given as a function of flux surface label \f$s\f$ or effective
radius \f$r_{eff}=a\sqrt{s}\f$ ( \a a is the minor plasma radius).
The coordinates transformation \f$s(\mathbf{X})\f$ from real
space coordinates to the flux coordinates is needed in order to
find plasma parameters along beam trajectory. This problem is common
for plasma diagnostics, for pellet injection modeling, ECRH ray/beam
tracing, NBI heating modeling.

\subsubsection  NewtonMethod Newton method

To do the transformation from cylindrical coordinates \f$(r_{0},\phi_{0},z_{0})\f$
to the flux coordinates system \f$\mathbf{u}=(u^{1},u^{2},u^{3})=(s,\theta,\varphi)\f$
the following system has to be solved
 
\f[
0=\mathbf{F}(\mathbf{u})\equiv
\left[
\begin{array}{r}
  \sum R_{mn}(s)\cos\left(m\theta-nN_{p}\varphi\right)\;\;-r_{0}\\
  \varphi-\frac{2\pi}{N_{p}}\sum\Phi_{mn}(s)\sin\left(m\theta-nN_{p}\varphi\right)\;-\phi_{0}\\
  \sum Z_{mn}(s)\sin\left(m\theta-nN_{p}\varphi\right)\;\;-z_{0}\end{array}
\right]
\f]

Using initial guess \f$\mathbf{u}_{k}\f$, the correction \f$\Delta\mathbf{u}_{k}\f$
is calculated from the system obtained by Taylor expansion

\f[
\frac{\partial\mathbf{F}(\mathbf{u}_{k})}{\partial\mathbf{u}_{k}}
\Delta\mathbf{u}_{k}=-\mathbf{F}(\mathbf{u}_{k})
\f]

The solution is then found using iterations 
\f[\mathbf{u}_{k+1}=\mathbf{u}_{k}+\gamma\!\Delta\mathbf{u}_{k}\f]
where 

\f[
\gamma=\Biggl\{\begin{array}{ll}
0.9\, s_{k}/\Delta s_{k} & \textrm{if}\,\,(s_{k}+\Delta s_{k})<0\\
1 & \textrm{otherwise}\end{array}
\f]

is the scale factor used to avoid negative flux label \a s.


This method is fast and accurate. Basis vectors \f$\partial\mathbf{X}/\partial u^{i}\f$
and thus other quantities are available right after transformation.
With a "good" guess only two or three iterations are needed to achieve accuracy 
of 0.1mm for W7-X magnetic configurations. The information from spatial points 
used in a previous coordinate transformations provides the good guess, for example 
during generating 3d-mesh or ray tracing. In cases when there is no information 
from nearest points available the library MConf uses rather slow root finder to 
compute the guess. The finder employs 2d-iterations on R-Z-plane.  

The figure below illustrates efficiency of the Newton method even in cases of "bad" 
guess. The iterations going from the guesses to solutions are shown.

<table border="0">
<tr>
	<td> \image html Newtonmethod02.png " "  </td>
</tr>
</table><br>

<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>


//************************************************************************************
\section  Gettingstarted Getting started

This documentation is build with Doxygen which generates 
a lot of information that is useful for developer, but makes usual user confused.  
Naturally the question arises how to start to use MConf and which classes are really needed. 

In most cases you need only MConf::CStconfig or MConf::C3dMesh

- MConf::CStconfig is the main class which implements geometric calculation. 
   In particular, this class provides forth and back coordinate transformation 
   between flux and real space coordinates. 
- MConf::C3dMesh is the class derived from CStconfig and all methods from CStconfig are available. 
   In addition, C3dMesh tabulates magnetic field, flux surface label, grad(s) on a 3d-mesh
   in cylindrical coordinates and provides functions for interpolation. Use this class if you need fast 
   coordinate mapping from real space to flux surface label. You can even do magnetic field line tracing.

You may also find the class MConf::CProfile to be usefull.   
   
  See also some <b><a class="el" href="examples.html">examples</a></b>.


<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

//************************************************************************************
\section  MinimalExample Minimal example

This is a minimal example of how to use the CStconfig class.

The steps are:

\code
#include "CStconfig.h"                          //1. include the header file
#include <iostream>

int main(int argc, char* argv[])
{
  using namespace MConf;                        //2. MConf functions are in MConf namespace
  
  const double degree = 3.1415926535898/180;
  CStconfig mc;                                 //3. Declare object

  if(false==mc.load("w7x-sc1.bc")) {            //4. load into it magnetic configuration
    std::cout << "Loading error"<<std::endl;
    return 1;   // exit if errors
  }
  //set point in cylindrical coord.
  Vector3d cyl(6,1*degree,0.5);                 //5. set point in real space
  // s is the normalized toroidal flux
  double s = mc.cyl2s(cyl);                     //6. find flux label that corresponds to cyl
  std::cout<<"cyl="<<cyl<<"  s="<<s<<std::endl; //7. print results
  return 0;
}
\endcode

Functions of MConf-library are defined in namespace  MConf, that is why 
<em>using namespace MConf</em> appears in step 2. 
However I would recommend not to use <em>using namespace</em> 
in your program. Otherwise you'll have problems in the future 
when someone (or even you) adds a symbol in other namespace 
and it will collide with one in yours. Using an individual class or 
function name (like <em>using Mconf::Vector3d;</em>) is okay, 
but don't do it at the file scope or a header file.

See also other <b><a class="el" href="examples.html">examples</a></b>.

<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

//************************************************************************************
\section  ObjectCanBeCopied Objects can be copied

MConf::CStconfig, MConf::CRayTrace, and MConf::C3dMesh  objects internally use data
sharing, so they can be cloned/copied without CPU time and memory overhead since 
only the pointers are copied (they are really smart).
Writing to the clone object may take some time because the data referenced 
by the corresponding pointer is automatically cloned to preserve original instance of the object.
The  method is called "copy-on-write." With copy-on-write, two or more objects can share the same data until 
the moment when one of those objects is changed, at which point the data is physically copied 
and changed in one of the objects.
This allows effective using of the MCONF in the multithreading programming.

see also MConf::pBase::ngArray description. 

<!--
A technique which blends the advantages of mutable and immutable objects, 
and is supported directly in almost all modern hardware, is copy-on-write (COW). 
Using this technique, when a user asks the system to copy an object, 
it will instead merely create a new reference which still points to the same object. 
As soon as a user modifies the object through a particular reference, 
the system makes a real copy and sets the reference to refer to the new copy. 
The other users are unaffected, because they still refer to the original object. 
Therefore, under COW, all users appear to have a mutable version of their objects, 
although in the case that users do not modify their objects, 
the space-saving and speed advantages of immutable objects are preserved. 
-->


\code
#include "C3dMesh.h"

int main(int argc, char* argv[])
{
  MConf::C3dMesh mc1, mc2, mc3; // C3dMesh is defined in MConf namespace
                                // the class is prefixed by namespace
  m1.load("filenameOfMesh.bin4");

  if(mc1.isMeshOK()==false) {  // create mesh if needed
    mc1.createMeshUsingSymmetry (0.02,0.02,2*degree);
    if(false==mc1.isMeshOK()) exit(1);    // exit if errors
  }

  mc2 = mc1;                  // This assignment is completely legal
  mc3 = mc2;

  MConf::Vector3d c1(5),c2(5.1),c3(5.4);  

  double s1 = mc1.cyl2s(c1);  // can be run in new execution thread, not shown here
  double s2 = mc2.cyl2s(c2);  // can be run in new execution thread
  double s3 = mc3.cyl2s(c3);  // can be run in new execution thread

  return 0;
}

\endcode

The advantage of the code above is that the mesh is generated only once and
then it can be used in different execution threads on a multi-processor systems
with shared memory. 

This technique is used in MConf::C3dMesh::createMeshUsingSymmetry(), MConf::CStconfig::epsEffCreate(), 
MConf::CStconfig::FbsCreate(), in flux-surface averaging procedures and 
in NBI-code to run several NBI beams through plasma concurrently, see \ref
RelatedProjects and \ref threads.h

<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

//************************************************************************************
\section  RelatedProjects Related Projects


<!--
\li <b><a class="el" href="../../../../../../MCviewer/development/implementation/docs/userManual/index.html">
  Magnetic Configuration viewer</a></b>
-->

\li <b><a class="el" href="../docs-mcviewer/index.html">Magnetic Configuration viewer</a></b> 
  displays a magnetic configuration and is written on a base of the 
  MConf library. The MCviewer (\ref MCviewer "sreenshoots") is written in C++,
  using OpenGL library and the FLTK project (a cross-platform GUI
  library http://www.fltk.org/ ). The FLTK library are distributed by its
  authors under the terms of the GNU Library General Public License (LGPL)
  with the exceptions, which do not require the author to provide source 
  code of MCviewer. 

\li \ref Travis (stands for TRAce and VISualize) is the multi-beam multi-pass
  ray-tracing code for electron cyclotron heating and current drive calculations;
  it uses MConf library.

\li \ref NBIviewer is the code for Neutral Beam Injection heating simulation with
  the graphical user interface (GUI) front-end to the
  NBI-code; it uses MConf library.

 <h2>Acknowledgments</h2>
 
The 
<a class="el" href="mailto:yuriy.turkinNOSPAMipp.mpg.de?body=Bitte NOSPAM in der Email-Adresse durch at-Zeichen ersetzen. 
- Please replace NOSPAM in the e-mail address by the at symbol.">author</a>
thanks
H. Maa&szlig;berg, J. Geiger, C.D. Beidler, A. Dinklage, H. Dreier, Y. Feng,
N. Marushchenko, M. Schmidt, J. Svensson, A. Werner
for suggestions, contributions, helpful comments, and bug reports.

\htmlonly
<div><table border="0">
<tr> 
<td>Special thanks go to Dimitri van Heesch for his program&nbsp;&nbsp;<br> 
that has been used for producing this documentation.</td> 
<td><a href="http://www.doxygen.org/index.html">
<img src="doxygen.png" alt="doxygen" align="right" border="0"></a></td> 
<td></td>
</tr>
</table></div>
\endhtmlonly

\anchor  BottomOfPage

*/

//***************************************************************************
//******other Pages**********************************************************
//***************************************************************************
//***************************************************************************
/** @page Screenshots Screenshots

\section  MCviewer

<table border="0">
<tr>
	<td> \image html mcviewer001.png " "  </td>
	<td> \image html mcviewer002.png " "  </td>
</tr>
</table><div class="no-print"><br><hr></div>

\section  Travis

<table border="0">
<tr>
	<td> \image html travis001.png " "  </td>
	<td> \image html travis002.png " "  </td>
</tr>
</table><div class="no-print"><br><hr></div>

\section  NBIviewer

<table border="0">
<tr>
	<td> \image html nbiviewer001.png " "  </td>
	<td> \image html nbiviewer002.png " "  </td>
</tr>
</table><br>

*/

//************************************************************************************
/** @page W7X_Format W7-X Format

\par Format of bc-file.

\verbatim
CC Boozer-coordinate data file
CC Version 0.1
CC Author J.Geiger
CC Created: 30.04.2004
 m0b  n0b nsurf nper flux/[Tm^2]     a/[m]     R/[m]
   36   36   97    5 -2.418619E+00   0.51092   5.52666
       s         iota  curr_pol/nper    curr_tor    pprime   sqrt g(0,0)
                            [A]            [A]   dp/ds,[Pa] (dV/ds)/nper
  1.5306E-02  8.5653E-01 -1.7885E+07 -1.2735E+00 -1.2474E+03 -5.7549E+00
    m    n      R_mn/[m]      Z_mn/[m]          F_mn           bmn/[T]
    0   36  1.66348028E-13 -6.08412908E-14  1.03164537E-13 -1.57826503E-11
    0   35  6.57460332E-14 -1.07438400E-13 -1.59783659E-14 -9.03979668E-12
    0   34 -1.88635922E-12  1.52797654E-12  1.68670282E-13  7.55025114E-12
 ...............   
\endverbatim

\par How to read bc-file, initial lines of which are shown above.

\include w7xFormat.f90

*/

//************************************************************************************
/** @page MCDB Magnetic Configuration Database

\par

The schema of Magnetic Configuration Database (MCDB) is implemented in a
relation database (ORACLE). Magnetic equilibria for the plasma core are
stored with their input parameters (e.g. plasma profiles) and code execution
parameters.

\sa

Dinklage, A., J. Geiger, G. Kuehner, M. Schmidt, Yu.A. Turkin and A. Werner: 
A magnetic configuration database for Wendelstein 7-X and its application 
programming interface. In: Fusion Engineering and Design 81, 1969-1973(2006).
http://dx.doi.org/10.1016/j.fusengdes.2006.04.047

<div><table border="0">
<tr> 
<td width="100"><a href="../MCDB/index.html">MCDB</a></td> 
<td> \image html undercon.gif " "  </td> </tr>
</table></div>

*/

//************************************************************************************
/** @page GEQDSK EFIT G EQDSK File

\par

This file describes tokamak equilibrium and contains information on P', FF', 
the poloidal flux psi \f$\psi\f$ on the rectangular grid used, q, the boundary, 
and the limiter contour. A right-handed cylindrical coordinate 
system \f$(r,\phi,z)\f$ is used. The G EQDSK provides information on the pressure, 
poloidal current function, q profile on a uniform flux grid from the magnetic axis 
to the plasma boundary and the poloidal flux function on the rectangular computation 
grid. Information on the plasma boundary and the surrounding limiter contour is also 
provided. http://fusion.gat.com/efit/g_eqdsk.html

\par  Format and  Variables 

\verbatim 
  character*10 case(6)
  dimension psirz(nw,nh),fpol(1),pres(1),ffprim(1), &
   pprime(1),qpsi(1),rbbbs(1),zbbbs(1),rlim(1),zlim(1)   
  read (neqdsk,'(6a8,3i4)') (case(i),i=1,6),idum,nw,nh
  read (neqdsk,'(5e16.9)') rdim,zdim,rcentr,rleft,zmid
  read (neqdsk,'(5e16.9)') rmaxis,zmaxis,simag,sibry,bcentr
  read (neqdsk,'(5e16.9)') current,simag,xdum,rmaxis,xdum
  read (neqdsk,'(5e16.9)') zmaxis,xdum,sibry,xdum,xdum
  read (neqdsk,'(5e16.9)') (fpol(i),i=1,nw)
  read (neqdsk,'(5e16.9)') (pres(i),i=1,nw)
  read (neqdsk,'(5e16.9)') (ffprim(i),i=1,nw)
  read (neqdsk,'(5e16.9)') (pprime(i),i=1,nw)
  read (neqdsk,'(5e16.9)') ((psirz(i,j),i=1,nw),j=1,nh)
  read (neqdsk,'(5e16.9)') (qpsi(i),i=1,nw)
  read (neqdsk,'(2i5)'   ) nbbbs,limitr
  read (neqdsk,'(5e16.9)') (rbbbs(i),zbbbs(i),i=1,nbbbs)
  read (neqdsk,'(5e16.9)') (rlim(i),zlim(i),i=1,limitr)
\endverbatim 

\arg \c case &nbsp;&nbsp;Identification character string
\arg \c nw &nbsp;&nbsp;Number of horizontal R grid points
\arg \c nh &nbsp;&nbsp;Number of vertical Z grid points
\arg \c rdim &nbsp;&nbsp;Horizontal dimension in meter of computational box
\arg \c zdim &nbsp;&nbsp;Vertical dimension in meter of computational box
\arg \c rleft &nbsp;&nbsp;Minimum R in meter of rectangular computational box
\arg \c zmid &nbsp;&nbsp;Z of center of computational box in meter
\arg \c rmaxis &nbsp;&nbsp;R of magnetic axis in meter
\arg \c zmaxis &nbsp;&nbsp;Z of magnetic axis in meter
\arg \c simag &nbsp;&nbsp;poloidal flux at magnetic axis in Weber/rad
\arg \c sibry &nbsp;&nbsp;poloidal flux at the plasma boundary in Weber/rad
\arg \c rcentr &nbsp;&nbsp;R in meter of vacuum toroidal magnetic field \c bcentr
\arg \c bcentr &nbsp;&nbsp;Vacuum toroidal magnetic field in Tesla at \c rcentr
\arg \c current &nbsp;&nbsp;Plasma current in Ampere
\arg \c fpol &nbsp;&nbsp;Poloidal current function in m T, F = RB_t on flux grid
\arg \c pres &nbsp;&nbsp;Plasma pressure in nt/m^2 on uniform flux grid
\arg \c ffprim &nbsp;&nbsp;FF'(psi) in (m T)^2/(Weber/rad) on uniform flux grid
\arg \c pprime &nbsp;&nbsp;P'(psi) in (nt/m^2)/(Weber/rad) on uniform flux grid
\arg \c psizr &nbsp;&nbsp;Poloidal flux in Weber/rad on the rectangular grid points
\arg \c qpsi &nbsp;&nbsp;q values on uniform flux grid from axis to boundary
\arg \c nbbbs &nbsp;&nbsp;Number of boundary points
\arg \c limitr &nbsp;&nbsp;Number of limiter points
\arg \c rbbbs &nbsp;&nbsp;R of boundary points in meter
\arg \c zbbbs &nbsp;&nbsp;Z of boundary points in meter
\arg \c rlim &nbsp;&nbsp;R of surrounding limiter contour in meter
\arg \c zlim &nbsp;&nbsp;Z of surrounding limiter contour in meter

 
The toroidal current J_t related to P'(psi) and FF'(psi) through

  J_t (Amp/m2) = R P'(psi) + FF'(psi) / R

*/

//************************************************************************************
/** @page MixedLanguage Mixed-language programming

\par

- \ref FortranCallsC 
- \ref CCallsFortran77
- \ref CCallsFortran90
- \ref Alignment  

The interface described below removes decoration of names and therefore simplifies
mixed-language programming. A decorated name is an internal representation
of a procedure name or variable name that contains information about where
it is declared; for procedures, the information includes how it is called.

Some examples. The name of the Fortran routine \c LOADVMEC with declaration
\code SUBROUTINE LOADVMEC(fname)\endcode 
in an object file becomes
\code loadvmec__ \endcode

The name of C++ method \c loadvmec with declaration
\code bool CStconfig::loadvmec(const char * fname) \endcode 
in an object file becomes
\code _ZN9CStconfig8loadvmecEPKc \endcode

Fortran can not call directly C++ methods also because of Fortran
knows nothing about C++ objects and the C++ object which the 
method belongs to must be created first.

\section FortranCallsC  Fortran calls C++

This is an example of how to write interface C-function in 
order to call C++ function from within Fortran programs.

Suppose we have C++ class \c Box declared in the header 
file \c \#include \c "Box.h"

\code 
class Box {
public:
  // Default constructor
  Box() { lb=bb=hb=0; }
  // Constructor 
  Box(double length, double breadth, double height) {
    set(length, breadth, height);
  }
  // Function to set dimensions of a box
  void set(double length, double breadth, double height) {
     lb = length;  // Set values of data members
     bb = breadth; 
     hb = height;
  }
  // Function to calculate the volume of a box
  double getVolume() {
     return lb*bb*hb;
  }
private:
  double lb; // Length of a box in inches
  double bb; // Breadth of a box in inches
  double hb; // Height of a box in inches
};
\endcode
 
and we want to call \c Box::getVolume(). The solution is to write C 
wrappers (a C-function that provides an interface to C++ function 
or object) as it follows 
 
\code 
#include "Box.h"

// the following preprocessor directives take care about
// different naming conventions in Windows and *nix world
#if (defined( __unix__ ) || (defined(linux)) )
  #define __stdcall
  #define BOX_LOAD    box_load_
  #define BOX_VOLUME  box_volume_
  #define BOX_FREE    box_free_
#elif defined( _WIN32 )
#endif


// Create and initialize an object Box. 
// The first call must be made with *box==0, 
// the following calls will reload data into the same object
extern "C" void __stdcall BOX_LOAD(Box **box,double *l,double *b,double *h)
{
  Box *bx;
  if(*box==0) bx = new Box;   // create object Box, 
  bx->set(*l, *b, *h);        //  initialize it
  *box = bx;                  //   and return to caller its address
}

// Function to get the volume of a box
extern "C" double __stdcall BOX_VOLUME(Box **box)
{
  Box &bx = **box;       // set reference to box
  return bx.getVolume(); // call object's method
}

// 'destructor' for Fortran
extern "C" void __stdcall BOX_FREE(Box **box)
{
  if (*box) delete (*box);
  *box=0;
}
\endcode

Finally the Fortran routine is

\code
  implicit none
  real(8) BOX_VOLUME, volume
  real(8) length/20./, breadth/10./, height/10./
  integer(8) toyBox/0/  ! this variable holds the address of C++ object
  external BOX_LOAD, BOX_VOLUME, BOX_FREE

! create  toyBox object of type Box and initialize it
  call BOX_LOAD(toyBox,length, breadth, height)
  
! get the volume of the toyBox
  volume = BOX_VOLUME(toyBox)
  
! destroy toyBox object and free memory
  call BOX_FREE(toyBox)       
  end
\endcode

\c toyBox is the pointer to the object of type Box; 
8-byte is used in order to compile and run this program on 
a CPU with 64bit architecture without modifications.
  
The interfacing described above is employed in cpp2for.cpp

<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

\section CCallsFortran77  C++ calls Fortran 77

This is an example of how to call Fortran subroutine from within C++ programs.

Suppose we have the Fortran 77 subroutine 

\code
      subroutine box_vol_f77(al,ab,ah,volume,arr,n)
\endcode

that we need to call from C++. The code below does this

\code
//**************************************************************
#if (defined( __unix__ ) || (defined(linux)) )
  #define SUBROUTINE        extern "C" void 
  #define REAL8_FUNCTION    extern "C" double 
  #define INT_FUNCTION      extern "C" int 
  
  #define BOX_VOL_f77  box_vol_f77_
#elif defined( _WIN32 )
  #define SUBROUTINE        extern "C" void   __stdcall
  #define REAL8_FUNCTION    extern "C" double __stdcall
  #define INT_FUNCTION      extern "C" int    __stdcall
#endif

#define REAL8 double
#define CALL    ;    // to mimic FORTRAN CALL

// prototype of Fortran routine
SUBROUTINE BOX_VOL_f77(
       const  REAL8 &l,       //[input]  scalar, [cm]
       const  REAL8 &b,       //[input]  scalar, [cm]
       const  REAL8 &h,       //[input]  scalar, [cm]
       const  REAL8 &volume   //[out] scalar, volume on output[cm^3]
       const  REAL8 *array,   // array(narr), not used, example only
       const  int   &narr,    // scalar, size of array
       );

// C function, that calls Fortran 77
double volume(double l,double b,double h)
{ 
  int n = 20;
  double arr[n], volume;      
  CALL BOX_VOL_f77(l,b,h,volume,arr,n);
  return volume;
}
\endcode

<div class="no-print" align="right"><table border="0">
<tr> 
<td>[ <a class="el" href="#TopOfPage">top</a> ]</td> 
<td>[ <a class="el" href="#BottomOfPage">bottom</a> ]</td>
<td>[ <a class="el" href="javascript:history.back(-1);"><b>back</b></a> ]</td> 
</tr>
</table></div><div class="no-print"><hr></div>

\section CCallsFortran90 C++ calls Fortran 90

Suppose we have the Fortran 90 module with subroutines \c box_set , \c box_volume

\code
  MODULE BOX_mod
    IMPLICIT NONE
    REAL(8), PRIVATE :: al,ab,ah
    CONTAINS
    
    SUBROUTINE box_set(l,b,h)
      REAL(8) l,b,h
      al = l
      ab = b
      ah = h
    END SUBROUTINE box_set
    
    SUBROUTINE box_volume(volume)
      REAL(8) volume
      volume = al*ab*ah
    END SUBROUTINE box_volume
    
   END MODULE BOX_mod
\endcode
 
At first we need Fortran 77 interface to remove decoration of names 

\code
SUBROUTINE BOX_VOL_f77(al,ab,ah,volume)
  USE BOX_mod, ONLY : box_set,box_volume
  implicit none
  real(8) l,ab,ah,volume
  call box_set(al,ab,ah)
  call box_volume(volume)
END SUBROUTINE BOX_f77
\endcode

then see section \ref CCallsFortran77

\section Alignment  Struct Member Alignment  

If you pass structures across the languages then the use of Alignment 
option is important
\code 
/Zp8  ---  Visual C, Intel C and Fortran ,  Intel/Digital/Compaq Fortran
\endcode 

The Struct Member Alignment (/Zpn) option controls how the members of 
a structure are packed into memory and specifies the same packing for 
all structures in a module. When you specify this option, each structure 
member after the first is stored on either the size of the member type 
or n-byte boundaries (where n is 1, 2, 4, 8, or 16), whichever is smaller.


<!--
 to be continued
<div>
<table border="1">
<tr> <td width="100">under construction</td> <td> \image html undercon.gif " "  </td> </tr>
</table>
</div>
-->

*/

//************************************************************************************
/** @page MatlabInterface Matlab (or Python) Interface

\par

The MConf can be used from within Matlab environment (or from Python) by 
calling MConf functions from a shared library: 
from dll under Windows, or from  so-library under Linux.
The only thing one have to provide is the 
C interface to MConf classes because Matlab (Python) 
can not call C++ methods directly.  
The technique to accomplish this task is similar to that
described in section \ref FortranCallsC. 
At first a C++ object is created by a C function and the resulting 
address is returned to a caller. Then the caller passes this address 
along with the parameters to another C function to invoke required object's method. 

Below is a working example of how to write interface in C in 
order to call C++ functions from within Matlab (Python) environment.

Header file; lines taken from mconf_matlab.h
\code
#ifndef __CPP2MATLAB_
#define __CPP2MATLAB_

#if (defined( __unix__ ) || (defined(linux)) || defined(__sun) )
  #define WINDLLEXPORT 
#elif defined( _WIN32 )
  #define WINDLLEXPORT  __declspec( dllexport )
#endif

#ifdef __GNUC__
 #define __int64 int64_t 
#endif

////class MConf::C3dMesh; typedef MConf::C3dMesh* MC_HANDLE;
////typedef int* MC_PTR;  
////typedef void* MC_PTR;

/// A platform-specific type that is used to represent a pointer.
//typedef __int64 MC_HANDLE;  // for 64-bit architecture
typedef unsigned int  MC_HANDLE;  // for 32-bit architecture

#ifdef __cplusplus
extern "C" {
#endif
  WINDLLEXPORT MC_HANDLE MCload      (char * fname);
  WINDLLEXPORT void   MCfree  (MC_HANDLE mConf);
  WINDLLEXPORT double MCgetB00(MC_HANDLE mConf);
  WINDLLEXPORT void   MCsetB00(MC_HANDLE mConf,double B00);
#ifdef __cplusplus
}
#endif

#endif  //  __CPP2MATLAB_

\endcode

Implemention; lines taken from  mconf_matlab.cpp

\code

#include "../include/C3dMesh.h"
#include <iostream>
#include "mconf_matlab.h"

using MConf::C3dMesh;
using MConf::CStconfig;
using MConf::Vector3d;

#undef True
#undef False

const static int True(1);
const static int False(0);

#ifdef __cplusplus
extern "C" {
#endif


//*********************************************************************
// constructor loads a file and returns the handle of the object 
// with the magnetic configuration
WINDLLEXPORT MC_HANDLE MCload(char * fname)
{
  C3dMesh * mc = new C3dMesh;  
  MC_HANDLE mConf = (MC_HANDLE)mc;
  bool ok = mc->load (fname);
  if(!ok) {
    std::cout << "MCLOAD: Loading error, file:'"<<fname<<"'\n";
    delete (mc);
    mConf = 0;    
  }
  return mConf;
}

//*********************************************************************
// destructor frees the memory occupied by the object with the handle mConf  
WINDLLEXPORT void MCfree(MC_HANDLE mConf)
{
  if (mConf) delete ((C3dMesh*)mConf);
}

//*********************************************************************
// function returns B_00 from the object with the handle mConf  
WINDLLEXPORT double MCgetB00(MC_HANDLE mConf)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  return mc.getB00();
}

//*********************************************************************
// function passes B_00 to the object with the handle mConf  
WINDLLEXPORT void MCsetB00(MC_HANDLE mConf, double B00)
{
  C3dMesh &mc = *((C3dMesh*)mConf);
  mc.setB00(B00);
}

#ifdef __cplusplus
}
#endif

\endcode

Finally the Matlab script is

\verbatim
loadlibrary('mconf_matlab.dll','mconf_matlab.h');
libfunctions mconf_matlab -full; 
       
fname='w7x-sc1(reduced).bc';
% load the magnetic configuration file
% @return -- if the function succeeds, the return value is 
% the handle of MConf object;  zero otherwise.
MC = calllib('mconf_matlab','MCload',fname); 

%  test the MC before next calls, it must be non-zero

B00 = calllib('mconf_matlab','MCgetB00',MC);
B00 = 2.3;
calllib('mconf_matlab','MCsetB00',MC,B00);

calllib('mconf_matlab','MCfree',MC);

unloadlibrary mconf_matlab

\endverbatim

Python script:

\verbatim

# This is working example of how to use the mconf.dll

# interfacing using numpy and ctypes.

import sys
import platform
import math
import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

os = platform.system()
is_64bits = sys.maxsize > 2**32

# load the library, using numpy mechanisms
if is_64bits:
  if os=='Windows':
    mconf = npct.load_library("mconf_matlab64.dll",".")
  elif os=='Linux':
    mconf = npct.load_library("mconf_matlab64.so",".")  
else:
  if os=='Windows':
    mconf = npct.load_library("mconf_matlab.dll",".")
  elif os=='Linux':
    mconf = npct.load_library("mconf_matlab.so",".")  

vec3 = npct.ndpointer(dtype=np.float64, ndim=1, flags='CONTIGUOUS')

# setup the return typs and argument types
mconf.MCload.restype = c_longlong  # for 64-bit architecture
mconf.MCload.argtypes = [c_char_p] 

mconf.MCgetRayIntersectionPoints.restype = c_int
mconf.MCgetRayIntersectionPoints.argtypes = [c_longlong,vec3,vec3,vec3,vec3]

mconf.MCgetB00.restype = c_double
mconf.MCgetB00.argtypes = [c_longlong]

mconf.MCsetB0.restype = c_double
mconf.MCsetB0.argtypes = [c_longlong,c_double,c_double]

mconf.MCgetBxyz.restype = c_double
mconf.MCgetBxyz.argtypes = [c_longlong,vec3,vec3]

mconf.MCVprime.restype = c_double
mconf.MCVprime.argtypes = [c_longlong,c_double]

mconf.MCVolume.restype = c_double
mconf.MCVolume.argtypes = [c_longlong,c_double]

mconf.MCtorFlux2polFlux.restype = c_double
mconf.MCtorFlux2polFlux.argtypes = [c_longlong,c_double]

fname=c_char_p("w7x-sc1beta=0.02.bc")
# load the magnetic configuration file
# @return -- if the function succeeds, the return value is 
# the address of C++ object;  zero otherwise.
mc = mconf.MCload(fname) # mc is like self in python
if mc == 0:
  print 'mconf: Could not load magnetic configuration'

B00 = mconf.MCgetB00(mc);
print B00 

# set magnetic field 2.5T at toroidal angle 0
##mconf.MCsetB0 (mc,c_double(2.5), c_double(0.));
B00 = mconf.MCgetB00(mc);
print B00 

# trace plasma along the ray through the port AEL41
# using Cartesian coordinates
r0 = np.array([-2.39133,-2.32718,-0.12071],dtype=np.float64)
r1 = np.array([-3.37847,-4.27681, 0.17038],dtype=np.float64)
rd = r1 - r0
rd = rd/np.sqrt(np.dot(rd,rd)) # normalize
print rd

# find entry point of the ray into plasma
entry= np.empty_like(r0) # ray entry
exit = np.empty_like(r0) # ray exit
retcode = mconf.MCgetRayIntersectionPoints(mc,r0,rd,entry,exit)
#  test retcode, it must be non-zero
if retcode == 0:
  print 'mconf: ray does not hit plasma'

print entry, exit, retcode

# ########################################################

maxPnt = 2000
dl = 0.01;     #  1cm
r0 = entry
          
# trace plasma along the ray
r = np.empty_like(r0)
rB = np.empty_like(r0)
for i in xrange(0,maxPnt):    
  lng  = i*dl          # length from r0 to r
  r = r0 + lng*rd      # move along the ray
  s = mconf.MCgetBxyz(mc,r,rB)
  if s>1:  break     # break if not inside plasma
  V =  mconf.MCVolume(mc,c_double(s))  # V  is the  volume inside the surface s 
  Vp = mconf.MCVprime(mc,c_double(s)) # Vp is the  dV/ds 
  #x = sqrt(s)       # x is the normalized plasma radius x=reff/a  
  #n = ne(x)         # density
  #t = Te(x)         # temperature
  sPol = mconf.MCtorFlux2polFlux(mc,c_double(s)) # sPol is the normalized poloidal flux, where s is the normalized toroidal flux.
  print lng, s, sPol, rB 

  
mconf.MCmix2xyz.argtypes = [c_longlong,vec3,vec3]
 
# ****************************************************
# plot flux surface s=0.5 at cyl. angle 2degree 
phi = 2*6.28318531/360    # 2 degree 
s   = 0.5         

m = np.empty_like(r)
maxPnt = 200
dth =6.28318531/(maxPnt-1)
for i in range(0,maxPnt):    
  th = i*dth  
  m[0] = s
  m[1] = th
  m[2] = phi
  
  mconf.MCmix2xyz(mc,m,r)
  R = math.sqrt(r[0]**2+r[1]**2)
  Z = r[2]
  print R, Z  # cyl. coordinates R,Z
\endverbatim
  
The interfacing described above is employed in mconf_matlab.cpp and mconf_matlab.h

*/


//************************************************************************************
/** @page References References

[1] W. D. D'haeseleer, W. N. G. Hitchon, W. I. van Rij, S. P. Hirshman, 
and J. L. Shohet, Flux Coordinates and Magnetic Field Structure 
(Springer-Verlag, New York, 1991). 

[2] Physics of Plasmas, Alen H.Boozer: What is a stellarator? 
Volume 5, 1647(1998) http://dx.doi.org/10.1063/1.872833

[3] A. Dinklage, J. Geiger, G. Kuehner, M. Schmidt, Yu.A. Turkin and A. Werner: 
A magnetic configuration database for Wendelstein 7-X and its application 
programming interface. In: Fusion Engineering and Design 81, 1969-1973(2006).
http://dx.doi.org/10.1016/j.fusengdes.2006.04.047

[4] S. P. Hirshman, W. I. van Rij and P. Merkel, Comput. Phys. Commun. \b 43, 143(1986).
<a href="http://dx.doi.org/10.1016/0010-4655(86)90058-5">http://dx.doi.org/10.1016/0010-4655(86)90058-5</a>

[5] N B. Marushchenko et. al. Ray Tracing Simulations of ECR Heating and ECE Diagnostic
at W7-X Stellarator. Plasma and Fusion Research, Vol.2, S1129(2007). DOI: 10.1585/pfr.2.S1129
http://dx.doi.org/10.1585/pfr.2.S1129

[6] Masanori NUNAMI, Tomo-Hiko WATANABE and Hideo SUGAMA. Plasma and Fusion Research
Vol. 5, 016(2010). DOI: 10.1585/pfr.5.016  http://dx.doi.org/10.1585/pfr.5.016

[7] R.C Grimm, R.L Dewar, J Manickam. Ideal MHD stability calculations in axisymmetric toroidal coordinate systems.
J. Comput. Phys. 49, 94 (1983). <a href="http://dx.doi.org/10.1016/0021-9991(83)90116-X">http://dx.doi.org/10.1016/0021-9991(83)90116-X</a>

<div class="no-print" align="center"><table border="0">
<tr> 
<td>\image html background.png</td>
<td>\image html background.png</td>
<td>\image html background.png</td>
</tr>
</table></div>

*/