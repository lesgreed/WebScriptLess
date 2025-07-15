#ifndef MC_VECTOR3D_H
#define MC_VECTOR3D_H

#include <stdio.h>
#include <math.h>
#include <iostream>

namespace MConf {

#ifdef _DEBUG
 #define  TESTZERO(x) if (!(x)) std::cerr<<"devide by zero in Vector3d operator/()"<<std::endl;
#else 
 #define  TESTZERO(x) 
#endif

//**********************************************************************
/// This class provides 3d-vector type.
/// Basically this is a cartesian vector.
/// In addition to the usual operation with vectors in Cartesian coordinates
/// the following operation are provided:
///  \li componentwise product,
///  \li componentwise division,
///  \li reading/writing from/to binary file,
///  \li seting/coping components from/to array.
///
/// This class can be used for storing the coordinates of a point or
/// vector in other coordinate system.
/// The transformations between cartesian and cylindrical
/// point or vector are provided;
/// Operations with vectors in cylindrical(or other) coordinates 
/// are not defined.
/// \par Usage:
/// \code
/// #include "CVector3d.h"
/// int main() {
///   Vector3d a(1,2,3),b(4,2,3),c;
///   a.x()= b.y();
///   a[0] = b[1];  //the same as above
///   double x=a*b; //dot product
///   c=2/b; //componentwise division
///   c=a/b; //componentwise division
///   c=a^b; //cross product
///   c=a&b; //componentwise product
///   Vector3d *d = new Vector3d[100]; // array of vectors
/// // do something with d, then save
///   FILE *fp=fopen("test","wb");
///   for(int i=0; i<100; i++) d[i].write(double(1),fp);
/// }
/// \endcode
/// \author Yuriy Turkin @ ipp mpg de
///
class Vector3d {
  double  xx; ///< first component
  double  yy; ///< second component
  double  zz; ///< third component
public:
  /// Default constructor
  Vector3d() :xx(0),yy(0),zz(0) {;}
  /// Copy constructor
  Vector3d(Vector3d const &v ) :xx(v.xx),yy(v.yy),zz(v.zz) {;}
  /// Constructor
  Vector3d(double x,double y=0.,double z=0.) :xx(x),yy(y),zz(z) {;}
  /// Constructor
  /// @param a is the array with 3 elements at least.
  template <class T> Vector3d(const T *a) {xx=double(a[0]);yy=double(a[1]);zz=double(a[2]);}
  /// Set components
  void  set    (double x,double y=0.,double z=0.) {xx=x;yy=y;zz=z;}
  /// Set components from T array
  template <class T> void setFrom(const T *a) {xx=double(a[0]);yy=double(a[1]);zz=double(a[2]);}
  /// Copy components to T array
  template <class T> void copyTo(T *a) {a[0]=T(xx);a[1]=T(yy);a[2]=T(zz);}
  /// Zero components
  void  zero() { xx = yy = zz = 0;}
  /// operator []
  /// @param row <em>=0,1,2</em> is the component number
  /// @return reference to the correspondent component
  double &operator [] (int row) {return (row==0)?xx:((row==1)?yy:zz);}
  /// operator []
  /// @param row <em>=0,1,2</em> is the component number
  /// @return reference to the correspondent component
  const double &operator [] (int row) const {return (row==0)?xx:((row==1)?yy:zz);}
  /// Get reference to the first component for reading/writing
  /// @return reference to the first component of the vector
  double &x() { return  xx;}
  /// Get the first component
  /// @return reference to the first component of the vector
  const double &x() const { return  xx;}
  /// Get reference to the second component for reading/writing
  /// @return reference to the second component of the vector
  double &y() { return  yy;}
  /// Get the second component
  /// @return reference to the second component of the vector
  const double &y() const { return  yy;}
  /// Get reference to the third component for reading/writing
  /// @return reference to the third component of the vector
  double &z() { return  zz;}
  /// Get the third component
  /// @return reference to the third component of the vector
  const double &z() const { return  zz;}
  /// Print, \b printf is used
  void print() const { printf("%g %g %g\n",xx,yy,zz); }
  /// Write to stream
  friend std::ostream& operator<< (std::ostream& o, const Vector3d & v) {
//        o.precision(4);
//        o.setf(ios::scientific, ios::floatfield);
	    o <<v.xx<<" "<<v.yy<<" "<<v.zz;
	    return o;
  }
  /// Read from stream
  friend std::istream& operator>> (std::istream& i, Vector3d & v) {
	    i >>v.xx>>v.yy>>v.zz;
	    return i;
  }
  /// Write into binary file
  /// @param x is a dummy parameter used to choose \b float
  ///    or \b double format for writing of vector components;
  ///    \e x must be \b float(1) or \b double(1).
  /// @param fp is the pointer to \b FILE structure
  /// @return 3 - the number of components actually written,
  ///    which may be less than 3 if an error occurs; see \b fwrite in <b><stdio.h></b>
  /// @note \b fwrite from <b><stdio.h></b> is used because it fast;
  ///  file must be opened before writing in binary mode: <tt>FILE *fp=fopen("test","wb");</tt>
  template <class T> size_t write(T x, FILE *fp) const {
    T w[3]; // double or float or int or char
    size_t wsize = sizeof(w[0]);
    w[0]=T(xx); w[1]=T(yy); w[2]=T(zz);
    return fwrite(w, wsize,3,fp);
  }
  /// Read from binary file
  /// @param x is a dummy parameter used to specify
  ///    the file format; \e x must be \b float(1) or \b double(1).
  /// @param fp is the pointer to \b FILE structure
  /// @return 3 - the number of components actually read,
  ///    which may be less than 3 if an error occurs or if
  ///    the end of the file is encountered; see \b fread in <b><stdio.h></b>
  /// @note \b fread from <b><stdio.h></b> is used because it fast
  ///  file must be opened before reading in binary mode: <tt>FILE *fp=fopen("test","rb");</tt>
  template <class T> size_t read(T x, FILE *fp) {
    T w[3]; // double or float or int or char
    size_t wsize = sizeof(w[0]);
    size_t cnt=fread(w, wsize,3,fp);
    xx=double(w[0]); yy=double(w[1]); zz=double(w[2]);
    return cnt;
  }
  /// Length of \b Vector3d
  /// @return x()+y()+z()
  double sum () const {return xx+yy+zz;}
  /// Length of \b Vector3d
  /// @return sqrt(x()*x()+y()*y()+z()*z())
  double abs () const {return sqrt(xx*xx+yy*yy+zz*zz);}
  /// Square of the length of \b Vector3d
  /// @return x()*x()+y()*y()+z()*z()
  double abs2() const {return (xx*xx+yy*yy+zz*zz);}
  // Normalize, set Length of \b Vector3d to 1
  //void   norm() {double a=abs(); xx/=a;yy/=a;zz/=a;}
  /// Normalize, set Length of \b Vector3d to 1
  /// @return  this vector
  Vector3d & norm() {
    double a=abs(); xx/=a;yy/=a;zz/=a;
    return( *this );
  }
  /// minimum value of component
  /// @return min (x(),y(),z())
  double mmin() const {
    double min = yy<xx?yy:xx;
    return zz<min?zz:min;
  }
  /// maximum value of component
  /// @return max (x(),y(),z())
  double mmax() const {
    double max = yy>xx?yy:xx;
    return zz>max?zz:max;
  }
  /// Find dominant coordinate
  /// @return the component number with maximum absolute value
  int  maxComponent() const {
    int i = 0;
    double max = fabs(xx);
    if(fabs(yy) > max) {max = fabs(yy); i=1;}
    if(fabs(zz) > max) i=2;
    return i;
  }
  /// Delete component.
  /// Method sets deleted component to zero and moves it to the end.
  /// <b> Vector3d(1,2,3).delComponent(2)</b> transforms to \b Vector3d(1,3,0)
  /// @param m is the component to delete, <em>m=0,1,2</em>
  void delComponent(int m) {
    switch(m) {
    case 0:
      xx = yy;
    case 1:
      yy = zz;
    case 2:
      zz = 0;
    default:
      break;
    }
  }
/// Simple assignment
  Vector3d &operator = ( Vector3d const &v ) {
    if( this != &v) {
      xx = v.xx;
      yy = v.yy;
      zz = v.zz;
    }
    return( *this );
  }
/// Simple assignment of \b double to the first component of \b Vector3d
  Vector3d &operator = ( double d ) {
      xx = d;
      yy = zz = 0.;
      return( *this );
  }
/// Addition assignment
  Vector3d &operator += ( Vector3d const &v ) {
      xx += v.xx;
      yy += v.yy;
      zz += v.zz;
      return( *this );
  }
/// Subtraction assignment
  Vector3d &operator -= ( Vector3d const &v ) {
      xx -= v.xx;
      yy -= v.yy;
      zz -= v.zz;
      return( *this );
  }
/// Multiplication assignment
  Vector3d &operator *= ( double d ) {
      xx *= d;
      yy *= d;
      zz *= d;
      return( *this );
  }
  /// Division assignment
  Vector3d &operator /= ( double d ) {
      xx /= d;
      yy /= d;
      zz /= d;
      return( *this );
  }
  /// Componentwise division assignment
  /// @return Vector3d
  Vector3d &operator /= ( Vector3d const &v ) {
      xx /= v.xx;
      yy /= v.yy;
      zz /= v.zz;
      return( *this );
  }
  /// Componentwise multiplication assignment
  Vector3d &operator &= ( Vector3d const &v ) {
      xx *= v.xx;
      yy *= v.yy;
      zz *= v.zz;
      return( *this );
  }
  /// Cross product assignment
  Vector3d operator ^= ( Vector3d const &v ) {
      double x,y;
      x  = yy*v.zz-zz*v.yy;
      y  = zz*v.xx-xx*v.zz;
      zz = xx*v.yy-yy*v.xx;
      xx = x;
      yy = y;
      return( *this );
  }
  /// Plus sign
  Vector3d operator + () const {
      return( *this );
  }
  /// Minus sign, change sign of \b Vector3d
  Vector3d operator - () const {
      return Vector3d( -xx, -yy, -zz );
  }
  /// Sum of two \b Vector3d
  friend Vector3d operator + ( Vector3d const &v1, Vector3d const &v2 ) {
      return Vector3d( v1.xx + v2.xx, v1.yy + v2.yy, v1.zz + v2.zz );
  }
  /// Subtraction of two \b Vector3d
  friend Vector3d operator - ( Vector3d const &v1, Vector3d const &v2 ) {
      return Vector3d( v1.xx - v2.xx, v1.yy - v2.yy, v1.zz - v2.zz );
  }
  /// Cross product of two \b Vector3d
  friend Vector3d operator ^ ( Vector3d const &v1, Vector3d const &v2 ) {
      return Vector3d(v1.yy*v2.zz-v1.zz*v2.yy,
                    v1.zz*v2.xx-v1.xx*v2.zz,
                    v1.xx*v2.yy-v1.yy*v2.xx );
  }
  /// Componentwise product of two \b Vector3d
  friend Vector3d operator & ( Vector3d const &v1, Vector3d const &v2 ) {
    return Vector3d(v1.xx*v2.xx, v1.yy*v2.yy, v1.zz*v2.zz);
  }
  /// Dot product of two \b Vector3d
  friend double operator * ( Vector3d const &v1, Vector3d const &v2 ) {
      return (v1.xx*v2.xx+v1.yy*v2.yy+v1.zz*v2.zz);
  }
  /// Product of \b Vector3d and \b double
  friend Vector3d operator * ( Vector3d const &v, double d ) {
      return Vector3d( v.xx*d, v.yy*d, v.zz*d );
  }
  /// Product of \b Vector3d and \b double
  friend Vector3d operator * ( double d, Vector3d const &v ) {
      return Vector3d( v.xx*d, v.yy*d, v.zz*d );
  }
  /// Devide \b Vector3d by \b double
  friend Vector3d operator / ( Vector3d const &v, double d ) {
      TESTZERO(d);
      return Vector3d( v.xx/d, v.yy/d, v.zz/d );
  }
  /// Componentwise division of two \b Vector3d
  /// @return Vector3d
  friend Vector3d operator / ( Vector3d const &v1, Vector3d const &v2 ) {
      TESTZERO(v2.xx*v2.yy*v2.zz);
      return Vector3d(v1.xx/v2.xx, v1.yy/v2.yy, v1.zz/v2.zz);
  }
  /// Componentwise division: \b Vector3d=double/Vector3d
  /// @return Vector3d
  friend Vector3d operator / ( double d, Vector3d const &v  ) {
      TESTZERO(v.xx*v.yy*v.zz);
      return Vector3d(d/v.xx, d/v.yy, d/v.zz);}
  /// Compare two \b Vector3d
  /// @return true if equal
  friend int operator == ( Vector3d const &v1, Vector3d const &v2 ) {
      return( v1.xx==v2.xx && v1.yy==v2.yy && v1.zz==v2.zz );
  }
  /// Compare the first component of \b Vector3d with \b double
  /// @return true if equal and two other components are zero
  friend int operator == ( Vector3d const &v, double d ) {
      return( v.xx==d && v.yy==0.0 && v.zz==0.0 );
  }
  /// Compare the first component of \b Vector3d with \b double
  /// @return true if equal and two other components are zero
  friend int operator == ( double d, Vector3d const &v ) {
      return( v.xx==d && v.yy==0.0 && v.zz==0.0 );
  }
  /// Compare two \b Vector3d
  /// @return true if not equal
  friend int operator != ( Vector3d const &v1, Vector3d const &v2 ) {
      return( v1.xx!=v2.xx || v1.yy!=v2.yy || v1.zz!=v2.zz );
  }
  /// Compare the first component of \b Vector3d with \b double
  /// @return true if not equal
  friend int operator != ( Vector3d const &v, double d ) {
      return( v.xx!=d || v.yy!=0.0 || v.zz!=0.0 );
  }
  /// Compare the first component of \b Vector3d with \b double
  /// @return true if not equal
  friend int operator != ( double d, Vector3d const &v ) {
      return( v.xx!=d || v.yy!=0.0 || v.zz!=0.0 );
  }
  /// Componentwise exponent of Vector3d
  inline static Vector3d exp( Vector3d const &v) {
    return Vector3d(::exp(v.x()),::exp(v.y()),::exp(v.z()));
  }
  /// Componentwise cosine of Vector3d
  inline static Vector3d cos( Vector3d const &v) {
    return Vector3d(::cos(v.x()),::cos(v.y()),::cos(v.z()));
  }
  /// Componentwise sine of Vector3d
  inline static Vector3d sin( Vector3d const &v) {
    return Vector3d(::sin(v.x()),::sin(v.y()),::sin(v.z()));
  }
  /// The method returns cylindrical coordinates of the point 
  /// assuming this object contains cartesian coordinates of the point. 
  Vector3d toCylindrical() const {
    return Vector3d(sqrt(xx*xx+yy*yy), atan2(yy,xx), zz);
  }
  /// The method returns cartesian coordinates of the point 
  /// assuming this object contains cylindrical coordinates of the point. 
  Vector3d toCartesian() const {
    return Vector3d(xx * ::cos(yy), xx * ::sin(yy), zz);
  }
  /// The method returns cartesian vector 
  /// assuming this object contains vector given in cylindrical coordinates.
  /// @param cyl is the cylindrical coordinates of the point 
  ///   at which the object vector is given
  Vector3d cylVector2Cartesian(const Vector3d &cyl) const {
    Vector3d cart;
    const Vector3d &c = *this;
    double fi = cyl[1]; 
    double cs = ::cos(fi);
    double sn = ::sin(fi);
    cart.x() = c[0]*cs-c[1]*sn;  // gx=cr*cos(fi)-cfi*sin(fi);
    cart.y() = c[0]*sn+c[1]*cs;  // gy=cr*sin(fi)+cfi*cos(fi);
    cart.z() = c[2];
    return cart;
  }
  /// The method returns vector in cylindrical coordinates 
  /// assuming this object contains vector given in cartesian coordinates.
  /// @param cyl is the cylindrical coordinates of the point 
  ///   at which the object vector is given
  Vector3d Cartesian2cylVector(const Vector3d &cyl) const {
    Vector3d v;
    const Vector3d &c = *this;
    double fi = cyl[1]; 
    double cs = ::cos(fi);
    double sn = ::sin(fi);
    v.x() =  c[0]*cs+c[1]*sn;  // vr  = cx*cos(fi)+cy*sin(fi);
    v.y() = -c[0]*sn+c[1]*cs;  // vphi=-cx*sin(fi)+cy*cos(fi);
    v.z() =  c[2];
    return v;
  }
  ///  The method returns dR^2, where dR is the distance 
  ///  between points c and this vector. 
  ///  Both vectors are given in cylindrical coordinates.
  double diffCyl2(const Vector3d &c) {
    // dR - distance between points: dR=(c-*this)
    // @return dR^2
    double dfi = c.yy-yy;
    double dR2 = xx*xx + c.xx*c.xx - 2*xx*c.xx*::cos(dfi) + (c.zz-zz)*(c.zz-zz);
    return dR2;
  }
  ///  The method returns dR.abs2(), where dR is the distance 
  ///  between points r and this vector. 
  ///  Both vectors are given in cylindrical coordinates.
  double diffCyl2x(const Vector3d &r) {
    Vector3d xyz =  this->toCartesian();
    Vector3d xyz2 =  r.toCartesian();
    xyz -= xyz2;
    return xyz.abs2();
  }

};

#undef TESTZERO

};   //namespace MConf

#endif  // MC_VECTOR3D_H
