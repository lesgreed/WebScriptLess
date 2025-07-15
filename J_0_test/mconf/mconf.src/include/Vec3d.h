#ifndef MC_VEC3D_H
#define MC_VEC3D_H

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "CVector3d.h"

namespace MConf {

#ifdef _DEBUG
 #define  TESTZERO(x) if (!(x)) std::cerr<<"devide by zero in Vec3d operator/()"<<std::endl;
#else 
 #define  TESTZERO(x) 
#endif

//**********************************************************************
/// 2015, Oct 28 Templated class.
///
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
/// #include "Vec3d.h"
/// int main() {
///   Vec3d<float> a(1,2,3),b(4,2,3),c;
///   a.x()= b.y();
///   a[0] = b[1];  //the same as above
///   float x=a*b; //dot product
///   c=2/b; //componentwise division
///   c=a/b; //componentwise division
///   c=a^b; //cross product
///   c=a&b; //componentwise product
///   Vec3d *d = new Vec3d[100]; // array of vectors
/// // do something with d, then save
///   FILE *fp=fopen("test","wb");
///   for(int i=0; i<100; i++) d[i].write(float(1),fp);
/// }
/// \endcode
/// \author Yuriy Turkin @ ipp mpg de
///
template < typename TV > class Vec3d {
  TV  xx; ///< first component
  TV  yy; ///< second component
  TV  zz; ///< third component
public:
  friend Vec3d<int>;
  friend Vec3d<long long>;
  friend Vec3d<double>;
  friend Vec3d<float>;
  friend Vector3d;
  /// Default constructor
  /// must always set zero
  Vec3d() :xx(0),yy(0),zz(0) {;}
  /// Copy constructor
  Vec3d(const Vector3d  &v ) :xx(TV(v.x())),yy(TV(v.y())),zz(TV(v.z())) {;}
  template <class T> Vec3d(Vec3d<T> const &v ) :xx(TV(v.xx)),yy(TV(v.yy)),zz(TV(v.zz)) {;}
  /// Constructor
  template <class T> Vec3d(T x,T y,T z) :xx((TV)x),yy((TV)y),zz((TV)z) {;}
  Vec3d(double x,double y,double z) :xx((TV)x),yy((TV)y),zz((TV)z) {;}
  /// Constructor
  /// @param a is the array with 3 elements at least.
  template <class T> Vec3d(const T *a) {xx=TV(a[0]);yy=TV(a[1]);zz=TV(a[2]);}
  /// Set components
  template <class T> void  set  (T x,T y,T z) {xx=(TV)x;yy=(TV)y;zz=(TV)z;}
  //void  set  (double x,double y,double z) {xx=(TV)x;yy=(TV)y;zz=(TV)z;}
  /// Set components from T array
  template <class T> void set    (const T *a) {xx=TV(a[0]);yy=TV(a[1]);zz=TV(a[2]);}
  template <class T> void setFrom(const T *a) {xx=TV(a[0]);yy=TV(a[1]);zz=TV(a[2]);}
  /// Copy components to T array
  template <class T> void copy  (T *a) const {a[0]=T(xx);a[1]=T(yy);a[2]=T(zz);}
  template <class T> void copyTo(T *a) const {a[0]=T(xx);a[1]=T(yy);a[2]=T(zz);}
  /// Zero components
  void  zero() { xx = yy = zz = 0;}
  /// operator []
  /// @param row <em>=0,1,2</em> is the component number
  /// @return reference to the correspondent component
  TV &operator [] (int row) {return (row==0)?xx:((row==1)?yy:zz);}
  /// operator []
  /// @param row <em>=0,1,2</em> is the component number
  /// @return reference to the correspondent component
  const TV &operator [] (int row) const {return (row==0)?xx:((row==1)?yy:zz);}
  /// Get reference to the first component for reading/writing
  /// @return reference to the first component of the vector
  TV &x() { return  xx;}
  /// Get the first component
  /// @return reference to the first component of the vector
  const TV &x() const { return  xx;}
  /// Get reference to the second component for reading/writing
  /// @return reference to the second component of the vector
  TV &y() { return  yy;}
  /// Get the second component
  /// @return reference to the second component of the vector
  const TV &y() const { return  yy;}
  /// Get reference to the third component for reading/writing
  /// @return reference to the third component of the vector
  TV &z() { return  zz;}
  /// Get the third component
  /// @return reference to the third component of the vector
  const TV &z() const { return  zz;}
  /// Print, \b printf is used
  void print() const { printf("%g %g %g\n",xx,yy,zz); }
  /// Write to stream
  friend std::ostream& operator<< (std::ostream& o, const Vec3d & v) {
//        o.precision(4);
//        o.setf(ios::scientific, ios::floatfield);
	    o <<v.xx<<" "<<v.yy<<" "<<v.zz;
	    return o;
  }
  /// Read from stream
  friend std::istream& operator>> (std::istream& i, Vec3d & v) {
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
    xx=TV(w[0]); yy=TV(w[1]); zz=TV(w[2]);
    return cnt;
  }
  /// Length of \b Vec3d
  /// @return x()+y()+z()
  TV sum () const {return xx+yy+zz;}
  /// Length of \b Vec3d
  /// @return sqrt(x()*x()+y()*y()+z()*z())
  TV abs () const {return sqrt(xx*xx+yy*yy+zz*zz);}
  /// Square of the length of \b Vec3d
  /// @return x()*x()+y()*y()+z()*z()
  TV abs2() const {return (xx*xx+yy*yy+zz*zz);}
  // Normalize, set Length of \b Vec3d to 1
  //void   norm() {TV a=abs(); xx/=a;yy/=a;zz/=a;}
  /// Normalize, set Length of \b Vec3d to 1
  /// @return  this vector
  Vec3d & norm() {
    TV a=abs(); xx/=a;yy/=a;zz/=a;
    return( *this );
  }
  /// minimum value of component
  /// @return min (x(),y(),z())
  TV mmin() const {
    TV min = yy<xx?yy:xx;
    return zz<min?zz:min;
  }
  /// maximum value of component
  /// @return max (x(),y(),z())
  TV mmax() const {
    TV max = yy>xx?yy:xx;
    return zz>max?zz:max;
  }
  /// Find dominant coordinate
  /// @return the component number with maximum absolute value
  int  maxComponent() const {
    int i = 0;
    TV max = fabs(xx);
    if(fabs(yy) > max) {max = fabs(yy); i=1;}
    if(fabs(zz) > max) i=2;
    return i;
  }
  /// Delete component.
  /// Method sets deleted component to zero and moves it to the end.
  /// <b> Vec3d(1,2,3).delComponent(2)</b> transforms to \b Vec3d(1,3,0)
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
  Vec3d &operator = ( Vec3d const &v ) {
    if( this != &v) {
      xx = v.xx;
      yy = v.yy;
      zz = v.zz;
    }
    return( *this );
  }
/// Simple assignment of \b TV to the first component of \b Vec3d
  Vec3d &operator = ( TV d ) {
      xx = d;
      yy = zz = 0.;
      return( *this );
  }
/// Addition assignment
  Vec3d &operator += ( Vec3d const &v ) {
      xx += v.xx;
      yy += v.yy;
      zz += v.zz;
      return( *this );
  }
/// Subtraction assignment
  Vec3d &operator -= ( Vec3d const &v ) {
      xx -= v.xx;
      yy -= v.yy;
      zz -= v.zz;
      return( *this );
  }
/// Multiplication assignment
  Vec3d &operator *= ( TV d ) {
      xx *= d;
      yy *= d;
      zz *= d;
      return( *this );
  }
  /// Division assignment
  Vec3d &operator /= ( TV d ) {
      xx /= d;
      yy /= d;
      zz /= d;
      return( *this );
  }
  /// Componentwise division assignment
  /// @return Vec3d
  Vec3d &operator /= ( Vec3d const &v ) {
      xx /= v.xx;
      yy /= v.yy;
      zz /= v.zz;
      return( *this );
  }
  /// Componentwise multiplication assignment
  Vec3d &operator &= ( Vec3d const &v ) {
      xx *= v.xx;
      yy *= v.yy;
      zz *= v.zz;
      return( *this );
  }
  /// Cross product assignment
  Vec3d operator ^= ( Vec3d const &v ) {
      TV x,y;
      x  = yy*v.zz-zz*v.yy;
      y  = zz*v.xx-xx*v.zz;
      zz = xx*v.yy-yy*v.xx;
      xx = x;
      yy = y;
      return( *this );
  }
  /// Plus sign
  Vec3d operator + () const {
      return( *this );
  }
  /// Minus sign, change sign of \b Vec3d
  Vec3d operator - () const {
      return Vec3d( -xx, -yy, -zz );
  }
  /// Sum of two \b Vec3d
  friend Vec3d operator + ( Vec3d const &v1, Vec3d const &v2 ) {
      return Vec3d( v1.xx + v2.xx, v1.yy + v2.yy, v1.zz + v2.zz );
  }
  /// Subtraction of two \b Vec3d
  friend Vec3d operator - ( Vec3d const &v1, Vec3d const &v2 ) {
      return Vec3d( v1.xx - v2.xx, v1.yy - v2.yy, v1.zz - v2.zz );
  }
  /// Cross product of two \b Vec3d
  friend Vec3d operator ^ ( Vec3d const &v1, Vec3d const &v2 ) {
      return Vec3d(v1.yy*v2.zz-v1.zz*v2.yy,
                    v1.zz*v2.xx-v1.xx*v2.zz,
                    v1.xx*v2.yy-v1.yy*v2.xx );
  }
  /// Componentwise product of two \b Vec3d
  friend Vec3d operator & ( Vec3d const &v1, Vec3d const &v2 ) {
    return Vec3d(v1.xx*v2.xx, v1.yy*v2.yy, v1.zz*v2.zz);
  }
  /// Dot product of two \b Vec3d
  friend TV operator * ( Vec3d const &v1, Vec3d const &v2 ) {
      return (v1.xx*v2.xx+v1.yy*v2.yy+v1.zz*v2.zz);
  }
  /// Product of \b Vec3d and \b TV
  friend Vec3d operator * ( Vec3d const &v, TV d ) {
      return Vec3d( v.xx*d, v.yy*d, v.zz*d );
  }
  /// Product of \b Vec3d and \b TV
  friend Vec3d operator * ( TV d, Vec3d const &v ) {
      return Vec3d( v.xx*d, v.yy*d, v.zz*d );
  }
  /// Devide \b Vec3d by \b TV
  friend Vec3d operator / ( Vec3d const &v, TV d ) {
      TESTZERO(d);
      return Vec3d( v.xx/d, v.yy/d, v.zz/d );
  }
  /// Componentwise division of two \b Vec3d
  /// @return Vec3d
  friend Vec3d operator / ( Vec3d const &v1, Vec3d const &v2 ) {
      TESTZERO(v2.xx*v2.yy*v2.zz);
      return Vec3d(v1.xx/v2.xx, v1.yy/v2.yy, v1.zz/v2.zz);
  }
  /// Componentwise division: \b Vec3d=TV/Vec3d
  /// @return Vec3d
  friend Vec3d operator / ( TV d, Vec3d const &v  ) {
      TESTZERO(v.xx*v.yy*v.zz);
      return Vec3d(d/v.xx, d/v.yy, d/v.zz);}
  /// Compare two \b Vec3d
  /// @return true if equal
  friend int operator == ( Vec3d const &v1, Vec3d const &v2 ) {
      return( v1.xx==v2.xx && v1.yy==v2.yy && v1.zz==v2.zz );
  }
  /// Compare the first component of \b Vec3d with \b TV
  /// @return true if equal and two other components are zero
  friend int operator == ( Vec3d const &v, TV d ) {
      return( v.xx==d && v.yy==0.0 && v.zz==0.0 );
  }
  /// Compare the first component of \b Vec3d with \b TV
  /// @return true if equal and two other components are zero
  friend int operator == ( TV d, Vec3d const &v ) {
      return( v.xx==d && v.yy==0.0 && v.zz==0.0 );
  }
  /// Compare two \b Vec3d
  /// @return true if not equal
  friend int operator != ( Vec3d const &v1, Vec3d const &v2 ) {
      return( v1.xx!=v2.xx || v1.yy!=v2.yy || v1.zz!=v2.zz );
  }
  /// Compare the first component of \b Vec3d with \b TV
  /// @return true if not equal
  friend int operator != ( Vec3d const &v, TV d ) {
      return( v.xx!=d || v.yy!=0.0 || v.zz!=0.0 );
  }
  /// Compare the first component of \b Vec3d with \b TV
  /// @return true if not equal
  friend int operator != ( TV d, Vec3d const &v ) {
      return( v.xx!=d || v.yy!=0.0 || v.zz!=0.0 );
  }
  /// Componentwise exponent of Vec3d
  inline static Vec3d exp( Vec3d const &v) {
    return Vec3d(::exp(v.x()),::exp(v.y()),::exp(v.z()));
  }
  /// Componentwise cosine of Vec3d
  inline static Vec3d cos( Vec3d const &v) {
    return Vec3d(::cos(v.x()),::cos(v.y()),::cos(v.z()));
  }
  /// Componentwise sine of Vec3d
  inline static Vec3d sin( Vec3d const &v) {
    return Vec3d(::sin(v.x()),::sin(v.y()),::sin(v.z()));
  }
  /// The method returns cylindrical coordinates of the point 
  /// assuming this object contains cartesian coordinates of the point. 
  Vec3d toCylindrical() const {
    return Vec3d(sqrt(xx*xx+yy*yy), atan2(yy,xx), zz);
  }
  /// The method returns cartesian coordinates of the point 
  /// assuming this object contains cylindrical coordinates of the point. 
  Vec3d toCartesian() const {
    return Vec3d(xx * ::cos(yy), xx * ::sin(yy), zz);
  }
  /// The method returns cartesian vector 
  /// assuming this object contains vector given in cylindrical coordinates.
  /// @param cyl is the cylindrical coordinates of the point 
  ///   at which the object vector is given
  Vec3d cylVector2Cartesian(const Vec3d &cyl) const {
    Vec3d cart;
    const Vec3d &c = *this;
    TV fi = cyl[1]; 
    TV cs = ::cos(fi);
    TV sn = ::sin(fi);
    cart.x() = c[0]*cs-c[1]*sn;  // gx=cr*cos(fi)-cfi*sin(fi);
    cart.y() = c[0]*sn+c[1]*cs;  // gy=cr*sin(fi)+cfi*cos(fi);
    cart.z() = c[2];
    return cart;
  }
  /// The method returns vector in cylindrical coordinates 
  /// assuming this object contains vector given in cartesian coordinates.
  /// @param cyl is the cylindrical coordinates of the point 
  ///   at which the object vector is given
  Vec3d Cartesian2cylVector(const Vec3d &cyl) const {
    Vec3d v;
    const Vec3d &c = *this;
    TV fi = cyl[1]; 
    TV cs = ::cos(fi);
    TV sn = ::sin(fi);
    v.x() =  c[0]*cs+c[1]*sn;  // vr  = cx*cos(fi)+cy*sin(fi);
    v.y() = -c[0]*sn+c[1]*cs;  // vphi=-cx*sin(fi)+cy*cos(fi);
    v.z() =  c[2];
    return v;
  }
  ///  The method returns dR^2, where dR is the distance 
  ///  between points c and this vector. 
  ///  Both vectors are given in cylindrical coordinates.
  TV diffCyl2(const Vec3d &c) {
    // dR - distance between points: dR=(c-*this)
    // @return dR^2
    TV dfi = c.yy-yy;
    Vec3d dR(c.xx-xx*::cos(dfi),xx*::sin(dfi),c.zz-zz);
    return dR.abs2();
 //   TV dR2 = xx*xx + c.xx*c.xx - 2*xx*c.xx*::cos(dfi) + (c.zz-zz)*(c.zz-zz);
 //   return dR2;
  }
  ///  The method returns dR.abs2(), where dR is the distance 
  ///  between points r and this vector. 
  ///  Both vectors are given in cylindrical coordinates.
  TV diffCyl2x(const Vec3d &r) {
    Vec3d xyz =  this->toCartesian();
    Vec3d xyz2 =  r.toCartesian();
    xyz -= xyz2;
    return xyz.abs2();
  }

};

//typedef Vec3d<double> Vector3d;

#undef TESTZERO

};   //namespace MConf

#endif  // MC_VEC3D_H
