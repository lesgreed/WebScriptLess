#ifndef MC_CArray2d_
#define MC_CArray2d_

#include <stdio.h>
#include "CVector3d.h"
#include "ngarray.h"

namespace MConf {

//*********************************************************************
/// This class provides 2d-array.
/// For internal use only.
/// \par Usage:
/// \code
/// #include "CArray2d.h"
/// int main() {
///   MConf::CArray2d<MConf::Vector3d> vertex;
///   vertex.init(100,200);
///   vertex.init(-10,20, -1,200); // resize to another shape
///   vertex(-10,200) = MConf::Vector3d(1,2,3);
///   MConf::Vector3d r = vertex(20,-1);
///   MConf::Vector3d r = vertex[20][-1];  // important! operator[] provides read only access
///   FILE *fp=fopen("test","wb");
///   vertex.write(float(1),fp);
/// }
/// \endcode
///
template <class TV> class CArray2d : public pBase::ngArray <TV> {
  int dim1;      ///< 1st dimention of the array
  int dim2;      ///< 2nd dimention of the array
  int low1;      ///< low bound of the array
  int low2;      ///< low bound of the array
public:
  /// Default constructor
  CArray2d() :pBase::ngArray <TV> () {dim1=dim2=low1=low2=0;}
  /// Constructor, creates array with dimentions \b dim1 and \b dim2.
  /// Use isOK() method to check whether the method succeeds.
  CArray2d(int dim1,int dim2) {
    resize(0, dim1-1,0,dim2-1);
  }
  /// Constructor.
  /// @param low1,upper1 are the bounds for the 1st index,
  /// @param low2,upper2 are the bounds for the 2nd index,
  /// Use isOK() method to check whether the method succeeds.
  CArray2d(int low1,int upper1,int low2,int upper2) { // array with 1st index from 'low1' to 'upper1',
    resize(low1,upper1,low2,upper2);
  }
  /// Copy constructor, creates array with the same content and shape as \b z.
  /// Use isOK() method to check whether the method succeeds.
//1  CArray2d(CArray2d const &z){ *this = z; }
  /// The method returns \b true if no error occurs for the last operation with array.
  bool isOK() const {return !this->empty();}
  /// get the the first dimention of the array.
  int size1() const {return dim1;}
  /// get the the second dimention of the array.
  int size2() const {return dim2;}
  /// get the first value of the first array index.
  int lowbound1() const {return low1;}
  /// get the first value of the second array index.
  int lowbound2() const {return low2;}
  /// get the last value of the first array index.
  int upperbound1() const {return low1+dim1-1;}
  /// get the last value of the second array index
  int upperbound2() const {return low2+dim2-1;}
  /// check whether indexes within bounds
  bool isIndexOK(int i,int j) const {
    return ((i<low1)||(j<low2)||(i>upperbound1())||(j>upperbound2()))?false:true;
  }

  /// The operator returns the value of array element.
  /// @param i,j -- element indexes
  /// @return reference to the correspondent element.
  TV &operator()      (int i, int j)       {return this->array()     [(i-low1)*dim2+(j-low2)];}
  const TV &operator()(int i, int j) const {return this->constArray()[(i-low1)*dim2+(j-low2)];}
  const TV *operator[](int i) const {return this->constArray()+(i-low1)*dim2-low2; }

  /// Addition assignment
  /// Both arrays must have the same shape
  CArray2d &operator += ( CArray2d const &v ) 
  {
    TV *This =  this->array();
    const TV *vcA = v.constArray();
    for(int i=0; i<dim1*dim2; i++) {
      This[i] += vcA[i];
    }
    return( *this );
  }
  /// Subtraction assignment
  /// Both arrays must have the same shape
  CArray2d &operator -= ( CArray2d const &v ) 
  {
    TV *This =  this->array();
    const TV *vcA = v.constArray();
    for(int i=0; i<dim1*dim2; i++) {
      This[i] -= vcA[i];
    }
    return( *this );
  }
  /// Multiplication assignment
  CArray2d &operator *= ( TV d ) 
  {
    TV *This =  this->array();
    for(int i=0; i<dim1*dim2; i++) {
      This[i] *= d;
    }
    return( *this );
  }
  /// Division assignment
  CArray2d &operator /= ( TV d ) 
  {
    TV *This =  this->array();
    for(int i=0; i<dim1*dim2; i++) {
      This[i] /= d;
    }
    return( *this );
  }
  /// Product of \b CArray2d and \b double
  CArray2d & addweighted(CArray2d const &v, TV w) 
  {
    TV *This =  this->array();
    const TV *vp  =  v.constArray();
    for(int i=0; i<dim1*dim2; i++) {
      This[i] += vp[i]*w;
    }
    return( *this );
  }

  // get direct write access to data in order to avoid copy on-write 
  TV &w(int i, int j) {return const_cast<TV &>(this->constArray()[(i-low1)*dim2+(j-low2)]);}
  
  int index(int i, int j) const { return (i-low1)*dim2+(j-low2);}
  // ptr()+ index(i,j) == &w(i,j);
  TV *ptr() {return const_cast<TV *>(this->constArray());}
  const TV *cptr() const {return this->constArray();} 
  /// The operator assigns  \b x to all elements.
  CArray2d &operator = (const TV & x) {
    this->fill(x);
    return( *this );
  }

  /// Assignment operator.
  CArray2d & operator = (const CArray2d & z) {
    if (this != &z) {
      this->dim1 = z.dim1;
      this->dim2 = z.dim2;
      this->low1 = z.low1;
      this->low2 = z.low2;
      pBase::ngArray <TV>  *ptr = this;
      ptr->operator=(z);
    }
    return *this;
  }
///////////////////////////////////////////////////////////////
  /// The method changes the dimentions of the array.
  /// Previous content of the array is not preserved.
  /// The method is the same as resize(int d1,int d2).
  /// @param d1,d2 are the new dimentions.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool init(int d1, int d2) 
  {
    return resize(d1, d2);
  }
  /// The method changes the shape of the array.
  /// Previous content of the array is not preserved.
  /// @param low1,upper1 are the bounds for the 1st index,
  /// @param low2,upper2 are the bounds for the 2nd index,
  /// indexes can take the values of bounds.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// You can also use isOK() method to check whether the method succeeds.
  /// The method is the same as resize(int low1,int upper1,int low2,int upper2).
  bool init(int low1,int upper1,int low2,int upper2) { // 1st index from 'low1' to 'upper1',
    return resize(low1,upper1,low2,upper2);
  }
  /// The method changes the shape of the array.
  /// Previous content of the array is not preserved.
  /// @param low_1,upper1 are the bounds for the 1st index,
  /// @param low_2,upper2 are the bounds for the 2nd index,
  /// indexes can take the values of bounds.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// You can also use isOK() method to check whether the method succeeds.
  bool resize(int low_1,int upper1,int low_2,int upper2) { // 1st index from 'low_1' to 'upper1',
    this->low1=low_1;                                     // 2nd index from 'low_2' to 'upper2
    this->low2=low_2;
    this->dim1 = upper1-low1+1;
    this->dim2 = upper2-low2+1;
    pBase::ngArray <TV>  *ptr = this;
    ptr->init(dim1*dim2);
    return isOK();
  }
  /// The method changes the dimentions of the array.
  /// Previous content of the array is not preserved.
  /// @param d1,d2 are the new dimentions.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool resize(int d1,int d2)
  {
    return resize(0, d1-1,0,d2-1);
  }

  /// The method frees memory
  /// This function will clear the array.
  /// Afterwards, you can't access any element or you
  /// might get a segmentation fault.
  /// \note All stored data will be forgotten.
  void clear() {
    dim1=dim2=low1=low2=0;
    (dynamic_cast<pBase::ngArray<TV>*>(this))->clear();
  }

  /// The method prints all elements of the array.
  void print() const {if(isOK()) {TV w=1; print(w);} }
  /// Write into binary file
  /// @param x is a dummy parameter used to choose \b char or \b float
  ///   or \b double format for writing of array elements;
  ///   \b x must be \b char(1) or \b float(1) or \b double(1).
  /// @param fp is the pointer to \b FILE structure.
  /// @return number of elements actually written,
  ///    which may be less than dim1*dim2(or dim1*dim2*3 for TV==Vector3d) if an error occurs;
  ///    see \b fwrite in <b><stdio.h></b>
  /// @note \b fwrite from <b><stdio.h></b> is used because it fast;
  ///  file must be opened in binary mode before writing:
  ///  <tt>FILE *fp=fopen("test","wb");</tt>
  template <class T> size_t write(T x, FILE *fp) const {
    if(!isOK()) return 0;
    fwrite(&dim1,sizeof(int),1,fp);
    fwrite(&dim2,sizeof(int),1,fp);
    fwrite(&low1,sizeof(int),1,fp);
    fwrite(&low2,sizeof(int),1,fp);
    TV w=1;
    return write(w,x,fp);
  }
  /// Read from binary file the array that was written by write() method.
  /// @param x is a dummy parameter used to specify
  ///    the file format; \b x must be \b char(1) or \b float(1) or \b double(1).
  /// @param fp is the pointer to \b FILE structure
  /// @return number of elements actually read.
  ///    which may be less than dim1*dim2(or dim1*dim2*3 for TV==Vector3d)
  ///    if an error occurs or if the end of the file is encountered;
  ///    see \b fread in <b><stdio.h></b>
  /// @note \b fread from <b><stdio.h></b> is used because it fast
  /// file must be opened in binary mode before reading: <tt>FILE *fp=fopen("test","rb");</tt>
  template <class T> size_t read(T x, FILE *fp) {
    fread(&dim1,sizeof(int),1,fp);
    fread(&dim2,sizeof(int),1,fp);
    fread(&low1,sizeof(int),1,fp);
    fread(&low2,sizeof(int),1,fp);
    int upper1 = dim1+low1-1;
    int upper2 = dim2+low2-1;
    if(!resize(low1,upper1,low2,upper2)) return 0;
    TV w=1;
    return read(w,x,fp);
  }
  /// The same as write() but the shape of the array are not saved, only dimensions are saved.
  /// @return number of elements actually written.
  template <class T> size_t write0(T x, FILE *fp) const {
    if(!isOK()) return 0;
    fwrite(&dim1,sizeof(int),1,fp);
    fwrite(&dim2,sizeof(int),1,fp);
    TV w=1;
    return write(w,x,fp);
  }
  /// The method reads the array that was written by write0() method.
  /// @return number of elements actually read.
  template <class T> size_t read0(T x, FILE *fp) {
    fread(&dim1,sizeof(int),1,fp);
    fread(&dim2,sizeof(int),1,fp);
    if(!resize(dim1,dim2)) return 0;
    TV w=1;
    return read(w,x,fp);
  }

private:
  TV &This      (int i, int j)       {return this->array()     [i*dim2+j]; }
  const TV &This(int i, int j) const {return this->constArray()[i*dim2+j]; }

  /// helper function for print()
  void print(double u) const { // TV is float or double
    int i,j;
    for(i=0; i<dim1; i++) {
      for(j=0; j<dim2; j++)
        printf("%g ", This(i,j));
      printf("\n");
    }
  }
  void print(int u) const { // TV is int or char
    int i,j;
    for(i=0; i<dim1; i++) {
      for(j=0; j<dim2; j++)  printf("%d ", This(i,j));
      printf("\n");
    }
  }
  /// helper function for print()
  void print(Vector3d &u) const {
    int i,j;
    for(i=0; i<dim1; i++)
      for(j=0; j<dim2; j++)  This(i,j).print();
  }
  /// helper function for write(T x, FILE *fp)
  template <class T> size_t write(double u, T x, FILE *fp) const {//TV is float,double or char
    int i,j;
    size_t cnt=0;
    T w;
    size_t wsize = sizeof(w);
    for(i=0; i<dim1; i++)
      for(j=0; j<dim2; j++) {
        w=T(This(i,j));
        cnt+=fwrite(&w, wsize,1,fp);
      }
    return cnt;
  }
  /// helper function for write(T x, FILE *fp)
  template <class T> size_t write(Vector3d &u, T x, FILE *fp) const {
    int i,j;
    size_t cnt=0;
    for(i=0; i<dim1; i++)
      for(j=0; j<dim2; j++) cnt+=This(i,j).write(x,fp);
    return cnt;
  }
  /// helper function for read(T x, FILE *fp)
  template <class T> size_t read(double u, T x, FILE *fp) {  // TV is float or double or char
    int i,j;
    size_t cnt=0;
    T w;
    size_t wsize = sizeof(w);
    for(i=0; i<dim1; i++)
      for(j=0; j<dim2; j++) {
        cnt+=fread(&w, wsize,1,fp);
        This(i,j)=w;
      }
    return cnt;
  }
  /// helper function for read(T x, FILE *fp)
  template <class T> size_t read(Vector3d &u, T x, FILE *fp) {
    int i,j;
    size_t cnt=0;
    for(i=0; i<dim1; i++)
      for(j=0; j<dim2; j++) cnt+=This(i,j).read(x,fp);
    return cnt;
  }

};

};   //namespace MConf

#endif
