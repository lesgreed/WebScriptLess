#ifndef MC_CArray1d_
#define MC_CArray1d_

#include <stdio.h>
#include "CVector3d.h"
#include "ngarray.h"

namespace MConf {

//*********************************************************************
/// This class provides 1d-array.
/// For internal use only.
/// \par Usage:
/// \code
/// #include "CArray1d.h"
/// int main() {
///   MConf::CArray1d<MConf::Vector3d> vertex;
///   vertex.resize(100);    // resize to array[0:99]
///   vertex.resize(-10,20); // resize to array[-10:20]
///   vertex.at(-10) = MConf::Vector3d(1,2,3);
///   vertex.rw()[-10] = MConf::Vector3d(1,2,3);
///   MConf::Vector3d r = vertex[20];  // important! operator[] provides read only access
///   FILE *fp=fopen("test","wb");
/// }
/// \endcode
///
template <class TV> class CArray1d : public pBase::ngArray <TV> {
  int dim1;      ///< 1st dimention of the array
  int low1;      ///< low bound of the array
public:
  /// Default constructor
  CArray1d() :pBase::ngArray <TV> () {dim1=low1=0;}
  /// Constructor, creates array with dimention \b dim1.
  /// Use isOK() method to check whether the method succeeds.
  CArray1d(int dim1) {
    resize(0, dim1-1);
  }
  /// Constructor.
  /// @param low1,upper1 are the bounds for the array index
  /// Use isOK() method to check whether the method succeeds.
  CArray1d(int low1,int upper1) { // array with 1st index from 'low1' to 'upper1',
    resize(low1,upper1);
  }
  /// Copy constructor, creates array with the same content and shape as \b z.
  /// Use isOK() method to check whether the method succeeds.
  CArray1d(CArray1d const &z) { *this = z; }
  /// Copy constructor, creates array1d with the same content \b z.
  /// Use isOK() method to check whether the method succeeds.
  CArray1d(TV * z, int size) {
    this->assign(z, size);
  }
  /// The method creates array from TV * z.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool assign(const TV * z, int size) {
    this->resize(size);
    for(int i=0; i<size; i++) this->at(i) = z[i];
    return this->isOK(); 
  }
  /// The method scales array.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool scale(const TV fact) {
    for(int i=0; i<size(); i++) this->at(i) *= fact;
    return this->isOK(); 
  }
  /// The method returns \b true if no error occurs for the last operation with array.
  bool isOK() const {return !this->empty();}
  /// get the the first dimention of the array.
  int size() const {return dim1;}
  /// get the first value of the first array index.
  int lowbound() const {return low1;}
  /// get the last value of the first array index.
  int upperbound() const {return low1+dim1-1;}
  /// check whether indexes within bounds
  bool isIndexOK(int i) const {
    return (i<low1||i>upperbound())?false:true;
  }
  /// Perform a binary search of the segment in which given value lies.
  /// @param u is the value to search for the segment
  /// @return left end of the segment in which u lies.
  /// \note This array must be in ascending order.
 int bsearch(double u) const
  {
    if(this->empty()) return 0;
    int L = 0;       // Left bracket
    int R = dim1-1;  // Right bracket
    const TV * x = this->constArray();

    if(u<=x[0]) return L;
    else if(u>=x[R]) L=R-1;
    else
      while(L+1 < R) {  // bisection
        int k = (L + R)/2;
        if(u < x[k]) R = k;
        else L = k;
      }
    return L+low1;
  }

  /// Linear interpolation of this array in x with Xa as the abscissa array.
  /// @param x is the abscissa point at which the interpolated value is needed.
  /// @param Xa is the abscissa array.
  /// @return the interpolated value of this array in the point x.
  /// \note The arrays (Xa and this) must have the same shape, 
  /// test for consistensy is not perfomed.
  TV interp1(double x, const CArray1d<double> &Xa) const
  {
    int L = Xa.bsearch(x);
    int R = L+1;
    const TV * y = this->constArray()-low1; 
    TV yp= (y[R]-y[L])/(Xa[R]-Xa[L]); // derivative in point x
    return(y[L] + yp*(x-Xa[L]));
  }

  /// Quadratic interpolation of this array in x with Xa as the abscissa array.
  /// @return the interpolated value of this array in the point x.
  /// \note The arrays (Xa and this) must have the same shape, 
  /// test for consistensy is not perfomed.
  TV interp2(double x, const CArray1d<double> &Xa) const
  {
    int N = Xa.size();
    int L = Xa.bsearch(x);
    if(L+2>=N) L=N-3;
    const TV * y = this->constArray()-low1; 
    double y0 = y[L  ];
    double y1 = y[L+1];
    double y2 = y[L+2];
    double x0 = Xa[L  ];
    double x1 = Xa[L+1];
    double x2 = Xa[L+2];
    double c1=(y1-y0)/(x1-x0);
    double c2=(y2-y1)/(x2-x1);
    c2=(c2-c1)/(x2-x0);
    return  y0 + (c1 + c2*(x-x1))*(x-x0);
  }

  /// The operator returns the reference to the array element.
  /// @param i is the element index
  /// @return reference to the correspondent array element.
  TV &              at(int i)       {return this->array()[i-low1]; }
  /// The operator returns the address to the memory allocated for the array.
  /// @return  (address of allocated memory) - (lowest array index),
  /// so that rw()[lowestIndex] is the first array element.
  TV *              rw()            {return this->array() - low1; }  // read/write access
  const TV &operator[](int i) const {return this->constArray()[i-low1]; }
  //TV &      operator[](int i)       {return this->array()[i-low1]; }
  // allow direct write access to data in order to avoid copy on-write 
  TV &              w(int i)        {return const_cast<TV &>(this->constArray()[i-low1]); } 
  TV *              w()             {return const_cast<TV *>(this->constArray() - low1); } 
  
  /// The operator assigns  \b x to all elements.
  CArray1d &operator = (const TV & x) {
    this->fill(x);
    return( *this );
  }

  /// Assignment operator.
  CArray1d & operator = (const CArray1d & z) {
    if (this != &z) {
      this->dim1 = z.dim1;
      this->low1 = z.low1;
      pBase::ngArray <TV>  *ptr = this;
      ptr->operator=(z);
    }
    return *this;
  }
  
///////////////////////////////////////////////////////////////
  /// The method changes the dimentions of the array.
  /// Previous content of the array is not preserved.
  /// The method is the same as resize(int dim1,int dim2).
  /// @param dim is the new dimention.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool init(int dim) 
  {
    return resize(dim);
  }
  /// The method changes the shape of the array.
  /// Previous content of the array is not preserved.
  /// @param low,upper are the bounds for the 1st index.
  /// @param val is the value to initialize array.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// You can also use isOK() method to check whether the method succeeds.
  /// The method is the same as resize(int low1,int upper1,const TV & val = TV()).
  bool init(int low,int upper, const TV & val = TV()) 
  {
    return resize(low,upper,val);
  }
  /// The method changes the shape of the array.
  /// Previous content of the array is not preserved.
  /// @param low,upper are the bounds for the 1st index.
  /// @param val is the value to initialize array.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// You can also use isOK() method to check whether the method succeeds.
  bool resize(int low,int upper, const TV & val = TV()) { // 1st index from 'low' to 'upper',
    this->low1 = low;                         
    this->dim1 = upper-low+1;
//  pBase::ngArray <TV>  *ptr = this;
//  ptr->init(dim1*dim2);
    (dynamic_cast<pBase::ngArray<TV>*>(this))->init(0,dim1-1,val);    
    return isOK();
  }
  /// The method changes the dimentions of the array.
  /// Previous content of the array is not preserved.
  /// @param dim is the new dimention.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  bool resize(int dim)
  {
    return resize(0, dim-1);
  }

  /// The method frees memory
  /// This function will clear the array.
  /// Afterwards, you can't access any element or you
  /// might get a segmentation fault.
  /// \note All stored data will be forgotten.
  void clear() {
    dim1=low1=0;
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
  ///    which may be less than dim1(or dim1*3 for TV==Vector3d) if an error occurs;
  ///    see \b fwrite in <b><stdio.h></b>
  /// @note \b fwrite from <b><stdio.h></b> is used because it fast;
  ///  file must be opened in binary mode before writing:
  ///  <tt>FILE *fp=fopen("test","wb");</tt>
  template <class T> size_t write(T x, FILE *fp) const {
    if(!isOK()) return 0;
    fwrite(&dim1,sizeof(int),1,fp);
    fwrite(&low1,sizeof(int),1,fp);
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
    fread(&low1,sizeof(int),1,fp);
    int upper1 = dim1+low1-1;
    if(!resize(low1,upper1)) return 0;
    TV w=1;
    return read(w,x,fp);
  }
  /// The same as write() but the shape of the array are not saved, only dimensions are saved.
  /// @return number of elements actually written.
  template <class T> size_t write0(T x, FILE *fp) const {
    if(!isOK()) return 0;
    fwrite(&dim1,sizeof(int),1,fp);
    TV w=1;
    return write(w,x,fp);
  }
  /// The method reads the array that was written by write0() method.
  /// @return number of elements actually read.
  template <class T> size_t read0(T x, FILE *fp) {
    fread(&dim1,sizeof(int),1,fp);
    if(!resize(dim1)) return 0;
    TV w=1;
    return read(w,x,fp);
  }

private:
  TV &This      (int i)       {return this->array()     [i]; }
  const TV &This(int i) const {return this->constArray()[i]; }

  /// helper function for print()
  void print(double u) const { // TV is float or double
    for(int i=0; i<dim1; i++) 
      printf("%g ", This(i));
    printf("\n");
    
  }
  void print(int u) const { // TV is int or char
    for(int i=0; i<dim1; i++) 
      printf("%d ", This(i));
    printf("\n");
  }
  /// helper function for print()
  void print(Vector3d &u) const {
    for(int i=0; i<dim1; i++)
      This(i).print();
  }
  /// helper function for write(T x, FILE *fp)
  template <class T> size_t write(double u, T x, FILE *fp) const {//TV is float,double or char
    size_t cnt=0;
    T w;
    size_t wsize = sizeof(w);
    for(int i=0; i<dim1; i++) {
      w=T(This(i));
      cnt+=fwrite(&w, wsize,1,fp);
    }
    return cnt;
  }
  /// helper function for write(T x, FILE *fp)
  template <class T> size_t write(Vector3d &u, T x, FILE *fp) const {
    size_t cnt=0;
    for(int i=0; i<dim1; i++)
      cnt+=This(i).write(x,fp);
    return cnt;
  }
  /// helper function for read(T x, FILE *fp)
  template <class T> size_t read(double u, T x, FILE *fp) {  // TV is float or double or char
    size_t cnt=0;
    T w;
    size_t wsize = sizeof(w);
    for(int i=0; i<dim1; i++) {
      cnt+=fread(&w, wsize,1,fp);
      This(i)=w;
    }
    return cnt;
  }
  /// helper function for read(T x, FILE *fp)
  template <class T> size_t read(Vector3d &u, T x, FILE *fp) {
    size_t cnt=0;
    for(int i=0; i<dim1; i++)
      cnt+=This(i).read(x,fp);
    return cnt;
  }

};

};   //namespace MConf

#endif
