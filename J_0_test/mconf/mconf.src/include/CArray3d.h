#ifndef MC_CArray3d_
#define MC_CArray3d_

#include <stdio.h>
#include "CArray2d.h"
#include "CVector3d.h"

namespace MConf {

//*********************************************************************
/// This class provides 3d-array.
/// For internal use only.
/// Strictly speaking this is not 3d-array, but 1d-array of CArray2d.
/// The class is designed for internal use in class C3dMesh.
///
///
template <class TV> class CArray3d : public pBase::ngArray< CArray2d <TV> > {
  int low1;  ///< low bound of the first index of the array
  int dim1;  ///< 1st dimention of the array
public:
  /// Default constructor
  CArray3d() :pBase::ngArray < CArray2d <TV> > () {low1=dim1=0;}
  /// Constructor, creates array with dimentions \b dim1,\b dim2,and \b dim2.
  /// Use isOK() method to check whether the method succeeds.
  CArray3d(int dim1,int dim2,int dim3) {low1=0;resize(dim1,dim2,dim3); }
  // Copy constructor, creates array with the same content and shape as \b z.
  // Use isOK() method to check whether the method succeeds.
  CArray3d(CArray3d const &z) { *this = z; }
  /// The method returns \b true if no error occurs in the last operation with array.
  bool isOK() const {return !this->empty();}
  /// The operator returns the reference to the CArray2d.
  /// @param i is the index of CArray2d.
  /// @return reference to the CArray2d.
  /// \note method ngArray::array() is very specific, see ngArray::detach() where the whole 
  ///  array may need to be copied if clones exist, but the clones are not affected! 
  CArray2d<TV> &operator[] (int i) {return this->array()[i-low1]; }
  /// The operator returns the constant reference to the CArray2d.
  const CArray2d<TV> &operator[] (int i) const {return this->constArray()[i-low1]; }
  /// The operator assigns  \b x to all element.
  CArray3d &operator = (const TV & x)
  {
    if(isOK())
      for(int i=0; i<dim1; i++) this->array()[i] = x;
    return( *this );
  }

  /// Assignment operator.
  CArray3d & operator = (const CArray3d & z) {
    if (this != &z) {
      // we need new container
      this->resize(z.dim1);
      this->low1 = z.low1;
      this->dim1 = z.dim1;
      for(int i=0; i<dim1; i++) 
        (this->array()[i]).operator=(z.constArray()[i]);
    }
    return *this;
  }

  /// The method changes the shape of the array.
  /// Previous content of the array is not preserved.
  /// @param low,upper are the bounds for the 1st index,
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// You can also use isOK() method to check whether the method succeeds.
  /// \note the memory for other dimentions must be allocated as follows:
  /// \code
  /// MConf::CArray3d<double> B;
  /// B.resize(-2,10);
  /// B[-1].resize(-10,20, -1,200); // resize the plane number -1
  /// \endcode
  bool resize(int low,int upper) { // array with 1st index from 'low1' to 'upper1',
    this->low1 = low;
    this->dim1 = upper-low+1;
    this->init(dim1);
    return isOK();
  }
  /// The method changes the fisrt dimention of the array.
  /// Previous content of the array is not preserved.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// You can also use isOK() method to check whether the method succeeds.
  /// \note the memory for other dimentions must be allocated as follows:
  /// \code
  /// MConf::CArray3d<double> B;
  /// B.resize(20);
  /// B[5].resize(-10,20, -1,200); // resize the plane number 5
  /// \endcode
  bool resize(int dim1) {
    return resize(0,dim1-1);
  }
  /// The method changes the dimentions of the array.
  /// Previous content of the array is not preserved.
  /// @param dim1,dim2,dim3 are the new dimentions.
  /// @return  if the method succeeds, the return value is \b true; \b false otherwise.
  /// You can also use isOK() method to check whether the method succeeds.
  bool resize(int dim1,int dim2,int dim3)
  {
    resize(dim1);
    if(dim2&&dim3)
      for(int i=0; i<dim1; i++)
        if(false==(this->array()[i]).resize(dim2,dim3)) {
          this->clear();
          break;
        }
    return isOK();
  }
  /// The method frees memory
  /// This function will clear the array.
  /// Afterwards, you can't access any element or you
  /// might get a segmentation fault.
  /// \note All stored data will be forgotten.
  void clear() {
    dim1=low1=0;
    (dynamic_cast< pBase::ngArray< CArray2d <TV> > *>(this))->clear();

  }
  /// The method prints all elements of the array.
  void print() const {
    if(isOK()) for(int i=0; i<dim1; i++) (this->constArray()[i]).print();
  }

  /// Write into binary file
  /// @param x is a dummy parameter used to choose \b char or \b float
  ///   or \b double format for writing of array elements;
  ///   \b x must be \b char(1) or \b float(1) or \b double(1).
  /// @param fp is the pointer to \b FILE structure.
  /// @param i1,i2 if i1<i2 then write planes number from i1 to i2, where the plane number
  ///   means the values of the 1st index of array.
  /// @return number of elements actually written,which may be less than
  ///   dim1*dim2*dim3(or dim1*dim2*dim3*3 for TV==Vector3d) if an error occurs;
  ///    see \b fwrite in <b><stdio.h></b>
  /// @note \b fwrite from <b><stdio.h></b> is used because it fast;
  ///  file must be opened in binary mode before writing:
  ///  <tt>FILE *fp=fopen("test","wb");</tt>
  template <class T> size_t write(T x, FILE *fp,int i1=0,int i2=0) const {
    if(!isOK()) return 0;
    size_t cnt=0;
    fwrite(&dim1,sizeof(int),1,fp);
    fwrite(&low1,sizeof(int),1,fp);
    fwrite(&i1,  sizeof(int),1,fp);
    fwrite(&i2,  sizeof(int),1,fp);
    if(i1>=i2)
      for(int i=0; i<dim1; i++) {
        size_t n = (this->constArray()[i]).write(x,fp);
        if(n==0) return 0;
        cnt+=n;
      }
    else
      for(int i=i1; i<=i2; i++) {
        size_t n = (*this)[i].write(x,fp);
        if(n==0) return 0;
        cnt+=n;
      }

    return cnt;
  }
  /// Read from binary file the array that was written by write() method.
  /// @param x is a dummy parameter used to specify
  ///    the file format; \b x must be \b char(1) or \b float(1) or \b double(1).
  /// @param fp is the pointer to \b FILE structure
  /// @return number of elements actually read.
  ///  which may be less than dim1*dim2*dim3(or dim1*dim2*dim3*3 for TV==Vector3d)
  ///  if an error occurs or if the end of the file is encountered;
  ///  see \b fread in <b><stdio.h></b>
  /// @note \b fread from <b><stdio.h></b> is used because it fast
  /// file must be opened in binary mode before reading:
  ///   <tt>FILE *fp=fopen("test","rb");</tt>
  template <class T> size_t read(T x, FILE *fp) {
    int i1=0,i2=0;
    fread(&dim1,sizeof(int),1,fp);
    fread(&low1,sizeof(int),1,fp);
    fread(&i1,  sizeof(int),1,fp);
    fread(&i2,  sizeof(int),1,fp);
    int upper1 = dim1+low1-1;
    if(!resize(low1,upper1)) return 0;
    size_t cnt=0;
    if(i1>=i2)
      for(int i=0; i<dim1; i++) {
        size_t n = (this->array()[i]).read(x,fp);
        if(n==0) return 0;
        cnt+=n;
      }
    else
      for(int i=i1; i<=i2; i++) {
        size_t n = (*this)[i].read(x,fp);
        if(n==0) return 0;
        cnt+=n;
      }
    return cnt;
  }
};

};   //namespace MConf

#endif
