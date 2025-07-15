/***************************************************************************
                          ngarray.h  -  description
                             -------------------
    begin                : Tue May 25 2004
    copyright            : (C) 2004 by Stefan Zegenhagen
    based on             : garray.h (C) 2000 by Andreas Werner (andreas.werner@ipp.mpg.de)
    email                : stefan.zegenhagen@ipp.mpg.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef MC__PBASE_NGARRAY_H__
#define MC__PBASE_NGARRAY_H__

#include <string>
#include <iostream>
#include <algorithm>
#include <limits.h>

#ifdef __INTEL_COMPILER
// These pragmas are to quiet icc about "virtual function override intended?"
#pragma warning( push )
#pragma warning( disable : 1125 )
#endif

namespace MConf {

/// Container namespace
namespace pBase {

/** \brief Array of objects.
 * \ingroup tools
  *
  * This template class is a generic array container for object of type \c T.
  * The objects stored should meet the requirements that usual containers oppose
  * to them. They must have a default constructor, a copy operator and copy constructor.
  * Furthermore, it is required that the objects are capable of writing and reading themselves to
  * and from a stream using the usual indirection operators \c << and \c >>.
  *
  * A great difference from normal C/C++ array types is that this class has an arbitrary index for the
  * first element. Where a normal C/C++ array type has it's indices always starting at zero, you can
  * construct or initialize a ngArray with an arbitrary start index. Using zero-based indices is
  * supported and as easy as usual. Examples:
  * \code
    pBase::ngArray<double> a(-10, 10)   // valid indices range from -10 to 10, size() = 21
    pBase::ngArray<double> b(20, 30)    // valid indices range from 20 to 30, size() = 11
    a.init(0, 9);                       // valid indices range from 0 to 9, size() = 10
    b.init(10);                         // same as above, but easier and more intuitive
    \endcode
  *
  * This class uses implicit data sharing, it is not expensive to copy it or pass it around as argument.
  * Writing to the array may, however, take some time because the data is automatically copied if it is
  * shared between different instances.
  *
  * To reshape or resize the array, the init() function is supplied. Note that this will destroy
  * the data stored. The array will be initialized with the new dimensions provided. When an array
  * created or reshaped, all elements are initialized with a default value. You may pass the initialization
  * value in the constructor or init() function.
  * \note If you create a ngArray<T> where T is an integer type, you should avoid passing the initialization
  * value. Use fill() instead to initialize the array. This is because the compiler can't destinguish between
  * ngArray<T>::init(size, init) and ngArray<T>::init(minIdx, maxIdx).
  * \note If you create a ngArray<double>, you should explicitely convert the initialization value to a \c double
  * when calling the constructors and init() functions. See the following example for explanation.
  \code
    int someVal = 1;
    ...
    ngArray<double> a(size, someVal);               // will use ngArray<double>::ngArray<double>(int, int, const double & = double())
    ngArray<double> b(size, 1);                     // will use ngArray<double>::ngArray<double>(int, int, const double & = double())
    ngArray<double> c(size, double(someVal));       // correct
  \endcode
  * @author Stefan Zegenhagen, Andreas Werner
  */

template < typename T > class ngArray {
public:
    /*! \brief Create ngArray.
     *
     * Generates a new empty ngArray. If you try
     * to access any element of this array, you will probably
     * get a segmentation fault.
     */
    ngArray() {
        p = new Rep();
    }

    /*! \brief Create and initialize ngArray.
     *
     * This creates a new ngArray with valid indices ranging from \c 0 to
     * \c size - 1 and fills all elements with a constant value.
     * If you don't specify a fill value, the default constructor for that type is called.
     */
    ngArray(size_t size, const T & val = T())
    {
        p = new Rep(0, int(size-1));
        fill(val);
    }

    /*! \brief Create and initialize ngArray.
     * This creates a new ngArray and fills all elements with a constant value.
     * If you don't specify a fill value, the default constructor for that type is called.
     */
    ngArray(int min, int max, const T &val = T())
    {
        p = new Rep(min, max);
        fill(val);
    }

    /**\brief Destructor.*/
    virtual ~ngArray() {
        p->deref();
    }

    /**\brief Copy constructor.*/
    ngArray(const ngArray & that) {
        p = that.p;
        p->ref();
    }

    /**\brief Assignment operator.*/
    ngArray & operator = (const ngArray & that) {
        if (p != that.p) {
            p->deref();
            p = that.p;
            p->ref();
        }
        return *this;
    }

    /** \brief Gives access to element \c i,
     *
     * Returns a reference of the element at \c i for reading and writing.
     * \c i should be limited to \c min ... \c max, otherwise it is made to
     * fit into that range.
     * \note You should check that the array is not empty. Accessing elements
     * from an empty ngArray will cause a segmentation fault.
     * \note This operator may take some time to return the element. It may be necessary
     * to copy the whole array.
     **/
    T & operator[] (int i)
    {
        detach();
        return (*p)[i];
    }

    /** \brief Gives direct access to element \c i, copy on write is not performed.
     *
     * Returns a reference of the element at \c i for reading and writing.
     * \c i should be limited to \c min ... \c max, otherwise it is made to
     * fit into that range.
     * \note You should check that the array is not empty. Accessing elements
     * from an empty ngArray will cause a segmentation fault.
     * \note This method doesn't detach this object memory from clones, 
     * all ngArray clones will see changes after assignment A.wd(i) = x, i.e.:
     *  ngArray A,B;
     *  A.init(10,1); // all A elements are 1
     *  B = A;
     *  A.wd(1) = 10; // B[1] will be also 10;
     * compare with  T & operator[] (int i)
     **/
    T & wd(int i)
    {
        return (*p)[i];
    }

    /**\brief Returns a constant reference to the element at \c i.
     *
     * The function returns a constant reference to the element at \c i.
     * \c i should be limited to \c min ... \c max, otherwise it is made to
     * fit into that range.
     * \note You should check that the array is not empty. Accessing elements
     * from an empty ngArray will cause a segmentation fault.
     * \note This operator usually returns quite fast, no copying will be
     * performed.
     **/
    const T & operator[] (int i) const
    {
        return (*p)[i];
    }

     /*! \brief clear the array
     *
     * This function will clear the array.
     * Afterwards, you can't access any element or you
     * might get a segmentation fault.
     * \note All stored data will be forgotten.
     */
    void clear()
    {
        p->deref();
        p = new Rep();
    }

    /*! \brief reshape and initialize the array
     *
     * This function will reshape the array to range from index \c 0 to
     * \c size - 1. All stored data will be forgotten and the new array
     * elements are initialized with a constant value. If you don't specify
     * an initialization value, the default constructor for the type will be
     * called.
     */
     void init(size_t size, const T &val = T())
     {
         p->deref();
         p = new Rep(0, int(size-1));
         fill(val);
     }

    /*! \brief reshape and initialize the array
     *
     * This function will reshape the array to the provided dimensions.
     * All stored data will be forgotten and the new array elements are
     * initialized with a constant value. If you don't specify an initialization
     * value, the default constructor for the type will be called.
     */
    void init(int min, int max, const T &val = T())
    {
        p->deref();
        p  = new Rep(min, max);
        fill(val);
    }

    /*! \brief initialize all elements with a constant value. */
    virtual void fill(const T &val)
    {
        detach();
        for (size_t i = 0; i < p->n; i++)
            p->data[i] = val;
    }

    /** \brief Returns actual array size.*/
    size_t size() const
    {
        return p->n;
    }

    /*! \brief return the minimum index.
     * This function returns the minimum valid index to the array. If the array is empty,
     * this value is \c INT_MAX. */
    int minIdx() const
    {
        return p->minIdx;
    }

    /*! \brief return the maximum index.
     * This function returns the maximum valid index to the array. If the array is empty,
     * the value \c INT_MAX-1 is returned - just to make constructs like
     * \code
     * for (int i = array.minIdx(); i < array.maxIdx(); i++)
     *      do_something(array[i]);
     * \endcode
     * work as expected. It's better you check that the array is not empty before trying
     * to access it.
     */
    int maxIdx() const
    {
        return p->maxIdx;
    }

    /*! \brief return whether the array is empty */
    bool empty() const
    {
        return (size() == 0);
    }

    /*! \brief Returns memory pointer
     *
     * This function may be expensive because the whole array may need to be copied.
     * This function may also be dangerous because the access to the elements is not
     * protected.
     */
    T *array()
    {
        detach();
        return p->data;
    }

    /*! \brief returns constant pointer to the elements
     *
     * This function returns a pointer to the memory where the elements
     * are stored. It is just like a normal C style array. If the array
     * is empty, you get a NULL pointer returned. This function may be
     * dangerous because the access to the elements is not protected.
     */
    const T * constArray() const
    {
        return p->data;
    }

protected:

    /*! \brief Protected class holds data
     *
     */
    struct Rep {
#ifdef NGARRAY_CHECK_RANGE
  template <class W> const W& min (const W& x, const W& y) { return (x<y)?x:y; }
  template <class W> const W& max (const W& x, const W& y) { return (x>y)?x:y; }
  #define SET_IN_RANGE(a,mi,ma)   a=max(mi,min(ma,a))
#else
  #define SET_IN_RANGE(a,mi,ma)
#endif
        Rep() {
            minIdx = INT_MAX;
            maxIdx = INT_MAX - 1;
            n = 0;
            data = 0;
            refcnt = 1;
        }
        Rep(int min, int max) {
            if (max < min)
                max = min;
            minIdx = min;
            maxIdx = max;
            n = max - min + 1;
            data = new T[n];
            refcnt = 1;
        }

        Rep(const Rep &that) {
            minIdx = that.minIdx;
            maxIdx = that.maxIdx;
            n = that.n;
            data = new T[n];
            for (size_t i = 0; i < n; i++)
                data[i] = that.data[i];
            refcnt = 1;
        }
        ~Rep()
        {
          if(data) {delete[] data; data=0;}
        }

        T& operator[](int idx) {
            SET_IN_RANGE(idx,minIdx,maxIdx);
            return data[idx - minIdx];
        }

        const T& operator[](int idx) const {
            SET_IN_RANGE(idx,minIdx,maxIdx);
            return data[idx - minIdx];
        }

        void ref() {
            ++refcnt;
        }

        void deref() {
          if(--refcnt==0) {
            minIdx = INT_MAX;
            maxIdx = INT_MAX-1;
            n = 0;
            delete this;
          }
        }
        int minIdx;
        int maxIdx;
        size_t  n;
        size_t  refcnt;
        T*  data;
    };

    void detach() {
        if (p->refcnt > 1) {
            Rep *old = p;
            p = new Rep (*old);
            old->deref();
        }
    }
    Rep *p;
};


/*! \brief a 2 dimensional matrix type
 * \ingroup tools
 *
 * This class is based on ngArray<T> and represents a 2-dimensional
 * matrix. The init() and fill() functions work as expected for a
 * matrix - they perform the initialization of all elements in the
 * array. Indexing is the same as for the C-style matrices - i.e.
 * \c matrix[i][j] will return a reference to the appropriate value.
 *
 * This class is optimized to be passed around as argument and handles
 * copy operations quite fast. It uses data sharing between multiple
 * instances as much as possible and so tries to minimize the memory
 * footprint.
 *
 * Because the data sharing happens transparently to the user, assignments
 * of values to an element might require deep copies of some part of the data
 * that is shared between multiple instances.
 *
 * It is possible for the matrix to be empty (to contain not a single element)
 * In this case, you must not use the access operators - they will fail and
 * cause a segmentation fault. Rather, you may check whether the matrix is empty
 * by calling empty(). If you want the matrix to contain at least a single
 * element, use a proper constructor or the init() function and pass a value
 * of 0 to each index.
 */
template < class T > class ng2Matrix : public ngArray < ngArray < T > >
{
    typedef typename ngArray < ngArray < T > >::Rep Rep;

public:
    /*! \brief constructor. Creates an empty matrix that contains no elements. */
    ng2Matrix() {}
    /*! \brief constructor. Creates a matrix that is shaped according to the dimensions supplied.
     * All elements are initialized with the supplied \c val. If no \c val is given, the default constructor
     * for the type is called. */
    ng2Matrix(int min1, int max1, int min2, int max2, const T & val = T())
    {
        init(min1, max1, min2, max2, val);
    }

    /*! \brief constructor. Creates a matrix that has the dimensions \c (0:s1-1, \c 0:s2-1).
     * All elements are initialized with the supplied \c val. If no \c val is given, the default constructor
     * for the type is called. */
    ng2Matrix(size_t s1, size_t s2, const T &val = T())
    {
        init(s1, s2, val);
    }

    /*! \brief initialize all elements with a constant value.
     *
     * This function will initialize all elements of the matrix with
     * a constant value. Data sharing is used as much as possible between
     * the rows of the matrix.
     */
    void fill(const T &val)
    {
        // unfortunately, all that "this"-ing is necessary for GCC-3.4 because
        // that does not find unqualified members in template base classes that
        // don't depend on a template parameter argument
        this->detach();
        this->p->data[0].fill(val);
        for (size_t i = 1; i < this->p->n; i++)
            this->p->data[i] = this->p->data[0];
    }


    /*! \brief reshape and initialize the matrix
     *
     * This function makes the matrix forget all stored data,
     * reshape itself to the new dimensions, and initialize all
     * elements with the constant value supplied. Data sharing is used
     * as much as possible to keep the memory requirements low.
     * If no initialization value is given, the default constructor
     * for the type is called.
     */
    void init(int min1, int max1, int min2, int max2, const T &val = T())
    {
        // unfortunately, all that "this"-ing is necessary for GCC-3.4 because
        // that does not find unqualified members in template base classes that
        // don't depend on a template parameter argument
        this->p->deref();
        this->p = new Rep(min1, max1);
        ngArray < T > tmp2(min2, max2, val);
        for (size_t i = 0; i < this->p->n; i++)
            this->p->data[i] = tmp2;
    }
    /*! \brief reshape and initialize the matrix
     *
     * This function makes the matrix forget all stored data,
     * reshape itself to the new dimensions, and initialize all
     * elements with the constant value supplied. Data sharing is used
     * as much as possible to keep the memory requirements low.
     * If no initialization value is given, the default constructor
     * for the type is called.
     * \note After calling this function, the array will have the dimensions
     * (0:s1-1,0:s2-1).
     */
    void init(size_t s1, size_t s2, const T &val = T())
    {
        init(0, int(s1-1), 0, int(s2-1), val);
    }
};

/*! \brief A 3-dimensional matrix type
 * \ingroup tools
 *
 * This class represents a 3-dimensional matrix type that uses implicit and
 * transparent data sharing. It is much like ng2Matrix<T>, except that it has
 * a third dimension.
 *
 * It is possible for the matrix to be empty (to contain not a single element)
 * In this case, you must not use the access operators - they will fail and
 * cause a segmentation fault. Rather, you may check whether the matrix is empty
 * by calling empty(). If you want the matrix to contain at least a single
 * element, use a proper constructor or the init() function and pass a value
 * of 0 to each index.
 */
template < class T > class ng3Matrix : public ngArray < ng2Matrix < T > >
{
    typedef typename ngArray < ng2Matrix < T > >::Rep Rep;
public:
    ng3Matrix() {}
    ng3Matrix(int min1, int max1, int min2, int max2, int min3, int max3, const T & val = T())
    {
        init(min1, max1, min2, max2, min3, max3, val);
    }
    ng3Matrix(size_t s1, size_t s2, size_t s3, const T &val = T())
    {
        init(s1, s2, s3, val);
    }

    void fill(const T & val)
    {
        // unfortunately, all that "this"-ing is necessary for GCC-3.4 because
        // that does not find unqualified members in template base classes that
        // don't depend on a template parameter argument
        this->detach();
        this->p->data[0].fill(val);
        for (size_t i = 1; i < this->p->n; i++)
            this->p->data[i] = this->p->data[0];
    }

    void init(int min1, int max1, int min2, int max2, int min3, int max3, const T & val = T())
    {
        // unfortunately, all that "this"-ing is necessary for GCC-3.4 because
        // that does not find unqualified members in template base classes that
        // don't depend on a template parameter argument
        this->p->deref();
        this->p = new Rep(min1, max1);
        ng2Matrix < T > tmp2(min2, max2, min3, max3, val);
        for (size_t i = 0; i < this->p->n; i++)
            this->p->data[i] = tmp2;
    }
    void init(size_t s1, size_t s2, size_t s3, const T &val = T())
    {
        init(0, s1-1, 0, s2-1, 0, s3-1, val);
    }
};


using std::ostream;
using std::istream;
using std::ws;
using std::string;
#if defined(__MWERKS__)  // Metrowerks CodeWarrior
  using std::strtol;
#endif

template <class T>
ostream & operator<< (ostream & o, const ngArray<T> & that)
{
    if (that.empty()) {
        o << "ngArray( ) ";
        return o;
    }
    o << "ngArray(" << that.minIdx() << "," << that.maxIdx() << ")=[ ";
    for (int i = that.minIdx(); i <= that.maxIdx(); i++) {
        o << that[i];
        if (i < that.maxIdx()) o << " , ";
    }
    o << " ]";
    return o;
}



template<>
inline ostream & operator<< (ostream & o, const ngArray<string> & that)
{
    if (that.empty()) {
        o << "ngArray( ) ";
        return o;
    }
    o << "ngArray(" << (that.minIdx()) << "," << (that.maxIdx()) << ")=[ ";
    for (int i = that.minIdx(); i <= that.maxIdx(); i++) {
        o << '"' << that[i] << '"';
        if (i < that.maxIdx()) o << " , ";
    }
    o << " ]";
    return o;
}



template <class T>
istream & operator>> (istream & i, ngArray<T> & po)
{
    /* if the input stream is invalid, return without modifying it */
    if (!i)
        return i;
    char *buf = new char[256];
    int min, max;
    i >> ws;
    i.read(buf, 8);
    buf[8] = 0;
    /* read the header and check it. If a mismatch is found, set the fail bit
     * of this stream */
    if (string(buf) != "ngArray(")
        i.setstate(i.failbit);
    if (i) {
        /* check whether the array is empty */
        i >> ws;
        char ch = i.get();
        if (ch == ')') {
            /* the array is empty */
            po.clear();
            delete[] buf;
            return i;
        }
        i.putback(ch);
        /* read the first number - it's the minimum index */
        i.getline(buf, 256, ',');
        char *tail;
        min = strtol(buf, &tail, 10);
        if (*tail != 0)
            i.setstate(i.failbit);
    }
    if (i) {
        /* read the next number - it is the maximum index and ends with a closing ')' */
        i.getline(buf, 256, ')');
        char *tail;
        max = strtol(buf, &tail, 10);
        if (*tail != 0)
            i.setstate(i.failbit);
    }
    /* read the next characters: they should be a '=' sign and an opening bracket */
    if (i) {
        i.read(buf, 2);
        buf[2] = 0;
        if (string(buf) != "=[")
            i.setstate(i.failbit);
    }
    if (i) {
        /* next come all the elements - we read them into a temporary array
         * if something goes wrong while reading, we don't modify the target array! */
        ngArray<T> that(min, max);
        for (int idx = min; idx <= max; idx++) {
            /* read the element */
            i >> that[idx];
            /* in case of failure, exit the loop */
            if (!i)
                break;
            i >> ws;
            /* get the next character */
            char ch = i.get();
            /* if it's not the last index, it should be a comma */
            if (idx < max) {
                if (ch != ',')
                    i.setstate(i.failbit);
            } else {
                /* otherwise, it must be a closing bracket */
                if (ch != ']')
                    i.setstate(i.failbit);
            }
            /* again, if the stream is not valid anymore, exit the loop */
            if (!i)
                break;
        }
        /* done reading the elements - if the stream is valid, we may copy it to
         * the target array. this is usually quite fast (i.e. only a pointer is copied. */
        if (i)
            po = that;
    }
    /* free the temporary buffer and return a reference to the input stream */
    delete[] buf;
    return i;
}



template<>
inline istream & operator>> (istream & i, ngArray<string> & po)
{
    /* if the input stream is invalid, return without modifying it */
    if (!i)
        return i;
    char *buf = new char[256];
    int min, max;
    i >> ws;
    i.read(buf, 8);
    buf[8] = 0;
    /* read the header and check it. If a mismatch is found, set the fail bit
     * of this stream */
    if (!(string(buf) == "ngArray("))
        i.setstate(i.failbit);
    if (i) {
        i >> ws;
        char ch = char(i.get());
        if (ch == ')') {
            /* the array is empty */
            po.clear();
            delete[] buf;
            return i;
        }
        i.putback(ch);
        /* read the first number - it's the minimum index */
        i.getline(buf, 256, ',');
        char *tail;
        min = strtol(buf, &tail, 10);
        if (*tail != 0)
            i.setstate(i.failbit);
    }
    if (i) {
        /* read the next number - it is the maximum index and ends with a closing ')' */
        i.getline(buf, 256, ')');
        char *tail;
        max = strtol(buf, &tail, 10);
        if (*tail != 0)
            i.setstate(i.failbit);
    }
    /* read the next characters: they should be a '=' sign and an opening bracket */
    if (i) {
        i.read(buf, 2);
        buf[2] = 0;
        if (string(buf) != "=[")
            i.setstate(i.failbit);
    }
    if (i) {
        /* next come all the elements - we read them into a temporary array
         * if something goes wrong while reading, we don't modify the target array! */
        ngArray<string> that(min, max);
        for (int idx = min; idx <= max; idx++) {
            /* read the element */
            i >> ws;
            if (i.get() == '"') {
                string temp_str;
                std::getline(i, temp_str, '"');
                that[idx] = temp_str;
            } else
                i.setstate(i.failbit);
            /* in case of failure, exit the loop */
            if (!i)
                break;
            i >> ws;
            /* get the next character */
            char ch = char(i.get());
            /* if it's not the last index, it should be a comma */
            if (idx < max) {
                if (ch != ',')
                    i.setstate(i.failbit);
            } else {
                /* otherwise, it must be a closing bracket */
                if (ch != ']')
                    i.setstate(i.failbit);
            }
            /* again, if the stream is not valid anymore, exit the loop */
            if (!i)
                break;
        }
        /* done reading the elements - if the stream is valid, we may copy it to
         * the target array. this is usually quite fast (i.e. only a pointer is copied. */
        if (i)
            po = that;
    }
    /* free the temporary buffer and return a reference to the input stream */
    delete[] buf;
    return i;
}

}

};   //namespace MConf

#ifdef __INTEL_COMPILER
#pragma warning( pop )
#endif

#endif // __PBASE_NGARRAY_H__
