#ifndef MC__THREADS_H
#define MC__THREADS_H

#include <vector>
#include <iostream>
#include "threadtypes.h"

// undefine this if you don't want debug print
//#define _TESTPRINTX

#if (defined( _DEBUGx ) || defined(_TESTPRINTX) )
 #define  TESTPRINT(x)  {(x);}
#else 
 #define  TESTPRINT(x) {;}
#endif

#ifdef _WIN32
  #include <windows.h>       // Windows threading
  #include <process.h>

  namespace MConf {

  ////typedef HANDLE  uiThread;                 // see threadtypes.h
  ////typedef HANDLE  uiMutex;                  // see threadtypes.h
  ////#define uiThreadFunc  unsigned __stdcall  // see threadtypes.h

//******************************************************************
  inline uiThread ui_beginthreadex(unsigned(__stdcall *f)(void *), void* p) {
    #ifdef _MT
      return (uiThread)_beginthreadex(NULL,0,(unsigned(__stdcall *)(void *))f, p, 0, 0);
    #else
//      #pragma message( "Singlethread application")
      f(p);
      return 0;  // single thread if th==0
    #endif
  };
  // Wait until thread terminates
  inline void ui_waitthread(uiThread thread)  {
    #ifdef _MT
      if(thread==0) return;
      WaitForSingleObject(thread, INFINITE);
      CloseHandle(thread); // Destroy the thread object
    #endif
  };
  #ifdef _MSC_VER
  #pragma warning( push )
  #pragma warning( disable : 4311  ) // warning C4311: 'Typumwandlung': Zeigerverkuerzung von 'void *' zu 'unsigned int'
  #endif
  inline void ui_endthreadex() {
    #ifdef _MT
     _endthread();
    #endif
  };
  #ifdef _MSC_VER
  #pragma warning( pop )
  #endif
//******************************************************************
  inline uiMutex ui_mutex_init() { return CreateMutex (NULL, FALSE, NULL);};
  inline int ui_mutex_lock  (uiMutex mutex) { return WaitForSingleObject(mutex,INFINITE);};
  inline int ui_mutex_unlock(uiMutex mutex) { return ReleaseMutex       (mutex); };
  inline void ui_mutex_destroy(uiMutex mutex) { CloseHandle(mutex); };

  };   //namespace MConf


#else                        // POSIX threading
  #include <pthread.h>

  namespace MConf {

  ////typedef pthread_t        uiThread;  // see threadtypes.h
  ////typedef pthread_mutex_t* uiMutex;   // see threadtypes.h
  ////#define uiThreadFunc  void *        // see threadtypes.h
//******************************************************************
  inline uiThread ui_beginthreadex(void *(*f)(void *), void* p)
  {
    uiThread th;
    //pthread_attr_t attr;
    //pthread_attr_init(&attr);
    //pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    //int rc = pthread_create((pthread_t*)&th, &attr, f, p);
    //pthread_attr_destroy(&attr);
    int rc = pthread_create((pthread_t*)&th, NULL, f, p);
    return rc==0?th:0;  // error if th==0
  };
  // Wait until thread terminates
  inline void ui_waitthread(uiThread thread)  {
    if(thread==0) return;
    pthread_join(thread, NULL);
  };
  inline void ui_endthreadex() { pthread_exit((void *)0); }
//******************************************************************
  inline uiMutex ui_mutex_init()
  {
    uiMutex t=new pthread_mutex_t;
    if(t) {
      if(pthread_mutex_init(t,NULL)!=0) { delete t; t=0;}
    }
    return t;
  };
  inline int ui_mutex_lock  (uiMutex mutex) { return pthread_mutex_lock  (mutex); };
  inline int ui_mutex_unlock(uiMutex mutex) { return pthread_mutex_unlock(mutex); };
  inline void ui_mutex_destroy(uiMutex mutex)
  {
    if(mutex==0) return;
    pthread_mutex_destroy(mutex);
    delete mutex;
  }

};   //namespace MConf

#endif

namespace MConf {

//****************************************************************************
//
// mutex is not provided; user must create it in caller routines if needed
template <class T> class Threads
{
  typedef  void(T::*TExeFun)(int i1, int i2, bool setup);
  T *obj;   
  TExeFun Exe;
  std::string thName;
  void setup__(int i1, int i2) {
    (obj->*Exe)(i1, i2,true); // call the member function of the obj
  }
  void execute__(int i1, int i2) {
    (obj->*Exe)(i1, i2,false);
  }
  struct Args {Threads *that; int i1, i2; };
  static uiThreadFunc EntryPoint(void *A)
  {
    Args *arg  =(Args*)A;
    Threads *that = arg->that;
    int i1 = arg->i1;
    int i2 = arg->i2;
    delete arg;
    that->execute__(i1, i2);
    ui_endthreadex();
    return 0;
  }
public:
  Threads(T *ob, TExeFun exe, const char* name=0) : obj(ob), Exe(exe) {if(name) thName=name;}   
  virtual ~Threads() {;}
  // 
  void run(int low, int upper, int nThreads) 
  { 
    //nThreads = nThreads<=0?2:nThreads;    // number of threads
#ifdef _WIN32
    #ifndef _MT
    nThreads = 1;   //must be singlethread application
    #endif
#endif
    const int NPPT = 5;  //min number of points per thread
    setup__(low, upper);  
    if(nThreads<=1) { execute__(low, upper); return; }
    int size = upper-low+1;
    int nppt = size/nThreads;         // number of points (per thread) to be handled by one thread
    if(nppt<=NPPT)  nThreads=size/NPPT;
    if(nThreads<=1) { execute__(low, upper); return; }
    nppt = size/nThreads;
    nppt = (nppt<=0)?1:nppt;
    //nppt = nThreads==1?size:nppt; 
    std::vector<uiThread> hThreads;   // to keep handles of threads
    int iTh = -1;
    for(int i1=low,i2=0; ;) {
      iTh++;
      i2 = i1+nppt;
      i2 = i2>upper?upper:i2;
      Args *arg = new Args;
      arg->that = this;
      arg->i1 = i1;
      arg->i2 = i2;
      uiThread th = ui_beginthreadex(EntryPoint,(void*)arg);      
      if(th==uiThread(0)) {    // if thread is not started 
        execute__(i1, i2);
        delete arg;        
      }
      else {
        hThreads.push_back(th); 
        TESTPRINT( std::cerr<<thName<<" starts the thread #"<<iTh<<" with handle="<<th<<":  for(i="<<i1<<";i<="<<i2<<";i++) {...}"<<std::endl<<std::flush ) ;
      }
      i1 = i2+1;
      if(i1>upper) break;
    }
    TESTPRINT( std::cerr<<std::endl<<std::flush ) ;
    for(int j=0; j<(int)hThreads.size(); j++)  
      ui_waitthread(hThreads[j]);    // Wait until threads terminate
  }
};

#if 0 
// 23.08.2021 must be revised 
// the following is the example of how to use class Threads

//****************************************************************************
// this class implements multithreading for coordinate transformation 
using MConf::Vector3d;
using std::vector;
class Mag2cyl {
  C3dMesh * mconf;
  vector<Vector3d> * cyl;
  const vector<Vector3d> * mag;
  uiMutex mutex;
  int numThreads;
public:
  Mag2cyl(C3dMesh *mc,const vector<Vector3d> *m,vector<Vector3d> *c) 
           :mconf(mc), cyl(c), mag(m) {};

  void calculate(int nThreads) {
    numThreads = nThreads;
    Threads<Mag2cyl> threads(this, &Mag2cyl::calculateEx);  
    int low = 0, upper = int(mag->size()-1);
    mutex = ui_mutex_init();
    threads.run(low, upper, nThreads); 
    ui_mutex_destroy(mutex);
  }

  void calculateEx(int low, int upper, bool setup)
  {
    if(setup) {  // setup part; it allocates arrays if needed
      int size = upper-low+1;
      cyl->resize(size);
    }
    else
    {  // exe part; it will be called from threads to fill arrays
      C3dMesh *mcf = mconf;
      if(numThreads>1) {
        mcf = new C3dMesh;
        ui_mutex_lock(mutex);  
        *mcf = *mconf;        // need a new copy of C3dMesh to be thread safe
        ui_mutex_unlock(mutex);
      }
      C3dMesh &mc = *mcf;
      
      vector<Vector3d> &cyl = *(this->cyl);
      const vector<Vector3d> &mag = *(this->mag);
      
      for(int i=low; i<=upper; i++) {
        cyl[i] = mc.mag2cyl(mag[i]);     
      }
      if(numThreads>1) delete mcf;
    }
  }

};

// in main programm

  C3dMesh mc;
  vector<Vector3d> m;
  vector<Vector3d> c;
  Mag2cyl m2c(&mc, &m, &c);
  m2c.calculate(numOfThreads);

#endif

template < class T, class TR > class Threads2
{
  typedef  void(T::*TExeFun)(TR *objR, int i1, int i2, bool setup);
  struct Args {Threads2 *This; TR *objR; int i1; int i2; uiMutex mutex; };  
public:
  Threads2(T *ob, TExeFun exe, TR *oRHS,const char* name=0) : obj(ob), Exe(exe), objR(oRHS) {if(name) thName=name;}  
  virtual ~Threads2() {;}
  /// @param[in] low
  /// @param[in] upper; low and upper are the first and last indexes in arrays to be filled, low<=upper.
  /// @param[in] nThreads  is the number of threads to run
  /// @param[in] NPPT is the min number of array indexes per thread to be filled
  /// \see  void(T::*TExeFun)(TR *objR, int i1, int i2, bool setup),
  ///      where [i1 - i2] is the portion interval of the full interval [low - upper] if setup is false
  void run(int low, int upper, int nThreads, int NPPT=5) // low<=upper
  { 
#ifdef _WIN32
    #ifndef _MT
    nThreads = 1;   //must be singlethread application
    #endif
#endif
    int size = upper-low+1;
    if(size<=0) return;  // indexes must be in ascending(increasing) order
    if(NPPT<1) NPPT = 1;
    if(size<=NPPT)  { setup__(objR,low, upper); execute__(objR,low, upper); return; }  //if1
    if(nThreads<=1) { setup__(objR,low, upper); execute__(objR,low, upper); return; }  //if3
    if(nThreads>size) nThreads=size;                                                   //if4
    int nppt = size/nThreads; // nppt>=1, number of points per thread, to be handled by one thread
    //here: NPPT >= 1, size > NPPT , 2 <= nThreads <= size, nppt >= 1 
    if(nppt<NPPT)  nThreads=size/NPPT+1; // here nThreads >= 2, see comment //if1
    nppt = size/nThreads;
    if(nppt==0) nppt=1;
    nppt--;
    int nTh = 0;
    for(int i1=low,i2=low-1; ;) {  // calculate nTh - number of theads
      i1 = i2+1;   
      if(i1>upper) break;
      i2 = i1+nppt; 
      i2 = i2>upper?upper:i2;
      nTh++;
    }
    setup__(objR,low, upper);
    TR &oR = *objR;
    std::vector<TR> oRc(nTh,oR);  // vector of clones
    std::vector<Args> args(nTh);  // vector of Args
    std::vector<uiThread> hThreads;   // to keep handles of threads
    mutex = ui_mutex_init();
    for(int i1=low,i2=low-1,ith=0; ith<nTh; ith++) {
      i1 = i2+1;   
      if(i1>upper) break;
      i2 = i1+nppt; 
      i2 = i2>upper?upper:i2;
      Args *arg = &args[ith];
      arg->objR = &oRc[ith];   
      arg->This = this;
      arg->i1 = i1;
      arg->i2 = i2;
      arg->mutex = mutex;
    }
    for(int j=0; j<nTh; j++) {
      Args *arg = &args[j];
      uiThread th = ui_beginthreadex(EntryPoint,(void*)(arg));      
      if(th==uiThread(0)) {    // if thread is not started 
        execute__(arg->objR,arg->i1, arg->i2);
      }
      else {
        hThreads.push_back(th); 
        TESTPRINT( std::cerr<<"'"<<thName<<"' starts the thread #"<<j<<" with handle="<<th<<":  for(i="<<arg->i1<<";i<="<<arg->i2<<";i++) {...}"<<std::endl<<std::flush );
      }
    }
    TESTPRINT( std::cerr<<std::endl<<std::flush ); 
    for(int j=0; j<(int)hThreads.size(); j++) ui_waitthread(hThreads[j]);    // Wait until threads terminate
    ui_mutex_destroy(mutex);
  }
private:
  T *obj;   
  TR *objR; // will be used to create clone for using at right hand side in expressions inside TExeFun()     
  TExeFun Exe;
  uiMutex mutex;
  std::string thName;
  void setup__(TR *objR,int i1, int i2) { (obj->*Exe)(objR,i1, i2,true); } // call the member function of the obj
  void execute__(TR *objR,int i1, int i2) { 
    ui_mutex_unlock(mutex);  //std::cerr<<"thread for "<<i1<<"<=i<="<<i2<<std::endl<<std::flush;
    (obj->*Exe)(objR,i1, i2,false); 
  }

  static uiThreadFunc EntryPoint(void *A)  {
    ui_mutex_lock(((Args*)A)->mutex);  
    ((Args*)A)->This->execute__(((Args*)A)->objR, ((Args*)A)->i1, ((Args*)A)->i2);
    ui_endthreadex();
    return 0;
  }
};

#if 0 
// the following is the example of how to use class Threads2
//****************************************************************************
// this class implements multithreading for coordinate transformation 
#include "../../MConf/include/C3dMesh.h"
using MConf::Vector3d;
using std::vector;
//class C3dMesh;
class Mag2cyl {
  C3dMesh * mconf;
  vector<Vector3d> * cyl;
  const vector<Vector3d> * mag;
public:
  Mag2cyl(C3dMesh *mc,const vector<Vector3d> *m,vector<Vector3d> *c) : mconf(mc), cyl(c), mag(m) {};

  void calculate(int nThreads) {
    Threads2<Mag2cyl,C3dMesh> threads(this, &Mag2cyl::calculateEx, mconf, "Mag2cyl");  
    int low = 0, upper = int(mag->size()-1);
    threads.run(low, upper, nThreads); 
  }

  template <class mconf> void calculateEx(mconf * mc,int low, int upper, bool setup)
  {
    if(setup) {  // setup part; it allocates arrays if needed
      int size = upper-low+1;
      cyl->resize(size);
    }
    else
    {  // multithreading exe part; it will be called from threads to fill portion [low - upper] of arrays
       // mc is one from mconf clones for this portion
      vector<Vector3d> &cyl = *(this->cyl);
      const vector<Vector3d> &mag = *(this->mag);    
      for(int i=low; i<=upper; i++) {
        cyl[i] = mc->mag2cyl(mag[i]);     
      }
    }
  }

};

// in main programm

  #include "../../MConf/include/C3dMesh.h"
  MConf::C3dMesh mc;
  using MConf::Vector3d;
  using std::vector;
  Vector3d magcoord(.25,1.,2.); 
  vector<Vector3d> m(100,magcoord); // vector with the same magcoord, just for example
  vector<Vector3d> c;               // empty vector, it will be filled
  Mag2cyl m2c(&mc, &m, &c);
  m2c.calculate(numOfThreads);      // transform nagnetic coordinates to cylindrical, c will be filled

#endif


};   //namespace MConf

#undef  TESTPRINT

#endif // MC__THREADS_H
