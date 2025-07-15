#ifndef C_mathf_
#define C_mathf_
#include <math.h>
#include <time.h>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>

//#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#if defined( _WIN32 )
#else  //#if (defined( __unix__ ) || (defined(linux)) || defined(__sun) || defined(__APPLE__) )
  #include <sys/times.h>
  #include <unistd.h>
#endif

namespace mathf {

#if defined( _WIN32 )
 static double wall_clock() {return double(clock())/CLOCKS_PER_SEC; }
#else 
  static double wall_clock() {struct tms buf; return double(times(&buf))/sysconf(_SC_CLK_TCK); }
#endif

  const double pi    = 3.14159265358979323846;
  const double twopi = 2*pi;
  const double degree = pi/180;

  template <class T> T cube  (T x)   { return( x*x*x );  };
  template <class T> T square  (T x)   { return( x*x );  };
  template <class T> T mmin (T x, T y) { return (x<y)?x:y; };
  template <class T> T mmax (T x, T y) { return (x>y)?x:y; };
  template <class T> T clamp (T x, T a, T b) { x=mmax(x,a); return mmin(x,b); };
  template <class T> void swap (T &v1, T &v2) {T w = v1;  v1 = v2;  v2 = w;};
  template <class T> int sign(T x) { return x>=0?1:-1; };
  // remove twopi and reduce x to range -pi..pi
  template <class T> double modpi(T x) { x = (x>=0)?fmod(double(x),twopi):(fmod(double(x),twopi)+twopi); return (x>pi)?(x-twopi):x; }; 
  template <class T> double modPeriod(T x,double mPeriod ) { return (x>=0)?fmod(double(x),mPeriod):(fmod(double(x),mPeriod)+mPeriod);}; // remove  period
  // C++ version std::string style "itoa":
  template <class T> std::string itoa(T value, int base) {
    enum { kMaxDigits = 35 };
    std::string buf;
    buf.reserve( kMaxDigits ); // Pre-allocate enough space.
    if (base<2||base>16) return buf;  // check that the base if valid
    int quotient = value;
    do {  // Translating number to string with base:
      buf += "0123456789abcdef"[abs(quotient%base)];
      quotient /= base;
    } while ( quotient );
    if(value<0&&base==10) buf += '-';  // Append the negative sign for base 10
    std::reverse(buf.begin(), buf.end());
    return buf;
  };
  //******************************************************************
  // at return: path with ending /, ext with the first symbol .
 static void getPathNameExt(std::string &fullname,  std::string &path, std::string &name, std::string &ext)
  {
    std::size_t  sz1 = fullname.size()-1;
    for (std::size_t  i=0; i<fullname.size(); i++) if(fullname[i]=='\\') fullname[i] = '/'; 
    std::size_t  j = fullname.find_last_of("/");
    if(j==std::string::npos) {
      name = fullname;
      path = "";
    }
    else {
      name = (j<sz1)?fullname.substr(j+1):"";         // filename;
      path = (j<sz1)?fullname.substr(0,j+1):fullname; // file path
    }
    std::size_t  i = name.find_last_of(".");
    ext = (i!=std::string::npos)?name.substr(i):"";   // file extention
    if(i!=std::string::npos) name.resize(i);          // discard file extention
  }   

//******************************************************************************
// Checks to see if a directory exists. Note: This method only checks the
// existence of the full path AND if path leaf is a dir.
// 
// @return  >0 if dir exists AND is a dir,
//           0 if dir does not exist OR exists but not a dir,
//          <0 if an error occurred (errno is also set)
static int isDir(const char* path)
{
  struct stat info;
  int statRC = stat( path, &info );
  if(statRC != 0 ){
    if(errno == ENOENT)  { return 0; } // something along the path does not exist
    if(errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
    return -1;
  }
  int k = (info.st_mode&S_IFDIR)?1:0;
  return k;
}

//******************************************************************************
// Create folder subdir in user %TEMP% forder
// @return full path to the %TEMP%\subdir with a trailing slash
static std::string createTempSubdir(std::string &subdir)
{
#ifdef _WIN32
  const char *key[] = {"TMP", "TEMP", "USERPROFILE", NULL};
#else
  const char *key[] = {"TMPDIR", "TMP", "TEMP", "TEMPDIR","USERPROFILE", NULL};
#endif
  char *val;
  for(int i=0; key[i]!=NULL; i++) {
    val = getenv(key[i]);
    if(val) break;
  }
#ifdef _WIN32
  std::string path = (val==NULL)?std::string(""):std::string(val);
  if(path[path.size()-1]=='\\') path.resize(path.size()-1);
  path += "\\";
  path += subdir;
  std::string mkdir("mkdir "); 
  mkdir += path;
  int k = isDir(path.c_str());
  path += "\\";
#else
  std::string path = (val==NULL)?std::string("/tmp"):std::string(val);
  if(path[path.size()-1]=='/') path.resize(path.size()-1);
  path += "/";
  path += subdir;
  std::string mkdir("mkdir -p "); 
  mkdir += path;
  int k = isDir(path.c_str());
  path += "/";
#endif
  if(k==0) system(mkdir.c_str());
  return path;
}

//******************************************************************************
// Create unique name for temporary file in current working directory
// name is the guess name
// name is the full filename with the path or name for file in the current working directory
// @return name for temporary file
static std::string createUniqueFilename(const std::string &name)
{
  std::string tmpname = name;
  for(int n=1;;++n) {
    FILE* fp = ::fopen(tmpname.c_str(),"r");   // try to open
    if(fp==0) {
      int k = isDir(tmpname.c_str());
      if(k<=0) break; // break if not found
    }
    if(fp) ::fclose(fp);
    std::ostringstream os;
    os<<name<<"_"<<n;
    tmpname=os.str();
  }
  return tmpname;
}

//******************************************************************************
// Create unique subdir at directory path.
// path is the directory path
// subdir is the guess name for subdir
// @return path to temporary subdir 
static std::string createUniqueSubdir(const std::string &path,const std::string &subdir)
{
  std::string subdir1 = subdir;
  std::string tmpname;
  for(int n=1;;++n) {
    tmpname = path+subdir1;
    int err = isDir(tmpname.c_str());
    if(!(err>0)) break; // break if not found
    std::ostringstream os;
    os<<subdir<<"_"<<n;
    subdir1=os.str();
  }
#ifdef _WIN32
  std::string mkdir("mkdir "); 
  mkdir += tmpname;
  tmpname += "\\";
#else
  std::string mkdir("mkdir -p "); 
  mkdir += tmpname;
  tmpname += "/";
#endif
  system(mkdir.c_str());
  return tmpname;
}

//Create directory for output files.
// @return the path of the directory created or empty std::string  
static std::string createOutputDir(std::string &outputFile)
{
 // std::string path = getPath(outputFile); //get absolute path from the full filename
  std::string path, name, ext;
  mathf::getPathNameExt(outputFile,path,name,ext);
  if(path[path.size()-1 ]=='/') path.resize(path.size()-1);
#ifdef _WIN32
  std::replace(path.begin(),path.end(), '/', '\\');
//    _mkdir(path.c_str());
  std::string mkdir("mkdir "); 
  mkdir += path;
  int k = isDir(path.c_str());
#else
  std::string mkdir("mkdir -p "); 
  mkdir += path;
  int k = isDir(path.c_str());
  //int status = mkdir("/home/fsk/TMP", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
  if(k==0) system(mkdir.c_str());
  if (isDir(path.c_str())<=0) path.clear();
  return path;
}
  
  
/// StopWatch provides functions for timings
///                      2006.12.28  turkin
///                      2007.09.12  
class StopWatch {
  clock_t tStart;
  clock_t tSplit;
  clock_t tElapsed;
  bool started;
  bool splited;
  bool suspended;
  bool doPrint;
public:
  /// Constructor
  StopWatch() {
    reset();
    doPrint=true;
  };
  /// Resets the stopwatch
  void setPrint(bool flag) {
    doPrint=flag;
  };
  /// Resets the stopwatch
  void reset() {
    started=splited=suspended=false; 
    tStart=tSplit=tElapsed=0;
  };
  /// Start the stopwatch
  void start()   {
    if(started) return;
    reset();
    started = true; 
    tStart = clock(); 
  };
  /// Restart the stopwatch
  void restart() {
    stop();
    start();
  };
  /// Stop the stopwatch
  void stop()    {
    if(!started) return; 
    clock_t tStop  = clock(); 
    if(suspended) tStop =tStart;
    tElapsed +=tStop-tStart; 
    started=false; 
  };
  /// Split the time
  void split()   {
    if(!started || splited) return;
    if(suspended) 
      tSplit = tElapsed;
    else
      tSplit = tElapsed + clock()-tStart; 
    splited = true;
  };
  /// Remove a split
  void unsplit() {
    if(!started || !splited) return;
    tSplit = 0;
    splited = false;
  };
  /// Suspend the stopwatch for later resumption
  void suspend() {
    if(!started || suspended) return; 
    tElapsed += clock()-tStart;
    suspended = true;
  };
  /// Resume the stopwatch after a suspend
  void resume()  {
    if(!started || !suspended) return; 
    tStart  = clock();
    suspended = false;
  };
  /// Get the time (in sec) on the stopwatch
  double getTime() {
    if(!started || suspended) return double(tElapsed)/CLOCKS_PER_SEC; 
    clock_t t = tElapsed + (clock()-tStart);
    return double(t)/CLOCKS_PER_SEC; 
  };
  void printTime(const char * comment,bool reStart=false) {
    if(doPrint) 
      std::cout<<comment<<" "<<getTime()<<"s"<<std::endl;
    if(reStart) restart();
  };
  /// Get the split time (in sec) on the stopwatch
  double getSplitTime() {
    if(!splited || !started) return 0; 
    return double(tSplit)/CLOCKS_PER_SEC; 
  };
};

/*
 /// remove initial and trailing blanks,LF,CR
  void strtrim(char* f)  
  {
    char *p,*p2;
    int i=strlen(f);
    for(p=f; p!=f+i-1; p++) if(*p!=' ') break;
    i=strlen(p);
    for(p2=p+i-1; p2!=p; p2--) 
      if(*p2==' '||*p2==0x0D||*p2==0x0A)
        continue;
      else 
      {*(p2+1)=0; break; }
    memmove(f,p,p2-p+2);
  };
*/

  /// The class provides vector of pointers to simplify dealocation.
  /// @parm pT is the pointer to some object, memory for which was allocated by new
  /// The memory is deallocated by clear() or automatically by destructor when this class is deleted.
  template < typename pT > class vectorOfPtr: public std::vector< pT >
  {
    protected:
      pT This(int i) {return std::vector<pT>::at(i);} // get element at position i
    public:
    ~vectorOfPtr() {
      clear();
    }
    void clear() {
      if(std::vector<pT>::empty()) return;
      for (unsigned int i = 0; i < std::vector<pT>::size(); i++) {
         delete std::vector<pT>::at(i); 
      }
      std::vector< pT >::clear();
    }
    void resize(int n) {
      clear();
      std::vector<pT>::resize(n);
    }
  };

// wrap an allocation by new to ensure destruction when control leaves a scope.
  template < typename T > class autoPtr  {
    T *data;
    autoPtr(int n) {data = new T[n];}
  public:
    autoPtr() {data = 0;}
    autoPtr(T* d) {data = d;}
    ~autoPtr() { delete[] data; }
    T* get() const {return data;}
    T *operator->() const {return (get());}
    T& operator*() const {return (*get());}
    autoPtr &operator = ( autoPtr &v ) {
      if( this != &v) {
        delete[] data;
        data = v.data;
        v.data = 0;
      }
      return( *this );
    }
  };

#define MATHFNEWARRAY(T, y, n) mathf::autoPtr<T> alloc__##y(new T[n]); T *y = alloc__##y.get();
#define MATHFNEWOBJECT(T, y) mathf::autoPtr<T> alloc__##y(new T); T *y = alloc__##y.get();

  //MATHFNEWARRAY(double, y, 10);
  //MATHFNEWOBJECT(MConf::C3dMesh, mc);


//_WIN32 is Defined for applications for Win32 and Win64. Always defined.
#if defined( _WIN32 ) 
    #ifdef _MSC_VER
      #if _MSC_VER > 1200       // Visual Studio higher than 6
/* class auto_ptr
 * - improved standard conforming implementation
 wrap an allocation by new to ensure destruction when control leaves a scope.

The following code example is taken from the book
The C++ Standard Library - A Tutorial and Reference
by Nicolai M. Josuttis, Addison-Wesley, 1999
© Copyright Nicolai M. Josuttis 1999

 */
    // auxiliary type to enable copies and assignments (now global)
    template<class Y>
    struct auto_ptr_ref {
        Y* yp;
        auto_ptr_ref (Y* rhs)
         : yp(rhs) {
        }
    };

    template<class T>
    class auto_ptr {
      private:
        T* ap;    // refers to the actual owned object (if any)
      public:
        typedef T element_type;
        // constructor
        explicit auto_ptr (T* ptr = 0) throw() : ap(ptr) { }
        // copy constructors (with implicit conversion)
        // - note: nonconstant parameter
        auto_ptr (auto_ptr& rhs) throw() : ap(rhs.release()) { }
        template<class Y>
        auto_ptr (auto_ptr<Y>& rhs) throw() : ap(rhs.release()) { }
        // assignments (with implicit conversion)
        // - note: nonconstant parameter
        auto_ptr& operator= (auto_ptr& rhs) throw() {
            reset(rhs.release());
            return *this;
        }
        template<class Y>
        auto_ptr& operator= (auto_ptr<Y>& rhs) throw() {
            reset(rhs.release());
            return *this;
        }
        // destructor
        ~auto_ptr() throw() { delete ap; }
        // value access
        T* get() const throw() { return ap; }
        T& operator*() const throw() { return *ap; }
        T* operator->() const throw() {  return ap; }
        // release ownership
        T* release() throw() {
            T* tmp(ap);
            ap = 0;
            return tmp;
        }
        // reset value
        void reset (T* ptr=0) throw() {
            if (ap != ptr) {
                delete ap;
                ap = ptr;
            }
        }
        // special conversions with auxiliary type to enable copies and assignments
        auto_ptr(auto_ptr_ref<T> rhs) throw() : ap(rhs.yp) {  }
        auto_ptr& operator= (auto_ptr_ref<T> rhs) throw() {  // new
             reset(rhs.yp);
             return *this;
        }
        template<class Y>
        operator auto_ptr_ref<Y>() throw() { return auto_ptr_ref<Y>(release()); }
        template<class Y>
        operator auto_ptr<Y>() throw() { return auto_ptr<Y>(release()); }
    };

#endif
#endif
#endif

};



#endif
