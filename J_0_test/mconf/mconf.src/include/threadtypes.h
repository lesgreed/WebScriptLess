#ifndef MC__THREADTYPES_H
#define MC__THREADTYPES_H

#ifdef _WIN32
  #include <windows.h>       // Windows threading
  #include <process.h>

  typedef HANDLE  uiThread;
  typedef HANDLE  uiMutex;
  #define uiThreadFunc  unsigned __stdcall

  
#else                        // POSIX threading
  #include <pthread.h>

  typedef pthread_t        uiThread;
  typedef pthread_mutex_t* uiMutex;
  #define uiThreadFunc  void *

#endif

#endif // MC__THREADTYPES_H
