#ifndef Export_h
#define Export_h

#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_DLL
    #ifdef __GNUC__
      #define EXPORT __attribute__ ((dllexport))
    #else
      #define EXPORT __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #else
    #ifdef __GNUC__
      #define EXPORT __attribute__ ((dllimport))
    #else
      #define EXPORT __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #endif
  #define HIDE
#else
  #if __GNUC__ >= 4
    #define EXPORT __attribute__ ((visibility ("default")))
    #define HIDE  __attribute__ ((visibility ("hidden")))
  #else
    #define EXPORT
    #define HIDE
  #endif
#endif

#endif
