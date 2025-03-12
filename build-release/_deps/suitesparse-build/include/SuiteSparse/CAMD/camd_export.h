
#ifndef CAMD_EXPORT_H
#define CAMD_EXPORT_H

#ifdef CAMD_STATIC_DEFINE
#  define CAMD_EXPORT
#  define CAMD_NO_EXPORT
#else
#  ifndef CAMD_EXPORT
#    ifdef camd_EXPORTS
        /* We are building this library */
#      define CAMD_EXPORT 
#    else
        /* We are using this library */
#      define CAMD_EXPORT 
#    endif
#  endif

#  ifndef CAMD_NO_EXPORT
#    define CAMD_NO_EXPORT 
#  endif
#endif

#ifndef CAMD_DEPRECATED
#  define CAMD_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CAMD_DEPRECATED_EXPORT
#  define CAMD_DEPRECATED_EXPORT CAMD_EXPORT CAMD_DEPRECATED
#endif

#ifndef CAMD_DEPRECATED_NO_EXPORT
#  define CAMD_DEPRECATED_NO_EXPORT CAMD_NO_EXPORT CAMD_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CAMD_NO_DEPRECATED
#    define CAMD_NO_DEPRECATED
#  endif
#endif

#endif /* CAMD_EXPORT_H */
