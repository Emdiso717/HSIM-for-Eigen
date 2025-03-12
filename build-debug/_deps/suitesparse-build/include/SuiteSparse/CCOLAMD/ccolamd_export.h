
#ifndef CCOLAMD_EXPORT_H
#define CCOLAMD_EXPORT_H

#ifdef CCOLAMD_STATIC_DEFINE
#  define CCOLAMD_EXPORT
#  define CCOLAMD_NO_EXPORT
#else
#  ifndef CCOLAMD_EXPORT
#    ifdef ccolamd_EXPORTS
        /* We are building this library */
#      define CCOLAMD_EXPORT 
#    else
        /* We are using this library */
#      define CCOLAMD_EXPORT 
#    endif
#  endif

#  ifndef CCOLAMD_NO_EXPORT
#    define CCOLAMD_NO_EXPORT 
#  endif
#endif

#ifndef CCOLAMD_DEPRECATED
#  define CCOLAMD_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CCOLAMD_DEPRECATED_EXPORT
#  define CCOLAMD_DEPRECATED_EXPORT CCOLAMD_EXPORT CCOLAMD_DEPRECATED
#endif

#ifndef CCOLAMD_DEPRECATED_NO_EXPORT
#  define CCOLAMD_DEPRECATED_NO_EXPORT CCOLAMD_NO_EXPORT CCOLAMD_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CCOLAMD_NO_DEPRECATED
#    define CCOLAMD_NO_DEPRECATED
#  endif
#endif

#endif /* CCOLAMD_EXPORT_H */
