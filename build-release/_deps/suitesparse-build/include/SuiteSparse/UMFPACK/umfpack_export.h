
#ifndef UMFPACK_EXPORT_H
#define UMFPACK_EXPORT_H

#ifdef UMFPACK_STATIC_DEFINE
#  define UMFPACK_EXPORT
#  define UMFPACK_NO_EXPORT
#else
#  ifndef UMFPACK_EXPORT
#    ifdef umfpack_EXPORTS
        /* We are building this library */
#      define UMFPACK_EXPORT 
#    else
        /* We are using this library */
#      define UMFPACK_EXPORT 
#    endif
#  endif

#  ifndef UMFPACK_NO_EXPORT
#    define UMFPACK_NO_EXPORT 
#  endif
#endif

#ifndef UMFPACK_DEPRECATED
#  define UMFPACK_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef UMFPACK_DEPRECATED_EXPORT
#  define UMFPACK_DEPRECATED_EXPORT UMFPACK_EXPORT UMFPACK_DEPRECATED
#endif

#ifndef UMFPACK_DEPRECATED_NO_EXPORT
#  define UMFPACK_DEPRECATED_NO_EXPORT UMFPACK_NO_EXPORT UMFPACK_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef UMFPACK_NO_DEPRECATED
#    define UMFPACK_NO_DEPRECATED
#  endif
#endif

#endif /* UMFPACK_EXPORT_H */
