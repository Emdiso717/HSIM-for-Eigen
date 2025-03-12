
#ifndef SPQR_EXPORT_H
#define SPQR_EXPORT_H

#ifdef SPQR_STATIC_DEFINE
#  define SPQR_EXPORT
#  define SPQR_NO_EXPORT
#else
#  ifndef SPQR_EXPORT
#    ifdef spqr_EXPORTS
        /* We are building this library */
#      define SPQR_EXPORT 
#    else
        /* We are using this library */
#      define SPQR_EXPORT 
#    endif
#  endif

#  ifndef SPQR_NO_EXPORT
#    define SPQR_NO_EXPORT 
#  endif
#endif

#ifndef SPQR_DEPRECATED
#  define SPQR_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef SPQR_DEPRECATED_EXPORT
#  define SPQR_DEPRECATED_EXPORT SPQR_EXPORT SPQR_DEPRECATED
#endif

#ifndef SPQR_DEPRECATED_NO_EXPORT
#  define SPQR_DEPRECATED_NO_EXPORT SPQR_NO_EXPORT SPQR_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef SPQR_NO_DEPRECATED
#    define SPQR_NO_DEPRECATED
#  endif
#endif

#endif /* SPQR_EXPORT_H */
