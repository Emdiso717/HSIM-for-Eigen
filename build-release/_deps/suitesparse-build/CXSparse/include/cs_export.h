
#ifndef CXSPARSE_EXPORT_H
#define CXSPARSE_EXPORT_H

#ifdef CXSPARSE_STATIC_DEFINE
#  define CXSPARSE_EXPORT
#  define CXSPARSE_NO_EXPORT
#else
#  ifndef CXSPARSE_EXPORT
#    ifdef cxsparse_EXPORTS
        /* We are building this library */
#      define CXSPARSE_EXPORT 
#    else
        /* We are using this library */
#      define CXSPARSE_EXPORT 
#    endif
#  endif

#  ifndef CXSPARSE_NO_EXPORT
#    define CXSPARSE_NO_EXPORT 
#  endif
#endif

#ifndef CXSPARSE_DEPRECATED
#  define CXSPARSE_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CXSPARSE_DEPRECATED_EXPORT
#  define CXSPARSE_DEPRECATED_EXPORT CXSPARSE_EXPORT CXSPARSE_DEPRECATED
#endif

#ifndef CXSPARSE_DEPRECATED_NO_EXPORT
#  define CXSPARSE_DEPRECATED_NO_EXPORT CXSPARSE_NO_EXPORT CXSPARSE_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CXSPARSE_NO_DEPRECATED
#    define CXSPARSE_NO_DEPRECATED
#  endif
#endif

#endif /* CXSPARSE_EXPORT_H */
