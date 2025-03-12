
#ifndef CHOLMOD_EXPORT_H
#define CHOLMOD_EXPORT_H

#ifdef CHOLMOD_STATIC_DEFINE
#  define CHOLMOD_EXPORT
#  define CHOLMOD_NO_EXPORT
#else
#  ifndef CHOLMOD_EXPORT
#    ifdef cholmod_EXPORTS
        /* We are building this library */
#      define CHOLMOD_EXPORT 
#    else
        /* We are using this library */
#      define CHOLMOD_EXPORT 
#    endif
#  endif

#  ifndef CHOLMOD_NO_EXPORT
#    define CHOLMOD_NO_EXPORT 
#  endif
#endif

#ifndef CHOLMOD_DEPRECATED
#  define CHOLMOD_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CHOLMOD_DEPRECATED_EXPORT
#  define CHOLMOD_DEPRECATED_EXPORT CHOLMOD_EXPORT CHOLMOD_DEPRECATED
#endif

#ifndef CHOLMOD_DEPRECATED_NO_EXPORT
#  define CHOLMOD_DEPRECATED_NO_EXPORT CHOLMOD_NO_EXPORT CHOLMOD_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CHOLMOD_NO_DEPRECATED
#    define CHOLMOD_NO_DEPRECATED
#  endif
#endif

#endif /* CHOLMOD_EXPORT_H */
