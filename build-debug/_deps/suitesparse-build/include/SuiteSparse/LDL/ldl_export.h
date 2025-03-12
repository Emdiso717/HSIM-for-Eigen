
#ifndef LDL_EXPORT_H
#define LDL_EXPORT_H

#ifdef LDL_STATIC_DEFINE
#  define LDL_EXPORT
#  define LDL_NO_EXPORT
#else
#  ifndef LDL_EXPORT
#    ifdef ldl_EXPORTS
        /* We are building this library */
#      define LDL_EXPORT 
#    else
        /* We are using this library */
#      define LDL_EXPORT 
#    endif
#  endif

#  ifndef LDL_NO_EXPORT
#    define LDL_NO_EXPORT 
#  endif
#endif

#ifndef LDL_DEPRECATED
#  define LDL_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef LDL_DEPRECATED_EXPORT
#  define LDL_DEPRECATED_EXPORT LDL_EXPORT LDL_DEPRECATED
#endif

#ifndef LDL_DEPRECATED_NO_EXPORT
#  define LDL_DEPRECATED_NO_EXPORT LDL_NO_EXPORT LDL_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef LDL_NO_DEPRECATED
#    define LDL_NO_DEPRECATED
#  endif
#endif

#endif /* LDL_EXPORT_H */
