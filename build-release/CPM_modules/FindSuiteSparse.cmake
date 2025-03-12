include("/home/emdiso/c/HSIM/HSIM-for-Eigen/cmake/CPM.cmake")
CPMAddPackage("NAME;SuiteSparse;GIT_REPOSITORY;https://github.com/sergiud/SuiteSparse.git;GIT_TAG;5.13.0-cmake.4;OPTIONS;WITH_CUDA OFF")
set(SuiteSparse_FOUND TRUE)