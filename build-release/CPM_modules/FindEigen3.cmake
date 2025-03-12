include("/home/emdiso/c/HSIM/HSIM-for-Eigen/cmake/CPM.cmake")
CPMAddPackage("NAME;Eigen3;GIT_REPOSITORY;https://gitlab.com/libeigen/eigen.git;GIT_TAG;3.4.0;PATCH_COMMAND;git;restore;.;&&;git;apply;/home/emdiso/c/HSIM/HSIM-for-Eigen/cmake/eigen3_skip_build_demo.patch;OPTIONS;BUILD_TESTING OFF;EIGEN_BUILD_DOC OFF;EIGEN_BUILD_PKGCONFIG OFF")
set(Eigen3_FOUND TRUE)