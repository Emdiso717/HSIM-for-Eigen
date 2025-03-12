include("/home/emdiso/c/HSIM/HSIM-for-Eigen/cmake/CPM.cmake")
CPMAddPackage("NAME;spectra;GIT_REPOSITORY;https://github.com/yixuan/spectra.git;GIT_TAG;v1.0.1;PATCH_COMMAND;git;restore;.;&&;git;apply;/home/emdiso/c/HSIM/HSIM-for-Eigen/cmake/spectra_fix_eigen.patch")
set(spectra_FOUND TRUE)