include("/home/emdiso/c/HSIM/HSIM-for-Eigen/cmake/CPM.cmake")
CPMAddPackage("NAME;fast_matrix_market;GIT_REPOSITORY;https://github.com/alugowski/fast_matrix_market.git;GIT_TAG;main;GIT_SHALLOW;TRUE")
set(fast_matrix_market_FOUND TRUE)