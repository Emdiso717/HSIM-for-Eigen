aux_source_directory(. source)
file(GLOB header *.h)

add_library(eigen_solver SHARED ${source})
target_link_libraries(eigen_solver
  Eigen3_WITH_BLAS
  Spectra
  SuiteSparse::CHOLMOD
  CLI11::CLI11
)


target_compile_features(eigen_solver PUBLIC cxx_std_20)