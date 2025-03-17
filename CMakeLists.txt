cmake_minimum_required(VERSION 3.5)
project(HSIM)

include(cmake/CPM.cmake)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/;"${CMAKE_MODULE_PATH})

message("cmake build type: ${CMAKE_BUILD_TYPE}")
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ")
endif()
#enable all instructions
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wno-c++11-narrowing -Wno-unused-variable -Wno-unused-parameter -Wno-array-bounds -Wno-dangling-reference -Wno-tautological-compare")

set(INCLUDES
  suitesparse
  openmp
  eigen 
  spectra
  fast_matrix_market
  cli11
)

foreach(INCLUDE IN LISTS INCLUDES)
    include(cmake/include_${INCLUDE}.cmake)
    message("Load ${INCLUDE} done.")
endforeach()


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
include_directories(src)
add_subdirectory(src)
add_executable(HSIM main.cpp)
add_executable(gen_U gen_U.cpp)
add_executable(testing testing.cpp)
target_link_libraries(HSIM eigen_solver)
target_link_libraries(gen_U eigen_solver)
target_link_libraries(testing eigen_solver)
target_link_libraries(gen_U ${VTK_LIBRARIES})

