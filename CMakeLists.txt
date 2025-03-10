cmake_minimum_required(VERSION 3.5)
project(HSIM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ffunction-sections -fdata-sections")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

# Include Eigen
include_directories(${PROJECT_SOURCE_DIR}/eigen-3.4.0)
include_directories(${PROJECT_SOURCE_DIR}/Spectra)

# Add source files
set(SOURCES
    main.cpp
    HSIM.cpp
    laplace_3d.cpp
)

# Add header files
set(HEADERS
    HSIM.h
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})