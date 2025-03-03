cmake_minimum_required(VERSION 3.5)
project(HSIM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include Eigen
include_directories(${PROJECT_SOURCE_DIR}/eigen-3.4.0)

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