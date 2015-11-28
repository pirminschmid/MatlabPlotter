# cmake template for Eigen3
# some modifications re path to make it work on my OSX machine
# 2015-10-20 Pirmin Schmid

project(ProblemSheet3)
cmake_minimum_required(VERSION 2.6.2)

# Tested with g++ and clang
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -O3 -Wextra -Wno-deprecated-declarations -pedantic" )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    # Use -O3 for strong optimization. We use stdc++11. Compile with Wall to get all warnings.
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -O3 -Wextra -Wno-deprecated-declarations -pedantic" )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

# If your compiler already finds and includes Eigen you can safely remove the following two lines
#find_package(Eigen3 REQUIRED)
#include_directories(${EIGEN3_INCLUDE_DIR})

# modified version since the package could not be found on OSX
find_package(Eigen3)

if( NOT EIGEN3_INCLUDE_DIR )
    message( "Package not found. Fallback: trying to use environment variable EIGEN3_INCLUDE_DIR.")
    set( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
endif()

# fallback failed
if( NOT EIGEN3_INCLUDE_DIR )
    message( FATAL_ERROR "Fallback failed. Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
endif()

# here we actually set the path
include_directories(${EIGEN3_INCLUDE_DIR})
message("Path to Eigen3 library found and set for the project.")

### Example for MatlabPlotter
# natural cubic spline interpolation
set(SRCS natcsi.cpp)
add_executable(natcsi ${SRCS}) 

