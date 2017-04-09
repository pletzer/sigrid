###################################################
#
# CMakeLists.txt for simplexgrid
#
###################################################

# Project information
project(sigrid)
set(VERSION_MAJOR "0")
set(VERSION_MINOR "1")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

#set(DEBUG_CMAKE TRUE)

###################################################
#
# Create build type Profile
#
###################################################

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS_PROFILE "-g -pg -O2")
  set(CMAKE_C_FLAGS_PROFILE "-g -pg -O2")
  #set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  #set(BUILD_SHARED_LIBRARIES OFF)
  #set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

###################################################

# Required version
cmake_minimum_required(VERSION 2.8.8)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(PythonInterp REQUIRED)

enable_testing()

add_subdirectory(cpp)
add_subdirectory(utils)
add_subdirectory(tests)
add_subdirectory(examples)

