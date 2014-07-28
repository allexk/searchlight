cmake_minimum_required(VERSION 2.6)
project(searchlight)

# default build type
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting the type of the build to: Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of the build: Release/Debug/RelWithDebInfo/MinSizeRel" FORCE)
endif()
message(STATUS "The type of build is: ${CMAKE_BUILD_TYPE}")

# additional modules and library paths
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# finding libraries
find_package(SciDB REQUIRED)
find_package(OrTools REQUIRED)
find_package(Log4CXX REQUIRED)
find_package(Protobuf REQUIRED)
find_package(GFlags REQUIRED)   # needed by or-tools
find_package(Protobuf REQUIRED) # needed by or-tools
find_package(Coin REQUIRED)     # needed by or-tools

# includes and libraries
include_directories(${SCIDB_INCLUDE_DIRS})
include_directories(${ORTOOLS_INCLUDE_DIRS})
include_directories(${GFLAGS_INCLUDE_DIRS})
include_directories(${LOG4CXX_INCLUDE_DIRS})

# add c++0x compatibility and turn-off STL deprecated warnings (needed for or-tools)
# other warnings are turned off because of SciDb and or-tools
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -Wextra -Wno-strict-aliasing -Wno-sign-compare -Wno-deprecated -Wno-unused-parameter -Wno-ignored-qualifiers")

# if we are 64-bit, define it (again, for or-tools)
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(STATUS "We are on a 64-bit platform...")
    message(STATUS "Setting ARCH_K8 CXX definition for or-tools...")
    add_definitions(-DARCH_K8)
endif()

# PROJECT_ROOT should be defined for SciDb error reporting.
# We will use the parent directory of this project. Assuming
# SciDb sits in the same parent directory, error reporting will
# be fine. Otherwise, there might be problems with error messages
# initiated from SciDb header files. Of course, if we include header
# files from something like /usr/include, we will have problems anyway.
get_filename_component(PROJECT_PARENT_DIR ${CMAKE_SOURCE_DIR} PATH)
message(STATUS "Setting PROJECT_ROOT for SciDb: ${PROJECT_PARENT_DIR}")
add_definitions("-DPROJECT_ROOT=\"\\\"${PROJECT_PARENT_DIR}\\\"\"")

# system library
if (UNIX AND NOT APPLE)
    set(SYS_LIBRARIES rt)
endif()

# Documentation
find_package(Doxygen)
if (DOXYGEN_FOUND)
    configure_file(${CMAKE_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_BINARY_DIR}/doc/Doxyfile @ONLY)
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/doc/Doxyfile
        COMMENT "Generating documentation using Doxygen..."
    )
endif()

# Output directories
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# Since we link with or-tools and we ourselves, which is not
# in the standard path, we need to modify
# RPATHs of the libraries.
set(CMAKE_INSTALL_RPATH "$ORIGIN")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the sources
include_directories(src)
add_subdirectory(src)
