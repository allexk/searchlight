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
find_package(OrTools REQUIRED)
find_package(GFlags REQUIRED)
find_package(Protobuf REQUIRED) # needed by or-tools
find_package(Coin REQUIRED)  # needed by or-tools

# includes and libraries
include_directories(${ORTOOLS_INCLUDE_DIRS})
include_directories(${GFLAGS_INCLUDE_DIRS})

# add c++0x compatibility and turn-off STL deprecated warnings (needed for or-tools)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wno-deprecated")

# if we are 64-bit, define it (again, for or-tools)
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(STATUS "We are on a 64-bit platform...")
    message(STATUS "Setting ARCH_K8 CXX definition for or-tools...")
    add_definitions(-DARCH_K8)
endif()

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

# the sources
include_directories(src)
add_subdirectory(src)
