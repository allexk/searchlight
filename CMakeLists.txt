cmake_minimum_required(VERSION 2.6)
project(searchlight)

# default build type
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting the type of the build to: Release")
    # use FORCE for caching since CMake creates an initial empty value by default, which prevents non-force caching
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
find_package(Boost REQUIRED)    # we need only headers because of SciDb
find_package(GFlags REQUIRED)   # needed by or-tools
find_package(Coin REQUIRED)     # needed by or-tools
find_package(FFTW REQUIRED)

# includes and libraries
include_directories(${SCIDB_INCLUDE_DIRS})
include_directories(${ORTOOLS_INCLUDE_DIRS})
include_directories(${GFLAGS_INCLUDE_DIRS})
include_directories(${LOG4CXX_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PROTOBUF_INCLUDE_DIRS})
include_directories(${FFTW_INCLUDE_DIRS})

# add c++0x compatibility and turn-off STL deprecated warnings (needed for or-tools)
# other warnings are turned off because of SciDb and or-tools
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Wno-strict-aliasing -Wno-sign-compare -Wno-deprecated -Wno-unused-parameter -Wno-ignored-qualifiers -Wno-unused-local-typedefs")

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

# Libraries install dir
if (NOT SL_LIB_DIR)
    message(STATUS "Libraries will be installed into: lib/scidb/plugins")
    set(SL_LIB_DIR "lib/scidb/plugins" CACHE STRING "Relatieve path to install the SL libraries to")
endif()

# Since we link with or-tools and we ourselves, which is not
# in the standard path, we need to modify
# RPATHs of the libraries.
set(CMAKE_INSTALL_RPATH "$ORIGIN:$ORIGIN/sl-aux")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

# Determine libraries to be installed in the plugin directory along with searchlight
foreach (lib IN LISTS ORTOOLS_LIBRARIES GFLAGS_LIBRARIES COIN_LIBRARIES FFTW_LIBRARIES)
    get_filename_component(lib_path ${lib} PATH)
    get_filename_component(lib_name ${lib} NAME)
    get_filename_component(lib_ext ${lib} EXT)
    if (NOT lib_ext STREQUAL ".a") # Don't install static libs
        list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES ${lib_path} is_system_dir)
        if (is_system_dir EQUAL -1)
            # A little bit of cheating: since these libraries are not handled
            # by CMake, we need to copy the corresponding symlink .so.* files
            # for the dynamic linker. We do it by simple GLOBing.
            file(GLOB lib_files_all "${lib}*")
            install(FILES ${lib_files_all} DESTINATION ${SL_LIB_DIR}/sl-aux)
            message(STATUS "Library ${lib_name} will be installed to the SciDB plugins directory...")
        endif()
    endif()
endforeach()

# the sources
include_directories(src)
add_subdirectory(src)
