#
# Try to find SciDb sources and libraries, and determine proper includes
#
# Once done this will define
#  SCIDB_FOUND        - TRUE if SciDb is found
#  SCIDB_SOURCE_DIR   - Where to find SciDb source directory
#  SCIDB_INCLUDE_DIRS - What directories to use for SciDb includes
#  SCIDB_LIBRARIES    - SciDB client libraries
#
# You can define SCIDB_DIR to hint the scidb's location.
#

# Establish a hint
set(TRY_SCIDB_VERSION "14.12")
if (NOT SCIDB_DIR)
    message(STATUS "You can use SCIDB_DIR variable to specify the SciDb sources directory. Using the defaut hint.")
    # assume we have a common source directory for projects
    set(SCIDB_DIR ${CMAKE_SOURCE_DIR}/../scidb-${TRY_SCIDB_VERSION})
endif()
message(STATUS "Using ${SCIDB_DIR} as the SciDb hint directory...")

# find includes and libraries
find_path(SCIDB_SOURCE_DIR include/SciDBAPI.h HINTS ${SCIDB_DIR} DOC "SciDb source directory")
mark_as_advanced(SCIDB_SOURCE_DIR)
find_library(SCIDB_CLIENT_LIB scidbclient HINTS "/opt/scidb/${TRY_SCIDB_VERSION}/lib" DOC "SciDB client library")
mark_as_advanced(SCIDB_CLIENT_LIB)

# standard handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SciDB DEFAULT_MSG SCIDB_SOURCE_DIR SCIDB_CLIENT_LIB)

# set export vars
set(SCIDB_INCLUDE_DIRS ${SCIDB_SOURCE_DIR}/include ${SCIDB_SOURCE_DIR}/src ${SCIDB_SOURCE_DIR}/build/src)
set(SCIDB_LIBRARIES ${SCIDB_CLIENT_LIB})
