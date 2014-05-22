#
# Try to find SciDb sources and determine proper includes
#
# Once done this will define
#  SCIDB_FOUND        - TRUE if SciDb is found
#  SCIDB_SOURCE_DIR   - Where to find SciDb source directory
#  SCIDB_INCLUDE_DIRS - What directories to use for SciDb includes
#
# You can define SCIDB_DIR to hint the scidb's location.
#

# Establish a hint
if (NOT SCIDB_DIR)
    message(STATUS "You can use SCIDB_DIR variable to specify the SciDb sources directory. Using the defaut hint.")
    # assume we have a common source directory for projects
    set(SCIDB_DIR ${CMAKE_SOURCE_DIR}/../scidb-13.6)
endif()
message(STATUS "Using ${SCIDB_DIR} as the SciDb hint directory...")

# find includes and libraries
find_path(SCIDB_SOURCE_DIR include/SciDBAPI.h HINTS ${SCIDB_DIR} DOC "SciDb source directory")
mark_as_advanced(SCIDB_SOURCE_DIR)

# standard handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SciDB DEFAULT_MSG SCIDB_SOURCE_DIR)

# set export vars
set(SCIDB_INCLUDE_DIRS ${SCIDB_SOURCE_DIR}/include ${SCIDB_SOURCE_DIR}/src)
