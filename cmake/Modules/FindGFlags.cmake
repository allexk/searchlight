#
# Try to find Google's GFLAGS
#
# Once done this will define
#  GFLAGS_FOUND        - TRUE if gflags is found
#  GFLAGS_INCLUDE_DIRS - Where to find gflags include sub-directory
#  GFLAGS_LIBRARIES    - List of libraries when using gflags
#
# You can define GFLAGS_DIR to hint the gflags's location.
#

# use a hint: if we have found or-tools, try to use its dependencies or
# the same root repo dir as our project uses
if (NOT GFLAGS_DIR)
    message(STATUS "You can use GFLAGS_DIR variable to specify the gflags directory. Using the defaut hint.")
    if (ORTOOLS_FOUND)
        set(GFLAGS_DIR ${ORTOOLS_INCLUDE_DIR}/../dependencies/install)
    else()
        set(GFLAGS_DIR ${CMAKE_SOURCE_DIR}/../gflags)
    endif()
endif()
message(STATUS "Using ${GFLAGS_DIR} as the gflags hint directory...")

# find includes and libraries
find_path(GFLAGS_INCLUDE_DIR gflags/gflags.h HINTS ${GFLAGS_DIR}/include DOC "gflags include directory")
find_library(GFLAGS_LIBRARY gflags HINTS ${GFLAGS_DIR}/lib DOC "gflags library")
mark_as_advanced(GFLAGS_LIBRARY GFLAGS_INCLUDE_DIR)

# standard handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GFlags DEFAULT_MSG GFLAGS_LIBRARY GFLAGS_INCLUDE_DIR)

# set export vars
set(GFLAGS_LIBRARIES ${GFLAGS_LIBRARY})
set(GFLAGS_INCLUDE_DIRS ${GFLAGS_INCLUDE_DIR})
