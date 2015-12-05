#
# Try to find Google's or-tools
#
# Once done this will define
#  ORTOOLS_FOUND        - TRUE if or-tools is found
#  ORTOOLS_INCLUDE_DIRS - Where to find ortools include sub-directory
#  ORTOOLS_LIBRARIES    - List of libraries when using or-tools
#
# You can define ORTOOLS_DIR to hint the or-tools's location.
#

# Establish a hint
if (NOT ORTOOLS_DIR)
    message(STATUS "You can use ORTOOLS_DIR variable to specify the or-tools directory. Using the defaut hint.")
    # assume we have a common source directory for projects
    set(ORTOOLS_DIR ${CMAKE_SOURCE_DIR}/../or-tools)
endif()
message(STATUS "Using ${ORTOOLS_DIR} as the or-tools hint directory...")

# find includes and libraries
#
# NOTE: we'd like to go static here. Unfortunately, we cannot do it reliably.
# We have several dynamic libs, like the main UDO and the task ones. Each
# might link with the or-tools static libs. Or-tools uses gflags options,
# so duplicating or-tools object files across multiple libs might result
# in gflags errors, like "something wrong with flag F". The same error might
# also occur if we open/close the same dynamic lib during the execution, and
# we do it with the task libs. Thus, WE Do THE DYNAMIC LINKING HERE.
find_path(ORTOOLS_INCLUDE_DIR constraint_solver/constraint_solver.h HINTS ${ORTOOLS_DIR}/src DOC "or-tools include directory")
find_path(ORTOOLS_INCLUDE_GEN_DIR constraint_solver/model.pb.h HINTS ${ORTOOLS_DIR}/src/gen DOC "or-tools generated files include directory")
find_library(ORTOOLS_CP_LIBRARY constraint_solver HINTS ${ORTOOLS_DIR}/lib DOC "or-tools constraint solver library")
find_library(ORTOOLS_LP_LIBRARY linear_solver HINTS ${ORTOOLS_DIR}/lib DOC "or-tools linear solver library")
find_library(ORTOOLS_BASE_LIBRARY base HINTS ${ORTOOLS_DIR}/lib DOC "or-tools base library")
find_library(ORTOOLS_UTIL_LIBRARY util HINTS ${ORTOOLS_DIR}/lib DOC "or-tools util library")
mark_as_advanced(ORTOOLS_LIBRARY ORTOOLS_INCLUDE_DIR ORTOOLS_INCLUDE_GEN_DIR ORTOOLS_BASE_LIBRARY ORTOOLS_UTIL_LIBRARY ORTOOLS_LP_LIBRARY ORTOOLS_CP_LIBRARY)

# standard handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OrTools DEFAULT_MSG ORTOOLS_CP_LIBRARY ORTOOLS_LP_LIBRARY ORTOOLS_BASE_LIBRARY
    ORTOOLS_UTIL_LIBRARY ORTOOLS_INCLUDE_DIR ORTOOLS_INCLUDE_GEN_DIR)

# set export vars
set(ORTOOLS_LIBRARIES ${ORTOOLS_CP_LIBRARY} ${ORTOOLS_BASE_LIBRARY} ${ORTOOLS_UTIL_LIBRARY} ${ORTOOLS_LP_LIBRARY})
set(ORTOOLS_INCLUDE_DIRS ${ORTOOLS_INCLUDE_DIR} ${ORTOOLS_INCLUDE_GEN_DIR})
