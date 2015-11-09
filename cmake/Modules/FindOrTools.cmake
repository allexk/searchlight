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

# find includes and libraries (prefer static)
find_path(ORTOOLS_INCLUDE_DIR constraint_solver/constraint_solver.h HINTS ${ORTOOLS_DIR}/src DOC "or-tools include directory")
find_path(ORTOOLS_INCLUDE_GEN_DIR constraint_solver/model.pb.h HINTS ${ORTOOLS_DIR}/src/gen DOC "or-tools generated files include directory")
find_library(ORTOOLS_CP_LIBRARY libconstraint_solver.a constraint_solver HINTS ${ORTOOLS_DIR}/lib DOC "or-tools constraint solver library")
find_library(ORTOOLS_LP_LIBRARY liblinear_solver.a linear_solver HINTS ${ORTOOLS_DIR}/lib DOC "or-tools linear solver library")
find_library(ORTOOLS_BASE_LIBRARY libbase.a base HINTS ${ORTOOLS_DIR}/lib DOC "or-tools base library")
find_library(ORTOOLS_UTIL_LIBRARY libutil.a util HINTS ${ORTOOLS_DIR}/lib DOC "or-tools util library")
mark_as_advanced(ORTOOLS_LIBRARY ORTOOLS_INCLUDE_DIR ORTOOLS_INCLUDE_GEN_DIR ORTOOLS_BASE_LIBRARY ORTOOLS_UTIL_LIBRARY ORTOOLS_LP_LIBRARY)

# standard handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OrTools DEFAULT_MSG ORTOOLS_CP_LIBRARY ORTOOLS_LP_LIBRARY ORTOOLS_BASE_LIBRARY
    ORTOOLS_UTIL_LIBRARY ORTOOLS_INCLUDE_DIR ORTOOLS_INCLUDE_GEN_DIR)

# set export vars
set(ORTOOLS_LIBRARIES ${ORTOOLS_CP_LIBRARY} ${ORTOOLS_BASE_LIBRARY} ${ORTOOLS_UTIL_LIBRARY} ${ORTOOLS_LP_LIBRARY})
set(ORTOOLS_INCLUDE_DIRS ${ORTOOLS_INCLUDE_DIR} ${ORTOOLS_INCLUDE_GEN_DIR})
