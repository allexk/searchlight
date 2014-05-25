#
# Try to find log4cxx
#
# Once done this will define
#  LOG4CXX_FOUND        - TRUE if log4cxx is found
#  LOG4CXX_INCLUDE_DIRS - Where to find log4cxx headers
#  LOG4CXX_LIBRARIES    - List of libraries when using log4cxx
#

# find includes and libraries
find_path(LOG4CXX_INCLUDE_DIR log4cxx/log4cxx.h DOC "log4cxx include directory")
find_library(LOG4CXX_LIBRARY log4cxx DOC "log4cxx library")
mark_as_advanced(LOG4CXX_LIBRARY LOG4CXX_INCLUDE_DIR)

# standard handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Log4CXX DEFAULT_MSG LOG4CXX_LIBRARY LOG4CXX_INCLUDE_DIR)

# set export vars
set(LOG4CXX_LIBRARIES ${LOG4CXX_LIBRARY})
set(LOG4CXX_INCLUDE_DIRS ${LOG4CXX_INCLUDE_DIR})
