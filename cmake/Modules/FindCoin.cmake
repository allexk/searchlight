#
# Try to find coin
# TODO: for now I need only clp, so that is what I am trying to find. Maybe add COMPONENTS later?
#
# Once done this will define
#  COIN_FOUND        - TRUE if coin is found
#  COIN_INCLUDE_DIRS - Where to find coin include sub-directory
#  COIN_LIBRARIES    - List of libraries when using coin
#
# You can define COIN_DIR to hint the coin's location.
#

# use a hint: if we have found or-tools, try to use its dependencies or
# the same root repo dir as our project uses
if (NOT COIN_DIR)
    message(STATUS "You can use COIN_DIR variable to specify the coin directory. Using the defaut hint.")
    if (ORTOOLS_FOUND)
        set(COIN_DIR ${ORTOOLS_INCLUDE_DIR}/../dependencies/install)
    else()
        set(COIN_DIR ${CMAKE_SOURCE_DIR}/../coin)
    endif()
endif()
message(STATUS "Using ${COIN_DIR} as the coin hint directory...")

# find includes and libraries (prefer static)
find_path(COIN_INCLUDE_DIR coin/CbcConfig.h HINTS ${COIN_DIR}/include DOC "coin include directory")
find_library(COIN_UTIL_LIBRARY libCoinUtils.a CoinUtils HINTS ${COIN_DIR}/lib DOC "coin utils library")
find_library(COIN_LP_LIBRARY libClp.a Clp HINTS ${COIN_DIR}/lib DOC "coin clp library")
mark_as_advanced(COIN_UTIL_LIBRARY COIN_LP_LIBRARY COIN_INCLUDE_DIR)

# standard handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Coin DEFAULT_MSG COIN_UTIL_LIBRARY COIN_LP_LIBRARY COIN_INCLUDE_DIR)

# set export vars
set(COIN_LIBRARIES ${COIN_UTIL_LIBRARY} ${COIN_LP_LIBRARY})
set(COIN_INCLUDE_DIRS ${COIN_INCLUDE_DIR})
