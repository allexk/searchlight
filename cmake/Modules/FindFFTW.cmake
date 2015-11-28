#
# Try to find FFTW (Fast Fourier Transform library)
#
# Once done this will define
#  FFTW_FOUND        - TRUE if FFTW is found
#  FFTW_INCLUDE_DIRS - Where to find FFTW include sub-directory
#  FFTW_LIBRARIES    - List of libraries when using FFTW
#

# find includes and libraries (prefer static)
find_path(FFTW_INCLUDE_DIR fftw3.h DOC "FFTW include directory")
find_library(FFTW_LIBRARY libfftw3.a fftw3 DOC "FFTW library")
mark_as_advanced(FFTW_LIBRARY FFTW_INCLUDE_DIR)

# standard handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARY FFTW_INCLUDE_DIR)

# set export vars
set(FFTW_LIBRARIES ${FFTW_LIBRARY})
set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
