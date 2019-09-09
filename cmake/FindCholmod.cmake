# Cholmod lib usually requires linking to a blas and lapack library.
# It is up to the user of this module to find a BLAS and link to it.

if (CHOLMOD_INCLUDES AND CHOLMOD_LIBRARIES)
  set(CHOLMOD_FIND_QUIETLY TRUE)
endif (CHOLMOD_INCLUDES AND CHOLMOD_LIBRARIES)

find_path(CHOLMOD_INCLUDES
  NAMES
  cholmod.h
  PATHS
  $ENV{CHOLMODDIR}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  suitesparse
  ufsparse
)

find_library(CHOLMOD_LIBRARIES cholmod PATHS $ENV{CHOLMODDIR} ${LIB_INSTALL_DIR})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CHOLMOD DEFAULT_MSG
                                  CHOLMOD_INCLUDES CHOLMOD_LIBRARIES)

mark_as_advanced(CHOLMOD_INCLUDES CHOLMOD_LIBRARIES)