################################################################################
# CMake download helpers
################################################################################

# download external dependencies
include(SCOFDownloadExternal)

################################################################################
# Required dependencies
################################################################################

find_package(OpenMesh QUIET)
if(NOT TARGET OpenMeshCore OR NOT TARGET OpenMeshTools)
    scof_download_openmesh()
endif()


find_package(OpenVolumeMesh QUIET)
if(NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    scof_download_openvolumemesh()
endif()


## Eigen
#find_package(EIGEN3)
#if(NOT EIGEN3_FOUND AND "${CMAKE_PROJECT_NAME}" STREQUAL "SCOF")
#    scof_download_eigen()
#endif()

# Eigen
find_package(EIGEN3)
if(NOT TARGET Eigen3::Eigen)
    scof_download_eigen()
endif()


# GMM
find_package(GMM)
if(NOT GMM_FOUND AND "${CMAKE_PROJECT_NAME}" STREQUAL "SCOF")
	scof_download_gmm()
endif()


# CHOLMOD
find_package(CHOLMOD)


