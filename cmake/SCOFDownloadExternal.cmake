################################################################################
include(DownloadProject)

if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(SCOF_DOWNLOAD_MISSINGS_DEPS_DEFAULT ON)
else()
    set(SCOF_DOWNLOAD_MISSINGS_DEPS_DEFAULT OFF)
endif()
option(SCOF_DOWNLOAD_MISSINGS_DEPS "Download missing SCOF dependencies" ${SCOF_DOWNLOAD_MISSINGS_DEPS_DEFAULT})

# Shortcut function
function(scof_download_project name)
    if (NOT ${SCOF_DOWNLOAD_MISSINGS_DEPS})
        message(WARNING "SCOF: not downloading missing dependency '${name}', because SCOF_DOWNLOAD_MISSINGS_DEPS is not set")
        return()
    endif()
    message(STATUS "SCOF: downloading missing dependency '${name}'")

    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${SCOF_EXTERNAL}/${name}
        DOWNLOAD_DIR ${SCOF_EXTERNAL}/.cache/${name}
        ${ARGN}
    )
endfunction()

################################################################################
## OpenMesh
function(scof_download_openmesh)
    scof_download_project(OpenMesh
        GIT_REPOSITORY       https://www.graphics.rwth-aachen.de:9000/OpenMesh/OpenMesh.git
        GIT_TAG              master
        )
    if (${SCOF_DOWNLOAD_MISSINGS_DEPS})
        add_subdirectory(${SCOF_EXTERNAL}/OpenMesh)
    endif()
endfunction()


## OpenVolumeMesh
function(scof_download_openvolumemesh)
    scof_download_project(OpenVolumeMesh
        GIT_REPOSITORY       https://www.graphics.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh.git
        GIT_TAG              v3.2.1
        )
    if (${SCOF_DOWNLOAD_MISSINGS_DEPS})
        add_subdirectory(${SCOF_EXTERNAL}/OpenVolumeMesh)
    endif()
endfunction()


## Eigen
function(scof_download_eigen)
    scof_download_project(eigen
        URL           http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
        URL_MD5       f2a417d083fe8ca4b8ed2bc613d20f07
    )
    if (${SCOF_DOWNLOAD_MISSINGS_DEPS})
        set(EIGEN3_INCLUDE_DIR $<BUILD_INTERFACE:${SCOF_EXTERNAL}/eigen> PARENT_SCOPE)
    endif()
endfunction()

## Eigen
function(scof_download_eigen)
    scof_download_project(eigen
            GIT_REPOSITORY      https://gitlab.com/libeigen/eigen.git
            GIT_TAG             master
            )
    add_library (Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories (Eigen3::Eigen SYSTEM INTERFACE $<BUILD_INTERFACE:${SCOF_EXTERNAL}/eigen>)
    #    target_include_directories (Eigen3::Eigen SYSTEM INTERFACE "${ALGOHEX_EXTERNAL}/eigen")
    target_compile_definitions(Eigen3::Eigen INTERFACE -DEIGEN_HAS_STD_RESULT_OF=0)
endfunction()


## GMM
function(scof_download_gmm)
    scof_download_project(gmm
            URL           http://download-mirror.savannah.gnu.org/releases/getfem/stable/gmm-5.4.tar.gz
            URL_HASH SHA256=7163d5080efbe6893d1950e4b331cd3e9160bb3dcf583d206920fba6af7b1e56
            )
    if (${SCOF_DOWNLOAD_MISSINGS_DEPS})
        set(GMM_INCLUDE_DIR $<BUILD_INTERFACE:${SCOF_EXTERNAL}/gmm/include> PARENT_SCOPE)
    endif()
endfunction()

