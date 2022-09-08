# This file is part of SCOF.
#
# SCOF is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# SCOF is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License
# along with SCOF.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required (VERSION 3.9)

# AutoMOC/AutoUIC
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.17")
    cmake_policy(SET CMP0100 NEW)
endif()

list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project (SCOF
        VERSION 0.0.1
        LANGUAGES C CXX
        )

set(CMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY ON)

option(BUILD_SHARED_LIBS "Build shared library (*.dll, *.so, *.dylib) instead of static library (*.a, *.lib)" ON)

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    if (WIN32)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build")
    else ()
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/bin")
    endif ()

    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")

    set(CMAKE_CXX_STANDARD 17)
endif()

set(SCOF_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set (ACG_COMMON_DO_NOT_COPY_POST_BUILD TRUE)
include (SCOFDependencies)

# add library
add_library (SCOF
        src/SCOF/AxisAlignment.hh
        src/SCOF/FrameFieldGeneratorT.hh
        src/SCOF/SingularityConstrainedOctahedralField.cc
        src/SCOF/SingularityGraphT.hh
        src/SCOF/TransitionQuaternionEigen.cc
        src/SCOF/CornerT.hh
        src/SCOF/MatchingsAndAlignmentExtractionT.cc
        src/SCOF/SingularityConstrainedOctahedralField.hh
        src/SCOF/StopWatch.hh
        src/SCOF/TransitionQuaternionEigen.hh
        src/SCOF/FrameFieldGeneratorT.cc
        src/SCOF/MatchingsAndAlignmentExtractionT.hh
        src/SCOF/SingularityGraphT.cc
        src/SCOF/TangentContinuityT.hh
        src/SCOF/Typedefs.hh
        )

add_library (SCOF::SCOF ALIAS SCOF)

target_compile_options(SCOF PUBLIC
        "$<$<CXX_COMPILER_ID:MSVC>:/bigobj>"
        # cmake by default adds /W3, that's too much for us:
        "$<$<CXX_COMPILER_ID:MSVC>:/W1>"
        # disable MSVC inline warnings that are enabled with /Ob2, which is in the default cmake release cxxflags:
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4710>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4711>"

        "$<$<CXX_COMPILER_ID:MSVC>:/wd4514>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4571>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd5045>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd5024>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd5025>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd5026>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd5027>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4625>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4626>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4820>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4464>"
        )
target_compile_definitions(SCOF
        PUBLIC
        "INCLUDE_TEMPLATES"
        "$<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_DEPRECATE>" # for GMM
        "$<$<CXX_COMPILER_ID:MSVC>:_USE_MATH_DEFINES>"
        )

set_target_properties(SCOF PROPERTIES
        VERSION ${SCOF_VERSION_MAJOR}.${SCOF_VERSION_MINOR}
        SOVERSION ${SCOF_VERSION_MAJOR}.${SCOF_VERSION_MINOR}
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        CXX_VISIBILITY_PRESET hidden
        )

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set_target_properties(SCOF PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
            ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
            )
endif()


configure_file (
        "${CMAKE_CURRENT_SOURCE_DIR}/src/SCOF/Config/Version.hh.in"
        "${CMAKE_CURRENT_BINARY_DIR}/src/SCOF/Config/Version.hh"
)

include(GenerateExportHeader)
generate_export_header(SCOF
        BASE_NAME SCOF
        EXPORT_FILE_NAME "src/Config/Export.hh"
        )

if(CHOLMOD_INCLUDES)
    add_definitions(-DENABLE_SUITESPARSE)
    target_link_libraries (SCOF PUBLIC
            OpenMeshCore
            OpenMeshTools
            OpenVolumeMesh::OpenVolumeMesh
            Eigen3::Eigen
            ${CHOLMOD_LIBRARIES}
            )
    #Set target properties
    target_include_directories(SCOF
            PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
            #    PRIVATE
            #"${CMAKE_CURRENT_BINARY_DIR}/src"
#            ${EIGEN3_INCLUDE_DIR}
            ${GMM_INCLUDE_DIR}
            ${CHOLMOD_INCLUDES}
            )
else(CHOLMOD_INCLUDES)
    target_link_libraries (SCOF PUBLIC
            OpenMeshCore
            OpenMeshTools
            OpenVolumeMesh::OpenVolumeMesh
            Eigen3::Eigen
            )
    #Set target properties
    target_include_directories(SCOF
            PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
            #    PRIVATE
            #"${CMAKE_CURRENT_BINARY_DIR}/src"
            ${GMM_INCLUDE_DIR}
#            ${EIGEN3_INCLUDE_DIR}
            )
endif()

#If set found by upper level, eigen needs it to compile
if(OpenMP_CXX_FOUND)
    target_link_libraries(SCOF PUBLIC OpenMP::OpenMP_CXX OpenMP::OpenMP_C)
endif()

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_subdirectory(demo/scof)
endif()


include(CMakePackageConfigHelpers)
include(GNUInstallDirs)


set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/SCOF)

configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/CMAKE/SCOFConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/SCOFConfig.cmake"
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

write_basic_package_version_file(
        SCOFConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/SCOFConfigVersion.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/SCOFConfig.cmake"
        DESTINATION ${INSTALL_CONFIGDIR})


install(DIRECTORY
        "${CMAKE_CURRENT_BINARY_DIR}/src/SCOF/Config"
        DESTINATION include/SCOF)

# Install Header Files
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/SCOF
        DESTINATION include
        FILES_MATCHING PATTERN "*.hh"
        )

install(TARGETS SCOF
        EXPORT SCOFTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )

install(EXPORT SCOFTargets
        FILE SCOFTargets.cmake
        NAMESPACE SCOF::
        DESTINATION ${INSTALL_CONFIGDIR}
        )

export(EXPORT SCOFTargets
        NAMESPACE SCOF::)
export(PACKAGE SCOF)

list(REMOVE_AT CMAKE_MODULE_PATH -1)
