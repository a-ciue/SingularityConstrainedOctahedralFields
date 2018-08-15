 #
# Copyright 2013 Computer Graphics Group, RWTH Aachen University
# Authors:
#
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

cmake_minimum_required (VERSION 2.6)

# Only set project name if not build from within another project.
if("${PROJECT_NAME}" STREQUAL "")
    project (SCOF)
endif()

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

add_definitions(-DINCLUDE_TEMPLATES)

set (SOURCES "")
set (INCLUDE_DIRS "")
set (LIBRARIES "")
set (LIBRARY_DIRS "")


file(GLOB SRCS "${CMAKE_CURRENT_LIST_DIR}/src/*.cc" "${CMAKE_CURRENT_LIST_DIR}/src/*.c")


find_package (OpenVolumeMesh REQUIRED)
find_package (OpenMesh REQUIRED)
find_package (EIGEN3 REQUIRED)
find_package (GMM REQUIRED)
find_package (CHOLMOD)

#enable suitesparse
if(CHOLMOD_INCLUDES)
	add_definitions(-DENABLE_SUITESPARSE)
	list (APPEND EIGEN3_INCLUDE_DIR ${CHOLMOD_INCLUDES})
	list (APPEND LIBRARIES ${CHOLMOD_LIBRARIES})
endif(CHOLMOD_INCLUDES)

list (APPEND SOURCES ${SRCS})
list (APPEND INCLUDE_DIRS ${OPENVOLUMEMESH_INCLUDE_DIR}
	${OPENMESH_INCLUDE_DIRS} 
	${EIGEN3_INCLUDE_DIR} 
	${GMM_INCLUDE_DIR})
list (APPEND LIBRARY_DIRS ${OPENVOLUMEMESH_LIBRARY_DIR})
list (APPEND LIBRARIES ${OPENVOLUMEMESH_LIBRARY}
	${OPENMESH_LIBRARIES})


include_directories (${INCLUDE_DIRS})

link_directories (${CMAKE_BINARY_DIR} ${LIBRARY_DIRS})

add_library (SCOF SHARED ${SOURCES})
add_library (SCOFStatic STATIC ${SOURCES})

target_link_libraries (SCOF ${LIBRARIES})
target_link_libraries (SCOFStatic ${LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -Wall -Wextra")
#
# Fake successful finder run if compiling as a dependent project.
#
if(NOT "${PROJECT_NAME}" MATCHES "SCOF")
    set (SCOF_FOUND true PARENT_SCOPE)
    set (SCOF_LIBRARIES SCOF PARENT_SCOPE)
    set (SCOF_LIBRARY SCOF PARENT_SCOPE)
    set (SCOF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src/" PARENT_SCOPE)
    set (SCOF_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE PATH "The directory where the SCOF libraries can be found.")
endif()

if("${PROJECT_NAME}" MATCHES "SCOF")
   add_subdirectory(demo/scof)
endif()
