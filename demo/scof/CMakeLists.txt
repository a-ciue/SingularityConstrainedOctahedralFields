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

cmake_minimum_required(VERSION 3.4)

find_package (OpenVolumeMesh REQUIRED)
find_package (OpenMesh REQUIRED)
find_package (GMM REQUIRED)
find_package (EIGEN3 REQUIRED)


# disable OpenVolumeMesh debug output
add_definitions(-DNDEBUG)

set (CMDLINE_TOOL_SOURCES "")
set (CMDLINE_TOOL_INCLUDE_DIRS "")
set (CMDLINE_TOOL_LIBRARIES "")
set (CMDLINE_TOOL_LIBRARY_DIRS "")

list (APPEND CMDLINE_TOOL_SOURCES main.cc)
list (APPEND CMDLINE_TOOL_INCLUDE_DIRS
    ${OPENVOLUMEMESH_INCLUDE_DIR}
    ${OPENMESH_INCLUDE_DIRS} 
    ${GMM_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src
    ${SCOF_INCLUDE_DIRS}
)
list (APPEND CMDLINE_TOOL_LIBRARIES
    ${OPENVOLUMEMESH_LIBRARY}
    ${OPENMESH_LIBRARIES}
    SCOFStatic
)

list (APPEND CMDLINE_TOOL_LIBRARY_DIRS ${OPENVOLUMEMESH_LIBRARY_DIR})


add_executable (scof ${CMDLINE_TOOL_SOURCES})

LINK_DIRECTORIES(scof ${CMDLINE_TOOL_LIBRARY_DIRS} )
target_link_libraries(scof ${CMDLINE_TOOL_LIBRARIES})

if (NOT WIN32 )
    set(SCOF_COMPILE_FLAGS "-Wall -Wextra -pedantic -std=c++11 -O3 -Wno-long-long")
else()
    set(SCOF_COMPILE_FLAGS "")
endif()

set_target_properties(scof
    PROPERTIES
    INCLUDE_DIRECTORIES "${CMDLINE_TOOL_INCLUDE_DIRS}"
    COMPILE_FLAGS "${SCOF_COMPILE_FLAGS}"
)

