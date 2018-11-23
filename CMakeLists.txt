#
# Copyright (C) 2015   Michael Lahnert
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

project (NUMSIM)

cmake_minimum_required (VERSION 2.8)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

#set (CMAKE_CC_COMPILER "/usr/bin/gcc")
#set (CMAKE_CXX_COMPILER "/usr/bin/g++")

set (GOOD_FLAGS "-Wall -Wextra -pedantic -std=c++11")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GOOD_FLAGS}")

# debug visu?
option (DEBUG_VISU
  "Enable debugging visualization" OFF)

# find sdl 2
INCLUDE(FindPkgConfig)
PKG_SEARCH_MODULE(SDL2 REQUIRED sdl2)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${SDL2_INCLUDE_DIRS})

# find mpi
find_package (MPI REQUIRED)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${MPI_INCLUDE_PATH})

# Using glob is discouraged in cmake.
set (HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/communicator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/compute.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/geometry.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/grid.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/iterator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/parameter.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/solver.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/typedef.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/vtk.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/zeitgeist.hpp
  )

set (SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/communicator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/compute.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/geometry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/grid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/iterator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/parameter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/solver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/vtk.cpp
  )

if (DEBUG_VISU)
  set (HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/visu.hpp
    )

  set (SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolver/src/visu.cpp
    )

  add_definitions(-DUSE_DEBUG_VISU)
endif (DEBUG_VISU)

add_executable (NumSim ${SOURCES})
target_link_libraries (NumSim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})

set (HEAD
  ${CMAKE_CURRENT_SOURCE_DIR}/Magrathea/src/tga.hpp
  )

set (SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/Magrathea/src/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Magrathea/src/tga.cpp
  )

add_executable (magrathea ${SRCS})

if (MPI_COMPILE_FLAGS)
  set_target_properties(NumSim PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if (MPI_LINK_FLAGS)
  set_target_properties (NumSim PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
