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

Import('*')

# get sources.
srcs = ['src/main.cpp',
        'src/communicator.cpp',
        'src/compute.cpp',
        'src/geometry.cpp',
        'src/grid.cpp',
        'src/iterator.cpp',
        'src/parameter.cpp',
        'src/solver.cpp',
        'src/vtk.cpp',
       ]

# check if debug-visualization should be build.
# if so, append its source file to sources and set preproc. define
if env['visu'] == 1:
    srcs.append('src/visu.cpp')
    env.Append(CPPDEFINES = ['USE_DEBUG_VISU'])

# give the program a name
name = 'NumSim'

# build it
env.Program(name, srcs)
