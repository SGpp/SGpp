# pyclass is a set of tools utilizing sparse grids to solve numerical problems
# Copyright (C) 2007  Joerg Blank (blankj@in.tum.de), Richard Roettger (roettger@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)
#               2008  Dirk Pflueger (pflueged@in.tum.de)
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

import os

###############################
# Copy this file to custom.py #
###############################


####### STLPort #######
# path to STLport
# global installation e.g.
STLPORT = '/usr/include/stlport'
# local installation (see documentation)
#STLPORT = os.path.expanduser('~/include/stlport')



####### GCC Settings #######
# additional compiler flags
CPPFLAGS = ['-O3','-g','-funroll-loops']

# MARCH variable for compiler optimization. See gcc man page
#MARCH = 'nocona'  # Core2
#MARCH = 'opteron'
#MARCH = 'pentium4'

# Set if Scons-Version = v0.96.1.D001
#OLDSCONS = 1

# Set to True to vectorize certain files
#VECTORIZE = True



####### ICC #######
# Uncomment to use Intels optimizing Compiler, please use version 11
# Take care that you have defined following env. variables for loading the shared libraries: LD_LIBRARY_PATH and LIBPATH
# both must contain the path to the intel shared libs
# for instance:
# LD_LIBRARY_PATH = /opt/intel/cce/default/lib/intel64:LD_LIBRARY_PATH
# LIBPATH = /opt/intel/cce/default/lib/intel64:LIBPATH

#ICC = 1
#INTELHOME = '/opt/intel/cce/default/bin/'
#CPPFLAGS = ['-axSSE3', '-O3', '-funroll-loops', '-ipo', '-intel-static', '-ip', '-fno-fnalias', '-no-alias-const', '-no-ansi-alias']


