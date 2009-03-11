#############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2007 Joerg Blank (blankj@in.tum.de)                         #
# Copyright (C) 2007 Richard Roettger (roettger@in.tum.de)                  #
# Copyright (C) 2008 Dirk Plueger (pflueged@in.tum.de)                      #
# Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################

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

####### OMP #######
# OpenMP parallelisation directives in version 3
#OMP = 1

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


