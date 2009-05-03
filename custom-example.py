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
# it under the terms of the GNU Lesser General Public License as published  #
# by the Free Software Foundation; either version 3 of the License, or      #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU Lesser General Public License  #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################

import os

###############################
# Copy this file to custom.py #
###############################

####### GCC Settings #######
# additional compiler flags
CPPFLAGS = ['-O3','-g','-funroll-loops', '-Wall', '-ansi', '-Wno-long-long', '-pedantic']

# MARCH variable for compiler optimization. See gcc man page
#MARCH = 'nocona'
#MARCH = 'opteron'
#MARCH = 'pentium4'

####### JSGPP #######
# Uncomment to compile jsgpp, a version of sgpp that works with java
# you have to define serveral paths
# JSGPP = 1
# JNI_CPPPATH = '/usr/java/jdk1.6.0_11/include'
# JNI_OS = 'linux'

####### OMP #######
# OpenMP parallelisation directives in version 2 or 3
#OMPTWO = 1        #tasks are not supported, badder performance of Laplace
#OMPTHREE = 1      #good parallel performance of Laplace due to tasks, is ignored if set with g++

####### ICC #######
# Uncomment to use Intels optimizing Compiler, please use version 11
# Take care that you have defined following env. variables for loading the shared libraries: LD_LIBRARY_PATH and LIBPATH
# both must contain the path to the intel shared libs
# for instance:
# LD_LIBRARY_PATH = /opt/intel/cce/default/lib/intel64:LD_LIBRARY_PATH
# LIBPATH = /opt/intel/cce/default/lib/intel64:LIBPATH
#
# FOR LRZ:
# lib: /lrz/sys/intel/icc_110_074/lib/ia64/
# bin: /lrz/sys/intel/icc_110_074/bin/ia64/ 

#ICC = 1

# FOR YOUR STANDARD X86 X86_64 System USE:
#INTELHOME = '/opt/intel/cce/default/bin/'
#CPPFLAGS = ['-axSSE3', '-O3', '-funroll-loops', '-ipo', '-intel-static', '-ip', '-fno-fnalias', '-no-alias-const', '-no-ansi-alias', '-Wall', '-ansi', '-wd981']

# FOR LRZ
#INTELHOME = '/lrz/sys/intel/icc_110_074/bin/ia64/'
#CPPFLAGS = ['-O3', '-alias-args', '-fno_alias', '-fno-fnalias', '-funroll-loops', '-no-alias-const', '-no-ansi-alias', '-i-static', '-gcc-version=400', '-unroll-aggressive', '-opt-jump-tables=never', '-Wall', '-ansi', '-wd981']
