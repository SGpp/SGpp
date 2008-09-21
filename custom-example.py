# pyclass is a set of tools utilizing sparse grids to solve numerical problems
# Copyright (C) 2007  Joerg Blank (blankj@in.tum.de), Richard Roettger (roettger@in.tum.de)
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



# Copy this file to custom.py


# additional compiler flags
CPPFLAGS = ['-O3','-ggdb','-funroll-loops']

# debug flags. You have to disable inline or you'll get nowhere
#CPPFLAGS = ['-O0','-ggdb']

# MARCH variable for compiler optimization. See gcc man page
MARCH = 'nocona'  # Core2
#MARCH = 'opteron'
#MARCH = 'pentium4'


####### ICC #######
# Uncomment to use Intels optimizing Compiler

#ICC = 1

#CPPFLAGS = ['-axN', '-O3', '-xN', '-funroll-loops', '-ipo']

###### JSGPP ######
# Uncomment to enable jsgpp

#JSGPP = True

#JNI_CPPPATH = ''
#JNI_LIBPATH = ''
#JNI_OS = 'linux'



