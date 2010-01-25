#############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2008-2009 Dirk Pflueger (dirk.pflueger@in.tum.de)           #
# Copyright (C) 2008 Joerg Blank (blankj@in.tum.de)                         #
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
import distutils.sysconfig

vars = Variables("custom.py")

# define the flags 
vars.Add('CPPFLAGS','Set additional Flags, they are compiler-depended','')
vars.Add('LINKFLAGS','Set additional Linker-flags, they are linker-depended','')

# define the target
vars.Add('MARCH','Sets the architecture if compiling with gcc, this is a pass-through option: just specify the gcc options!', None)
vars.Add('TARGETCPU',"Sets the processor you are compiling for. 'default' means using gcc with standard configuration. Also available are: 'opteronICC', 'core2ICC', 'ia64ICC'; here Intel Compiler in version 11 must be used", 'default')
vars.Add('OMP', "Sets if OpenMP should be used; with gcc OpenMP 2 is used, with all icc configurations OpenMP 3 is used!", False)
vars.Add('TRONE', "Sets if the tr1/unordered_map should be uesed", False)

# for building the the jsgpp lib
vars.Add('JSGPP', 'Build jsgpp if set to True', False)
vars.Add('JNI_CPPPATH', 'Path to JNI includes', None)
vars.Add('JNI_OS', 'JNI os path', None)

env = Environment(variables = vars, ENV = os.environ)

# Specifying the target
# there are several target avialable:
# 	- default: using the gcc toolchain with OpenMP 2
#	- opteronICC: using the ICC 11.x toolchain with OpenMP 3 with standard x86_64 options
#	- core2ICC: using the ICC 11.x toolchain with OpenMP 3 with Intel x86_64 options
#	- ia64ICC: using the ICC 11.x toolchain with OpenMP 3 with Itanium options
#
# Take care that you have defined following env. variables for loading the 
# shared libraries: LD_LIBRARY_PATH and LIBPATH
# both must contain the path to the intel shared libs
# for instance:
# LD_LIBRARY_PATH = /opt/intel/cce/default/lib:LD_LIBRARY_PATH
# LIBPATH = /opt/intel/cce/default/lib:LIBPATH
#
# FOR LRZ:
# lib: /lrz/sys/intel/icc_110_074/lib/ia64/
# bin: /lrz/sys/intel/icc_110_074/bin/ia64/

if env['TRONE']:
    env.Append(CPPDEFINES=['USETRONE'])
    env.Append(CPPFLAGS=['-std=c++0x'])

if env['TARGETCPU'] == 'default':
    print "Using default gcc"
    # -Wno-long-long as swig uses long long
    # -fno-strict-aliasing: http://www.swig.org/Doc1.3/Java.html or http://www.swig.org/Release/CHANGES, 03/02/2006
    #    "If you are going to use optimisations turned on with gcc > 4.0 (for example -O2), 
    #     ensure you also compile with -fno-strict-aliasing"
    env.Append(CPPFLAGS=['-Wall', '-ansi', '-pedantic', '-Wno-long-long', 
                         '-fno-strict-aliasing', '-fopenmp', '-O3', '-g',
                         '-funroll-loops', '-pthread', '-ffloat-store'])
    if env['OMP']:
    	env.Append(CPPDEFINES=['USEOMP'])
    	env.Append(LINKFLAGS=['-fopenmp'])
    	
elif env['TARGETCPU'] == 'ia64ICC':
    print "Using icc 11.0 for Itanium systems"
    # ICC doesn't know '-pedantic'
    # ICC has different options on ia64
    env.Append(CPPFLAGS = ['-O3', '-fno-fnalias', '-funroll-loops', '-no-alias-const', 
                           '-no-ansi-alias', '-i-static', '-gcc-version=400', 
                           '-unroll-aggressive', '-opt-jump-tables=large', '-Wall', 
                           '-ansi', '-wd981', '-fno-strict-aliasing', '-openmp', '-pthread']) 
elif env['TARGETCPU'] == 'opteronICC':
    print "Using icc 11.x for Opteron systems"
    env.Append(CPPFLAGS = ['-axSSE3', '-O3', '-funroll-loops', '-ipo', '-ip', '-fno-fnalias', 
                           '-no-alias-const', '-no-ansi-alias', '-Wall', '-ansi', '-wd981', 
                           '-fno-strict-aliasing', '-openmp', '-pthread'])
elif env['TARGETCPU'] == 'core2ICC':
    print "Using icc 11.x for Core2 systems"
    env.Append(CPPFLAGS = ['-axSSE3', '-O3', '-funroll-loops', '-ipo', '-ip', '-fno-fnalias', 
                           '-no-alias-const', '-no-ansi-alias', '-Wall', '-ansi', '-wd981', 
                           '-fno-strict-aliasing', '-openmp', '-pthread'])
else:
    print "You must specify a valid value for TARGETCPU."
    print "Available configurations are: default, core2ICC, opteronICC, ia64ICC"
    Exit(1)
    
# sets ICC-wide commen options and the tool chain   
if env['TARGETCPU'] in ['ia64ICC', 'opteronICC', 'core2ICC']:
    env['CC'] = ('icc')
    env['LINK'] = ('icpc')
    env['CXX'] = ('icpc')	    
    if env['OMP']:
        env.Append(LINKFLAGS=['-openmp']) 
        env.Append(CPPDEFINES=['USEOMP', 'USEOMPTHREE'])
    
# sets the architecture option for gcc
if env.has_key('MARCH'):
    if env['TARGETCPU'] == 'default':
        env.Append(CPPFLAGS=('-march=' + env['MARCH']))
    else:
        print "Warning: Ignoring option MARCH"
         
# add path of python includes
env.Append(CPPPATH=[distutils.sysconfig.get_python_inc()])

if not env.GetOption('clean'):	
    config = env.Configure()
	
    # check if the intel omp lib is available
    if env['TARGETCPU'] in ['ia64ICC', 'opteronICC', 'core2ICC'] and env['OMP']:
        if not config.CheckLib('iomp5'):
            print "Error: Intel omp library iomp5 is missing."
            Exit(1)
               
    # check if the the intel vector lib is available
    if env['TARGETCPU'] in ['ia64ICC', 'opteronICC', 'core2ICC']:
        if not config.CheckLib('svml'):
            print "SVML should be available when using intelc. Consider runnning scons --config=force!"

    # check if the math header is available
    if not config.CheckLibWithHeader('m', 'math.h', 'c++'):
        print "Error: Math headers are missing."
        Exit(1)
        
    # check if the Python headers are available
    # @todo (heinecke) some old things that not work, should be fixed
#    if not config.CheckCHeader('Python.h'):
#        print "Error: Python.h not found. Check path to Python include files."
#        print distutils.sysconfig.get_python_inc()
#        Exit(1)

    env = config.Finish()

Export('env')

#start build of pysgpp and jsgpp
SConscript('src/sgpp/SConscript', build_dir='tmp/build_sg', duplicate=0)
SConscript('src/pysgpp/SConscript', build_dir='tmp/build_pysgpp', duplicate=0)
if env['JSGPP']:
    SConscript('src/jsgpp/SConscript', build_dir='tmp/build_jsgpp', duplicate=0)


# Copy required files
cpy = []
cpy += Command("#lib/pysgpp/_pysgpp.dylib", "#/tmp/build_pysgpp/_pysgpp.dylib", Copy("$TARGET", "$SOURCE"))
cpy += Command("#lib/pysgpp/pysgpp.py", "#/tmp/build_pysgpp/pysgpp.py", Copy("$TARGET", "$SOURCE"))
cpy += Command("#bin/_pysgpp.dylib", "#/tmp/build_pysgpp/_pysgpp.dylib", Copy("$TARGET", "$SOURCE"))
cpy += Command("#bin/pysgpp.py", "#/tmp/build_pysgpp/pysgpp.py", Copy("$TARGET", "$SOURCE"))
cpy += Command("#lib/sgpp/libsgpp.a", "#/tmp/build_sg/libsgpp.a", Copy("$TARGET", "$SOURCE"))
cpy += Command("#bin/sgpp.a", "#/tmp/build_sg/libsgpp.a", Copy("$TARGET", "$SOURCE"))

# Execute Unit Tests
SConscript('tests/SConscript')


Help(vars.GenerateHelpText(env))
