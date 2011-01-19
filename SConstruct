# Copyright (C) 2009 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp

# author Dirk Pflueger (Dirk.Pflueger@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)


import os
import distutils.sysconfig

vars = Variables("custom.py")

# define the flags 
vars.Add('CPPFLAGS','Set additional Flags, they are compiler-depended','-Wno-deprecated')
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

# for compiling on LRZ without errors: omit unit tests
vars.Add('NO_UNIT_TESTS', 'Omit UnitTests if set to True', False)


env = Environment(variables = vars, ENV = os.environ)

# Specifying the target
# there are several targets avialable:
# 	- default: using the gcc toolchain with OpenMP 2
#	- opteronICC: using the ICC 11.x toolchain with OpenMP 3 with standard x86_64 options
#	- core2ICC: using the ICC 11.x toolchain with OpenMP 3 with Intel x86_64 options (core architecture)
#   - nehalemICC: using the ICC 11.x toolchain with OpenMP 3 with Intel x86_64 options (nehalem architecture)
#	- ia64ICC: using the ICC 11.x toolchain with OpenMP 3 with Itanium options
#
# FOR LRZ please execute:
# module load python
# module load gcc/4.5
#
# FOR LRZ and when using intel compiler:
#
# execute:
# export LIBPATH=$LD_LIBRARY_PATH


# scons usually adds double quotes around the command-line arguments containing 
# white spaces  this whould produce compilation error, therefore replace string 
# with corresponding list of parameters
opt_flags = Split(env['CPPFLAGS'])
env['CPPFLAGS'] = []

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
                         '-fno-strict-aliasing', '-O3',
                         '-funroll-loops', '-ffloat-store'])
    if env['OMP']:
	env.Append(CPPFLAGS=['-fopenmp'])
    	env.Append(CPPDEFINES=['USEOMP'])
    	env.Append(LINKFLAGS=['-fopenmp'])
    	
elif env['TARGETCPU'] == 'ia64ICC':
    print "Using icc 11.0 for Itanium systems"
    # ICC doesn't know '-pedantic'
    # ICC has different options on ia64
    env.Append(CPPFLAGS = ['-O3', '-funroll-loops', 
                           '-no-alias', '-i-static', '-gcc-version=400', 
                           '-unroll-aggressive', '-opt-jump-tables=large', '-Wall', 
                           '-ansi', '-wd981', '-fno-strict-aliasing']) 
elif env['TARGETCPU'] == 'opteronICC':
    print "Using icc 11.x/12.0 for Opteron systems"
    env.Append(CPPFLAGS = ['-axSSE3', '-O3', '-funroll-loops', '-ipo', '-ip', '-ansi-alias', 
                           '-Wall', '-ansi', '-wd981', 
                           '-fno-strict-aliasing'])
elif env['TARGETCPU'] == 'core2ICC':
    print "Using icc 11.x/12.0 for Core2 systems"
    env.Append(CPPFLAGS = ['-axSSE3', '-O3', '-funroll-loops', '-ipo', '-ip', '-ansi-alias', 
                           '-Wall', '-ansi', '-wd981', 
                           '-fno-strict-aliasing'])
elif env['TARGETCPU'] == 'nehalemICC':
    print "Using icc 11.x/12.0 for Nehalem/Westmere systems"
    env.Append(CPPFLAGS = ['-axSSE4.1', '-O3', '-funroll-loops', '-ipo', '-ip', '-ansi-alias', 
                           '-Wall', '-ansi', '-wd981', 
                           '-fno-strict-aliasing'])
else:
    print "You must specify a valid value for TARGETCPU."
    print "Available configurations are: default, core2ICC, opteronICC, ia64ICC"
    Exit(1)
    
# sets ICC-wide commen options and the tool chain   
if env['TARGETCPU'] in ['ia64ICC', 'opteronICC', 'core2ICC', 'nehalemICC']:
    env['CC'] = ('icc')
    env['LINK'] = ('icpc')
    env['CXX'] = ('icpc')	    
    if env['OMP']:
	env.Append(CPPFLAGS=['-openmp'])
        env.Append(LINKFLAGS=['-openmp']) 
        env.Append(CPPDEFINES=['USEOMP', 'USEOMPTHREE', 'USEICCINTRINSICS'])
    
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
    if env['TARGETCPU'] in ['ia64ICC', 'opteronICC', 'core2ICC', 'nehalemICC'] and env['OMP']:
        if not config.CheckLib('iomp5'):
            print "Error: Intel omp library iomp5 is missing."
            Exit(1)
               
    # check if the the intel vector lib is available
    if env['TARGETCPU'] in ['ia64ICC', 'opteronICC', 'core2ICC', 'nehalemICC']:
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


# the optional CPPFLAGS at the end will override the previous flags
env['CPPFLAGS'] = env['CPPFLAGS'] + opt_flags

Export('env')

#start build of pysgpp and jsgpp
SConscript('src/sgpp/SConscript', build_dir='tmp/build_sg', duplicate=0)
SConscript('src/pysgpp/SConscript', build_dir='tmp/build_pysgpp', duplicate=0)
if env['JSGPP']:
    SConscript('src/jsgpp/SConscript', build_dir='tmp/build_jsgpp', duplicate=0)
    SConscript('src/jsgpp_weka/SConscript', build_dir='tmp/build_jsgpp_weka', duplicate=0)

# Copy required files
cpy = []
cpy += Command("#lib/pysgpp/_pysgpp.so", "#/tmp/build_pysgpp/_pysgpp.so", Copy("$TARGET", "$SOURCE"))
cpy += Command("#lib/pysgpp/pysgpp.py", "#/tmp/build_pysgpp/pysgpp.py", Copy("$TARGET", "$SOURCE"))
cpy += Command("#bin/_pysgpp.so", "#/tmp/build_pysgpp/_pysgpp.so", Copy("$TARGET", "$SOURCE"))
cpy += Command("#bin/pysgpp.py", "#/tmp/build_pysgpp/pysgpp.py", Copy("$TARGET", "$SOURCE"))
cpy += Command("#lib/sgpp/libsgpp.a", "#/tmp/build_sg/libsgpp.a", Copy("$TARGET", "$SOURCE"))
cpy += Command("#bin/sgpp.a", "#/tmp/build_sg/libsgpp.a", Copy("$TARGET", "$SOURCE"))

# Execute Unit Tests
if not env['NO_UNIT_TESTS']:
    SConscript('tests/SConscript')


# Help Text
Help("""Type: 'scons [parameters]' to build the libraries

There are compiler optimizations for different platforms which can be
specified via parameters.

Parameters can be set either by setting the corresponding environment
variables, or directly via the commandline, e.g.,
> scons OMP=True
to enable OpenMP support.

---------------------------------------------------------------------

Parameters are:
""" +
vars.GenerateHelpText(env))
