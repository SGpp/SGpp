# Copyright (C) 2009 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp

# author Dirk Pflueger (Dirk.Pflueger@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), David Pfander (David.Pfander@ipvs.uni-stuttgart.de)

import os, sys, subprocess
import distutils.sysconfig
import glob
import SCons
import fnmatch
import os
import ConfigureExtend

from SConsHelper import *

# Check for versions of Scons and Python
EnsurePythonVersion(2, 7)

# languageWrapperList = ['SG_PYTHON', 'SG_JAVA']
allLanguageWrapperList = ['SG_PYTHON']

# find all modules
allModules = getModules('src/sgpp')

allModuleNames = []
for name in allModules:
    allModuleNames.append('SG_' + name.upper())

vars = Variables("custom.py")

# define the flags 
vars.Add('CPPFLAGS', 'Set additional Flags, they are compiler-depended (multiple flags combined with comma, e.g. -lpython,-lm)', '', converter=multiParamConverter)
vars.Add('LINKFLAGS', 'Set additional Linker-flags, they are linker-depended (multiple flags combined with comma, e.g. -lpython,-lm)', '', converter=multiParamConverter)
# define the target
vars.Add('MARCH', 'Sets the architecture if compiling with gcc, this is a pass-through option: just specify the gcc options!', None)
vars.Add('TARGETCPU', "Sets the processor you are compiling for. 'default' means using gcc with standard configuration. Also available are: 'ICC', here Intel Compiler in version 11 or higher must be used", 'default')
vars.Add('OPT', "Sets optimization on and off", False)
# for compiling on LRZ without errors: omit unit tests
vars.Add(BoolVariable('NO_UNIT_TESTS', 'Omit UnitTests if set to True', False))
# for compiling different modules
for moduleName in allModuleNames:
    vars.Add(BoolVariable(moduleName, 'Build  Module: ' + moduleName, False))
vars.Add(BoolVariable('SG_ALL', 'Build all modules', False))
vars.Add(BoolVariable('SG_PYTHON', 'Build Python Support', False))
vars.Add(BoolVariable('SG_JAVA', 'Build Java Support', False))
vars.Add('OUTPUT_PATH', 'Path where built libraries are installed. Needs a trailing slash!', '')
vars.Add(BoolVariable('VERBOSE', 'Set output verbosity', False))
vars.Add('CMD_LOGFILE', 'Specifies a file to capture the build log', 'build.log')

# initialize environment
env = Environment(variables=vars, ENV=os.environ)

# Help Text
Help("""---------------------------------------------------------------------

Type: 'scons [parameters]' to build the libraries

There are compiler optimizations for different platforms which can be
specified via parameters.

Parameters can be set either by setting the corresponding environment
variables, or directly via the commandline, e.g.,
> scons VERBOSE=True
to enable verbose compilation.


Specifying the target, the following options are available:
    - default: using the gcc toolchain with OpenMP 2
    - ICC: using the ICC toolchain with OpenMP 3 with standard x86_64 options

For LRZ, please execute:
module load python
module load gcc/4.5

FOR LRZ and when using intel compiler, execute:
export LIBPATH=$LD_LIBRARY_PATH

---------------------------------------------------------------------

Parameters are:
""" + 
vars.GenerateHelpText(env))

# adds trailing slashes were required and if not present
env['OUTPUT_PATH'] = os.path.join(env['OUTPUT_PATH'], '')

BUILD_DIR = Dir(env['OUTPUT_PATH'] + 'lib/sgpp')
Export('BUILD_DIR')
PYSGPP_BUILD_DIR = Dir(env['OUTPUT_PATH'] + 'lib/pysgpp')
Export('PYSGPP_BUILD_DIR')
JAVASGPP_BUILD_DIR = Dir(env['OUTPUT_PATH'] + 'lib/javasgpp')
Export('JAVASGPP_BUILD_DIR')
TEST_DIR = Dir(env['OUTPUT_PATH'] + 'tests')
Export('TEST_DIR')

# find out which modules to compile and the folders they reside in
modulesToBuild, modulesToBuildFolder = setupSelectedModules(env, allModules, allModuleNames, allLanguageWrapperList)
setupCompilerDefines(env, modulesToBuild)

# no checks if clean:
if not env.GetOption('clean'):
    print "Checking programs and libraries: "
  
    config = env.Configure(custom_tests={ 'CheckExec' : ConfigureExtend.CheckExec,
                                            'CheckJNI' : ConfigureExtend.CheckJNI,
                                            'CheckFlag' : ConfigureExtend.CheckFlag })
  
    # check scons
    EnsureSConsVersion(2, 0)
    print "Using SCons", SCons.__version__
      
    # check C++11 support    
    if not config.CheckFlag("-std=c++11"):
        sys.stderr.write("Error: compiler doesn't seem to support the C++11 standard. Abort!\n")
        Exit(0)
    config.env.AppendUnique(CPPFLAGS = "-std=c++11")
    
    # check avx support    
    if not config.CheckFlag("-mavx"):
        sys.stderr.write("Error: compiler doesn't seem to support AVX. Abort!\n")
        Exit(0)
    config.env.AppendUnique(CPPFLAGS = "-mavx")

    # check whether swig installed
    if not config.CheckExec('doxygen'):
        sys.stderr.write("Warning: doxygen cannot be found.\n  You will not be able to generate the documentation.\n  Check PATH environment variable!\n")
  
    # check whether dot installed
    if not config.CheckExec('dot'):
        sys.stderr.write("Warning: dot (Graphviz) cannot be found.\n  The documentation might lack diagrams.\n  Check PATH environment variable!\n")
    

  

    if env["SG_PYTHON"]:
        # check whether swig installed
        if not config.CheckExec('swig'):
            sys.stderr.write("Error: swig cannot be found, but required for SG_PYTHON. Check PATH environment variable!\n")
            Exit(0)
        config.env.AppendUnique(CPPPATH=[distutils.sysconfig.get_python_inc()])
        print "pythonpath: ", distutils.sysconfig.get_python_inc()
        if not config.CheckCXXHeader('Python.h'):
            sys.stderr.write("Error: Python.h not found, but required for SG_PYTHON. Check path to Python include files: "
                         + distutils.sysconfig.get_python_inc() + "\n")
            sys.stderr.write("Hint: You might have to install package python-dev\n")
            Exit(0)
           
        if not config.CheckCXXHeader('pyconfig.h'):
            sys.stderr.write("Error: pyconfig.h not found, but required for SG_PYTHON. Check path to Python include files: "
                         + distutils.sysconfig.get_python_inc() + "\n")
            sys.stderr.write("Hint: You might have to install package python-dev\n")
            Exit(0) 
            
        import numpy
        numpy_path = os.path.join(os.path.split(numpy.__file__)[0], "core", "include")
        config.env.AppendUnique(CPPPATH = [numpy_path])
        if not config.CheckCXXHeader(['Python.h', 'pyconfig.h', 'numpy/arrayobject.h']):
            print config.env['CPPPATH']
            sys.stderr.write('Error: Cannot find NumPy header files in: "' + str(numpy_path) + '", required for SG_PYTHON\n')
            Exit(0)
    else:
        print 'Warning: Python extension ("SG_PYTHON") not enabled, skipping python unit tests'
  
    if env['SG_JAVA']:
        # check for $JAVA_HOME; prepend to search path
        if os.environ.get('JAVA_HOME'):
            config.env.PrependENVPath('PATH', os.path.join(os.environ.get('JAVA_HOME'), 'bin'))
        # check whether javac installed
        if not config.CheckExec('javac'):
            sys.stderr.write("Error: javac cannot be found, but required by SG_JAVA. Check PATH environment variable!\n")
            Exit(0)
        # check whether javac installed
        if not config.CheckExec('java'):
            sys.stderr.write("Warning: java cannot be found, but required by SG_JAVA. Check PATH environment variable!\n")
            Exit(0)
      
        # check for JNI headers
        if os.environ.get('JNI_CPPINCLUDE'):
            config.env.AppendUnique(CPPPATH=[os.environ.get('JNI_CPPINCLUDE')])
        if not config.CheckCXXHeader('jni.h'):
            # not found; try to find
            if not config.CheckJNI():
                sys.stderr.write("Error: jni.h not found.\n"
                                 + "Please set JAVA_HOME environment variable "
                                 + "with $JAVA_HOME/bin/javac, $JAVA_HOME/include/jni.h\n"
                                 + "or directly $JNI_CPPINCLUDE with $JNI_CPPINCLUDE/jni.h\n")
                Exit(0)
    else:
        print "Info: Compiling without java support"
        
    # now set up all further environment settings that should never fail    
    # compiler setup should be always after checking headers and flags, as they can make the checks invalid
    # e.g. by setting "-Werror"
    
    # TODO check
    if env['OPT'] == True:
       env.Append(CPPFLAGS=['-O3'])
    else:    
       env.Append(CPPFLAGS=['-g', '-O0'])
    
    if env['TARGETCPU'] == 'default':
        print "Using default gcc"
    
        allWarnings = "-Wall -pedantic -pedantic-errors -Wextra \
            -Wcast-align -Wcast-qual -Wconversion -Wdisabled-optimization -Wformat=2 \
            -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport  -Winit-self  \
            -Winvalid-pch -Wlong-long -Wmissing-field-initializers -Wmissing-format-attribute -Wmissing-include-dirs \
             -Wpacked   -Wpointer-arith -Wredundant-decls -Wstack-protector \
            -Wstrict-aliasing=2 -Wswitch-default -Wswitch-enum -Wunreachable-code -Wunused -Wunused-parameter \
            -Wvariadic-macros -Wwrite-strings".split(" ")
        # rather uninteresting:  -Wpadded -Wshadow -Wfloat-equal -Waggregate-return -Wimplicit -Wmissing-noreturn -Weffc++
        # cannot really use: 
        
        # -Wno-long-long as swig uses long long
        # -fno-strict-aliasing: http://www.swig.org/Doc1.3/Java.html or http://www.swig.org/Release/CHANGES, 03/02/2006
        #    "If you are going to use optimisations turned on with gcc > 4.0 (for example -O2), 
        #     ensure you also compile with -fno-strict-aliasing"
        env.Append(CPPFLAGS=allWarnings + [
                             # '-Wall', '-Wextra', 
                             '-std=c++11',  # '-Wno-long-long', '-Wno-deprecated', 
                             # '-Werror',
                             '-Wno-unused-parameter',
                             # '-Wconversion',
                             '-fno-strict-aliasing',
                             '-funroll-loops', '-mfpmath=sse', '-msse3',
                             '-DDEFAULT_RES_THRESHOLD=-1.0', '-DTASKS_PARALLEL_UPDOWN=4'])
        env.Append(CPPFLAGS=['-fopenmp'])
        env.Append(LINKFLAGS=['-fopenmp'])
        
        if env.has_key('MARCH'):
            env.Append(CPPFLAGS=('-march=' + env['MARCH']))
    
    elif env['TARGETCPU'] == 'ICC':
        print "Using icc"
        env.Append(CPPFLAGS=['-Wall', '-ansi', '-Werror', '-Wno-deprecated', '-wd1125',
                               '-fno-strict-aliasing', '-O3',
                               '-ip', '-ipo', '-funroll-loops', '-msse3',
                               '-ansi-alias', '-fp-speculation=safe',
                               '-DDEFAULT_RES_THRESHOLD=-1.0', '-DTASKS_PARALLEL_UPDOWN=4', '-no-offload'])
    
        env['CC'] = ('icc')
        env['LINK'] = ('icpc')
        env['CXX'] = ('icpc')        
        env.Append(CPPFLAGS=['-openmp'])
        env.Append(LINKFLAGS=['-openmp']) 
    
    else:
        print "You must specify a valid value for TARGETCPU."
        print "Available configurations are: ICC"
        Exit(1)
    
    # special treatment for different platforms
    if env['PLATFORM'] == 'darwin':
        # the "-undefined dynamic_lookup"-switch is required to actually build a shared library 
        # in OSX. "-dynamiclib" alone results in static linking of all further dependent shared libraries
        # beware: if symbols are missing that are actually required (because the symbols don't reside in a shared library), there will be no error during compilation
        # the python binding (pysgpp) requires lpython and a flat namespace
        # also for the python binding, the library must be suffixed with '*.so' even though it is a dynamiclib and not a bundle (see SConscript in src/pysgpp)
        env.Append(LINKFLAGS=['-flat_namespace', '-undefined', 'dynamic_lookup', '-lpython'])
        env['SHLIBSUFFIX'] = '.dylib'
    elif env['PLATFORM'] == 'cygwin':
        # required to find the static libraries compiled before the shared libraries
        # the static libraries are required as the linker on windows cannot ignore undefined symbols
        # (as is done on linux automatically and is done on OSX with the settings above)
        env.Append(LIBPATH=[BUILD_DIR])
    
    # will lead to a warning on cygwin (and we have -Werror enabled)
    # is enabled by default on cygwin
    if env['PLATFORM'] != 'cygwin':
        env.Append(CPPFLAGS=['-fPIC'])
        
    # setup the include base folder
    env.Append(CPPPATH=['#/src/sgpp'])
  
    env = config.Finish()

# clear build_log file
logfile = open(env['CMD_LOGFILE'], 'a')
logfile.seek(0)
logfile.truncate()

# detour compiler output
env['PRINT_CMD_LINE_FUNC'] = print_cmd_line

# environement setup finished, export environment
Export('env')
  
# compile selected modules
moduleDependencies = {}
all_objs = []
all_srcs = []
for name in modulesToBuildFolder:
    if name == "test":
        print "Preparing to build module: ", "test"
        env.SConscript('src/sgpp/' + "test" + '/SConscript', {'env': env, 'moduleName': name})
        continue
    else:
        print "Preparing to build module: ", name
        # SConscript('src/sgpp/SConscript' + name, variant_dir='#/tmp/build/', duplicate=0)
        env.SConscript('src/sgpp/' + name + '/SConscript', {'env': env, 'moduleName': name})
    
        Import('dependencies')
        moduleDependencies['SG_' + name.upper()] = dependencies
    
        print 'Module SG_' + name.upper() + ' depends on:'
        for dep in dependencies:
            print '\t' + dep
            if not dep in modulesToBuildFolder:
                print "Error!"
                print name + " depends on non-existent module " + dep
                Exit(1)
        Import('srcs')
        Import('objs')
        all_objs += objs
        all_srcs += srcs
        # src_objs[name] = objs
        # src_files[name] = srcs
     
Export('moduleDependencies')
Export('all_objs')
Export('all_srcs')

#   
# build python lib
if env['SG_PYTHON']:
    libpysgpp = SConscript('src/pysgpp/SConscript')

# build java lib
if env['SG_JAVA']:
    libjsgpp = env.SConscript('src/jsgpp/SConscript',
                              variant_dir='tmp/build_jsgpp', duplicate=0)
    # install
    jinst = env.Install(env['OUTPUT_PATH'] + 'lib/jsgpp', [libjsgpp])
    
  
# Unit tests
#########################################################################
  
if not env['NO_UNIT_TESTS'] and env['SG_PYTHON']:
    testdep = env.SConscript('tests/SConscript')
    # execute after all installations (even where not necessary)
    if env['SG_JAVA']:
        Depends(testdep, [jinst, pyinst])
    else:
        Depends(testdep, [pyinst])
else:
    print "Warning: Skipping python unit tests"

