# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import distutils.sysconfig
import os
import sys
import commands
import re

import SGppConfigureExtend
import Helper

def doConfigure(env, moduleFolders, languageWrapperFolders):
    print "Checking programs and libraries: "

    config = env.Configure(custom_tests={ 'CheckExec' : SGppConfigureExtend.CheckExec,
                                            'CheckJNI' : SGppConfigureExtend.CheckJNI,
                                            'CheckFlag' : SGppConfigureExtend.CheckFlag })    
    # boost library
    # TODO add check and error
    config.env.AppendUnique(CPPPATH=config.env['BOOST_INCLUDE_PATH'])

    # check C++11 support
    if not config.CheckFlag("-std=c++11"):
        sys.stderr.write("Error: compiler doesn't seem to support the C++11 standard. Abort!\n")
        sys.exit(1) #TODO: exist undefined, fix
    config.env.AppendUnique(CPPFLAGS="-std=c++11")

#    if "-msse3" in config.env['CPPFLAGS'] or "-avx" in config.env['CPPFLAGS'] or env['TARGETCPU'] == "MIC":
#      print "architecture set from outside"
#    else:
#      if not env['SSE3_FALLBACK']:
#          # check avx support
#          if not config.CheckFlag("-mavx"):
#              sys.stderr.write("Error: compiler doesn't seem to support AVX. Abort! Fallin\n")
#              sys.exit(1)
#          config.env.AppendUnique(CPPFLAGS="-mavx")
#      else:
#          config.env.AppendUnique(CPPFLAGS="-msse3")

    # check whether swig installed
    if not config.CheckExec('doxygen'):
        sys.stderr.write("Warning: doxygen cannot be found.\n  You will not be able to generate the documentation.\n  Check PATH environment variable!\n")
    else:
        print "Using doxygen " + commands.getoutput('doxygen --version')

    # check whether dot installed
    if not config.CheckExec('dot'):
        sys.stderr.write("Warning: dot (Graphviz) cannot be found.\n  The documentation might lack diagrams.\n  Check PATH environment variable!\n")
    else:
        print "Using " + commands.getoutput('dot -V').splitlines()[0]
    
    if config.env['USE_OCL']:
      if 'OCL_INCLUDE_PATH' in config.env['ENV']:
        config.env.AppendUnique(CPPPATH=config.env['ENV']['OCL_INCLUDE_PATH'])
      elif 'OCL_INCLUDE_PATH' in config.env:
        config.env.AppendUnique(CPPPATH=config.env['OCL_INCLUDE_PATH'])
      else:
        sys.stderr.write("Info: Trying to find the OpenCL without the variable \"OCL_INCLUDE_PATH\"\n")
         
      if not config.CheckCXXHeader('CL/cl.h'):
        sys.stderr.write("Error: \"CL/cl.h\" not found, but required for OpenCL\n")
        sys.exit(1)
        
      if 'OCL_LIBRARY_PATH' in config.env['ENV']:
        config.env.AppendUnique(LIBPATH=config.env['ENV']['OCL_LIBRARY_PATH'])
      elif 'OCL_LIBRARY_PATH' in config.env:
        config.env.AppendUnique(LIBPATH=config.env['OCL_LIBRARY_PATH'])
      else:
        sys.stderr.write("Info: Trying to find the OpenCL library \"libOpenCL\" without the variable \"OCL_LIBRARY_PATH\"\n")
        
      if not config.CheckLib('OpenCL', language="C++"):
        sys.stderr.write("Error: \"libOpenCL\" not found, but required for OpenCL\n")
        sys.exit(1)
        
      config.env.AppendUnique(CPPDEFINES="USE_OCL")
    else:
      print "Info: OpenCL is not enabled"
    
    # Check the availability of the boost unit test dependencies
    if env['COMPILE_BOOST_TESTS']:
        if not config.CheckHeader("boost/test/unit_test.hpp", language="c++"):
            env['COMPILE_BOOST_TESTS'] = False
            print """****************************************************
No Boost Unit Test Headers found. Omitting Boost unit tests. 
Please install the corresponding package, e.g. using command on Ubuntu
> sudo apt-get install libboost-test-dev
****************************************************
"""

        if not config.CheckLib('boost_unit_test_framework', autoadd=0):
            env['COMPILE_BOOST_TESTS'] = False
            print """****************************************************
No Boost Unit Test library found. Omitting Boost unit tests. 
Please install the corresponding package, e.g. using command on Ubuntu
> sudo apt-get install libboost-test-dev
****************************************************
"""
    
    if env["SG_PYTHON"]:
        # check whether swig installed
        if not config.CheckExec('swig'):
            sys.stderr.write("Error: swig cannot be found, but required for SG_PYTHON. Check PATH environment variable!\n")
            sys.exit(1)

        print "Using SWIG " + re.findall(r"[0-9.]*[0-9]+", commands.getoutput('swig -version'))[0]
        config.env.AppendUnique(CPPPATH=[distutils.sysconfig.get_python_inc()])
        print "pythonpath: ", distutils.sysconfig.get_python_inc()

        if not config.CheckCXXHeader('Python.h'):
            sys.stderr.write("Error: Python.h not found, but required for SG_PYTHON. Check path to Python include files: "
                         + distutils.sysconfig.get_python_inc() + "\n")
            sys.stderr.write("Hint: You might have to install package python-dev\n")
            sys.exit(1)

        if not config.CheckCXXHeader('pyconfig.h'):
            sys.stderr.write("Error: pyconfig.h not found, but required for SG_PYTHON. Check path to Python include files: "
                         + distutils.sysconfig.get_python_inc() + "\n")
            sys.stderr.write("Hint: You might have to install package python-dev\n")
            sys.exit(1)

        try:
            import numpy
            numpy_path = os.path.join(os.path.split(numpy.__file__)[0], "core", "include")
            config.env.AppendUnique(CPPPATH=[numpy_path])
            if not config.CheckCXXHeader(['Python.h', 'pyconfig.h', 'numpy/arrayobject.h']):
                sys.stderr.write('Warning: Cannot find NumPy header files in: "' + str(numpy_path) + '", disabling unit tests\n')
                env['NO_UNIT_TESTS'] = True
        except:
            sys.stderr.write('Warning: Numpy doesn\'t seem to be installed. Disabling unit tests\n')
            env['NO_UNIT_TESTS'] = True
    else:
        print 'Warning: Python extension ("SG_PYTHON") not enabled, skipping python unit tests'

    if env['SG_JAVA']:
        # check for $JAVA_HOME; prepend to search path
        if os.environ.get('JAVA_HOME'):
            config.env.PrependENVPath('PATH', os.path.join(os.environ.get('JAVA_HOME'), 'bin'))
        # check whether javac installed
        if not config.CheckExec('javac'):
            sys.stderr.write("Error: javac cannot be found, but required by SG_JAVA. Check PATH environment variable!\n")
            sys.exit(1)
        # check whether javac installed
        if not config.CheckExec('java'):
            sys.stderr.write("Warning: java cannot be found, but required by SG_JAVA. Check PATH environment variable!\n")
            sys.exit(1)

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
                sys.exit(1)
    else:
        print "Info: Compiling without java support"

    # now set up all further environment settings that should never fail
    # compiler setup should be always after checking headers and flags, as they can make the checks invalid
    # e.g. by setting "-Werror"

    # TODO check
    if env['OPT'] == True:
       env.Append(CPPFLAGS=['-O3'])
    else:
       env.Append(CPPFLAGS=['-g']) #, '-O0'

    if env['COMPILER'].upper() == 'GNU':
        print "Using default gcc " + commands.getoutput(env['CXX'] + ' -dumpversion')

        allWarnings = "-Wall -pedantic -pedantic-errors -Wextra \
            -Wcast-align -Wcast-qual -Wconversion -Wdisabled-optimization -Wformat=2 \
            -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport  -Winit-self  \
            -Winvalid-pch -Wmissing-field-initializers -Wmissing-format-attribute -Wmissing-include-dirs \
             -Wpacked   -Wpointer-arith -Wredundant-decls -Wstack-protector \
            -Wstrict-aliasing=2 -Wswitch-default -Wswitch-enum -Wunreachable-code -Wunused -Wunused-parameter \
            -Wvariadic-macros -Wwrite-strings -Wuninitialized".split(" ")
        # rather uninteresting: -Wlong-long -Wpadded -Wshadow -Wfloat-equal -Waggregate-return -Wimplicit -Wmissing-noreturn -Weffc++
        # cannot really use:

        # -fno-strict-aliasing: http://www.swig.org/Doc1.3/Java.html or http://www.swig.org/Release/CHANGES, 03/02/2006
        #    "If you are going to use optimisations turned on with gcc > 4.0 (for example -O2),
        #     ensure you also compile with -fno-strict-aliasing"
        env.Append(CPPFLAGS=allWarnings + [
                             # '-Wall', '-Wextra',
                             #'-std=c++11',  # '-Wno-long-long', '-Wno-deprecated',
                             # '-Werror',
                             '-Wno-unused-parameter',
                             # '-Wconversion',
                             '-fno-strict-aliasing',
                             '-funroll-loops', '-mfpmath=sse',
                             '-DDEFAULT_RES_THRESHOLD=-1.0', '-DTASKS_PARALLEL_UPDOWN=4'])
        env.Append(CPPFLAGS=['-fopenmp'])
        env.Append(LINKFLAGS=['-fopenmp'])
        
        if env['ARCH'].upper() == 'SSE3':
            config.env.AppendUnique(CPPFLAGS="-msse3")
        elif env['ARCH'].upper() == 'SSE42':
            config.env.AppendUnique(CPPFLAGS="-msse4.2")
        elif env['ARCH'].upper() == 'AVX':
            config.env.AppendUnique(CPPFLAGS="-mavx")
        elif env['ARCH'].upper() == 'FMA4':
            config.env.AppendUnique(CPPFLAGS="-mavx")
            config.env.AppendUnique(CPPFLAGS="-mfma4")
        elif env['ARCH'].upper() == 'AVX2':
            config.env.AppendUnique(CPPFLAGS="-mavx2")
            config.env.AppendUnique(CPPFLAGS="-mfma")
        elif env['ARCH'].upper() == 'AVX512':
            config.env.AppendUnique(CPPFLAGS="-mavx512f")
            config.env.AppendUnique(CPPFLAGS="-mavx512cd")
            config.env.AppendUnique(CPPFLAGS="-mfma")
        else:
            print "You must specify a valid ARCH value for gnu."
            print "Available configurations are: sse3, sse4.2, avx, fma4, avx2, avx512"
            sys.exit(1)            
    elif env['COMPILER'].upper() == 'INTEL':
        print "Using icc"
        env.Append(CPPFLAGS=['-Wall', '-ansi', '-Wno-deprecated', '-wd1125',
                               '-fno-strict-aliasing',
                               '-ip', '-ipo', '-funroll-loops',
                               '-ansi-alias', '-fp-speculation=safe',
                               '-DDEFAULT_RES_THRESHOLD=-1.0', '-DTASKS_PARALLEL_UPDOWN=4', '-no-offload'])

        env['CC'] = ('icc')
        env['LINK'] = ('icpc')
        env['CXX'] = ('icpc')
        env.Append(CPPFLAGS=['-openmp'])
        env.Append(LINKFLAGS=['-openmp'])

        if env['ARCH'].upper() == 'SSE3':
            config.env.AppendUnique(CPPFLAGS="-msse3")
        elif env['ARCH'].upper() == 'SSE42':
            config.env.AppendUnique(CPPFLAGS="-msse4.2")
        elif env['ARCH'].upper() == 'AVX':
            config.env.AppendUnique(CPPFLAGS="-mavx")
        elif env['ARCH'].upper() == 'AVX2':
            config.env.AppendUnique(CPPFLAGS="-xCORE-AVX2")
            config.env.AppendUnique(CPPFLAGS="-fma")
        elif env['ARCH'].upper() == 'AVX512':
            config.env.AppendUnique(CPPFLAGS="-xCOMMON-AVX512")
            config.env.AppendUnique(CPPFLAGS="-fma")
        elif env['ARCH'].upper() == 'MIC':
            config.env.AppendUnique(CPPFLAGS="-mmic")
            config.env.Append(LINKFLAGS=['-mmic'])
        else:
            print "You must specify a valid ARCH value for intel."
            print "Available configurations are: sse3, sse4.2, avx, avx2, avx512, mic"
            sys.exit(1)  

        env.Append(CPPPATH=[distutils.sysconfig.get_python_inc()])  
    else:
        print "You must specify a valid value for Compiler."
        print "Available configurations are: gnu and intel"
        sys.exit(1)

    # special treatment for different platforms
    if env['PLATFORM'] == 'darwin':
        # the "-undefined dynamic_lookup"-switch is required to actually build a shared library
        # in OSX. "-dynamiclib" alone results in static linking of all further dependent shared libraries
        # beware: if symbols are missing that are actually required (because the symbols don't reside in a shared library), there will be no error during compilation
        # the python binding (pysgpp) requires lpython and a flat namespace
        # also for the python binding, the library must be suffixed with '*.so' even though it is a dynamiclib and not a bundle (see SConscript in src/pysgpp)
        env.Append(LINKFLAGS=['-flat_namespace', '-undefined', 'dynamic_lookup', '-lpython'])
        #The GNU assembler (GAS) is not supported in Mac OS X. A solution that fixed this problem is by adding -Wa,-q to the compiler flags.
        #From the man pages for as (version 1.38): -q Use the clang(1) integrated assembler instead of the GNU based system assembler.
        #Note that the CPPFLAG is exactly "-Wa,-q", where -Wa passes flags to the assembler and -q is the relevant flag to make it use integrated assembler
        env.Append(CPPFLAGS=['-Wa,-q'])
        env.AppendUnique(CPPPATH="/usr/local/include")
        env.AppendUnique(LIBPATH="/usr/local/lib")
        env['SHLIBSUFFIX'] = '.dylib'
    elif env['PLATFORM'] == 'cygwin':
        # required to find the static libraries compiled before the shared libraries
        # the static libraries are required as the linker on windows cannot ignore undefined symbols
        # (as is done on linux automatically and is done on OSX with the settings above)
        #env.Append(LIBPATH=[BUILD_DIR])
        pass

    # will lead to a warning on cygwin (and we have -Werror enabled)
    # is enabled by default on cygwin
    if env['PLATFORM'] != 'cygwin':
        env.Append(CPPFLAGS=['-fPIC'])

    # setup the include base folder
    # env.Append(CPPPATH=['#/src/sgpp'])
    for moduleFolder in moduleFolders:
      if moduleFolder in languageWrapperFolders:
        continue
      env.Append(CPPPATH=['#/' + moduleFolder + '/src/'])

    # detour compiler output
    env['PRINT_CMD_LINE_FUNC'] = Helper.print_cmd_line

    env = config.Finish()

    # clear build_log file
    with open(env['CMD_LOGFILE'], 'a') as logFile:
        logFile.seek(0)
        logFile.truncate()
