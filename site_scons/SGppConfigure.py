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

def doConfigure(env, moduleFolders):
    print "Checking programs and libraries: "

    config = env.Configure(custom_tests={ 'CheckExec' : SGppConfigureExtend.CheckExec,
                                            'CheckJNI' : SGppConfigureExtend.CheckJNI,
                                            'CheckFlag' : SGppConfigureExtend.CheckFlag })


    # check C++11 support
    if not config.CheckFlag("-std=c++11"):
        sys.stderr.write("Error: compiler doesn't seem to support the C++11 standard. Abort!\n")
        Exit(1)
    config.env.AppendUnique(CPPFLAGS="-std=c++11")

    if not env['SSE3_FALLBACK']:
        # check avx support
        if not config.CheckFlag("-mavx"):
            sys.stderr.write("Error: compiler doesn't seem to support AVX. Abort! Fallin\n")
            Exit(1)
        config.env.AppendUnique(CPPFLAGS="-mavx")
    else:
        config.env.AppendUnique(CPPFLAGS="-msse3")

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

    if env["SG_PYTHON"]:
        # check whether swig installed
        if not config.CheckExec('swig'):
            sys.stderr.write("Error: swig cannot be found, but required for SG_PYTHON. Check PATH environment variable!\n")
            Exit(1)

        print "Using SWIG " + re.findall(r"[0-9.]*[0-9]+", commands.getoutput('swig -version'))[0]
        config.env.AppendUnique(CPPPATH=[distutils.sysconfig.get_python_inc()])
        print "pythonpath: ", distutils.sysconfig.get_python_inc()

        if not config.CheckCXXHeader('Python.h'):
            sys.stderr.write("Error: Python.h not found, but required for SG_PYTHON. Check path to Python include files: "
                         + distutils.sysconfig.get_python_inc() + "\n")
            sys.stderr.write("Hint: You might have to install package python-dev\n")
            Exit(1)

        if not config.CheckCXXHeader('pyconfig.h'):
            sys.stderr.write("Error: pyconfig.h not found, but required for SG_PYTHON. Check path to Python include files: "
                         + distutils.sysconfig.get_python_inc() + "\n")
            sys.stderr.write("Hint: You might have to install package python-dev\n")
            Exit(1)

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
            Exit(1)
        # check whether javac installed
        if not config.CheckExec('java'):
            sys.stderr.write("Warning: java cannot be found, but required by SG_JAVA. Check PATH environment variable!\n")
            Exit(1)

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
                Exit(1)
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
    # env.Append(CPPPATH=['#/src/sgpp'])
    for moduleFolder in moduleFolders:
        env.Append(CPPPATH=['#/' + moduleFolder + '/src/'])

    # detour compiler output
    env['PRINT_CMD_LINE_FUNC'] = Helper.print_cmd_line

    env = config.Finish()

    # clear build_log file
    logfile = open(env['CMD_LOGFILE'], 'a')
    logfile.seek(0)
    logfile.truncate()
