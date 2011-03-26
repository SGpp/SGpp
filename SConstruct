# Copyright (C) 2009 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp

# author Dirk Pflueger (Dirk.Pflueger@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)


import os, sys
import distutils.sysconfig

# Check for versions of Scons and Python
EnsureSConsVersion(1, 0)
EnsurePythonVersion(2, 5)

# Custom test for executables used during configuration
def CheckExec(context, cmd):
    context.Message( 'Checking for %s...' % (cmd) )
    ret = context.env.WhereIs(cmd)
    if ret == None:
        ret = ''
    context.Result(ret)
    return ret

# get all subdirs of path, required by CheckJNI
def getSubdirs(path):
    pathlist = []
    for f in os.listdir(path):
        if os.path.isdir(os.path.join(path, f)):
            pathlist.append(os.path.join(path, f))
    return pathlist

# Check for jni header file
# if found, additionally add all subdirs to CPPPATH (platform dependent files)
def CheckJNI(context):
    found = False
    print "Trying to locate jni.h..."
    # message if JNI_CPPINCLUDE not set
    if not os.environ.get('JNI_CPPINCLUDE'):
        print "... JNI_CPPINCLUDE not set"
    # check for JAVA_HOME first
    if os.environ.get('JAVA_HOME'):
        pname = os.path.join(os.environ.get('JAVA_HOME'), 'include')
        if os.path.exists(os.path.join(pname, 'jni.h')):
            context.env.Append(CPPPATH = [pname]+getSubdirs(pname))
            res = "... found in "+pname
            context.Result(res)
            return res
        else:
            print "... not found in $JAVA_HOME/include"
    else:
        print "... JAVA_HOME not set"
        # not found, try guessing:
        # look, where java and javac are located: include/ directory might be 1 or 2 dirs below
        print "... trying to guess"
        for f in ['java', 'javacc']:
            fdir = context.env.WhereIs(f)
            if not fdir:
                continue
            # os.path.realpath to resolve links
            basedir = os.path.dirname(os.path.realpath(fdir))
            for subdir in ['..', os.path.join('..','..')]:
                pname = os.path.join(basedir, subdir, 'include')
                if os.path.exists(os.path.join(pname, 'jni.h')):
                    context.env.Append(CPPPATH = [pname]+getSubdirs(pname))
                    res = "... found in "+pname
                    context.Result(res)
                    return res
    ret = 0
    context.Result('... nothing found!')
    return ret






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
#vars.Add('JNI_CPPPATH', 'Path to JNI includes', None)
#vars.Add('JNI_OS', 'JNI os path', None)


# for compiling on LRZ without errors: omit unit tests
vars.Add('NO_UNIT_TESTS', 'Omit UnitTests if set to True', False)


env = Environment(variables = vars, ENV = os.environ)


# Help Text
Help("""---------------------------------------------------------------------

Type: 'scons [parameters]' to build the libraries

There are compiler optimizations for different platforms which can be
specified via parameters.

Parameters can be set either by setting the corresponding environment
variables, or directly via the commandline, e.g.,
> scons OMP=True
to enable OpenMP support.


Specifying the target, the following options are available:
    - default: using the gcc toolchain with OpenMP 2
    - opteronICC: using the ICC 11.x toolchain with OpenMP 3 with standard x86_64 options
    - core2ICC: using the ICC 11.x toolchain with OpenMP 3 with Intel x86_64 options (core architecture)
    - nehalemICC: using the ICC 11.x toolchain with OpenMP 3 with Intel x86_64 options (nehalem architecture)
    - snbICC: using the ICC 12.x toolchain with OpenMP 3 with Intel x86_64 options (sandy bridge architecture)
    - ia64ICC: using the ICC 11.x toolchain with OpenMP 3 with Itanium options

For LRZ, please execute:
module load python
module load gcc/4.5

FOR LRZ and when using intel compiler, execute:
export LIBPATH=$LD_LIBRARY_PATH

---------------------------------------------------------------------

Parameters are:
""" +
vars.GenerateHelpText(env))





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
elif env['TARGETCPU'] == 'snbICC':
    print "Using icc 12.0 for Sandy Bridge systems"
    env.Append(CPPFLAGS = ['-axAVX', '-O3', '-funroll-loops', '-ipo', '-ip', '-ansi-alias', 
                           '-Wall', '-ansi', '-wd981', 
                           '-fno-strict-aliasing'])
else:
    print "You must specify a valid value for TARGETCPU."
    print "Available configurations are: default, core2ICC, opteronICC, ia64ICC"
    Exit(1)
    
# sets ICC-wide commen options and the tool chain   
if env['TARGETCPU'] in ['ia64ICC', 'opteronICC', 'core2ICC', 'nehalemICC', 'snbICC']:
    env['CC'] = ('icc')
    env['LINK'] = ('icpc')
    env['CXX'] = ('icpc')	    
    if env['OMP']:
	env.Append(CPPFLAGS=['-openmp'])
        env.Append(LINKFLAGS=['-openmp']) 
    
# sets the architecture option for gcc
if env.has_key('MARCH'):
    if env['TARGETCPU'] == 'default':
        env.Append(CPPFLAGS=('-march=' + env['MARCH']))
    else:
        print "Warning: Ignoring option MARCH"
         


# boolean variables for environment 
pyAvail = True
swigAvail = True
javaAvail = True

# configure environment
# ---------------------
if not env.GetOption('clean'):

    config = env.Configure(custom_tests = { 'CheckExec' : CheckExec,
                                            'CheckJNI' : CheckJNI })

    # check if the math header is available
    if not config.CheckLibWithHeader('m', 'math.h', 'c++'):
        sys.stderr.write("Error: Math headers are missing.\n")
        Exit(1)

    # check whether swig installed
    if not config.CheckExec('swig'):
        sys.stderr.write("Error: swig cannot be found. Check PATH environment variable!\n")
        swigAvail = False

    # check for Python headers
    config.env.AppendUnique(CPPPATH = distutils.sysconfig.get_python_inc())
    if not config.CheckCXXHeader('Python.h'):
        sys.stderr.write("Error: Python.h not found. Check path to Python include files: "
                         + distutils.sysconfig.get_python_inc() + "\n")
        sys.stderr.write("Warning: You might have to install package python-dev\n")
        sys.stderr.write("... skipping Python support and unit tests")
        pyAvail = False

    # check for $JAVA_HOME; prepend to search path
    if os.environ.get('JAVA_HOME'):
        config.env.PrependENVPath('PATH', os.path.join(os.environ.get('JAVA_HOME'), 'bin'))

    # check whether javac installed
    if not config.CheckExec('javac'):
        sys.stderr.write("Error: javac cannot be found. Check PATH environment variable!\n")
        javaAvail = False
    # check whether javac installed
    if javaAvail and not config.CheckExec('java'):
        sys.stderr.write("Warning: java cannot be found. Check PATH environment variable!\n")

    # check for JNI headers
    if javaAvail and os.environ.get('JNI_CPPINCLUDE'):
        config.env.AppendUnique(CPPPATH = [os.environ.get('JNI_CPPINCLUDE')])
    if javaAvail and not config.CheckCXXHeader('jni.h'):
        # not found; try to find
        if not config.CheckJNI():
            sys.stderr.write("Error: jni.h not found.\n"
                             +"Please set JAVA_HOME environment variable "
                             +"with $JAVA_HOME/bin/javac, $JAVA_HOME/include/jni.h\n"
                             +"or directly $JNI_CPPINCLUDE with $JNI_CPPINCLUDE/jni.h\n")
            javaAvail = False
    if not javaAvail:
        sys.stderr.write("No Java support...\n")

    # check if the intel omp lib is available
    if env['TARGETCPU'] in ['ia64ICC', 'opteronICC', 'core2ICC', 'nehalemICC', 'snbICC'] and env['OMP']:
        if not config.CheckLib('iomp5'):
            print "Error: Intel omp library iomp5 is missing."
            Exit(1)
               
    # check if the the intel vector lib is available
    if env['TARGETCPU'] in ['ia64ICC', 'opteronICC', 'core2ICC', 'nehalemICC', 'snbICC']:
        if not config.CheckLib('svml'):
            print "SVML should be available when using intelc. Consider runnning scons --config=force!"
            
    env = config.Finish()


# the optional CPPFLAGS at the end will override the previous flags
env['CPPFLAGS'] = env['CPPFLAGS'] + opt_flags

Export('env')

print "finished configuration"


# Now start for compilation...
# ----------------------------

# build c++ lib
(libsgpp, libsgppa) = env.SConscript('src/sgpp/SConscript',
                                     build_dir='tmp/build_sg', duplicate=0)
# install
env.Install('lib/sgpp', [libsgpp, libsgppa])

# build python lib
if swigAvail and pyAvail:
    libpysgpp = env.SConscript('src/pysgpp/SConscript',
                               build_dir='tmp/build_pysgpp', duplicate=0)
    # install
    pyinst = env.Install('lib/pysgpp', [libpysgpp, 'tmp/build_pysgpp/pysgpp.py'])
    Depends(pyinst, libpysgpp)
    dep = env.Install('bin', [libpysgpp, 'tmp/build_pysgpp/pysgpp.py'])
    Depends(dep, libpysgpp)
    
# build java lib
if swigAvail and javaAvail and env['JSGPP']:
    libjsgpp = env.SConscript('src/jsgpp/SConscript',
                              build_dir='tmp/build_jsgpp', duplicate=0)
    libweka = env.SConscript('src/jsgpp_weka/SConscript',
                             build_dir='tmp/build_jsgpp_weka', duplicate=0)
    # install
    jinst = env.Install('lib/jsgpp', [libjsgpp])

# Execute Unit Tests
if not env['NO_UNIT_TESTS'] and pyAvail:
    dep = env.SConscript('tests/SConscript')
    # execute after all installations (even where not necessary)
    if javaAvail and env['JSGPP']:
        Depends(dep, [jinst, pyinst])
    else:
        Depends(dep, [pyinst])
else:
    sys.stderr.write("Warning!! Skipping unit tests!!\n")

