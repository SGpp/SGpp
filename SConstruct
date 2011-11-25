# Copyright (C) 2009 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp

# author Dirk Pflueger (Dirk.Pflueger@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)


import os, sys
import distutils.sysconfig

# Check for versions of Scons and Python
EnsureSConsVersion(1, 0)
EnsurePythonVersion(2, 5)

# Definitions and functions
#########################################################################

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
        # look, where java and javac are located:
        # include/ directory might be 1 or 2 dirs below
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


# Definition of flags / command line parameters for SCons
#########################################################################

vars = Variables("custom.py")

# define the flags 
vars.Add('CPPFLAGS','Set additional Flags, they are compiler-depended','-Wno-deprecated')
vars.Add('LINKFLAGS','Set additional Linker-flags, they are linker-depended','')

# define the target
vars.Add('MARCH','Sets the architecture if compiling with gcc, this is a pass-through option: just specify the gcc options!', None)
vars.Add('TARGETCPU',"Sets the processor you are compiling for. 'default' means using gcc with standard configuration. Also available are: 'opteronICC', 'core2ICC', 'ia64ICC'; here Intel Compiler in version 11 must be used", 'default')
vars.Add('OMP', "Sets if OpenMP should be used; with gcc OpenMP 2 is used, with all icc configurations OpenMP 3 is used!", False)
vars.Add('TRONE', "Sets if the tr1/unordered_map should be uesed", False)

# for compiling on LRZ without errors: omit unit tests
vars.Add('NO_UNIT_TESTS', 'Omit UnitTests if set to True', False)

# for compiling different modules
vars.Add('SG_ALL', 'Build all modules', False)
vars.Add('SG_BASE', 'Build Basis Module', False)
vars.Add('SG_DATADRIVEN', 'Build Datadriven Module', False)
vars.Add('SG_SOLVER', 'Build Solver Module', False)
vars.Add('SG_FINANCE', 'Build Finance Module', False)
vars.Add('SG_PDE', 'Build PDE Module', False)
vars.Add('SG_PARALLEL', 'Build Parallel Module', False)
vars.Add('SG_COMBIGRID', 'Build Combigrid Module', False)
vars.Add('SG_PYTHON', 'Build Python Support', False)
vars.Add('SG_JAVA', 'Build Java Support', False)
# modules and dependencies
moduleList = {'SG_BASE': (), 
              'SG_DATADRIVEN': ('SG_BASE'), 
              'SG_SOLVER': ('SG_BASE'), 
              'SG_FINANCE': ('SG_BASE', 'SG_PDE'),
              'SG_PDE': ('SG_BASE'), 
              'SG_PARALLEL': ('SG_BASE', 'SG_PDE', 'SG_DATADRIVEN', 'SG_FINANCE'), 
              'SG_COMBIGRID': ('SG_BASE')}
supportList = ['SG_PYTHON', 'SG_JAVA']

# initialize environment
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



# Set compiler switches and check architectures
#########################################################################

# scons usually adds double quotes around the command-line arguments containing 
# white spaces. As this whould produce compilation error, replace string 
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
                         '-funroll-loops', '-mfpmath=sse', '-msse3'])
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
        env.Append(CPPDEFINES=['USEOMP', 'USEOMPTHREE', 'USEICCINTRINSICS'])
    
# sets the architecture option for gcc
if env.has_key('MARCH'):
    if env['TARGETCPU'] == 'default':
        env.Append(CPPFLAGS=('-march=' + env['MARCH']))
    else:
        print "Warning: Ignoring option MARCH"

# special treatment for MAC OS-X
if env['PLATFORM']=='darwin':
    env.Append(LINKFLAGS=['-lpython'])
    env['SHLIBSUFFIX'] = '.so'


# the optional CPPFLAGS at the end will override the previous flags
env['CPPFLAGS'] = env['CPPFLAGS'] + opt_flags



# Decide what to compile
#########################################################################

# for clean enable everything:
if env.GetOption('clean'):
    for entry in moduleList.keys()+supportList:
        env[entry] = True

# if neither module nor support language set, do all
anySet = False
for entry in moduleList.keys() + supportList:
    if env[entry]:
        anySet = True
if not anySet:
    env['SG_ALL'] = True

# SG_ALL activates all modules and Python 
if env['SG_ALL']:
    print "Compiling all modules..."
    print "Compiling all support..."
    for entry in moduleList.keys() + supportList:
        env[entry] = True

# check dependencies
for modl in moduleList.keys():
    if env[modl]:
        for dep in moduleList[modl]:
            env[dep] = True
for modl in moduleList.keys():
    if env[modl]:
        print "Compiling module", modl
# support for non-C++
for sup in supportList:
    if env[sup]:
        print "Compiling support for", sup

# add C++ defines for all modules
cppdefines = []
for modl in moduleList.keys():
    if env[modl]:
        cppdefines.append(modl)
env.Append(CPPDEFINES=cppdefines)


# Initialize environment + support for Python and Java
#########################################################################
# boolean variables for environment 
pyAvail = True
swigAvail = True
javaAvail = True

# no checks if clean:
if not env.GetOption('clean'):
    print """
******************************************
* Configuring system                     *
******************************************"""

    config = env.Configure(custom_tests = { 'CheckExec' : CheckExec,
                                            'CheckJNI' : CheckJNI })
    # check whether swig installed
    if not config.CheckExec('doxygen'):
        sys.stderr.write("Warning: doxygen cannot be found.\n  You will not be able to generate the documentation.\n  Check PATH environment variable!\n")

    # check whether dot installed
    if not config.CheckExec('dot'):
        sys.stderr.write("Warning: dot (Graphviz) cannot be found.\n  The documentation might lack diagrams.\n  Check PATH environment variable!\n")

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

    print """******************************************
* Finished configuring system            *
******************************************"""

# End of configuration
#########################################################################
Export('env')


# Now compile
#########################################################################
lib_sgpp_targets = []

if env['SG_BASE']:
	SConscript('src/sgpp/SConscriptBase', build_dir='tmp/build_sgbase', duplicate=0)
	Import('libsgppbase')
	Import('libsgppbasestatic')
	lib_sgpp_targets.append(libsgppbase)
	lib_sgpp_targets.append(libsgppbasestatic)
	
if env['SG_PDE']:
	SConscript('src/sgpp/SConscriptPde', build_dir='tmp/build_sgpde', duplicate=0)
	Import('libsgpppde')
	Import('libsgpppdestatic')
	lib_sgpp_targets.append(libsgpppde)
	lib_sgpp_targets.append(libsgpppdestatic)
	
if env['SG_DATADRIVEN']:
	SConscript('src/sgpp/SConscriptDatadriven', build_dir='tmp/build_sgdatadriven', duplicate=0)
	Import('libsgppdatadriven')
	Import('libsgppdatadrivenstatic')
	lib_sgpp_targets.append(libsgppdatadriven)
	lib_sgpp_targets.append(libsgppdatadrivenstatic)
	
if env['SG_SOLVER']:
	SConscript('src/sgpp/SConscriptSolver', build_dir='tmp/build_sgsolver', duplicate=0)
	Import('libsgppsolver')
	Import('libsgppsolverstatic')
	lib_sgpp_targets.append(libsgppsolver)
	lib_sgpp_targets.append(libsgppsolverstatic)
	
if env['SG_FINANCE']:
	SConscript('src/sgpp/SConscriptFinance', build_dir='tmp/build_sgfinance', duplicate=0)
	Import('libsgppfinance')
	Import('libsgppfinancestatic')
	lib_sgpp_targets.append(libsgppfinance)
	lib_sgpp_targets.append(libsgppfinancestatic)
	
if env['SG_PARALLEL']:
	SConscript('src/sgpp/SConscriptParallel', build_dir='tmp/build_sgparallel', duplicate=0)
	Import('libsgppparallel')
	Import('libsgppparallelstatic')
	lib_sgpp_targets.append(libsgppparallel)
	lib_sgpp_targets.append(libsgppparallelstatic)

if env['SG_COMBIGRID']:
	SConscript('src/sgpp/SConscriptCombigrid', build_dir='tmp/build_sgcombigrid', duplicate=0)
	Import('libsgppcombigrid')
	Import('libsgppcombigridstatic')
	lib_sgpp_targets.append(libsgppcombigrid)
	lib_sgpp_targets.append(libsgppcombigridstatic)

# build python lib
if env['SG_PYTHON'] and swigAvail and pyAvail:
	libpysgpp = SConscript('src/pysgpp/SConscript', build_dir='tmp/build_pysgpp', duplicate=0)
	pyinst = env.Install('lib/pysgpp', [libpysgpp, 'tmp/build_pysgpp/pysgpp.py'])
	Depends(pyinst, libpysgpp)
	pybin = env.Install('bin', [libpysgpp, 'tmp/build_pysgpp/pysgpp.py'])
	Depends(pybin, libpysgpp)
    
# build java lib
if swigAvail and javaAvail and env['SG_JAVA']:
    libjsgpp = env.SConscript('src/jsgpp/SConscript',
                              build_dir='tmp/build_jsgpp', duplicate=0)
#    libweka = env.SConscript('src/jsgpp_weka/SConscript',
#                             build_dir='tmp/build_jsgpp_weka', duplicate=0)
    # install
    jinst = env.Install('lib/jsgpp', [libjsgpp])
	
env.Install('#lib/sgpp', lib_sgpp_targets)


# Unit tests
#########################################################################

if not env['NO_UNIT_TESTS'] and env['SG_PYTHON'] and pyAvail and swigAvail:
    testdep = env.SConscript('tests/SConscript')
    # execute after all installations (even where not necessary)
    if javaAvail and env['SG_JAVA']:
        Depends(testdep, [jinst, pyinst])
    else:
        Depends(testdep, [pyinst])
else:
    sys.stderr.write("Warning!! Skipping unit tests!!\n\n\n")

