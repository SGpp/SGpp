# Copyright (C) 2009 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp

# author Dirk Pflueger (Dirk.Pflueger@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), David Pfander (David.Pfander@ipvs.uni-stuttgart.de)


import os, sys
import distutils.sysconfig
import glob
import SCons

# Check for versions of Scons and Python
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
            
    context.Result('... nothing found!')
    return 0

# get all files in a folder matching "SConscript*"
def getModules(path):
    modules = glob.glob(path + '/SConscript*')
    for i in range(len(modules)):            
        modules[i] = modules[i].split('SConscript', 1)[1]
    
    # filter out backup files        
    for module in modules:
        if '~' in module:
            modules.remove(module)
        
    return modules

# Definition of flags / command line parameters for SCons
#########################################################################

def multiParamConverter(s):
    print s
    return s.split(',')

vars = Variables("custom.py")

# define the flags 
vars.Add('CPPFLAGS','Set additional Flags, they are compiler-depended (multiple flags combined with comma, e.g. -lpython,-lm)', '', converter=multiParamConverter)
vars.Add('LINKFLAGS','Set additional Linker-flags, they are linker-depended (multiple flags combined with comma, e.g. -lpython,-lm)', '', converter=multiParamConverter)

# define the target
vars.Add('MARCH','Sets the architecture if compiling with gcc, this is a pass-through option: just specify the gcc options!', None)
vars.Add('TARGETCPU',"Sets the processor you are compiling for. 'default' means using gcc with standard configuration. Also available are: 'ICC', here Intel Compiler in version 11 or higher must be used", 'default')
vars.Add(BoolVariable('OMP', "Sets if OpenMP should be used; with gcc OpenMP 2 is used, with all icc configurations OpenMP 3 is used!", False))
vars.Add(BoolVariable('TRONE', "Sets if the tr1/unordered_map should be uesed", False))
vars.Add('OPT', "Sets optimization on and off", False)

# for compiling on LRZ without errors: omit unit tests
vars.Add(BoolVariable('NO_UNIT_TESTS', 'Omit UnitTests if set to True', False))

# modules and dependencies
moduleList = {}
src_files = {}
supportList = ['SG_PYTHON', 'SG_JAVA']

# find all modules
modules = getModules('src/sgpp')

# check dependencies and import src file locations
for name in modules:
    SConscript('src/sgpp/SConscript' + name, variant_dir='tmp/build_sg' + name.lower(), duplicate=0)
    Import('srcs')
    Import('dependencies')
    moduleList['SG_' + name.upper()] = dependencies
    print 'Module SG_' + name.upper() + ' depends on:'
    for dep in dependencies:
        print '\t' + dep
        if not dep in modules:
            print "Error!"
            print name + " depends on non-existent module " + dep
            Exit(1)
    src_files[name.lower()] = srcs

# for compiling different modules
for module in moduleList:
    vars.Add(BoolVariable(module, 'Build  Module: ' + module, False))

vars.Add(BoolVariable('SG_ALL', 'Build all modules', False))
vars.Add(BoolVariable('SG_PYTHON', 'Build Python Support', False))
vars.Add(BoolVariable('SG_JAVA', 'Build Java Support', False))
vars.Add('OUTPUT_PATH', 'Path where built libraries are installed. Needs a trailing slash!', '')


# verbosity options
vars.Add(BoolVariable('VERBOSE', 'Set output verbosity', False))
vars.Add('CMD_LOGFILE','Specifies a file to capture the build log','build_log.txt')

# initialize environment
env = Environment(variables = vars, ENV = os.environ)

# sanity check in case user didn't read the variable description
if not env['OUTPUT_PATH'] == '':
    if not env['OUTPUT_PATH'][-1] == '/':
        env['OUTPUT_PATH'] = env['OUTPUT_PATH'] + '/'

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

# clear build_log file
logfile = open(env['CMD_LOGFILE'], 'a')
logfile.seek(0)
logfile.truncate()

# detour compiler output
def print_cmd_line(s, target, src, env):
    if env['VERBOSE']:
        sys.stdout.write(u'%s\n'%s)
    else:
        sys.stdout.write(u'.')
        sys.stdout.flush()
    if env['CMD_LOGFILE']:
        open(env['CMD_LOGFILE'], 'a').write('%s\n'%s);


env['PRINT_CMD_LINE_FUNC'] = print_cmd_line


# Set compiler switches and check architectures
#########################################################################

# scons usually adds double quotes around the command-line arguments containing 
# white spaces. As this whould produce compilation error, replace string 
# with corresponding list of parameters
opt_flags = Split(env['CPPFLAGS'])
env['CPPFLAGS'] = []

if env['TRONE']:
    env.Append(CPPDEFINES=['USETRONE'])
    env.Append(CPPFLAGS=['-std=c++11'])

if env['OPT']:
   env.Append(CPPFLAGS=['-O3'])
else:	
   env.Append(CPPFLAGS=['-g'])	

if env['TARGETCPU'] == 'default':
    print "Using default gcc"
    # -Wno-long-long as swig uses long long
    # -fno-strict-aliasing: http://www.swig.org/Doc1.3/Java.html or http://www.swig.org/Release/CHANGES, 03/02/2006
    #    "If you are going to use optimisations turned on with gcc > 4.0 (for example -O2), 
    #     ensure you also compile with -fno-strict-aliasing"
    env.Append(CPPFLAGS=['-Wall', '-ansi', '-pedantic', '-Wno-long-long', '-Werror', '-Wno-deprecated', 
                         '-fno-strict-aliasing', '-O3', '-Wconversion',
                         '-funroll-loops', '-mfpmath=sse', '-msse3', 
                         '-DDEFAULT_RES_THRESHOLD=-1.0', '-DTASKS_PARALLEL_UPDOWN=4'])
    if env['OMP']:
        env.Append(CPPFLAGS=['-fopenmp'])
        env.Append(LINKFLAGS=['-fopenmp'])
    else:
        # do not stop for unknown pragmas (due to #pragma omp ... )
        env.AppendUnique(CPPFLAGS=['-Wno-unknown-pragmas'])

elif env['TARGETCPU'] == 'ICC':
    print "Using icc"
    env.Append(CPPFLAGS = ['-Wall', '-ansi', '-Werror', '-Wno-deprecated', '-wd1125',  
                           '-fno-strict-aliasing', '-O3',
                           '-ip', '-ipo', '-funroll-loops', '-msse3',
                           '-ansi-alias', '-fp-speculation=safe', 
                           '-DDEFAULT_RES_THRESHOLD=-1.0', '-DTASKS_PARALLEL_UPDOWN=4', '-no-offload'])

    env['CC'] = ('icc')
    env['LINK'] = ('icpc')
    env['CXX'] = ('icpc')        

    if env['OMP']:
        env.Append(CPPFLAGS=['-openmp'])
        env.Append(LINKFLAGS=['-openmp']) 
    else:
        # do not stop for unknown pragmas (due to #pragma omp ... )
        env.AppendUnique(CPPFLAGS=['-Wno-unknown-pragmas'])


else:
    print "You must specify a valid value for TARGETCPU."
    print "Available configurations are: ICC"
    Exit(1)
    
# sets the architecture option for gcc
if env.has_key('MARCH'):
    if env['TARGETCPU'] == 'default':
        env.Append(CPPFLAGS=('-march=' + env['MARCH']))
    else:
        print "Warning: Ignoring option MARCH"

# special treatment for different platforms
if env['PLATFORM']=='darwin':
    # the "-undefined dynamic_lookup"-switch is required to actually build a shared library 
    # in OSX. "-dynamiclib" alone results in static linking of all further dependent shared libraries
    # beware: if symbols are missing that are actually required (because the symbols don't reside in a shared library), there will be no error during compilation
    # the python binding (pysgpp) requires lpython and a flat namespace
    # also for the python binding, the library must be suffixed with '*.so' even though it is a dynamiclib and not a bundle (see SConscript in src/pysgpp)
    env.Append(LINKFLAGS=['-flat_namespace', '-undefined', 'dynamic_lookup', '-lpython'])
    env['SHLIBSUFFIX'] = '.dylib'
elif env['PLATFORM']=='cygwin':
    # required to find the static libraries compiled before the shared libraries
    # the static libraries are required as the linker on windows cannot ignore undefined symbols
    # (as is done on linux automatically and is done on OSX with the settings above)
    env.Append(LIBPATH=['tmp/build'])
    # required because of usage of obsolete "rand_r()"
    # (should be changed to "-std=c++11" after refactoring)
    env.Append(CPPFLAGS=['-std=gnu++03'])

# will lead to a warning on cygwin (and we have -Werror enabled)
# is enabled by default on cygwin
if env['PLATFORM'] != 'cygwin':
    env.Append(CPPFLAGS=['-fPIC'])

# the optional CPPFLAGS at the end will override the previous flags
env['CPPFLAGS'] = env['CPPFLAGS'] + opt_flags
env.Append(CPPFLAGS=['-std=c++11'])

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
else:
    # if no module set (but at least one support language), select all modules
    anyModule = False
    for entry in moduleList.keys():
        if env[entry]:
            anyModule = True
    if not anyModule:
        print "Compiling all modules..."
        for entry in moduleList.keys():
            env[entry] = True

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
            env['SG_' + dep.upper()] = True

for modl in moduleList.keys():
    if env[modl]:
        print "Compiling module", modl
# support for non-C++
for sup in supportList:
    if env[sup]:
        print "Compiling support for", sup

# include Parallel and dependent modules only if OpenMP support is activated
# FIXME: it is actually a work around, the proper solution would involve change of source files. Afterwards this fix should be removed.
if env['SG_PARALLEL'] and not env['OMP']:
            print 'Warning: Building module Parallel requires OpenMP support. Please compile with OMP=1.'
            print 'Skipping modules Parallel, Misc, and Java support'
            env['SG_PARALLEL'] = False
            env['SG_MISC'] = False
            env['SG_JAVA'] = False

# add C++ defines for all modules
cppdefines = []
print moduleList
for modl in moduleList.keys():
    if env[modl]:
        cppdefines.append(modl)
print cppdefines
env.Append(CPPDEFINES=cppdefines)



# Initialize environment + support for Python and Java
#########################################################################

# no checks if clean:
if not env.GetOption('clean'):
    print ""
    print "******************************************"
    print "* Configuring system                     *"
    print "******************************************"
    log_file_orig = env['CMD_LOGFILE']
    env['CMD_LOGFILE'] = "config_commandline.log"

    config = env.Configure(custom_tests = { 'CheckExec' : CheckExec,
                                            'CheckJNI' : CheckJNI })

    # check scons
    EnsureSConsVersion(1, 0)
    print "Using SCons", SCons.__version__

    # check whether swig installed
    if not config.CheckExec('doxygen'):
        sys.stderr.write("Warning: doxygen cannot be found.\n  You will not be able to generate the documentation.\n  Check PATH environment variable!\n")

    # check whether dot installed
    if not config.CheckExec('dot'):
        sys.stderr.write("Warning: dot (Graphviz) cannot be found.\n  The documentation might lack diagrams.\n  Check PATH environment variable!\n")

    # check if the math header is available
    #if not config.CheckCXXHeader('cmath'):
    #    sys.stderr.write("Error: c++ math header cmath.h is missing.\n")
    #    Exit(1)

    # check whether swig installed
    swigAvail = True
    if not config.CheckExec('swig'):
        sys.stderr.write("Error: swig cannot be found. Check PATH environment variable!\n")
        swigAvail = False

    # check for Python headers
    pyAvail = True
    config.env.AppendUnique(CPPPATH = [distutils.sysconfig.get_python_inc()])
    if not config.CheckCXXHeader('Python.h'):
        sys.stderr.write("Error: Python.h not found. Check path to Python include files: "
                         + distutils.sysconfig.get_python_inc() + "\n")
        sys.stderr.write("Warning: You might have to install package python-dev\n")
        sys.stderr.write("... skipping Python support and unit tests")
        pyAvail = False
    else:
        numPyAvail = True
        # remove -Werror, if set. Elsewise, test will fail
        flagErrorRemoved = False
        if '-Werror' in config.env.get('CPPFLAGS'):
            config.env['CPPFLAGS'].remove('-Werror')
            flagErrorRemoved = True
        if not config.CheckCXXHeader(['pyconfig.h','Python.h','numpy/arrayobject.h']):
            try:
                print "... trying to extend path:"
                # get path to numpy header files
                import numpy
                numpy_path = os.path.join(os.path.split(numpy.__file__)[0],"core","include")
                if os.path.exists(numpy_path):
                    config.env.AppendUnique(CPPPATH = [numpy_path])
                    if not config.CheckCXXHeader(['pyconfig.h','Python.h','numpy/arrayobject.h']):
                        numPyAvail = False
                else:
                    sys.stderr.write("   Cannot find NumPy header files in:", numpy_path, "\n")
            except Exception, e:
                sys.stderr.write("   NumPy not available!\nException: %s\n" % e)
                numPyAvail = False
        if not numPyAvail:
            sys.stderr.write("   No NumPy support.\n   Corresponding unit tests and extended functionality are missing!\n")
        else:
            config.env.Append(NUMPY_AVAIL=1)
        # reappend -Werror if removed
        if flagErrorRemoved:
            config.env.Append(CPPFLAGS=['-Werror'])
    

    # check for $JAVA_HOME; prepend to search path
    javaAvail = True
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

    env = config.Finish()

    env['CMD_LOGFILE'] = log_file_orig
    print "******************************************"
    print "* Finished configuring system            *"
    print "******************************************"
else:
    swigAvail = True
    javaAvail = True
    pyAvail = True

# End of configuration
#########################################################################
Export('env')
Export('moduleList')


# Now compile
#########################################################################
lib_sgpp_targets = []
src_objs = {}                

# compile libraries
for name in modules:
    if env['SG_' + name.upper()]:
        print 'Building: ' + name
        name = name.lower()
        env.Append(CPPPATH=['#/src/sgpp'])
        
        # there is probably a more elgant way to do this
        for index in range(0, len(src_files[name])):
           src_files[name][index] = 'src/sgpp/' + src_files[name][index]

        src_objs[name] = env.SharedObject(src_files[name])

        lib = None
        
        # symbols cannot be undefined on windows (cygwin). Therefore static libraries with the 
        # symbols interdependent libraries have to be provided. Will not statically include the 
        # content of the static libraries. The static libraries work as import libraries.
        if env["PLATFORM"] == "cygwin":
            libdependencies = []
            # lib dependencies have to be ordered correctly!
            if name == "combigrid":
                libdependencies = ['sgppbasestatic']
            if name == "datadriven":
                libdependencies = ['sgppbasestatic', 'sgppsolverstatic', 'sgppmiscstatic', 'sgpppdestatic']
            elif name == "parallel":
                libdependencies = ['sgppdatadrivenstatic', 'sgppsolverstatic', 'sgppmiscstatic', 'sgppbasestatic']
            elif name == "pde":
                libdependencies = ['sgppbasestatic', 'sgppdatadrivenstatic', 'sgppsolverstatic']
            elif name == "solver":
                libdependencies = ['sgppbasestatic', 'sgppdatadrivenstatic', 'sgpppdestatic']
            elif name == "finance":
                libdependencies = ['sgppbasestatic', 'sgppdatadrivenstatic', 'sgpppdestatic', 'sgppsolverstatic']
            elif name == "misc":
                libdependencies = ['sgppbasestatic', 'sgppdatadrivenstatic']

            lib = env.SharedLibrary(target="tmp/build/sgpp" + name, source = src_objs[name], SHLIBPREFIX = 'lib', LIBS = libdependencies)
        else:
            lib = env.SharedLibrary(target="tmp/build/sgpp" + name, source = src_objs[name], SHLIBPREFIX = 'lib')
        # static libraries get the suffix "static" which allos scons to correctly resolve the dependencies 
        # of the shared libaries on the static libraries on windows
        libstatic = env.StaticLibrary(target="tmp/build/sgpp" + name + "static", source = src_objs[name], SHLIBPREFIX = 'lib')
        lib_sgpp_targets.append(lib)
        lib_sgpp_targets.append(libstatic)

Export('src_objs')

# build python lib
if env['SG_PYTHON'] and swigAvail and pyAvail:

    libpysgpp = SConscript('src/pysgpp/SConscript', variant_dir='tmp/build_pysgpp', duplicate=0)
    pyinst = env.Install(env['OUTPUT_PATH'] + 'lib/pysgpp', [libpysgpp, 'tmp/build_pysgpp/pysgpp.py'])
    Depends(pyinst, libpysgpp)
    pybin = env.Install(env['OUTPUT_PATH'] + 'bin', [libpysgpp, 'tmp/build_pysgpp/pysgpp.py'])
    Depends(pybin, libpysgpp)
    
# build java lib
if swigAvail and javaAvail and env['SG_JAVA']:
    libjsgpp = env.SConscript('src/jsgpp/SConscript',
                              variant_dir='tmp/build_jsgpp', duplicate=0)
#    libweka = env.SConscript('src/jsgpp_weka/SConscript',
#                             variant_dir='tmp/build_jsgpp_weka', duplicate=0)
    # install
    jinst = env.Install(env['OUTPUT_PATH'] + 'lib/jsgpp', [libjsgpp])
    
env.Install(env['OUTPUT_PATH'] + 'lib/sgpp', lib_sgpp_targets)

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
    sys.stderr.write("WARNING!! Skipping unit tests!!\n\n\n")
