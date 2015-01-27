# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import glob
import SCons
import fnmatch
import os
import SGppConfigure

from Helper import *

# Check for versions of Scons and Python
EnsurePythonVersion(2, 7)
# check scons
EnsureSConsVersion(2, 0)
print "Using SCons", SCons.__version__

# languageWrapperList = ['SG_PYTHON', 'SG_JAVA']
allLanguageWrapperList = ['SG_PYTHON', 'SG_JAVA']

ignoreFolders = ['tests', 'jsgpp', 'pysgpp']

# find all modules
moduleFolders = getModules(ignoreFolders)
languageSupport = ['pysgpp']

prepareDoxyfile(moduleFolders)

moduleNames = []
for name in moduleFolders:
    moduleNames.append('SG_' + name.upper())

vars = Variables("custom.py")

# define the flags 
vars.Add('CPPFLAGS', 'Set additional Flags, they are compiler-depended (multiple flags combined with comma, e.g. -lpython,-lm)', '', converter=multiParamConverter)
vars.Add('LINKFLAGS', 'Set additional Linker-flags, they are linker-depended (multiple flags combined with comma, e.g. -lpython,-lm)', '', converter=multiParamConverter)
# define the target
vars.Add('MARCH', 'Sets the architecture if compiling with gcc, this is a pass-through option: just specify the gcc options!', None)
vars.Add('TARGETCPU', "Sets the processor you are compiling for. 'default' means using gcc with standard configuration. Also available are: 'ICC', here Intel Compiler in version 11 or higher must be used", 'default')
vars.Add(BoolVariable('OPT', "Sets optimization on and off", False))
# for compiling on LRZ without errors: omit unit tests
vars.Add(BoolVariable('NO_UNIT_TESTS', 'Omit UnitTests if set to True', False))
vars.Add(BoolVariable('SG_PYTHON', 'Build with python Support', True))
vars.Add(BoolVariable('SG_JAVA', 'Build with java Support', False))
vars.Add(BoolVariable('SSE3_FALLBACK', 'Tries to build as much as possible with SSE3 instead of AVX (intrinsics based functions won\'t work)', False))
vars.Add('OUTPUT_PATH', 'Path where built libraries are installed. Needs a trailing slash!', '')
vars.Add(BoolVariable('VERBOSE', 'Set output verbosity', False))
vars.Add('CMD_LOGFILE', 'Specifies a file to capture the build log', 'build.log')

# initialize environment
env = Environment(variables=vars, ENV=os.environ)
env.Export('moduleNames')
env.Export('moduleFolders')

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

# no checks if clean:
if not env.GetOption('clean'):
    SGppConfigure.doConfigure(env, moduleFolders)

# add C++ defines for all modules
cppdefines = []
for module in moduleNames:
    cppdefines.append(module)
env.Append(CPPDEFINES=cppdefines)

# environement setup finished, export environment
Export('env')

#Install alglib

libalglib, alglibstatic  = env.SConscript('tools/SConscriptAlglib', variant_dir='tmp/build_alglib', duplicate=0)
alglibinst = env.Install(env['OUTPUT_PATH'] + 'lib/alglib', [libalglib, alglibstatic])
# make base depend on alglib

env.Depends("#/" + BUILD_DIR.path + "/libsgppbase.so", alglibinst)
  
env.Append(CPPPATH=['#/tools'])
  
libraryTargetList = []
installTargetList = []
env.Export('libraryTargetList')
env.Export('installTargetList')
  
# compile selected modules
for name in moduleFolders:
    print "Preparing to build module: ", name
    # SConscript('src/sgpp/SConscript' + name, variant_dir='#/tmp/build/', duplicate=0)
    env.SConscript('#/' + name + '/SConscript', {'env': env, 'moduleName': name})
       
# build python lib
if env['SG_PYTHON']:
    SConscript('#/pysgpp/SConscript')
    Import('pysgppInstall')

# build java lib
if env['SG_JAVA']:
    libjsgpp = env.SConscript('#/jsgpp/SConscript',
                              variant_dir='tmp/build_jsgpp', duplicate=0)
    # install
    jinst = env.Install(env['OUTPUT_PATH'] + 'lib/jsgpp', [libjsgpp])
    
  
# Unit tests
#########################################################################
  
if not env['NO_UNIT_TESTS'] and env['SG_PYTHON']:
    
    env.Import('pysgppSimpleImportTest')
    
    # run tests
    builder = Builder(action = "python $SOURCE.file", chdir=1)
    env.Append(BUILDERS = {'Test' : builder})
    pysgppTestTargets = []
    dependency = None
    for moduleFolder in moduleFolders:
        if moduleFolder == "parallel":
        #    # these modules don't currently have tests
            continue
        moduleTest = env.Test('#/' + moduleFolder + '/tests/test_' + moduleFolder + '.py')
        #env.Requires(moduleTest, installTargetList)
        env.Requires(moduleTest, pysgppInstall)
        env.Depends(moduleTest, pysgppSimpleImportTest)
        env.AlwaysBuild(moduleTest) 
        if dependency == None:
            dependency = moduleTest
        else:
            env.Depends(moduleTest, dependency)
            dependency = moduleTest
            
        pysgppTestTargets.append(moduleTest)
else:
    print "Warning: Skipping python unit tests"
    
# not strictly necessary, used to enforce an order on the final steps of the building of the wrapper
# used to execute the unittests at the very end
if env['SG_PYTHON'] and env['SG_JAVA']:
    env.Requires(pysgppInstall, jsgppInstall)

