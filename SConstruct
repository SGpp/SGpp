# Copyright (C) 2009 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp

# author Dirk Pflueger (Dirk.Pflueger@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), David Pfander (David.Pfander@ipvs.uni-stuttgart.de)

import glob
import SCons
import fnmatch
import os
import SGppConfigure

from SConsHelper import *

# Check for versions of Scons and Python
EnsurePythonVersion(2, 7)
# check scons
EnsureSConsVersion(2, 0)
print "Using SCons", SCons.__version__

# languageWrapperList = ['SG_PYTHON', 'SG_JAVA']
allLanguageWrapperList = ['SG_PYTHON']

# find all modules
moduleFolders = getModules('src/sgpp')

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

# no checks if clean:
if not env.GetOption('clean'):
    SGppConfigure.doConfigure(env)

# add C++ defines for all modules
cppdefines = []
for module in moduleNames:
    cppdefines.append(module)
env.Append(CPPDEFINES=cppdefines)

# environement setup finished, export environment
Export('env')
  
# compile selected modules
for name in moduleFolders:
        print "Preparing to build module: ", name
        # SConscript('src/sgpp/SConscript' + name, variant_dir='#/tmp/build/', duplicate=0)
        env.SConscript('src/sgpp/' + name + '/SConscript', {'env': env, 'moduleName': name})
       
# build python lib
if env['SG_PYTHON']:
    SConscript('src/pysgpp/SConscript')
    Import('pysgppInstall')

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
        Depends(testdep, [jinst, pysgppInstall])
    else:
        Depends(testdep, [pysgppInstall])
else:
    print "Warning: Skipping python unit tests"

